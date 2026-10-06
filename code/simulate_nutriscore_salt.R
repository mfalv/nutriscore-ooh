if (!requireNamespace("truncnorm", quietly = TRUE)) {
  personal_lib <- Sys.getenv("R_LIBS_USER", file.path(Sys.getenv("HOME"), "R", "library"))
  dir.create(personal_lib, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(personal_lib, .libPaths()))
  install.packages("truncnorm", lib = personal_lib, repos = "https://cloud.r-project.org")
}

library(truncnorm)
library(fst)
source("./global.R")

# Mozaffarian et al. NEJM (2014) - salt -> SBP effect function
# salt_difference: g NaCl/day (negative = reduction); 5.85 g NaCl ~ 100 mmol Na
# age, sbp, black_race: per-row vectors (data.table columns); n = .N; stochastic = TRUE in MC runs
salt_sbp_effect <- function(salt_difference, age, sbp, black_race, n, stochastic) {
  if (stochastic) {
    y <- rnorm(n, -3.735113,  0.7303861) +
      rnorm(n, -0.1052782, 0.0294371) * (age - 50) +
      rnorm(n, -1.873587,  0.8841412) * (sbp > 140) +
      rnorm(n, -2.489173,  1.188258)  * black_race
    effect <- -salt_difference * y / 5.85
    # Salt reduction -> SBP reduction only (clamp spurious sign-flips)
    effect <- pmin(effect, 0)
  } else {
    y      <- -3.735113 - 0.1052782 * (age - 50) -
      1.873587 * (sbp > 140) - 2.489173 * black_race
    effect <- -salt_difference * y / 5.85
  }
  return(effect)
}

IMPACTncd <- Simulation$new("sim_design_nutriscore_salt.yaml")

# Generate the 3 custom OOH exposures on-demand 
# ooh_freq_day/salt_eo_ooh/energykcal_eo_ooh
gen_ooh_exposures <- function(sp) {
  d <- IMPACTncd$design

  rs <- .Random.seed
  dqrs <- dqrng_get_state()
  set.seed(20250724L + sp$mc_aggr)
  dqset.seed(20250724L + sp$mc_aggr)
  sp$pop[, `:=`(
    rank_ooh_freq_day = dqrunif(1),
    rank_salt_eo_ooh = dqrunif(1),
    rank_energykcal_eo_ooh = dqrunif(1)
  ), by = pid]
  set.seed(rs)
  dqrng_set_state(dqrs)

  sp$pop[, `:=`(
    age_acc = age,
    age = pmin(pmax(age, 30L), 99L),
    bmi = as.integer(round(clamp(bmi_curr_xps, 14, 70), 0))
  )]
  d$exposures$ooh_freq_day$generate(sp$pop, d)
  d$exposures$salt_eo_ooh$generate(sp$pop, d)
  d$exposures$energykcal_eo_ooh$generate(sp$pop, d)
  sp$pop[, `:=`(age = age_acc, age_acc = NULL, bmi = NULL)]
  invisible(sp)
}
set_scn <- function(fn) {
  .current_scn_fn <<- fn
  IMPACTncd$update_primary_prevention_scn(function(sp) {
    gen_ooh_exposures(sp)
    .current_scn_fn(sp)
  })
}

# Load salt-energy Tweedie GLM (elasticity) model 
salt_energy_model <- readRDS("./inputs/salt_energy_model_final.rds")
tweedie_coefs     <- salt_energy_model$coefs

# Load Nutri-Score / calorie-labelling percentage effect sizes
ns_kcal_pct_reg <- readRDS("./inputs/ns_kcal_pct_reg.rds")

tweedie_predict_salt <- function(k) {
  unname(exp(tweedie_coefs["(Intercept)"] + tweedie_coefs["log_kcal"] * log(k)))
}

compute_delta_salt <- function(kcal0, delta_kcal) {
  kcal0     <- pmax(kcal0, 1e-6)
  kcal_post <- pmax(kcal0 + delta_kcal, 1e-6)  # numerical safety floor
  tweedie_predict_salt(kcal_post) - tweedie_predict_salt(kcal0)
}

cat("Salt-energy Tweedie GLM (elasticity): b1 =",
    round(unname(tweedie_coefs["log_kcal"]), 4),
    "(1% change in kcal -> ~", round(unname(tweedie_coefs["log_kcal"]), 4), "% change in salt)\n")
cat("  For -31.2 kcal/meal at mean baseline (", round(salt_energy_model$kcal_mean, 1), "): Δsalt ≈",
    round(compute_delta_salt(salt_energy_model$kcal_mean, -31.2), 5), "g/meal\n")

# Load Hall BMI lookup table (once, shared across all 8 scenarios)
change_bmi_Hall <- fst::read_fst(
  "./tables/Data_Hall_bmi_13_94_kcal_1_200.fst",
  as.data.table = TRUE
)
cat("Hall BMI lookup table loaded:", nrow(change_bmi_Hall),
    "rows (age 20-99, bmi 13-94, calories -1 to -200)\n")

n_runs <- 200L  # 50L test; change to 200L for publication run

# =============================================================================
# SC0 - Counterfactual (no OOH labelling, no structural salt policy)
# =============================================================================
IMPACTncd$
  del_logs()$
  del_outputs()$
  run(1:n_runs, multicore = TRUE, "sc0")


# =============================================================================
# SC1 - Nutri-Score | Large businesses (18%) | No structural salt policy
# =============================================================================
set_scn(
  function(sp) {
    set.seed(as.integer(sp$mc) + 20250L)

    # NS kcal draw (truncated normal, reductions only), no dietary compensation
    ns_pct   <- truncnorm::rtruncnorm(1, b = 0, mean = ns_kcal_pct_reg$ns$mean, sd = ns_kcal_pct_reg$ns$se)
    ns_kcal_meal   <- sp$pop$energykcal_eo_ooh * ns_pct
    coverage       <- 0.18
    # delta_kcal_meal is kcal per meal (coverage-adjusted).
    # ooh_freq_day is applied column-wise so each person uses their own GAMLSS-modelled daily OOH frequency.
    delta_kcal_meal <- ns_kcal_meal * coverage  # kcal/meal
    sc_year <- 26L

    # Pathway 1 - BMI
    # calories: person-specific daily kcal reduction = per-meal effect x meals/day
    # bmi_hall: baseline integer BMI
    sp$pop[, `:=`(
      bmi_hall = pmax(13L, pmin(94L, as.integer(bmi_curr_xps))),
      calories = pmax(-200L, as.integer(delta_kcal_meal * ooh_freq_day))
    )]
    sp$pop[is.na(calories) | calories > 0L, calories := 0L]
    for (yr in sc_year:max(sp$pop$year)) {
      syn_yr <- copy(sp$pop[year == yr,
                            .(age, sex, bmi = bmi_hall, calories, pid, bmi_curr_xps)])
      syn_yr <- merge(syn_yr, change_bmi_Hall,
                      by = c("age", "sex", "bmi", "calories"), all.x = TRUE)
      setorder(syn_yr, pid)
      syn_yr[is.na(change_bmi_last) | calories == 0L, change_bmi_last := 0]
      syn_yr[, bmi_curr_xps := bmi_curr_xps + change_bmi_last]
      sp$pop[year == yr, bmi_curr_xps := syn_yr$bmi_curr_xps]
    }
    sp$pop[, c("bmi_hall", "calories") := NULL]

    # Pathway 2 - SBP
    sp$pop[, delta_salt_meal := ifelse(energykcal_eo_ooh > 0,
      compute_delta_salt(energykcal_eo_ooh, ns_kcal_meal) * coverage, 0)]
    sp$pop[, black_race := as.integer(ethnicity %in% c("black african", "black caribbean"))]
    sp$pop[year >= sc_year,
           sbp_curr_xps := sbp_curr_xps + salt_sbp_effect(
             salt_difference = delta_salt_meal * ooh_freq_day,
             age        = age,
             sbp        = sbp_curr_xps,
             black_race = black_race,
             n          = .N,
             stochastic = TRUE
           )]

    sp$pop[, c("black_race", "delta_salt_meal") := NULL]
  }
)
IMPACTncd$run(1:n_runs, multicore = TRUE, "sc1_ns_large")


# =============================================================================
# SC2 - Nutri-Score | All businesses (100%) | No structural salt policy
# =============================================================================
set_scn(
  function(sp) {
    set.seed(as.integer(sp$mc) + 20250L)

    ns_pct   <- truncnorm::rtruncnorm(1, b = 0, mean = ns_kcal_pct_reg$ns$mean, sd = ns_kcal_pct_reg$ns$se)
    ns_kcal_meal   <- sp$pop$energykcal_eo_ooh * ns_pct
    coverage       <- 1.0
    delta_kcal_meal <- ns_kcal_meal * coverage  # kcal/meal
    sc_year <- 26L

    # Pathway 1 - BMI
    sp$pop[, `:=`(
      bmi_hall = pmax(13L, pmin(94L, as.integer(bmi_curr_xps))),
      calories = pmax(-200L, as.integer(delta_kcal_meal * ooh_freq_day))
    )]
    sp$pop[is.na(calories) | calories > 0L, calories := 0L]
    for (yr in sc_year:max(sp$pop$year)) {
      syn_yr <- copy(sp$pop[year == yr,
                            .(age, sex, bmi = bmi_hall, calories, pid, bmi_curr_xps)])
      syn_yr <- merge(syn_yr, change_bmi_Hall,
                      by = c("age", "sex", "bmi", "calories"), all.x = TRUE)
      setorder(syn_yr, pid)
      syn_yr[is.na(change_bmi_last) | calories == 0L, change_bmi_last := 0]
      syn_yr[, bmi_curr_xps := bmi_curr_xps + change_bmi_last]
      sp$pop[year == yr, bmi_curr_xps := syn_yr$bmi_curr_xps]
    }
    sp$pop[, c("bmi_hall", "calories") := NULL]

    # Pathway 2 - SBP 
    sp$pop[, delta_salt_meal := ifelse(energykcal_eo_ooh > 0,
      compute_delta_salt(energykcal_eo_ooh, ns_kcal_meal) * coverage, 0)]
    sp$pop[, black_race := as.integer(ethnicity %in% c("black african", "black caribbean"))]
    sp$pop[year >= sc_year,
           sbp_curr_xps := sbp_curr_xps + salt_sbp_effect(
             salt_difference = delta_salt_meal * ooh_freq_day,
             age        = age,
             sbp        = sbp_curr_xps,
             black_race = black_race,
             n          = .N,
             stochastic = TRUE
           )]

    sp$pop[, c("black_race", "delta_salt_meal") := NULL]
  }
)
IMPACTncd$run(1:n_runs, multicore = TRUE, "sc2_ns_all")


# =============================================================================
# SC3 - Energy + Nutri-Score | Large businesses (18%) | No structural salt policy
# =============================================================================
# Energy label adds a stochastic draw
set_scn(
  function(sp) {
    set.seed(as.integer(sp$mc) + 20250L)

    ns_pct     <- truncnorm::rtruncnorm(1, b = 0, mean = ns_kcal_pct_reg$ns$mean, sd = ns_kcal_pct_reg$ns$se)
    ns_kcal_meal     <- sp$pop$energykcal_eo_ooh * ns_pct
    energy_pct <- truncnorm::rtruncnorm(1, b = 0, mean = ns_kcal_pct_reg$energy_label$pct_mean, sd = ns_kcal_pct_reg$energy_label_se)
    energy_kcal_meal <- sp$pop$energykcal_eo_ooh * energy_pct
    kcal_meal        <- ns_kcal_meal + energy_kcal_meal
    coverage         <- 0.18
    delta_kcal_meal   <- kcal_meal * coverage  # kcal/meal
    sc_year <- 26L

    # Pathway 1 - BMI
    sp$pop[, `:=`(
      bmi_hall = pmax(13L, pmin(94L, as.integer(bmi_curr_xps))),
      calories = pmax(-200L, as.integer(delta_kcal_meal * ooh_freq_day))
    )]
    sp$pop[is.na(calories) | calories > 0L, calories := 0L]
    for (yr in sc_year:max(sp$pop$year)) {
      syn_yr <- copy(sp$pop[year == yr,
                            .(age, sex, bmi = bmi_hall, calories, pid, bmi_curr_xps)])
      syn_yr <- merge(syn_yr, change_bmi_Hall,
                      by = c("age", "sex", "bmi", "calories"), all.x = TRUE)
      setorder(syn_yr, pid)
      syn_yr[is.na(change_bmi_last) | calories == 0L, change_bmi_last := 0]
      syn_yr[, bmi_curr_xps := bmi_curr_xps + change_bmi_last]
      sp$pop[year == yr, bmi_curr_xps := syn_yr$bmi_curr_xps]
    }
    sp$pop[, c("bmi_hall", "calories") := NULL]

    # Pathway 2 - SBP 
    sp$pop[, delta_salt_meal := ifelse(energykcal_eo_ooh > 0,
      compute_delta_salt(energykcal_eo_ooh, kcal_meal) * coverage, 0)]
    sp$pop[, black_race := as.integer(ethnicity %in% c("black african", "black caribbean"))]
    sp$pop[year >= sc_year,
           sbp_curr_xps := sbp_curr_xps + salt_sbp_effect(
             salt_difference = delta_salt_meal * ooh_freq_day,
             age        = age,
             sbp        = sbp_curr_xps,
             black_race = black_race,
             n          = .N,
             stochastic = TRUE
           )]

    sp$pop[, c("black_race", "delta_salt_meal") := NULL]
  }
)
IMPACTncd$run(1:n_runs, multicore = TRUE, "sc3_energy_ns_large")


# =============================================================================
# SC4 - Energy + Nutri-Score | All businesses (100%) | No structural salt policy
# =============================================================================
set_scn(
  function(sp) {
    set.seed(as.integer(sp$mc) + 20250L)

    ns_pct     <- truncnorm::rtruncnorm(1, b = 0, mean = ns_kcal_pct_reg$ns$mean, sd = ns_kcal_pct_reg$ns$se)
    ns_kcal_meal     <- sp$pop$energykcal_eo_ooh * ns_pct
    energy_pct <- truncnorm::rtruncnorm(1, b = 0, mean = ns_kcal_pct_reg$energy_label$pct_mean, sd = ns_kcal_pct_reg$energy_label_se)
    energy_kcal_meal <- sp$pop$energykcal_eo_ooh * energy_pct
    kcal_meal        <- ns_kcal_meal + energy_kcal_meal
    coverage         <- 1.0
    delta_kcal_meal   <- kcal_meal * coverage  # kcal/meal
    sc_year <- 26L

    # Pathway 1 - BMI
    sp$pop[, `:=`(
      bmi_hall = pmax(13L, pmin(94L, as.integer(bmi_curr_xps))),
      calories = pmax(-200L, as.integer(delta_kcal_meal * ooh_freq_day))
    )]
    sp$pop[is.na(calories) | calories > 0L, calories := 0L]
    for (yr in sc_year:max(sp$pop$year)) {
      syn_yr <- copy(sp$pop[year == yr,
                            .(age, sex, bmi = bmi_hall, calories, pid, bmi_curr_xps)])
      syn_yr <- merge(syn_yr, change_bmi_Hall,
                      by = c("age", "sex", "bmi", "calories"), all.x = TRUE)
      setorder(syn_yr, pid)
      syn_yr[is.na(change_bmi_last) | calories == 0L, change_bmi_last := 0]
      syn_yr[, bmi_curr_xps := bmi_curr_xps + change_bmi_last]
      sp$pop[year == yr, bmi_curr_xps := syn_yr$bmi_curr_xps]
    }
    sp$pop[, c("bmi_hall", "calories") := NULL]

    # Pathway 2 - SBP
    sp$pop[, delta_salt_meal := ifelse(energykcal_eo_ooh > 0,
      compute_delta_salt(energykcal_eo_ooh, kcal_meal) * coverage, 0)]
    sp$pop[, black_race := as.integer(ethnicity %in% c("black african", "black caribbean"))]
    sp$pop[year >= sc_year,
           sbp_curr_xps := sbp_curr_xps + salt_sbp_effect(
             salt_difference = delta_salt_meal * ooh_freq_day,
             age        = age,
             sbp        = sbp_curr_xps,
             black_race = black_race,
             n          = .N,
             stochastic = TRUE
           )]

    sp$pop[, c("black_race", "delta_salt_meal") := NULL]
  }
)
IMPACTncd$run(1:n_runs, multicore = TRUE, "sc4_energy_ns_all")


# =============================================================================
# SC5 - Nutri-Score | Large businesses (18%) | FSA structural salt reduction
# =============================================================================
# Structural policy: 1.5%/yr cumulative reduction of OOH salt (always 100% sector).
#   OOH salt baseline: salt_eo_ooh (g/meal) x ooh_freq_day g/day - person-specific.
#   At policy year t (1-20): structural delta = -(salt_eo_ooh * ooh_freq_day) * 0.015 * t.
set_scn(
  function(sp) {
    set.seed(as.integer(sp$mc) + 20250L)

    ns_pct         <- truncnorm::rtruncnorm(1, b = 0, mean = ns_kcal_pct_reg$ns$mean, sd = ns_kcal_pct_reg$ns$se)
    ns_kcal_meal         <- sp$pop$energykcal_eo_ooh * ns_pct
    coverage             <- 0.18
    delta_kcal_meal       <- ns_kcal_meal * coverage  # kcal/meal
    sc_year <- 26L

    # Pathway 1 - BMI
    sp$pop[, `:=`(
      bmi_hall = pmax(13L, pmin(94L, as.integer(bmi_curr_xps))),
      calories = pmax(-200L, as.integer(delta_kcal_meal * ooh_freq_day))
    )]
    sp$pop[is.na(calories) | calories > 0L, calories := 0L]
    for (yr in sc_year:max(sp$pop$year)) {
      syn_yr <- copy(sp$pop[year == yr,
                            .(age, sex, bmi = bmi_hall, calories, pid, bmi_curr_xps)])
      syn_yr <- merge(syn_yr, change_bmi_Hall,
                      by = c("age", "sex", "bmi", "calories"), all.x = TRUE)
      setorder(syn_yr, pid)
      syn_yr[is.na(change_bmi_last) | calories == 0L, change_bmi_last := 0]
      syn_yr[, bmi_curr_xps := bmi_curr_xps + change_bmi_last]
      sp$pop[year == yr, bmi_curr_xps := syn_yr$bmi_curr_xps]
    }
    sp$pop[, c("bmi_hall", "calories") := NULL]

    # Pathway 2 - SBP
    sp$pop[, delta_salt_consumer := ifelse(energykcal_eo_ooh > 0,
      compute_delta_salt(energykcal_eo_ooh, ns_kcal_meal) * coverage, 0)]
    sp$pop[, black_race := as.integer(ethnicity %in% c("black african", "black caribbean"))]
    sp$pop[year >= sc_year,
           sbp_curr_xps := sbp_curr_xps + salt_sbp_effect(
             salt_difference = (delta_salt_consumer * ooh_freq_day) +
                               (-(salt_eo_ooh * ooh_freq_day) * 0.015 *
                                pmin(year - sc_year + 1L, 20L)),
             age        = age,
             sbp        = sbp_curr_xps,
             black_race = black_race,
             n          = .N,
             stochastic = TRUE
           )]

    sp$pop[, c("black_race", "delta_salt_consumer") := NULL]
  }
)
IMPACTncd$run(1:n_runs, multicore = TRUE, "sc5_ns_fsa_large")


# =============================================================================
# SC6 - Nutri-Score | All businesses (100%) | FSA structural salt reduction
# =============================================================================
set_scn(
  function(sp) {
    set.seed(as.integer(sp$mc) + 20250L)

    ns_pct        <- truncnorm::rtruncnorm(1, b = 0, mean = ns_kcal_pct_reg$ns$mean, sd = ns_kcal_pct_reg$ns$se)
    ns_kcal_meal        <- sp$pop$energykcal_eo_ooh * ns_pct
    coverage            <- 1.0
    delta_kcal_meal      <- ns_kcal_meal * coverage  # kcal/meal
    sc_year <- 26L

    # Pathway 1 - BMI
    sp$pop[, `:=`(
      bmi_hall = pmax(13L, pmin(94L, as.integer(bmi_curr_xps))),
      calories = pmax(-200L, as.integer(delta_kcal_meal * ooh_freq_day))
    )]
    sp$pop[is.na(calories) | calories > 0L, calories := 0L]
    for (yr in sc_year:max(sp$pop$year)) {
      syn_yr <- copy(sp$pop[year == yr,
                            .(age, sex, bmi = bmi_hall, calories, pid, bmi_curr_xps)])
      syn_yr <- merge(syn_yr, change_bmi_Hall,
                      by = c("age", "sex", "bmi", "calories"), all.x = TRUE)
      setorder(syn_yr, pid)
      syn_yr[is.na(change_bmi_last) | calories == 0L, change_bmi_last := 0]
      syn_yr[, bmi_curr_xps := bmi_curr_xps + change_bmi_last]
      sp$pop[year == yr, bmi_curr_xps := syn_yr$bmi_curr_xps]
    }
    sp$pop[, c("bmi_hall", "calories") := NULL]

    # Pathway 2 - SBP
    sp$pop[, delta_salt_consumer := ifelse(energykcal_eo_ooh > 0,
      compute_delta_salt(energykcal_eo_ooh, ns_kcal_meal) * coverage, 0)]
    sp$pop[, black_race := as.integer(ethnicity %in% c("black african", "black caribbean"))]
    sp$pop[year >= sc_year,
           sbp_curr_xps := sbp_curr_xps + salt_sbp_effect(
             salt_difference = (delta_salt_consumer * ooh_freq_day) +
                               (-(salt_eo_ooh * ooh_freq_day) * 0.015 *
                                pmin(year - sc_year + 1L, 20L)),
             age        = age,
             sbp        = sbp_curr_xps,
             black_race = black_race,
             n          = .N,
             stochastic = TRUE
           )]

    sp$pop[, c("black_race", "delta_salt_consumer") := NULL]
  }
)
IMPACTncd$run(1:n_runs, multicore = TRUE, "sc6_ns_fsa_all")


# =============================================================================
# SC7 - Energy + Nutri-Score | Large businesses (18%) | FSA structural salt
# =============================================================================
set_scn(
  function(sp) {
    set.seed(as.integer(sp$mc) + 20250L)

    ns_pct        <- truncnorm::rtruncnorm(1, b = 0, mean = ns_kcal_pct_reg$ns$mean, sd = ns_kcal_pct_reg$ns$se)
    ns_kcal_meal        <- sp$pop$energykcal_eo_ooh * ns_pct
    energy_pct    <- truncnorm::rtruncnorm(1, b = 0, mean = ns_kcal_pct_reg$energy_label$pct_mean, sd = ns_kcal_pct_reg$energy_label_se)
    energy_kcal_meal    <- sp$pop$energykcal_eo_ooh * energy_pct
    kcal_meal           <- ns_kcal_meal + energy_kcal_meal
    coverage            <- 0.18
    delta_kcal_meal      <- kcal_meal * coverage  # kcal/meal
    sc_year <- 26L

    # Pathway 1 - BMI
    sp$pop[, `:=`(
      bmi_hall = pmax(13L, pmin(94L, as.integer(bmi_curr_xps))),
      calories = pmax(-200L, as.integer(delta_kcal_meal * ooh_freq_day))
    )]
    sp$pop[is.na(calories) | calories > 0L, calories := 0L]
    for (yr in sc_year:max(sp$pop$year)) {
      syn_yr <- copy(sp$pop[year == yr,
                            .(age, sex, bmi = bmi_hall, calories, pid, bmi_curr_xps)])
      syn_yr <- merge(syn_yr, change_bmi_Hall,
                      by = c("age", "sex", "bmi", "calories"), all.x = TRUE)
      setorder(syn_yr, pid)
      syn_yr[is.na(change_bmi_last) | calories == 0L, change_bmi_last := 0]
      syn_yr[, bmi_curr_xps := bmi_curr_xps + change_bmi_last]
      sp$pop[year == yr, bmi_curr_xps := syn_yr$bmi_curr_xps]
    }
    sp$pop[, c("bmi_hall", "calories") := NULL]

    # Pathway 2 - SBP
    sp$pop[, delta_salt_consumer := ifelse(energykcal_eo_ooh > 0,
      compute_delta_salt(energykcal_eo_ooh, kcal_meal) * coverage, 0)]
    sp$pop[, black_race := as.integer(ethnicity %in% c("black african", "black caribbean"))]
    sp$pop[year >= sc_year,
           sbp_curr_xps := sbp_curr_xps + salt_sbp_effect(
             salt_difference = (delta_salt_consumer * ooh_freq_day) +
                               (-(salt_eo_ooh * ooh_freq_day) * 0.015 *
                                pmin(year - sc_year + 1L, 20L)),
             age        = age,
             sbp        = sbp_curr_xps,
             black_race = black_race,
             n          = .N,
             stochastic = TRUE
           )]

    sp$pop[, c("black_race", "delta_salt_consumer") := NULL]
  }
)
IMPACTncd$run(1:n_runs, multicore = TRUE, "sc7_energy_ns_fsa_large")


# =============================================================================
# SC8 - Energy + Nutri-Score | All businesses (100%) | FSA structural salt
# =============================================================================
set_scn(
  function(sp) {
    set.seed(as.integer(sp$mc) + 20250L)

    ns_pct        <- truncnorm::rtruncnorm(1, b = 0, mean = ns_kcal_pct_reg$ns$mean, sd = ns_kcal_pct_reg$ns$se)
    ns_kcal_meal        <- sp$pop$energykcal_eo_ooh * ns_pct
    energy_pct    <- truncnorm::rtruncnorm(1, b = 0, mean = ns_kcal_pct_reg$energy_label$pct_mean, sd = ns_kcal_pct_reg$energy_label_se)
    energy_kcal_meal    <- sp$pop$energykcal_eo_ooh * energy_pct
    kcal_meal           <- ns_kcal_meal + energy_kcal_meal
    coverage            <- 1.0
    delta_kcal_meal      <- kcal_meal * coverage  # kcal/meal
    sc_year <- 26L

    # Pathway 1 - BMI
    sp$pop[, `:=`(
      bmi_hall = pmax(13L, pmin(94L, as.integer(bmi_curr_xps))),
      calories = pmax(-200L, as.integer(delta_kcal_meal * ooh_freq_day))
    )]
    sp$pop[is.na(calories) | calories > 0L, calories := 0L]
    for (yr in sc_year:max(sp$pop$year)) {
      syn_yr <- copy(sp$pop[year == yr,
                            .(age, sex, bmi = bmi_hall, calories, pid, bmi_curr_xps)])
      syn_yr <- merge(syn_yr, change_bmi_Hall,
                      by = c("age", "sex", "bmi", "calories"), all.x = TRUE)
      setorder(syn_yr, pid)
      syn_yr[is.na(change_bmi_last) | calories == 0L, change_bmi_last := 0]
      syn_yr[, bmi_curr_xps := bmi_curr_xps + change_bmi_last]
      sp$pop[year == yr, bmi_curr_xps := syn_yr$bmi_curr_xps]
    }
    sp$pop[, c("bmi_hall", "calories") := NULL]

    # Pathway 2 - SBP
    sp$pop[, delta_salt_consumer := ifelse(energykcal_eo_ooh > 0,
      compute_delta_salt(energykcal_eo_ooh, kcal_meal) * coverage, 0)]
    sp$pop[, black_race := as.integer(ethnicity %in% c("black african", "black caribbean"))]
    sp$pop[year >= sc_year,
           sbp_curr_xps := sbp_curr_xps + salt_sbp_effect(
             salt_difference = (delta_salt_consumer * ooh_freq_day) +
                               (-(salt_eo_ooh * ooh_freq_day) * 0.015 *
                                pmin(year - sc_year + 1L, 20L)),
             age        = age,
             sbp        = sbp_curr_xps,
             black_race = black_race,
             n          = .N,
             stochastic = TRUE
           )]

    sp$pop[, c("black_race", "delta_salt_consumer") := NULL]
  }
)
IMPACTncd$run(1:n_runs, multicore = TRUE, "sc8_energy_ns_fsa_all")


# =============================================================================
# SENSITIVITY ANALYSES
# =============================================================================
# SA1  Dietary compensation (11%, 26.5%, 42%)            - SC2, SC4, SC6, SC8
# SA2  Industry response to mandatory calorie labelling  - SC3, SC4, SC7, SC8
# SA3  Calorie + traffic-light labelling                 - SC3, SC4, SC7, SC8
# SA4  Nutri-Score coverage at 47% turnover              - SC1, SC3, SC5, SC7
#
# All SA scenarios share the same seed as their base scenario so that the NS
# and energy-label draws are identical, enabling clean paired comparisons.
# =============================================================================


# SA1: Dietary compensation
# Compensation fraction comp reduces the net calorie deficit by (1 - comp). 
# Salt-reduction effect and FSA structural salt is NOT compensated.
# Coverage = 1.0 (all-business scenarios only: SC2, SC4, SC6, SC8).
.build_comp_scn <- function(energy_label, fsa, comp) {
  force(energy_label); force(fsa); force(comp)
  function(sp) {
    set.seed(as.integer(sp$mc) + 20250L)

    ns_pct <- truncnorm::rtruncnorm(1, b = 0, mean = ns_kcal_pct_reg$ns$mean, sd = ns_kcal_pct_reg$ns$se)
    ns_kcal_meal <- sp$pop$energykcal_eo_ooh * ns_pct
    if (energy_label) {
      energy_pct <- truncnorm::rtruncnorm(1, b = 0, mean = ns_kcal_pct_reg$energy_label$pct_mean, sd = ns_kcal_pct_reg$energy_label_se)
      energy_kcal_meal <- sp$pop$energykcal_eo_ooh * energy_pct
      kcal_meal <- ns_kcal_meal + energy_kcal_meal
    } else {
      kcal_meal <- ns_kcal_meal
    }

    delta_kcal_meal <- kcal_meal * (1 - comp)   # BMI path: compensated
    sc_year <- 26L

    sp$pop[, `:=`(
      bmi_hall = pmax(13L, pmin(94L, as.integer(bmi_curr_xps))),
      calories = pmax(-200L, as.integer(delta_kcal_meal * ooh_freq_day))
    )]
    sp$pop[is.na(calories) | calories > 0L, calories := 0L]
    for (yr in sc_year:max(sp$pop$year)) {
      syn_yr <- copy(sp$pop[year == yr,
                            .(age, sex, bmi = bmi_hall, calories, pid, bmi_curr_xps)])
      syn_yr <- merge(syn_yr, change_bmi_Hall,
                      by = c("age", "sex", "bmi", "calories"), all.x = TRUE)
      setorder(syn_yr, pid)
      syn_yr[is.na(change_bmi_last) | calories == 0L, change_bmi_last := 0]
      syn_yr[, bmi_curr_xps := bmi_curr_xps + change_bmi_last]
      sp$pop[year == yr, bmi_curr_xps := syn_yr$bmi_curr_xps]
    }
    sp$pop[, c("bmi_hall", "calories") := NULL]

    # Salt path:
    sp$pop[, delta_salt_consumer := ifelse(energykcal_eo_ooh > 0,
      compute_delta_salt(energykcal_eo_ooh, kcal_meal), 0)]
    sp$pop[, black_race := as.integer(ethnicity %in% c("black african", "black caribbean"))]
    if (fsa) {
      sp$pop[year >= sc_year,
             sbp_curr_xps := sbp_curr_xps + salt_sbp_effect(
               salt_difference = (delta_salt_consumer * ooh_freq_day) +
                 (-(salt_eo_ooh * ooh_freq_day) * 0.015 * pmin(year - sc_year + 1L, 20L)),
               age = age, sbp = sbp_curr_xps, black_race = black_race,
               n = .N, stochastic = TRUE)]
    } else {
      sp$pop[year >= sc_year,
             sbp_curr_xps := sbp_curr_xps + salt_sbp_effect(
               salt_difference = delta_salt_consumer * ooh_freq_day,
               age = age, sbp = sbp_curr_xps, black_race = black_race,
               n = .N, stochastic = TRUE)]
    }
    sp$pop[, c("black_race", "delta_salt_consumer") := NULL]
  }
}

.comp_levels <- list(
  list(val = 0.11,  str = "11"),
  list(val = 0.265, str = "265"),
  list(val = 0.42,  str = "42")
)

for (.cl in .comp_levels) {
  # SC2 analogue: NS only, all businesses, no FSA
  set_scn(.build_comp_scn(FALSE, FALSE, .cl$val))
  IMPACTncd$run(1:n_runs, multicore = TRUE, paste0("sc2_ns_all_sa_comp",             .cl$str))
  # SC4 analogue: energy+NS, all businesses, no FSA
  set_scn(.build_comp_scn(TRUE,  FALSE, .cl$val))
  IMPACTncd$run(1:n_runs, multicore = TRUE, paste0("sc4_energy_ns_all_sa_comp",      .cl$str))
  # SC6 analogue: NS only, all businesses, FSA
  set_scn(.build_comp_scn(FALSE, TRUE,  .cl$val))
  IMPACTncd$run(1:n_runs, multicore = TRUE, paste0("sc6_ns_fsa_all_sa_comp",         .cl$str))
  # SC8 analogue: energy+NS, all businesses, FSA
  set_scn(.build_comp_scn(TRUE,  TRUE,  .cl$val))
  IMPACTncd$run(1:n_runs, multicore = TRUE, paste0("sc8_energy_ns_fsa_all_sa_comp",  .cl$str))
}
rm(.build_comp_scn, .comp_levels, .cl)


# SA2: Industry response to mandatory calorie labelling (Essman et al. 2025)
# SA2-SC3: energy+NS, large (18%), no FSA
set_scn(
  function(sp) {
    set.seed(as.integer(sp$mc) + 20250L)

    ns_pct       <- truncnorm::rtruncnorm(1, b = 0, mean = ns_kcal_pct_reg$ns$mean, sd = ns_kcal_pct_reg$ns$se)
    ns_kcal_meal       <- sp$pop$energykcal_eo_ooh * ns_pct
    energy_pct   <- truncnorm::rtruncnorm(1, b = 0, mean = ns_kcal_pct_reg$energy_label$pct_mean, sd = ns_kcal_pct_reg$energy_label_se)
    energy_kcal_meal   <- sp$pop$energykcal_eo_ooh * energy_pct
    industry_pct <- truncnorm::rtruncnorm(1, b = 0, mean = ns_kcal_pct_reg$industry$pct_mean, sd = ns_kcal_pct_reg$industry_se)
    industry_kcal_meal <- sp$pop$energykcal_eo_ooh * industry_pct
    coverage           <- 0.18
    delta_kcal_meal     <- (ns_kcal_meal + energy_kcal_meal + industry_kcal_meal) * coverage
    sc_year <- 26L

    sp$pop[, `:=`(
      bmi_hall = pmax(13L, pmin(94L, as.integer(bmi_curr_xps))),
      calories = pmax(-200L, as.integer(delta_kcal_meal * ooh_freq_day))
    )]
    sp$pop[is.na(calories) | calories > 0L, calories := 0L]
    for (yr in sc_year:max(sp$pop$year)) {
      syn_yr <- copy(sp$pop[year == yr,
                            .(age, sex, bmi = bmi_hall, calories, pid, bmi_curr_xps)])
      syn_yr <- merge(syn_yr, change_bmi_Hall,
                      by = c("age", "sex", "bmi", "calories"), all.x = TRUE)
      setorder(syn_yr, pid)
      syn_yr[is.na(change_bmi_last) | calories == 0L, change_bmi_last := 0]
      syn_yr[, bmi_curr_xps := bmi_curr_xps + change_bmi_last]
      sp$pop[year == yr, bmi_curr_xps := syn_yr$bmi_curr_xps]
    }
    sp$pop[, c("bmi_hall", "calories") := NULL]

    sp$pop[, delta_salt_meal := ifelse(energykcal_eo_ooh > 0,
      compute_delta_salt(energykcal_eo_ooh, (ns_kcal_meal + energy_kcal_meal + industry_kcal_meal)) * coverage, 0)]
    sp$pop[, black_race := as.integer(ethnicity %in% c("black african", "black caribbean"))]
    sp$pop[year >= sc_year,
           sbp_curr_xps := sbp_curr_xps + salt_sbp_effect(
             salt_difference = delta_salt_meal * ooh_freq_day,
             age = age, sbp = sbp_curr_xps, black_race = black_race,
             n = .N, stochastic = TRUE)]
    sp$pop[, c("black_race", "delta_salt_meal") := NULL]
  }
)
IMPACTncd$run(1:n_runs, multicore = TRUE, "sc3_energy_ns_large_sa_essman")

# SA2-SC4: energy+NS, all (100%), no FSA
set_scn(
  function(sp) {
    set.seed(as.integer(sp$mc) + 20250L)

    ns_pct       <- truncnorm::rtruncnorm(1, b = 0, mean = ns_kcal_pct_reg$ns$mean, sd = ns_kcal_pct_reg$ns$se)
    ns_kcal_meal       <- sp$pop$energykcal_eo_ooh * ns_pct
    energy_pct   <- truncnorm::rtruncnorm(1, b = 0, mean = ns_kcal_pct_reg$energy_label$pct_mean, sd = ns_kcal_pct_reg$energy_label_se)
    energy_kcal_meal   <- sp$pop$energykcal_eo_ooh * energy_pct
    industry_pct <- truncnorm::rtruncnorm(1, b = 0, mean = ns_kcal_pct_reg$industry$pct_mean, sd = ns_kcal_pct_reg$industry_se)
    industry_kcal_meal <- sp$pop$energykcal_eo_ooh * industry_pct
    coverage           <- 1.0
    delta_kcal_meal     <- (ns_kcal_meal + energy_kcal_meal + industry_kcal_meal) * coverage
    sc_year <- 26L

    sp$pop[, `:=`(
      bmi_hall = pmax(13L, pmin(94L, as.integer(bmi_curr_xps))),
      calories = pmax(-200L, as.integer(delta_kcal_meal * ooh_freq_day))
    )]
    sp$pop[is.na(calories) | calories > 0L, calories := 0L]
    for (yr in sc_year:max(sp$pop$year)) {
      syn_yr <- copy(sp$pop[year == yr,
                            .(age, sex, bmi = bmi_hall, calories, pid, bmi_curr_xps)])
      syn_yr <- merge(syn_yr, change_bmi_Hall,
                      by = c("age", "sex", "bmi", "calories"), all.x = TRUE)
      setorder(syn_yr, pid)
      syn_yr[is.na(change_bmi_last) | calories == 0L, change_bmi_last := 0]
      syn_yr[, bmi_curr_xps := bmi_curr_xps + change_bmi_last]
      sp$pop[year == yr, bmi_curr_xps := syn_yr$bmi_curr_xps]
    }
    sp$pop[, c("bmi_hall", "calories") := NULL]

    sp$pop[, delta_salt_meal := ifelse(energykcal_eo_ooh > 0,
      compute_delta_salt(energykcal_eo_ooh, (ns_kcal_meal + energy_kcal_meal + industry_kcal_meal)) * coverage, 0)]
    sp$pop[, black_race := as.integer(ethnicity %in% c("black african", "black caribbean"))]
    sp$pop[year >= sc_year,
           sbp_curr_xps := sbp_curr_xps + salt_sbp_effect(
             salt_difference = delta_salt_meal * ooh_freq_day,
             age = age, sbp = sbp_curr_xps, black_race = black_race,
             n = .N, stochastic = TRUE)]
    sp$pop[, c("black_race", "delta_salt_meal") := NULL]
  }
)
IMPACTncd$run(1:n_runs, multicore = TRUE, "sc4_energy_ns_all_sa_essman")

# SA2-SC7: energy+NS, large (18%), FSA structural salt
set_scn(
  function(sp) {
    set.seed(as.integer(sp$mc) + 20250L)

    ns_pct       <- truncnorm::rtruncnorm(1, b = 0, mean = ns_kcal_pct_reg$ns$mean, sd = ns_kcal_pct_reg$ns$se)
    ns_kcal_meal       <- sp$pop$energykcal_eo_ooh * ns_pct
    energy_pct   <- truncnorm::rtruncnorm(1, b = 0, mean = ns_kcal_pct_reg$energy_label$pct_mean, sd = ns_kcal_pct_reg$energy_label_se)
    energy_kcal_meal   <- sp$pop$energykcal_eo_ooh * energy_pct
    industry_pct <- truncnorm::rtruncnorm(1, b = 0, mean = ns_kcal_pct_reg$industry$pct_mean, sd = ns_kcal_pct_reg$industry_se)
    industry_kcal_meal <- sp$pop$energykcal_eo_ooh * industry_pct
    coverage           <- 0.18
    delta_kcal_meal     <- (ns_kcal_meal + energy_kcal_meal + industry_kcal_meal) * coverage
    sc_year <- 26L

    sp$pop[, `:=`(
      bmi_hall = pmax(13L, pmin(94L, as.integer(bmi_curr_xps))),
      calories = pmax(-200L, as.integer(delta_kcal_meal * ooh_freq_day))
    )]
    sp$pop[is.na(calories) | calories > 0L, calories := 0L]
    for (yr in sc_year:max(sp$pop$year)) {
      syn_yr <- copy(sp$pop[year == yr,
                            .(age, sex, bmi = bmi_hall, calories, pid, bmi_curr_xps)])
      syn_yr <- merge(syn_yr, change_bmi_Hall,
                      by = c("age", "sex", "bmi", "calories"), all.x = TRUE)
      setorder(syn_yr, pid)
      syn_yr[is.na(change_bmi_last) | calories == 0L, change_bmi_last := 0]
      syn_yr[, bmi_curr_xps := bmi_curr_xps + change_bmi_last]
      sp$pop[year == yr, bmi_curr_xps := syn_yr$bmi_curr_xps]
    }
    sp$pop[, c("bmi_hall", "calories") := NULL]

    sp$pop[, delta_salt_consumer := ifelse(energykcal_eo_ooh > 0,
      compute_delta_salt(energykcal_eo_ooh, (ns_kcal_meal + energy_kcal_meal + industry_kcal_meal)) * coverage, 0)]
    sp$pop[, black_race := as.integer(ethnicity %in% c("black african", "black caribbean"))]
    sp$pop[year >= sc_year,
           sbp_curr_xps := sbp_curr_xps + salt_sbp_effect(
             salt_difference = (delta_salt_consumer * ooh_freq_day) +
               (-(salt_eo_ooh * ooh_freq_day) * 0.015 * pmin(year - sc_year + 1L, 20L)),
             age = age, sbp = sbp_curr_xps, black_race = black_race,
             n = .N, stochastic = TRUE)]
    sp$pop[, c("black_race", "delta_salt_consumer") := NULL]
  }
)
IMPACTncd$run(1:n_runs, multicore = TRUE, "sc7_energy_ns_fsa_large_sa_essman")

# SA2-SC8: energy+NS, all (100%), FSA structural salt
set_scn(
  function(sp) {
    set.seed(as.integer(sp$mc) + 20250L)

    ns_pct       <- truncnorm::rtruncnorm(1, b = 0, mean = ns_kcal_pct_reg$ns$mean, sd = ns_kcal_pct_reg$ns$se)
    ns_kcal_meal       <- sp$pop$energykcal_eo_ooh * ns_pct
    energy_pct   <- truncnorm::rtruncnorm(1, b = 0, mean = ns_kcal_pct_reg$energy_label$pct_mean, sd = ns_kcal_pct_reg$energy_label_se)
    energy_kcal_meal   <- sp$pop$energykcal_eo_ooh * energy_pct
    industry_pct <- truncnorm::rtruncnorm(1, b = 0, mean = ns_kcal_pct_reg$industry$pct_mean, sd = ns_kcal_pct_reg$industry_se)
    industry_kcal_meal <- sp$pop$energykcal_eo_ooh * industry_pct
    coverage           <- 1.0
    delta_kcal_meal     <- (ns_kcal_meal + energy_kcal_meal + industry_kcal_meal) * coverage
    sc_year <- 26L

    sp$pop[, `:=`(
      bmi_hall = pmax(13L, pmin(94L, as.integer(bmi_curr_xps))),
      calories = pmax(-200L, as.integer(delta_kcal_meal * ooh_freq_day))
    )]
    sp$pop[is.na(calories) | calories > 0L, calories := 0L]
    for (yr in sc_year:max(sp$pop$year)) {
      syn_yr <- copy(sp$pop[year == yr,
                            .(age, sex, bmi = bmi_hall, calories, pid, bmi_curr_xps)])
      syn_yr <- merge(syn_yr, change_bmi_Hall,
                      by = c("age", "sex", "bmi", "calories"), all.x = TRUE)
      setorder(syn_yr, pid)
      syn_yr[is.na(change_bmi_last) | calories == 0L, change_bmi_last := 0]
      syn_yr[, bmi_curr_xps := bmi_curr_xps + change_bmi_last]
      sp$pop[year == yr, bmi_curr_xps := syn_yr$bmi_curr_xps]
    }
    sp$pop[, c("bmi_hall", "calories") := NULL]

    sp$pop[, delta_salt_consumer := ifelse(energykcal_eo_ooh > 0,
      compute_delta_salt(energykcal_eo_ooh, (ns_kcal_meal + energy_kcal_meal + industry_kcal_meal)) * coverage, 0)]
    sp$pop[, black_race := as.integer(ethnicity %in% c("black african", "black caribbean"))]
    sp$pop[year >= sc_year,
           sbp_curr_xps := sbp_curr_xps + salt_sbp_effect(
             salt_difference = (delta_salt_consumer * ooh_freq_day) +
               (-(salt_eo_ooh * ooh_freq_day) * 0.015 * pmin(year - sc_year + 1L, 20L)),
             age = age, sbp = sbp_curr_xps, black_race = black_race,
             n = .N, stochastic = TRUE)]
    sp$pop[, c("black_race", "delta_salt_consumer") := NULL]
  }
)
IMPACTncd$run(1:n_runs, multicore = TRUE, "sc8_energy_ns_fsa_all_sa_essman")


# SA3: Calorie + traffic-light labelling (Ellison et al. 2014)
# SA3-SC3: energy+NS, large (18%), no FSA - Ellison replaces consumer draw
set_scn(
  function(sp) {
    set.seed(as.integer(sp$mc) + 20250L)

    ellison_pct <- truncnorm::rtruncnorm(1, b = 0, mean = ns_kcal_pct_reg$ellison$pct_mean, sd = ns_kcal_pct_reg$ellison_se)
    ellison_kcal_meal <- sp$pop$energykcal_eo_ooh * ellison_pct
    coverage          <- 0.18
    delta_kcal_meal    <- ellison_kcal_meal * coverage
    sc_year <- 26L

    sp$pop[, `:=`(
      bmi_hall = pmax(13L, pmin(94L, as.integer(bmi_curr_xps))),
      calories = pmax(-200L, as.integer(delta_kcal_meal * ooh_freq_day))
    )]
    sp$pop[is.na(calories) | calories > 0L, calories := 0L]
    for (yr in sc_year:max(sp$pop$year)) {
      syn_yr <- copy(sp$pop[year == yr,
                            .(age, sex, bmi = bmi_hall, calories, pid, bmi_curr_xps)])
      syn_yr <- merge(syn_yr, change_bmi_Hall,
                      by = c("age", "sex", "bmi", "calories"), all.x = TRUE)
      setorder(syn_yr, pid)
      syn_yr[is.na(change_bmi_last) | calories == 0L, change_bmi_last := 0]
      syn_yr[, bmi_curr_xps := bmi_curr_xps + change_bmi_last]
      sp$pop[year == yr, bmi_curr_xps := syn_yr$bmi_curr_xps]
    }
    sp$pop[, c("bmi_hall", "calories") := NULL]

    sp$pop[, delta_salt_meal := ifelse(energykcal_eo_ooh > 0,
      compute_delta_salt(energykcal_eo_ooh, ellison_kcal_meal) * coverage, 0)]
    sp$pop[, black_race := as.integer(ethnicity %in% c("black african", "black caribbean"))]
    sp$pop[year >= sc_year,
           sbp_curr_xps := sbp_curr_xps + salt_sbp_effect(
             salt_difference = delta_salt_meal * ooh_freq_day,
             age = age, sbp = sbp_curr_xps, black_race = black_race,
             n = .N, stochastic = TRUE)]
    sp$pop[, c("black_race", "delta_salt_meal") := NULL]
  }
)
IMPACTncd$run(1:n_runs, multicore = TRUE, "sc3_energy_ns_large_sa_ellison")

# SA3-SC4: energy+NS, all (100%), no FSA
set_scn(
  function(sp) {
    set.seed(as.integer(sp$mc) + 20250L)

    ellison_pct <- truncnorm::rtruncnorm(1, b = 0, mean = ns_kcal_pct_reg$ellison$pct_mean, sd = ns_kcal_pct_reg$ellison_se)
    ellison_kcal_meal <- sp$pop$energykcal_eo_ooh * ellison_pct
    coverage          <- 1.0
    delta_kcal_meal    <- ellison_kcal_meal * coverage
    sc_year <- 26L

    sp$pop[, `:=`(
      bmi_hall = pmax(13L, pmin(94L, as.integer(bmi_curr_xps))),
      calories = pmax(-200L, as.integer(delta_kcal_meal * ooh_freq_day))
    )]
    sp$pop[is.na(calories) | calories > 0L, calories := 0L]
    for (yr in sc_year:max(sp$pop$year)) {
      syn_yr <- copy(sp$pop[year == yr,
                            .(age, sex, bmi = bmi_hall, calories, pid, bmi_curr_xps)])
      syn_yr <- merge(syn_yr, change_bmi_Hall,
                      by = c("age", "sex", "bmi", "calories"), all.x = TRUE)
      setorder(syn_yr, pid)
      syn_yr[is.na(change_bmi_last) | calories == 0L, change_bmi_last := 0]
      syn_yr[, bmi_curr_xps := bmi_curr_xps + change_bmi_last]
      sp$pop[year == yr, bmi_curr_xps := syn_yr$bmi_curr_xps]
    }
    sp$pop[, c("bmi_hall", "calories") := NULL]

    sp$pop[, delta_salt_meal := ifelse(energykcal_eo_ooh > 0,
      compute_delta_salt(energykcal_eo_ooh, ellison_kcal_meal) * coverage, 0)]
    sp$pop[, black_race := as.integer(ethnicity %in% c("black african", "black caribbean"))]
    sp$pop[year >= sc_year,
           sbp_curr_xps := sbp_curr_xps + salt_sbp_effect(
             salt_difference = delta_salt_meal * ooh_freq_day,
             age = age, sbp = sbp_curr_xps, black_race = black_race,
             n = .N, stochastic = TRUE)]
    sp$pop[, c("black_race", "delta_salt_meal") := NULL]
  }
)
IMPACTncd$run(1:n_runs, multicore = TRUE, "sc4_energy_ns_all_sa_ellison")

# SA3-SC7: energy+NS, large (18%), FSA structural salt
set_scn(
  function(sp) {
    set.seed(as.integer(sp$mc) + 20250L)

    ellison_pct   <- truncnorm::rtruncnorm(1, b = 0, mean = ns_kcal_pct_reg$ellison$pct_mean, sd = ns_kcal_pct_reg$ellison_se)
    ellison_kcal_meal   <- sp$pop$energykcal_eo_ooh * ellison_pct
    coverage            <- 0.18
    delta_kcal_meal      <- ellison_kcal_meal * coverage
    sc_year <- 26L

    sp$pop[, `:=`(
      bmi_hall = pmax(13L, pmin(94L, as.integer(bmi_curr_xps))),
      calories = pmax(-200L, as.integer(delta_kcal_meal * ooh_freq_day))
    )]
    sp$pop[is.na(calories) | calories > 0L, calories := 0L]
    for (yr in sc_year:max(sp$pop$year)) {
      syn_yr <- copy(sp$pop[year == yr,
                            .(age, sex, bmi = bmi_hall, calories, pid, bmi_curr_xps)])
      syn_yr <- merge(syn_yr, change_bmi_Hall,
                      by = c("age", "sex", "bmi", "calories"), all.x = TRUE)
      setorder(syn_yr, pid)
      syn_yr[is.na(change_bmi_last) | calories == 0L, change_bmi_last := 0]
      syn_yr[, bmi_curr_xps := bmi_curr_xps + change_bmi_last]
      sp$pop[year == yr, bmi_curr_xps := syn_yr$bmi_curr_xps]
    }
    sp$pop[, c("bmi_hall", "calories") := NULL]

    sp$pop[, delta_salt_consumer := ifelse(energykcal_eo_ooh > 0,
      compute_delta_salt(energykcal_eo_ooh, ellison_kcal_meal) * coverage, 0)]
    sp$pop[, black_race := as.integer(ethnicity %in% c("black african", "black caribbean"))]
    sp$pop[year >= sc_year,
           sbp_curr_xps := sbp_curr_xps + salt_sbp_effect(
             salt_difference = (delta_salt_consumer * ooh_freq_day) +
               (-(salt_eo_ooh * ooh_freq_day) * 0.015 * pmin(year - sc_year + 1L, 20L)),
             age = age, sbp = sbp_curr_xps, black_race = black_race,
             n = .N, stochastic = TRUE)]
    sp$pop[, c("black_race", "delta_salt_consumer") := NULL]
  }
)
IMPACTncd$run(1:n_runs, multicore = TRUE, "sc7_energy_ns_fsa_large_sa_ellison")

# SA3-SC8: energy+NS, all (100%), FSA structural salt
set_scn(
  function(sp) {
    set.seed(as.integer(sp$mc) + 20250L)

    ellison_pct   <- truncnorm::rtruncnorm(1, b = 0, mean = ns_kcal_pct_reg$ellison$pct_mean, sd = ns_kcal_pct_reg$ellison_se)
    ellison_kcal_meal   <- sp$pop$energykcal_eo_ooh * ellison_pct
    coverage            <- 1.0
    delta_kcal_meal      <- ellison_kcal_meal * coverage
    sc_year <- 26L

    sp$pop[, `:=`(
      bmi_hall = pmax(13L, pmin(94L, as.integer(bmi_curr_xps))),
      calories = pmax(-200L, as.integer(delta_kcal_meal * ooh_freq_day))
    )]
    sp$pop[is.na(calories) | calories > 0L, calories := 0L]
    for (yr in sc_year:max(sp$pop$year)) {
      syn_yr <- copy(sp$pop[year == yr,
                            .(age, sex, bmi = bmi_hall, calories, pid, bmi_curr_xps)])
      syn_yr <- merge(syn_yr, change_bmi_Hall,
                      by = c("age", "sex", "bmi", "calories"), all.x = TRUE)
      setorder(syn_yr, pid)
      syn_yr[is.na(change_bmi_last) | calories == 0L, change_bmi_last := 0]
      syn_yr[, bmi_curr_xps := bmi_curr_xps + change_bmi_last]
      sp$pop[year == yr, bmi_curr_xps := syn_yr$bmi_curr_xps]
    }
    sp$pop[, c("bmi_hall", "calories") := NULL]

    sp$pop[, delta_salt_consumer := ifelse(energykcal_eo_ooh > 0,
      compute_delta_salt(energykcal_eo_ooh, ellison_kcal_meal) * coverage, 0)]
    sp$pop[, black_race := as.integer(ethnicity %in% c("black african", "black caribbean"))]
    sp$pop[year >= sc_year,
           sbp_curr_xps := sbp_curr_xps + salt_sbp_effect(
             salt_difference = (delta_salt_consumer * ooh_freq_day) +
               (-(salt_eo_ooh * ooh_freq_day) * 0.015 * pmin(year - sc_year + 1L, 20L)),
             age = age, sbp = sbp_curr_xps, black_race = black_race,
             n = .N, stochastic = TRUE)]
    sp$pop[, c("black_race", "delta_salt_consumer") := NULL]
  }
)
IMPACTncd$run(1:n_runs, multicore = TRUE, "sc8_energy_ns_fsa_all_sa_ellison")


# SA4: Nutri-Score coverage at 47% large-business turnover
# SA4-SC1: NS only, large (47% turnover), no FSA
set_scn(
  function(sp) {
    set.seed(as.integer(sp$mc) + 20250L)

    ns_pct   <- truncnorm::rtruncnorm(1, b = 0, mean = ns_kcal_pct_reg$ns$mean, sd = ns_kcal_pct_reg$ns$se)
    ns_kcal_meal   <- sp$pop$energykcal_eo_ooh * ns_pct
    coverage       <- 0.47
    delta_kcal_meal <- ns_kcal_meal * coverage
    sc_year <- 26L

    sp$pop[, `:=`(
      bmi_hall = pmax(13L, pmin(94L, as.integer(bmi_curr_xps))),
      calories = pmax(-200L, as.integer(delta_kcal_meal * ooh_freq_day))
    )]
    sp$pop[is.na(calories) | calories > 0L, calories := 0L]
    for (yr in sc_year:max(sp$pop$year)) {
      syn_yr <- copy(sp$pop[year == yr,
                            .(age, sex, bmi = bmi_hall, calories, pid, bmi_curr_xps)])
      syn_yr <- merge(syn_yr, change_bmi_Hall,
                      by = c("age", "sex", "bmi", "calories"), all.x = TRUE)
      setorder(syn_yr, pid)
      syn_yr[is.na(change_bmi_last) | calories == 0L, change_bmi_last := 0]
      syn_yr[, bmi_curr_xps := bmi_curr_xps + change_bmi_last]
      sp$pop[year == yr, bmi_curr_xps := syn_yr$bmi_curr_xps]
    }
    sp$pop[, c("bmi_hall", "calories") := NULL]

    sp$pop[, delta_salt_meal := ifelse(energykcal_eo_ooh > 0,
      compute_delta_salt(energykcal_eo_ooh, ns_kcal_meal) * coverage, 0)]
    sp$pop[, black_race := as.integer(ethnicity %in% c("black african", "black caribbean"))]
    sp$pop[year >= sc_year,
           sbp_curr_xps := sbp_curr_xps + salt_sbp_effect(
             salt_difference = delta_salt_meal * ooh_freq_day,
             age = age, sbp = sbp_curr_xps, black_race = black_race,
             n = .N, stochastic = TRUE)]
    sp$pop[, c("black_race", "delta_salt_meal") := NULL]
  }
)
IMPACTncd$run(1:n_runs, multicore = TRUE, "sc1_ns_large_sa_turnover47")

# SA4-SC3: energy+NS, large (47% turnover), no FSA
set_scn(
  function(sp) {
    set.seed(as.integer(sp$mc) + 20250L)

    ns_pct     <- truncnorm::rtruncnorm(1, b = 0, mean = ns_kcal_pct_reg$ns$mean, sd = ns_kcal_pct_reg$ns$se)
    ns_kcal_meal     <- sp$pop$energykcal_eo_ooh * ns_pct
    energy_pct <- truncnorm::rtruncnorm(1, b = 0, mean = ns_kcal_pct_reg$energy_label$pct_mean, sd = ns_kcal_pct_reg$energy_label_se)
    energy_kcal_meal <- sp$pop$energykcal_eo_ooh * energy_pct
    coverage         <- 0.47
    delta_kcal_meal   <- (ns_kcal_meal + energy_kcal_meal) * coverage
    sc_year <- 26L

    sp$pop[, `:=`(
      bmi_hall = pmax(13L, pmin(94L, as.integer(bmi_curr_xps))),
      calories = pmax(-200L, as.integer(delta_kcal_meal * ooh_freq_day))
    )]
    sp$pop[is.na(calories) | calories > 0L, calories := 0L]
    for (yr in sc_year:max(sp$pop$year)) {
      syn_yr <- copy(sp$pop[year == yr,
                            .(age, sex, bmi = bmi_hall, calories, pid, bmi_curr_xps)])
      syn_yr <- merge(syn_yr, change_bmi_Hall,
                      by = c("age", "sex", "bmi", "calories"), all.x = TRUE)
      setorder(syn_yr, pid)
      syn_yr[is.na(change_bmi_last) | calories == 0L, change_bmi_last := 0]
      syn_yr[, bmi_curr_xps := bmi_curr_xps + change_bmi_last]
      sp$pop[year == yr, bmi_curr_xps := syn_yr$bmi_curr_xps]
    }
    sp$pop[, c("bmi_hall", "calories") := NULL]

    sp$pop[, delta_salt_meal := ifelse(energykcal_eo_ooh > 0,
      compute_delta_salt(energykcal_eo_ooh, (ns_kcal_meal + energy_kcal_meal)) * coverage, 0)]
    sp$pop[, black_race := as.integer(ethnicity %in% c("black african", "black caribbean"))]
    sp$pop[year >= sc_year,
           sbp_curr_xps := sbp_curr_xps + salt_sbp_effect(
             salt_difference = delta_salt_meal * ooh_freq_day,
             age = age, sbp = sbp_curr_xps, black_race = black_race,
             n = .N, stochastic = TRUE)]
    sp$pop[, c("black_race", "delta_salt_meal") := NULL]
  }
)
IMPACTncd$run(1:n_runs, multicore = TRUE, "sc3_energy_ns_large_sa_turnover47")

# SA4-SC5: NS only, large (47% turnover), FSA structural salt
set_scn(
  function(sp) {
    set.seed(as.integer(sp$mc) + 20250L)

    ns_pct        <- truncnorm::rtruncnorm(1, b = 0, mean = ns_kcal_pct_reg$ns$mean, sd = ns_kcal_pct_reg$ns$se)
    ns_kcal_meal        <- sp$pop$energykcal_eo_ooh * ns_pct
    coverage            <- 0.47
    delta_kcal_meal      <- ns_kcal_meal * coverage
    sc_year <- 26L

    sp$pop[, `:=`(
      bmi_hall = pmax(13L, pmin(94L, as.integer(bmi_curr_xps))),
      calories = pmax(-200L, as.integer(delta_kcal_meal * ooh_freq_day))
    )]
    sp$pop[is.na(calories) | calories > 0L, calories := 0L]
    for (yr in sc_year:max(sp$pop$year)) {
      syn_yr <- copy(sp$pop[year == yr,
                            .(age, sex, bmi = bmi_hall, calories, pid, bmi_curr_xps)])
      syn_yr <- merge(syn_yr, change_bmi_Hall,
                      by = c("age", "sex", "bmi", "calories"), all.x = TRUE)
      setorder(syn_yr, pid)
      syn_yr[is.na(change_bmi_last) | calories == 0L, change_bmi_last := 0]
      syn_yr[, bmi_curr_xps := bmi_curr_xps + change_bmi_last]
      sp$pop[year == yr, bmi_curr_xps := syn_yr$bmi_curr_xps]
    }
    sp$pop[, c("bmi_hall", "calories") := NULL]

    sp$pop[, delta_salt_consumer := ifelse(energykcal_eo_ooh > 0,
      compute_delta_salt(energykcal_eo_ooh, ns_kcal_meal) * coverage, 0)]
    sp$pop[, black_race := as.integer(ethnicity %in% c("black african", "black caribbean"))]
    sp$pop[year >= sc_year,
           sbp_curr_xps := sbp_curr_xps + salt_sbp_effect(
             salt_difference = (delta_salt_consumer * ooh_freq_day) +
               (-(salt_eo_ooh * ooh_freq_day) * 0.015 * pmin(year - sc_year + 1L, 20L)),
             age = age, sbp = sbp_curr_xps, black_race = black_race,
             n = .N, stochastic = TRUE)]
    sp$pop[, c("black_race", "delta_salt_consumer") := NULL]
  }
)
IMPACTncd$run(1:n_runs, multicore = TRUE, "sc5_ns_fsa_large_sa_turnover47")

# SA4-SC7: energy+NS, large (47% turnover), FSA structural salt
set_scn(
  function(sp) {
    set.seed(as.integer(sp$mc) + 20250L)

    ns_pct        <- truncnorm::rtruncnorm(1, b = 0, mean = ns_kcal_pct_reg$ns$mean, sd = ns_kcal_pct_reg$ns$se)
    ns_kcal_meal        <- sp$pop$energykcal_eo_ooh * ns_pct
    energy_pct    <- truncnorm::rtruncnorm(1, b = 0, mean = ns_kcal_pct_reg$energy_label$pct_mean, sd = ns_kcal_pct_reg$energy_label_se)
    energy_kcal_meal    <- sp$pop$energykcal_eo_ooh * energy_pct
    coverage            <- 0.47
    delta_kcal_meal      <- (ns_kcal_meal + energy_kcal_meal) * coverage
    sc_year <- 26L

    sp$pop[, `:=`(
      bmi_hall = pmax(13L, pmin(94L, as.integer(bmi_curr_xps))),
      calories = pmax(-200L, as.integer(delta_kcal_meal * ooh_freq_day))
    )]
    sp$pop[is.na(calories) | calories > 0L, calories := 0L]
    for (yr in sc_year:max(sp$pop$year)) {
      syn_yr <- copy(sp$pop[year == yr,
                            .(age, sex, bmi = bmi_hall, calories, pid, bmi_curr_xps)])
      syn_yr <- merge(syn_yr, change_bmi_Hall,
                      by = c("age", "sex", "bmi", "calories"), all.x = TRUE)
      setorder(syn_yr, pid)
      syn_yr[is.na(change_bmi_last) | calories == 0L, change_bmi_last := 0]
      syn_yr[, bmi_curr_xps := bmi_curr_xps + change_bmi_last]
      sp$pop[year == yr, bmi_curr_xps := syn_yr$bmi_curr_xps]
    }
    sp$pop[, c("bmi_hall", "calories") := NULL]

    sp$pop[, delta_salt_consumer := ifelse(energykcal_eo_ooh > 0,
      compute_delta_salt(energykcal_eo_ooh, (ns_kcal_meal + energy_kcal_meal)) * coverage, 0)]
    sp$pop[, black_race := as.integer(ethnicity %in% c("black african", "black caribbean"))]
    sp$pop[year >= sc_year,
           sbp_curr_xps := sbp_curr_xps + salt_sbp_effect(
             salt_difference = (delta_salt_consumer * ooh_freq_day) +
               (-(salt_eo_ooh * ooh_freq_day) * 0.015 * pmin(year - sc_year + 1L, 20L)),
             age = age, sbp = sbp_curr_xps, black_race = black_race,
             n = .N, stochastic = TRUE)]
    sp$pop[, c("black_race", "delta_salt_consumer") := NULL]
  }
)
IMPACTncd$run(1:n_runs, multicore = TRUE, "sc7_energy_ns_fsa_large_sa_turnover47")


# =============================================================================
# EXPORT
# =============================================================================
IMPACTncd$export_summaries(multicore = TRUE)
IMPACTncd$export_tables(multicore = TRUE)
message("All scenarios complete.")
