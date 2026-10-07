# Nutri-Score labelling and salt reformulation in England’s out-of-home food sector

This repository contains the code and data used to estimate the health and economic effects of introducing Nutri-Score labelling and a salt reformulation programme in England’s out-of-home (OOH) food sector.

## Model framework and code development

This study uses the [IMPACTncd microsimulation framework for England](https://github.com/ChristK/IMPACTncd_Engl), developed by Chris Kypridemos and contributors.

All study-specific code provided in this repository was developed by **Flavia Alves**.

## Repository contents

| File | Description |
| --- | --- |
| [`code/sim_design_nutriscore_salt.yaml`](code/sim_design_nutriscore_salt.yaml) | Simulation design and configuration for the Nutri-Score and salt reformulation scenarios. |
| [`code/simulate_nutriscore_salt.R`](code/simulate_nutriscore_salt.R) | R script for running the policy simulations. |
| [`data_models/ns_kcal_pct_reg.rds`](data_models/ns_kcal_pct_reg.rds) | Saved regression model used to estimate Nutri-Score-related percentage changes in energy intake. |
| [`data_models/salt_energy_model_final.rds`](data_models/salt_energy_model_final.rds) | Saved model describing the relationship between salt and energy intake. |
| [`data_tables/energykcal_eo_ooh_table.csv`](data_tables/energykcal_eo_ooh_table.csv) | OOH energy intake per eating occasion. |
| [`data_tables/ooh_freq_day_table.csv`](data_tables/ooh_freq_day_table.csv) | Daily frequency of OOH eating occasions. |
| [`data_tables/salt_eo_ooh_table.csv`](data_tables/salt_eo_ooh_table.csv) | OOH salt intake per eating occasion. |
| [`LICENSE`](LICENSE) | GNU General Public License v3. |

## Using the code

The study-specific files in this repository are intended for use with the IMPACTncd framework for England. Please refer to the [IMPACTncd repository](https://github.com/ChristK/IMPACTncd_Engl) for model installation instructions, dependencies and documentation.

The study configuration is provided in `code/sim_design_nutriscore_salt.yaml`, and the simulation script is `code/simulate_nutriscore_salt.R`. Supporting models and input tables are provided in `data_models/` and `data_tables/`, respectively.

## Associated manuscript

**Estimated health and economic effects of introducing Nutri-Score labelling and a salt reformulation programme in England’s out-of-home food sector: a microsimulation study**

Flavia Alves¹, Martin O’Flaherty¹, Amy Finlay², Zoi Toumpakari³, Nick Townsend³, Andrew Jones⁴, James Garbutt⁵, Chris Kypridemos¹, Eric Robinson², Zoé Colombet¹

1. Department of Public Health, Policy and Systems, University of Liverpool, Liverpool, UK
2. Department of Psychology, University of Liverpool, Liverpool, UK
3. School for Policy Studies, University of Bristol, Bristol, UK
4. School of Psychology, Liverpool John Moores University, Liverpool, UK
5. Department for Health, University of Bath, Bath, UK

## Citation

If you use these study-specific materials, please cite the associated manuscript and acknowledge the IMPACTncd framework and its developers. Full publication details and a DOI will be added when available.

## Licence

This repository is licensed under the GNU General Public License v3. See the [LICENSE](LICENSE) file for details.

## Contact

For questions about the study-specific code or analysis, please contact [Flavia Alves](mailto:Flavia.Alves@liverpool.ac.uk).
