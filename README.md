# Elastic Net Modeling of AMR data.

This repository is to document all code used for the analysis of antimicrobial resistance in Escherichia coli isolates taken from veterinary animals. 
Analysis comprises descriptive statistcs, figure generation, and elastic net modeling of isolate genotype and phenotype data. Analysis is intended to run in sequential order from 1-4. 

Data for analysis is available on figshare (link)

Shield: [![CC BY 4.0][cc-by-shield]][cc-by]

This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg

## Scripts

- 1_descriptive: generates figures and tables calculating descriptive statistics for sample genotype and phenotype data.

- 2_elasticnet: fits six different elastic nets to sample data to predict isolate resistance. 

- 3_power: estimates power for candidate predictors used in best fit elastic net model.

- 3_animal : calculates prevalence of antibiotic resistance in animals and determines significant differences in resistance proportion. 

- 4_outputs: calculates statistics and generates tables combining the important predictors and estimated power from previous scripts.

