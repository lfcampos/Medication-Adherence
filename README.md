# Introduction

This repository contains code needed to analyze data as in the following two papers:

'Measuring Effects of Medication Adherence on Time-Varying Health Outcome'
https://academic.oup.com/biostatistics/advance-article-abstract/doi/10.1093/biostatistics/kxz059/5687018

'Inferring medication adherence from time-varying health outcomes'
URL TBD.

# Inference on model parameters

To run the analysis in 'Measuring Effects of Medication Adherence on Time-Varying Health Outcome', the main file is run_analysis.R. Sourcing this file will generate a simulated dataset and produce posterior draws of the state space model and save them out in a folder called 'data'.

Before running this file, change the code to point to the proper working directory.

The other files for inference on the model parameters are:

- stan_functions.R
- state_space_functions.R
- state_space_adherence.stan

# Predicting adherence

To run the analysis in 'Inferring medication adherence from time-varying health outcomes', the files are in the 'predict' folder.

The main file is predict.R. Sourcing this file will read in the data and posterior draws of the model parameters from the 'data' folder, produce posterior draws of adherence values, and save them in a folder with the timestamp of the runtime in the 'out' folder.

Before running this file, change the code to point to the proper working directory.

The file run_params.R allows the user to easily adjust the tuning parameters, such as the number of particles or the number of test patients.

After running predict.R, the user should run summarize.R. At the top of the file, the user provides the specific analysis run date of interest and summarizes the posterior draws in plots.

The other files for predicting adherence behavior are:

- setup.R
- functions.R
- functions_pf.R
- plots.R
- ran_ef.stan

# References

Particle filter code in functions_pf.R is from the [CoupledCPF](https://github.com/pierrejacob/CoupledCPF) repo by pierrejacob.


