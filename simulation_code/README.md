# simulation_code
Code to Run Simulation Study on Multivariate Longitduinal Models

Below are descriptions of each document included.

* gen_data.R: R script to generate data under a random slopes, random intercept mixed model with normally distributed variance components
* sim_function.R: R script including functions to simulate a scenario where data is generated, models are fit, and model information is recorded
* run_sims.R: R script that runs the simulations under different scenarios.
* organize_results.R: R script to throw out simulations that crashed due to computational issues and organize the results into a useful dataframe.
* analyze_results.R: R script that makes figures and tables of the simulation results.
