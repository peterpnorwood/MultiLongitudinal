## ---------------------------------------------------- ##
## run_sims.R ----------------------------------------- ##
## Purpose: start simulation for a certain rhos ------- ##
## Author: Peter Norwood, NCSU, Janssen R&D ----------- ##
## ---------------------------------------------------- ##

setwd("~/Summer Project/Multivariate Simulation/Precision")
source("additional_analysis.R")

## vector of sample sizes and correlations
Ns <- c(40,100,200,300,400)
rhos <- c(0,0.25,0.5,0.75)
## number of monte carlo reps
b <- 1e3


## Note: in practice I made four different files, one for
## each rho simulation and ran them simultaneously
## I also could have used mclapply or parLapply instead of 
## lapply to add additional parallelization (the mc reps inside
## of rep_sim are paralleized already


rho_0_sim <- lapply(Ns,rep_sim,B=b,r=rhos[1])
save(rho_0_sim,file="rho_0.RData")

rho_25_sim <- lapply(Ns,rep_sim,B=b,r=rhos[2])
save(rho_25_sim,file="rho_25.RData")

rho_50_sim <- lapply(Ns,rep_sim,B=b,r=rhos[3])
save(rho_50_sim,file="rho_50.RData")

rho_75_sim <- lapply(Ns,rep_sim,B=b,r=rhos[4])
save(rho_75_sim,file="rho_75.RData")
