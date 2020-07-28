## ---------------------------------------------------- ##
## organize_results.R --------------------------------- ##
## Purpose: organize simulation results --------------- ##
## Author: Peter Norwood, NCSU, Janssen R&D ----------- ##
## ---------------------------------------------------- ##

## set wd
setwd("<LIBRARY>")

## load packages
library(dplyr)
library(tidyr)
library(Matrix)

## covariance matricies
G <- rbind(c(0.959,0.049,0.730,-0.043),
           c(0.049,0.031,0.033,0.010),
           c(0.730,0.033,1.826,-0.025),
           c(-0.043,0.010,-0.025,0.008))

R_diag <- c(sqrt(1.55),sqrt(0.571))

Z <- rbind(c(1,0,0,0),
           c(0,0,1,0),
           c(1,2,0,0),
           c(0,0,1,2),
           c(1,3,0,0),
           c(0,0,1,3),
           c(1,4,0,0),
           c(0,0,1,4))

## true mean parameters
t_beta10=7.044        
t_beta20=1.243
t_beta11=0.068
t_beta21=0.081
t_beta12=-0.056
t_beta22=-0.004
t_beta13=-0.107
t_beta23=-0.081

t_beta <- c(t_beta10,t_beta20,t_beta11,t_beta21,
            t_beta12,t_beta22,t_beta13,t_beta23)



## ---------------------------------------------------------------------------------- ##
## ---------------------------------------------------------------------------------- ##


## extract the non-broken simulations for rho_0
load("rho_0.RData")
good_set <- list()
good_set[[1]] <- list()
good_set[[2]] <- list()
good_set[[3]] <- list()
good_set[[4]] <- list()
good_set[[5]] <- list()

## generate the good stes to loop through
for(i in 1:5){
  for(j in 1:1000){
    if(is.list(rho_0_sim[[i]][[j]])){
      good_set[[i]][j] <- j
    }else{
      ## do nothin
    }
  }
}

## extract the non-broken rho_0 information
rho_0_dat <- data.frame()
rho_0_pvals <- data.frame()
for(i in 1:5){
  tick=1
  for(j in unlist(good_set[[i]])){
    
    lst <- rho_0_sim[[i]][[j]]
    
    ## generate R matrix
    rho=lst$rho
    R <- diag(R_diag) %*% rbind(c(1,rho),c(rho,1)) %*% diag(R_diag)
    
    full_R <- bdiag(R,R,R,R)
    
    V <- Z %*% G %*% t(Z) + full_R
    
    V_norm_gls <- norm(V-lst$gls_V,"F")
    
    beta_norm_lme <- sum((t_beta-lst$lme_beta)**2)
    beta_norm_lme2 <- sum((t_beta-lst$lme_big_beta)**2)
    beta_norm_gls <- sum((t_beta-lst$gls_beta)**2)
    
    
    temp <- data.frame(
              N=rep(lst$N,3),
              rho=rep(lst$rho,3),
              rep=rep(tick,3),
              model=c("mixed_diag_R","mixed_unstr_R","gls"),
              V_norm=c(lst$model_info$V_norm[1],lst$model_info$V_norm[2],V_norm_gls),
              G_norm=c(lst$model_info$G_norm),
              beta_norm=c(beta_norm_lme,beta_norm_lme2,beta_norm_gls),
              AIC=c(lst$model_info$AIC),
              BIC=c(lst$model_info$BIC),
              var_y11=c(lst$lme_V[1,1],lst$lme_big_V[1,1],lst$gls_V[1,1]),
              cov_y11_y22=c(lst$lme_V[1,4],lst$lme_big_V[1,4],lst$gls_V[1,4]),
              cov_y11_y21=c(lst$lme_V[1,2],lst$lme_big_V[1,2],lst$gls_V[1,2]),
              var_e1=c(lst$lme_R[1,1],lst$lme_big_R[1,1],NA),
              var_e2=c(lst$lme_R[2,2],lst$lme_big_R[2,2],NA),
              cor_e11_e21=c(NA,lst$lme_big_R[1,2]/(sqrt(lst$lme_big_R[1,1])*sqrt(lst$lme_big_R[2,2])),NA),
              beta10=c(lst$lme_beta[1],lst$lme_big_beta[1],lst$gls_beta[1]),
              beta20=c(lst$lme_beta[2],lst$lme_big_beta[2],lst$gls_beta[2]),
              beta11=c(lst$lme_beta[3],lst$lme_big_beta[3],lst$gls_beta[3]),
              beta21=c(lst$lme_beta[4],lst$lme_big_beta[4],lst$gls_beta[4]),
              beta12=c(lst$lme_beta[5],lst$lme_big_beta[5],lst$gls_beta[5]),
              beta22=c(lst$lme_beta[6],lst$lme_big_beta[6],lst$gls_beta[6]),
              beta13=c(lst$lme_beta[7],lst$lme_big_beta[7],lst$gls_beta[7]),
              beta23=c(lst$lme_beta[8],lst$lme_big_beta[8],lst$gls_beta[8])
    )
    
    temp_pvals <- lst$p_vals
    
    rho_0_dat <- rbind(rho_0_dat,temp)
    rho_0_pvals <- rbind(rho_0_pvals,temp_pvals)
    tick=tick+1

    
  }
}

## save the rho_0 information
save(rho_0_dat,file="rho_0_dat.RData")


## ---------------------------------------------------------------------------------- ##
## ---------------------------------------------------------------------------------- ##

## same process for rho_25

load("rho_25.RData")
good_set <- list()
good_set[[1]] <- list()
good_set[[2]] <- list()
good_set[[3]] <- list()
good_set[[4]] <- list()
good_set[[5]] <- list()


## since some of the rho_25 data for certain sample sizes
## never failed to converge, we have different list structure
## we can figure this out by just looking at the list with View(rho_25)
## and seeing which have certain structures
## so we split it up into part1 and part2

## generate the good stes to loop through
for(i in c(1,5)){
  for(j in 1:1000){
    if(is.list(rho_25_sim[[i]][[j]])){
      good_set[[i]][j] <- j
    }else{
      ## do nothin
    }
  }
}

rho_25_dat <- data.frame()
rho_25_pvals <- data.frame()
for(i in c(1,5)){
  tick=1
  for(j in unlist(good_set[[i]])){
    
    lst <- rho_25_sim[[i]][[j]]
    
    ## generate R matrix
    rho=lst$rho
    R <- diag(R_diag) %*% rbind(c(1,rho),c(rho,1)) %*% diag(R_diag)
    
    full_R <- bdiag(R,R,R,R)
    
    V <- Z %*% G %*% t(Z) + full_R
    
    V_norm_gls <- norm(V-lst$gls_V,"F")
    
    beta_norm_lme <- sum((t_beta-lst$lme_beta)**2)
    beta_norm_lme2 <- sum((t_beta-lst$lme_big_beta)**2)
    beta_norm_gls <- sum((t_beta-lst$gls_beta)**2)
    
    
    temp <- data.frame(
      N=rep(lst$N,3),
      rho=rep(lst$rho,3),
      rep=rep(tick,3),
      model=c("mixed_diag_R","mixed_unstr_R","gls"),
      V_norm=c(lst$model_info$V_norm[1],lst$model_info$V_norm[2],V_norm_gls),
      G_norm=c(lst$model_info$G_norm),
      beta_norm=c(beta_norm_lme,beta_norm_lme2,beta_norm_gls),
      AIC=c(lst$model_info$AIC),
      BIC=c(lst$model_info$BIC),
      var_y11=c(lst$lme_V[1,1],lst$lme_big_V[1,1],lst$gls_V[1,1]),
      cov_y11_y22=c(lst$lme_V[1,4],lst$lme_big_V[1,4],lst$gls_V[1,4]),
      cov_y11_y21=c(lst$lme_V[1,2],lst$lme_big_V[1,2],lst$gls_V[1,2]),
      var_e1=c(lst$lme_R[1,1],lst$lme_big_R[1,1],NA),
      var_e2=c(lst$lme_R[2,2],lst$lme_big_R[2,2],NA),
      beta10=c(lst$lme_beta[1],lst$lme_big_beta[1],lst$gls_beta[1]),
      cor_e11_e21=c(NA,lst$lme_big_R[1,2]/(sqrt(lst$lme_big_R[1,1])*sqrt(lst$lme_big_R[2,2])),NA),
      beta20=c(lst$lme_beta[2],lst$lme_big_beta[2],lst$gls_beta[2]),
      beta11=c(lst$lme_beta[3],lst$lme_big_beta[3],lst$gls_beta[3]),
      beta21=c(lst$lme_beta[4],lst$lme_big_beta[4],lst$gls_beta[4]),
      beta12=c(lst$lme_beta[5],lst$lme_big_beta[5],lst$gls_beta[5]),
      beta22=c(lst$lme_beta[6],lst$lme_big_beta[6],lst$gls_beta[6]),
      beta13=c(lst$lme_beta[7],lst$lme_big_beta[7],lst$gls_beta[7]),
      beta23=c(lst$lme_beta[8],lst$lme_big_beta[8],lst$gls_beta[8])
    )
    
    temp_pvals <- lst$p_vals
    
    rho_25_dat <- rbind(rho_25_dat,temp)
    rho_25_pvals <- rbind(rho_25_pvals,temp_pvals)
    tick=tick+1
    
    
  }
}

rho_25_dat_part1 <- rho_25_dat
rho_25_pvals_part1 <- rho_25_pvals

## PART TWO
rho_25_dat <- data.frame()
rho_25_pvals <- data.frame()
for(i in c(2,3,4)){
  tick=1
  for(j in 1:1000){
    
    lst <- rho_25_sim[[i]][,j]
    
    ## generate R matrix
    rho=lst$rho
    R <- diag(R_diag) %*% rbind(c(1,rho),c(rho,1)) %*% diag(R_diag)
    
    full_R <- bdiag(R,R,R,R)
    
    V <- Z %*% G %*% t(Z) + full_R
    
    V_norm_gls <- norm(V-lst$gls_V,"F")
    
    beta_norm_lme <- sum((t_beta-lst$lme_beta)**2)
    beta_norm_lme2 <- sum((t_beta-lst$lme_big_beta)**2)
    beta_norm_gls <- sum((t_beta-lst$gls_beta)**2)
    
    
    temp <- data.frame(
      N=rep(lst$N,3),
      rho=rep(lst$rho,3),
      rep=rep(tick,3),
      model=c("mixed_diag_R","mixed_unstr_R","gls"),
      V_norm=c(lst$model_info$V_norm[1],lst$model_info$V_norm[2],V_norm_gls),
      G_norm=c(lst$model_info$G_norm),
      beta_norm=c(beta_norm_lme,beta_norm_lme2,beta_norm_gls),
      AIC=c(lst$model_info$AIC),
      BIC=c(lst$model_info$BIC),
      var_y11=c(lst$lme_V[1,1],lst$lme_big_V[1,1],lst$gls_V[1,1]),
      cov_y11_y22=c(lst$lme_V[1,4],lst$lme_big_V[1,4],lst$gls_V[1,4]),
      cov_y11_y21=c(lst$lme_V[1,2],lst$lme_big_V[1,2],lst$gls_V[1,2]),
      var_e1=c(lst$lme_R[1,1],lst$lme_big_R[1,1],NA),
      var_e2=c(lst$lme_R[2,2],lst$lme_big_R[2,2],NA),
      cor_e11_e21=c(NA,lst$lme_big_R[1,2]/(sqrt(lst$lme_big_R[1,1])*sqrt(lst$lme_big_R[2,2])),NA),
      beta10=c(lst$lme_beta[1],lst$lme_big_beta[1],lst$gls_beta[1]),
      beta20=c(lst$lme_beta[2],lst$lme_big_beta[2],lst$gls_beta[2]),
      beta11=c(lst$lme_beta[3],lst$lme_big_beta[3],lst$gls_beta[3]),
      beta21=c(lst$lme_beta[4],lst$lme_big_beta[4],lst$gls_beta[4]),
      beta12=c(lst$lme_beta[5],lst$lme_big_beta[5],lst$gls_beta[5]),
      beta22=c(lst$lme_beta[6],lst$lme_big_beta[6],lst$gls_beta[6]),
      beta13=c(lst$lme_beta[7],lst$lme_big_beta[7],lst$gls_beta[7]),
      beta23=c(lst$lme_beta[8],lst$lme_big_beta[8],lst$gls_beta[8])
    )
    
    temp_pvals <- lst$p_vals
    
    rho_25_dat <- rbind(rho_25_dat,temp)
    rho_25_pvals <- rbind(rho_25_pvals,temp_pvals)
    tick=tick+1
    
    
  }
}

rho_25_dat_part2 <- rho_25_dat
rho_25_pvals_part2 <- rho_25_pvals

rho_25_dat <- rbind(rho_25_dat_part1,rho_25_dat_part2)
rho_25_pvals <- rbind(rho_25_pvals_part1,rho_25_pvals_part2)

save(rho_25_dat,file="rho_25_dat.RData")
save(rho_25_pvals,file="rho_25_pvals.RData")


## ---------------------------------------------------------------------------------- ##
## ---------------------------------------------------------------------------------- ##

## same process for rho_50

load("rho_50.RData")
good_set <- list()
good_set[[1]] <- list()
good_set[[2]] <- list()
good_set[[3]] <- list()
good_set[[4]] <- list()
good_set[[5]] <- list()

## generate the good stes to loop through
for(i in 1:5){
  for(j in 1:1000){
    if(is.list(rho_50_sim[[i]][[j]])){
      good_set[[i]][j] <- j
    }else{
      ## do nothin
    }
  }
}

rho_50_dat <- data.frame()
rho_50_pvals <- data.frame()
for(i in 1:5){
  tick=1
  for(j in unlist(good_set[[i]])){
    
    lst <- rho_50_sim[[i]][[j]]
    
    ## generate R matrix
    rho=lst$rho
    R <- diag(R_diag) %*% rbind(c(1,rho),c(rho,1)) %*% diag(R_diag)
    
    full_R <- bdiag(R,R,R,R)
    
    V <- Z %*% G %*% t(Z) + full_R
    
    V_norm_gls <- norm(V-lst$gls_V,"F")
    
    beta_norm_lme <- sum((t_beta-lst$lme_beta)**2)
    beta_norm_lme2 <- sum((t_beta-lst$lme_big_beta)**2)
    beta_norm_gls <- sum((t_beta-lst$gls_beta)**2)
    
    
    temp <- data.frame(
      N=rep(lst$N,3),
      rho=rep(lst$rho,3),
      rep=rep(tick,3),
      model=c("mixed_diag_R","mixed_unstr_R","gls"),
      V_norm=c(lst$model_info$V_norm[1],lst$model_info$V_norm[2],V_norm_gls),
      G_norm=c(lst$model_info$G_norm),
      beta_norm=c(beta_norm_lme,beta_norm_lme2,beta_norm_gls),
      AIC=c(lst$model_info$AIC),
      BIC=c(lst$model_info$BIC),
      var_y11=c(lst$lme_V[1,1],lst$lme_big_V[1,1],lst$gls_V[1,1]),
      cov_y11_y22=c(lst$lme_V[1,4],lst$lme_big_V[1,4],lst$gls_V[1,4]),
      cov_y11_y21=c(lst$lme_V[1,2],lst$lme_big_V[1,2],lst$gls_V[1,2]),
      var_e1=c(lst$lme_R[1,1],lst$lme_big_R[1,1],NA),
      var_e2=c(lst$lme_R[2,2],lst$lme_big_R[2,2],NA),
      cor_e11_e21=c(NA,lst$lme_big_R[1,2]/(sqrt(lst$lme_big_R[1,1])*sqrt(lst$lme_big_R[2,2])),NA),
      beta10=c(lst$lme_beta[1],lst$lme_big_beta[1],lst$gls_beta[1]),
      beta20=c(lst$lme_beta[2],lst$lme_big_beta[2],lst$gls_beta[2]),
      beta11=c(lst$lme_beta[3],lst$lme_big_beta[3],lst$gls_beta[3]),
      beta21=c(lst$lme_beta[4],lst$lme_big_beta[4],lst$gls_beta[4]),
      beta12=c(lst$lme_beta[5],lst$lme_big_beta[5],lst$gls_beta[5]),
      beta22=c(lst$lme_beta[6],lst$lme_big_beta[6],lst$gls_beta[6]),
      beta13=c(lst$lme_beta[7],lst$lme_big_beta[7],lst$gls_beta[7]),
      beta23=c(lst$lme_beta[8],lst$lme_big_beta[8],lst$gls_beta[8])
    )
    
    temp_pvals <- lst$p_vals
    
    rho_50_dat <- rbind(rho_50_dat,temp)
    rho_50_pvals <- rbind(rho_50_pvals,temp_pvals)
    tick=tick+1
    
    
  }
}


save(rho_50_dat,file="rho_50_dat.RData")
save(rho_50_pvals,file="rho_50_pvals.RData")


## ---------------------------------------------------------------------------------- ##
## ---------------------------------------------------------------------------------- ##

## same process for rho_75
load("7_7_rho_75.RData")
good_set <- list()
good_set[[1]] <- list()
good_set[[2]] <- list()
good_set[[3]] <- list()
good_set[[4]] <- list()
good_set[[5]] <- list()

## generate the good stes to loop through
for(i in 1:5){
  for(j in 1:1000){
    if(is.list(rho_75_sim[[i]][[j]])){
      good_set[[i]][j] <- j
    }else{
      ## do nothin
    }
  }
}

rho_75_dat <- data.frame()
rho_75_pvals <- data.frame()
for(i in 1:5){
  tick=1
  for(j in unlist(good_set[[i]])){
    
    lst <- rho_75_sim[[i]][[j]]
    
    ## generate R matrix
    rho=lst$rho
    R <- diag(R_diag) %*% rbind(c(1,rho),c(rho,1)) %*% diag(R_diag)
    
    full_R <- bdiag(R,R,R,R)
    
    V <- Z %*% G %*% t(Z) + full_R
    
    V_norm_gls <- norm(V-lst$gls_V,"F")
    
    beta_norm_lme <- sum((t_beta-lst$lme_beta)**2)
    beta_norm_lme2 <- sum((t_beta-lst$lme_big_beta)**2)
    beta_norm_gls <- sum((t_beta-lst$gls_beta)**2)
    
    
    temp <- data.frame(
      N=rep(lst$N,3),
      rho=rep(lst$rho,3),
      rep=rep(tick,3),
      model=c("mixed_diag_R","mixed_unstr_R","gls"),
      V_norm=c(lst$model_info$V_norm[1],lst$model_info$V_norm[2],V_norm_gls),
      G_norm=c(lst$model_info$G_norm),
      beta_norm=c(beta_norm_lme,beta_norm_lme2,beta_norm_gls),
      AIC=c(lst$model_info$AIC),
      BIC=c(lst$model_info$BIC),
      var_y11=c(lst$lme_V[1,1],lst$lme_big_V[1,1],lst$gls_V[1,1]),
      cov_y11_y22=c(lst$lme_V[1,4],lst$lme_big_V[1,4],lst$gls_V[1,4]),
      cov_y11_y21=c(lst$lme_V[1,2],lst$lme_big_V[1,2],lst$gls_V[1,2]),
      var_e1=c(lst$lme_R[1,1],lst$lme_big_R[1,1],NA),
      var_e2=c(lst$lme_R[2,2],lst$lme_big_R[2,2],NA),
      cor_e11_e21=c(NA,lst$lme_big_R[1,2]/(sqrt(lst$lme_big_R[1,1])*sqrt(lst$lme_big_R[2,2])),NA),
      beta10=c(lst$lme_beta[1],lst$lme_big_beta[1],lst$gls_beta[1]),
      beta20=c(lst$lme_beta[2],lst$lme_big_beta[2],lst$gls_beta[2]),
      beta11=c(lst$lme_beta[3],lst$lme_big_beta[3],lst$gls_beta[3]),
      beta21=c(lst$lme_beta[4],lst$lme_big_beta[4],lst$gls_beta[4]),
      beta12=c(lst$lme_beta[5],lst$lme_big_beta[5],lst$gls_beta[5]),
      beta22=c(lst$lme_beta[6],lst$lme_big_beta[6],lst$gls_beta[6]),
      beta13=c(lst$lme_beta[7],lst$lme_big_beta[7],lst$gls_beta[7]),
      beta23=c(lst$lme_beta[8],lst$lme_big_beta[8],lst$gls_beta[8])
    )
    
    temp_pvals <- lst$p_vals
    
    rho_75_dat <- rbind(rho_75_dat,temp)
    rho_75_pvals <- rbind(rho_75_pvals,temp_pvals)
    tick=tick+1
    
    
  }
}


save(rho_75_dat,file="rho_75_dat.RData")
save(rho_75_pvals,file="rho_75_pvals.RData")


## ---------------------------------------------------------------------------------- ##
## ---------------------------------------------------------------------------------- ##

## save overall simulation results

simulation_dat <- rbind(rho_0_dat,rho_25_dat,rho_50_dat,rho_75_dat)
simulation_pvals <- rbind(rho_0_pvals,rho_25_pvals,rho_50_pvals,rho_75_pvals)

save(simulation_dat,file="simulation_dat.RData")
save(simulation_pvals,file="simulation_pvals.RData")
