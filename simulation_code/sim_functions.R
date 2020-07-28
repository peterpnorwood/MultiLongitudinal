## ---------------------------------------------------- ##
## sim_function.R ------------------------------------- ##
## Purpose: simulate data under a random slope, ------- ##
## random intercept model with normal errors ---------- ##
## Author: Peter Norwood, NCSU, Janssen --------------- ##
## ---------------------------------------------------- ##


setwd("<LIBRARY>")
source("gen_data.R")

library(dplyr)
library(nlme)

## covariance matricies
G <- rbind(c(0.959,0.049,0.730,-0.043),
           c(0.049,0.031,0.033,0.010),
           c(0.730,0.033,1.826,-0.025),
           c(-0.043,0.010,-0.025,0.008))

R_diag <- c(sqrt(1.55),sqrt(0.571))

## mean parameters
beta1 <- c(7.044,0.068,-0.056,-0.107)
beta2 <- c(1.234,0.081,0.004,-0.081)

## design for one individual with treatment 0
X0 <- rbind(c(1,0,0,0),
            c(1,0,2,0),
            c(1,0,4,0),
            c(1,0,8,0))

## design for one individual with treatment 1
X1 <- rbind(c(1,1,0,0),
            c(1,1,2,2),
            c(1,1,4,4),
            c(1,1,8,8))

cntrl <- lmeControl(maxIter=1000,
                    msMaxIter=1000,
                    niterEM=1000,
                    opt="optim")

## sim
## Purpose: generate random slopes, random intercept data,
##          fit models, and extract model information
## param G: random effects cov matrix
## param R_diag: error variances
## param rho: correlation between random errors
## param beta1: mean parameters for y1
## param beta2: mean parameters for y2
## param X0: design matrix for one individual with trt=0
## param X1: design matrix for on eindividual with trt=1
## param cntrl: nlme control
## return lst: a list containing model information for lme, gls models
##             plus other simulation information like rho and N
sim <- function(G,R_diag,rho,
                beta1,beta2,
                X0,X1,samples,cntrl){
  
  ## generate R matrix
  R <- diag(R_diag) %*% rbind(c(1,rho),c(rho,1)) %*% diag(R_diag)
  
  ## mean vectors
  mu0 <- t(cbind(X0 %*% beta1, X0 %*% beta2))
  mu1 <- t(cbind(X1 %*% beta1, X1 %*% beta2))
  
  ## simulate data
  n0 = samples/2
  n1 = samples/2
  dat0 <- gen_data_mixed(fixed_mean=mu0,G=G,R=R,
                         time_points=c(0,2,4,8),num_responses=2,
                         samples=n0)
  dat1 <- gen_data_mixed(fixed_mean=mu1,G=G,R=R,
                         time_points=c(0,2,4,8),num_responses=2,
                         samples=n1)
  
  ## put data in nice format
  dat <- rbind(dat0,dat1) %>% as.data.frame()
  colnames(dat) <- c("id","y_type","time","y")
  dat <- dat %>%
         mutate(id=c(rep(1:samples,each=8)),
                trt=c(rep(0,4*samples),rep(1,4*samples)),
                y_type=paste0("y",y_type),
                y_type_time = paste0(y_type,time)) %>%
         arrange(id,time,y_type)
  
  ## fit models 
  fit_lme <- lme(y~-1+y_type+y_type:(trt+time+trt:time),
              ## random slope and intercept for each y_type
              random=~-1+(y_type+y_type:time)|id,
              ## different error variances for different y_type
              weights=varIdent(form=~1|y_type),
              control=cntrl,
              data=dat)
  
  ## extract lme information
  V_lme <- getVarCov(fit_lme,type="marginal",individual=1)$`1`
  R_lme <- getVarCov(fit_lme,type="conditional",individual=1)$`1` 
  G_star <- getVarCov(fit_lme,type="random.effects")
  G_lme <- matrix(c(G_star[1,1],G_star[1,3],G_star[1,2],G_star[1,4],
                    G_star[3,1],G_star[3,3],G_star[3,2],G_star[3,4],
                    G_star[2,1],G_star[2,3],G_star[2,2],G_star[2,4],
                    G_star[4,1],G_star[4,3],G_star[4,2],G_star[4,4]),4,4,byrow=TRUE)
  beta_lme <- fit_lme$coefficients$fixed
  
  ## gls model
  fit_gls <- gls(y~-1+y_type+y_type:(trt+time+trt:time),
                 ## different error variances for different y_type_time
                 correlation=corSymm(form=~1|id),
                 weights=varIdent(form=~1|y_type_time),
                 control=cntrl,
                 data=dat)
  
  V_gls <- getVarCov(fit_gls,type="marginal",individual=1)
  beta_gls <- fit_gls$coefficients
  
  ## save outputs
  lst <- list(lme=list(V=V_lme,R=R_lme,G=G_lme,beta=beta_lme,rho=rho,N=samples),
              gls=list(V=V_gls,beta=beta_gls,rho=rho,N=samples))
  
  return(lst)

}


## rep_sim
## Purpose: wrapper function to replicate the simulation multiple times
##          with general across-simulation information inside
## param rho: correlation between random errors
## param N: number of subjects
## return a: a list of lists, each inner list containing a single simulation
rep_sim <- function(rho,N,reps=5){
  
  ## covariance matricies
  G <- rbind(c(0.959,0.049,0.730,-0.043),
             c(0.049,0.031,0.033,0.010),
             c(0.730,0.033,1.826,-0.025),
             c(-0.043,0.010,-0.025,0.008))
  
  R_diag <- c(sqrt(1.55),sqrt(0.571))
  
  ## mean parameters
  beta1 <- c(7.044,0.068,-0.056,-0.107)
  beta2 <- c(1.234,0.081,0.004,-0.081)
  
  ## design for one individual with treatment 0
  X0 <- rbind(c(1,0,0,0),
              c(1,0,2,0),
              c(1,0,4,0),
              c(1,0,8,0))
  
  ## design for one individual with treatment 1
  X1 <- rbind(c(1,1,0,0),
              c(1,1,2,2),
              c(1,1,4,4),
              c(1,1,8,8))
  
  cntrl <- lmeControl(maxIter=1000,
                      msMaxIter=1000,
                      niterEM=1000,
                      opt="optim")
  
  a <- replicate(try(sim(G=G,R_diag=R_diag,rho=rho,beta1=beta1,beta2=beta2,
                     X0=X0,X1=X1,samples=N,cntrl=cntrl)),n=reps)
  
  return(a)
  
}
