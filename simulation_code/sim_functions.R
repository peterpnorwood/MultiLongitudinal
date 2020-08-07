## ---------------------------------------------------- ##
## sim_function.R ------------------------------------- ##
## Purpose: simulate data under a random slope, ------- ##
## random intercept model with normal errors ---------- ##
## Author: Peter Norwood, NCSU, Janssen --------------- ##
## ---------------------------------------------------- ##


library(Matrix)
library(nlme)
library(dplyr)
library(parallel)

setwd("<LIBRARY>")
source("gen_data.R")

## info_analysis
## Purpose: generate random slopes, random intercept data,
##          fit models, and extract model information
## param rho: correlation between random errors
## param N: number of subject sin a trial
## return lst: a list containing model information for lme, gls models
##             plus other simulation information like rho and N and hypothesis
##             testing results for three hypothesis tests
info_analysis <- function(rho,N){

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

  ## generate R matrix
  R <- diag(R_diag) %*% rbind(c(1,rho),c(rho,1)) %*% diag(R_diag)

  Z <- rbind(c(1,0,0,0),
             c(0,0,1,0),
             c(1,2,0,0),
             c(0,0,1,2),
             c(1,4,0,0),
             c(0,0,1,4),
             c(1,8,0,0),
             c(0,0,1,8))


  full_R <- bdiag(R,R,R,R)

  V <- Z %*% G %*% t(Z) + full_R

  ## mean vectors
  mu0 <- rbind(c(beta1[1],beta1[3]),
               c(beta2[1],beta2[3]))
  
  mu1 <- rbind(c(beta1[1]+beta1[2],beta1[3]+beta1[4]),
               c(beta2[1]+beta2[2],beta2[3]+beta2[4]))

  ## simulate data
  n0 = N/2
  n1 = N/2
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
    mutate(id=c(rep(1:N,each=8)),
           trt=c(rep(0,4*N),rep(1,4*N)),
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

  fit_lme2 <- lme(y~-1+y_type+y_type:(trt+time+trt:time),
                ## random slope and intercept for each y_type
                random=~-1+(y_type+y_type:time)|id,
                ## unstructured correlation matrix
                correlation=corSymm(form=~1|id),
                ## different error variances for different y_type
                weights=varIdent(form=~1|y_type),
                control=cntrl,
                data=dat)

  ## extract lme information
  V_lme2 <- getVarCov(fit_lme2,type="marginal",individual=1)$`1`
  R_lme2 <- getVarCov(fit_lme2,type="conditional",individual=1)$`1` 
  G_star2 <- getVarCov(fit_lme2,type="random.effects")
  G_lme2 <- matrix(c(G_star2[1,1],G_star2[1,3],G_star2[1,2],G_star2[1,4],
                  G_star2[3,1],G_star2[3,3],G_star2[3,2],G_star2[3,4],
                  G_star2[2,1],G_star2[2,3],G_star2[2,2],G_star2[2,4],
                  G_star2[4,1],G_star2[4,3],G_star2[4,2],G_star2[4,4]),4,4,byrow=TRUE)
  beta_lme2 <- fit_lme2$coefficients$fixed

  ## gls model
  fit_gls <- gls(y~-1+y_type+y_type:(trt+time+trt:time),
               ## different error variances for different y_type_time
               correlation=corSymm(form=~1|id),
               weights=varIdent(form=~1|y_type_time),
               control=cntrl,
               data=dat)
  
  ## extract gls information
  V_gls <- getVarCov(fit_gls,type="marginal",individual=1)
  beta_gls <- fit_gls$coefficients

  ## G norms
  lme_G_norm <- norm(G_lme-G,"F")
  lme2_G_norm <-norm(G_lme2-G,"F")

  ## V norms
  lme_V_norm <- norm(V_lme-V,"F")
  lme2_V_norm <- norm(V_lme2-V,"F")
  gls_V_norm <- norm(V_gls-V,"F")

  ## information criteria and LRT
  aov1 <- anova(fit_lme,fit_lme2)
  aov2 <- anova(fit_lme,fit_gls)
  aov3 <- anova(fit_gls,fit_lme2)
  
  df1 <- data.frame(N=rep(N,3),rho=rep(rho,3),model=c("mixed_diag_R","mixed_unstr_R","gls"),
                    V_norm=c(lme_V_norm,lme2_V_norm,gls_V_norm),
                    G_norm=c(lme_G_norm,lme2_G_norm,NA),
                    loglik=c(aov1$logLik[1],aov1$logLik[2],aov2$logLik[2]),
                    AIC=c(aov1$logLik[1],aov1$logLik[2],aov2$logLik[2]),
                    BIC=c(aov1$AIC[1],aov1$AIC[2],aov2$AIC[2]))
  
  df2 <- data.frame(N,rho,
                    p_lme_small_large=aov1$`p-value`[2],
                    p_lme_small_gls=aov2$`p-value`[2],
                    p_gls_lme_large=aov3$`p-value`[2])
  
  lst <- list(model_info=df1,
              p_vals=df2,
              lme_big_R=R_lme2,
              lme_R=R_lme,
              lme_G=G_lme,
              lme_big_G=G_lme2,
              lme_big_V=V_lme2,
              lme_V=V_lme,
              gls_V=V_gls,
              lme_big_beta=beta_lme2,
              lme_beta=beta_lme,
              gls_beta=beta_gls,
              rho=rho,
              N=N)
  
  return(lst)
}


## rep_sim
## Purpose: wrapper function to replicate the simulation multiple times
##          with general across-simulation information inside
## param r: correlation between random errors
## param n: number of subjects
## return a: a list of lists, each inner list containing a single simulation
rep_sim <- function(B,r,n){
  
  cl <- makeCluster(detectCores()-1)  
  # get library support needed to run the code
  clusterExport(cl,c('info_analysis','bdiag','gen_data_mixed','gen_data','%>%',
                     'lme','gls','getVarCov','mutate','arrange','norm','lmeControl','varIdent','corSymm'))
  # Set a different seed on each member of the cluster (just in case)
  clusterSetRNGStream(cl)
  #... then parallel replicate...
  a <- parSapply(cl, 1:B, function(i,r,n){try(info_analysis(rho=r,N=n))},r=r,n=n)
  #stop the cluster
  stopCluster(cl)
  return(a)
  
}

