## ---------------------------------------------------- ##
## fit_models.R --------------------------------------- ##
## Purpose: fit multivariate longitduinal models ------ ##
## on real clinical data ------------------------------ ##
## Author: Peter Norwood, NCSU, Janssen --------------- ##
## ---------------------------------------------------- ##

## Note: study and drug names redacted

setwd("~/Summer Project/Multivariate Simulation")

library(tidyr)
library(dplyr)
library(ggplot2)
library(taRifx)
library(nlme)

## read in the data
dat <- read.csv("<DATASET>.csv")

## function that will indicate if an entire column is null
not_all_na <- function(x) any(!is.na(x))

## transform data
dat <- dat %>% 
  ## select only certain colums
  select(ID,Treatment,Week,mMAYO,MAYO,pMAYO,bmkCRP,bmkCAL,bmkLAC) %>%
  ## filter to only weeks 0,2,4,8
  filter(Week<16) %>%
  ## transform biomarkers to be numeric and/or normally distributed
  mutate(mMayo=as.numeric(mMAYO),
         MAYO=as.numeric(MAYO),
         pMAYO=as.numeric(pMAYO),
         log_bmkCRP=log(bmkCRP),
         log_bmkCAL=log(bmkCAL),
         sqrt_bmkLAC=sqrt(bmkLAC)) 

## function that will indicate if an entire column is null
not_all_na <- function(x) any(!is.na(x))

## transform data
clean <- dat %>% 
  ## select only certain colums
  select(ID,Treatment,Week,mMAYO,MAYO,pMAYO,bmkCRP,bmkCAL,bmkLAC) %>%
  ## filter to only weeks 0,2,4,8
  filter(Week<16) %>%
  ## transform biomarkers to be numeric and/or normally distributed
  mutate(mMayo=as.numeric(mMAYO),
         MAYO=as.numeric(MAYO),
         pMAYO=as.numeric(pMAYO),
         log_bmkCRP=log(bmkCRP),
         log_bmkCAL=log(bmkCAL),
         sqrt_bmkLAC=sqrt(bmkLAC)) %>%
  ## widen the dataset
  pivot_wider(id_cols=ID,
              names_from=Week,
              values_from=c(mMAYO,MAYO,pMAYO,log_bmkCRP,log_bmkCAL,sqrt_bmkLAC),
              values_fn = list(MAYO = max, ##,
                               mMAYO=max,
                               pMAYO=max,
                               log_bmkCAL=max,
                               log_bmkCRP=max,
                               sqrt_bmkLAC=max)) %>%
  select_if(not_all_na)


## means of variables
means <- apply(clean[,2:ncol(clean)],2,mean,na.rm=TRUE)

## covariances and correlations
cov_mat <- cov(clean[,2:ncol(clean)],use="na.or.complete")
cor_mat <- cov2cor(cov_mat)


## MAYO 
MAYO_dat <- dat %>% mutate(y=MAYO,
                           y_type="MAYO",
                           y_type_time=paste0("MAYO_",Week)) %>%
  group_by(ID) %>%
  mutate(delta=y-first(y)) %>%
  select(ID,Treatment,Week,y,delta,y_type,y_type_time)

## pMAYO
pMAYO_dat <- dat %>% mutate(y=pMAYO,
                            y_type="pMAYO",
                            y_type_time=paste0("pMAYO_",Week)) %>%
  group_by(ID) %>%
  mutate(delta=y-first(y)) %>%
  select(ID,Treatment,Week,y,delta,y_type,y_type_time)



## combine all the data
clean_dat <- bind_rows(MAYO_dat,pMAYO_dat) %>%
  mutate(Treatment=case_when(Treatment=="Placebo"~"Placebo",
                             Treatment=="<Treatment One>"~"A",
                             Treatment=="<Treatment Two>"~"B"))


## control
cntrl <- lmeControl(maxIter=5000,
                    msMaxIter=5000,
                    niterEM=5000,
                    opt="optim")

cntrl2 <- lmeControl(maxIter=10000,
                     msMaxIter=10000,
                     niterEM=10000)

## fit pMAYO and MAYO

## small lme model
fit1_lme <- lme(y~-1+y_type+y_type:(Treatment+Week+Treatment:Week),
                ## random slope and intercept for each y_type
                random=~-1+(y_type+y_type:Week)|ID,
                control=cntrl,
                na.action = "na.omit",
                data= clean_dat %>% filter(y_type_time %in% c("MAYO_0","MAYO_8",
                                                              "pMAYO_0","pMAYO_2","pMAYO_4",
                                                              "pMAYO_8")) %>%
                  arrange(ID,Week,y_type))


## extract lme information
fit1_lme_beta <- fit1_lme$coefficients$fixed
fit1_lme_V <- getVarCov(fit1_lme,type="marginal",individual=1)$`CNTO1275UCO3001-100001`
fit1_lme_cor <- cov2cor(fit1_lme_V)
fit1_lme_G <- getVarCov(fit1_lme,type="random.effects")
fit1_lme_R <- getVarCov(fit1_lme,type="conditional")
sum1_lme <- summary(fit1_lme)

## compare lme covariance estimate to the sample covariance estimate
abs(cov_mat[c(3,5,6,7,4,8),c(3,5,6,7,4,8)]-fit1_lme_V)
abs(cor_mat[c(3,5,6,7,4,8),c(3,5,6,7,4,8)]-fit1_lme_cor)


## generalized least squares fit
fit1_gls <- gls(y~-1+y_type+y_type:(Treatment+Week+Treatment:Week),
                correlation=corSymm(form=~1|ID),
                weights=varIdent(form=~1|y_type_time),
                control=cntrl,
                na.action = "na.omit",
                data= clean_dat %>% filter(y_type_time %in% c("MAYO_0","MAYO_8",
                                                              "pMAYO_0","pMAYO_2","pMAYO_4",
                                                              "pMAYO_8")) %>%
                  arrange(ID,Week,y_type))

## extract gls model information
Gamma1 <- corMatrix(fit1_gls$modelStruct$corStruct)[[1]]
parms1 <- coef(fit1_gls$model,unconstrained=FALSE)   
vars1 <- ( c(1,parms1[29:33])*fit1_gls$sigma)^2
fit1_gls_V <- diag(sqrt(vars1))%*%Gamma1%*%diag(sqrt(vars1))
fit1_gls_cor <- cov2cor(fit1_gls_V)
colnames(fit1_gls_cor) <- c("MAYO_0","pMAYO_0","pMAYO_2","pMAYO_4","MAYO_8","pMAYO_8")
rownames(fit1_gls_cor) <- colnames(fit1_gls_cor)
sum1_gls <- summary(fit1_gls)

## estimates for fixed effects
## mixed model
round(sum1_lme$tTable[,-3],3)
## gls model
round(sum1_gls$tTable,3)

## variance of fixed effects
norm(fit1_lme$varFix,"F")
norm(fit1_gls$varBeta,"F")

norm(fit1_lme$varFix,"2")
norm(fit1_gls$varBeta,"2")

## covariance and correlation matricies for random effects
fit1_lme_G
cov2cor(fit1_lme_G)

## conditional covariance
fit1_lme_R

## marginal covariance
fit1_lme_V
fit1_gls_V
## marginal correlation
cov2cor(fit1_lme_V)
cov2cor(fit1_gls_V)

## gather correlation matrix
Gamma1 <- corMatrix(fit1_gls$modelStruct$corStruct)[[1]]
parms1 <- coef(fit1_gls$model,unconstrained=FALSE)   
vars1 <- ( c(1,parms1[29:33])*fit1_gls$sigma)^2
fit1_gls_V <- diag(sqrt(vars1))%*%Gamma1%*%diag(sqrt(vars1))
fit1_gls_cor <- cov2cor(fit1_gls_V)
colnames(fit1_gls_cor) <- c("MAYO_0","pMAYO_0","pMAYO_2","pMAYO_4","MAYO_8","pMAYO_8")
rownames(fit1_gls_cor) <- colnames(fit1_gls_cor)
## print correlation matrix
fit1_gls_cor


## gathering information criteria and the LRT
anova(fit1_lme,fit1_gls)


## hypothesis testing on MAYO
trt_130_delta <- c(1,0,1,0,0,0,8,0,8,0,0,0) - c(1,0,1,0,0,0,0,0,0,0,0,0)
trt_6_delta <-   c(1,0,0,0,1,0,8,0,0,0,8,0) - c(1,0,0,0,1,0,0,0,0,0,0,0)
trt_placebo_delta <- c(1,0,0,0,0,0,8,0,0,0,0,0) - c(1,0,0,0,0,0,0,0,0,0,0,0)

## contrasts
c1 <- trt_130_delta - trt_placebo_delta
c2 <- trt_6_delta - trt_placebo_delta
c3 <- trt_130_delta - trt_6_delta

## mean and covariance estimates
est_coef <- fit1_gls$coefficients
var_coef <- fit1_gls$varBeta

## test one
est1 <- c1 %*% est_coef
se1 <- sqrt(t(c1) %*% var_coef %*% c1)
z1 <- est1/se1
p1 <- pnorm(z1) 

## test two
est2 <- c2 %*% est_coef
se2 <- sqrt(t(c2) %*% var_coef %*% c2)
z2 <- est2/se2
p2 <- pnorm(z2) 

## test three
est3 <- c3 %*% est_coef
se3 <- sqrt(t(c3) %*% var_coef %*% c3)
z3 <- est3/se3
p3 <- 2*(1-pnorm(z3)) 

## hypothesis testing on pmayo at week 4
pMAYO_trt_130_delta <- c(0,1,0,1,0,0,0,4,0,4,0,0) - c(0,1,0,1,0,0,0,0,0,0,0,0)
pMAYO_trt_6_delta <-   c(0,1,0,0,0,1,0,4,0,0,0,4) - c(0,1,0,0,0,1,0,0,0,0,0,0)
pMAYO_trt_placebo_delta <- c(0,1,0,0,0,0,0,4,0,0,0,0) - c(0,1,0,0,0,0,0,0,0,0,0,0)

## contrasts
pMAYO_c1 <- pMAYO_trt_130_delta - pMAYO_trt_placebo_delta
pMAYO_c2 <- pMAYO_trt_6_delta - pMAYO_trt_placebo_delta
pMAYO_c3 <- pMAYO_trt_130_delta - pMAYO_trt_6_delta

## test one
pMAYO_est1 <- pMAYO_c1 %*% est_coef
pMAYO_se1 <- sqrt(t(pMAYO_c1) %*% var_coef %*% pMAYO_c1)
pMAYO_z1 <- pMAYO_est1/pMAYO_se1
pMAYO_p1 <- pnorm(pMAYO_z1) 

## test two
pMAYO_est2 <- pMAYO_c2 %*% est_coef
pMAYO_se2 <- sqrt(t(pMAYO_c2) %*% var_coef %*% pMAYO_c2)
pMAYO_z2 <- pMAYO_est2/pMAYO_se2
pMAYO_p2 <- pnorm(pMAYO_z2) 

## test three
pMAYO_est3 <- pMAYO_c3 %*% est_coef
pMAYO_se3 <- sqrt(t(pMAYO_c3) %*% var_coef %*% pMAYO_c3)
pMAYO_z3 <- pMAYO_est3/pMAYO_se3
pMAYO_p3 <- 2*(1-pnorm(pMAYO_z3))
