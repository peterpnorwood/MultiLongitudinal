## ---------------------------------------------------- ##
## analyze_results.R ---------------------------------- ##
## Purpose: make plots and tables analyzing ----------- ##
## the results ---------------------------------------- ##
## Author: Peter Norwood, NCSU, Janssen R&D ----------- ##
## ---------------------------------------------------- ##

## set wd and bring in data
setwd("<LIBRARY>")
load("simulation_dat.RData")
load("simulation_pvals.RData")

## load packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(xtable)

## change/add variables to make them more plot-friendly 
big_dat <- simulation_dat %>%
           mutate(cor_char=paste0("Correlation=",rho),
                  N_char=paste0("N=",N),
                  model_char=case_when(model=="gls"~"gls",
                                        model=="mixed_diag_R"~"mixed_small",
                                        model=="mixed_unstr_R"~"mixed_large"))


big_dat$N_char <- factor(big_dat$N_char,levels=c("N=40",
                                    "N=100",
                                    "N=200",
                                    "N=300",
                                    "N=400"))

## true covariance parameters
t_var_y11=2.509
t_cov_y11_y21=0.730 ## plus covariance
t_cov_y11_y12=1.057
t_cov_y11_y22=0.645  
t_cov_y11_y13=1.155   
t_cov_y11_y23=0.559
t_cov_y11_y14=1.351
t_cov_y11_y24=0.387
t_var_y21=2.398       
  
## var(y11) sampling distribution plot
ggplot(data=big_dat) +
     geom_density(aes(x=var_y11,color=model_char),fill="white",alpha=0.5) +
  scale_color_manual(values=c("#012169","#CC0000","#4B9CD3")) +
     geom_vline(xintercept = t_var_y11) +
     facet_grid(cols=vars(cor_char),rows=vars(N_char)) +
     labs(x="Estimate",y="",color="Model") +
     ggtitle("Estimates of Var(y11)") +
     theme_bw() +
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())

## simulation data for y11
y11_dat <- big_dat %>%
           group_by(N,rho,model_char) %>%
           mutate(abs_error=abs(var_y11-t_var_y11),
                  sq_error=(var_y11-t_var_y11)**2) %>%
           summarise(mean_y11=mean(var_y11),
                     sd_y11=sd(var_y11),
                     mean_abs_error=mean(abs_error),
                     mse=mean(sq_error)) %>%
           arrange(N,rho,model_char)

## print for latex tables
print(xtable(y11_dat),include.rownames = FALSE)



## cov(y11,y22) sampling distribution plot
ggplot(data=big_dat) +
  geom_density(aes(x=cov_y11_y22,color=model_char),fill="white",alpha=0.5) +
  scale_color_manual(values=c("#012169","#CC0000","#4B9CD3")) +
  geom_vline(xintercept = t_cov_y11_y22) +
  facet_grid(cols=vars(cor_char),rows=vars(N_char)) +
  labs(x="Estimate",y="",color="Model") +
  ggtitle("Estimates of Cov(y11,y22)") +
  theme_bw() +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

## cov(y11,y21) sampling distribution plot
ggplot(data=big_dat) +
  geom_density(aes(x=cov_y11_y21,color=model_char),fill="white",alpha=0.5) +
  scale_color_manual(values=c("#012169","#CC0000","#4B9CD3")) +
  geom_vline(data=big_dat%>%filter(rho==0), aes(xintercept = t_cov_y11_y21)) +
  geom_vline(data=big_dat%>%filter(rho==0.25), aes(xintercept = t_cov_y11_y21+0.25*sqrt(t_var_e1)*sqrt(t_var_e2))) +
  geom_vline(data=big_dat%>%filter(rho==0.5), aes(xintercept = t_cov_y11_y21+0.5*sqrt(t_var_e1)*sqrt(t_var_e2))) +
  geom_vline(data=big_dat%>%filter(rho==0.75), aes(xintercept = t_cov_y11_y21+0.75*sqrt(t_var_e1)*sqrt(t_var_e2))) +
  facet_grid(cols=vars(cor_char),rows=vars(N_char)) +
  labs(x="Estimate",y="",color="Model") +
  ggtitle("Estimates of Cov(y11,y21)") +
  theme_bw() +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())


## summary of norms
summary_dat <- big_dat %>% group_by(N,cor_char,model_char) %>% 
  summarise(mean_V_norm=mean(V_norm),
  sd_V_norm=sd(V_norm),
  mean_G_norm=mean(G_norm,na.rm = TRUE),
  sd_G_norm=sd(G_norm,na.rm=TRUE),
  mean_beta_norm=mean(beta_norm),
  sd_beta_norm=sd(beta_norm))

## V norm plot
ggplot(data=summary_dat) +
  geom_line(aes(x=N,y=mean_V_norm,color=model_char)) +
  facet_grid(cols=vars(cor_char)) +
  scale_color_manual(values=c("#012169","#CC0000","#4B9CD3")) +
  labs(x="N",y="F Norm",color="Model") +
  ggtitle("Model Performance on Overall Covariance Matrix V") +
  theme_bw()

## beta norm plot
ggplot(data=summary_dat) +
  geom_line(aes(x=N,y=mean_beta_norm,color=model_char)) +
  facet_grid(cols=vars(cor_char)) +
  scale_color_manual(values=c("#012169","#CC0000","#4B9CD3")) +
  labs(x="N",y="2 Norm (SSE)",color="Model") +
  ggtitle("Model Performance on Fixed Effects") +
  theme_bw()

## G norm plot 
ggplot(data=summary_dat %>% filter(model_char != "gls")) +
  geom_line(aes(x=N,y=mean_G_norm,color=model_char)) +
  facet_grid(cols=vars(cor_char)) +
  scale_color_manual(values=c("#CC0000","#4B9CD3")) +
  labs(x="N",y="F Norm",color="Model") +
  ggtitle("Model Performance on Random Effects Covariance Matrix G") +
  theme_bw()


## estimation of conditional error
conditional <- big_dat %>% filter(model_char!="gls") %>%
               group_by(N,rho) %>%
               summarise(mean_e1=mean(var_e1),
                         sd_e1=sd(var_e1),
                         mean_e2=mean(var_e2),
                         sd_e2=sd(var_e2),
                         mean_cor_e11_e21=mean(cor_e11_e21,na.rm=TRUE),
                         sd_cor_e11_e21=sd(cor_e11_e21,na.rm = TRUE))

## true values for error variances
t_var_e1=1.55
t_var_e2=0.571

## var(e1) sampling distribution plot
ggplot(data=big_dat %>% filter(model_char!="gls")) +
  geom_density(aes(x=var_e1,color=model_char)) +
  scale_color_manual(values=c("#CC0000","#4B9CD3")) +
  geom_vline(xintercept = t_var_e1) +
  facet_grid(cols=vars(cor_char),rows=vars(N_char)) +
  labs(x="Estimate",y="",color="Model") +
  ggtitle("Estimates of Var(e1)") +
  theme_bw() +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

## var(e2) sampling distribution plot
ggplot(data=big_dat %>% filter(model!="gls")) +
  geom_density(aes(x=var_e2,color=model_char)) +
  scale_color_manual(values=c("#CC0000","#4B9CD3")) +
  geom_vline(xintercept = t_var_e2) +
  facet_grid(cols=vars(cor_char),rows=vars(N_char)) +
  labs(x="Estimate",y="",color="Model") +
  ggtitle("Estimates of Var(e2)") +
  xlim(c(0,1.25)) +
  theme_bw() +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

## sampling distribution of cor(e11,e21)
ggplot(data=big_dat %>% filter(model_char=="mixed_large")) +
  geom_density(aes(x=cor_e11_e21,color=model_char)) +
  scale_color_manual(values=c("#CC0000","#4B9CD3")) +
  facet_grid(cols=vars(cor_char),rows=vars(N_char)) +
  labs(x="Estimate",y="",color="Model") +
  ggtitle("Estimates of Cor(e11,e21)") +
  geom_vline(data=big_dat%>%filter(rho==0), aes(xintercept = rho)) +
  geom_vline(data=big_dat%>%filter(rho==0.25), aes(xintercept = rho)) +
  geom_vline(data=big_dat%>%filter(rho==0.5), aes(xintercept = rho)) +
  geom_vline(data=big_dat%>%filter(rho==0.75), aes(xintercept = rho)) +
  theme_bw() +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())


## more focus on mean parameters

## fixed effects parameters
t_beta10=7.044        
t_beta20=1.243
t_beta11=0.068
t_beta21=0.081
t_beta12=-0.056
t_beta22=-0.004
t_beta13=-0.107
t_beta23=-0.081

means <- big_dat %>%
  group_by(model_char,N,rho) %>%
  mutate(abs_error_beta12=abs(beta12-t_beta12),
         mse_beta12=(beta12-t_beta12)**2) %>%
  summarise(mean_beta12=mean(beta12),
            sd_beta12=sd(beta12),
            mean_abs_beta12=mean(abs_error_beta12),
            mean_mse_beta12=mean(mse_beta12)) %>%
  arrange(N,rho,model_char)


## beta12 sampling distribution plot
ggplot(data=big_dat) +
  geom_density(aes(x=beta12,color=model),fill="white",alpha=0.5) +
  geom_vline(xintercept = t_beta12) +
  facet_grid(cols=vars(cor_char),rows=vars(N_char)) +
  labs(x="Estimate",y="",color="Model") +
  ggtitle("Estimates of beta12") +
  theme_bw() +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

## table of beta12 estimates
xtable(means %>% select(N,rho,model_char,mean_beta12,sd_beta12,mean_abs_beta12,mean_mse_beta12) %>% filter(N==100))


## looking at AIC and BIC further
info_summary <- big_dat %>% group_by(N_char,rho,cor_char,model_char) %>%
                summarise(across(ends_with("IC"),list(mean=mean,sd=sd),.names="{fn}_{col}"))

info_summary$N_char <- factor(info_summary$N_char,levels=c("N=40",
                                                           "N=100",
                                                           "N=200",
                                                           "N=300",
                                                           "N=400"))
## AIC plot                 
ggplot(data=info_summary) +
  geom_line(aes(x=rho,y=mean_AIC,color=model_char)) +
  facet_wrap(vars(N_char),nrow=1,ncol=5,scales = "free_y") +
  scale_color_manual(values=c("#012169","#CC0000","#4B9CD3")) +
  labs(x="Error Correlation",y="AIC",color="Model") +
  ggtitle("AIC for Different Models") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

## BIC plot
ggplot(data=info_summary) +
  geom_line(aes(x=rho,y=mean_BIC,color=model_char)) +
  facet_wrap(vars(N_char),nrow=1,ncol=5,scales = "free_y") +
  scale_color_manual(values=c("#012169","#CC0000","#4B9CD3")) +
  labs(x="Error Correlation",y="BIC",color="Model") +
  ggtitle("BIC for Different Models") +
  theme_bw() 


## look at probability of rejecting h0 in likelihood ratio tests
pvals_summary <- simulation_pvals %>%
                 mutate(N_char = paste0("N=",N),
                        reject_lme_small_large=I(p_lme_small_large<0.05),
                        reject_lme_small_gls=I(p_lme_small_gls<0.05),
                        reject_gls_lme_large=I(p_gls_lme_large<0.05)) %>%
                 group_by(N_char,rho) %>%
                 summarise(across(starts_with("reject"),mean,.names="p_{col}"))


pvals_summary$N_char <- factor(pvals_summary$N_char,levels=c("N=40",
                                                 "N=100",
                                                 "N=200",
                                                 "N=300",
                                                 "N=400"))

## LRT plot one
ggplot(data=pvals_summary) +
  geom_line(aes(x=rho,y=p_reject_lme_small_large)) +
  geom_hline(yintercept = 0.05,linetype="dashed",color="red") +
  facet_wrap(vars(N_char),nrow=1,ncol=5) +
  labs(x="Error Correlation",y="Probability of Rejecting H0",color="Model") +
  ggtitle("Probability of Rejecting H0 in LRT Small vs Large Mixed Models") +
  theme_bw() #+
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


## LRT plot two
ggplot(data=pvals_summary) +
  geom_line(aes(x=rho,y=p_reject_lme_small_gls)) +
  geom_hline(yintercept = 0.05,linetype="dashed",color="red") +
  facet_wrap(vars(N_char),nrow=1,ncol=5) +
  labs(x="Error Correlation",y="Probability of Rejecting H0",color="Model") +
  ggtitle("Probability of Rejecting H0 in LRT Small Mixed vs Generalized Least Squares") +
  theme_bw()

## LRT plot three
ggplot(data=pvals_summary) +
  geom_line(aes(x=rho,y=p_reject_gls_lme_large)) +
  geom_hline(yintercept = 0.05,linetype="dashed",color="red") +
  facet_wrap(vars(N_char),nrow=1,ncol=5) +
  ylim(c(0,1)) +
  labs(x="Error Correlation",y="Probability of Rejecting H0",color="Model") +
  ggtitle("Probability of Rejecting H0 in LRT Generalized Least Squares vs Large Mixed") +
  theme_bw()
