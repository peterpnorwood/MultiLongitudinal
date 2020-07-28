## ---------------------------------------------------- ##
## gen_data.R ----------------------------------------- ##
## Purpose: simulate data under a random slope, ------- ##
## random intercept model with normal errors ---------- ##
## Author: Peter Norwood, NCSU, Janssen --------------- ##
## ---------------------------------------------------- ##



## gen_data
## Purpose: generate correlated longitudinal data
## param mean_vec: vector of means
## param cov: covariance matrix
## param samples: number of samples to take
## return Y: a vector of responses for an individual
gen_data <- function(mean_vec,cov_mat,samples){
  
  ## generate correlation matrix
  cor_mat <- cov2cor(cov_mat)
  ## generate diag of cov matrix
  cov_diag <- diag(diag(cov_mat))
  ## cholesky factorization of corr
  A <- t(chol(cor_mat))
  ## multiplying matrix
  B <- cov_diag**(0.5) %*% A
  Binv <- solve(B)
  
  ## mean x vector
  mu_x <- Binv %*% mean_vec
  
  ## generate random x vector 
  X <- sapply(mu_x,rnorm,n=samples,sd=1)
  
  ## generate random y vector
  Y <- t(B %*% t(X))
  
  return(Y)
}

## gen_data_mixed
## Purpose: generate correlated longitudinal data with
## random slopes and random intercepts
## param fixed_mean: vector of means parameters, 
## rows indicate different responses
## param G: covariance matrix for the random effects
## param R: covariance matrix for the random error
## param time_points: number of time points the responses have
## num_responses: number of responses
## param samples: number of samples to take
## return dat: a matrix with the following columns: ID, y_type, time, y 
gen_data_mixed <- function(fixed_mean,
                           G,R,
                           time_points,
                           num_responses,
                           samples){
  
  
  random_mean <- c(rep(0,ncol(G)))
  random_coef <- gen_data(mean_vec=random_mean,
                          cov_mat=G,samples=samples)
  
  rand_list <- list()
  splits <- ncol(G)/num_responses
  start <- 1
  for(r in 1:splits){
    #print(start:(start+(num_responses-1)))
    rand_list[[r]] <- random_coef[,start:(start+(num_responses-1))]
    start <- start+num_responses
  }
  
  error_mean <- c(rep(0,ncol(R)))
  error_list <- list()
  for(r in 1:length(time_points)){
    error_list[[r]] <- gen_data(mean_vec=error_mean,
                                cov_mat=R,samples=samples)
  }
  
  ## add the slopes to the correct columns
  dat <- matrix(NA,nrow=length(time_points)*num_responses*samples,ncol=4)
  tick <- 1
  ## loop through all individuals
  for(i in 1:samples){
    ## loop through different responses
    for(j in 1:num_responses){
      ## loop through time points
      for(k in 1:length(time_points)){
        t <- time_points[k]
        ## int j-th response    ## random int j-th response, i-th individual
        y <- (fixed_mean[j,1] + rand_list[[j]][i,1]) + 
          ## slope j-th response  ## random slope j-th response, i-th individual
          (fixed_mean[j,2] + rand_list[[j]][i,2])*t +
          ## random error k-th time point, j-th response, i-th individual
          error_list[[k]][i,j]
        
        dat[tick,] <- c(i,j,t,y)
        tick=tick+1
        
      }
    }
  }
  
  return(dat)
}
