## HMM Fitting ## 
## Fitting to each stock individually


if(!dir.exists("./output/models/hmm-ss/")) dir.create("./figures/stat/pub/")
if(!dir.exists("./figures/hmm-ss/")) dir.create("./figures/stat/pub/")


# Parameters of interest
params_out <- c("beta1", "gamma", "alpha", "beta", "log_lik")



## hmm_ss1 [Single Stock version 1] : 2 states, original transition matrix (2,1,1,2)
## --------------------------------------------------------------------------------

# list to store outputs
hmm_out_ss1 <- list()

#loop over all stocks
for ( i in 1:nlevels(sock$Stock) ) {
  
  # Number of states
  nlvl = 2 
  
  # Filter to one stock
  hmm.dat <- dplyr::filter(sock, Stock == levels(sock$Stock)[i]) 
  
  #output array 
  hmm_outputs <- array(NA, dim=c(length(params_out), nlvl, nrow(hmm.dat)))
  dimnames(hmm_outputs) <- list( params_out, c("State 1", "State 2"), hmm.dat$BY)
  
  
  # Define data
  hmm.stan.dat <- list( 
       N = nrow(hmm.dat),
       R_S = hmm.dat$lnRS,
       S = hmm.dat$S,
       K = nlvl,
       X1 = hmm.dat$early_sst_stnd, # covariate = SST
       alpha_dirichlet = matrix(c(2,1,1,2), nrow=2) )
  
  
  # Kept all the specs from hmm.fit script ??
  hmm_ss1 <- rstan::stan(file = "./stan/hmm_test.stan",
                       data = hmm.stan.dat,
                       warmup = 1000,
                       iter = 4000,
                       cores = 4,
                       chains = 4,
                       seed = 123,
                       control = list(adapt_delta = 0.90,
                                      max_treedepth = 10))
  save(hmm_ss1, file = paste0("./output/models/hmm-ss/hmm_ss1_", levels(sock$Stock)[i], ".Rdata"))
  
  ## Extract parameters of interest
  post=rstan::extract(hmm_ss1)
  
  
  # beta1 : dim (6000, 2=K)
  hmm_outputs[1,,] <- rep(apply(post$beta1, 2, median), times=hmm.stan.dat$N) 
  # gamma : dim (6000, N years, 2=K)
  hmm_outputs[2,,] <- t( apply(post$gamma, c(2,3), median) )
  # alpha : dim (6000, N years, 2=K)
  hmm_outputs[3,,] <- t( apply(post$alpha, c(2,3), median) )
  # beta : dim (6000, N years, 2=K) #  
  hmm_outputs[4,,] <- t( apply(post$beta, c(2,3), median) )
  # log_lik : dim(6000, N years) 
  hmm_outputs[5,,] <- rep( apply(post$log_lik, c(2), median), each=nlvl )
  
  
  hmm_out_ss1[[i]] <- hmm_outputs


}


## hmm_ss2 [Single Stock version 2] : 3 states, 2,1,1 transition matrix
## --------------------------------------------------------------------------------

# list to store outputs
hmm_out_ss2 <- list()

#loop over all stocks
for ( i in 1:nlevels(sock$Stock) ) {
  
  # Number of states
  nlvl = 3 
  
  # Filter to one stock
  hmm.dat <- dplyr::filter(sock, Stock == levels(sock$Stock)[i]) 
  
  #output array 
  hmm_outputs <- array(NA, dim=c(length(params_out), nlvl, nrow(hmm.dat)))
  dimnames(hmm_outputs) <- list( params_out, c("State 1", "State 2", "State 3"), hmm.dat$BY)
  
  
  # Define data
  hmm.stan.dat <- list( 
    N = nrow(hmm.dat),
    R_S = hmm.dat$lnRS,
    S = hmm.dat$S,
    K = nlvl,
    X1 = hmm.dat$early_sst_stnd, # covariate = SST
    alpha_dirichlet = matrix(c(2,1,1,1,2,1,1,1,2), nrow=3) )
  
  
  # Kept all the specs from hmm.fit script ??
  hmm_ss2 <- rstan::stan(file = "./stan/hmm_test.stan",
                          data = hmm.stan.dat,
                          warmup = 1000,
                          iter = 4000,
                          cores = 4,
                          chains = 4,
                          seed = 123,
                          control = list(adapt_delta = 0.90,
                                         max_treedepth = 10))
  save(hmm_ss2, file = paste0("./output/models/hmm-ss/hmm_ss2_", levels(sock$Stock)[i], ".Rdata"))
  
  ## Extract parameters of interest
  post=rstan::extract(hmm_ss2)
  
  
  # beta1 : dim (6000, 2=K)
  hmm_outputs[1,,] <- rep(apply(post$beta1, 2, median), times=hmm.stan.dat$N) 
  # gamma : dim (6000, N years, 2=K)
  hmm_outputs[2,,] <- t( apply(post$gamma, c(2,3), median) )
  # alpha : dim (6000, N years, 2=K)
  hmm_outputs[3,,] <- t( apply(post$alpha, c(2,3), median) )
  # beta : dim (6000, N years, 2=K) #  
  hmm_outputs[4,,] <- t( apply(post$beta, c(2,3), median) )
  # log_lik : dim(6000, N years) 
  hmm_outputs[5,,] <- rep( apply(post$log_lik, c(2), median), each=nlvl )
  
  
  hmm_out_ss2[[i]] <- hmm_outputs
  
  
}



## hmm_ss3 [Single Stock version 3] : 2 states, conservative transition matrix (4,1,1,4)
## --------------------------------------------------------------------------------

# list to store outputs
hmm_out_ss3 <- list()

#loop over all stocks
for ( i in 1:nlevels(sock$Stock) ) {
  
  # Number of states
  nlvl = 2 
  
  # Filter to one stock
  hmm.dat <- dplyr::filter(sock, Stock == levels(sock$Stock)[i]) 
  
  #output array 
  hmm_outputs <- array(NA, dim=c(length(params_out), nlvl, nrow(hmm.dat)))
  dimnames(hmm_outputs) <- list( params_out, c("State 1", "State 2"), hmm.dat$BY)
  
  
  # Define data
  hmm.stan.dat <- list( 
    N = nrow(hmm.dat),
    R_S = hmm.dat$lnRS,
    S = hmm.dat$S,
    K = nlvl,
    X1 = hmm.dat$early_sst_stnd, # covariate = SST
    alpha_dirichlet = matrix(c(4,1,1,4), nrow=2) )
  
  
  # Kept all the specs from hmm.fit script ??
  hmm_ss3 <- rstan::stan(file = "./stan/hmm_test.stan",
                         data = hmm.stan.dat,
                         warmup = 1000,
                         iter = 4000,
                         cores = 4,
                         chains = 4,
                         seed = 123,
                         control = list(adapt_delta = 0.90,
                                        max_treedepth = 10))
  save(hmm_ss3, file = paste0("./output/models/hmm-ss/hmm_ss3_", levels(sock$Stock)[i], ".Rdata"))
  
  ## Extract parameters of interest
  post=rstan::extract(hmm_ss3)
  
  
  # beta1 : dim (6000, 2=K)
  hmm_outputs[1,,] <- rep(apply(post$beta1, 2, median), times=hmm.stan.dat$N) 
  # gamma : dim (6000, N years, 2=K)
  hmm_outputs[2,,] <- t( apply(post$gamma, c(2,3), median) )
  # alpha : dim (6000, N years, 2=K)
  hmm_outputs[3,,] <- t( apply(post$alpha, c(2,3), median) )
  # beta : dim (6000, N years, 2=K) #  
  hmm_outputs[4,,] <- t( apply(post$beta, c(2,3), median) )
  # log_lik : dim(6000, N years) 
  hmm_outputs[5,,] <- rep( apply(post$log_lik, c(2), median), each=nlvl )
  
  
  hmm_out_ss3[[i]] <- hmm_outputs
  
  
}



## Input-output HMM (IOHMM) test : varying transition matrix
## ------------------------------------------------------------ ##

iohmm_out_ss1 <- list() # Master list length = N stocks


for(i in 1:nlevels(sock$Stock)){
  
  nlvl = 2 # 2 state model 
    
  hmm.dat <- dplyr::filter(sock, Stock == levels(sock$Stock)[i]) # extract data for 1 stock
  
  
  hmm.stan.dat <- list( #Get stan data
    N = nrow(hmm.dat),
    R_S = hmm.dat$lnRS,
    S = hmm.dat$S,
    K = nlvl,
    X1 = hmm.dat$early_sst_stnd, # covariate = SST
    alpha_dirichlet = matrix(c(2,1,1,2), nrow=2) )
  
  
  # Run the stan model
  iohmm <- rstan::stan(file = "./stan/iohmm_test.stan",
                         data = hmm.stan.dat,
                         warmup = 1000,
                         iter = 4000,
                         cores = 4,
                         chains = 4,
                         seed = 123,
                         control = list(adapt_delta = 0.90,
                                        max_treedepth = 10))
  save(iohmm, file = paste0("./output/models/hmm-ss/iohmm_ss1_", levels(sock$Stock)[i], ".Rdata"))
  
  
  ## Extract parameters of interest
  post <- rstan::extract(iohmm)
  params_out <- c("beta1", "gamma", "alpha", "beta", "log_lik")
  
  ## Store parameters of interest
  hmm_outputs <- array(NA, dim=c(length(params_out), nlvl, nrow(hmm.dat))) # Stock-specific array of results - overwritten each loop
  dimnames(hmm_outputs) <- list( params_out, c("State 1", "State 2"), hmm.dat$BY)
  
  # beta1 : dim (6000, 2=K)
  hmm_outputs[1,,] <- rep( apply(post$beta1, 2, median), times=hmm.stan.dat$N) 
  # gamma : dim (6000, N years, 2=K)
  hmm_outputs[2,,] <- t( apply(post$gamma, c(2,3), median) )
  # alpha : dim (6000, N years, 2=K)
  hmm_outputs[3,,] <- t( apply(post$alpha, c(2,3), median) )
  # beta : dim (6000, N years, 2=K) 
  hmm_outputs[4,,] <- t( apply(post$beta, c(2,3), median) )
  # log_lik : dim(6000, N years) 
  hmm_outputs[5,,] <- rep( apply(post$log_lik, c(2), median), each=nlvl )
  
  
  iohmm_out_ss1[[i]] <- hmm_outputs

}


