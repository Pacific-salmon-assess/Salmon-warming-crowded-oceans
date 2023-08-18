##
## Single Stock HMM Fitting ## 
## -------------------------------

if(!dir.exists("./output/models/hmm-ss/")) dir.create("./output/models/hmm-ss/")
if(!dir.exists("./figures/hmm-ss/")) dir.create("./figures/hmm-ss/")


## single-covariate HMMs
## ------------------------------------------------------------ ##

params_out <- c("beta1", "gamma", "alpha", "beta", "log_lik", "rho") # Parameters of interest
nlvl = 2 # K=2 states


## SST ---

hmm_ac_out_sst <- list()  # Master list to store outputs

# Loop over stocks:
for(i in 1:nlevels(sock$Stock)){
  
  hmm.dat <- dplyr::filter(sock, Stock == levels(sock$Stock)[i]) # extract data for 1 stock
  
  
  hmm.stan.dat <- list( #Get stan data
    N = nrow(hmm.dat),
    R_S = hmm.dat$lnRS,
    S = hmm.dat$S,
    K = nlvl,
    X1 = hmm.dat$early_sst_stnd, # covariate = SST
    alpha_dirichlet = matrix(c(2,1,1,2), nrow=2) )
  
  
  # Run the stan model
  hmm_ac <- rstan::stan(file = "./stan/ss_hmm_c1.stan",
                       data = hmm.stan.dat,
                       warmup = 1000,
                       iter = 4000,
                       cores = 4,
                       chains = 4,
                       seed = 123,
                       control = list(adapt_delta = 0.95,
                                      max_treedepth = 10))
  save(hmm_ac, file = paste0("./output/models/hmm-ss/hmm_ac_sst_", levels(sock$Stock)[i], ".Rdata"))
  
  
  ## Extract parameters of interest
  
  hmm_ac_out_sst[[i]] <- rstan::summary(hmm_ac, pars = params_out, probs=c(0.025, 0.2, 0.5, 0.8, 0.975))$summary
  
}



## Comp ---

hmm_ac_out_comp <- list() # Master list to store outputs

# Loop over stocks
for(i in 1:nlevels(sock$Stock)){
  
  nlvl = 2 # 2 state model 
  
  hmm.dat <- dplyr::filter(sock, Stock == levels(sock$Stock)[i]) # extract data for 1 stock
  
  hmm.stan.dat <- list( #Get stan data
    N = nrow(hmm.dat),
    R_S = hmm.dat$lnRS,
    S = hmm.dat$S,
    K = nlvl,
    X1 = hmm.dat$np_pinks_sec_stnd, # covariate = Comp
    alpha_dirichlet = matrix(c(2,1,1,2), nrow=2) )
  
  
  # Run the stan model
  hmm_ac <- rstan::stan(file = "./stan/ss_hmm_c1.stan",
                        data = hmm.stan.dat,
                        warmup = 1000,
                        iter = 4000,
                        cores = 4,
                        chains = 4,
                        seed = 123,
                        control = list(adapt_delta = 0.95,
                                       max_treedepth = 10))
  save(hmm_ac, file = paste0("./output/models/hmm-ss/hmm_ac_comp_", levels(sock$Stock)[i], ".Rdata"))


  hmm_ac_out_comp[[i]] <- rstan::summary(hmm_ac, pars = params_out, probs=c(0.025, 0.2, 0.5, 0.8, 0.975))$summary
  
}



## Two- covariate HMM: SST and Comp
## -----------------------------------------------------

# AC 2-covar model

params_out <- c("beta1", "beta2", "gamma", "alpha", "beta", "log_lik")
nlvl = 2 # 2 state model
hmm_ac_out_c2 <- list() # Master list to store outputs

for(i in 1:nlevels(sock$Stock)) {
  
hmm.dat <- dplyr::filter(sock, Stock == levels(sock$Stock)[i])

hmm.stan.dat <- list( 
  N = nrow(hmm.dat),
  R_S = hmm.dat$lnRS,
  S = hmm.dat$S,
  K = nlvl,
  X1 = hmm.dat$early_sst_stnd, # covariate 1 = SST # determines states
  X2 =  hmm.dat$np_pinks_sec_stnd, # covariate 2 = Comp
  alpha_dirichlet = matrix(c(2,1,1,2), nrow=2) )

hmm_ac_c2 <- rstan::stan(file = "./stan/ss_hmm_c2.stan",
                      data = hmm.stan.dat,
                      warmup = 1000,
                      iter = 4000,
                      cores = 4,
                      chains = 4,
                      seed = 123,
                      control = list(adapt_delta = 0.95,
                                     max_treedepth = 10))
save(hmm_ac_c2, file = paste0("./output/models/hmm-ss/hmm_ac_c2_", levels(sock$Stock)[21], ".Rdata"))

post <- rstan::extract(hmm_ac_c2)

# beta1 : dim (6000, 2=K)
hmm_outputs[1,,] <- rep(apply(post$beta1, 2, median), times=hmm.stan.dat$N)
# beta2 : dim (6000, 2=K)
hmm_outputs[2,,] <- rep(apply(post$beta2, 2, median), times=hmm.stan.dat$N) 
# gamma : dim (6000, N years, 2=K)
hmm_outputs[3,,] <- t( apply(post$gamma, c(2,3), median) )
# alpha : dim (6000, N years, 2=K)
hmm_outputs[4,,] <- t( apply(post$alpha, c(2,3), median) )
# beta : dim (6000, N years, 2=K) #  
hmm_outputs[5,,] <- t( apply(post$beta, c(2,3), median) )
# log_lik : dim(6000, N years) 
hmm_outputs[6,,] <- rep( apply(post$log_lik, c(2), median), each=nlvl )


hmm_ac_out_c2[[i]] <- hmm_outputs

}


## Input-output HMM (IOHMM) : time-varying transition matrix
## ------------------------------------------------------------ ##

## SST as input vector, Comp as covariate ---

iohmm_out <- list() # Master list length = N stocks
nlvl=2

for(i in 1:nlevels(sock$Stock)){
  
  hmm.dat <- dplyr::filter(sock, Stock == levels(sock$Stock)[i]) # extract data for 1 stock
  
  hmm.stan.dat <- list( #Get stan data
    N = nrow(hmm.dat),
    R_S = hmm.dat$lnRS,
    S = hmm.dat$S,
    K = nlvl,
    u = hmm.dat$early_sst_stnd, # input vector = sst
    X1 = hmm.dat$np_pinks_sec_stnd, # covariate = Comp
    alpha_dirichlet = matrix(c(2,1,1,2), nrow=2) )
  
  
  # Run the stan model
  iohmm <- rstan::stan(file = "./stan/ss_iohmm.stan",
                       data = hmm.stan.dat,
                       warmup = 1000,
                       iter = 4000,
                       cores = 4,
                       chains = 4,
                       seed = 123,
                       control = list(adapt_delta = 0.95,
                                      max_treedepth = 10))
  save(iohmm, file = paste0("./output/models/hmm-ss/iohmm_c2", levels(sock$Stock)[i], ".Rdata"))
  
  
  ## Extract parameters of interest
  post <- rstan::extract(iohmm)
  params_out <- c("beta1", "gamma", "alpha", "beta", "log_lik")
  
  ## Store parameters of interest
  hmm_outputs <- array(NA, dim=c(length(params_out), nlvl, nrow(hmm.dat))) 
  
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
  
  
  iohmm_out[[i]] <- hmm_outputs
  
}


