##
## Single Stock HMM Fitting ## 
## -------------------------------
if(!dir.exists("./output/models/")) dir.create("./output/models/")
if(!dir.exists("./output/models/hmm-ss/")) dir.create("./output/models/hmm-ss/")


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
  hmm_ac <- rstan::stan(file = "./stan/ss_hmm_1c.stan",
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
save(hmm_ac_out_sst, file="./output/hmm_ac_out_sst.Rdata")


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
  hmm_ac <- rstan::stan(file = "./stan/ss_hmm_1c.stan",
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
hmm_ac_out_2c <- list() # Master list to store outputs

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

hmm_ac_2c <- rstan::stan(file = "./stan/ss_hmm_2c.stan",
                      data = hmm.stan.dat,
                      warmup = 1000,
                      iter = 4000,
                      cores = 4,
                      chains = 4,
                      seed = 123,
                      control = list(adapt_delta = 0.95,
                                     max_treedepth = 10))
save(hmm_ac_2c, file = paste0("./output/models/hmm-ss/hmm_ac_2c_", levels(sock$Stock)[i], ".Rdata"))

hmm_ac_out_2c[[i]] <- rstan::summary(hmm_ac_2c, pars = params_out, probs=c(0.025, 0.2, 0.5, 0.8, 0.975))$summary

print(paste(i, "of", nlevels(sock$Stock)))

}


## Input-output HMM (IOHMM) : time-varying transition matrix
## ------------------------------------------------------------ ##

## SST as input vector, Comp as covariate ---

params_out <- c("beta1", "gamma", "alpha", "beta", "log_lik")
nlvl=2
iohmm_out <- list() # Master list length = N stocks


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
iohmm_out[[i]] <-  rstan::summary(iohmm, pars = params_out, probs=c(0.025, 0.2, 0.5, 0.8, 0.975))$summary
  
}


