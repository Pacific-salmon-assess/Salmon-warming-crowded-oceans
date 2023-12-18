## Fit Nonstationary hierarchical Bayesian models ("Era" and Random walk)
# This is adapted from: https://github.com/michaelmalick/sockeye-nonstationary

# Species
if(speciesFlag=="pink") 
  data_master <- pink else if(speciesFlag=="chum") 
    data_master <- chum else if(speciesFlag=="sockeye")
      data_master <- sock
    
# Set paths to output locations - dependent on species 
fig.dir <- here("figures", "dyn", speciesFlag, "hbm_fit") # place to store figures generated in this script
fit.dir <- here("output", "models", "dyn", speciesFlag) # place to store model fits
diag.dir <- here("output", "diagnostics", "dyn", speciesFlag) # place to store diagnostics
    

# Make them if they don't exist
if(!dir.exists(fig.dir))
  dir.create(fig.dir, recursive = T)

if(!dir.exists(fit.dir))
    dir.create(fit.dir, recursive = T)

if(!dir.exists(diag.dir))
    dir.create(diag.dir, recursive = T)


## Get data for Stan ---------------------------------------

stan.dat.2c <- stan_data_dyn(data_master, 
                             var.x2 = "early_sst_stnd",
                             var.x3 = "np_pinks_sec_stnd",
                             breakpoint1 = 1977,
                             breakpoint2 = 1989,
                             var.region="Ocean.Region2", 
                             scale.x1 = TRUE,
                             alpha.group = FALSE)


## Set pars to monitor -------------------------------------

pars_era_2c <- c("alpha", "beta", "sigma", "phi", "mu_alpha", "sigma_alpha",
                 "gamma1", "gamma2", "gamma3", 
                 "mu_gamma1", "mu_gamma2", "mu_gamma3", "sigma_gamma", 
                 "kappa1", "kappa2", "kappa3", 
                 "mu_kappa1", "mu_kappa2", "mu_kappa3", "sigma_kappa" )
pars_dyn_2c <- c("alpha", "beta", "sigma", "phi", "mu_alpha", "sigma_alpha",
                 "gamma", "sigma_gamma", "signal_noise_g",
                 "kappa", "sigma_kappa", "signal_noise_k")

pars.gen.quant <- c("log_lik", "yrep") ## Generated quantities to monitor


save(pars_era_2c, file = "./output/pars_era_2c.RData")
save(pars_dyn_2c, file = "./output/pars_dyn_2c.RData")

save(pars.gen.quant, file = "./output/pars_gen_quant.RData")



## Run MCMC  -----------------------------------------

## 2-covariate

era.2c <- rstan::stan(file = "./stan/hbm_era_2c.stan",
                       data = stan.dat.2c,
                       pars = c(pars_era_2c, pars.gen.quant),
                       warmup = 1000,
                       iter = 2000,
                       cores = 4,
                       chains = 4,
                       seed = 123,
                       control = list(adapt_delta = 0.99,
                                      max_treedepth = 20))
save(era.2c, file = here(fit.dir, "hbm_era_2c.RData"))

dyn.2c <- rstan::stan(file = "./stan/hbm_dyn_2c.stan",
                      data = stan.dat.2c,
                      pars = c(pars_dyn_2c, pars.gen.quant),
                      warmup = 1000,
                      iter = 4000, # try more iterations to fix Rhat & Neff
                      cores = 4,
                      chains = 4,
                      seed = 123,
                      control = list(adapt_delta = 0.99,
                                     max_treedepth = 20))
save(dyn.2c, file = here(fit.dir, "hbm_dyn_2c.Rdata"))



