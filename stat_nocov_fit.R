## Run a version of stationary models with no covariates - primarily for looking at residuals, in other scripts


# Species
if(speciesFlag=="pink") {
  data_master <- pink
  info_master <- pink.info} else if (speciesFlag=="pinkeven"){
    data_master <- pinkeven
    info_master <- pinkeven.info } else if (speciesFlag=="pinkodd"){
      data_master <- pinkodd
      info_master <- pinkodd.info } else if (speciesFlag=="chum") {
        data_master <- chum
        info_master <- chum.info } else if(speciesFlag=="sockeye"){
          data_master <- sock
          info_master <- sock.info }


# Set paths to output locations - dependent on species
fit.dir <- here("output", "models", "stat", speciesFlag) # place to store model fits

# Make them if they don't exist
if(!dir.exists(fit.dir))
  dir.create(fit.dir, recursive = TRUE)


## Monitor params
pars.stat <- c("alpha", "beta", "sigma", "phi", "mu_alpha", "sigma_alpha")
pars.gen.quant <- c("log_lik", "yhat", "yrep", "yresid") ## Generated quantities to monitor


## Run MCMC
stan.dat.nocov <- stan_data_stat(data_master,
                               scale.x1 = TRUE,
                               var.x2 = "early_sst_stnd",
                               var.x3 = "np_pinks_sec_stnd", # comp = pink abundance
                               var.region = "Ocean.Region2",
                               alpha.group = ifelse(speciesFlag=="sockeye", TRUE, FALSE)) # set to TRUE for sockeye
stat_nocov <- rstan::stan(file = "./stan/hbm_stat_nocov.stan",
                      data = stan.dat.nocov,
                      pars = c(pars.stat, pars.gen.quant),
                      warmup = 1000,
                      iter = 2000,
                      cores = 4,
                      chains = 4,
                      seed = 123,
                      control = list(adapt_delta = 0.99,
                                     max_treedepth = 20)) # increased treedepth
save(stat_nocov, file = here(fit.dir, "stat_nocov.RData"))
