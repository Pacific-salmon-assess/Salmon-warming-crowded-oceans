## Stationary model fits (new) ##


# Species
if(speciesFlag=="pink") {
  data_master <- pink
  info_master <- pink.info} else if(speciesFlag=="chum") {
    data_master <- chum
    info_master <- chum.info } else if(speciesFlag=="sockeye"){
      data_master <- sock
      info_master <- sock.info }


# Set paths to output locations - dependent on species 
fig.dir <- here("figures", "stat", speciesFlag, "hbm_fit") # place to store figures generated in this script
fit.dir <- here("output", "models", "stat", speciesFlag) # place to store model fits
diag.dir <- here("output", "diagnostics", "stat", speciesFlag) # place to store diagnostics

# Make them if they don't exist 
if(!dir.exists(fig.dir))
  dir.create(fig.dir, recursive = TRUE) 
if(!dir.exists(fit.dir))
  dir.create(fit.dir, recursive = TRUE) 
if(!dir.exists(diag.dir))
  dir.create(diag.dir, recursive = TRUE) 


## Stationary model with SUBSETTED stocks [Sockeye ONLY for now] ------------------##

# subset data    
stk.sub <- info_master$Stock.ID[info_master$yr_start < 1985 & info_master$yr_end >= 2014]
data_sub <- data_master %>% filter(Stock.ID %in% stk.sub)
info_sub <- info_master %>% filter(Stock.ID %in% stk.sub)
data_sub <- geographic.order(data_sub)
info_sub <- geographic.order(info_sub)

# set up MCMC
pars.stat <- c("alpha", "beta", "sigma", "phi", "mu_alpha", "sigma_alpha",
               "gamma", "mu_gamma", "sigma_gamma",
               "kappa", "mu_kappa", "sigma_kappa")
pars.gen.quant <- c("log_lik", "yhat", "yrep", "yresid") ## Generated quantities to monitor


## Run MCMC
stan.dat.sub <- stan_data_stat(data_sub,
                               scale.x1 = TRUE,
                               var.x2 = "early_sst_stnd",
                               var.x3 = "np_pinks_sec_stnd", # comp = pink abundance
                               var.region = "Ocean.Region2",
                               alpha.group = ifelse(speciesFlag=="sockeye", TRUE, FALSE)) # set to TRUE for sockeye
stat_sub <- rstan::stan(file = "./stan/hbm_stat_2c.stan",
                      data = stan.dat.sub,
                      pars = c(pars.stat, pars.gen.quant),
                      warmup = 1000,
                      iter = 2000,
                      cores = 4,
                      chains = 1,
                      seed = 123,
                      control = list(adapt_delta = 0.99,
                                     max_treedepth = 20)) # increased treedepth
save(stat_sub, file = here(fit.dir, "stat_sub.RData"))


