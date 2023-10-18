## Fit Nonstationary hierarchical Bayesian models ("Era" and Random walk)
# This is adapted from: https://github.com/michaelmalick/sockeye-nonstationary


## Set Species -----
#speciesFlag = "pink"
speciesFlag = "chum"
#speciesFlag = "sockeye"

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
stan.dat.sst  <- stan_data_dyn(data_master, "early_sst_stnd",
                           breakpoint1 = 1977,
                           breakpoint2 = 1989,
                           var.region = "Ocean.Region2",
                           scale.x1 = TRUE)

stan.dat.comp <- stan_data_dyn(data_master, "np_pinks_sec_stnd",
                               breakpoint1 = 1977,
                               breakpoint2 = 1989,
                               var.region = "Ocean.Region2",
                               scale.x1 = TRUE)

stan.dat.2c <- stan_data_dyn(data_master, 
                             var.x2 = "early_sst_stnd",
                             var.x3 = "np_pinks_sec_stnd",
                             breakpoint1 = 1977,
                             breakpoint2 = 1989,
                             var.region="Ocean.Region2", 
                             scale.x1 = TRUE)


## Set pars to monitor -------------------------------------

pars_era_1c <- c("alpha", "beta", "sigma", "phi", "mu_alpha", "sigma_alpha",
                 "gamma1", "gamma2", "gamma3",
                 "mu_gamma1", "mu_gamma2", "mu_gamma3", "sigma_gamma")
pars_dyn_1c <- c("alpha", "beta", "sigma", "phi", "mu_alpha", "sigma_alpha",
                 "gamma", "sigma_gamma", "signal_noise")
pars_era_2c <- c("alpha", "beta", "sigma", "phi", "mu_alpha", "sigma_alpha",
                 "gamma1", "gamma2", "gamma3", 
                 "mu_gamma1", "mu_gamma2", "mu_gamma3", "sigma_gamma", 
                 "kappa1", "kappa2", "kappa3", 
                 "mu_kappa1", "mu_kappa2", "mu_kappa3", "sigma_kappa" )
pars_dyn_2c <- c("alpha", "beta", "sigma", "phi", "mu_alpha", "sigma_alpha",
                 "gamma", "sigma_gamma", "signal_noise_g",
                 "kappa", "sigma_kappa", "signal_noise_k")

pars.gen.quant <- c("log_lik", "yrep") ## Generated quantities to monitor


save(pars_era_1c, file = "./output/pars_era_1c.RData")
save(pars_dyn_1c, file = "./output/pars_dyn_1c.RData")
save(pars_era_2c, file = "./output/pars_era_2c.RData")
save(pars_dyn_2c, file = "./output/pars_dyn_2c.RData")

save(pars.gen.quant, file = "./output/pars_gen_quant.RData")



## Run MCMC  -----------------------------------------

## SST

era.sst <- rstan::stan(file = "./stan/hbm_era_1c.stan",
                        data = stan.dat.sst,
                        pars = c(pars_era_1c, pars.gen.quant),
                        warmup = 1000,
                        iter = 4000,
                        cores = 4,
                        chains = 4,
                        thin = 3,
                        seed = sample(1:1e6, 1),
                        control = list(adapt_delta = 0.98,
                                       max_treedepth = 10))
save(era.sst, file = here(fit.dir, "hbm_era_sst.RData"))

dyn.sst <- rstan::stan(file = "./stan/hbm_dyn_1c.stan",
                        data = stan.dat.sst,
                        pars = c(pars_dyn_1c, pars.gen.quant),
                        warmup = 1000,
                        iter = 4000,
                        cores = 4,
                        chains = 4,
                        thin = 3,
                        seed = sample(1:1e6, 1),
                        control = list(adapt_delta = 0.98,
                                       max_treedepth = 10))
save(dyn.sst, file = here(fit.dir, "hbm_dyn_sst.RData"))


## Comp

era.comp <- rstan::stan(file = "./stan/hbm_era_1c.stan",
                        data = stan.dat.comp,
                        pars = c(pars_era_1c, pars.gen.quant),
                        warmup = 1000,
                        iter = 4000,
                        cores = 4,
                        chains = 4,
                        thin = 3,
                        seed = sample(1:1e6, 1),
                        control = list(adapt_delta = 0.98,
                                       max_treedepth = 10))
save(era.comp, file = here(fit.dir, "hbm_era_comp.RData"))

dyn.comp <- rstan::stan(file = "./stan/hbm_dyn_1c.stan",
                        data = stan.dat.comp,
                        pars = c(pars_dyn_1c, pars.gen.quant),
                        warmup = 1000,
                        iter = 4000,
                        cores = 4,
                        chains = 4,
                        thin = 3,
                        seed = sample(1:1e6, 1),
                        control = list(adapt_delta = 0.98,
                                       max_treedepth = 10))
save(dyn.comp, file = here(fit.dir, "hbm_dyn_comp.RData"))


## 2-covariate

era.2c <- rstan::stan(file = "./stan/hbm_era_2c.stan",
                       data = stan.dat.2c,
                       pars = c(pars_era_2c, pars.gen.quant),
                       warmup = 1000,
                       iter = 4000,
                       cores = 4,
                       chains = 4,
                       thin = 3,
                       seed = sample(1:1e6, 1),
                       control = list(adapt_delta = 0.98,
                                      max_treedepth = 10))
save(era.2c, file = here(fit.dir, "hbm_era_2c.RData"))

dyn.2c <- rstan::stan(file = "./stan/hbm_dyn_2c.stan",
                      data = stan.dat.2c,
                      pars = c(pars_dyn_2c, pars.gen.quant),
                      warmup = 1000,
                      iter = 4000,
                      cores = 4,
                      chains = 4,
                      thin = 3,
                      seed = sample(1:1e6, 1),
                      control = list(adapt_delta = 0.98,
                                     max_treedepth = 10))
save(dyn.2c, file = here(fit.dir, "hbm_dyn_2c.Rdata"))

## Check pathology -----------------------------------------
# does not include 2-covar models yet

## NUTS diag

rstan::check_hmc_diagnostics(era.sst)
rstan::check_hmc_diagnostics(dyn.sst)

rstan::check_hmc_diagnostics(era.comp)
rstan::check_hmc_diagnostics(dyn.comp)

## Neff

neff_lowest(era.sst, pars = pars_era_1c)
neff_lowest(dyn.sst, pars = pars_dyn_1c)


neff_lowest(era.comp, pars = pars_era_1c)
neff_lowest(dyn.comp, pars = pars_dyn_1c)

## Rhat

rhat_highest(era.sst, pars = pars_era_1c)
rhat_highest(dyn.sst, pars = pars_dyn_1c)

rhat_highest(era.comp, pars = pars_era_1c)
rhat_highest(dyn.comp, pars = pars_dyn_1c)


## Pairs
#
# pairs_lowest(era.sst, pars = pars_era_c1)
# pairs_lowest(dyn.sst, pars = pars_dyn_c1)
## pairs_lowest(era.comp, pars = pars_era_c1)
# pairs_lowest(dyn.comp, pars = pars_dyn_c1)

## Run time

rstan::get_elapsed_time(era.sst)
rstan::get_elapsed_time(dyn.sst)

rstan::get_elapsed_time(era.comp)
rstan::get_elapsed_time(dyn.comp)


## MCMC diagnostics ----------------------------------------

## SST

pdf(here(fig.dir, "era_sst_diag.pdf"), width = 7, height = 5)
    coda_neff(get_neff(era.sst, pars = pars_era_1c), total_draws(era.sst))
    coda_rhat(get_rhat(era.sst, pars = pars_era_1c))
    coda_diag(As.mcmc.list(era.sst, pars = pars_era_1c))
dev.off()


pdf(here(fig.dir, "dyn_sst_diag.pdf"), width = 7, height = 5)
    coda_neff(get_neff(dyn.sst, pars = pars_dyn_1c), total_draws(dyn.sst))
    coda_rhat(get_rhat(dyn.sst, pars = pars_dyn_1c))
    coda_diag(As.mcmc.list(dyn.sst, pars = pars_dyn_1c))
dev.off()

## Comp

pdf(here(fig.dir, "era_comp_diag.pdf"), width = 7, height = 5)
  coda_neff(get_neff(era.comp, pars = pars_era_1c), total_draws(era.comp))
  coda_rhat(get_rhat(era.comp, pars = pars_era_1c))
  coda_diag(As.mcmc.list(era.comp, pars = pars_era_1c))
dev.off()


pdf(here(fig.dir, "dyn_comp_diag.pdf"), width = 7, height = 5)
  coda_neff(get_neff(dyn.comp, pars = pars_dyn_1c), total_draws(dyn.comp))
  coda_rhat(get_rhat(dyn.comp, pars = pars_dyn_1c))
  coda_diag(As.mcmc.list(dyn.comp, pars = pars_dyn_1c))
dev.off()


## Posterior predictive checks -----------------------------
plot_post_pc(era.sst, stan.dat.sst$y, pdf.path = here(fig.dir, "era_sst_yrep.pdf"))
plot_post_pc(dyn.sst, stan.dat.sst$y, pdf.path = here(fig.dir, "dyn_sst_yrep.pdf"))

plot_post_pc(era.comp, stan.dat.comp$y, pdf.path = here(fig.dir, "era_comp_yrep.pdf"))
plot_post_pc(dyn.comp, stan.dat.comp$y, pdf.path = here(fig.dir, "dyn_comp_yrep.pdf"))



## LOOIC ---------------------------------------------------

loo.era.sst <- rstan::loo(era.sst, cores = 4)
loo.dyn.sst <- rstan::loo(dyn.sst, cores = 4)


loo.era.comp <- rstan::loo(era.comp, cores = 4)
loo.dyn.comp <- rstan::loo(dyn.comp, cores = 4)



save(loo.era.sst, file = here(diag.dir, "loo_era_sst.RData"))
save(loo.dyn.sst, file = here(diag.dir, "loo_dyn_sst.RData"))

save(loo.era.comp, file = here(diag.dir, "loo_era_comp.RData"))
save(loo.dyn.comp, file = here(diag.dir, "loo_dyn_comp.RData"))



sum(pareto_k_values(loo.era.sst) > 0.7)
sum(pareto_k_values(loo.dyn.sst) > 0.7)

sum(pareto_k_values(loo.era.comp) > 0.7)
sum(pareto_k_values(loo.dyn.comp) > 0.7)



## SST

pdf(here(fig.dir, "era_sst_loo.pdf"), width = 7, height = 5)
    plot(loo.era.sst, label_points = TRUE)
dev.off()

pdf(here(fig.dir, "dyn_sst_loo.pdf"), width = 7, height = 5)
    plot(loo.dyn.sst, label_points = TRUE)
dev.off()

## Comp

pdf(here(fig.dir, "era_comp_loo.pdf"), width = 7, height = 5)
plot(loo.era.comp, label_points = TRUE)
dev.off()

pdf(here("dyn_comp_loo.pdf"), width = 7, height = 5)
plot(loo.dyn.comp, label_points = TRUE)
dev.off()

