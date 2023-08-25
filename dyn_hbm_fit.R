## Fit Nonstationary hierarchical Bayesian models ("Era" and Random walk)
# This is adapted from: https://github.com/michaelmalick/sockeye-nonstationary

if(!dir.exists("./figures/dyn/hbm_fit/"))
    dir.create("./figures/dyn/hbm_fit/")

if(!dir.exists("./output/models/dyn/"))
    dir.create("./output/models/dyn/")

if(!dir.exists("./output/diagnostics/dyn/"))
    dir.create("./output/diagnostics/dyn/")


## Get data for Stan ---------------------------------------
stan.dat.sst  <- stan_data_dyn(sock, "early_sst_stnd",
                           breakpoint1 = 1977,
                           breakpoint2 = 1989,
                           var.region = "Ocean.Region2",
                           scale.x1 = TRUE)
save(stan.dat.sst, file = "./output/stan_dat_sst.RData")

stan.dat.comp <- stan_data_dyn(sock, "np_pinks_sec_stnd",
                               breakpoint1 = 1977,
                               breakpoint2 = 1989,
                               var.region = "Ocean.Region2",
                               scale.x1 = TRUE)
save(stan.dat.comp, file = "./output/stan_dat_comp.RData")

stan.dat.2c <- stan_data_dyn(sock, 
                             var.x2 = "early_sst_stnd",
                             var.x3 = "np_pinks_sec_stnd",
                             breakpoint1 = 1977,
                             breakpoint2 = 1989,
                             var.region="Ocean.Region2", 
                             scale.x1 = TRUE)

save(stan.dat.2c, file = "./output/stan_dat_2c.RData")


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
                 "gamma", "sigma_gamma", "signal_noise_gamma",
                 "kappa", "sigma_kappa", "signal_noise_kappa")

pars.gen.quant <- c("log_lik", "yrep") ## Generated quantities to monitor


save(pars_era_1c, file = "./output/pars_era_1c.RData")
save(pars_dyn_1c, file = "./output/pars_dyn_1c.RData")
save(pars_era_2c, file = "./output/pars_era_2c.RData")
save(pars_dyn_2c, file = "./output/pars_dyn_2c.RData")

save(pars.gen.quant, file = "./output/pars_gen_quant.RData")



## Run MCMC  -----------------------------------------

## SST

era.sst <- rstan::stan(file = "./stan/hbm_era_1c",
                        data = stan.dat.sst,
                        pars = c(pars_era_1c, pars.gen.quant),
                        warmup = 1000,
                        iter = 4000,
                        cores = 4,
                        chains = 4,
                        thin = 3,
                        seed = sample(1:1e6, 1),
                        control = list(adapt_delta = 0.90,
                                       max_treedepth = 10))
save(era.sst, file = "./output/models/dyn/hbm_era_sst.RData")

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
save(dyn.sst, file = "./output/models/dyn/hbm_dyn_sst.RData")


## Comp

era.comp <- rstan::stan(file = "./stan/hbm_era_1c",
                        data = stan.dat.comp,
                        pars = c(pars_era_1c, pars.gen.quant),
                        warmup = 1000,
                        iter = 4000,
                        cores = 4,
                        chains = 4,
                        thin = 3,
                        seed = sample(1:1e6, 1),
                        control = list(adapt_delta = 0.90,
                                       max_treedepth = 10))
save(era.comp, file = "./output/models/dyn/hbm_era_comp.RData")

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
save(dyn.comp, file = "./output/models/dyn/hbm_dyn_comp.RData")


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
                       control = list(adapt_delta = 0.90,
                                      max_treedepth = 10))
save(era.2c, file = "./output/models/dyn/hbm_era_2c.RData")

## Check pathology -----------------------------------------

## NUTS diag

rstan::check_hmc_diagnostics(era.sst)
rstan::check_hmc_diagnostics(dyn.sst)

rstan::check_hmc_diagnostics(era.comp)
rstan::check_hmc_diagnostics(dyn.comp)

## Neff

neff_lowest(era.sst, pars = pars_era_c1)
neff_lowest(dyn.sst, pars = pars_dyn_c1)


neff_lowest(era.comp, pars = pars_era_c1)
neff_lowest(dyn.comp, pars = pars_dyn_c1)

## Rhat

rhat_highest(era.sst, pars = pars_era_c1)
rhat_highest(dyn.sst, pars = pars_dyn_c1)

rhat_highest(era.comp, pars = pars_era_c1)
rhat_highest(dyn.comp, pars = pars_dyn_c1)


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

pdf("./figures/dyn/hbm_fit/era_sst_diag.pdf", width = 7, height = 5)
    coda_neff(get_neff(era.sst, pars = pars_era_c1), total_draws(era.sst))
    coda_rhat(get_rhat(era.sst, pars = pars_era_c1))
    coda_diag(As.mcmc.list(era.sst, pars = pars_era_c1))
dev.off()


pdf("./figures/dyn/hbm_fit/dyn_sst_diag.pdf", width = 7, height = 5)
    coda_neff(get_neff(dyn.sst, pars = pars_dyn_c1), total_draws(dyn.sst))
    coda_rhat(get_rhat(dyn.sst, pars = pars_dyn_c1))
    coda_diag(As.mcmc.list(dyn.sst, pars = pars_dyn_c1))
dev.off()

## Comp

pdf("./figures/dyn/hbm_fit/era_comp_diag.pdf", width = 7, height = 5)
coda_neff(get_neff(era.comp, pars = pars_era_c1), total_draws(era.comp))
coda_rhat(get_rhat(era.comp, pars = pars_era_c1))
coda_diag(As.mcmc.list(era.comp, pars = pars_era_c1))
dev.off()


pdf("./figures/dyn/hbm_fit/dyn_comp_diag.pdf", width = 7, height = 5)
coda_neff(get_neff(dyn.comp, pars = pars_dyn_c1), total_draws(dyn.comp))
coda_rhat(get_rhat(dyn.comp, pars = pars_dyn_c1))
coda_diag(As.mcmc.list(dyn.comp, pars = pars_dyn_c1))
dev.off()



## Posterior predictive checks -----------------------------
plot_post_pc(era.sst, stan.dat.sst$y, pdf.path = "./figures/dyn/hbm_fit/era_sst_yrep.pdf")
plot_post_pc(dyn.sst, stan.dat.sst$y, pdf.path = "./figures/dyn/hbm_fit/dyn_sst_yrep.pdf")

plot_post_pc(era.comp, stan.dat.comp$y, pdf.path = "./figures/dyn/hbm_fit/era_comp_yrep.pdf")
plot_post_pc(dyn.comp, stan.dat.comp$y, pdf.path = "./figures/dyn/hbm_fit/dyn_comp_yrep.pdf")



## LOOIC ---------------------------------------------------

loo.era.sst <- rstan::loo(era.sst, cores = 4)
loo.dyn.sst <- rstan::loo(dyn.sst, cores = 4)


loo.era.comp <- rstan::loo(era.comp, cores = 4)
loo.dyn.comp <- rstan::loo(dyn.comp, cores = 4)



save(loo.era.sst, file = "./output/diagnostics/dyn/loo_era_sst.RData")
save(loo.dyn.sst, file = "./output/diagnostics/dyn/loo_dyn_sst.RData")

save(loo.era.comp, file = "./output/diagnostics/dyn/loo_era_comp.RData")
save(loo.dyn.comp, file = "./output/diagnostics/dyn/loo_dyn_comp.RData")



sum(pareto_k_values(loo.era.sst) > 0.7)
sum(pareto_k_values(loo.dyn.sst) > 0.7)

sum(pareto_k_values(loo.era.comp) > 0.7)
sum(pareto_k_values(loo.dyn.comp) > 0.7)



## SST

pdf("./figures/dyn/hbm_fit/era_sst_loo.pdf", width = 7, height = 5)
    plot(loo.era.sst, label_points = TRUE)
dev.off()

pdf("./figures/dyn/hbm_fit/dyn_sst_loo.pdf", width = 7, height = 5)
    plot(loo.dyn.sst, label_points = TRUE)
dev.off()

## Comp

pdf("./figures/dyn/hbm_fit/era_comp_loo.pdf", width = 7, height = 5)
plot(loo.era.comp, label_points = TRUE)
dev.off()

pdf("./figures/dyn/hbm_fit/dyn_comp_loo.pdf", width = 7, height = 5)
plot(loo.dyn.comp, label_points = TRUE)
dev.off()

