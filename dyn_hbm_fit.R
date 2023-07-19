## Fit NONSTATIONARY hierarchical Bayesian models
# This is adapted from: https://github.com/michaelmalick/sockeye-nonstationary

if(!dir.exists("./figures/dyn/hbm_fit/"))
    dir.create("./figures/dyn/hbm_fit/")

if(!dir.exists("./output/models/dyn/"))
    dir.create("./output/models/dyn/")

if(!dir.exists("./output/diagnostics/dyn/"))
    dir.create("./output/diagnostics/dyn/")


## Get data for Stan ---------------------------------------
stan.dat.sst  <- stan_data_dyn(sock, "sst_anom_stnd",
                           breakpoint1 = 1977,
                           breakpoint2 = 1989,
                           scale.x1 = TRUE)
save(stan.dat.sst, file = "./output/stan_dat_sst.RData")

stan.dat.comp <- stan_data_dyn(sock, "np_pinks_sec_stnd",
                               breakpoint1 = 1977,
                               breakpoint2 = 1989,
                               scale.x1 = TRUE)
save(stan.dat.comp, file = "./output/stan_dat_comp.RData")


## Set pars to monitor -------------------------------------

pars.hbm3 <- c("alpha", "beta", "sigma", "phi", "mu_alpha", "sigma_alpha",
               "gamma1", "gamma2", "gamma3",
               "mu_gamma1", "mu_gamma2", "mu_gamma3", "sigma_gamma")
pars.hbm5 <- c("alpha", "beta", "sigma", "phi", "mu_alpha", "sigma_alpha",
               "gamma", "sigma_gamma", "signal_noise")

pars.gen.quant <- c("log_lik", "yrep") ## Generated quantities to monitor


save(pars.hbm3, file = "./output/pars_hbm3.RData")
save(pars.hbm5, file = "./output/pars_hbm5.RData")
save(pars.gen.quant, file = "./output/pars_gen_quant.RData")



## Run MCMC  -----------------------------------------

## SST

hbm3.sst <- rstan::stan(file = "./stan/3_hbm_ar1_x2era3.stan",
                        data = stan.dat.sst,
                        pars = c(pars.hbm3, pars.gen.quant),
                        warmup = 1000,
                        iter = 4000,
                        cores = 4,
                        chains = 4,
                        thin = 3,
                        seed = sample(1:1e6, 1),
                        control = list(adapt_delta = 0.90,
                                       max_treedepth = 10))
save(hbm3.sst, file = "./output/models/dyn/hbm3_sst.RData")

hbm5.sst <- rstan::stan(file = "./stan/5_hbm_ar1_x2nonstat_sigma_diff.stan",
                        data = stan.dat.sst,
                        pars = c(pars.hbm5, pars.gen.quant),
                        warmup = 1000,
                        iter = 4000,
                        cores = 4,
                        chains = 4,
                        thin = 3,
                        seed = sample(1:1e6, 1),
                        control = list(adapt_delta = 0.98,
                                       max_treedepth = 10))
save(hbm5.sst, file = "./output/models/dyn/hbm5_sst.RData")


## Comp

hbm3.comp <- rstan::stan(file = "./stan/3_hbm_ar1_x2era3.stan",
                        data = stan.dat.comp,
                        pars = c(pars.hbm3, pars.gen.quant),
                        warmup = 1000,
                        iter = 4000,
                        cores = 4,
                        chains = 4,
                        thin = 3,
                        seed = sample(1:1e6, 1),
                        control = list(adapt_delta = 0.90,
                                       max_treedepth = 10))
save(hbm3.comp, file = "./output/models/dyn/hbm3_comp.RData")

hbm5.comp <- rstan::stan(file = "./stan/5_hbm_ar1_x2nonstat_sigma_diff.stan",
                        data = stan.dat.comp,
                        pars = c(pars.hbm5, pars.gen.quant),
                        warmup = 1000,
                        iter = 4000,
                        cores = 4,
                        chains = 4,
                        thin = 3,
                        seed = sample(1:1e6, 1),
                        control = list(adapt_delta = 0.98,
                                       max_treedepth = 10))
save(hbm5.comp, file = "./output/models/dyn/hbm5_comp.RData")


## Check pathology -----------------------------------------

## NUTS diag

rstan::check_hmc_diagnostics(hbm3.sst)
rstan::check_hmc_diagnostics(hbm5.sst)

rstan::check_hmc_diagnostics(hbm3.comp)
rstan::check_hmc_diagnostics(hbm5.comp)

## Neff

neff_lowest(hbm3.sst, pars = pars.hbm3)
neff_lowest(hbm5.sst, pars = pars.hbm5)


neff_lowest(hbm3.comp, pars = pars.hbm3)
neff_lowest(hbm5.comp, pars = pars.hbm5)

## Rhat

rhat_highest(hbm3.sst, pars = pars.hbm3)
rhat_highest(hbm5.sst, pars = pars.hbm5)

rhat_highest(hbm3.comp, pars = pars.hbm3)
rhat_highest(hbm5.comp, pars = pars.hbm5)


## Pairs
#
# pairs_lowest(hbm3.sst, pars = pars.hbm3)
# pairs_lowest(hbm5.sst, pars = pars.hbm5)
## pairs_lowest(hbm3.comp, pars = pars.hbm3)
# pairs_lowest(hbm5.comp, pars = pars.hbm5)

## Run time

rstan::get_elapsed_time(hbm3.sst)
rstan::get_elapsed_time(hbm5.sst)

rstan::get_elapsed_time(hbm3.comp)
rstan::get_elapsed_time(hbm5.comp)


## MCMC diagnostics ----------------------------------------

## SST

pdf("./figures/dyn/hbm_fit/hbm3_sst_diag.pdf", width = 7, height = 5)
    coda_neff(get_neff(hbm3.sst, pars = pars.hbm3), total_draws(hbm3.sst))
    coda_rhat(get_rhat(hbm3.sst, pars = pars.hbm3))
    coda_diag(As.mcmc.list(hbm3.sst, pars = pars.hbm3))
dev.off()


pdf("./figures/dyn/hbm_fit/hbm5_sst_diag.pdf", width = 7, height = 5)
    coda_neff(get_neff(hbm5.sst, pars = pars.hbm5), total_draws(hbm5.sst))
    coda_rhat(get_rhat(hbm5.sst, pars = pars.hbm5))
    coda_diag(As.mcmc.list(hbm5.sst, pars = pars.hbm5))
dev.off()

## Comp

pdf("./figures/dyn/hbm_fit/hbm3_comp_diag.pdf", width = 7, height = 5)
coda_neff(get_neff(hbm3.comp, pars = pars.hbm3), total_draws(hbm3.comp))
coda_rhat(get_rhat(hbm3.comp, pars = pars.hbm3))
coda_diag(As.mcmc.list(hbm3.comp, pars = pars.hbm3))
dev.off()


pdf("./figures/dyn/hbm_fit/hbm5_comp_diag.pdf", width = 7, height = 5)
coda_neff(get_neff(hbm5.comp, pars = pars.hbm5), total_draws(hbm5.comp))
coda_rhat(get_rhat(hbm5.comp, pars = pars.hbm5))
coda_diag(As.mcmc.list(hbm5.comp, pars = pars.hbm5))
dev.off()



## Posterior predictive checks -----------------------------
plot_post_pc(hbm3.sst, stan.dat.sst$y, pdf.path = "./figures/dyn/hbm_fit/hbm3_sst_yrep.pdf")
plot_post_pc(hbm5.sst, stan.dat.sst$y, pdf.path = "./figures/dyn/hbm_fit/hbm5_sst_yrep.pdf")

plot_post_pc(hbm3.comp, stan.dat.comp$y, pdf.path = "./figures/dyn/hbm_fit/hbm3_comp_yrep.pdf")
plot_post_pc(hbm5.comp, stan.dat.comp$y, pdf.path = "./figures/dyn/hbm_fit/hbm5_comp_yrep.pdf")



## LOOIC ---------------------------------------------------

loo.hbm3.sst <- rstan::loo(hbm3.sst, cores = 4)
loo.hbm5.sst <- rstan::loo(hbm5.sst, cores = 4)


loo.hbm3.comp <- rstan::loo(hbm3.comp, cores = 4)
loo.hbm5.comp <- rstan::loo(hbm5.comp, cores = 4)



save(loo.hbm3.sst, file = "./output/diagnostics/dyn/loo_hbm3_sst.RData")
save(loo.hbm5.sst, file = "./output/diagnostics/dyn/loo_hbm5_sst.RData")

save(loo.hbm3.comp, file = "./output/diagnostics/dyn/loo_hbm3_comp.RData")
save(loo.hbm5.comp, file = "./output/diagnostics/dyn/loo_hbm5_comp.RData")



sum(pareto_k_values(loo.hbm3.sst) > 0.7)
sum(pareto_k_values(loo.hbm5.sst) > 0.7)

sum(pareto_k_values(loo.hbm3.comp) > 0.7)
sum(pareto_k_values(loo.hbm5.comp) > 0.7)



## SST

pdf("./figures/dyn/hbm_fit/hbm3_sst_loo.pdf", width = 7, height = 5)
    plot(loo.hbm3.sst, label_points = TRUE)
dev.off()

pdf("./figures/dyn/hbm_fit/hbm5_sst_loo.pdf", width = 7, height = 5)
    plot(loo.hbm5.sst, label_points = TRUE)
dev.off()

## Comp

pdf("./figures/dyn/hbm_fit/hbm3_comp_loo.pdf", width = 7, height = 5)
plot(loo.hbm3.comp, label_points = TRUE)
dev.off()

pdf("./figures/dyn/hbm_fit/hbm5_comp_loo.pdf", width = 7, height = 5)
plot(loo.hbm5.comp, label_points = TRUE)
dev.off()

