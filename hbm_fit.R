## (7) Fit hierarchical Bayesian models ##
## -------------------------------- ##

## A simplified version including only the "hb07a" model: np pinks, early sst as covariates


if(!dir.exists("./figures/hbm_fit/"))
    dir.create("./figures/hbm_fit/")
if(!dir.exists("./output/models/"))
    dir.create("./output/models/")

pars.gen.quant <- c("log_lik", "yhat", "yrep", "yresid") ## Generated quantities to monitor


## hb07a.pr1 
## Total pink North Pacific

## Monitor params
pars.hb07 <- c("alpha", "beta", "sigma", "phi", "mu_alpha", "sigma_alpha",
               "gamma", "mu_gamma", "sigma_gamma",
               "kappa", "mu_kappa", "sigma_kappa",
               "chi", "mu_chi", "sigma_chi")
save(pars.hb07, file = "./output/pars_hb07.RData")


## Run MCMC
stan.dat.hb07a <- stan_data(sock,
                            scale.x1 = TRUE,
                            var.x2 = "early_sst_stnd",
                            var.x3 = "np_pinks_sec_stnd")
hb07a <- rstan::stan(file = "./stan/hb07_pr1.stan",
                     data = stan.dat.hb07a,
                     pars = c(pars.hb07, pars.gen.quant),
                     warmup = 1000,
                     iter = 4000,
                     cores = 4,
                     chains = 4,
                     thin = 2,
                     seed = 123,
                     control = list(adapt_delta = 0.90,
                                    max_treedepth = 10))
save(hb07a, file = "./output/models/hb07a.RData")

pdf("./figures/hbm_fit/hb07a_diag.pdf", width = 7, height = 5)
    coda_neff(get_neff(hb07a, pars = pars.hb07), total_draws(hb07a))
    coda_rhat(get_rhat(hb07a, pars = pars.hb07))
    coda_diag(As.mcmc.list(hb07a, pars = pars.hb07))
dev.off()

plot_post_pc(hb07a, stan.dat.hb07a$y, "./figures/hbm_fit/hb07a_yrep.pdf")

loo.hb07a <- rstan::loo(hb07a, cores = 4)
save(loo.hb07a, file = "./output/loo_hb07a.RData")
waic.hb07a <- loo::waic(loo::extract_log_lik(hb07a, "log_lik"))
save(waic.hb07a, file = "./output/waic_hb07a.RData")
pdf("./figures/hbm_fit/hb07a_loo.pdf", width = 7, height = 5)
    plot(loo.hb07a, label_points = TRUE)
dev.off()

r2.hb07a <- bayes_R2(sock$lnRS, as.matrix(hb07a, pars = "yhat"))
save(r2.hb07a, file = "./output/r2_hb07a.RData")

pdf("./figures/hbm_fit/hb07a_resid.pdf", width = 8, height = 8)
    plot_hbm_resids(hb07a, sock)
dev.off()



## hb07r2.pr1 ----------------------------------------------
## regime: only brood years post 76/77

## Monitor params
pars.hb07 <- c("alpha", "beta", "sigma", "phi", "mu_alpha", "sigma_alpha",
               "gamma", "mu_gamma", "sigma_gamma",
               "kappa", "mu_kappa", "sigma_kappa",
               "chi", "mu_chi", "sigma_chi")
save(pars.hb07, file = "./output/pars_hb07.RData")


## Run MCMC
stan.dat.hb07r2 <- stan_data(sock[sock$BY >= 1975, ],
                             scale.x1 = TRUE,
                             var.x2 = "early_sst_stnd",
                             var.x3 = "np_pinks_sec_stnd")
hb07r2 <- rstan::stan(file = "./stan/hb07_pr1.stan",
                      data = stan.dat.hb07r2,
                      pars = c(pars.hb07, pars.gen.quant),
                      warmup = 1000,
                      iter = 4000,
                      cores = 4,
                      chains = 4,
                      thin = 2,
                      seed = 123,
                      control = list(adapt_delta = 0.90,
                                     max_treedepth = 10))
save(hb07r2, file = "./output/models/hb07r2.RData")

pdf("./figures/hbm_fit/hb07r2_diag.pdf", width = 7, height = 5)
coda_neff(get_neff(hb07r2, pars = pars.hb07), total_draws(hb07r2))
coda_rhat(get_rhat(hb07r2, pars = pars.hb07))
coda_diag(As.mcmc.list(hb07r2, pars = pars.hb07))
dev.off()

plot_post_pc(hb07r2, stan.dat.hb07r2$y, data = sock[sock$BY >= 1975,],
             pdf.path = "./figures/hbm_fit/hb07r2_yrep.pdf")

loo.hb07r2 <- rstan::loo(hb07r2, cores = 4)
save(loo.hb07r2, file = "./output/loo_hb07r2.RData")
waic.hb07r2 <- loo::waic(loo::extract_log_lik(hb07r2, "log_lik"))
save(waic.hb07r2, file = "./output/waic_hb07r2.RData")
pdf("./figures/hbm_fit/hb07r2_loo.pdf", width = 7, height = 5)
plot(loo.hb07r2, label_points = TRUE)
dev.off()

r2.hb07r2 <- bayes_R2(sock$lnRS[sock$BY >= 1975], as.matrix(hb07r2, pars = "yhat"))
save(r2.hb07r2, file = "./output/r2_hb07r2.RData")

pdf("./figures/hbm_fit/hb07r2_resid.pdf", width = 8, height = 8)
plot_hbm_resids(hb07r2, sock[sock$BY >= 1975,])
dev.off()




## Check pathology -----------------------------------------
rstan::check_hmc_diagnostics(hb07a)
rstan::check_hmc_diagnostics(hb07ar2)

rstan::get_elapsed_time(hb07a)
rstan::get_elapsed_time(hb07r2)

summary(hb07a, pars = pars.hb07)
summary(hb07r2, pars = pars.hb07)

neff_lowest(hb07a, pars = pars.hb07)
neff_lowest(hb07r2, pars = pars.hb07)

rhat_highest(hb07a, pars = pars.hb07)
rhat_highest(hb07r2, pars = pars.hb07)

pairs_lowest(hb07a, pars = pars.hb07)
pairs_lowest(hb07r2, pars = pars.hb07)

