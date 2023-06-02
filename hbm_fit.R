## (7) Fit hierarchical Bayesian models ##
## -------------------------------- ##

## A simplified version including only the "hb05" model: np pinks + early sst (NO interaction)


if(!dir.exists("./figures/hbm_fit/"))
    dir.create("./figures/hbm_fit/")
if(!dir.exists("./output/models/"))
    dir.create("./output/models/")

pars.gen.quant <- c("log_lik", "yhat", "yrep", "yresid") ## Generated quantities to monitor


## hb05a.pr1 
## Total pink North Pacific

## Monitor params
pars.hb05 <- c("alpha", "beta", "sigma", "phi", "mu_alpha", "sigma_alpha",
               "gamma", "mu_gamma", "sigma_gamma",
               "kappa", "mu_kappa", "sigma_kappa")
save(pars.hb05, file = "./output/pars_hb05.RData")


## Run MCMC
stan.dat.hb05a <- stan_data(sock,
                            scale.x1 = TRUE,
                            var.x2 = "early_sst_stnd",
                            var.x3 = "np_pinks_sec_stnd")
hb05a <- rstan::stan(file = "./stan/hb05_pr1.stan",
                     data = stan.dat.hb05a,
                     pars = c(pars.hb05, pars.gen.quant),
                     warmup = 1000,
                     iter = 4000,
                     cores = 4,
                     chains = 4,
                     thin = 2,
                     seed = 123,
                     control = list(adapt_delta = 0.90,
                                    max_treedepth = 10))
save(hb05a, file = "./output/models/hb05a.RData")

pdf("./figures/hbm_fit/hb05a_diag.pdf", width = 7, height = 5)
    coda_neff(get_neff(hb05a, pars = pars.hb05), total_draws(hb05a))
    coda_rhat(get_rhat(hb05a, pars = pars.hb05))
    coda_diag(As.mcmc.list(hb05a, pars = pars.hb05))
dev.off()

plot_post_pc(hb05a, stan.dat.hb05a$y, "./figures/hbm_fit/hb05a_yrep.pdf")

loo.hb05a <- rstan::loo(hb05a, cores = 4)
save(loo.hb05a, file = "./output/loo_hb05a.RData")
waic.hb05a <- loo::waic(loo::extract_log_lik(hb05a, "log_lik"))
save(waic.hb05a, file = "./output/waic_hb05a.RData")
pdf("./figures/hbm_fit/hb05a_loo.pdf", width = 7, height = 5)
    plot(loo.hb05a, label_points = TRUE)
dev.off()

r2.hb05a <- bayes_R2(sock$lnRS, as.matrix(hb05a, pars = "yhat"))
save(r2.hb05a, file = "./output/r2_hb05a.RData")

pdf("./figures/hbm_fit/hb05a_resid.pdf", width = 8, height = 8)
    plot_hbm_resids(hb05a, sock)
dev.off()



## hb05r2.pr1 ----------------------------------------------
## regime: only brood years post 76/77

## Monitor params
pars.hb05 <- c("alpha", "beta", "sigma", "phi", "mu_alpha", "sigma_alpha",
               "gamma", "mu_gamma", "sigma_gamma",
               "kappa", "mu_kappa", "sigma_kappa"
               )
save(pars.hb05, file = "./output/pars_hb05.RData")


## Run MCMC
stan.dat.hb05r2 <- stan_data(sock[sock$BY >= 1975, ],
                             scale.x1 = TRUE,
                             var.x2 = "early_sst_stnd",
                             var.x3 = "np_pinks_sec_stnd")
hb05r2 <- rstan::stan(file = "./stan/hb05_pr1.stan",
                      data = stan.dat.hb05r2,
                      pars = c(pars.hb05, pars.gen.quant),
                      warmup = 1000,
                      iter = 4000,
                      cores = 4,
                      chains = 4,
                      thin = 2,
                      seed = 123,
                      control = list(adapt_delta = 0.90,
                                     max_treedepth = 10))
save(hb05r2, file = "./output/models/hb05r2.RData")

pdf("./figures/hbm_fit/hb05r2_diag.pdf", width = 7, height = 5)
coda_neff(get_neff(hb05r2, pars = pars.hb05), total_draws(hb05r2))
coda_rhat(get_rhat(hb05r2, pars = pars.hb05))
coda_diag(As.mcmc.list(hb05r2, pars = pars.hb05))
dev.off()

plot_post_pc(hb05r2, stan.dat.hb05r2$y, data = sock[sock$BY >= 1975,],
             pdf.path = "./figures/hbm_fit/hb05r2_yrep.pdf")

loo.hb05r2 <- rstan::loo(hb05r2, cores = 4)
save(loo.hb05r2, file = "./output/loo_hb05r2.RData")
waic.hb05r2 <- loo::waic(loo::extract_log_lik(hb05r2, "log_lik"))
save(waic.hb05r2, file = "./output/waic_hb05r2.RData")
pdf("./figures/hbm_fit/hb05r2_loo.pdf", width = 7, height = 5)
plot(loo.hb05r2, label_points = TRUE)
dev.off()

r2.hb05r2 <- bayes_R2(sock$lnRS[sock$BY >= 1975], as.matrix(hb05r2, pars = "yhat"))
save(r2.hb05r2, file = "./output/r2_hb05r2.RData")

pdf("./figures/hbm_fit/hb05r2_resid.pdf", width = 8, height = 8)
plot_hbm_resids(hb05r2, sock[sock$BY >= 1975,])
dev.off()




## Check pathology -----------------------------------------
rstan::check_hmc_diagnostics(hb05a)
rstan::check_hmc_diagnostics(hb05r2)

rstan::get_elapsed_time(hb05a)
rstan::get_elapsed_time(hb05r2)

summary(hb05a, pars = pars.hb05)
summary(hb05r2, pars = pars.hb05)

neff_lowest(hb05a, pars = pars.hb05)
neff_lowest(hb05r2, pars = pars.hb05)

rhat_highest(hb05a, pars = pars.hb05)
rhat_highest(hb05r2, pars = pars.hb05)

pairs_lowest(hb05a, pars = pars.hb05)
pairs_lowest(hb05r2, pars = pars.hb05)

