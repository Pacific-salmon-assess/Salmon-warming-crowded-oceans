## (7) Fit hierarchical Bayesian models ##
## -------------------------------- ##

## A simplified version including only the hb05 model: np pinks + early sst (NO interaction)


if(!dir.exists("./figures/hbm_fit/"))
    dir.create("./figures/hbm_fit/")
if(!dir.exists("./output/models/"))
    dir.create("./output/models/")

pars.gen.quant <- c("log_lik", "yhat", "yrep", "yresid") ## Generated quantities to monitor


## hb05a.pr1 
## Total pink North Pacific, all years, three ocean regions

## Monitor params
pars.hb05 <- c("alpha", "beta", "sigma", "phi", "mu_alpha", "sigma_alpha",
               "gamma", "mu_gamma", "sigma_gamma",
               "kappa", "mu_kappa", "sigma_kappa")
save(pars.hb05, file = "./output/pars_hb05.RData")


## Run MCMC
stan.dat.hb05a <- stan_data_stat(sock,
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

## Diagnostic plots
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
## regime: only brood years post 76/77, total pink North Pacific, three ocean regions

## Monitor params
pars.hb05 <- c("alpha", "beta", "sigma", "phi", "mu_alpha", "sigma_alpha",
               "gamma", "mu_gamma", "sigma_gamma",
               "kappa", "mu_kappa", "sigma_kappa"
               )
save(pars.hb05, file = "./output/pars_hb05.RData")


## Run MCMC
stan.dat.hb05r2 <- stan_data_stat(sock[sock$BY >= 1975, ],
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

## Diagnostic plots
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


## hb05c.pr1 ----------------------------------------------
## 'control': only brood years post 76/77, total pink North Pacific, three ocean regions, only years & stocks in 2020 analysis 

## Monitor params
pars.hb05 <- c("alpha", "beta", "sigma", "phi", "mu_alpha", "sigma_alpha",
               "gamma", "mu_gamma", "sigma_gamma",
               "kappa", "mu_kappa", "sigma_kappa")
save(pars.hb05, file = "./output/pars_hb05.RData")

dat_2020 <- read.csv("./data/master_brood_table_2020.csv", header=T)
c.by <- unique(dat_2020$BY)
c.stk <- unique(dat_2020$Stock)

## Run MCMC
stan.dat.hb05c <- stan_data_stat(sock[sock$BY %in% c.by & sock$Stock %in% c.stk, ],
                            scale.x1 = TRUE,
                            var.x2 = "early_sst_stnd",
                            var.x3 = "np_pinks_sec_stnd")
hb05c <- rstan::stan(file = "./stan/hb05_pr1.stan",
                     data = stan.dat.hb05c,
                     pars = c(pars.hb05, pars.gen.quant),
                     warmup = 1000,
                     iter = 4000,
                     cores = 4,
                     chains = 4,
                     thin = 2,
                     seed = 123,
                     control = list(adapt_delta = 0.90,
                                    max_treedepth = 10))
save(hb05c, file = "./output/models/hb05c.RData")

## Diagnostic plots 
pdf("./figures/hbm_fit/hb05c_diag.pdf", width = 7, height = 5)
coda_neff(get_neff(hb05c, pars = pars.hb05), total_draws(hb05c))
coda_rhat(get_rhat(hb05c, pars = pars.hb05))
coda_diag(As.mcmc.list(hb05c, pars = pars.hb05))
dev.off()

plot_post_pc(hb05c, stan.dat.hb05c$y, data = sock[sock$BY %in% c.by & sock$Stock %in% c.stk, ],
             pdf.path = "./figures/hbm_fit/hb05c_yrep.pdf")

loo.hb05c <- rstan::loo(hb05c, cores = 4)
save(loo.hb05c, file = "./output/loo_hb05c.RData")
waic.hb05c <- loo::waic(loo::extract_log_lik(hb05c, "log_lik"))
save(waic.hb05c, file = "./output/waic_hb05c.RData")
pdf("./figures/hbm_fit/hb05c_loo.pdf", width = 7, height = 5)
plot(loo.hb05c, label_points = TRUE)
dev.off()

r2.hb05c <- bayes_R2(sock$lnRS[sock$BY %in% c.by & sock$Stock %in% c.stk, ], as.matrix(hb05r2, pars = "yhat"))
save(r2.hb05c, file = "./output/r2_hb05c.RData")

pdf("./figures/hbm_fit/hb05c_resid.pdf", width = 8, height = 8)
plot_hbm_resids(hb05c, sock[sock$BY %in% c.by & sock$Stock %in% c.stk, ])
dev.off()


## hb05oc.pr1 ----------------------------------------------
## Total pink North Pacific, all years, four ocean regions

pars.hb05 <- c("alpha", "beta", "sigma", "phi", "mu_alpha", "sigma_alpha",
               "gamma", "mu_gamma", "sigma_gamma",
               "kappa", "mu_kappa", "sigma_kappa")
save(pars.hb05, file = "./output/pars_hb05.RData")

## Run MCMC
stan.dat.hb05oc <- stan_data_stat(sock,
                             scale.x1 = TRUE,
                             var.region = "Ocean.Region2",
                             var.x2 = "early_sst_stnd",
                             var.x3 = "np_pinks_sec_stnd")
hb05oc <- rstan::stan(file = "./stan/hb05_pr1.stan",
                      data = stan.dat.hb05oc,
                      pars = c(pars.hb05, pars.gen.quant),
                      warmup = 1000,
                      iter = 4000,
                      cores = 4,
                      chains = 4,
                      thin = 2,
                      seed = 123,
                      control = list(adapt_delta = 0.90,
                                     max_treedepth = 10))
save(hb05oc, file = "./output/models/hb05oc.RData")

## Diagnostic plots
pdf("./figures/hbm_fit/hb05oc_diag.pdf", width = 7, height = 5)
coda_neff(get_neff(hb05oc, pars = pars.hb05), total_draws(hb05oc))
coda_rhat(get_rhat(hb05oc, pars = pars.hb05))
coda_diag(As.mcmc.list(hb05oc, pars = pars.hb05))
dev.off()

plot_post_pc(hb05oc, stan.dat.hb05oc$y, data = sock[sock$BY >= 1975,],
             pdf.path = "./figures/hbm_fit/hb05oc_yrep.pdf")

loo.hb05oc <- rstan::loo(hb05oc, cores = 4)
save(loo.hb05oc, file = "./output/loo_hb05oc.RData")
waic.hb05oc <- loo::waic(loo::extract_log_lik(hb05oc, "log_lik"))
save(waic.hb05oc, file = "./output/waic_hb05oc.RData")
pdf("./figures/hbm_fit/hb05oc_loo.pdf", width = 7, height = 5)
plot(loo.hb05oc, label_points = TRUE)
dev.off()

r2.hb05oc <- bayes_R2(sock$lnRS[sock$BY >= 1975], as.matrix(hb05oc, pars = "yhat"))
save(r2.hb05oc, file = "./output/r2_hb05oc.RData")

pdf("./figures/hbm_fit/hb05oc_resid.pdf", width = 8, height = 8)
plot_hbm_resids(hb05oc, sock[sock$BY >= 1975,])
dev.off()

## hb05ocr2.pr1 ----------------------------------------------
## Total pink North Pacific, only brood years post 76/77, four ocean regions

## Monitor params
pars.hb05 <- c("alpha", "beta", "sigma", "phi", "mu_alpha", "sigma_alpha",
               "gamma", "mu_gamma", "sigma_gamma",
               "kappa", "mu_kappa", "sigma_kappa")
save(pars.hb05, file = "./output/pars_hb05.RData")


## Run MCMC
stan.dat.hb05ocr2 <- stan_data_stat(sock[sock$BY >= 1975, ],
                             scale.x1 = TRUE,
                             var.region = "Ocean.Region2",
                             var.x2 = "early_sst_stnd",
                             var.x3 = "np_pinks_sec_stnd")
hb05ocr2 <- rstan::stan(file = "./stan/hb05_pr1.stan",
                      data = stan.dat.hb05ocr2,
                      pars = c(pars.hb05, pars.gen.quant),
                      warmup = 1000,
                      iter = 4000,
                      cores = 4,
                      chains = 4,
                      thin = 2,
                      seed = 123,
                      control = list(adapt_delta = 0.90,
                                     max_treedepth = 10))
save(hb05ocr2, file = "./output/models/hb05ocr2.RData")

## Diagnostic plots
pdf("./figures/hbm_fit/hb05ocr2_diag.pdf", width = 7, height = 5)
coda_neff(get_neff(hb05ocr2, pars = pars.hb05), total_draws(hb05ocr2))
coda_rhat(get_rhat(hb05ocr2, pars = pars.hb05))
coda_diag(As.mcmc.list(hb05ocr2, pars = pars.hb05))
dev.off()

plot_post_pc(hb05ocr2, stan.dat.hb05ocr2$y, data = sock[sock$BY >= 1975,],
             pdf.path = "./figures/hbm_fit/hb05ocr2_yrep.pdf")

loo.hb05ocr2 <- rstan::loo(hb05ocr2, cores = 4)
save(loo.hb05ocr2, file = "./output/loo_hb05ocr2.RData")
waic.hb05ocr2 <- loo::waic(loo::extract_log_lik(hb05ocr2, "log_lik"))
save(waic.hb05ocr2, file = "./output/waic_hb05ocr2.RData")
pdf("./figures/hbm_fit/hb05ocr2_loo.pdf", width = 7, height = 5)
plot(loo.hb05ocr2, label_points = TRUE)
dev.off()

r2.hb05ocr2 <- bayes_R2(sock$lnRS[sock$BY >= 1975], as.matrix(hb05ocr2, pars = "yhat"))
save(r2.hb05ocr2, file = "./output/r2_hb05ocr2.RData")

pdf("./figures/hbm_fit/hb05ocr2_resid.pdf", width = 8, height = 8)
plot_hbm_resids(hb05ocr2, sock[sock$BY >= 1975,])
dev.off()


## Check pathology ----------------------------------------- 
rstan::check_hmc_diagnostics(hb05a)
rstan::check_hmc_diagnostics(hb05r2)
rstan::check_hmc_diagnostics(hb05c)
rstan::check_hmc_diagnostics(hb05oc)
rstan::check_hmc_diagnostics(hb05ocr2)

rstan::get_elapsed_time(hb05a)
rstan::get_elapsed_time(hb05r2)
rstan::get_elapsed_time(hb05c)
rstan::get_elapsed_time(hb05oc)
rstan::get_elapsed_time(hb05ocr2)

summary(hb05a, pars = pars.hb05)
summary(hb05r2, pars = pars.hb05)
summary(hb05c, pars = pars.hb05)
summary(hb05oc, pars = pars.hb05)
summary(hb05ocr2, pars = pars.hb05)

neff_lowest(hb05a, pars = pars.hb05)
neff_lowest(hb05r2, pars = pars.hb05)
neff_lowest(hb05c, pars = pars.hb05)
neff_lowest(hb05oc, pars = pars.hb05)
neff_lowest(hb05ocr2, pars = pars.hb05)

rhat_highest(hb05a, pars = pars.hb05)
rhat_highest(hb05r2, pars = pars.hb05)
rhat_highest(hb05c, pars = pars.hb05)
rhat_highest(hb05oc, pars = pars.hb05)
rhat_highest(hb05ocr2, pars = pars.hb05)

pairs_lowest(hb05a, pars = pars.hb05)
pairs_lowest(hb05r2, pars = pars.hb05)
pairs_lowest(hb05c, pars = pars.hb05)
pairs_lowest(hb05oc, pars = pars.hb05)
pairs_lowest(hb05ocr2, pars = pars.hb05)

