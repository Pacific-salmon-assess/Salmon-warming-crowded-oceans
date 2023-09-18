## (7) Fit hierarchical Bayesian models ##
## -------------------------------- ##

## A simplified version including only the model: np pinks + early sst (NO interaction)


if(!dir.exists("./figures/stat/hbm_fit/"))
    dir.create("./figures/stat/hbm_fit/")
if(!dir.exists("./output/models/stat/"))
    dir.create("./output/models/stat/")
if(!dir.exists("./output/diagnostics/"))
  dir.create("./output/diagnostics/")

## Monitor params
pars.stat <- c("alpha", "beta", "sigma", "phi", "mu_alpha", "sigma_alpha",
               "gamma", "mu_gamma", "sigma_gamma",
               "kappa", "mu_kappa", "sigma_kappa")
save(pars.stat, file = "./output/pars_stat.RData")
pars.gen.quant <- c("log_lik", "yhat", "yrep", "yresid") ## Generated quantities to monitor


## stat_a
## Total pink North Pacific, all years, four ocean regions

## Run MCMC
stan.dat.all <- stan_data_stat(sock,
                            scale.x1 = TRUE,
                            var.x2 = "early_sst_stnd",
                            var.x3 = "np_pinks_sec_stnd",
                            var.region = "Ocean.Region2")
stat_a <- rstan::stan(file = "./stan/hbm_stat_2c.stan",
                     data = stan.dat.all,
                     pars = c(pars.stat, pars.gen.quant),
                     warmup = 1000,
                     iter = 4000,
                     cores = 4,
                     chains = 4,
                     thin = 2,
                     seed = 123,
                     control = list(adapt_delta = 0.95,
                                    max_treedepth = 10))
save(stat_a, file = "./output/models/stat/stat_a.RData")

## Diagnostic plots
pdf("./figures/stat/hbm_fit/stat_a_diag.pdf", width = 7, height = 5)
    coda_neff(get_neff(stat_a, pars = pars.stat), total_draws(stat_a))
    coda_rhat(get_rhat(stat_a, pars = pars.stat))
    coda_diag(As.mcmc.list(stat_a, pars = pars.stat))
dev.off()

#plot_post_pc(stat_a, stan.dat.all$y, "./figures/stat/hbm_fit/stat_a_yrep.pdf") # Bugging

loo.stat_a <- rstan::loo(stat_a, cores = 4)
save(loo.stat_a, file = "./output/diagnostics/loo_stat_a.RData")
waic.stat_a <- loo::waic(loo::extract_log_lik(stat_a, "log_lik"))
save(waic.stat_a, file = "./output/diagnostics/waic_stat_a.RData")
pdf("./figures/stat/hbm_fit/stat_a_loo.pdf", width = 7, height = 5)
    plot(loo.stat_a, label_points = TRUE)
dev.off()

r2.stat_a <- bayes_R2(sock$lnRS, as.matrix(stat_a, pars = "yhat"))
save(r2.stat_a, file = "./output/diagnostics/r2_stat_a.RData")

pdf("./figures/stat/hbm_fit/stat_a_resid.pdf", width = 8, height = 8)
    plot_hbm_resids(stat_a, sock)
dev.off()



## stat_tr ----------------------------------------------
## truncated timeseries: only brood years post 76/77, total pink North Pacific, four ocean regions


## Run MCMC
stan.dat.tr <- stan_data_stat(sock[sock$BY >= 1975, ],
                             scale.x1 = TRUE,
                             var.x2 = "early_sst_stnd",
                             var.x3 = "np_pinks_sec_stnd",
                             var.region = "Ocean.Region2")
stat_tr <- rstan::stan(file = "./stan/hbm_stat_2c.stan",
                      data = stan.dat.tr,
                      pars = c(pars.stat, pars.gen.quant),
                      warmup = 1000,
                      iter = 4000,
                      cores = 4,
                      chains = 4,
                      thin = 2,
                      seed = 123,
                      control = list(adapt_delta = 0.95,
                                     max_treedepth = 10))
save(stat_tr, file = "./output/models/stat/stat_tr.RData")

## Diagnostic plots
pdf("./figures/stat/hbm_fit/stat_tr_diag.pdf", width = 7, height = 5)
coda_neff(get_neff(stat_tr, pars = pars.stat), total_draws(stat_tr))
coda_rhat(get_rhat(stat_tr, pars = pars.stat))
coda_diag(As.mcmc.list(stat_tr, pars = pars.stat))
dev.off()

#plot_post_pc(stat_tr, stan.dat.tr$y, data = sock[sock$BY >= 1975,],
         #    pdf.path = "./figures/stat/hbm_fit/stat_tr_yrep.pdf")
# Bugging 

loo.stat_tr <- rstan::loo(stat_tr, cores = 4)
save(loo.stat_tr, file = "./output/diagnostics/loo_stat_tr.RData")
waic.stat_tr <- loo::waic(loo::extract_log_lik(stat_tr, "log_lik"))
save(waic.stat_tr, file = "./output/diagnostics/waic_stat_tr.RData")
pdf("./figures/stat/hbm_fit/stat_tr_loo.pdf", width = 7, height = 5)
plot(loo.stat_tr, label_points = TRUE)
dev.off()

r2.stat_tr <- bayes_R2(sock$lnRS[sock$BY >= 1975], as.matrix(stat_tr, pars = "yhat"))
save(r2.stat_tr, file = "./output/diagnostics/r2_stat_tr.RData")

pdf("./figures/stat/hbm_fit/stat_tr_resid.pdf", width = 8, height = 8)
plot_hbm_resids(stat_tr, sock[sock$BY >= 1975,])
dev.off()


## stat_ctrl ----------------------------------------------
## 'control': only brood years post 76/77, total pink North Pacific, THREE ocean regions, only years & stocks in 2020 analysis 


info_2020 <- read.csv("./data/master_stock_info_2020.csv", header=T)
ctrl_dat <- dplyr::filter(sock, Stock %in% info_2020$Stock)
ctrl_dat <- ddply(ctrl_dat, .(Stock), function(x) {
     df <- dplyr::filter(x, 
                         BY >= info_2020$yr_start[info_2020$Stock == unique(x$Stock)],
                         BY <= info_2020$yr_end[info_2020$Stock == unique(x$Stock)])
            return(df)
            } )

## Run MCMC
stan.dat.ctrl <- stan_data_stat(ctrl_dat,
                            scale.x1 = TRUE,
                            var.x2 = "early_sst_stnd",
                            var.x3 = "np_pinks_sec_stnd",
                            var.region = "Ocean.Region")
stat_ctrl <- rstan::stan(file = "./stan/hbm_stat_2c.stan",
                     data = stan.dat.ctrl,
                     pars = c(pars.stat, pars.gen.quant),
                     warmup = 1000,
                     iter = 4000,
                     cores = 4,
                     chains = 4,
                     thin = 2,
                     seed = 123,
                     control = list(adapt_delta = 0.95,
                                    max_treedepth = 10))
save(stat_ctrl, file = "./output/models/stat/stat_ctrl.RData")

## Diagnostic plots 
pdf("./figures/stat/hbm_fit/stat_ctrl_diag.pdf", width = 7, height = 5)
coda_neff(get_neff(stat_ctrl, pars = pars.stat), total_draws(stat_ctrl))
coda_rhat(get_rhat(stat_ctrl, pars = pars.stat))
coda_diag(As.mcmc.list(stat_ctrl, pars = pars.stat))
dev.off()

#plot_post_pc(stat_ctrl, stan.dat.ctrl$y, data = sock[sock$BY %in% c.by & sock$Stock %in% c.stk, ],
 #            pdf.path = "./figures/stat/hbm_fit/stat_ctrl_yrep.pdf")
# BUgging

loo.stat_ctrl <- rstan::loo(stat_ctrl, cores = 4)
save(loo.stat_ctrl, file = "./output/diagnostics/loo_stat_ctrl.RData")
waic.stat_ctrl <- loo::waic(loo::extract_log_lik(stat_ctrl, "log_lik"))
save(waic.stat_ctrl, file = "./output/diagnostics/waic_stat_ctrl.RData")
pdf("./figures/stat/hbm_fit/stat_ctrl_loo.pdf", width = 7, height = 5)
plot(loo.stat_ctrl, label_points = TRUE)
dev.off()

# r2.stat_ctrl <- bayes_R2(sock$lnRS[sock$BY %in% c.by & sock$Stock %in% c.stk, ], as.matrix(stat_ctrl, pars = "yhat")) # BC: not working, incorrect dimensions error

#save(r2.stat_ctrl, file = "./output/diagnostics/r2_stat_ctrl.RData")

pdf("./figures/stat/hbm_fit/stat_ctrl_resid.pdf", width = 8, height = 8)
plot_hbm_resids(stat_ctrl, ctrl_dat)
dev.off()



## Check pathology ----------------------------------------- 
rstan::check_hmc_diagnostics(stat_a)
rstan::check_hmc_diagnostics(stat_tr)
rstan::check_hmc_diagnostics(stat_ctrl)

rstan::get_elapsed_time(stat_a)
rstan::get_elapsed_time(stat_tr)
rstan::get_elapsed_time(stat_ctrl)

summary(stat_a, pars = pars.stat)$summary
summary(stat_tr, pars = pars.stat)$summary
summary(stat_ctrl, pars = pars.stat)$summary

neff_lowest(stat_a, pars = pars.stat)
neff_lowest(stat_tr, pars = pars.stat)
neff_lowest(stat_ctrl, pars = pars.stat)

rhat_highest(stat_a, pars = pars.stat)
rhat_highest(stat_tr, pars = pars.stat)
rhat_highest(stat_ctrl, pars = pars.stat)

pairs_lowest(stat_a, pars = pars.stat) # can ignore 'warning: not a graphical parameter'
pairs_lowest(stat_tr, pars = pars.stat)
pairs_lowest(stat_ctrl, pars = pars.stat)

