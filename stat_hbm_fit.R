## Fit stationary hierarchical Bayesian models ##
## -------------------------------- ##
## A simplified version including only the model: np pinks + early sst (NO interaction)

## Set Species -----
#speciesFlag = "pink"
speciesFlag = "chum"
#speciesFlag = "sockeye"

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

## Monitor params
pars.stat <- c("alpha", "beta", "sigma", "phi", "mu_alpha", "sigma_alpha",
               "gamma", "mu_gamma", "sigma_gamma",
               "kappa", "mu_kappa", "sigma_kappa")
save(pars.stat, file = "./output/pars_stat.RData")
pars.gen.quant <- c("log_lik", "yhat", "yrep", "yresid") ## Generated quantities to monitor


## stat_a - all yrs of data

## Run MCMC
stan.dat.all <- stan_data_stat(data_master,
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
save(stat_a, file = here(fit.dir, "stat_a.RData"))

## Diagnostic plots
pdf(here(fig.dir, "stat_a_diag.pdf"), width = 7, height = 5)
    coda_neff(get_neff(stat_a, pars = pars.stat), total_draws(stat_a))
    coda_rhat(get_rhat(stat_a, pars = pars.stat))
    coda_diag(As.mcmc.list(stat_a, pars = pars.stat))
dev.off()

#plot_post_pc(stat_a, stan.dat.all$y, here(fig.dir, "stat_a_yrep.pdf")) # Bugging

loo.stat_a <- rstan::loo(stat_a, cores = 4)
save(loo.stat_a, file = here(diag.dir, "loo_stat_a.RData"))
waic.stat_a <- loo::waic(loo::extract_log_lik(stat_a, "log_lik"))
save(waic.stat_a, file = here(diag.dir, "waic_stat_a.RData"))
pdf(here(fig.dir, "stat_a_loo.pdf"), width = 7, height = 5)
    plot(loo.stat_a, label_points = TRUE)
dev.off()

r2.stat_a <- bayes_R2(data_master$lnRS, as.matrix(stat_a, pars = "yhat"))
save(r2.stat_a, file = here(diag.dir, "r2_stat_a.Rdata"))

pdf( here(fig.dir, "stat_a_resid.pdf"), width = 8, height = 8)
    plot_hbm_resids(stat_a, data_master)
dev.off()



## stat_tr ----------------------------------------------
## truncated timeseries: only brood years post 76/77, total pink North Pacific, four ocean regions


## Run MCMC
stan.dat.tr <- stan_data_stat(data_master[data_master$BY >= 1975, ],
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
save(stat_tr, file = here(fit.dir, "stat_tr.RData"))

## Diagnostic plots
pdf(here(fig.dir, "stat_tr_diag.pdf"), width = 7, height = 5)
coda_neff(get_neff(stat_tr, pars = pars.stat), total_draws(stat_tr))
coda_rhat(get_rhat(stat_tr, pars = pars.stat))
coda_diag(As.mcmc.list(stat_tr, pars = pars.stat))
dev.off()

#plot_post_pc(stat_tr, stan.dat.tr$y, data = data_master[data_master$BY >= 1975,],
         #    pdf.path = here(fig.dir, "stat_tr_yrep.pdf"))
# Bugging 

loo.stat_tr <- rstan::loo(stat_tr, cores = 4)
save(loo.stat_tr, file = here(diag.dir, "loo_stat_tr.RData"))
waic.stat_tr <- loo::waic(loo::extract_log_lik(stat_tr, "log_lik"))
save(waic.stat_tr, file = here(diag.dir, "waic_stat_tr.RData"))
pdf(here(fig.dir, "stat_tr_loo.pdf"), width = 7, height = 5)
plot(loo.stat_tr, label_points = TRUE)
dev.off()

r2.stat_tr <- bayes_R2(data_master$lnRS[data_master$BY >= 1975], as.matrix(stat_tr, pars = "yhat"))
save(r2.stat_tr, file = here(diag.dir, "r2_stat_tr.RData"))

pdf(here(fig.dir, "stat_tr_resid.pdf"), width = 8, height = 8)
plot_hbm_resids(stat_tr, data_master[data_master$BY >= 1975,])
dev.off()


## stat_ctrl ----------------------------------------------
## 'control': replicates the main analysis from Connors et al. 2020
# ONLY for sockeye.

if(speciesFlag == "Sockeye") {
    
  info_2020 <- read.csv("./data/master_stock_info_2020.csv", header=T)
  ctrl_dat <- dplyr::filter(data_master, Stock %in% info_2020$Stock)
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
  save(stat_ctrl, file = here(fit.dir, "stat_ctrl.RData"))
  
  ## Diagnostic plots 
  pdf(here(fig.dir, "stat_ctrl_diag.pdf"), width = 7, height = 5)
  coda_neff(get_neff(stat_ctrl, pars = pars.stat), total_draws(stat_ctrl))
  coda_rhat(get_rhat(stat_ctrl, pars = pars.stat))
  coda_diag(As.mcmc.list(stat_ctrl, pars = pars.stat))
  dev.off()
  
  #plot_post_pc(stat_ctrl, stan.dat.ctrl$y, data = data_master[data_master$BY %in% c.by & data_master$Stock %in% c.stk, ],
   #            pdf.path = here(fig.dir, "stat_ctrl_yrep.pdf"))
  # BUgging
  
  loo.stat_ctrl <- rstan::loo(stat_ctrl, cores = 4)
  save(loo.stat_ctrl, file = here(diag.dir, "loo_stat_ctrl.RData"))
  waic.stat_ctrl <- loo::waic(loo::extract_log_lik(stat_ctrl, "log_lik"))
  save(waic.stat_ctrl, file = here(diag.dir, "waic_stat_ctrl.RData"))
  pdf(here(fig.dir, "stat_ctrl_loo.pdf"), width = 7, height = 5)
  plot(loo.stat_ctrl, label_points = TRUE)
  dev.off()
  
  # r2.stat_ctrl <- bayes_R2(data_master$lnRS[data_master$BY %in% c.by & data_master$Stock %in% c.stk, ], as.matrix(stat_ctrl, pars = "yhat")) # BC: not working, incorrect dimensions error
  
  #save(r2.stat_ctrl, file = here(diag.dir, "r2_stat_ctrl.RData"))
  
  pdf(here(fig.dir, "stat_ctrl_resid.pdf"), width = 8, height = 8)
  plot_hbm_resids(stat_ctrl, ctrl_dat)
  dev.off()

} # end if(sockeye)


## Check pathology ----------------------------------------- 
rstan::check_hmc_diagnostics(stat_a) # not outputting anything?
rstan::check_hmc_diagnostics(stat_tr)
if(speciesFlag=="sockeye") rstan::check_hmc_diagnostics(stat_ctrl)

rstan::get_elapsed_time(stat_a)
rstan::get_elapsed_time(stat_tr)
if(speciesFlag=="sockeye") rstan::get_elapsed_time(stat_ctrl)

summary(stat_a, pars = pars.stat)$summary
summary(stat_tr, pars = pars.stat)$summary
if(speciesFlag=="sockeye") summary(stat_ctrl, pars = pars.stat)$summary

neff_lowest(stat_a, pars = pars.stat)
neff_lowest(stat_tr, pars = pars.stat)
if(speciesFlag=="sockeye") neff_lowest(stat_ctrl, pars = pars.stat)

rhat_highest(stat_a, pars = pars.stat)
rhat_highest(stat_tr, pars = pars.stat)
if(speciesFlag=="sockeye") rhat_highest(stat_ctrl, pars = pars.stat)

pairs_lowest(stat_a, pars = pars.stat) # can ignore 'warning: not a graphical parameter'
pairs_lowest(stat_tr, pars = pars.stat)
if(speciesFlag=="sockeye") pairs_lowest(stat_ctrl, pars = pars.stat)

