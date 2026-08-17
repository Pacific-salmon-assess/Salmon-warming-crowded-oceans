## Fit Nonstationary hierarchical Bayesian models ('Era' and Random walk)
# This is adapted from: https://github.com/michaelmalick/sockeye-nonstationary

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

stan.dat.2c <- stan_data_dyn(data_master,
                             var.x2 = "early_sst_stnd",
                             var.x3 = "np_pinks_sec_stnd",
                             breakpoint1 = 1989,
                             breakpoint2 = 2011,
                             var.region="Ocean.Region2",
                             scale.x1 = TRUE,
                             alpha.group = ifelse(speciesFlag=="sockeye", TRUE, FALSE))


## Set pars to monitor -------------------------------------

pars_era_2c <- c("alpha", "beta", "sigma", "phi", "mu_alpha", "sigma_alpha",
                 "gamma1", "gamma2", "gamma3",
                 "mu_gamma1", "mu_gamma2", "mu_gamma3", "sigma_gamma",
                 "kappa",                                # CHANGED: kappa is now stationary -> single parameter, no era suffix
                 "mu_kappa", "sigma_kappa" )              # CHANGED: single mu_kappa, no era suffix


pars.gen.quant <- c("log_lik", "yrep", "yhat") ## Generated quantities to monitor

save(pars_era_2c, file = "./output/pars_era_2c.RData")

save(pars.gen.quant, file = "./output/pars_gen_quant.RData")



## Run MCMC  -----------------------------------------

era.2c <- rstan::stan(file = "./stan/hbm_era_gamma_stat_kappa_2c.stan",
                       data = stan.dat.2c,
                       pars = c(pars_era_2c, pars.gen.quant),
                       warmup = 1000,
                       iter = ifelse(speciesFlag=="pink", 3000, 2000), # 2000 not enough for ESS - pinks
                       cores = 4,
                       chains = 4,
                       seed = 123,
                       control = list(adapt_delta = 0.99,
                                      max_treedepth = 20))
save(era.2c, file = here(fit.dir, "hbm_era_2c.RData"))


## Diagnostics -----------------------------------------

## Check pathology

rstan::check_hmc_diagnostics(era.2c)

neff_lowest(era.2c, pars = pars_era_2c)

rhat_highest(era.2c, pars = pars_era_2c)

# pairs_lowest(era.2c, pars = pars_era_2c)

rstan::get_elapsed_time(era.2c)

## MCMC diagnostics

pdf(here(diag.dir, "era_2c_diag.pdf"), width = 7, height = 5)
coda_neff(get_neff(era.2c, pars = pars_era_2c), total_draws(era.2c))
coda_rhat(get_rhat(era.2c, pars = pars_era_2c))
coda_diag(As.mcmc.list(era.2c, pars = pars_era_2c))
dev.off()




## Posterior predictive checks
plot_post_pc(era.2c, stan.dat.2c$y, pdf.path = here(diag.dir, "era_2c_yrep.pdf"))

## LOOIC

loo.era.2c <- rstan::loo(era.2c, cores = 4)

save(loo.era.2c, file = here(diag.dir, "loo_era_2c.RData"))

sum(pareto_k_values(loo.era.2c) > 0.7)

pdf(here(diag.dir, "era_2c_loo.pdf"), width = 7, height = 5)
plot(loo.era.2c, label_points = TRUE)
dev.off()

