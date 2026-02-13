## Sensitivity analysis of including NPGO, PDO in models


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


# Specify output directories
sens.fig.dir <- here("sensitivity-analyses", "pdo-npgo", speciesFlag)
sens.fit.dir <- here("sensitivity-analyses", "fits", speciesFlag)

# Set up
brood <- if(speciesFlag=="sockeye"){
  read.csv("./data/sockeye/master_sockeye_brood_table.csv", header=T)} else if(speciesFlag=="pink") {
    read.csv("./data/pink/master_pink_brood_table.csv", header=T)
     } else if (speciesFlag=="chum"){
       read.csv("./data/chum/master_chum_brood_table.csv", header=T)}

brood <- brood[complete.cases(brood),]

#### Data - NPGO ####

## NPGO index as in Malick et al 2017: avg Dec-Mar *before* ocean entry
# Upload data
npgo <- read.csv("./data-downloaded/climate-data/npgo.csv", header=T)
# Average over Dec-Mar
npgo <- npgo %>% filter(Month %in% month.abb[c(1:3, 12)]) %>%
  mutate(yr_eff = if_else(Month=="Dec", Year+1, Year)) %>%  # Yr for Dec becomes Yr+1
  dplyr::summarize(npgo_index=mean(NPGO), .by=yr_eff) %>%
  rename(Year=yr_eff)
# Join with Stock.IDs
stk <- expand.grid(min(data_master$BY):(max(data_master$BY)+5), info_master$Stock.ID)
names(stk) <- c("Year", "Stock.ID")
npgo <- left_join(stk, npgo, by=c("Year"))
# Need brood tbl with DetailFlag column:
brood <- geographic.order(brood)
# Get weighted average index for Sockeye
if(speciesFlag=="sockeye"){
  npgo.dat <- clim.wgt.avg(brood, npgo, env.covar="npgo_index", type="first_year", out.covar = "npgo_index")
} else {
  npgo$BY <- npgo$Year - 1
  npgo.dat <- npgo[,names(npgo) != "Year"]
}
# Get standardized NPGO index
npgo.dat <- ddply(npgo.dat, .(Stock.ID), transform,
                  npgo_stnd = scale(npgo_index))
br.npgo <- left_join(data_master, npgo.dat, by=c("BY", "Stock.ID"))

# Combine species dataframes for all-species analyses
all_spp_br <- bind_rows(sock, pink, chum)
all_spp_info <- bind_rows(sock.info, pink.info, chum.info)


#### Data - PDO ####

## PDO index as in Malick et al 2017: avg Dec-Mar *before* ocean entry

# Upload data
pdo <- read.csv("./data-downloaded/climate-data/pdo.csv", header=T)
# Average over Dec-Mar
pdo <- pdo %>% filter(Month %in% month.abb[c(1:3, 12)]) %>%
  mutate(yr_eff = if_else(Month=="Dec", Year+1, Year)) %>%  # Yr for Dec becomes Yr+1
  dplyr::summarize(pdo_index=mean(PDO), .by=yr_eff) %>%
  rename(Year=yr_eff)

# Join with Stock.IDs
stk <- expand.grid(min(data_master$BY):(max(data_master$BY)+5), info_master$Stock.ID)
names(stk) <- c("Year", "Stock.ID")
pdo <- left_join(stk, pdo, by=c("Year"))
# Need brood tbl with DetailFlag column:
brood <- geographic.order(brood)
# Get weighted average index
if(speciesFlag=="sockeye"){
  pdo.dat <- clim.wgt.avg(brood, pdo, env.covar="pdo_index", type="first_year", out.covar = "pdo_index")
} else {
  pdo$BY <- npgo$Year - 1
  pdo.dat <- pdo[,names(pdo) != "Year"]
}
# Get standardized PDO index
pdo.dat <- ddply(pdo.dat, .(Stock.ID), transform,
                 pdo_stnd = scale(pdo_index))
br.pdo <- left_join(data_master, pdo.dat, by=c("BY", "Stock.ID"))

# Combine pdo & npgo
br.pdo.npgo <- cbind(br.pdo, br.npgo[, c("npgo_index", "npgo_stnd")])
write.csv(br.pdo.npgo, here(sens.fit.dir, 'data_covar_pdo_npgo.csv'), row.names = F)



#### -- MCMC - NPGO -- ####

## Run MCMC, stationary model
stan.dat.npgo <- stan_data_stat(data_master,
                               scale.x1 = TRUE,
                               var.x2 = "early_sst_stnd",
                               var.x3 = "np_pinks_sec_stnd",
                               var.region = "Ocean.Region2",
                               alpha.group = ifelse(speciesFlag=="sockeye", TRUE, FALSE))
stan.dat.npgo$x4 <- br.npgo$npgo_stnd

stat_npgo <- rstan::stan(file = "./stan/hbm_stat_3c_inter.stan",
                      data = stan.dat.npgo,
                      #pars = c(pars.stat, pars.gen.quant),
                      warmup = 1000,
                      iter = ifelse(speciesFlag=="pink", 3000, 2000),
                      cores = 4,
                      chains = 4,
                      seed = 123,
                      control = list(adapt_delta = ifelse(speciesFlag=="pink", 0.999, 0.99),
                                     max_treedepth = 20)) # increased treedepth
save(stat_npgo, file = here(sens.fit.dir, "stat_npgo.RData"))


## Run MCMC, Eras model
stan.dat.era.npgo <- stan_data_dyn(data_master,
                                   var.x2 = "early_sst_stnd",
                                   var.x3 = "np_pinks_sec_stnd",
                                   breakpoint1 = 1989,
                                   breakpoint2 = 2011,
                                   var.region="Ocean.Region2",
                                   scale.x1 = TRUE,
                                   alpha.group = ifelse(speciesFlag=="sockeye", TRUE, FALSE))
stan.dat.era.npgo$x4 <- br.npgo$npgo_stnd


era_npgo <- rstan::stan(file = "./stan/hbm_era_3c.stan",
                       data = stan.dat.era.npgo,
                       warmup = 1000,
                       iter = 2000,
                       cores = 4,
                       chains = 4,
                       seed = 123,
                       control = list(adapt_delta = 0.9))
save(era_npgo, file = here(sens.fit.dir, "era_npgo.RData"))


## Diagnostic plots -- NPGO
pdf(here(sens.fig.dir, "stat_npgo_diag.pdf"), width = 7, height = 5)
coda_neff(get_neff(stat_npgo, pars = pars.stat), total_draws(stat_npgo))
coda_rhat(get_rhat(stat_npgo, pars = pars.stat))
coda_diag(As.mcmc.list(stat_npgo, pars = pars.stat))
dev.off()

plot_post_pc(stat_npgo, stan.dat.npgo$y, data = data_master,
             pdf.path = here(sens.fig.dir, "stat_npgo_yrep.pdf"))

loo.stat_npgo <- rstan::loo(stat_npgo, cores = 4)
save(loo.stat_npgo, file = here(sens.fig.dir, "loo_stat_npgo.RData"))
waic.stat_npgo <- loo::waic(loo::extract_log_lik(stat_npgo, "log_lik"))
save(waic.stat_npgo, file = here(sens.fig.dir, "waic_stat_npgo.RData"))

pdf(here(sens.fig.dir, "stat_npgo_loo.pdf"), width = 7, height = 5)
plot(loo.stat_npgo, label_points = TRUE)
dev.off()

#r2.stat_npgo <- bayes_R2(data_master$lnRS, as.matrix(stat_npgo, pars = "yhat"))
#save(r2.stat_npgo, file = here(sens.fig.dir, "r2_stat_npgo.RData"))

pdf(here(sens.fig.dir, "stat_npgo_resid.pdf"), width = 8, height = 8)
plot_hbm_resids(stat_npgo, data_master)
dev.off()

## Table: coefficients ----

#NPGO
gamma <- rstan::summary(stat_npgo, pars = "mu_gamma")$summary
kappa <- rstan::summary(stat_npgo, pars = "mu_kappa")$summary
chi <- rstan::summary(stat_npgo, pars = "mu_chi")$summary
lambda <- rstan::summary(stat_npgo, pars = "mu_lambda")$summary

reg   <- unique(info_master$ocean_region_lab)

tab.g <- data.frame(reg = reg,
                    coef = "SST",
                    lower = gamma[ , "2.5%"],
                    mean = gamma[ , "mean"],
                    upper = gamma[ , "97.5%"])
tab.k <- data.frame(reg = reg,
                    coef = "Comp",
                    lower = kappa[ , "2.5%"],
                    mean = kappa[ , "mean"],
                    upper = kappa[ , "97.5%"])
tab.c <- data.frame(reg = reg,
                    coef = "SST x Comp",
                    lower = chi[ , "2.5%"],
                    mean = chi[ , "mean"],
                    upper = chi[ , "97.5%"])
tab.l <- data.frame(reg = reg,
                    coef = "NPGO",
                    lower = lambda[ , "2.5%"],
                    mean = lambda[ , "mean"],
                    upper = lambda[ , "97.5%"])

tab.coef <- rbind(tab.g, tab.k, tab.c, tab.l)
tab.coef$perc <- (exp(tab.coef$mean) - 1) * 100
row.names(tab.coef) <- NULL
names(tab.coef) <- c("Ocean Region", "Coefficient", "Lower_95_CI", "Mean",
                     "Upper_95_CI", "perc_change_RS")

write.csv(tab.coef, file = here(sens.fig.dir, paste0("model_coefficients_stat_npgo_", speciesFlag, ".csv")))






#### MCMC - PDO -- ####

## Run MCMC, stationary model
stan.dat.pdo <- stan_data_stat(data_master,
                                scale.x1 = TRUE,
                                var.x2 = "early_sst_stnd",
                                var.x3 = "np_pinks_sec_stnd",
                                var.region = "Ocean.Region2",
                                alpha.group = ifelse(speciesFlag=="sockeye", TRUE, FALSE))
stan.dat.pdo$x4 <- br.pdo$pdo_stnd

stat_pdo <- rstan::stan(file = "./stan/hbm_stat_3c_inter.stan",
                         data = stan.dat.pdo,
                         #pars = c(pars.stat, pars.gen.quant),
                         warmup = 1000,
                         iter = ifelse(speciesFlag=="pink", 3000, 2000),
                         cores = 4,
                         chains = 4,
                         seed = 123,
                         control = list(adapt_delta = ifelse(speciesFlag=="pink", 0.999, 0.99),
                                        max_treedepth = 20)) # increased treedepth
save(stat_pdo, file = here(sens.fit.dir, "stat_pdo.RData"))


## Run MCMC, Eras model

stan.dat.era.pdo <- stan_data_dyn(data_master,
                                   var.x2 = "early_sst_stnd",
                                   var.x3 = "np_pinks_sec_stnd",
                                   breakpoint1 = 1989,
                                   breakpoint2 = 2011,
                                   var.region="Ocean.Region2",
                                   scale.x1 = TRUE,
                                   alpha.group = ifelse(speciesFlag=="sockeye", TRUE, FALSE))
stan.dat.era.pdo$x4 <- br.pdo$pdo_stnd


era_pdo <- rstan::stan(file = "./stan/hbm_era_3c.stan",
                        data = stan.dat.era.pdo,
                        warmup = 1000,
                        iter = 2000,
                        cores = 4,
                        chains = 4,
                        seed = 123,
                        control = list(adapt_delta = 0.9))
save(era_pdo, file = here(sens.fit.dir, "era_pdo.RData"))


#### -- Visualize fits -- ####

## Diagnostic plots -- PDO
pdf(here(sens.fig.dir, "stat_pdo_diag.pdf"), width = 7, height = 5)
coda_neff(get_neff(stat_pdo, pars = pars.stat), total_draws(stat_pdo))
coda_rhat(get_rhat(stat_pdo, pars = pars.stat))
coda_diag(As.mcmc.list(stat_pdo, pars = pars.stat))
dev.off()

plot_post_pc(stat_pdo, stan.dat.pdo$y, data = data_master,
             pdf.path = here(sens.fig.dir, "stat_pdo_yrep.pdf"))

loo.stat_pdo <- rstan::loo(stat_pdo, cores = 4)
save(loo.stat_pdo, file = here(sens.fig.dir, "loo_stat_pdo.RData"))
waic.stat_pdo <- loo::waic(loo::extract_log_lik(stat_pdo, "log_lik"))
save(waic.stat_pdo, file = here(sens.fig.dir, "waic_stat_pdo.RData"))

pdf(here(sens.fig.dir, "stat_pdo_loo.pdf"), width = 7, height = 5)
plot(loo.stat_pdo, label_points = TRUE)
dev.off()

#r2.stat_pdo <- bayes_R2(data_master$lnRS, as.matrix(stat_pdo, pars = "yhat"))
#save(r2.stat_pdo, file = here(sens.fig.dir, "r2_stat_pdo.RData"))

pdf(here(sens.fig.dir, "stat_pdo_resid.pdf"), width = 8, height = 8)
plot_hbm_resids(stat_pdo, data_master)
dev.off()


## Table: coefficients ----

# PDO
gamma <- rstan::summary(stat_pdo, pars = "mu_gamma")$summary
kappa <- rstan::summary(stat_pdo, pars = "mu_kappa")$summary
lambda <- rstan::summary(stat_pdo, pars = "mu_lambda")$summary

reg   <- unique(info_master$ocean_region_lab)

tab.g <- data.frame(reg = reg,
                    coef = "SST",
                    lower = gamma[ , "2.5%"],
                    mean = gamma[ , "mean"],
                    upper = gamma[ , "97.5%"])
tab.k <- data.frame(reg = reg,
                    coef = "Comp",
                    lower = kappa[ , "2.5%"],
                    mean = kappa[ , "mean"],
                    upper = kappa[ , "97.5%"])
tab.l <- data.frame(reg = reg,
                    coef = "PDO",
                    lower = lambda[ , "2.5%"],
                    mean = lambda[ , "mean"],
                    upper = lambda[ , "97.5%"])

tab.coef <- rbind(tab.g, tab.k, tab.l) # add tab.c if exists
tab.coef$perc <- (exp(tab.coef$mean) - 1) * 100
row.names(tab.coef) <- NULL
names(tab.coef) <- c("Ecosystem", "Coefficient", "Lower 95% CI", "Mean",
                     "Upper 95% CI", "Mean % change in R/S")

write.csv(tab.coef, file = here(sens.fig.dir, paste0("model_coefficients_stat_pdo_", speciesFlag, ".csv")))

