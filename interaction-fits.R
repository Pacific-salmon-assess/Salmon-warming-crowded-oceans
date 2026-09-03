## Script that investigates including an interaction term in the model ##


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
fig.dir <- here("figures", "stat", speciesFlag, "hbm_inf") # place to store figures generated in this script
diag.fig.dir <- here("figures", "stat", speciesFlag, "hbm_fit") # place to store figures generated in this script
fit.dir <- here("output", "models", "stat", speciesFlag) # place to store model fits
diag.dir <- here("output", "diagnostics", "stat", speciesFlag) # place to store diagnostics

# Make them if they don't exist
if(!dir.exists(fig.dir))
  dir.create(fig.dir, recursive = TRUE)
if(!dir.exists(diag.fig.dir))
  dir.create(diag.fig.dir, recursive = TRUE)
if(!dir.exists(fit.dir))
  dir.create(fit.dir, recursive = TRUE)
if(!dir.exists(diag.dir))
  dir.create(diag.dir, recursive = TRUE)

## Monitor params
pars.stat <- c("alpha", "beta", "sigma", "phi", "mu_alpha", "sigma_alpha",
               "gamma", "mu_gamma", "sigma_gamma",
               "kappa", "mu_kappa", "sigma_kappa",
               "chi", "mu_chi", "sigma_chi")
save(pars.stat, file = "./output/pars_inter_stat.RData")
pars.gen.quant <- c("log_lik", "yhat", "yrep", "yresid") ## Generated quantities to monitor


## Run MCMC for trended competitor covariate
stan.dat.inter <- stan_data_stat(data_master,
                               scale.x1 = TRUE,
                               var.x2 = "early_sst_stnd",
                               var.x3 = "np_all_spp_sec_stnd", # comp = all spp abundance
                               var.region = "Ocean.Region2",
                               alpha.group = ifelse(speciesFlag=="sockeye", TRUE, FALSE)) # set to TRUE for sockeye
stat_inter <- rstan::stan(file = "./stan/hbm_stat_inter.stan",
                      data = stan.dat.inter,
                      pars = c(pars.stat, pars.gen.quant),
                      warmup = 1000,
                      iter = 3000,
                      cores = 4,
                      chains = 4,
                      seed = 123,
                      control = list(adapt_delta = 0.99,
                                     max_treedepth = 20)) # increased treedepth
save(stat_inter, file = here(fit.dir, "stat_inter.RData"))


## Diagnostic plots

pdf(here(diag.fig.dir, "stat_inter_diag.pdf"), width = 7, height = 5)
coda_neff(get_neff(stat_inter, pars = pars.stat), total_draws(stat_inter))
coda_rhat(get_rhat(stat_inter, pars = pars.stat))
coda_diag(As.mcmc.list(stat_inter, pars = pars.stat))
dev.off()

plot_post_pc(stat_inter, stan.dat.inter$y, pdf.path = here(diag.fig.dir, "stat_inter_yrep.pdf")) # Working again?

loo.stat_inter <- rstan::loo(stat_inter, cores = 4)
save(loo.stat_inter, file = here(diag.dir, "loo_stat_inter.RData"))
waic.stat_inter <- loo::waic(loo::extract_log_lik(stat_inter, "log_lik"))
save(waic.stat_inter, file = here(diag.dir, "waic_stat_inter.RData"))
pdf(here(diag.fig.dir, "stat_inter_loo.pdf"), width = 7, height = 5)
plot(loo.stat_inter, label_points = TRUE)
dev.off()

r2.stat_inter <- bayes_R2(data_master$lnRS, as.matrix(stat_inter, pars = "yhat"))
save(r2.stat_inter, file = here(diag.dir, "r2_stat_inter.Rdata"))

pdf( here(diag.fig.dir, "stat_inter_resid.pdf"), width = 8, height = 8)
plot_hbm_resids(stat_inter, data_master)
dev.off()




## ----------------------------- Inference and plots


# Coefficient table -- Regional (hyperparameters)

gamma <- rstan::summary(stat_inter, pars = "mu_gamma")$summary
kappa <- rstan::summary(stat_inter, pars = "mu_kappa")$summary
chi <- rstan::summary(stat_inter, pars = "mu_chi")$summary

reg <- c("West Coast", "Southeast Alaska", "Gulf of Alaska", "Bering Sea")
if(speciesFlag=="chum") reg <- reg[-2]

tab.g <- data.frame(species = speciesFlag,
                    reg = reg,
                    coef = "SST",
                    lower = gamma[ , "2.5%"],
                    mean = gamma[ , "mean"],
                    upper = gamma[ , "97.5%"])
tab.k <- data.frame(species = speciesFlag,
                    reg = reg,
                    coef = "Comp",
                    lower = kappa[ , "2.5%"],
                    mean = kappa[ , "mean"],
                    upper = kappa[ , "97.5%"])
tab.c <- data.frame(species = speciesFlag,
                    reg = reg,
                    coef = "SST x Comp",
                    lower = chi[ , "2.5%"],
                    mean = chi[ , "mean"],
                    upper = chi[ , "97.5%"])

tab.coef <- rbind(tab.g, tab.k, tab.c)
tab.coef$perc <- (exp(tab.coef$mean) - 1) * 100
row.names(tab.coef) <- NULL
names(tab.coef) <- c("Species", "Ecosystem", "Coefficient", "Lower 95% CI", "Mean",
                     "Upper 95% CI", "Mean % change in R/S")

write.csv(tab.coef, file = here(fig.dir, paste0("reg_coefficients_", "stat_inter_", speciesFlag, ".csv")))


# Coefficient table -- Population lvl
gamma <- rstan::summary(stat_inter, pars = "gamma")$summary
kappa <- rstan::summary(stat_inter, pars = "kappa")$summary
chi <- rstan::summary(stat_inter, pars = "chi")$summary

tab.g <- data.frame(species = speciesFlag,
                    stock = info_master$Stock,
                    reg = info_master$ocean_region_lab,
                    coef = "SST",
                    lower = gamma[ , "2.5%"],
                    mean = gamma[ , "mean"],
                    upper = gamma[ , "97.5%"])
tab.k <- data.frame(species = speciesFlag,
                    stock = info_master$Stock,
                    reg = info_master$ocean_region_lab,
                    coef = "Comp",
                    lower = kappa[ , "2.5%"],
                    mean = kappa[ , "mean"],
                    upper = kappa[ , "97.5%"])
tab.c <- data.frame(species = speciesFlag,
                    stock = info_master$Stock,
                    reg = info_master$ocean_region_lab,
                    coef = "SST x Comp",
                    lower = chi[ , "2.5%"],
                    mean = chi[ , "mean"],
                    upper = chi[ , "97.5%"])

tab.coef <- rbind(tab.g, tab.k, tab.c)
tab.coef$perc <- (exp(tab.coef$mean) - 1) * 100
row.names(tab.coef) <- NULL
names(tab.coef) <- c("Species", "Stock", "Region", "Coefficient", "Lower 95% CI", "Mean",
                     "Upper 95% CI", "Mean % change in R/S")

write.csv(tab.coef, file = here(fig.dir, paste0("stk_coefficients_", "stat_inter_", speciesFlag, ".csv")))



# Density plot

lst <- hb07_density_df(stat_inter, ocean.regions = ifelse(speciesFlag=="chum", 3, 4))
s.df <- lst$stock
m.df <- lst$region
m.df$region <- factor(m.df$region, levels = c("West Coast", "Gulf of Alaska", "Southeast Alaska", "Bering Sea"))

## Covariate labels
vars <- data.frame(var = levels(m.df$var))
vars$lab <- paste0("(", letters[1:nrow(vars)], ") ", vars$var)
vars$var <- factor(vars$var, levels = c("SST", "Comp", "SST x Comp", "SST + Comp"))

g <- ggplot(m.df) +
  geom_vline(xintercept = 0, color = "grey50", linetype = 2, linewidth = 0.25) +
  geom_path(data = s.df[s.df$region == "West Coast", ],
            aes(x = x, y = y, group = stock), color = col.region[["West Coast"]], alpha=0.3,
            na.rm = TRUE) +
  geom_path(data = s.df[s.df$region == "Gulf of Alaska", ],
            aes(x = x, y = y, group = stock), color = col.region[["Gulf of Alaska"]], alpha=0.3,
            na.rm = TRUE) +
  geom_path(data = s.df[s.df$region == "Bering Sea", ],
            aes(x = x, y = y, group = stock), color = col.region[["Bering Sea"]], alpha=0.3,
            na.rm = TRUE) +
  geom_path(aes(x = x, y = y, color = region), linewidth = 1, alpha=1,
            na.rm = TRUE) +
  scale_colour_manual(values=col.region) +
  labs(x = "Percent change in R/S",
       y = "Posterior density",
       color = "") +
  scale_x_continuous(limits = c(-50, 50), expand = c(0, 0)) +
  scale_y_continuous(breaks=NULL) +
  geom_text(data = vars,
            aes(x = -48.1,
                y = max(m.df$y) - 0.008,
                label = lab),
            hjust = 0,
            size = 2.7,
            color = "grey30") +
  facet_wrap( ~ var, ncol = 1) +
  theme_sleek(base_size = 9) +
  theme(legend.justification = c(0, 0),
        legend.position = c(0.7, 0.91),
        legend.key.size = unit(10, "pt"),
        legend.background = element_blank(),
        legend.text = element_text(size = 8),
        panel.spacing.y = unit(-0.5, "pt"),
        strip.background = element_blank(),
        strip.text.x = element_blank())
if(speciesFlag != "chum"){
  g <- g + geom_path(data = s.df[s.df$region == "Southeast Alaska", ],
                     aes(x = x, y = y, group = stock), color = col.region[["Southeast Alaska"]], alpha=0.3, na.rm = TRUE) }

pdf(here(fig.dir, "dens_stat_inter.pdf"), width = 4, height = 6)
print(g)
dev.off()

## Run MCMC for trended competitor covariate
stan.dat.inter.detp <- stan_data_stat(data_master,
                                      scale.x1 = TRUE,
                                      var.x2 = "early_sst_stnd",
                                      var.x3 = "det_np_all_spp_sec_stnd", # comp = all spp abundance detrended
                                      var.region = "Ocean.Region2",
                                      alpha.group = ifelse(speciesFlag=="sockeye", TRUE, FALSE)) # set to TRUE for sockeye
stat_inter.detp <- rstan::stan(file = "./stan/hbm_stat_inter.stan",
                               data = stan.dat.inter.detp,
                               pars = c(pars.stat, pars.gen.quant),
                               warmup = 1000,
                               iter = 3000,
                               cores = 4,
                               chains = 4,
                               seed = 123,
                               control = list(adapt_delta = 0.99,
                                              max_treedepth = 20)) # increased treedepth
save(stat_inter.detp, file = here(fit.dir, "stat_inter_detr.RData"))
