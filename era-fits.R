## Re-fitting era models w/ different eras

# Set species
#speciesFlag = "pink"
speciesFlag = "chum"
#speciesFlag = "sockeye"

# Species
if(speciesFlag=="pink") 
  data_master <- pink else if(speciesFlag=="chum") 
    data_master <- chum else if(speciesFlag=="sockeye")
      data_master <- sock

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

## ---- 2 eras : pre- and post- 88/89 --------------------------------------
## Get data for Stan
stan.dat.2c.2 <- stan_data_dyn(data_master, 
                             var.x2 = "early_sst_stnd",
                             var.x3 = "np_pinks_sec_stnd",
                             breakpoint1 = 1989,
                             breakpoint2 = NULL,
                             var.region="Ocean.Region2", 
                             scale.x1 = TRUE,
                             alpha.group = FALSE)


pars_era_2c_2 <- c("alpha", "beta", "sigma", "phi", "mu_alpha", "sigma_alpha",
                 "gamma1", "gamma2",
                 "mu_gamma1", "mu_gamma2", "sigma_gamma", 
                 "kappa1", "kappa2",
                 "mu_kappa1", "mu_kappa2", "sigma_kappa" )
pars.gen.quant <- c("log_lik", "yrep") ## Generated quantities to monitor


# Run MCMC 
era.2c.2 <- rstan::stan(file = "./stan/hbm_era_2c.stan",
                      data = stan.dat.2c.2,
                      pars = c(pars_era_2c_2, pars.gen.quant),
                      warmup = 1000,
                      iter = 2000,
                      cores = 4,
                      chains = 4,
                      seed = 123,
                      control = list(adapt_delta = 0.99,
                                     max_treedepth = 20))
save(era.2c.2, file = here(fit.dir, "hbm_era_2c_2.RData"))


### --- Data

# Density dataframe - by stock
dens.df.st.2c.2 <- era_density_df(era.2c.2, par=c("gamma", "kappa"), neras=2)
dens.df.st.2c.2 <- ocean_region_lab(dens.df.st.2c.2)
# Density dataframe (regional lvl)
dens.df.reg.2c.2 <- era_density_df(era.2c.2, par=c("gamma", "kappa"), mu=T, neras=2)
dens.df.reg.2c.2 <- ocean_region_lab(dens.df.reg.2c.2)

### --- Figures

## Density plot
g <- ggplot(dens.df.st.2c.2) + 
  geom_vline(xintercept = 0, color = "grey50", linetype = 2, linewidth = 0.25) +
  geom_path(aes(x=x, y=dens, group=stock, col=ocean_region_lab), alpha=0.2) + 
  geom_path(data=dens.df.reg.2c.3, aes(x=x, y=dens, col=ocean_region_lab), alpha=0.85, linewidth=1) +
  scale_colour_manual(values=col.region) +
  facet_grid(rows=vars(era), cols=vars(varnam)) + 
  coord_cartesian(xlim=c(-1, 1)) +
  theme_minimal() + labs(x="covariate effect", y="", col="Ocean Region") +
  theme(axis.text.y=element_blank(), 
        legend.position="bottom")
pdf(here(fig.dir, "2_era_dens.pdf"))
print(g)
dev.off()


# Time series length - show overlap
prod_dat <- ocean_region_lab(data_master)
g <- ggplot(prod_dat) + 
  geom_vline(xintercept=c(1989), color = "grey50", linetype = 2, linewidth = 0.25) +
  geom_line(aes(x=BY, y=lnRS, col=ocean_region_lab)) + 
  facet_grid(rows=vars(Stock), switch ="y", scales="free_y", as.table=F) + 
  scale_colour_manual(values=col.region) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.y.left = element_text(angle=0, hjust=0, margin=margin(l=0, r=0)),
        strip.background = element_rect(fill="transparent", colour="transparent"),
        strip.text = element_text(size=7, ),
        panel.spacing.y = unit(0, unit="cm"),
        panel.background = element_rect(fill="white"),
        legend.position = "none") +
  labs(y="", x="Brood Year")
pdf(here(fig.dir, "2_era_tslength.pdf"))
print(g)
dev.off()


## 3 eras : pre- 88/89, 89-2010, post- 2010/11 ---------------------------------
## Get data for Stan

stan.dat.2c.3 <- stan_data_dyn(data_master, 
                             var.x2 = "early_sst_stnd",
                             var.x3 = "np_pinks_sec_stnd",
                             breakpoint1 = 1989,
                             breakpoint2 = 2011,
                             var.region="Ocean.Region2", 
                             scale.x1 = TRUE,
                             alpha.group = ifelse(speciesFlag=="sockeye", TRUE, FALSE))

pars_era_2c_3 <- c("alpha", "beta", "sigma", "phi", "mu_alpha", "sigma_alpha",
                 "gamma1", "gamma2", "gamma3", 
                 "mu_gamma1", "mu_gamma2", "mu_gamma3", "sigma_gamma", 
                 "kappa1", "kappa2", "kappa3", 
                 "mu_kappa1", "mu_kappa2", "mu_kappa3", "sigma_kappa" )
pars_dyn_2c_3 <- c("alpha", "beta", "sigma", "phi", "mu_alpha", "sigma_alpha",
                 "gamma", "sigma_gamma", "signal_noise_g",
                 "kappa", "sigma_kappa", "signal_noise_k")

pars.gen.quant <- c("log_lik", "yrep") ## Generated quantities to monitor

era.2c.3 <- rstan::stan(file = "./stan/hbm_era_2c.stan",
                      data = stan.dat.2c.3,
                      pars = c(pars_era_2c_3, pars.gen.quant),
                      warmup = 1000,
                      iter = 2000,
                      cores = 4,
                      chains = 4,
                      seed = 123,
                      control = list(adapt_delta = 0.99,
                                     max_treedepth = 20))
save(era.2c.3, file = here(fit.dir, "hbm_era_2c_3.RData"))

### --- Data

# Density dataframe - by stock
dens.df.st.2c.3 <- era_density_df(era.2c.3, par=c("gamma", "kappa"), neras=3)
dens.df.st.2c.3 <- ocean_region_lab(dens.df.st.2c.3)
# Density dataframe (regional lvl)
dens.df.reg.2c.3 <- era_density_df(era.2c.3, par=c("gamma", "kappa"), mu=T, neras=3)
dens.df.reg.2c.3 <- ocean_region_lab(dens.df.reg.2c.3)

### --- Figures

## Density plot
g <- ggplot(dens.df.st.2c.3) + 
  geom_vline(xintercept = 0, color = "grey50", linetype = 2, linewidth = 0.25) +
  geom_path(aes(x=x, y=dens, group=stock, col=ocean_region_lab), alpha=0.2) + 
  geom_path(data=dens.df.reg.2c.3, aes(x=x, y=dens, col=ocean_region_lab), alpha=0.85, linewidth=1) +
  scale_colour_manual(values=col.region) +
  facet_grid(rows=vars(era), cols=vars(varnam)) + 
  coord_cartesian(xlim=c(-1, 1)) +
  theme_minimal() + labs(x="covariate effect", y="", col="Ocean Region") +
  theme(axis.text.y=element_blank(), 
        legend.position="bottom")
pdf(here(fig.dir, "3_era_dens.pdf"))
print(g)
dev.off()


# Time series length - show overlap
prod_dat <- ocean_region_lab(data_master)
g <- ggplot(prod_dat) + 
  geom_vline(xintercept=c(1989, 2011), color = "grey50", linetype = 2, linewidth = 0.25) +
  geom_line(aes(x=BY, y=lnRS, col=ocean_region_lab)) + 
  facet_grid(rows=vars(Stock), switch ="y", scales="free_y", as.table=F) + 
  scale_colour_manual(values=col.region) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.y.left = element_text(angle=0, hjust=0, margin=margin(l=0, r=0)),
        strip.background = element_rect(fill="transparent", colour="transparent"),
        strip.text = element_text(size=7, ),
        panel.spacing.y = unit(0, unit="cm"),
        panel.background = element_rect(fill="white"),
        legend.position = "none") +
  labs(y="", x="Brood Year")
pdf(here(fig.dir, "3_era_tslength.pdf"))
print(g)
dev.off()



## 3 eras : pre- 88/89, 89-2010, post- 2010/11 w/ TRIMMED data ------------------
## Get data for Stan

## # Data trimming function - to test sensitivity to short 'exposure' to eras
data_trimmed <- trim.era.ts(data_master, breakpoint1 = 1989, breakpoint2 = 2011)

stan.dat.2c.3t <- stan_data_dyn(data_trimmed, 
                               var.x2 = "early_sst_stnd",
                               var.x3 = "np_pinks_sec_stnd",
                               breakpoint1 = 1989,
                               breakpoint2 = 2011,
                               var.region="Ocean.Region2", 
                               scale.x1 = TRUE,
                               alpha.group = ifelse(speciesFlag=="sockeye", TRUE, FALSE))

pars_era_2c_3t <- c("alpha", "beta", "sigma", "phi", "mu_alpha", "sigma_alpha",
                   "gamma1", "gamma2", "gamma3", 
                   "mu_gamma1", "mu_gamma2", "mu_gamma3", "sigma_gamma", 
                   "kappa1", "kappa2", "kappa3", 
                   "mu_kappa1", "mu_kappa2", "mu_kappa3", "sigma_kappa" )
pars_dyn_2c_3t <- c("alpha", "beta", "sigma", "phi", "mu_alpha", "sigma_alpha",
                   "gamma", "sigma_gamma", "signal_noise_g",
                   "kappa", "sigma_kappa", "signal_noise_k")

pars.gen.quant <- c("log_lik", "yrep") ## Generated quantities to monitor

era.2c.3t <- rstan::stan(file = "./stan/hbm_era_2c.stan",
                        data = stan.dat.2c.3t,
                        pars = c(pars_era_2c_3t, pars.gen.quant),
                        warmup = 1000,
                        iter = 2000,
                        cores = 4,
                        chains = 4,
                        seed = 123,
                        control = list(adapt_delta = 0.99,
                                       max_treedepth = 20))
save(era.2c.3t, file = here(fit.dir, "hbm_era_2c_3t.RData"))

### --- Data

# Density dataframe - by stock
dens.df.st.2c.3t <- era_density_df(era.2c.3t, par=c("gamma", "kappa"), neras=3)
dens.df.st.2c.3t <- ocean_region_lab(dens.df.st.2c.3t)
# Density dataframe (regional lvl)
dens.df.reg.2c.3t <- era_density_df(era.2c.3t, par=c("gamma", "kappa"), mu=T, neras=3)
dens.df.reg.2c.3t <- ocean_region_lab(dens.df.reg.2c.3t)

### --- Figures

## Density plot
g <- ggplot(dens.df.st.2c.3t) + 
  geom_vline(xintercept = 0, color = "grey50", linetype = 2, linewidth = 0.25) +
  geom_path(aes(x=x, y=dens, group=stock, col=ocean_region_lab), alpha=0.2) + 
  geom_path(data=dens.df.reg.2c.3t, aes(x=x, y=dens, col=ocean_region_lab), alpha=0.85, linewidth=1) +
  scale_colour_manual(values=col.region) +
  facet_grid(rows=vars(era), cols=vars(varnam)) + 
  coord_cartesian(xlim=c(-1, 1)) +
  theme_minimal() + labs(x="covariate effect", y="", col="Ocean Region") +
  theme(axis.text.y=element_blank(),
        legend.position="bottom")
pdf(here(fig.dir, "3t_era_dens.pdf"))
print(g)
dev.off()


# Time series length - show overlap
prod_dat_t <- ocean_region_lab(data_trimmed)
g <- ggplot(prod_dat_t) + 
  geom_vline(xintercept=c(1989, 2011), color = "grey50", linetype = 2, linewidth = 0.25) +
  geom_line(aes(x=BY, y=lnRS, col=ocean_region_lab)) + 
  facet_grid(rows=vars(Stock), switch ="y", scales="free_y", as.table=F) + 
  scale_colour_manual(values=col.region) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.y.left = element_text(angle=0, hjust=0, margin=margin(l=0, r=0)),
        strip.background = element_rect(fill="transparent", colour="transparent"),
        strip.text = element_text(size=7, ),
        panel.spacing.y = unit(0, unit="cm"),
        panel.background = element_rect(fill="white"),
        legend.position = "none") +
  labs(y="", x="Brood Year")
pdf(here(fig.dir, "3t_era_tslength.pdf"))
print(g)
dev.off()


