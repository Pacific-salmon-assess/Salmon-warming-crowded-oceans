## Re-fitting era models w/ different eras

# Set species
#speciesFlag = "pink"
speciesFlag = "chum"
#speciesFlag = "sockeye"

if(speciesFlag=="pink") 
  data_master <- pink else if(speciesFlag=="chum") 
    data_master <- chum else if(speciesFlag=="sockeye")
      data_master <- sock
    

dir <- here('hannahland', 'era-fits', speciesFlag)

if(!exists(dir)) dir.create(dir, recursive=T)


## Get data for Stan
stan.dat.2c <- stan_data_dyn(data_master, 
                             var.x2 = "early_sst_stnd",
                             var.x3 = "np_pinks_sec_stnd",
                             breakpoint1 = 1989,
                             breakpoint2 = NULL,
                             var.region="Ocean.Region2", 
                             scale.x1 = TRUE,
                             alpha.group = FALSE)


pars_era_2c <- c("alpha", "beta", "sigma", "phi", "mu_alpha", "sigma_alpha",
                 "gamma1", "gamma2",
                 "mu_gamma1", "mu_gamma2", "sigma_gamma", 
                 "kappa1", "kappa2",
                 "mu_kappa1", "mu_kappa2", "sigma_kappa" )
pars.gen.quant <- c("log_lik", "yrep") ## Generated quantities to monitor


# Run MCMC 
era.2c <- rstan::stan(file = "./stan/hbm_era_2c.stan",
                      data = stan.dat.2c,
                      pars = c(pars_era_2c, pars.gen.quant),
                      warmup = 1000,
                      iter = 2000,
                      cores = 4,
                      chains = 4,
                      seed = 123,
                      control = list(adapt_delta = 0.99,
                                     max_treedepth = 20))
save(era.2c, file = here(dir, "hbm_era_2c_2.RData"))


### --- Eras model: Data

# Density dataframe - by stock
dens.df.st.2c <- era_density_df(era.2c, par=c("gamma", "kappa"), neras=2)
dens.df.st.2c <- ocean_region_lab(dens.df.st.2c)
# Density dataframe (regional lvl)
dens.df.reg.2c <- era_density_df(era.2c, par=c("gamma", "kappa"), mu=T, neras=2)
dens.df.reg.2c <- ocean_region_lab(dens.df.reg.2c)

### --- Eras model: Figures

## Density plot
g <- ggplot(dens.df.st.2c) + 
  geom_vline(xintercept = 0, color = "grey50", linetype = 2, linewidth = 0.25) +
  geom_path(aes(x=x, y=dens, group=stock, col=ocean_region_lab), alpha=0.2) + 
  geom_path(data=dens.df.reg.2c, aes(x=x, y=dens, col=ocean_region_lab), alpha=0.85, linewidth=1) +
  scale_colour_manual(values=col.region) +
  facet_grid(rows=vars(era), cols=vars(varnam)) + 
  coord_cartesian(xlim=c(-1, 1)) +
  theme_minimal() + labs(x="covariate effect", y="", col="Ocean Region") +
  theme(axis.text.y=element_blank())
pdf(here(dir, "2_era_dens.pdf"))
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
pdf(here(dir, "2_era_tslength.pdf"))
print(g)
dev.off()

# Data trimming function - to test sensitivity to short 'exposure' to eras
data_trimmed <- trim.era.ts()


