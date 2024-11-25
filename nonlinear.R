# Nonlinearities analysis

# Set up
fig.dir <- here('hannahland', 'nonlin-stnd') #here('hannahland')

raw.clim <- if(speciesFlag=="sockeye"){
  raw.clim.sock} else if(speciesFlag=="pink") {
    raw.clim.pink } else if (speciesFlag=="chum"){
     raw.clim.chum  }

raw.sst.index <- clim.wgt.avg(brood.table = data_master,
                              env.data = raw.clim,
                              env.covar = "sst_raw",
                              type = "first_year",
                              out.covar = "early_sst")

cov.dat <- full_join(raw.sst.index, data_master[,c("Stock.ID", "Ocean.Region2","early_sst_stnd", "np_pinks_sec", "BY", "lnRS")] , by=c("Stock.ID", "BY"))

# Show the relationship between raw SST and stnd sst
ggplot(cov.dat) + geom_point(aes(x=early_sst, y=early_sst_stnd, col=Ocean.Region2))

# Plot SST density by region # more different than by era
g <- ggplot(cov.dat) + 
  geom_density(aes(x=early_sst, col=Ocean.Region2, fill=Ocean.Region2), alpha=0.1, linewidth=1) +
  scale_colour_manual(values=rev(unname(col.region)), aesthetics=c("colour", "fill")) + 
  theme_sleek() + labs(x="SST Index (°C)", y="Density", col="Region", fill="Region")
png(filename=here(fig.dir, 'sst_dens.png'), width=801*2, height=486*2, res=72*2)
print(g)
dev.off()



# Load fit with no covariates
load(file=here('output', 'models', 'stat', speciesFlag, 'stat_nocov.RData'), verbose=T)
summ_nc <- rstan::summary(nocov, pars="yresid", probs=c(0.025, 0.5, 0.975))$summary

# Also load residuals from fit with covariates
load(file=here('output', 'models', 'stat', speciesFlag, 'stat_sub.RData'), verbose=T) ## This should be changed eventually to 'stat_a'
summ_c <- rstan::summary(stat_sub, pars="yresid", probs=NULL)$summary

# Bind with data
cov.dat.resid <- cbind(cov.dat, summ_nc[,"mean"], summ_c[,"mean"])
names(cov.dat.resid) <- c(names(cov.dat), "resids", "cov_resids")
cov.dat.resid <- ocean_region_lab(cov.dat.resid)


## Standardize residuals
cov.dat.resid <- ddply(cov.dat.resid, .(Stock.ID), transform,
                          stnd_resids = scale(resids),
                          stnd_cov_resids = scale(cov_resids))

cov.dat.resid %>% dplyr::summarize(mean_raw = mean(resids), sd_raw = sd(resids), mean_stnd=mean(stnd_resids), sd_stnd=sd(stnd_resids), .by=c(Ocean.Region2))

#### Residuals versus time ####

# 1. Residuals from stationary model without covariates
g <- ggplot(cov.dat.resid) + 
    geom_point(aes(x=BY, y=stnd_resids, col=ocean_region_lab), alpha=0.2, shape=16) + 
    hline_0(alpha=0.2) + geom_vline(xintercept=c(1989, 2011), alpha=0.2, lty="dashed") +
    geom_smooth(aes(x=BY, y=stnd_resids, col=ocean_region_lab)) +
    scale_colour_manual(values=col.region, guide="none") +
    theme_sleek() + labs(x="Brood Year", y="Residuals", title="Residuals (no covariates) vs. time") + 
    facet_wrap(vars(ocean_region_lab), scales = "free_y") 
png(filename=here(fig.dir, 'no_cov_stnd_residsxtime.png'), width=801*2, height=486*2, res=72*3)
print(g)
dev.off()

# 2. Residuals from stationary model with covariates
g <- ggplot(cov.dat.resid) + 
    geom_point(aes(x=BY, y=stnd_cov_resids, col=ocean_region_lab), alpha=0.2, shape=16) + 
    hline_0(alpha=0.2) + geom_vline(xintercept=c(1989, 2011), alpha=0.2, lty="dashed") +
    geom_smooth(aes(x=BY, y=stnd_resids), col="gray50", alpha=0.7, se=F) +
    geom_smooth(aes(x=BY, y=stnd_cov_resids, col=ocean_region_lab)) +
    scale_colour_manual(values=col.region, guide="none") +
    theme_sleek() + labs(x="Brood Year", y="Residuals", title="Residuals (w/ covariates) vs. time") + 
    facet_wrap(vars(ocean_region_lab), scales = "free_y")
png(filename=here(fig.dir, 'cov_residsxtime.png'), width=801*2, height=486*2, res=72*3)
print(g)
dev.off()


#### Residuals versus covariates ####

# 1. Residuals vs SST by region
g <- ggplot(cov.dat.resid) + geom_point(aes(x=early_sst, y=stnd_resids, col=ocean_region_lab), alpha=0.2, shape=16) +
      hline_0(alpha=0.2) +
      geom_smooth(aes(x=early_sst, y=stnd_resids, col=ocean_region_lab)) + 
      theme_sleek() + labs(x="SST (degrees C)", y="Residuals", title="Residuals (no covariates) vs. SST") + 
      scale_colour_manual(values=col.region, guide="none") +
      facet_wrap(vars(ocean_region_lab), scales="free_y")
png(filename=here(fig.dir, 'no_cov_residsxsst.png'), width=801*2, height=486*2, res=72*3)
print(g)
dev.off()

# 2. Residuals vs SST - 'global'
g <- ggplot(cov.dat.resid) + geom_point(aes(x=early_sst, y=stnd_resids, col=ocean_region_lab), alpha=0.2, shape=16) +
      hline_0(alpha=0.2) +
      geom_smooth(aes(x=early_sst, y=stnd_resids), col="grey50") +
      theme_sleek() + labs(x="SST (degrees C)", y="Residuals", title="Residuals (no covariates) vs. SST", col="Region") + 
      scale_colour_manual(values=col.region) 
png(filename=here(fig.dir, 'no_cov_residsxsst_global.png'), width=801*2, height=486*2, res=72*3)
print(g)
dev.off()

# 3. Residuals vs Pinks by region
g <- ggplot(cov.dat.resid) + geom_point(aes(x=np_pinks_sec, y=stnd_resids, col=ocean_region_lab), alpha=0.2, shape=16) + 
      hline_0(alpha=0.2) +
      geom_smooth(aes(x=np_pinks_sec, y=stnd_resids, col=ocean_region_lab), alpha=0.5) +
      theme_sleek() + labs(x="NP Pink Abundance", y="Residuals", title="Residuals (no covariates) vs. NP Pinks") + 
      scale_colour_manual(values=col.region, guide="none") +
      facet_wrap(vars(ocean_region_lab), scales="free_y") 
png(filename=here(fig.dir, 'no_cov_residsxpinks.png'), width=801*2, height=486*2, res=72*3)
print(g)
dev.off()

# 4. 2-D covariates (SST and pinks) with residual size
g <- ggplot(cov.dat.resid) + 
      geom_point(aes(x=early_sst, y=BY, col=stnd_resids, size=abs(stnd_resids)), alpha=0.7) + 
      theme_sleek() +
      guides(size="none") +
      scale_colour_distiller(palette = "RdBu") +
      labs(x="SST (degrees C)", y="NP Pink abundance", title="Residuals vs. covariates", col="Residual", size="Residual")  
png(filename=here(fig.dir, 'residsxcovars.png'), width=801*2, height=486*2, res=72*3)
print(g)
dev.off()


#### Resids vs covars by era ####

# Fitting a simple lm for the lines - NOT hbm estimated rels.
g <- cov.dat.resid %>%
  mutate(era = case_when(BY<1989 ~ "Early",
                          BY >= 1989 & BY <2011 ~ "Middle",
                          BY >= 2011 ~ "Late")) %>%
  mutate(regxera = paste(ocean_region_lab, era, sep=".")) %>%
  ggplot() + hline_0(col="gray70") + 
  geom_point(aes(x=early_sst, y=stnd_resids, col=factor(regxera, levels=unique(regxera))), alpha=0.2, shape=16) +
  geom_smooth(aes(x=early_sst, y=stnd_resids, col=factor(regxera, levels=unique(regxera))), method="lm", se=F) +
  scale_colour_manual(values=col.eras) +
  facet_wrap(vars(ocean_region_lab), scales="free_y") +
  labs(x="SST (degrees C)", y="Residuals", col="Region x Era") +
  theme_sleek()
png(filename=here(fig.dir, 'eras_residsxsst.png'), width=801*2, height=486*2, res=72*3)
print(g)
dev.off()

g <- cov.dat.resid %>%
  mutate(era = case_when(BY<1989 ~ "Early",
                         BY >= 1989 & BY <2011 ~ "Middle",
                         BY >= 2011 ~ "Late")) %>%
  mutate(regxera = paste(ocean_region_lab, era, sep=".")) %>%
  ggplot() + hline_0(col="gray70") + 
  geom_point(aes(x=np_pinks_sec, y=stnd_resids, col=factor(regxera, levels=unique(regxera))), alpha=0.2, shape=16) +
  geom_smooth(aes(x=np_pinks_sec, y=stnd_resids, col=factor(regxera, levels=unique(regxera))), method="lm", se=F) +
  scale_colour_manual(values=col.eras) +
  facet_wrap(vars(ocean_region_lab), scales="free_y") +
  labs(x="NP Pink abundance", y="Residuals", col="Region x Era") +
  theme_sleek()
png(filename=here(fig.dir, 'eras_residsxcomp.png'), width=801*2, height=486*2, res=72*3)
print(g)
dev.off()


#### Detrended covariates ####
fig.dir <- here('hannahland')
#source("detrend_covars.R")

# Load detrended model fit
load(file = "./output/models/stat_detrend.RData", verbose=T)
df.reg.dt <- era_hb_param_df(era_detrend, par=c("gamma", "kappa"), mu=T, neras=3)


# load *not detrended* era model fits
load(file=here('output', 'models', 'dyn', speciesFlag, 'hbm_era_2c_3sub.RData'), verbose=T)
era.reg <- era_hb_param_df(era.2c.3sub, par=c("gamma", "kappa"), mu = TRUE)

### --- Figures

## Dot-and-whiskers to compare detrended to not detrended era model
g <- ggplot(ocean_region_lab(df.reg.dt)) + 
  geom_hline(yintercept = 0, color = "grey50", linetype = 2, linewidth = 0.25) +
  geom_point(aes(y=reg_mean, x=era, col=ocean_region_lab), shape=18, size=3) +
  geom_segment(aes(y=lower_10, yend=upper_90, x=era, xend=era, col=ocean_region_lab), 
               linewidth=.8) +
  geom_point(data=ocean_region_lab(era.reg), aes(y=reg_mean, x=era, col=ocean_region_lab), 
             position=position_nudge(x=.1), shape=15, size=2) + 
  geom_segment(data=ocean_region_lab(era.reg), aes(y=lower_10, yend=upper_90, x=era, xend=era, col=ocean_region_lab), 
               position=position_nudge(x=.1), linewidth=.8) + # squares = original era estimates!
  facet_grid(cols=vars(varnam), rows=vars(ocean_region_lab), scales="free_y") +
  scale_colour_manual(values=col.region, guide="none") +
  labs(x="Era", y="Covariate effect") +  theme_sleek()
png(filename=here(fig.dir, 'eras_detrend_compare.png'), width=600*2, height=500*2, res=72*3)
print(g)
dev.off()



# Get detrended model residuals
summ_dt <- rstan::summary(stat_detrend, pars="yresid", probs=NULL)$summary

# Bind with data
detrend.resid <- cbind(detrend, summ_dt[,"mean"])
names(detrend.resid) <- c(names(detrend.resid)[1:(ncol(detrend.resid)-1)], "cov_resids")
detrend.resid <- ocean_region_lab(detrend.resid)

### Plot detrended covariates
ggplot(detrend.resid) + geom_line(aes(x=BY, y=np_pinks_gam_sec, group=Stock)) + 
  facet_wrap(vars(ocean_region_lab))


#### Residuals versus time ####

g <- ggplot(detrend.resid) + 
  geom_point(aes(x=BY, y=cov_resids, col=ocean_region_lab), alpha=0.2, shape=16) + 
  hline_0(alpha=0.2) + geom_vline(xintercept=c(1989, 2011), alpha=0.2, lty="dashed") +
  geom_smooth(aes(x=BY, y=cov_resids, col=ocean_region_lab)) +
  scale_colour_manual(values=col.region, guide="none") +
  theme_sleek() + labs(x="Brood Year", y="Residuals", title="Residuals vs. time (detrended covars)") + 
  facet_wrap(vars(ocean_region_lab), scales = "free_y") 
png(filename=here(fig.dir, 'detrended_cov_residsxtime.png'), width=801*2, height=486*2, res=72*3)
print(g)
dev.off()

#### Residuals vs covariates ####

# 1. Residuals vs SST by region
g <- ggplot(detrend.resid) + geom_point(aes(x=early_sst_lm, y=cov_resids, col=ocean_region_lab), alpha=0.2, shape=16) +
  hline_0(alpha=0.2) +
  geom_smooth(aes(x=early_sst_lm, y=cov_resids, col=ocean_region_lab)) +
  theme_sleek() + labs(x="Detrended SST", y="Residuals", title="Residuals (detrended covariates) vs. SST") + 
  scale_colour_manual(values=col.region, guide="none") +
  facet_wrap(vars(ocean_region_lab), scales="free_y")
png(filename=here(fig.dir, 'detrended_cov_residsxsst.png'), width=801*2, height=486*2, res=72*3)
print(g)
dev.off()

# 2. Residuals vs SST - 'global' # NA for now; unsure if matching to 'raw' sst is appropriate

# 3. Residuals vs Pinks by region
g <- ggplot(detrend.resid) + geom_point(aes(x=np_pinks_gam_sec, y=cov_resids, col=ocean_region_lab), alpha=0.2, shape=16) + 
  hline_0(alpha=0.2) +
  geom_smooth(aes(x=np_pinks_gam_sec, y=cov_resids, col=ocean_region_lab), alpha=0.5) +
  theme_sleek() + labs(x="Detrended NP Pink Abundance", y="Residuals", title="Residuals (detrended covariates) vs. NP Pinks") + 
  scale_colour_manual(values=col.region, guide="none") +
  facet_wrap(vars(ocean_region_lab), scales="free_y") 
png(filename=here(fig.dir, 'detrended_cov_residsxcomp.png'), width=801*2, height=486*2, res=72*3)
print(g)
dev.off()

#### Detrended eras model ####
load( file = "./output/models/era_detrend.RData", verbose=T)
summ_era_dt <- rstan::summary(era_detrend, pars="yresid", probs=NULL)$summary

# Bind with data
detrend.resid.era <- cbind(detrend, summ_era_dt[,"mean"])
names(detrend.resid.era) <- c(names(detrend.resid.era)[1:(ncol(detrend.resid.era)-1)], "cov_resids")
detrend.resid.era <- ocean_region_lab(detrend.resid.era)

# Residuals over time
ggplot(detrend.resid.era) + 
  geom_point(aes(x=BY, y=cov_resids, col=ocean_region_lab), alpha=0.2, shape=16) + 
  hline_0(alpha=0.2) + geom_vline(xintercept=c(1989, 2011), alpha=0.2, lty="dashed") +
  geom_smooth(aes(x=BY, y=cov_resids, col=ocean_region_lab)) +
  scale_colour_manual(values=col.region, guide="none") +
  theme_sleek() + labs(x="Brood Year", y="Residuals", title="Residuals vs. time (detrended covars)") + 
  facet_wrap(vars(ocean_region_lab), scales = "free_y")

# Residuals vs covariate by era
detrend.resid.era %>% mutate(era = case_when(BY<1989 ~ "Early",
                                       BY >= 1989 & BY <2011 ~ "Middle",
                                       BY >= 2011 ~ "Late")) %>%
ggplot() + 
  geom_point(aes(x=early_sst_lm, y=cov_resids, col=ocean_region_lab), alpha=0.2, shape=16) + 
  hline_0(alpha=0.2) + 
  #geom_vline(xintercept=c(1989, 2011), alpha=0.2, lty="dashed") +
  geom_smooth(aes(x=early_sst_lm, y=cov_resids, group=era, col=ocean_region_lab), method="lm") +
  scale_colour_manual(values=col.region, guide="none") +
  theme_sleek() + labs(x="Brood Year", y="Residuals", title="Residuals vs. time (detrended covars)") + 
  facet_wrap(vars(ocean_region_lab), scales = "free_y")
# is this interesting?

