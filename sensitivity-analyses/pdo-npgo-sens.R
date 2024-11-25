## Sensitivity analysis of including NPGO, PDO in models

# This script could be organized more efficiently - e.g. put plotting sections
# together so base model fits are only loaded once.


# Specify output directories
sens.fig.dir <- here("sensitivity-analyses", "pdo-npgo", speciesFlag)
sens.fit.dir <- here("sensitivity-analyses", "fits", speciesFlag)

# Set up
brood <- if(speciesFlag=="sockeye"){
  bt.complete} else if(speciesFlag=="pink") {
    bt.complete.pink } else if (speciesFlag=="chum"){
      bt.complete.chum  }


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


#### Exploratory plots ####

# Explore correlation between PDO and other covariates
# By region and by time period

pdf(file=here(sens.fig.dir, "pdo-npgo-covar-corr.pdf"))

## create empty 3D array
array.cor <- array(NA, dim = c(4,4,nrow(info_master)))

## calculate stock specific covar correlations
for(i in seq_along(info_master$Stock.ID)) { 
  stk.i <- br.pdo.npgo[br.pdo.npgo$Stock.ID == info_master$Stock.ID[i], ]
  covar.i <- subset(stk.i, select = c("early_sst", "np_pinks_sec", "pdo_index", "npgo_index"))
  cor.i <- cor(covar.i, use = "pairwise.complete.obs")
  array.cor[ , , i] <- cor.i
}

## Average across stocks
cor.covars <- apply(array.cor, c(1, 2), mean)
row.names(cor.covars) <- c("early_sst", "np_pinks_sec", "pdo_index", "npgo_index")
colnames(cor.covars) <- c("early_sst", "np_pinks_sec", "pdo_index", "npgo_index")

corrplot::corrplot(cor.covars, method="color", type="upper",
         addCoef.col = "black", tl.col="black", 
         diag=F,
         col=rev(corrplot::COL2('RdBu', 200)),
         title="Average stock-specific covariate correlations (All years)",
         mar=c(0,0,2,0))

# Split by era
era.yrs <- list(start=c(1950,1989,2011),
                end=c(1988,2010,2019))
for(e in 1:3){

  ## create empty 3D array
  array.cor <- array(NA, dim = c(4,4,nrow(info_master)))
  ## calculate stock specific covar correlations
  for(i in seq_along(info_master$Stock.ID)) { 
    stk.i <- filter(br.pdo.npgo, Stock.ID==info_master$Stock.ID[i], BY %in% era.yrs$start[e]:era.yrs$end[e])
    if(empty(stk.i)) next
    covar.i <- subset(stk.i, select = c("early_sst", "np_pinks_sec", "pdo_index", "npgo_index"))
    cor.i <- cor(covar.i, use = "pairwise.complete.obs")
    array.cor[ , , i] <- cor.i
  }
  
  ## Average across stocks
  cor.covars <- apply(array.cor, c(1, 2), mean, na.rm=T)
  row.names(cor.covars) <- c("early_sst", "np_pinks_sec", "pdo_index", "npgo_index")
  colnames(cor.covars) <- c("early_sst", "np_pinks_sec", "pdo_index", "npgo_index")
  
  g <- corrplot::corrplot(cor.covars, method="color", type="upper",
                addCoef.col = "black", tl.col="black", 
                diag=F,
                col=rev(corrplot::COL2('RdBu', 200)),
                title=paste("Average stock-specific covariate correlations", 
                         era.yrs$start[e], "-", era.yrs$end[e]),
                mar=c(0,0,2,0))
  print(g)
}
dev.off()


## Explore correlation between PDO, NPGO and lnRS
era.yrs <- list(start=c(1950,1989,2011),
                end=c(1988,2010,2019))

pdf(file= here(sens.fig.dir, "pdo-npgo-prod-corr.pdf"))

cor.stock <- plyr::ddply(br.pdo.npgo, .(Stock.ID), plyr::summarize,
                         Ocean.Region2 = unique(Ocean.Region2),
                         early_sst = cor(lnRS, early_sst, use = "pairwise.complete.obs"),
                         np_pinks_sec = cor(lnRS, np_pinks_sec, use = "pairwise.complete.obs"),
                         npgo = cor(lnRS, npgo_index, use="pairwise.complete.obs"),
                         pdo = cor(lnRS, pdo_index, use="pairwise.complete.obs"))

cor.stock$Stock.ID <- NULL
cor.stock <- reshape2::melt(cor.stock, id.vars = "Ocean.Region2")

cor.ocean <- plyr::ddply(cor.stock, .(Ocean.Region2, variable), summarize,
                         cor.avg = mean(value))

g <- ggplot(ocean_region_lab(cor.ocean)) + 
  geom_col(aes(x=variable, y=cor.avg, fill=ocean_region_lab), position=position_dodge(), alpha=0.7) +
  scale_fill_manual(values=col.region) +
  theme_sleek() + labs(x="", y="Average Correlation", fill="Region",
                       title="All Years")
print(g)

for(e in 1:3) { # loop over eras
    dat <-br.pdo.npgo %>% filter(BY %in% era.yrs$start[e]:era.yrs$end[e]) 
    cor.stock <- plyr::ddply(dat, .(Stock.ID), plyr::summarize,
                             Ocean.Region2 = unique(Ocean.Region2),
                             early_sst = cor(lnRS, early_sst, use = "pairwise.complete.obs"),
                             np_pinks_sec = cor(lnRS, np_pinks_sec, use = "pairwise.complete.obs"),
                             npgo = cor(lnRS, npgo_index, use="pairwise.complete.obs"),
                             pdo = cor(lnRS, pdo_index, use="pairwise.complete.obs"))
    
    cor.stock$Stock.ID <- NULL
    cor.stock <- reshape2::melt(cor.stock, id.vars = "Ocean.Region2")
    
    cor.ocean <- plyr::ddply(cor.stock, .(Ocean.Region2, variable), summarize,
                             cor.avg = mean(value))
    g <- ggplot(ocean_region_lab(cor.ocean)) + 
      geom_col(aes(x=variable, y=cor.avg, fill=ocean_region_lab), position=position_dodge(), alpha=0.7) +
      scale_fill_manual(values=col.region) +
      theme_sleek() + labs(x="", y="Average Correlation", fill="Region",
                           title=paste("Era:", era.yrs$start[e], "-", era.yrs$end[e]))
    
    print(g)
}
dev.off()


## No-covariate model residuals with PDO overlaid 

# Base model w/o covariates may not exist for all species currently, skip if it doesn't
if(file.exists(here('output', 'models', 'stat', speciesFlag, 'stat_nocov.RData'))) {
  # Load fit with no covariates
  load(file=here('output', 'models', 'stat', speciesFlag, 'stat_nocov.RData'), verbose=T)
  summ_nc <- rstan::summary(nocov, pars="yresid", probs=NULL)$summary
  
  # Also load residuals from fit with covariates
  load(file=here('output', 'models', 'stat', speciesFlag, 'stat_a.RData'), verbose=T)
  summ_c <- rstan::summary(stat_a, pars="yresid", probs=NULL)$summary
  
  # Bind residuals with brood table
  resids <- cbind(br.pdo, br.npgo[,"npgo_index"], summ_c[,"mean"], summ_nc[,"mean"])
  names(resids) <- c(names(br.pdo), "npgo_index", "2c_resids", "base_resids")
  
  
  pdf(file=here(sens.fig.dir, "pdo-npgo-rec-resids.pdf"))
  # Plot PDO index over recruitment residuals 
g <-  ggplot(ocean_region_lab(resids)) + 
    geom_point(aes(x=BY, y=base_resids, col=ocean_region_lab), alpha=0.2, shape=16) + 
    hline_0(alpha=0.2) + geom_vline(xintercept=c(1989, 2011), alpha=0.2, lty="dashed") +
    geom_smooth(aes(x=BY, y=pdo_index), alpha=0.2) +
    scale_colour_manual(values=col.region, guide="none") +
    theme_sleek() + labs(x="Brood Year", y="Residuals", 
                         title="Recruitment Residuals vs. time w/ Dec-Mar PDO Index") + 
    facet_wrap(vars(Ocean.Region2), scales = "free_y") 
print(g)
  
  # Plot NPGO over recruitment residuals
g <- ggplot(ocean_region_lab(resids)) + 
    geom_point(aes(x=BY, y=base_resids, col=ocean_region_lab), alpha=0.2, shape=16) + 
    hline_0(alpha=0.2) + geom_vline(xintercept=c(1989, 2011), alpha=0.2, lty="dashed") +
    geom_smooth(aes(x=BY, y=npgo_index), alpha=0.2) +
    scale_colour_manual(values=col.region, guide="none") +
    theme_sleek() + labs(x="Brood Year", y="Residuals", 
                         title="Recruitment Residuals vs. time w/ Dec-Mar NPGO Index") + 
    facet_wrap(vars(Ocean.Region2), scales = "free_y") 
print(g)

  dev.off()
}


#### -- MCMC - NPGO -- ####

## Run MCMC, stationary model
stan.dat.npgo <- stan_data_stat(data_master,
                               scale.x1 = TRUE,
                               var.x2 = "early_sst_stnd",
                               var.x3 = "np_pinks_sec_stnd",
                               var.region = "Ocean.Region2",
                               alpha.group = ifelse(speciesFlag=="sockeye", TRUE, FALSE))
stan.dat.npgo$x4 <- br.npgo$npgo_stnd

stat_npgo <- rstan::stan(file = "./stan/hbm_stat_3c.stan",
                      data = stan.dat.npgo,
                      #pars = c(pars.stat, pars.gen.quant),
                      warmup = 1000,
                      iter = 2000,
                      cores = 4,
                      chains = 4,
                      seed = 123,
                      control = list(adapt_delta = 0.99,
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


#### ---- Visualize fits ---- #####

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
# This should be a function

#NPGO
gamma <- rstan::summary(stat_npgo, pars = "mu_gamma")$summary
kappa <- rstan::summary(stat_npgo, pars = "mu_kappa")$summary
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
tab.l <- data.frame(reg = reg,
                    coef = "NPGO",
                    lower = lambda[ , "2.5%"],
                    mean = lambda[ , "mean"],
                    upper = lambda[ , "97.5%"])

tab.coef <- rbind(tab.g, tab.k, tab.l) # add tab.c if exists
tab.coef$perc <- (exp(tab.coef$mean) - 1) * 100
row.names(tab.coef) <- NULL
names(tab.coef) <- c("Ecosystem", "Coefficient", "Lower 95% CI", "Mean",
                     "Upper 95% CI", "Mean % change in R/S")

write.csv(tab.coef, file = here(sens.fig.dir, "model_coefficients_stat_npgo.csv"))



## Visualize effects -- NPGO ## 

# Load base model results for comparison
load(here('output', 'models', 'stat', speciesFlag, 'stat_a.RData'), verbose=T)
gamma.stock <- hb_param_df(stat_a, "gamma", "Ocean.Region2", "SST", info=info_master)
kappa.stock <- hb_param_df(stat_a, "kappa", "Ocean.Region2", "Competitors", info=info_master)
df.dot <- rbind(gamma.stock, kappa.stock )
df.dot <- ocean_region_lab(df.dot, "region", FALSE)
df.dot$Stock <- factor(df.dot$Stock, levels = levels(data_master$Stock))
df.dot$var <- factor(df.dot$var, levels = c("SST", "Competitors" )) # ,"SST x Comp"))
df.dot <- dplyr::arrange(df.dot, Stock)
df.mu <- plyr::ddply(df.dot, .(region, var), summarize,
                     mu_mean = unique(mu_mean),
                     mu_2.5 = unique(`mu_2.5%`),
                     mu_97.5 = unique(`mu_97.5%`),
                     ocean_region_lab = unique(ocean_region_lab),
                     ystart = Stock[1],
                     yend = Stock[length(Stock)])

# NPGO stationary model
npgo.dot <- rbind(hb_param_df(stat_npgo, par=c("gamma"), region.var="Ocean.Region2"),
                  hb_param_df(stat_npgo, par=c("kappa"), region.var="Ocean.Region2"),
                  hb_param_df(stat_npgo, par=c("lambda"), region.var="Ocean.Region2"))
npgo.dot <- ocean_region_lab(npgo.dot, "region", FALSE)
npgo.dot$Stock <- factor(npgo.dot$Stock, levels = levels(data_master$Stock))
npgo.dot <- mutate(npgo.dot, var = case_when(par=="gamma" ~ "SST",
                                             par=="kappa" ~ "Competitors",
                                             par=="lambda" ~ "NPGO"))
npgo.dot <- dplyr::arrange(npgo.dot, Stock)
npgo.mu <- plyr::ddply(npgo.dot, .(region, var), summarize,
                       mu_mean = unique(mu_mean),
                       mu_2.5 = unique(`mu_2.5%`),
                       mu_97.5 = unique(`mu_97.5%`),
                       ocean_region_lab = unique(ocean_region_lab),
                       ystart = Stock[1],
                       yend = Stock[length(Stock)])


# Dot and whiskers plot by stock
pdf(file=here(sens.fig.dir, "npgo-stat-dot.pdf"))
g <- ggplot(npgo.dot) +
  geom_vline(xintercept = 0, color = "grey50", linetype = 2, linewidth = 0.25) +
  geom_point(aes(x = mean, y = Stock, color = ocean_region_lab, shape = ocean_region_lab)) +
  geom_segment(aes(y = Stock, yend = Stock, x = `2.5%`, xend = `97.5%`,
                   color = ocean_region_lab), linewidth = 0.25) +
  geom_segment(data = npgo.mu, aes(y = ystart, yend = yend, x = mu_mean, xend = mu_mean,
                                   color = ocean_region_lab), linewidth = 0.25) +
  geom_segment(data=df.mu, aes(y = ystart, yend = yend, x = mu_mean, xend = mu_mean,
                               color = ocean_region_lab), linewidth = 0.25, col="grey50",
               alpha=0.6) + 
  geom_rect(data = npgo.mu, aes(xmin = mu_2.5, xmax = mu_97.5, ymin = ystart,
                                ymax = yend, fill = ocean_region_lab),
            alpha = 0.2) +
  col.scale.reg +
  scale_shape_manual(values = shp.reg) +
  scale_fill_manual(values = col.region, guide="none") +
  labs(x = "Coefficient",
       y = "",
       color = "",
       shape = "") +
  facet_wrap( ~ factor(var, levels=c("SST", "Competitors", "NPGO"))) +
  scale_x_continuous(breaks=c(-.5, -.25,0,.25, .5), limits=c(-0.7, 0.7))+
  theme_sleek(base_size = 10) +
  theme(legend.justification = c(0, 0),
        legend.position = c(0.01, 0.87),
        legend.key.size = unit(10, "pt"),
        legend.background = element_blank(),
        legend.text = element_text(size = 8),
        panel.spacing.x = unit(-0.5, "pt"))
print(g)
dev.off()

# Era comparison plot with NPGO

df.reg.npgo <- era_hb_param_df(era_npgo, par=c("gamma", "kappa", "lambda"), mu=T, neras=3)
df.reg.npgo$varnam[is.na(df.reg.npgo$varnam)] <- "NPGO"

# load basic era model fits
load(here('output', 'models', 'dyn', speciesFlag, 'hbm_era_2c.RData'), verbose=T)
df.reg <- era_hb_param_df(era.2c, par=c("gamma", "kappa"), mu=T, neras=3, region.var="Ocean.Region2")

pdf(file=here(sens.fig.dir, "npgo-era-dot.pdf"))
g <- ggplot(ocean_region_lab(df.reg.npgo)) + 
  geom_hline(yintercept = 0, color = "grey50", linetype = 2, linewidth = 0.25) +
  geom_point(data=ocean_region_lab(df.reg), aes(y=reg_mean, x=era), col="grey50", shape=18, size=3, position=position_nudge(x=.1)) +
  geom_segment(data=ocean_region_lab(df.reg), aes(y=lower_10, yend=upper_90, x=era, xend=era), col="grey50",  
               linewidth=.8, position=position_nudge(x=.1)) + # plot base models for comparison
  geom_point(aes(y=reg_mean, x=era, col=ocean_region_lab), shape=18, size=3) +
  geom_segment(aes(y=lower_10, yend=upper_90, x=era, xend=era, col=ocean_region_lab), 
               linewidth=.8) +
  facet_grid(cols=vars(factor(varnam, levels=c("SST", "Competitors", "NPGO"))), rows=vars(ocean_region_lab), scales="free_y") +
  scale_colour_manual(values=col.region, guide="none") +
  labs(x="Era", y="Covariate effect") +  theme_sleek() + lims(y=c(-1,1))
print (g)
dev.off()





#### MCMC - PDO -- ####

## Run MCMC, stationary model
stan.dat.pdo <- stan_data_stat(data_master,
                                scale.x1 = TRUE,
                                var.x2 = "early_sst_stnd",
                                var.x3 = "np_pinks_sec_stnd",
                                var.region = "Ocean.Region2",
                                alpha.group = ifelse(speciesFlag=="sockeye", TRUE, FALSE))
stan.dat.pdo$x4 <- br.pdo$pdo_stnd

stat_pdo <- rstan::stan(file = "./stan/hbm_stat_3c.stan",
                         data = stan.dat.pdo,
                         #pars = c(pars.stat, pars.gen.quant),
                         warmup = 1000,
                         iter = 2000,
                         cores = 4,
                         chains = 4,
                         seed = 123,
                         control = list(adapt_delta = 0.99,
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
# This should be a function

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

write.csv(tab.coef, file = here(sens.fig.dir, "model_coefficients_stat_pdo.csv"))


## Visualize effects -- PDO ## 

# Load base model results for comparison
load(here('output', 'models', 'stat', speciesFlag, 'stat_a.RData'), verbose=T)
gamma.stock <- hb_param_df(stat_a, "gamma", "Ocean.Region2", "SST", info=info_master)
kappa.stock <- hb_param_df(stat_a, "kappa", "Ocean.Region2", "Competitors", info=info_master)
df.dot <- rbind(gamma.stock, kappa.stock )
df.dot <- ocean_region_lab(df.dot, "region", FALSE)
df.dot$Stock <- factor(df.dot$Stock, levels = levels(data_master$Stock))
df.dot$var <- factor(df.dot$var, levels = c("SST", "Competitors" )) # ,"SST x Comp"))
df.dot <- dplyr::arrange(df.dot, Stock)
df.mu <- plyr::ddply(df.dot, .(region, var), summarize,
                     mu_mean = unique(mu_mean),
                     mu_2.5 = unique(`mu_2.5%`),
                     mu_97.5 = unique(`mu_97.5%`),
                     ocean_region_lab = unique(ocean_region_lab),
                     ystart = Stock[1],
                     yend = Stock[length(Stock)])

# PDO stationary model
pdo.dot <- rbind(hb_param_df(stat_pdo, par=c("gamma"), region.var="Ocean.Region2"),
                  hb_param_df(stat_pdo, par=c("kappa"), region.var="Ocean.Region2"),
                  hb_param_df(stat_pdo, par=c("lambda"), region.var="Ocean.Region2"))
pdo.dot <- ocean_region_lab(pdo.dot, "region", FALSE)
pdo.dot$Stock <- factor(pdo.dot$Stock, levels = levels(data_master$Stock))
pdo.dot <- mutate(pdo.dot, var = case_when(par=="gamma" ~ "SST",
                                             par=="kappa" ~ "Competitors",
                                             par=="lambda" ~ "PDO"))
pdo.dot <- dplyr::arrange(pdo.dot, Stock)
pdo.mu <- plyr::ddply(pdo.dot, .(region, var), summarize,
                       mu_mean = unique(mu_mean),
                       mu_2.5 = unique(`mu_2.5%`),
                       mu_97.5 = unique(`mu_97.5%`),
                       ocean_region_lab = unique(ocean_region_lab),
                       ystart = Stock[1],
                       yend = Stock[length(Stock)])


# Dot and whiskers plot by stock
pdf(file=here(sens.fig.dir, "pdo-stat-dot.pdf"))
g <- ggplot(pdo.dot) +
  geom_vline(xintercept = 0, color = "grey50", linetype = 2, linewidth = 0.25) +
  geom_point(aes(x = mean, y = Stock, color = ocean_region_lab, shape = ocean_region_lab)) +
  geom_segment(aes(y = Stock, yend = Stock, x = `2.5%`, xend = `97.5%`,
                   color = ocean_region_lab), linewidth = 0.25) +
  geom_segment(data = pdo.mu, aes(y = ystart, yend = yend, x = mu_mean, xend = mu_mean,
                                   color = ocean_region_lab), linewidth = 0.25) +
  geom_segment(data=df.mu, aes(y = ystart, yend = yend, x = mu_mean, xend = mu_mean,
                               color = ocean_region_lab), linewidth = 0.25, col="grey50",
               alpha=0.6) + 
  geom_rect(data = pdo.mu, aes(xmin = mu_2.5, xmax = mu_97.5, ymin = ystart,
                                ymax = yend, fill = ocean_region_lab),
            alpha = 0.2) +
  col.scale.reg +
  scale_shape_manual(values = shp.reg) +
  scale_fill_manual(values = col.region, guide="none") +
  labs(x = "Coefficient",
       y = "",
       color = "",
       shape = "") +
  facet_wrap( ~ factor(var, levels=c("SST", "Competitors", "PDO"))) +
  scale_x_continuous(breaks=c(-.5, -.25,0,.25, .5), limits=c(-0.7, 0.7))+
  theme_sleek(base_size = 10) +
  theme(legend.justification = c(0, 0),
        legend.position = c(0.01, 0.87),
        legend.key.size = unit(10, "pt"),
        legend.background = element_blank(),
        legend.text = element_text(size = 8),
        panel.spacing.x = unit(-0.5, "pt"))
print(g)
dev.off()

# Era comparison plot with PDO

df.reg.pdo <- era_hb_param_df(era_pdo, par=c("gamma", "kappa", "lambda"), mu=T, neras=3)
df.reg.pdo$varnam[is.na(df.reg.pdo$varnam)] <- "PDO"

# load basic era model fits
load(here('output', 'models', 'dyn', speciesFlag, 'hbm_era_2c.RData'), verbose=T)
df.reg <- era_hb_param_df(era.2c, par=c("gamma", "kappa"), mu=T, neras=3, region.var="Ocean.Region2")

pdf(file=here(sens.fig.dir, "pdo-era-dot.pdf"))
g <- ggplot(ocean_region_lab(df.reg.pdo)) + 
  geom_hline(yintercept = 0, color = "grey50", linetype = 2, linewidth = 0.25) +
  geom_point(data=ocean_region_lab(df.reg), aes(y=reg_mean, x=era), col="grey50", shape=18, size=3, position=position_nudge(x=.1)) +
  geom_segment(data=ocean_region_lab(df.reg), aes(y=lower_10, yend=upper_90, x=era, xend=era), col="grey50",  
               linewidth=.8, position=position_nudge(x=.1)) + # plot base models for comparison
  geom_point(aes(y=reg_mean, x=era, col=ocean_region_lab), shape=18, size=3) +
  geom_segment(aes(y=lower_10, yend=upper_90, x=era, xend=era, col=ocean_region_lab), 
               linewidth=.8) +
  facet_grid(cols=vars(factor(varnam, levels=c("SST", "Competitors", "PDO"))), rows=vars(ocean_region_lab), scales="free_y") +
  scale_colour_manual(values=col.region, guide="none") +
  labs(x="Era", y="Covariate effect") +  theme_sleek() + lims(y=c(-1,1))
print (g)
dev.off()

