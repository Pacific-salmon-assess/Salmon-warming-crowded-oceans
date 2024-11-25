## De-trend covariates

if(!dir.exists("./figures/detrend-covars/"))
    dir.create("./figures/detrend-covars/")

brood <- if(speciesFlag=="sockeye"){
    bt.complete} else if(speciesFlag=="pink") {
      bt.complete.pink } else if (speciesFlag=="chum"){
        bt.complete.chum  }


## Pink abundance ------------------------------------------

## Get pink abundance data
comp.dat <- read.table(file = "./data-downloaded/competitor_indices_2024.csv",
                       stringsAsFactors = FALSE, sep = ",",
                       header=TRUE)
comp.dat <- subset(comp.dat, select = c(Year, pink_numbers_np))
names(comp.dat) <- c("Year", "pink_raw")


## Fit time-trend models
pink.lm <- lm(pink_raw ~ Year, data = comp.dat)

pink.exp <- nls(pink_raw ~ I(exp(1)^(a + b * Year)),
               data = comp.dat, start = c(a = 0, b = 0.001), trace = FALSE)

pink.gam <- mgcv::gam(pink_raw ~ s(Year, k = 3), data = comp.dat, gamma = 2)


## Save de-trended data
comp.dat$pink_lm  <- residuals(pink.lm)
comp.dat$pink_exp <- residuals(pink.exp)
comp.dat$pink_gam <- residuals(pink.gam)
write.csv(comp.dat, "./data/pink_abundance_detrend.csv", row.names = FALSE)


## MSE
sum(comp.dat$pink_lm^2)
sum(comp.dat$pink_exp^2)
sum(comp.dat$pink_gam^2)


## Vis model fits
pdf("./figures/detrend-covars/pink_fits.pdf", width = 8, height = 9)
par(mfrow = c(2, 1))
plot(comp.dat$Year, comp.dat$pink_raw,
     type = "o",
     lty = 2,
     pch = 19,
     xlab = "Year",
     ylab = "Pink salmon abundance",
     axes = FALSE)
box(col = "grey50")
axis(side = 1, lwd = 0, lwd.tick = 1, col = "grey65")
axis(side = 2, lwd = 0, lwd.tick = 1, las = 1, col = "grey65")
lines(comp.dat$Year, predict(pink.lm), lwd = 2, col = "blue3")
lines(comp.dat$Year, predict(pink.exp),lwd = 2,  col = "red3")
lines(comp.dat$Year, predict(pink.gam),lwd = 2,  col = "green3")
legend("topleft", legend = c("data", "lm", "exp", "gam"), lwd = 2,
       col = c("black", "blue3", "red3", "green3"),
       bty = 'n')
##
plot(comp.dat$Year, resid(pink.lm),
     ylim = c(-200, 200),
     xlab = "Year",
     ylab = "Residual pink index",
     axes = FALSE,
     type = "o",
     pch = 19,
     col = "blue3",
     panel.first = abline(h = 0, col = "grey50", lty = 2))
box(col = "grey50")
axis(side = 1, lwd = 0, lwd.tick = 1, col = "grey65")
axis(side = 2, lwd = 0, lwd.tick = 1, las = 1, col = "grey65")
lines(comp.dat$Year, resid(pink.exp), type = "o", col = "red3", pch = 19)
lines(comp.dat$Year, resid(pink.gam), type = "o", col = "green3", pch = 19)
dev.off()


## Vis residual series correlations
pdf("./figures/detrend-covars/pink_cor.pdf", width = 8, height = 8)
g <- splom(comp.dat[ , names(comp.dat) != "Year"],
           par.settings = theme.mjm(),
           pch = 19, las = 1,
           panel = function(x, y, ...) {
               panel.splom(x, y, ...)
                panel.text(min(x), max(y),
                           paste("r =", round(cor(x, y, use = "pairwise.complete.obs"), 2)),
                           cex = 0.8, adj = 0, col = "black")
           })
print(g)
dev.off()


## Vis ACF of residual series
pdf("./figures/detrend-covars/pink_acf.pdf", width = 9, height = 8)
par(mfrow = c(2, 2), las = 1)
acf(comp.dat$pink_raw, main = "Pink data")
acf(resid(pink.lm), main = "Linear model residuals")
acf(resid(pink.exp), main = "Exponential model residuals")
acf(resid(pink.gam), main = "GAM model residuals")
dev.off()

# acf(comp.dat$pink_raw, main = "Pink data")$acf
# acf(resid(pink.lm), main = "Linear model residuals")$acf
# acf(resid(pink.exp), main = "Exponential model residuals")$acf
# acf(resid(pink.gam), main = "GAM model residuals")$acf



## SST 1st year --------------------------------------------
## Here, I use the raw SST values, rather than the anomalies

## Get SST data
sst.dat <- read.table("./data/sst_yr_1_stock_anomalies.csv",
                      stringsAsFactors = FALSE, sep = ",",
                      header = TRUE)
sst.dat <- sst.dat %>% filter(Stock.ID %in% info_master$Stock.ID)

## Fit models
sst.lst  <- vector("list", length(unique(sst.dat$Stock.ID)))
sst.pred <- vector("list", length(unique(sst.dat$Stock.ID)))
ind <- 1
for(i in unique(sst.dat$Stock.ID)) {
    sst.dat.i  <- sst.dat[sst.dat$Stock.ID == i, ]
    sst.pred.i <- sst.dat[sst.dat$Stock.ID == i, ]
    sst.lm <- lm(sst_raw ~ Year, data = sst.dat.i)
    sst.exp <- nls(sst_raw ~ I(exp(1)^(a + b * Year)),
                   control = list(maxiter = 500),
                   data = sst.dat.i, start = c(a = -1, b = 0.001), trace = FALSE)
    sst.gam <- mgcv::gam(sst_raw ~ s(Year, k = 3), data = sst.dat.i, gamma = 2)

    sst.pred.i$sst_lm  <- predict(sst.lm)
    sst.pred.i$sst_exp <- predict(sst.exp)
    sst.pred.i$sst_gam <- predict(sst.gam)
    sst.pred[[ind]] <- sst.pred.i

    sst.dat.i$sst_lm  <- residuals(sst.lm)
    sst.dat.i$sst_exp <- residuals(sst.exp)
    sst.dat.i$sst_gam <- residuals(sst.gam)
    sst.lst[[ind]] <- sst.dat.i
    ind <- ind + 1
}
sst.detrend <- plyr::rbind.fill(sst.lst) # detrended data = residuals of models
sst.predict <- plyr::rbind.fill(sst.pred)
write.csv(sst.detrend, "./data/sst_yr_1_stock_anomalies_detrend.csv", row.names = FALSE)


## Vis model fits
pdf("./figures/detrend-covars/sst_fits.pdf", width = 8, height = 9)
for(i in unique(sst.detrend$Stock.ID)) {
    sst.i <- sst.detrend[sst.detrend$Stock.ID == i, ]
    pred.i <- sst.predict[sst.predict$Stock.ID == i, ]
    par(mfrow = c(2, 1))
    plot(sst.i$Year, sst.i$sst_raw,
         type = "o",
         main = i,
         lty = 2,
         pch = 19,
         xlab = "Year",
         ylab = "SST raw",
         axes = FALSE,
         panel.first = abline(h = 0, col = "grey50", lty = 2))
    box(col = "grey50")
    axis(side = 1, lwd = 0, lwd.tick = 1, col = "grey65")
    axis(side = 2, lwd = 0, lwd.tick = 1, las = 1, col = "grey65")
    lines(pred.i$Year, pred.i$sst_lm, lwd = 2, col = "blue3")
    lines(pred.i$Year, pred.i$sst_exp, lwd = 2,  col = "red3")
    lines(pred.i$Year, pred.i$sst_gam, lwd = 2,  col = "green3")
    legend("topleft", legend = c("data", "lm", "exp", "gam"), lwd = 2,
           col = c("black", "blue3", "red3", "green3"),
           bty = 'n')
    ##
    plot(sst.i$Year, sst.i$sst_lm,
         ylim = c(-2, 2),
         xlab = "Year",
         ylab = "Residual sst index",
         axes = FALSE,
         type = "o",
         pch = 19,
         col = "blue3",
         panel.first = abline(h = 0, col = "grey50", lty = 2))
    box(col = "grey50")
    axis(side = 1, lwd = 0, lwd.tick = 1, col = "grey65")
    axis(side = 2, lwd = 0, lwd.tick = 1, las = 1, col = "grey65")
    lines(sst.i$Year, sst.i$sst_exp, type = "o", col = "red3", pch = 19)
    lines(sst.i$Year, sst.i$sst_gam, type = "o", col = "green3", pch = 19)
}
dev.off()


## Vis residual series correlations
pdf("./figures/detrend-covars/sst_cor.pdf", width = 8, height = 8)
for(i in unique(sst.detrend$Stock.ID)) {
    sst.i <- sst.detrend[sst.detrend$Stock.ID == i, ]
    g <- splom(sst.i[ , !names(sst.i) %in% c("Year", "sst_anomaly", "Stock.ID", "Species")],
               par.settings = theme.mjm(),
               main = as.character(i),
               pch = 19, las = 1,
               panel = function(x, y, ...) {
                   panel.splom(x, y, ...)
                    panel.text(min(x), max(y),
                               paste("r =", round(cor(x, y, use = "pairwise.complete.obs"), 2)),
                               cex = 0.8, adj = 0, col = "black")
               })
    print(g)
}
dev.off()


## Vis ACF of residual series
pdf("./figures/detrend-covars/sst_acf.pdf", width = 9, height = 8)
for(i in unique(sst.detrend$Stock.ID)) {
    sst.i <- sst.detrend[sst.detrend$Stock.ID == i, ]
    par(mfrow = c(2, 2), las = 1)
    acf(sst.i$sst_raw, main = "SST raw")
    acf(sst.i$sst_lm, main = "Linear model residuals")
    acf(sst.i$sst_exp, main = "Exponential model residuals")
    acf(sst.i$sst_gam, main = "GAM model residuals")
    title(as.character(i), line = -1.5, outer = TRUE)
}
dev.off()



#### Join covariate data ####

## SST during early marine life 
sst.dt <- clim.wgt.avg(brood.table = brood,
                          env.data = sst.detrend,
                          env.covar = "sst_lm",
                          type = "first_year",
                          out.covar = "early_sst_lm")


## Age weighted competitor index: Pinks in 2nd yr marine life

## competitors in second year of marine life
comp.dt <- pink.wgt.avg(brood.table = brood,
                            pink.data = comp.dat,
                            pink.covar = "pink_gam",
                            type = "second_year",
                            out.covar = "np_pinks_gam_sec")


## Merge datasets 
detrend <- dplyr::left_join(brood, sst.dt, by=c("BY","Stock.ID"))
detrend <- dplyr::left_join(detrend, comp.dt, by=c("BY","Stock.ID"))
detrend <- geographic.order(detrend) 


## Add derived columns
detrend <- ddply(detrend, .(Stock), transform,
                      RS = R/S,
                      RS_stnd = scale(R/S)[ , 1],
                      lnRS = log(R/S),
                      S_stnd = scale(S)[ , 1],
                      early_sst_lm_stnd = scale(early_sst_lm)[ , 1],
                      np_pinks_gam_sec_stnd = scale(np_pinks_gam_sec)[ , 1])


#### Fit HBMs with detrended data ####
## stationary model detrended ---------------------------------------

## Run MCMC
stan.dat <- stan_data_stat(detrend,
                          scale.x1 = TRUE,
                          var.x2 = "early_sst_lm_stnd", # use lm detrend for sst
                          var.x3 = "np_pinks_gam_sec_stnd", # use gam detrend for comp
                          var.region = "Ocean.Region2",
                          alpha.group = ifelse(speciesFlag=="sockeye", TRUE, FALSE)) 
stat_detrend <- rstan::stan(file = "./stan/hbm_stat_2c.stan",
                           data = stan.dat,
                           warmup = 1000,
                           iter = 3000,
                           cores = 4,
                           chains = 4,
                           thin = 2,
                           seed = 123,
                           control = list(adapt_delta = 0.90,
                                          max_treedepth = 10))
save(stat_detrend, file = "./output/models/stat_detrend.RData")




## eras model detrended -----------------------------------------------
## Run MCMC
stan.dat <- stan_data_dyn(detrend,
                          breakpoint2 = 2011,
                          scale.x1 = TRUE,
                          var.x2 = "early_sst_lm_stnd", # use lm detrend for sst
                          var.x3 = "np_pinks_gam_sec_stnd") # use gam detrend for comp
era_detrend <- rstan::stan(file = "./stan/hbm_era_2c.stan",
                     data = stan.dat,
                     warmup = 1000,
                     iter = 3000,
                     cores = 4,
                     chains = 4,
                     thin = 2,
                     seed = 123,
                     control = list(adapt_delta = 0.90,
                                    max_treedepth = 10))
save(era_detrend, file = "./output/models/era_detrend.RData")

# Density dataframes - detrended
#load(file = "./output/models/era_detrend.RData", verbose=T)
df.s.dt <- era_hb_param_df(era_detrend, par=c("gamma", "kappa"))
df.reg.dt <- era_hb_param_df(era_detrend, par=c("gamma", "kappa"), mu=T, neras=3)


# load *NOT DETRENDED* era model fits
load(file=here('output','models', 'dyn', speciesFlag, 'hbm_era_2c_3sub.RData', verbose=T)) 
# Eras df by stock
era.stk <- era_hb_param_df(era.2c.3sub, par=c("gamma", "kappa"))
era.reg <- era_hb_param_df(era.2c.3sub, par=c("gamma", "kappa"), mu = TRUE)


### --- Figures

## Dot-and-whiskers to compare detrended to 'regular'
ggplot(ocean_region_lab(df.reg.dt)) + 
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


## Density plot
ggplot(dens.df.s.dt) + 
  geom_vline(xintercept = 0, color = "grey50", linetype = 2, linewidth = 0.25) +
  geom_path(aes(x=x, y=dens, group=stock, col=ocean_region_lab), alpha=0.2) + 
  geom_path(data=dens.df.reg.dt, aes(x=x, y=dens, col=ocean_region_lab), alpha=0.85, linewidth=1) +
  scale_colour_manual(values=col.region) +
  facet_grid(rows=vars(era), cols=vars(varnam)) + 
  coord_cartesian(xlim=c(-1, 1)) +
  theme_minimal() + labs(x="covariate effect", y="", col="Ocean Region") +
  theme(axis.text.y=element_blank(), 
        legend.position="bottom")
pdf(here(fig.dir, "era_detrend_dens.pdf"))
print(g)
dev.off()

## Detrended time series
g1 <- ggplot(ocean_region_lab(detrend)) + 
    hline_0(lty="dashed", col="gray50", alpha=0.5) +
    geom_line(data=data_master, aes(x=BY, y=early_sst, group=Stock), col="orange", alpha=0.7) +
    geom_line(aes(x=BY, y=early_sst_lm, group=Stock), col="navyblue", alpha=0.7) +
    facet_wrap(vars(Ocean.Region2)) + labs(x="Brood Year", y="SST Index (anomalies)") +
    theme_sleek()

g2 <- ggplot(detrend) + 
  hline_0(lty="dashed", col="gray50", alpha=0.5) +
  geom_line(data=data_master, aes(x=BY, y=np_pinks_sec, group=Stock), col="orange", alpha=0.7) +
  geom_line(aes(x=BY, y=np_pinks_gam_sec, group=Stock), col="navyblue", alpha=.7) +
  theme_sleek() + labs(x="Brood Year", y="NP Pink Abundance")

png(filename=here(fig.dir, 'detrended_covars.png'), width=900*2, height=400*2, res=72*3)
cowplot::plot_grid(g1,g2)
dev.off()


