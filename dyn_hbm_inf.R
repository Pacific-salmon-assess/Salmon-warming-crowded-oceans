## Bayesian model inference

if(!dir.exists("./figures/hbm_inf/"))
    dir.create("./figures/hbm_inf/")

## Load model fits
for(i in list.files(path = "./output/hbm_fit/", pattern = "*.RData$")) {
    load(paste("./output/hbm_fit/", i, sep = ""))
}

## Extract model outputs 
## ----------------------------------------------------------


## Eras model : hbm3 ----------------------------------------

mcmc.hbm3.sst  <- As.mcmc.list(hbm3.sst,  pars = pars.hbm3)
lst.hbm3 <- list(SST  = mcmc.hbm3.sst)


# Summarized dataframe (across stocks)

dfl.hbm3 <- lapply(seq_along(lst.hbm3), function(i) {
    mcmc <- lst.hbm3[[i]]
    df.wide <- coda_df(mcmc, grep("mu_gamma", coda::varnames(mcmc),
                                  value = TRUE))
    df.long <- reshape2::melt(df.wide, id.vars = c("chain", "iter"),
                              variable.name = "par")
    info <- data.frame(par = as.character(unique(df.long$par)),
                       ocean_region = rep(c("WC", "GOA", "BS"), 3),
                       era = c(rep("early", 3),
                               rep("middle", 3),
                               rep("late", 3)))
    df <- plyr::join(df.long, info, by = "par")
    df$var <- names(lst.hbm3)[i]
    return(df)
})


# Stock-specific dataframe

dfl.hbm3.st <- lapply(seq_along(lst.hbm3), function(i) {
    mcmc <- lst.hbm3[[i]]
    df.wide <- coda_df(mcmc, grep("^gamma", coda::varnames(mcmc),
                                  value = TRUE))
    df.long <- reshape2::melt(df.wide, id.vars = c("chain", "iter"),
                              variable.name = "par")
    info <- sock.info[ , c("Stock.ID", "Stock", "Ocean.Region2")]
    info1 <- info; info2 <- info; info3 <- info
    info1$par <- paste0("gamma1[", 1:nlevels(info$Stock), "]")
    info2$par <- paste0("gamma2[", 1:nlevels(info$Stock), "]")
    info3$par <- paste0("gamma3[", 1:nlevels(info$Stock), "]")
    info1$era <- "early"
    info2$era <- "middle"
    info3$era <- "late"
    info <- rbind(info1, info2, info3)
    info <- ocean_region_lab(info, var="Ocean.Region2")
    df <- plyr::join(info, df.long, by = "par")
    df$var <- names(lst.hbm3)[i]
    return(df)
})


## Dynamic model: hbm5  ------------------------------------

lst.hbm5 <- list(SST  = hbm5.sst)
                 
dfl.hbm5 <- lapply(seq_along(lst.hbm5), function(i) {
    probs <- c(0.025, 0.05, 0.10, 0.50, 0.90, 0.95, 0.975)
    summ <- summary(lst.hbm5[[i]], pars = "gamma", probs = probs)[[1]]
    df <- data.frame(Stock = sock$Stock,
                     Ocean.Region = sock$Ocean.Region,
                     BY = sock$BY,
                     gamma = summ[ , "50%"],
                     lower = summ[ , "10%"],
                     upper = summ[ , "90%"],
                     # lower = summ[ , "5%"],
                     # upper = summ[ , "95%"],
                     var = names(lst.hbm5)[i])
    return(df)
})

df.hbm5 <- plyr::rbind.fill(dfl.hbm5)
df.hbm5 <- ocean_region_lab(df.hbm5)
df.hbm5 <- plyr::ddply(df.hbm5, .(ocean_region_lab, var, Stock), transform,
                       sig = (lower < 0 & upper < 0) |
                             (lower > 0 & upper > 0))
df.hbm5$var <- factor(df.hbm5$var, levels = c("SST")) 

hbm.gamma.diff <- df.hbm5
save(hbm.gamma.diff, file = "./output/hbm_gamma_diff.RData")



## Graphics 
## ------------------------------------------------

## Dot-plots of alpha, beta, sigma 
plot_hbm_dot(hbm3.sst, pdf.file = "./figures/dyn/hbm_inf/hbm3_sst_dot.pdf")
plot_hbm_dot(hbm5.sst, pdf.file = "./figures/dyn/hbm_inf/hbm5_sst_dot.pdf")

## Posterior density of alpha
plot_hbm_dens(hbm3.sst, pdf.file = "./figures/dyn/hbm_inf/hbm3_sst_dens.pdf")
plot_hbm_dens(hbm5.sst, pdf.file = "./figures/dyn/hbm_inf/hbm5_sst_dens.pdf")


## Dot-plots of covariate effects: Eras (hbm3)

# hb_param_df function copy with change allowing multiple eras: 
#TODO modify function code instead of copying here
  s.par <- as.data.frame(summary(hbm3.sst, pars = c("gamma1", "gamma2", "gamma3"), probs = c(0.025, 0.05, 0.10, 0.50, 0.90, 0.95, 0.975))[[1]])
  s.par$par <- rep(c("gamma1", "gamma2", "gamma3"), each=nlevels(sock$Stock))
  s.par$Stock <- sock.info$Stock
  s.par$Ocean.Region <- sock.info$Ocean.Region
  
  s.mu <- as.data.frame(summary(hbm3.sst, pars = paste0("mu_", c("gamma1", "gamma2", "gamma3")),
                                probs = c(0.025, 0.05, 0.10, 0.50, 0.90, 0.95, 0.975))[[1]])
  names(s.mu) <- paste0("mu_", names(s.mu))
  s.mu[["Ocean.Region"]] <- unique(sock.info[["Ocean.Region"]])
  s.mu$par <- rep(c("gamma1", "gamma2", "gamma3"), each=3)
  
  df.dot.era <- dplyr::left_join(s.par, s.mu, by = c('Ocean.Region', 'par'))

  df.dot.era <- ocean_region_lab(df.dot.era, factor=FALSE)
  df.dot.era$Stock <- factor(df.dot.era$Stock, levels = levels(sock$Stock)) 
  df.dot.era$era <- df.dot.era$par
#df.dot.era$var <- factor(df.dot.era$var, levels = c("SST")) #, "Comp" )) # ,"SST x Comp"))

df.mu.era <- plyr::ddply(df.dot.era, .(Ocean.Region), dplyr::summarize,
                     mu_mean = unique(mu_mean),
                     mu_2.5 = unique(`mu_2.5%`),
                     mu_97.5 = unique(`mu_97.5%`),
                     ocean_region_lab = unique(ocean_region_lab),
                     ystart = Stock[1],
                     yend = Stock[length(Stock)],
                     era = unique(par))

# Set colours 
col.region <- chroma::qpal(7, luminance = 40)[c(1, 3, 5)]

# plot and save
pdf("./figures/dyn/hbm_inf/eras_coef_dot.pdf")
g <- ggplot(df.dot.era) +
  geom_vline(xintercept = 0, color = "grey50", linetype = 2, linewidth = 0.25) +
  geom_point(aes(x = mean, y = Stock, color = ocean_region_lab, shape = ocean_region_lab, fill=ocean_region_lab)) +
  geom_segment(aes(y = Stock, yend = Stock, x = `2.5%`, xend = `97.5%`,
                   color = ocean_region_lab), linewidth = 0.25) +
  geom_segment(data = df.mu.era, aes(y = ystart, yend = yend, x = mu_mean, xend = mu_mean,
                                 color = ocean_region_lab), linewidth = 0.25) +
  geom_rect(data = df.mu.era, aes(xmin = mu_2.5, xmax = mu_97.5, ymin = ystart,
                              ymax = yend, fill = ocean_region_lab),
            alpha = 0.2) +
  scale_color_manual(values = rev(col.region)) +
  scale_shape_manual(values = c(15:17), guide = "legend") +
  scale_fill_manual(values = rev(col.region), guide="none") +
  guides(shape = guide_legend(override.aes = list(shape=c(15:17)))) +
  labs(x = "Coefficient",
       y = "",
       color = "",
       shape = "") +
  facet_wrap( ~ era, labeller = labeller(era = c("gamma1" = "< 1976", "gamma2" = "1977 - 1988", "gamma3" = "1989 - 2015"))) +
  scale_x_continuous(breaks=c(-0.25,0,0.25))+
  theme_sleek(base_size = 10) +
  theme(legend.justification = c(0, 0),
        legend.position = c(0.80, 0.52),
        legend.key.size = unit(10, "pt"),
        legend.background = element_blank(),
        legend.text = element_text(size = 8),
        panel.spacing.x = unit(-0.5, "pt"))
print(g)
dev.off()


## Stock-specific covariate effects: random walk # Stacked gamma timeseries

g <- ggplot(hbm.gamma.diff) + geom_line(aes(x=BY, y=gamma, col=Ocean.Region)) + 
  facet_grid(rows=vars(Stock), switch ="y", scales="free_y", as.table=F) + 
  scale_colour_manual(values=col.region) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.y.left = element_text(angle=0),
        strip.background = element_rect(fill="transparent", colour="transparent"),
        strip.text = element_text(size=7, ),
        panel.spacing.y = unit(0, unit="cm"),
        panel.background = element_rect(fill="white"),
        legend.position = "none") +
  labs(y="", x="Brood Year")
pdf("./figures/dyn/hbm_inf/hbm5_stacked_gamma_ts.pdf", width=5, height=10)
print(g)
dev.off()



## Regional average timeseries : hbm 5

hbm5.gamma.st <- ddply(hbm.gamma.diff, "Stock", moving_average_df, value="gamma")
hbm5.gamma.avg <- hbm5.gamma.st %>% 
  group_by(Ocean.Region, BY) %>%
  dplyr::summarize(mean_ma= mean(ma, na.rm=T), mean_ma_sd = sd(ma, na.rm=T))

g <- ggplot(hbm5.gamma.avg) +
  geom_line(data=hbm5.gamma.st, aes(x=BY, y=gamma, group=Stock, col=Ocean.Region), alpha=0.2) +
  geom_line(aes(x=BY, y=mean_ma, col=Ocean.Region), linewidth=1) +
  geom_ribbon(aes(x=BY, y=mean_ma, ymin=mean_ma-mean_ma_sd, ymax=mean_ma+mean_ma_sd, fill=Ocean.Region), alpha=0.2) +
  facet_wrap(vars(Ocean.Region), nrow=2) + ylim(c(-1,1)) + 
  scale_colour_manual(values=col.region, aesthetics=c("colour", "fill")) +
  theme_minimal()
pdf("./figures/dyn/hbm_inf/hbm5_grouped_gamma.pdf")  
print(g)
dev.off()
  
