## Bayesian model inference

if(!dir.exists("./figures/dyn/hbm_inf/"))
    dir.create("./figures/dyn/hbm_inf/")

## Load model fits
for(i in list.files(path = "./output/models/dyn/", pattern = "*.RData$")) {
    load(paste("./output/models/dyn/", i, sep = ""))
}

## Extract model outputs 
## ----------------------------------------------------------


## Eras model  ----------------------------------------

mcmc.era.sst  <- As.mcmc.list(era.sst, pars = pars.era)
mcmc.era.comp  <- As.mcmc.list(era.comp, pars = pars.era)
lst.era <- list(SST  = mcmc.era.sst,
                 Comp = mcmc.era.comp)


# Stock-specific dataframe

dfl.era.st <- lapply(seq_along(lst.era), function(i) {
    mcmc <- lst.era[[i]]
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
    df$var <- names(lst.era)[i]
    df$era <- factor(df$era, levels=c("early", "middle", "late"))
    df.summ <- dplyr::summarize(df, gamma_mu=mean(value), 
                                gamma_10=quantile(value, 0.1), 
                                gamma_90=quantile(value, 0.9), 
                                .by=c(Stock, par, ocean_region_lab, era))
    return(df.summ)
})
names(dfl.era.st) <- c("SST", "Comp")


# Summarized dataframe (across stocks)

dfl.era.reg <- lapply(seq_along(lst.era), function(i) {
  mcmc <- lst.era[[i]]
  df.wide <- coda_df(mcmc, grep("mu_gamma", coda::varnames(mcmc),
                                value = TRUE))
  df.long <- reshape2::melt(df.wide, id.vars = c("chain", "iter"),
                            variable.name = "par")
  info <- data.frame(par = as.character(unique(df.long$par)),
                     ocean_region = factor(rep(c("West Coast", "SEAK", "Gulf of Alaska", "Bering Sea"), 3), 
                                           levels=c("West Coast", "SEAK", "Gulf of Alaska", "Bering Sea") ), # This hard coding could be problematic 
                     era = factor( rep(c("early", "middle", "late"), each=4), 
                                   levels=c("early", "middle", "late") )
                    )
  df <- plyr::join(df.long, info, by = "par")
  df$var <- names(lst.era)[i]
  df.summ <- dplyr::summarize(df, gamma_mu=mean(value), 
                              gamma_10=quantile(value, 0.1), 
                              gamma_90=quantile(value, 0.9), 
                              .by=c(par, ocean_region, era))
  df.summ$ystart <- rep(c("Early Stuart", "Alastair", "Copper", "Nelson"), 3)
  df.summ$yend <- rep(c("Atnarko", "Chilkoot", "Chignik Lake", "Goodnews"), 3)
  return(df.summ)
})
names(dfl.era.reg) <- c("SST", "Comp")


# Density dataframe
post <- rstan::extract(era.sst, pars=c(paste0("gamma", c(1:3)), paste0("mu_gamma", c(1:3))))
dens.l <- lapply(post, function(x){
                        dens.out <- col_density(x, plot.it=F)
                        dens.list <- lapply(dens.out, adply, .margins=c(1,2))
                        dens.df <- join(dens.list$x, dens.list$y, by=c("X1", "X2"))
                        names(dens.df) <- c("n", "stock", "gamma.x", "gamma.dens")
                        dens.df <- mutate(dens.df, stock=levels(sock$Stock)[as.numeric(stock)] ) 
                        return(dens.df)
                        } )
dens.df <- bind_rows(dens.l, .id="par") # rename - this is confusing with the above naming
dens.df$era <- case_when(dens.df$par %in% c("gamma1", "mu_gamma1") ~ "early",
                         dens.df$par %in% c("gamma2", "mu_gamma2") ~ "middle",
                         dens.df$par %in% c("gamma3", "mu_gamma3") ~ "late") 
dens.df$era <- factor(dens.df$era, levels=c("early", "middle", "late"))
dens.df <- ocean_region_col(dens.df, stock.col = "stock")

# mean density dataframe # should summarize then apply col_density
dens.df.summ <- dplyr::summarize(dens.df, mu.x = mean(gamma.x), mu.dens = mean(gamma.dens), 
                          .by=c("era", "region", "n"))

## Graphics : Eras model
## ------------------------------------------------

# Set colours 
col.region <- chroma::qpal(7, luminance = 40)[c(1, 3, 5, 7)]
col.eras <- c("#4db8ff", "#b3a100", "#ff80d7", "#00b39e", 
              "#0070BDFF", "#6D6200FF", "#BC007FFF", "#008070",
              "#002e4d", "#4d4500", "#4d0034",  "#00332d")



# SST dot plot - Eras
pdf("./figures/dyn/hbm_inf/eras_coef_dot_sst_comb.pdf", width=5, height=8)
g <- ggplot(dfl.era.st$SST) +
  geom_vline(xintercept = 0, color = "grey50", linetype = 2, linewidth = 0.25) +
  geom_point(aes(x = gamma_mu, y = Stock, color = interaction(ocean_region_lab, era), shape = ocean_region_lab)) +
  #geom_segment(aes(y = Stock, yend = Stock, x = `2.5%`, xend = `97.5%`,
  #                 color = ocean_region_lab), linewidth = 0.25) +
  geom_segment(data = dfl.era.reg$SST, aes(y = ystart, yend = yend, x = gamma_mu, xend=gamma_mu,
                                 color = interaction(ocean_region, era)), linewidth = 0.25) +
  geom_rect(data = dfl.era.reg$SST, aes(xmin = gamma_10, xmax = gamma_90, ymin = ystart,
                         ymax = yend, fill = interaction(ocean_region, era)), alpha=0) +
  scale_color_manual(values = col.eras) +
  scale_shape_manual(values = c(15:18), guide = "legend") +
  scale_fill_manual(values = col.eras, guide="legend") +
  labs(x = "Coefficient",
       y = "",
       color = "",
       shape = "") +
  scale_x_continuous(breaks=c(-0.25,0,0.25))+
  theme_sleek(base_size = 10) +
  theme(legend.position = "none", 
        legend.justification = c(0, 0),
        legend.key.size = unit(10, "pt"),
        legend.background = element_blank(),
        legend.text = element_text(size = 8)
        )
print(g)
dev.off()


# SST Density plot - eras 
g <- ggplot(dens.df) + 
  geom_vline(xintercept = 0, color = "grey50", linetype = 2, linewidth = 0.25) +
  geom_path(aes(x=gamma.x, y=gamma.dens, group=stock, col=region), alpha=0.2) + 
  geom_path(data=dens.df.summ, aes(x=mu.x, y=mu.dens, col=region ), alpha=0.85, linewidth=1.5) +
  scale_colour_manual(values=col.region) +
  facet_wrap(~era, nrow=3) + coord_cartesian(xlim=c(-1, 1)) +
  theme_minimal() + labs(x="SST covariate effect", y="", col="Ocean Region") +
  theme(axis.text.y=element_blank())
pdf("./figures/dyn/hbm_inf/eras_sst_postdens.pdf")
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
pdf("./figures/dyn/hbm_inf/dyn_stacked_gamma_ts.pdf", width=5, height=10)
print(g)
dev.off()


## Dynamic model - Data ------------------------------------

lst.dyn <- list(SST  = dyn.sst,
                 Comp = dyn.comp)

dfl.dyn <- lapply(seq_along(lst.dyn), function(i) {
  probs <- c(0.025, 0.05, 0.10, 0.50, 0.90, 0.95, 0.975)
  summ <- summary(lst.dyn[[i]], pars = "gamma", probs = probs)[[1]]
  df <- data.frame(Stock = sock$Stock,
                   Ocean.Region = sock$Ocean.Region,
                   BY = sock$BY,
                   gamma = summ[ , "50%"],
                   lower = summ[ , "10%"],
                   upper = summ[ , "90%"],
                   # lower = summ[ , "5%"],
                   # upper = summ[ , "95%"],
                   var = names(lst.dyn)[i])
  return(df)
})

dfl.dyn <- lapply(dfl.dyn, ocean_region_lab)
names(dfl.dyn) <- c("SST", "Comp")


## Regional average timeseries : hbm 5

dyn.gamma.st <- ddply(hbm.gamma.diff, "Stock", moving_average_df, value="gamma")
dyn.gamma.avg <- dyn.gamma.st %>% 
  group_by(Ocean.Region, BY) %>%
  dplyr::summarize(mean_ma= mean(ma, na.rm=T), mean_ma_sd = sd(ma, na.rm=T))

g <- ggplot(dyn.gamma.avg) +
  geom_line(data=dyn.gamma.st, aes(x=BY, y=gamma, group=Stock, col=Ocean.Region), alpha=0.2) +
  geom_line(aes(x=BY, y=mean_ma, col=Ocean.Region), linewidth=1) +
  geom_ribbon(aes(x=BY, y=mean_ma, ymin=mean_ma-mean_ma_sd, ymax=mean_ma+mean_ma_sd, fill=Ocean.Region), alpha=0.2) +
  facet_wrap(vars(Ocean.Region), nrow=2) + ylim(c(-1,1)) + 
  scale_colour_manual(values=col.region, aesthetics=c("colour", "fill")) +
  theme_minimal()
pdf("./figures/dyn/hbm_inf/dyn_grouped_gamma.pdf")  
print(g)
dev.off()
  


## Dot-plots of alpha, beta, sigma 
plot_hbm_dot(era.sst, pdf.file = "./figures/dyn/hbm_inf/era_sst_dot.pdf")
plot_hbm_dot(dyn.sst, pdf.file = "./figures/dyn/hbm_inf/dyn_sst_dot.pdf")

## Posterior density of alpha
plot_hbm_dens(era.sst, pdf.file = "./figures/dyn/hbm_inf/era_sst_dens.pdf")
plot_hbm_dens(dyn.sst, pdf.file = "./figures/dyn/hbm_inf/dyn_sst_dens.pdf")


