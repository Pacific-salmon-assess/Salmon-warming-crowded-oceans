## Stationary Bayesian Model inference
## Produces figures

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


# Set directory paths
fit.dir <- here("output", "models", "stat", speciesFlag) # place to store model fits
fig.dir <- here("Rmd", "figures") # place to store figures generated in this script

# Create destination folder if it doesn't exist
if(!dir.exists(fig.dir)) dir.create(fig.dir, recursive = T)


# Load model fits
load(here(fit.dir, "stat_a.Rdata"))

load(here(fit.dir, "single-stock", "single_stock_lms.Rdata"), verbose=T) # load single stock fits for comparison if we have them

# Define shape for figs
shp.reg <- c(18, 16, 17, 15)
names(shp.reg) <-  c("West Coast", "Southeast Alaska", "Gulf of Alaska", "Bering Sea")


gamma.stock <- hb_param_df(stat_a, "gamma", "Ocean.Region2", "SST", info=info_master)
kappa.stock <- hb_param_df(stat_a, "kappa", "Ocean.Region2", "Comp", info=info_master)
df.dot <- rbind(gamma.stock, kappa.stock )
df.dot <- ocean_region_lab(df.dot, "region", FALSE)
df.dot$Stock <- factor(df.dot$Stock, levels = levels(data_master$Stock))
df.dot$var <- factor(df.dot$var, levels = c("SST", "Comp" )) # ,"SST x Comp"))
df.mu <- plyr::ddply(df.dot, .(region, var), summarize,
                     mu_mean = unique(mu_mean),
                     mu_2.5 = unique(`mu_2.5%`),
                     mu_97.5 = unique(`mu_97.5%`),
                     ocean_region_lab = unique(ocean_region_lab),
                     ystart = Stock[1],
                     yend = Stock[length(Stock)])


## Table: Regional coefficients (hyperparameters) ----

gamma <- rstan::summary(stat_a, pars = "mu_gamma")$summary
kappa <- rstan::summary(stat_a, pars = "mu_kappa")$summary
reg <- c("West Coast", "Southeast Alaska", "Gulf of Alaska", "Bering Sea")
if(speciesFlag=="chum") reg <- reg[-2]

tab.g <- data.frame(species = speciesFlag,
                    reg = reg,
                    coef = "SST",
                    lower = gamma[ , "2.5%"],
                    mean = gamma[ , "mean"],
                    median = gamma[ , "50%"],
                    upper = gamma[ , "97.5%"])
tab.k <- data.frame(species = speciesFlag,
                    reg = reg,
                    coef = "Comp",
                    lower = kappa[ , "2.5%"],
                    mean = kappa[ , "mean"],
                    median = kappa[ , "50%"],
                    upper = kappa[ , "97.5%"])

tab.coef <- rbind(tab.g, tab.k)
tab.coef$perc <- (exp(tab.coef$mean) - 1) * 100
row.names(tab.coef) <- NULL
names(tab.coef) <- c("Species", "Ecosystem", "Coefficient", "Lower 95% CI", "Mean", "Median",
                     "Upper 95% CI", "Mean % change in R/S")

write.csv(tab.coef, file = here(fig.dir, paste0("reg_coefficients_stat_a_", speciesFlag, ".csv")))


# Coefficient table -- Population lvl
gamma <- rstan::summary(stat_a, pars = "gamma")$summary
kappa <- rstan::summary(stat_a, pars = "kappa")$summary

tab.g <- data.frame(species = speciesFlag,
                    stock = info_master$Stock,
                    reg = info_master$ocean_region_lab,
                    coef = "SST",
                    lower = gamma[ , "2.5%"],
                    mean = gamma[ , "mean"],
                    median = gamma[ , "50%"],
                    upper = gamma[ , "97.5%"])
tab.k <- data.frame(species = speciesFlag,
                    stock = info_master$Stock,
                    reg = info_master$ocean_region_lab,
                    coef = "Comp",
                    lower = kappa[ , "2.5%"],
                    mean = kappa[ , "mean"],
                    median = kappa[ , "50%"],
                    upper = kappa[ , "97.5%"])

tab.coef <- rbind(tab.g, tab.k)
tab.coef$perc <- (exp(tab.coef$mean) - 1) * 100
row.names(tab.coef) <- NULL
names(tab.coef) <- c("Species", "Stock", "Region", "Coefficient", "Lower 95% CI", "Mean", "Median",
                     "Upper 95% CI", "Mean % change in R/S")

write.csv(tab.coef, file = here(fig.dir, paste0("stk_coefficients_", "stat_a_", speciesFlag, ".csv")))



## --  Plot timeseries length (R/S)

prod_dat <- fill.time.series(data_master)
prod_dat <- ocean_region_lab(prod_dat)

g <- ggplot(prod_dat) +
  geom_vline(xintercept=c(1988,2011), color = "grey50", linetype = 2, linewidth = 0.25) +
  geom_line(data=na.omit(prod_dat), aes(x=BY, y=lnRS), col="grey75") +
  geom_line(aes(x=BY, y=lnRS, col=ocean_region_lab)) +
  facet_grid(rows=vars(Stock), switch ="y", scales="free_y", as.table=F) +
  scale_colour_manual(values=col.region) +
  xlim(1960,2020) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.y.left = element_text(angle=0, hjust=0, margin=margin(l=0, r=0)),
        strip.background = element_rect(fill="transparent", colour="transparent"),
        strip.text = element_text(size=7, ),
        panel.spacing.y = unit(0, unit="cm"),
        panel.background = element_rect(fill="white"),
        legend.position = "none") +
  labs(y="Population", x="Brood Year")


if(speciesFlag=="chum") {
  png(here(fig.dir, paste0(speciesFlag, "_ts_length.png")), height = 5.5, width = 6,  units = "in", res = 72*4) }
if(speciesFlag=="pink") {
  png(here(fig.dir, paste0(speciesFlag, "_ts_length.png")), height = 6, width = 6,  units = "in", res = 72*4) }

print(g)
dev.off()

if(speciesFlag=="sockeye") {
  prod_dat1 <- prod_dat |>
      group_by(Stock) %>%
      mutate(group_id = cur_group_id()) %>%
      ungroup()|>
    filter(group_id %in% c(1:33))

  prod_dat2 <- prod_dat |>
    group_by(Stock) %>%
    mutate(group_id = cur_group_id()) %>%
    ungroup()|>
    filter(group_id %in% c(34:66))

  g1 <- ggplot(prod_dat1) +
    geom_vline(xintercept=c(1988,2011), color = "grey50", linetype = 2, linewidth = 0.25) +
    geom_line(data=na.omit(prod_dat1), aes(x=BY, y=lnRS), col="grey75") +
    geom_line(aes(x=BY, y=lnRS, col=ocean_region_lab)) +
    facet_grid(rows=vars(Stock), switch ="y", scales="free_y", as.table=F) +
    scale_colour_manual(values=col.region) +
    xlim(1960,2020) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          strip.text.y.left = element_text(angle=0, hjust=0, margin=margin(l=0, r=0)),
          strip.background = element_rect(fill="transparent", colour="transparent"),
          strip.text = element_text(size=7, ),
          panel.spacing.y = unit(0, unit="cm"),
          panel.background = element_rect(fill="white"),
          legend.position = "none") +
    labs(y="Population", x="Brood Year")

  g2 <- ggplot(prod_dat2) +
    geom_vline(xintercept=c(1988,2011), color = "grey50", linetype = 2, linewidth = 0.25) +
    geom_line(data=na.omit(prod_dat2), aes(x=BY, y=lnRS), col="grey75") +
    geom_line(aes(x=BY, y=lnRS, col=ocean_region_lab)) +
    facet_grid(rows=vars(Stock), switch ="y", scales="free_y", as.table=F) +
    scale_colour_manual(values=col.region) +
    xlim(1960,2020) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          strip.text.y.left = element_text(angle=0, hjust=0, margin=margin(l=0, r=0)),
          strip.background = element_rect(fill="transparent", colour="transparent"),
          strip.text = element_text(size=7, ),
          panel.spacing.y = unit(0, unit="cm"),
          panel.background = element_rect(fill="white"),
          legend.position = "none") +
    labs(y="Population", x="Brood Year")

  png(here(fig.dir, paste0(speciesFlag, "2_ts_length.png")), height = 6, width = 6,  units = "in", res = 72*4)
  print(g1)
  dev.off()

  png(here(fig.dir, paste0(speciesFlag, "1_ts_length.png")), height = 6, width = 6,  units = "in", res = 72*4)
  print(g2)
  dev.off()
}

if(exists("ss.all.yrs") & speciesFlag == "sockeye") {

  ## Fig: Dot + density main with single stock estimates overlaid
  ss.dat <- ss.all.yrs$coef$model4a %>%
    dplyr::filter(variable %in% c("early_sst_stnd", "np_all_spp_sec_stnd")) %>%
    dplyr::mutate(var = ifelse(variable == "early_sst_stnd", "SST", "Comp"))
  ss.dat$Stock <- factor(ss.dat$Stock, levels=levels(data_master$Stock))
  ss.dat$var <- factor(ss.dat$var, levels=c("SST", "Comp"))
  df.dot.ss <- dplyr::left_join(df.dot, ss.dat, by=c("Stock", "var"))

  g <- ggplot(df.dot.ss) +
    geom_vline(xintercept = 0, color = "grey50", linetype = 2, linewidth = 0.25) +
    geom_point(aes(x = mean, y = Stock, color = ocean_region_lab, shape = ocean_region_lab, fill=ocean_region_lab)) +
    geom_segment(aes(y = Stock, yend = Stock, x = `2.5%`, xend = `97.5%`,
                     color = ocean_region_lab), linewidth = 0.25) +
    geom_segment(data = df.mu, aes(y = ystart, yend = yend, x = mu_mean, xend = mu_mean,
                                   color = ocean_region_lab), linewidth = 0.25) +
    geom_rect(data = df.mu, aes(xmin = mu_2.5, xmax = mu_97.5, ymin = ystart,
                                ymax = yend, fill = ocean_region_lab),
              alpha = 0.2) +
    geom_point(aes(x=value, y=Stock, colour=ocean_region_lab, shape=ocean_region_lab), fill="transparent") +
    col.scale.reg +
    scale_shape_manual(values = rev(c(23, 21, 24, 22)), guide = "legend") +
    scale_fill_manual(values = col.region, guide="none") +
    labs(x = "Coefficient",
         y = "Population",
         color = "",
         shape = "") +
    guides(shape = "none") +
    facet_wrap( ~ var) +
    scale_x_continuous(breaks=c(-0.5,0,0.5)) +
    coord_cartesian(xlim=c(-1,1), clip="off") +
    theme_sleek(base_size = 10) +
    theme(legend.justification = c(0, 0),
          legend.position = c(0.01, 0.87),
          legend.key.size = unit(10, "pt"),
          legend.background = element_blank(),
          legend.text = element_text(size = 8),
          panel.spacing.x = unit(-0.5, "pt"))

  png(here(fig.dir, paste0(speciesFlag, "_coef_catepillar_heir_indv.png")), height = 9, width = 7,  units = "in", res = 72*4)
  print(g)
  dev.off()

}

if(exists("ss.all.yrs") & speciesFlag == "chum") {

  ## Fig: Dot + density main with single stock estimates overlaid
  ss.dat <- ss.all.yrs$coef$model4a %>%
    dplyr::filter(variable %in% c("early_sst_stnd", "np_all_spp_sec_stnd")) %>%
    dplyr::mutate(var = ifelse(variable == "early_sst_stnd", "SST", "Comp"))
  ss.dat$Stock <- factor(ss.dat$Stock, levels=levels(data_master$Stock))
  ss.dat$var <- factor(ss.dat$var, levels=c("SST", "Comp"))
  df.dot.ss <- dplyr::left_join(df.dot, ss.dat, by=c("Stock", "var"))

  g <- ggplot(df.dot.ss) +
    geom_vline(xintercept = 0, color = "grey50", linetype = 2, linewidth = 0.25) +
    geom_point(aes(x = mean, y = Stock, color = ocean_region_lab, shape = ocean_region_lab, fill=ocean_region_lab)) +
    geom_segment(aes(y = Stock, yend = Stock, x = `2.5%`, xend = `97.5%`,
                     color = ocean_region_lab), linewidth = 0.25) +
    geom_segment(data = df.mu, aes(y = ystart, yend = yend, x = mu_mean, xend = mu_mean,
                                   color = ocean_region_lab), linewidth = 0.25) +
    geom_rect(data = df.mu, aes(xmin = mu_2.5, xmax = mu_97.5, ymin = ystart,
                                ymax = yend, fill = ocean_region_lab),
              alpha = 0.2) +
    geom_point(aes(x=value, y=Stock, colour=ocean_region_lab, shape=ocean_region_lab), fill="transparent") +
    col.scale.reg +
    scale_shape_manual(values = rev(c(23, 21, 24, 22)), guide = "legend") +
    scale_fill_manual(values = col.region, guide="none") +
    labs(x = "Coefficient",
         y = "Population",
         color = "",
         shape = "") +
    guides(shape = "none") +
    facet_wrap( ~ var) +
    scale_x_continuous(breaks=c(-0.5,0,0.5)) +
    coord_cartesian(xlim=c(-1,1), clip="off") +
    theme_sleek(base_size = 10) +
    theme(legend.justification = c(0, 0),
          legend.position = c(0.01, 0.8),
          legend.key.size = unit(10, "pt"),
          legend.background = element_blank(),
          legend.text = element_text(size = 8),
          panel.spacing.x = unit(-0.5, "pt"))

  png(here(fig.dir, paste0(speciesFlag, "_coef_catepillar_heir_indv.png")), height = 4.5, width = 7,  units = "in", res = 72*4)
  print(g)
  dev.off()

}

if(exists("ss.all.yrs") & speciesFlag == "pink") {

  ## Fig: Dot + density main with single stock estimates overlaid
  ss.dat <- ss.all.yrs$coef$model4a %>%
    dplyr::filter(variable %in% c("early_sst_stnd", "np_all_spp_sec_stnd")) %>%
    dplyr::mutate(var = ifelse(variable == "early_sst_stnd", "SST", "Comp"))
  ss.dat$Stock <- factor(ss.dat$Stock, levels=levels(data_master$Stock))
  ss.dat$var <- factor(ss.dat$var, levels=c("SST", "Comp"))
  df.dot.ss <- dplyr::left_join(df.dot, ss.dat, by=c("Stock", "var"))

  g <- ggplot(df.dot.ss) +
    geom_vline(xintercept = 0, color = "grey50", linetype = 2, linewidth = 0.25) +
    geom_point(aes(x = mean, y = Stock, color = ocean_region_lab, shape = ocean_region_lab, fill=ocean_region_lab)) +
    geom_segment(aes(y = Stock, yend = Stock, x = `2.5%`, xend = `97.5%`,
                     color = ocean_region_lab), linewidth = 0.25) +
    geom_segment(data = df.mu, aes(y = ystart, yend = yend, x = mu_mean, xend = mu_mean,
                                   color = ocean_region_lab), linewidth = 0.25) +
    geom_rect(data = df.mu, aes(xmin = mu_2.5, xmax = mu_97.5, ymin = ystart,
                                ymax = yend, fill = ocean_region_lab),
              alpha = 0.2) +
    geom_point(aes(x=value, y=Stock, colour=ocean_region_lab, shape=ocean_region_lab), fill="transparent") +
    col.scale.reg +
    scale_shape_manual(values = rev(c(23, 21, 24, 22)), guide = "legend") +
    scale_fill_manual(values = col.region, guide="none") +
    labs(x = "Coefficient",
         y = "Population",
         color = "",
         shape = "") +
    guides(shape = "none") +
    facet_wrap( ~ var) +
    scale_x_continuous(breaks=c(-0.5,0,0.5)) +
    coord_cartesian(xlim=c(-1,1), clip="off") +
    theme_sleek(base_size = 10) +
    theme(legend.justification = c(0, 0),
          legend.position = c(0.01, 0.77),
          legend.key.size = unit(10, "pt"),
          legend.background = element_blank(),
          legend.text = element_text(size = 8),
          panel.spacing.x = unit(-0.5, "pt"))

  png(here(fig.dir, paste0(speciesFlag, "_coef_catepillar_heir_indv.png")), height = 4.8, width = 7,  units = "in", res = 72*4)
  print(g)
  dev.off()

}
## --- Remove large model fits (saved in stat_hbm_fit)
rm(list = c("stat_a"))
