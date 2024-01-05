## Stationary Bayesian Model inference
## Produces figures

# Set directory paths
fit.dir <- here("output", "models", "stat", speciesFlag) # place to store model fits
fig.dir <- here("figures", "stat", speciesFlag, "hbm_inf") # place to store figures generated in this script

# Create destination folder if it doesn't exist
if(!dir.exists(fig.dir)) dir.create(fig.dir, recursive = T)


# Load model fits if not just run
if(!exists("stat_a")){
  for(i in list.files(path = fit.dir, pattern = "*.RData$")) {
    load(here(fit.dir, i), verbose=T)
        }
}

load(here(fit.dir, "single-stock", "single_stock_lms.Rdata"), verbose=T) # add back in when we have these fits


## Define colors
col.region <- rev(chroma::qpal(7, luminance = 40)[c(1, 3, 5, 7)])
names(col.region) <- unique(info_master$ocean_region_lab)
col.scale.reg <- scale_colour_manual(name = "Ocean Region", values=col.region)
col.region.3 <- chroma::qpal(7, luminance = 40)[c(1, 4, 6)]
# Define shape
shp.reg <- c(18, 16, 17, 15)
names(shp.reg) <- unique(info_master$ocean_region_lab)


# Coefficient table
fitnam <- strsplit(list.files(path = fit.dir, pattern = "*.RData$"), ".RData")
fit.list <- if("stat_ctrl" %in% fitnam) list(stat_a, stat_tr, stat_ctrl) else list(stat_a, stat_tr) 

for (n in 1:length(fit.list)){
  
  ## Table: coefficients ----
  
  gamma <- rstan::summary(fit.list[[n]], pars = "mu_gamma")$summary
  kappa <- rstan::summary(fit.list[[n]], pars = "mu_kappa")$summary
  reg   <- if(fitnam[[n]] == "stat_ctrl") c("West Coast", "Gulf of Alaska", "Bering Sea") else
    c("West Coast", "Southeast Alaska", "Gulf of Alaska", "Bering Sea")
  
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
  
  tab.coef <- rbind(tab.g, tab.k) # add tab.c if exists
  tab.coef$perc <- (exp(tab.coef$mean) - 1) * 100
  row.names(tab.coef) <- NULL
  names(tab.coef) <- c("Ecosystem", "Coefficient", "Lower 95% CI", "Mean",
                       "Upper 95% CI", "Mean % change in R/S")
  
  write.csv(tab.coef, file = here::here(fig.dir, paste0("model_coefficients_", fitnam[[n]], ".csv")))
  
}

## Plot timeseries length (R/S) 
prod_dat <- ocean_region_lab(data_master)
g <- ggplot(prod_dat) + 
  geom_vline(xintercept=c(1976,1988), color = "grey50", linetype = 2, linewidth = 0.25) +
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

pdf(here(fig.dir, "ts_length.pdf"), width=4, height=6)
print(g)
dev.off()

## stat_a -----------------------------------------------------

## Fig: Posterior percent change density ------------------- 
lst <- hb05_density_df(stat_a, ocean.regions = 4)
s.df <- lst$stock
m.df <- lst$region
m.df$region <- factor(m.df$region, levels = c("West Coast", "Gulf of Alaska", "Southeast Alaska", "Bering Sea"))

## Covariate labels
vars <- data.frame(var = levels(m.df$var))
vars$lab <- paste0("(", letters[1:nrow(vars)], ") ", vars$var)
vars$var <- factor(vars$var, levels = c("SST", "Comp", "SST + Comp"))

g <- ggplot(m.df) +
  geom_vline(xintercept = 0, color = "grey50", linetype = 2, linewidth = 0.25) +
  geom_path(data = s.df[s.df$region == "West Coast", ],
            aes(x = x, y = y, group = stock), color = col.region[["West Coast"]], alpha=0.3,
            na.rm = TRUE) +
  geom_path(data = s.df[s.df$region == "Gulf of Alaska", ],
            aes(x = x, y = y, group = stock), color = col.region[["Gulf of Alaska"]], alpha=0.3,
            na.rm = TRUE) +
  geom_path(data = s.df[s.df$region == "Southeast Alaska", ],
            aes(x = x, y = y, group = stock), color = col.region[["Southeast Alaska"]], alpha=0.3,
            na.rm = TRUE) +
  geom_path(data = s.df[s.df$region == "Bering Sea", ],
            aes(x = x, y = y, group = stock), color = col.region[["Bering Sea"]], alpha=0.3,
            na.rm = TRUE) +
  geom_path(aes(x = x, y = y, color = region), linewidth = 1, alpha=1, 
            na.rm = TRUE) +
  col.scale.reg +
  labs(x = "Percent change in R/S",
       y = "Posterior density",
       color = "") +
  scale_x_continuous(limits = c(-100, 100), expand = c(0, 0)) +
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

pdf(here::here(fig.dir, "dens_stat_a.pdf"), width = 4, height = 6)
print(g)
dev.off()



## Fig: dot + density main --------------------------------- 
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

g <- ggplot(df.dot) +
  geom_vline(xintercept = 0, color = "grey50", linetype = 2, linewidth = 0.25) +
  geom_point(aes(x = mean, y = Stock, color = ocean_region_lab, shape = ocean_region_lab)) +
  geom_segment(aes(y = Stock, yend = Stock, x = `2.5%`, xend = `97.5%`,
                   color = ocean_region_lab), linewidth = 0.25) +
  geom_segment(data = df.mu, aes(y = ystart, yend = yend, x = mu_mean, xend = mu_mean,
                                 color = ocean_region_lab), linewidth = 0.25) +
  geom_rect(data = df.mu, aes(xmin = mu_2.5, xmax = mu_97.5, ymin = ystart,
                              ymax = yend, fill = ocean_region_lab),
            alpha = 0.2) +
  col.scale.reg +
  scale_shape_manual(values = shp.reg) +
  scale_fill_manual(values = col.region, guide="none") +
  labs(x = "Coefficient",
       y = "",
       color = "",
       shape = "") +
  facet_wrap( ~ var) +
  scale_x_continuous(breaks=c(-0.25,0,0.25), limits=c(-0.7, 0.7))+
  theme_sleek(base_size = 10) +
  theme(legend.justification = c(0, 0),
        legend.position = c(0.01, 0.87),
        legend.key.size = unit(10, "pt"),
        legend.background = element_blank(),
        legend.text = element_text(size = 8),
        panel.spacing.x = unit(-0.5, "pt"))

pdf(here::here(fig.dir, "coef_dot_stat_a.pdf"), width = 6.5, height = 6.0)
print(g)
dev.off()

if(exists("ss.all.yrs")) {
  
  ## Fig: Dot + density main with single stock estimates overlaid
  ss.dat <- ss.all.yrs$coef$model4a %>% 
    dplyr::filter(variable %in% c("early_sst_stnd", "np_pinks_sec_stnd")) %>% 
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
         y = "",
         color = "",
         shape = "") +
    facet_wrap( ~ var) +
    scale_x_continuous(breaks=c(-0.25,0,0.25), limits=c(-0.7, 0.7))+
    theme_sleek(base_size = 10) +
    theme(legend.justification = c(0, 0),
          legend.position = c(0.01, 0.87),
          legend.key.size = unit(10, "pt"),
          legend.background = element_blank(),
          legend.text = element_text(size = 8),
          panel.spacing.x = unit(-0.5, "pt"))
  
  pdf(here::here(fig.dir, "coef_dot_ss_stat_a.pdf"), width = 6.5, height = 6.0)
  print(g)
  dev.off()
  
}

## stat_tr -----------------------------------------------------

## Fig: Posterior percent change density ------------------- 
lst <- hb05_density_df(stat_tr, ocean.regions = 4)
s.df <- lst$stock
m.df <- lst$region
m.df$region <- factor(m.df$region, levels = c("West Coast", "Gulf of Alaska", "Southeast Alaska", "Bering Sea"))

## Covariate labels
vars <- data.frame(var = levels(m.df$var))
vars$lab <- paste0("(", letters[1:nrow(vars)], ") ", vars$var)
vars$var <- factor(vars$var, levels = c("SST", "Comp", "SST + Comp"))

g <- ggplot(m.df) +
  geom_vline(xintercept = 0, color = "grey50", linetype = 2, linewidth = 0.25) +
  geom_path(data = s.df[s.df$region == "West Coast", ],
            aes(x = x, y = y, group = stock), color = col.region["West Coast"], alpha=0.3,
            na.rm = TRUE) +
  geom_path(data = s.df[s.df$region == "Gulf of Alaska", ],
            aes(x = x, y = y, group = stock), color = col.region["Gulf of Alaska"], alpha=0.3,
            na.rm = TRUE) +
  geom_path(data = s.df[s.df$region == "Southeast Alaska", ],
            aes(x = x, y = y, group = stock), color = col.region["Southeast Alaska"], alpha=0.3,
            na.rm = TRUE) +
  geom_path(data = s.df[s.df$region == "Bering Sea", ],
            aes(x = x, y = y, group = stock), color = col.region["Bering Sea"], alpha=0.3,
            na.rm = TRUE) +
  geom_path(aes(x = x, y = y, color = region), linewidth = 1, alpha=1, 
            na.rm = TRUE) +
  col.scale.reg +
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

pdf(here::here(fig.dir, "dens_stat_tr.pdf"), width = 4, height = 6)
print(g)
dev.off()



## Fig: dot + density main --------------------------------- 
gamma.stock <- hb_param_df(stat_tr, "gamma", "Ocean.Region2", "SST")
kappa.stock <- hb_param_df(stat_tr, "kappa", "Ocean.Region2", "Comp")
df.dot <- rbind(gamma.stock, kappa.stock ) # , chi.stock)
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

g <- ggplot(df.dot) +
  geom_vline(xintercept = 0, color = "grey50", linetype = 2, linewidth = 0.25) +
  geom_point(aes(x = mean, y = Stock, color = ocean_region_lab, shape = ocean_region_lab)) +
  geom_segment(aes(y = Stock, yend = Stock, x = `2.5%`, xend = `97.5%`,
                   color = ocean_region_lab), linewidth = 0.25) +
  geom_segment(data = df.mu, aes(y = ystart, yend = yend, x = mu_mean, xend = mu_mean,
                                 color = ocean_region_lab), linewidth = 0.25) +
  geom_rect(data = df.mu, aes(xmin = mu_2.5, xmax = mu_97.5, ymin = ystart,
                              ymax = yend, fill = ocean_region_lab),
            alpha = 0.2) +
  col.scale.reg +
  scale_shape_manual(values = shp.reg) +
  scale_fill_manual(values = col.region, guide="none") +
  labs(x = "Coefficient",
       y = "",
       color = "",
       shape = "") +
  facet_wrap( ~ var) +
  scale_x_continuous(breaks=c(-0.25,0,0.25))+
  theme_sleek(base_size = 10) +
  theme(legend.justification = c(0, 0),
        legend.position = c(0.01, 0.87),
        legend.key.size = unit(10, "pt"),
        legend.background = element_blank(),
        legend.text = element_text(size = 8),
        panel.spacing.x = unit(-0.5, "pt"))

pdf(here::here(fig.dir, "coef_dot_stat_tr.pdf"), width = 6.5, height = 6.0)
print(g)
dev.off()



## stat_ctrl -----------------------------------------------------
  if("stat_ctrl" %in% fitnam){
  ## Fig: Posterior percent change density ------------------- 
  lst <- hb05_density_df(stat_ctrl, ocean.regions = 3)
  s.df <- lst$stock
  m.df <- lst$region
  m.df$region <- factor(m.df$region, levels = c("West Coast", "Gulf of Alaska", "Bering Sea"))
  
  ## Covariate labels
  vars <- data.frame(var = levels(m.df$var))
  vars$lab <- paste0("(", letters[1:nrow(vars)], ") ", vars$var)
  vars$var <- factor(vars$var, levels = c("SST", "Comp", "SST + Comp"))
  
  g <- ggplot(m.df) +
    geom_vline(xintercept = 0, color = "grey50", linetype = 2, linewidth = 0.25) +
    geom_path(data = s.df[s.df$region == "West Coast", ],
              aes(x = x, y = y, group = stock), color = col.region.3[1], alpha=0.3,
              na.rm = TRUE) +
    geom_path(data = s.df[s.df$region == "Gulf of Alaska", ],
              aes(x = x, y = y, group = stock), color = col.region.3[2], alpha=0.3,
              na.rm = TRUE) +
    geom_path(data = s.df[s.df$region == "Bering Sea", ],
              aes(x = x, y = y, group = stock), color = col.region.3[3], alpha=0.3,
              na.rm = TRUE) +
    geom_path(aes(x = x, y = y, color = region), linewidth = 1, alpha=1, 
              na.rm = TRUE) +
    scale_color_manual(values = col.region.3) +
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
  
  pdf(here::here(fig.dir, "dens_stat_ctrl.pdf"), width = 4, height = 6)
  print(g)
  dev.off()
  
  
  
  ## Fig: dot + density main --------------------------------- 
  gamma.stock <- hb_param_df(stat_ctrl, "gamma", "Ocean.Region", "SST", info = info_master[info_master$Stock %in% ctrl_dat$Stock,])
  kappa.stock <- hb_param_df(stat_ctrl, "kappa", "Ocean.Region", "Comp", info = info_master[info_master$Stock %in% ctrl_dat$Stock,])
  df.dot <- rbind(gamma.stock, kappa.stock ) # , chi.stock)
  df.dot <- ocean_region_lab(df.dot, "region", FALSE)
  df.dot$Stock <- factor(df.dot$Stock, levels = unique(ctrl_dat$Stock))
  df.dot$var <- factor(df.dot$var, levels = c("SST", "Comp" )) # ,"SST x Comp"))
  df.mu <- plyr::ddply(df.dot, .(region, var), summarize,
                       mu_mean = unique(mu_mean),
                       mu_2.5 = unique(`mu_2.5%`),
                       mu_97.5 = unique(`mu_97.5%`),
                       ocean_region_lab = unique(ocean_region_lab),
                       ystart = Stock[1],
                       yend = Stock[length(Stock)])
  
  g <- ggplot(df.dot) +
    geom_vline(xintercept = 0, color = "grey50", linetype = 2, size = 0.25) +
    geom_point(aes(x = mean, y = Stock, color = ocean_region_lab, shape = ocean_region_lab)) +
    geom_segment(aes(y = Stock, yend = Stock, x = `2.5%`, xend = `97.5%`,
                     color = ocean_region_lab), linewidth = 0.25) +
    geom_segment(data = df.mu, aes(y = ystart, yend = yend, x = mu_mean, xend = mu_mean,
                                   color = ocean_region_lab), linewidth = 0.25) +
    geom_rect(data = df.mu, aes(xmin = mu_2.5, xmax = mu_97.5, ymin = ystart,
                                ymax = yend, fill = ocean_region_lab),
              alpha = 0.2) +
    scale_color_manual(values = col.region.3) +
    scale_shape_manual(values = shp.reg) +
    scale_fill_manual(values = col.region.3, guide="none") +
    labs(x = "Coefficient",
         y = "",
         color = "",
         shape = "") +
    facet_wrap( ~ var) +
    scale_x_continuous(breaks=c(-0.25,0,0.25))+
    theme_sleek(base_size = 10) +
    theme(legend.justification = c(0, 0),
          legend.position = c(0.01, 0.87),
          legend.key.size = unit(10, "pt"),
          legend.background = element_blank(),
          legend.text = element_text(size = 8),
          panel.spacing.x = unit(-0.5, "pt"))
  
  pdf(here::here(fig.dir, "coef_dot_stat_ctrl.pdf"), width = 6.5, height = 6.0)
  print(g)
  dev.off()
  
  ## Fig: Dot + density main with single stock estimates overlaid
  # N/A

} # end if (stat_ctrl) loop


## Fig: Map + covars (Hannah's version) --------------------

# Downlad map and convert to sp
cl <- rnaturalearth::ne_states(country = c("United States of America", "Canada"))
na_map <- sf::st_as_sf(cl)

axes <- list( xlims=c(-165.5, -121), 
              ylims=c(47, 61),
              xbreaks=seq(-160,-120,10), 
              xlabels=as.character(seq(-160,-120,10)),
              seq(45, 65, 5), 
              ybreaks=seq(50, 60, 5),
              ylabels=as.character(seq(50,60,5)))

#make map data
map.info <- info_master %>% select(Stock, lon, lat, ocean_region_lab) %>% 
  mutate(stock.no = 1:nrow(info_master)) %>% 
  dplyr::summarize(n.stk = n_distinct(Stock),
                   first.stk = first(stock.no),
                   last.stk = last(stock.no),
                   .by= c("lat", "lon", "ocean_region_lab")) #%>%
#filter(!(first.stk %in% c(54:56))) # remove some overlapping points
map.info[map.info$first.stk %in% c(52:56, 41:45), c("first.stk", "last.stk")] <- NA
map.info <- map.info %>% mutate(num = ifelse(first.stk == last.stk, first.stk, paste(first.stk, last.stk, sep="-")))
misc.lab <- data.frame(x=c(-157.9, -152.5), y=c(59.5, 57.25), label=c("52-56", "41-45"))

map <- ggplot(map.info) + 
  geom_sf(data=na_map, color="grey30", fill="white", linewidth=0.1, ) + 
  ggspatial::geom_spatial_point(aes(x=lon, y=lat, col=ocean_region_lab, shape=ocean_region_lab, fill=ocean_region_lab), 
                                crs=4326, size=1, stroke=1.75, alpha=0.7) +
  #geom_text(aes(x=lon, y=lat, label=num), vjust=1.4, col="gray20", size=3) +
  ggrepel::geom_text_repel(aes(x=lon, y=lat, label=num), col="gray20", size=2.5, min.segment.length = 0.25, box.padding=0.1) +
  geom_text(data=misc.lab, aes(x, y, label=label), col="gray20", size=2.5) +
  coord_sf(xlim=axes$xlims, ylim=axes$ylims) +
  scale_x_continuous(breaks=axes$xbreaks, labels=axes$xlabels) +
  scale_y_continuous(breaks=axes$ybreaks, labels=axes$ylabels) +
  scale_colour_manual(values=col.region, name="Ocean Region Grouping") + 
  scale_fill_manual(values=col.dk, name="Ocean Region Grouping") + 
  scale_shape_manual(values=c(22, 24, 21, 23), name="Ocean Region Grouping") +
  labs(x="Longitude (°E)", y="Latitude (°N)") +
  theme_sleek() + 
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust=0.5),
        legend.position = c(0.25,0.25)
  )

## Make dataframe of covariate info
covar.dat.st <- data_master %>% select(Stock, BY, Ocean.Region2, early_sst, np_pinks_sec) %>% tidyr::pivot_longer(cols=c(early_sst, np_pinks_sec), names_to = "covar") %>% mutate(covar_nam = ifelse(covar=="early_sst", "SST Index", "Competitor Index"))
covar.dat.st <- ocean_region_lab(covar.dat.st)
covar.dat.reg <- dplyr::summarize(.data=covar.dat.st, mean_covar = mean(value), .by = c("BY", "ocean_region_lab", "covar_nam"))

covar <- ggplot(covar.dat.st) + 
  geom_vline(xintercept=c(1976,1988), color = "grey50", linetype = 2, linewidth = 0.25) +
  geom_line(aes(x=BY, y=value, group=Stock, colour=ocean_region_lab), linewidth=0.5, alpha=0.4) +
  #geom_line(data=covar.dat.reg, aes(x=BY, y=mean_covar, col=ocean_region_lab), linewidth=0.5) +
  facet_grid(rows=vars(covar_nam), cols=vars(ocean_region_lab), scales="free_y") +
  scale_colour_manual(values=col.region) +
  scale_y_continuous(n.breaks=4) +
  theme_sleek() + theme(legend.position="none") +
  labs(x= "Brood Year", y="")

cowplot::plot_grid(map, covar, nrow=2, rel_heights=c(3,2), rel_widths = c(2,2), labels="auto")



## --- Remove large model fits (saved in stat_hbm_fit)
rm(list = c("stat_a", "stat_tr"))
if(exists("stat_ctrl")) rm("stat_ctrl")
