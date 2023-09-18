## (8) Recreate the tables and figures in the 2020 publication ##
## ------------------------------------------------------ ##

if(!dir.exists("./figures/stat/pub/")) dir.create("./figures/stat/pub/")


## Model fits with 3 Ocean Regions
# ---------------------------------------------------------------------

# Load fits if not just run
if(!exists("stat_a")){
  for(i in list.files(path = "./output/models/stat/", pattern = "*.RData$")) {
            load(paste("./output/models/stat/", i, sep = ""), verbose=T)
        }
}
load("./output/models/stat/single-stock/single_stock_lms.Rdata", verbose=T)

fit.list <- list(stat_a, stat_tr, stat_ctrl)
fitnam <- list("stat_a", "stat_tr", "stat_ctrl")

## Define colors
col.region <- chroma::qpal(7, luminance = 40)[c(7, 5, 3, 1)]
names(col.region) <- unique(sock.info$ocean_label2)
col.scale.reg <- scale_colour_manual(name = "Ocean Region", values=col.region)
col.region.3 <- chroma::qpal(7, luminance = 40)[c(1, 4, 6)]

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
  
  write.csv(tab.coef, file = paste0("./figures/stat/pub/model_coefficients_", fitnam[[n]],  ".csv"))
  
}

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

pdf("./figures/stat/pub/dens_stat_a.pdf", width = 4, height = 6)
print(g)
dev.off()



## Fig: dot + density main --------------------------------- 
gamma.stock <- hb_param_df(stat_a, "gamma", "Ocean.Region2", "SST")
kappa.stock <- hb_param_df(stat_a, "kappa", "Ocean.Region2", "Comp")
df.dot <- rbind(gamma.stock, kappa.stock ) # , chi.stock)
df.dot <- ocean_region_lab(df.dot, "region", FALSE)
df.dot$Stock <- factor(df.dot$Stock, levels = levels(sock$Stock))
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
  scale_shape_manual(values = c(15, 16, 17, 18)) +
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

pdf("./figures/stat/pub/coef_dot_stat_a.pdf", width = 6.5, height = 6.0)
print(g)
dev.off()

## Fig: Dot + density main with single stock estimates overlaid
ss.dat <- ss.all.yrs$coef$model4a %>% 
  dplyr::filter(variable %in% c("early_sst_stnd", "np_pinks_sec_stnd")) %>% 
  dplyr::mutate(var = ifelse(variable == "early_sst_stnd", "SST", "Comp"))
ss.dat$Stock <- factor(ss.dat$Stock, levels=levels(sock$Stock), ordered = T)
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
  scale_shape_manual(values = c(22, 21, 24, 23), guide = "legend") +
  scale_fill_manual(values = col.region, guide="none") +
  guides(shape = guide_legend(override.aes = list(shape=c(15, 16, 17, 18)))) +
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

pdf("./figures/stat/pub/coef_dot_ss_stat_a.pdf", width = 6.5, height = 6.0)
print(g)
dev.off()


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

pdf("./figures/stat/pub/dens_stat_tr.pdf", width = 4, height = 6)
print(g)
dev.off()



## Fig: dot + density main --------------------------------- 
gamma.stock <- hb_param_df(stat_tr, "gamma", "Ocean.Region2", "SST")
kappa.stock <- hb_param_df(stat_tr, "kappa", "Ocean.Region2", "Comp")
df.dot <- rbind(gamma.stock, kappa.stock ) # , chi.stock)
df.dot <- ocean_region_lab(df.dot, "region", FALSE)
df.dot$Stock <- factor(df.dot$Stock, levels = levels(sock$Stock))
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
  scale_shape_manual(values = c(15, 16, 17, 18)) +
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

pdf("./figures/stat/pub/coef_dot_stat_tr.pdf", width = 6.5, height = 6.0)
print(g)
dev.off()

## Fig: Dot + density main with single stock estimates overlaid
# N/A because single stock uses all yrs of data


## stat_ctrl -----------------------------------------------------

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

pdf("./figures/stat/pub/dens_stat_ctrl.pdf", width = 4, height = 6)
print(g)
dev.off()



## Fig: dot + density main --------------------------------- 
gamma.stock <- hb_param_df(stat_ctrl, "gamma", "Ocean.Region", "SST", info = sock.info[sock.info$Stock %in% ctrl_dat$Stock,])
kappa.stock <- hb_param_df(stat_ctrl, "kappa", "Ocean.Region", "Comp", info = sock.info[sock.info$Stock %in% ctrl_dat$Stock,])
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
  scale_shape_manual(values = c(15, 16, 17, 18)) +
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

pdf("./figures/stat/pub/coef_dot_stat_ctrl.pdf", width = 6.5, height = 6.0)
print(g)
dev.off()

## Fig: Dot + density main with single stock estimates overlaid
# N/A



## Fig: Map + covars ---------------------------------------

# 4 ocean regions
pdf("./figures/stat/pub/fig_map_covars.pdf",
    width = 5, height = 6.0)

## hi-res map
cl <- rnaturalearth::ne_states(country = c("United States of America",
                                           "Mexico",
                                           "Canada"))
namerica.state.sp <- sp::spTransform(cl, sp::CRS("EPSG:4326"))


## Define graphic layout
m <- rbind(c(1, 1, 1, 1, 1, 1, 1, 1),
           c(1, 1, 1, 1, 1, 1, 1, 1),
           c(1, 1, 1, 1, 1, 1, 1, 1),
           c(2, 2, 3, 3, 4, 4, 5, 5),
           c(6, 6, 7, 7, 8, 8, 9, 9),
           c(10, 10, 11, 11, 12, 12, 13, 13))
layout(m)
# print(m)
#layout.show(13)


## Set colors
col.lt <- rev(chroma::qpal(7)[c(1, 3, 5, 7)])
col.dk <- rev(chroma::qpal(7, luminance = 30)[c(1, 3, 5, 7)])
col.stk <- rev(chroma::qpal(7, alpha=0.4)[c(1, 3, 5, 7)])


oc.reg     <- unique(sock.info$Ocean.Region2)
rep.stocks <- c("Chilko", "Copper", "Wood")
p.title    <- c("West Coast", "Southeast Alaska", "Gulf of Alaska", "Bering Sea")
f.ind      <- 1

par(ps = 8,
    mar = c(0,0,0,0),
    oma = c(4, 4, 0.8, 0.5))


## ~~~~~~~~~~~~~~~~~~~~~~
## Map
## ~~~~~~~~~~~~~~~~~~~~~~

# Ignoring stock # labels for now

## Plot map
pch <- c(21, 22, 24, 25)
par(mar = c(5, 0, 0, 0))
plot(sock.info$lon, sock.info$lat,
     type = "n",
     xlab = "",
     ylab = "",
     xlim = c(-165, -120),
     ylim = c(47, 61),
     col = "red",
     axes = FALSE)
box(col = "grey65", lwd = 0.6)
axis(side = 1, lwd = 0, lwd.tick = 1, col = "grey65",
     padj = -1.9, tck = -0.02)
axis(side = 2, lwd = 0, lwd.tick = 1, las = 1, col = "grey65",
     hadj = 0.3, tck = -0.02)
mtext(expression(paste("Longitude (", degree, "E)")),
      side = 1, line = 1.8, cex = 1.0, col = "grey25")
mtext(expression(paste("Latitude (", degree, "N)")),
      side = 2, line = 1.8, cex = 1.0, col = "grey25")
sp::plot(namerica.state.sp, add = TRUE, border = "grey50", lwd = 0.5)
#add.label("a", font = 2, cex = 2, xfrac = -0.1, yfrac = 0, xpd = NA)


## Plot stock points
regs <- unique(sock.info$Ocean.Region2)
for(i in seq_along(regs)) {
  points(sock.info$lon[sock.info$Ocean.Region2 == regs[i]],
         sock.info$lat[sock.info$Ocean.Region2 == regs[i]],
         pch = pch[i], col = col.dk[i], bg = col.lt[i],
         cex = 1.5, lwd = 2)
}

#text(x = labs$x, y = labs$y, cex = 0.9, labels = labs$lab)
## text(x = labs$lon + labs$x.off, y = labs$lat + labs$y.off, cex = 0.9,
##      labels = labs$lab)

legend(x = -159, y = 53,
       legend = rev(p.title),
       pch = rev(pch), col = rev(col.dk), pt.bg = rev(col.lt),
       bty = "n", cex = 1.4)

## reset mar
par(mar = c(0, 0, 0, 0))
f.ind <- f.ind + 1


## ~~~~~~~~~~~~~~~~~~~~~~
## Smolt Age-1
## ~~~~~~~~~~~~~~~~~~~~~~
for(i in rev(seq_along(oc.reg))) {
  
  dat.i <- sock[sock$Ocean.Region2 == oc.reg[i], ]
  age.mean <- plyr::ddply(dat.i, .(BY), summarize,
                          mean = mean(ocean_1, na.rm = TRUE))
  plot(0, type = "n",
       ylim = c(0, 1),
       xlim = c(1950, 2015),
       ylab = "",
       xlab = "",
       axes = FALSE)
  mtext(p.title[i], col = "grey25")
  abline(v=1975,lty=2,col=" dark grey",lwd=1)
  if(i == 4) {
    mtext("Smolt age-1",
          side = 2, line = 2, cex = 1.0, col = "grey25")
    axis(side = 2, lwd = 0, lwd.tick = 1, col = "grey65", las = 1, hadj = 0.3, tck = -0.04)
    #add.label("b", font = 2, cex = 2, xfrac = -0.3, yfrac = -0.1, xpd = NA)
  }
  box(col = "grey65", lwd = 0.6)
  
  x <- lapply(split(dat.i, dat.i$Stock), function(x) {
    lines(x$BY, x$ocean_1, col = col.stk[i])
  })
  lines(age.mean$BY, age.mean$mean, col = col.dk[i], lwd = 2)
  f.ind <- f.ind + 1
}


## ~~~~~~~~~~~~~~~~~~~~~~
## Covariates
## ~~~~~~~~~~~~~~~~~~~~~~
covars <- c("np_pinks_sec", "early_sst")
scale <- c(100, 1)
at <- list(c(0, 2, 4, 6), -1:1)
ind <- 1
for(k in 1:2 ) {  ## loop over covars
  
  for(i in rev(seq_along(oc.reg))) {  ## loop over regions
    
    dat.i <- sock[sock$Ocean.Region2 == oc.reg[i], ]
    mm <- plyr::ddply(dat.i, .(BY), function(x) {
      min  <- min(x[ , covars[k]] / scale[k], na.rm = TRUE)
      max  <- max(x[ , covars[k]] / scale[k], na.rm = TRUE)
      mean <- mean(x[ , covars[k]] / scale[k], na.rm = TRUE)
      data.frame(min = min, max = max, mean = mean)
    })
    if(k == 1)
      ylim <- c(0, max(sock[ , covars[k]] / scale[k], na.rm = TRUE))
    if(k == 2)
      ylim <- range(sock[ , covars[k]] / scale[k], na.rm = TRUE)
    plot(0, type = "n",
         xlim = c(1950, 2015),
         ylim = ylim,
         axes = FALSE)
    abline(v=1975,lty=2,col="dark grey",lwd=1)
    if(i == 4) {
      axis(side = 2, lwd = 0, lwd.tick = 1, las = 1,
           col = "grey65", at = at[[k]], hadj = 0, tck = -0.04)
    }
    if(k == 2) {
      axis(side = 1, lwd = 0, lwd.tick = 1, las = 1, col = "grey65",
           at = c(1960, 1980, 2000), padj = -1.9, tck = -0.04)
    }
    if(ind == 6) {
      mtext("Brood year", side = 1, line = 1.8, cex = 1.0, col = "grey25")
    }
    if(ind == 1) {
      mtext("Pink index", side = 2, line = 2, cex = 1.0, col = "grey25")
    }
    if(ind == 5) {
      mtext("SST index", side = 2, line = 2, cex = 1.0, col = "grey25")
    }
    box(col = "grey65", lwd = 0.6)
    
    x <- lapply(split(dat.i, dat.i$Stock), function(x) {
      lines(x$BY, x[ , covars[k]] / scale[k], col = col.stk[i])
    })
    if(k == 1) {
      x <- lapply(split(dat.i, dat.i$Stock), function(x) {
        lines(x$BY, x[ , "np_pinks_sec"] / scale[k], col = col.stk[i], lty = 2)
      })
      ph <- plyr::ddply(dat.i, .(BY), function(x) {
        min  <- min(x[ , "np_pinks_sec"] / scale[k], na.rm = TRUE)
        max  <- max(x[ , "np_pinks_sec"] / scale[k], na.rm = TRUE)
        mean <- mean(x[ , "np_pinks_sec"] / scale[k], na.rm = TRUE)
        data.frame(min = min, max = max, mean = mean)
      })
      lines(ph$BY, ph$mean, col = col.dk[i], lwd = 1, lty = 2)
    }
    lines(mm$BY, mm$mean, col = col.dk[i], lwd = 1)
    ind <- ind + 1
    f.ind <- f.ind + 1
  }
}



## ~~~~~~~~~~~~~~~~~~~~~~
## Inset map
## ~~~~~~~~~~~~~~~~~~~~~~
# c(x_left, x_right, y_bottom, y_top)
par(fig=c(0.78, 0.99, 0.85, 0.99), mar = c(0,0,0,0), new = TRUE)
plot(1, 1,
     xlim = c(-180, -57),
     ylim = c(12.5, 80),
     col  = "white",
     xlab = "",
     ylab = "",
     axes = FALSE)
inc.cont <- c("United States of America",
              "Canada",
              "Mexico",
              "Guatemala",
              "Honduras",
              "Nicaragua",
              "Costa Rica",
              "Panama")
plot(countriesLow[countriesLow$ADMIN %in% inc.cont, ],
     col = "grey55", border = "grey55", add = TRUE)
rect(xleft = -165, ybottom = 47, xright = -120, ytop = 61,
     lwd = 2, border = 1)


dev.off()


## --- Remove large model fits (saved in stat_hbm_fit)
rm(list = c("stat_a", "stat_tr", "stat_ctrl"))
