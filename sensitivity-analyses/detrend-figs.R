# plot of detrended north pacific competitor index
library("mgcv")
library("cowplot")

# 1. Fit gam to pink time series
comp_gam <- gam(all_spp_numbers_np  ~ s(Year), data = raw.comp)
raw.comp$comp_detrend_np <- residuals(comp_gam)

# 2. Create a prediction grid (e.g., 200 points for a smooth line)
p_data <- data.frame(Year = seq(min(raw.comp$Year), max(raw.comp$Year), length.out = 200))

# 3. Predict the trend and standard errors
preds <- predict(comp_gam, newdata = p_data, se.fit = TRUE)
p_data$fit <- preds$fit
p_data$se  <- preds$se.fit

# 4. Plot using ggplot2
a <- ggplot() +
  geom_line(data = raw.comp, aes(x = Year, y = all_spp_numbers_np), alpha = 0.4) +
  geom_line(data = p_data, aes(x = Year, y = fit), color = "blue", size = 1) +
  ylab("Salmon abundance (millions)") +
  theme_sleek(base_size = 9)

b <- ggplot() +
  geom_line(data = raw.comp, aes(x = Year, y = comp_detrend_np), alpha = 0.4) +
  ylab("Detrended salmon abundance (millions)") +
  theme_sleek(base_size = 9) +
  scale_y_continuous(position = "right") # Moves axis and label to the right



png(here('figures/manuscript/SI/detrend-comp.png'), width=800*3, height=500*3, res=72*5)
plot_grid(a, b, labels = "auto", nrow=1)
dev.off()

## All spp. interaction density plot -------------------------------------------
s.df <- m.df <- NULL
for(sp in c("sockeye", "pink", "chum")){
  print(sp)
  # Species
  if(sp=="pink") {
    data_tmp <- pink
    info_tmp <- pink.info} else if (sp=="chum") {
      data_tmp <- chum
      info_tmp <- chum.info } else if(sp=="sockeye"){
        data_tmp <- sock
        info_tmp <- sock.info }
  load(here("output", "models", "stat", sp, "stat_inter_detr.RData"), verbose=T)
  lst <- hb07_density_df(stat_inter.detp, data=data_tmp, info = info_tmp,
                         ocean.regions = ifelse(sp=="chum", 3, 4))
  lst$stock$species <- sp
  lst$region$species <- sp
  s.df <- rbind(s.df, lst$stock)
  m.df <- rbind(m.df, lst$region)
  m.df$region <- factor(m.df$region, levels = c("West Coast", "Gulf of Alaska", "Southeast Alaska", "Bering Sea"))
}

## Covariate and species labels
m.df <- filter(m.df, var != "SST + Comp")
s.df <- filter(s.df, var != "SST + Comp")
vars <- data.frame(var = unique(m.df$var))
vars$lab <- paste0("(", letters[1:nrow(vars)], ") ", vars$var)
vars$var <- factor(vars$var, levels = c("SST", "Comp", "SST x Comp"))
m.df$species <- factor(str_to_sentence(m.df$species), levels=c("Sockeye","Chum","Pink"))
s.df$species <- factor(str_to_sentence(s.df$species), levels=c("Sockeye","Chum","Pink"))

m.df.plot <- m.df
s.df.plot <- s.df

g<-ggplot(m.df.plot) +
  geom_vline(xintercept = 0, color = "grey50", linetype = 2, linewidth = 0.25) +
  geom_path(data=s.df.plot, aes(x=x, y=y, group=stock, color = region), alpha=0.3, na.rm=T) +
  geom_path(aes(x = x, y = y, color = region), linewidth = 1, alpha=1,
            na.rm = TRUE) +
  scale_colour_manual(values=col.region) +
  labs(x = "Percent change in R/S",
       y = "Posterior density",
       color = "") +
  scale_y_continuous(breaks=NULL) +
  facet_grid(rows=vars(var), cols=vars(species)) +
  xlim(-50,50) +
  theme_sleek(base_size = 9) +
  theme(legend.justification = c(0, 0),
        legend.position = c(0.83, 0.76),
        legend.key.size = unit(10, "pt"),
        legend.background = element_blank(),
        legend.text = element_text(size = 6),
        panel.spacing.y = unit(-0.5, "pt"),
        strip.background = element_blank(),
        strip.text.x.top = element_text(size=9),
        legend.key.spacing = unit(-0.5, "pt"))


png(here('figures', 'manuscript/SI/', "detrend-posterior-percent-change.png"), width = 800*2, height = 500*2, res = 72*4)
print(g)
dev.off()
