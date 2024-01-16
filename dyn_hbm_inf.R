## Dynamic (time-varying) Bayesian model inference
## Produces figures


# Set paths to output locations - dependent on species 
fit.dir <- here("output", "models", "dyn", speciesFlag) # where model fits are stored
fig.dir <- here("figures", "dyn", speciesFlag, "hbm_inf") # place to store figures generated in this script

# Create destination folder if it doesn't exist
if(!dir.exists(fig.dir))
  dir.create(fig.dir, recursive = T)


# Load fits
if(!exists("era.2c")) {
    for(i in list.files(path = fit.dir, pattern = "*.RData$")) {
    load(here(fit.dir, i), verbose=T)
  }
}

# Set colours 
col.region <- rev(chroma::qpal(7, luminance = 40)[c(1, 3, 5, 7)])
names(col.region) <- unique(info_master$ocean_region_lab)

col.eras <- c("#00b39e", "#b3a100", "#ff80d7", "#4db8ff",
              "#008070", "#6D6200FF", "#BC007FFF",  "#0070BDFF",
              "#00332d", "#4d4500", "#4d0034", "#002e4d")
names(col.eras) <- paste0(rep(names(col.region), 3), ".", rep(c("Early", "Middle", "Late"), each=4))

shp.reg <- c(18, 16, 17, 15)
names(shp.reg) <- unique(info_master$ocean_region_lab)


## Two-covariate models ---------------------------------------------------

### --- Eras model: Data

# Stock-specific dataframe
df.era.st.2c <- era_hb_param_df(era.2c, par=c("gamma", "kappa")) 
df.era.st.2c <- ocean_region_lab(df.era.st.2c)

# Summarized dataframe (regional-level)
df.era.reg.2c <- era_hb_param_df(era.2c, par=c("gamma", "kappa"), mu = TRUE) 
df.era.reg.2c <- ocean_region_lab(df.era.reg.2c)

# Density dataframe - by stock
dens.df.st.2c <- era_density_df(era.2c, par=c("gamma", "kappa"))
dens.df.st.2c <- ocean_region_lab(dens.df.st.2c)
# Density dataframe (regional lvl)
dens.df.reg.2c <- era_density_df(era.2c, par=c("gamma", "kappa"), mu=T)
dens.df.reg.2c <- ocean_region_lab(dens.df.reg.2c)


### --- Eras model: Figures

# Caterpillar plot (combined) - 2-covar
g <- ggplot(df.era.st.2c) +
  geom_vline(xintercept = 0, color = "grey50", linetype = 2, linewidth = 0.25) +
  geom_point(aes(x = mu, y = Stock, color = interaction(ocean_region_lab, era), shape = ocean_region_lab)) +
  geom_segment(data = df.era.reg.2c, aes(y = ystart, yend = yend, x = reg_mean, xend=reg_mean,
                                           color = interaction(ocean_region_lab, era)), linewidth = 0.25) +
  geom_rect(data = df.era.reg.2c, aes(xmin = lower_10, xmax = upper_90, ymin = ystart,
                                        ymax = yend, fill = interaction(ocean_region_lab, era)), alpha=0) +
  facet_wrap(vars(varnam), ncol=2) +
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
pdf(here(fig.dir, "eras_coef_dot_2c_comb.pdf"), width=5, height=8)
print(g)
dev.off()

# Caterpillar plot (panel) - 2-covar
g <- ggplot(df.era.st.2c) +
  geom_vline(xintercept = 0, color = "grey50", linetype = 2, linewidth = 0.25) +
  geom_point(aes(x = mu, y = Stock, color = ocean_region_lab, shape = ocean_region_lab)) +
  geom_segment(aes(y = Stock, yend = Stock, x = lower_10, xend = upper_90,
                   color = ocean_region_lab), linewidth = 0.25) +
  geom_segment(data = df.era.reg.2c, aes(y = ystart, yend = yend, x = reg_mean, xend=reg_mean,
                                           color = ocean_region_lab), linewidth = 0.25) +
  geom_rect(data = df.era.reg.2c, aes(xmin = lower_10, xmax = upper_90, ymin = ystart,
                                        ymax = yend, fill = ocean_region_lab), alpha=0.2) +
  facet_grid(cols=vars(era), rows=vars(varnam)) +
  scale_color_manual(values = col.region) +
  scale_shape_manual(values = c(15:18), guide = "legend") +
  scale_fill_manual(values = col.region, guide="legend") +
  labs(x = "SST Coefficient",
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
pdf(here(fig.dir, "eras_coef_dot_2c_panel.pdf"), width=8, height=8)
print(g)
dev.off()

# Posterior density plot - 2-covar
g <- ggplot(dens.df.st.2c) + 
  geom_vline(xintercept = 0, color = "grey50", linetype = 2, linewidth = 0.25) +
  geom_path(aes(x=x, y=dens, group=stock, col=ocean_region_lab), alpha=0.2) + 
  geom_path(data=dens.df.reg.2c, aes(x=x, y=dens, col=ocean_region_lab), alpha=0.85, linewidth=1) +
  scale_colour_manual(values=col.region) +
  facet_grid(rows=vars(era), cols=vars(varnam)) + 
  coord_cartesian(xlim=c(-1, 1)) +
  theme_minimal() + labs(x="covariate effect", y="", col="Ocean Region") +
  theme(axis.text.y=element_blank())
pdf(here(fig.dir, "eras_2c_dens.pdf"))
print(g)
dev.off()



## NEW scatterplot by eras



### --- Dynamic model: Data

# Stock-specific dataframe
probs <- c(0.025, 0.05, 0.10, 0.50, 0.90, 0.95, 0.975)
summ <- rstan::summary(dyn.2c, pars = c("gamma", "kappa"), probs = probs)[[1]]
df.dyn.st.2c <- data.frame(Stock = data_master$Stock,
                   Ocean.Region2 = data_master$Ocean.Region2,
                   BY = data_master$BY,
                   mu = summ[, "mean"],
                   se = summ[, "se_mean"],
                   lower_10 = summ[, "10%"],
                   upper_90 = summ[ , "90%"],
                   var = str_extract(rownames(summ), "[a-z]+"),
                   varnam = case_when(grepl("^gamma", rownames(summ)) ~ "SST",
                                      grepl("^kappa", rownames(summ)) ~ "Competitors")
                   )
df.dyn.st.2c <- ocean_region_lab(df.dyn.st.2c)

# Summarized dataframe (regional-level)
# gamma/kappa are series-specific; no mu output. Summarize stocks instead
if(speciesFlag=="pink"){
df.dyn.reg.2c <- df.dyn.st.2c %>% 
  mutate(even_odd = ifelse(grepl("Even$", Stock), "Even", "Odd")) %>% 
  dplyr::summarize(reg_mean=mean(mu, na.rm=T), 
                                  n_stk=n_distinct(Stock),
                                  lower_10=quantile(mu, 0.1), 
                                  upper_90=quantile(mu, 0.9), 
                                  .by=c(Ocean.Region2, BY, varnam, even_odd))
} else {
  df.dyn.reg.2c <- dplyr::summarize(df.dyn.st.2c, 
                     reg_mean=mean(mu, na.rm=T), 
                     n_stk=n_distinct(Stock),
                     lower_10=quantile(mu, 0.1), 
                     upper_90=quantile(mu, 0.9), 
                     .by=c(Ocean.Region2, BY, varnam))
}
df.dyn.reg.2c <- ddply(df.dyn.reg.2c, .(Ocean.Region2), dplyr::filter, n_stk >= max(n_stk)*0.1) # remove years with less than 10% of stocks observed
df.dyn.reg.2c <- ocean_region_lab(df.dyn.reg.2c)

### --- Dynamic model: Figures 

# Grouped gamma timeseries : 2-covar
g <- ggplot(df.dyn.reg.2c) +
  geom_line(data= df.dyn.st.2c, aes(x=BY, y=mu, group=Stock, col=ocean_region_lab), alpha=0.2) +
  facet_grid(rows=vars(ocean_region_lab), cols=vars(varnam)) + 
  ylim(c(-1,1)) + 
  scale_colour_manual(values=col.region, aesthetics=c("colour", "fill")) +
  theme_minimal()

if(speciesFlag=="pink"){
  g <- g + geom_line(aes(x=BY, y=reg_mean, group=even_odd, col=ocean_region_lab), linewidth=1) + geom_ribbon(aes(x=BY, y=reg_mean, ymin=lower_10, ymax=upper_90, fill=ocean_region_lab, group=even_odd), alpha=0.2)

} else {
  g <- g + geom_line(aes(x=BY, y=reg_mean, col=ocean_region_lab), linewidth=1) + geom_ribbon(aes(x=BY, y=reg_mean, ymin=lower_10, ymax=upper_90, fill=ocean_region_lab), alpha=0.2)
}

pdf(here(fig.dir, "dyn_2c_grouped.pdf"))
print(g)
dev.off()

# Grouped gamma with Eras coefficients overlaid
df.era.reg.2c <- df.era.reg.2c %>% mutate(xstart = case_when(era=="Early" ~ 1950,
                                                             era=="Middle" ~ 1979,
                                                             era=="Late" ~ 1989),
                                          xend = case_when(era=="Early" ~ 1978,
                                                           era=="Middle" ~ 1988,
                                                           era=="Late" ~ 2019))
g <- ggplot(df.dyn.reg.2c) +
  geom_line(data= df.dyn.st.2c, aes(x=BY, y=mu, group=Stock, col=ocean_region_lab), alpha=0.2) +
  geom_segment(data=df.era.reg.2c, aes(x=xstart, xend=xend, y=reg_mean, yend=reg_mean), col="gray40", linetype="dashed") +
  facet_grid(rows=vars(ocean_region_lab), cols=vars(varnam)) + 
  ylim(c(-1,1)) + 
  scale_colour_manual(values=rev(col.region), aesthetics=c("colour", "fill")) +
  scale_y_continuous(limits=c(-.75,.75), breaks=c(-0.5, 0, 0.5), oob=scales::squish) +
  theme_sleek() +
  theme(legend.position = "none") + labs(x="Brood Year", y="Mean covariate effects")

if(speciesFlag=="pink"){
  g <- g + geom_line(aes(x=BY, y=reg_mean, group=even_odd, col=ocean_region_lab), linewidth=1) + geom_ribbon(aes(x=BY, y=reg_mean, ymin=lower_10, ymax=upper_90, fill=ocean_region_lab, group=even_odd), alpha=0.2)
} else {
  g <- g + geom_line(aes(x=BY, y=reg_mean, col=ocean_region_lab), linewidth=1) + + geom_ribbon(aes(x=BY, y=reg_mean, ymin=lower_10, ymax=upper_90, fill=ocean_region_lab), alpha=0.2)
}

pdf(here(fig.dir, "dyn_era_2c.pdf"))  
print(g)
dev.off()



