## Dynamic (time-varying) Bayesian model inference
## Produces figures


# Set paths to output locations - dependent on species
fit.dir <- here("output", "models", "dyn", speciesFlag) # where model fits are stored
fig.dir <- here("figures", "dyn", speciesFlag, "hbm_inf") # place to store figures generated in this script

# Create destination folder if it doesn't exist
if(!dir.exists(fig.dir))
  dir.create(fig.dir, recursive = T)


# Load fits
load(here(fit.dir, "hbm_era_2c.Rdata"))


# Set colours for era plots
# CHANGED: added a 4th ("Stationary") colour block, one per region, since
# kappa is now stationary and its rows carry era = "Stationary" instead of
# Early/Middle/Late. Without these, points/segments for kappa would have no
# matching colour in the combined caterpillar plot below.
col.eras <- c("#00b39e", "#b3a100", "#ff80d7", "#4db8ff",
              "#008070", "#6D6200FF", "#BC007FFF",  "#0070BDFF",
              "#00332d", "#4d4500", "#4d0034", "#002e4d",
              "#595959", "#8c8c8c", "#404040", "#bfbfbf")   # CHANGED: colours for the "Stationary" era
names(col.eras) <- paste0(rep(names(col.region), 4), ".", rep(c("Early", "Middle", "Late", "Stationary"), each=4)) # CHANGED: 4 eras now, not 3


## Save coefficients in tables ------------------------------------------------

# Stock-specific dataframe
df.stk.out <- era_hb_param_df(era.2c, par=c("gamma", "kappa"), lower_CI=2.5, upper_CI=97.5)
df.stk.out <- df.stk.out |> mutate(Species = speciesFlag) |> ocean_region_lab() |>
  rename(Ecosystem = ocean_region_lab, Coefficient = varnam, Era = era, 'Lower 95% CI' = lower,
         Mean = mu, Median = med, 'Upper 95% CI' = upper) |>
  select(Species, Stock, Ecosystem, Coefficient, Era, 'Lower 95% CI', Mean, Median, 'Upper 95% CI')

write.csv(df.stk.out, file = here(fig.dir, paste0("stk_coefficients_era_", speciesFlag, ".csv")), row.names = FALSE)

# Summarized dataframe (regional-level)
df.reg.out <- era_hb_param_df(era.2c, par=c("gamma", "kappa"), mu = TRUE, lower_CI=2.5, upper_CI=97.5)
df.reg.out <- df.reg.out |> mutate(Species = speciesFlag) |> ocean_region_lab() |>
  rename(Ecosystem = ocean_region_lab, Coefficient = varnam, Era = era, 'Lower 95% CI' = lower,
         Mean = reg_mean, Median = med, 'Upper 95% CI' = upper) |>
  select(Species, Ecosystem, Coefficient, Era, 'Lower 95% CI', Mean, Median, 'Upper 95% CI')

write.csv(df.reg.out, file = here(fig.dir, paste0("reg_coefficients_era_", speciesFlag, ".csv")), row.names = FALSE)


## Figures ------------------------------------------------------------

## Get dataframes for plots

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



# Caterpillar plot (combined) - 2-covar
g <- ggplot(df.era.st.2c) +
  geom_vline(xintercept = 0, color = "grey50", linetype = 2, linewidth = 0.25) +
  geom_point(aes(x = mu, y = Stock, color = interaction(ocean_region_lab, era), shape = ocean_region_lab)) +
  geom_segment(data = df.era.reg.2c, aes(y = ystart, yend = yend, x = reg_mean, xend=reg_mean,
                                           color = interaction(ocean_region_lab, era)), linewidth = 0.25) +
  geom_rect(data = df.era.reg.2c, aes(xmin = lower, xmax = upper, ymin = ystart,
                                        ymax = yend, fill = interaction(ocean_region_lab, era)), alpha=0) +
  facet_wrap(vars(varnam), ncol=2) +
  scale_color_manual(values = col.eras) +
  scale_shape_manual(values = c(15:18), guide = "legend") +
  scale_fill_manual(values = col.eras, guide="legend") +
  labs(x = "Coefficient",
       y = "",
       color = "",
       shape = "") +
  scale_x_continuous(breaks=c(-0.5,0,0.5))+
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
  geom_segment(aes(y = Stock, yend = Stock, x = lower, xend = upper,
                   color = ocean_region_lab), linewidth = 0.25) +
  geom_segment(data = df.era.reg.2c, aes(y = ystart, yend = yend, x = reg_mean, xend=reg_mean,
                                           color = ocean_region_lab), linewidth = 0.25) +
  geom_rect(data = df.era.reg.2c, aes(xmin = lower, xmax = upper, ymin = ystart,
                                        ymax = yend, fill = ocean_region_lab), alpha=0.2) +
  facet_grid(rows=vars(era), cols=vars(varnam)) +
  scale_color_manual(values = col.region) +
  scale_shape_manual(values = c(15:18), guide = "legend") +
  scale_fill_manual(values = col.region, guide="legend") +
  labs(x = "SST Coefficient",
       y = "",
       color = "",
       shape = "") +
  scale_x_continuous(breaks=c(-0.5,0,0.5))+
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
  theme(axis.text.y=element_blank(),
        legend.position = "bottom")
pdf(here(fig.dir, "eras_2c_dens.pdf"))
print(g)
dev.off()


# Remove large stan fits
rm(list = c("era.2c"))

