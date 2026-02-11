
# Set paths to output locations - dependent on species
fit.dir <- here("output", "models", "dyn", speciesFlag) # where model fits are stored
fig.dir <- here("figures", "dyn", speciesFlag, "hbm_inf") # place to store figures generated in this script

# Create destination folder if it doesn't exist
if(!dir.exists(fig.dir))
  dir.create(fig.dir, recursive = T)


# Load fits
load(here(fit.dir, "hbm_dyn_2c.Rdata"))

# Set colours
col.region <- rev(chroma::qpal(7, luminance = 40)[c(1, 3, 5, 7)])
names(col.region) <- c("West Coast", "Southeast Alaska", "Gulf of Alaska", "Bering Sea")

col.eras <- c("#00b39e", "#b3a100", "#ff80d7", "#4db8ff",
              "#008070", "#6D6200FF", "#BC007FFF",  "#0070BDFF",
              "#00332d", "#4d4500", "#4d0034", "#002e4d")
names(col.eras) <- paste0(rep(names(col.region), 3), ".", rep(c("Early", "Middle", "Late"), each=4))

shp.reg <- c(18, 16, 17, 15)
names(shp.reg) <- c("West Coast", "Southeast Alaska", "Gulf of Alaska", "Bering Sea")



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
  g <- g + geom_line(aes(x=BY, y=reg_mean, linetype=even_odd, col=ocean_region_lab), linewidth=1) + geom_ribbon(aes(x=BY, y=reg_mean, ymin=lower_10, ymax=upper_90, fill=ocean_region_lab, group=even_odd), alpha=0.2)

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
  g <- g + geom_line(aes(x=BY, y=reg_mean, linetype=even_odd, col=ocean_region_lab), linewidth=1) + geom_ribbon(aes(x=BY, y=reg_mean, ymin=lower_10, ymax=upper_90, fill=ocean_region_lab, group=even_odd), alpha=0.2)
} else {
  g <- g + geom_line(aes(x=BY, y=reg_mean, col=ocean_region_lab), linewidth=1) + geom_ribbon(aes(x=BY, y=reg_mean, ymin=lower_10, ymax=upper_90, fill=ocean_region_lab), alpha=0.2)
}

pdf(here(fig.dir, "dyn_era_2c.pdf"))
print(g)
dev.off()
