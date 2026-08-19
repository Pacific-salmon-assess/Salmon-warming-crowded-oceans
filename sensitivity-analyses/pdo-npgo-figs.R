## Plots and inference for NPGO sensitivity analysis

sens.fig.dir <- here('sensitivity-analyses', 'pdo-npgo')
#### Exploratory plots ####

# Explore correlation between PDO and other covariates
# By region and by time period, for ALL species together

# Load data and combine species
sock.br.pdo.npgo <- read.csv(file=here("sensitivity-analyses", "fits", "sockeye", "data_covar_pdo_npgo.csv")) # sockeye
chum.br.pdo.npgo <- read.csv(file=here("sensitivity-analyses", "fits", "chum", "data_covar_pdo_npgo.csv"))
 # chum
pink.br.pdo.npgo <- read.csv(file=here("sensitivity-analyses", "fits", "pink", "data_covar_pdo_npgo.csv"))
 # pink
all.br.pdo.npgo <- bind_rows(sock.br.pdo.npgo, chum.br.pdo.npgo, pink.br.pdo.npgo)

all_spp_info <- bind_rows(sock.info, chum.info, pink.info)


pdf(file=here(sens.fig.dir, "pdo-npgo-covar-corr.pdf"))

## create empty 3D array
array.cor <- array(NA, dim = c(4,4,nrow(all_spp_info)))

## calculate stock specific covar correlations
for(i in seq_along(all_spp_info$Stock.ID)) {
  stk.i <- all.br.pdo.npgo[all.br.pdo.npgo$Stock.ID == all_spp_info$Stock.ID[i], ]
  covar.i <- subset(stk.i, select = c("early_sst", "np_all_spp_sec", "pdo_index", "npgo_index"))
  cor.i <- cor(covar.i, use = "pairwise.complete.obs")
  array.cor[ , , i] <- cor.i
}

## Average across stocks
cor.covars <- apply(array.cor, c(1, 2), mean)
row.names(cor.covars) <- c("early_sst", "np_all_spp_sec", "pdo_index", "npgo_index")
colnames(cor.covars) <- c("early_sst", "np_all_spp_sec", "pdo_index", "npgo_index")

corrplot::corrplot(cor.covars, method="color", type="upper",
                   addCoef.col = "black", tl.col="black",
                   diag=F,
                   col=rev(corrplot::COL2('RdBu', 200)),
                   title="Average stock-specific covariate correlations (All years)",
                   mar=c(0,0,2,0))

# Split by era
era.yrs <- list(start=c(1950,1989,2011),
                end=c(1988,2010,2019))
for(e in 1:3){

  ## create empty 3D array
  array.cor <- array(NA, dim = c(4,4,nrow(all_spp_info)))
  ## calculate stock specific covar correlations
  for(i in seq_along(all_spp_info$Stock.ID)) {
    stk.i <- filter(all.br.pdo.npgo, Stock.ID==all_spp_info$Stock.ID[i], BY %in% era.yrs$start[e]:era.yrs$end[e])
    if(empty(stk.i)) next
    covar.i <- subset(stk.i, select = c("early_sst", "np_all_spp_sec", "pdo_index", "npgo_index"))
    cor.i <- cor(covar.i, use = "pairwise.complete.obs")
    array.cor[ , , i] <- cor.i
  }

  ## Average across stocks
  cor.covars <- apply(array.cor, c(1, 2), mean, na.rm=T)
  row.names(cor.covars) <- c("early_sst", "np_all_spp_sec", "pdo_index", "npgo_index")
  colnames(cor.covars) <- c("early_sst", "np_all_spp_sec", "pdo_index", "npgo_index")

  g <- corrplot::corrplot(cor.covars, method="color", type="upper",
                          addCoef.col = "black", tl.col="black",
                          diag=F,
                          col=rev(corrplot::COL2('RdBu', 200)),
                          title=paste("Average stock-specific covariate correlations",
                                      era.yrs$start[e], "-", era.yrs$end[e]),
                          mar=c(0,0,2,0))
  print(g)
}
dev.off()


## Explore correlation between PDO, NPGO and lnRS
# For all species, by region and time period
era.yrs <- list(start=c(1950,1989,2011),
                end=c(1988,2010,2019))

pdf(file= here(sens.fig.dir, "pdo-npgo-prod-corr.pdf"))

cor.stock <- plyr::ddply(all.br.pdo.npgo, .(Stock.ID), plyr::summarize,
                         Ocean.Region2 = unique(Ocean.Region2),
                         early_sst = cor(lnRS, early_sst, use = "pairwise.complete.obs"),
                         np_all_spp_sec = cor(lnRS, np_all_spp_sec, use = "pairwise.complete.obs"),
                         npgo = cor(lnRS, npgo_index, use="pairwise.complete.obs"),
                         pdo = cor(lnRS, pdo_index, use="pairwise.complete.obs"))

cor.stock$Stock.ID <- NULL
cor.stock <- reshape2::melt(cor.stock, id.vars = "Ocean.Region2")

cor.ocean <- plyr::ddply(cor.stock, .(Ocean.Region2, variable), summarize,
                         cor.avg = mean(value))

g <- ggplot(ocean_region_lab(cor.ocean)) +
  geom_col(aes(x=variable, y=cor.avg, fill=ocean_region_lab), position=position_dodge(), alpha=0.7) +
  scale_fill_manual(values=col.region) +
  theme_sleek() + labs(x="", y="Average Correlation", fill="Region",
                       title="All Years")
print(g)

for(e in 1:3) { # loop over eras
  dat <-all.br.pdo.npgo %>% filter(BY %in% era.yrs$start[e]:era.yrs$end[e])
  cor.stock <- plyr::ddply(dat, .(Stock.ID), plyr::summarize,
                           Ocean.Region2 = unique(Ocean.Region2),
                           early_sst = cor(lnRS, early_sst, use = "pairwise.complete.obs"),
                           np_all_spp_sec = cor(lnRS, np_all_spp_sec, use = "pairwise.complete.obs"),
                           npgo = cor(lnRS, npgo_index, use="pairwise.complete.obs"),
                           pdo = cor(lnRS, pdo_index, use="pairwise.complete.obs"))

  cor.stock$Stock.ID <- NULL
  cor.stock <- reshape2::melt(cor.stock, id.vars = "Ocean.Region2")

  cor.ocean <- plyr::ddply(cor.stock, .(Ocean.Region2, variable), summarize,
                           cor.avg = mean(value, na.rm=T))
  g <- ggplot(ocean_region_lab(cor.ocean)) +
    geom_col(aes(x=variable, y=cor.avg, fill=ocean_region_lab), position=position_dodge(), alpha=0.7) +
    scale_fill_manual(values=col.region) +
    theme_sleek() + labs(x="", y="Average Correlation", fill="Region",
                         title=paste("Era:", era.yrs$start[e], "-", era.yrs$end[e]))

  print(g)
}
dev.off()



## Visualize NPGO fits
# Could repeat with PDO, but the above shows it is highly correlated with SST.

# Load stationary NPGO model coefficient summaries and combine them
sock.npgo.coefs <- read.csv(here(sens.fig.dir, "model_coefficients_stat_npgo_sockeye.csv"))
chum.npgo.coefs <- read.csv(here(sens.fig.dir, "model_coefficients_stat_npgo_chum.csv"))
pink.npgo.coefs <- read.csv(here(sens.fig.dir, "model_coefficients_stat_npgo_pink.csv"))

# bind them
npgo.coefs <- bind_rows(sock.npgo.coefs, chum.npgo.coefs, pink.npgo.coefs, .id="Species")
npgo.coefs <- npgo.coefs |> mutate(Species = case_when(Species == 1 ~ "Sockeye",
                                                       Species == 2 ~ "Chum",
                                                       Species == 3 ~ "Pink"))



# Load base (interaction) stationary model results for comparison

# sockeye
load(here('output', 'models', 'stat', 'sockeye', 'stat_inter.RData'), verbose=T)
summ.sock <- rstan::summary(stat_inter, pars=c("mu_gamma", "mu_kappa", "mu_chi"), probs=c(0.025, 0.5, 0.975))$summary
sock.stat.df <- data.frame(region = unique(sock$Ocean.Region2),
                           var = rep(c("SST", "Competitors", "SST x Competitors"), each=n_distinct(sock$Ocean.Region2)),
                           lower_2.5 = summ.sock[,"2.5%"],
                           mean = summ.sock[,"mean"],
                           upper_97.5 = summ.sock[,"97.5%"])

# chum
load(here('output', 'models', 'stat', 'chum', 'stat_inter.RData'), verbose=T)
summ.chum <- rstan::summary(stat_inter, pars=c("mu_gamma", "mu_kappa", "mu_chi"), probs=c(0.025, 0.5, 0.975))$summary
chum.stat.df <- data.frame(region = unique(chum$Ocean.Region2),
                           var = rep(c("SST", "Competitors", "SST x Competitors"), each=n_distinct(chum$Ocean.Region2)),
                           lower_2.5 = summ.chum[,"2.5%"],
                           mean = summ.chum[,"mean"],
                           upper_97.5 = summ.chum[,"97.5%"])

# pink
load(here('output', 'models', 'stat', 'pink', 'stat_inter.RData'), verbose=T)
summ.pink <- rstan::summary(stat_inter, pars=c("mu_gamma", "mu_kappa", "mu_chi"), probs=c(0.025, 0.5, 0.975))$summary
pink.stat.df <- data.frame(region = unique(pink$Ocean.Region2),
                           var = rep(c("SST", "Competitors", "SST x Competitors"), each=n_distinct(pink$Ocean.Region2)),
                           lower_2.5 = summ.pink[,"2.5%"],
                           mean = summ.pink[,"mean"],
                           upper_97.5 = summ.pink[,"97.5%"])

# bind them
stat.df <- bind_rows(sock.stat.df, chum.stat.df, pink.stat.df, .id="Species")
stat.df <- stat.df |> mutate(Species = case_when(Species == 1 ~ "Sockeye",
                                                 Species == 2 ~ "Chum",
                                                 Species == 3 ~ "Pink"))



# Load era model w/ NPGO fits

# sockeye
load(file=here('sensitivity-analyses', 'fits', 'sockeye', 'era_npgo.RData'), verbose=T)
df.era.npgo.sock <- era_hb_param_df(era_npgo, par=c("gamma", "kappa", "lambda"), mu=T, lower_CI=2.5, upper_CI=97.5, info=sock.info)
df.era.npgo.sock$varnam[is.na(df.era.npgo.sock$varnam)] <- "NPGO"

# chum
load(file=here('sensitivity-analyses', 'fits', 'chum', 'era_npgo.RData'), verbose=T)
df.era.npgo.chum <- era_hb_param_df(era_npgo, par=c("gamma", "kappa", "lambda"), mu=T, lower_CI=2.5, upper_CI=97.5, info=chum.info)
df.era.npgo.chum$varnam[is.na(df.era.npgo.chum$varnam)] <- "NPGO"

# pink
load(file=here('sensitivity-analyses', 'fits', 'pink', 'era_npgo.RData'), verbose=T)
df.era.npgo.pink <- era_hb_param_df(era_npgo, par=c("gamma", "kappa", "lambda"), mu=T, lower_CI=2.5, upper_CI=97.5, info=pink.info)
df.era.npgo.pink$varnam[is.na(df.era.npgo.pink$varnam)] <- "NPGO"

# bind them
df.era.npgo <- bind_rows(df.era.npgo.sock, df.era.npgo.chum, df.era.npgo.pink, .id="Species")
df.era.npgo <- df.era.npgo |> mutate(Species = case_when(Species == 1 ~ "Sockeye",
                                                 Species == 2 ~ "Chum",
                                                 Species == 3 ~ "Pink"))


# Load base (SST + Comp) era model fits

# sockeye
load(here('output', 'models', 'dyn', 'sockeye', 'hbm_era_2c.RData'), verbose=T)
df.era.sock <- era_hb_param_df(era.2c, par=c("gamma", "kappa"), mu=T, lower_CI=2.5, upper_CI=97.5, info=sock.info)

# chum
load(here('output', 'models', 'dyn', 'chum', 'hbm_era_2c.RData'), verbose=T)
df.era.chum <- era_hb_param_df(era.2c, par=c("gamma", "kappa"), mu=T, lower_CI=2.5, upper_CI=97.5, info=chum.info)

# pink
load(here('output', 'models', 'dyn', 'pink', 'hbm_era_2c.RData'), verbose=T)
df.era.pink <- era_hb_param_df(era.2c, par=c("gamma", "kappa"), mu=T, lower_CI=2.5, upper_CI=97.5, info=pink.info)

# bind species
df.era <- bind_rows(df.era.sock, df.era.chum, df.era.pink, .id="Species")
df.era <- df.era |> mutate(Species = case_when(Species == 1 ~ "Sockeye",
                                                 Species == 2 ~ "Chum",
                                                 Species == 3 ~ "Pink"))



# Combine era and stationary model coefficients -- NPGO model
df.era.npgo <- ocean_region_lab(df.era.npgo) |> rename(Lower_95_CI = lower, Upper_95_CI = upper, Mean=reg_mean) |> mutate(Era = factor(era, levels=c("Early", "Middle", "Late", "All"))) |> select(!c("ystart", "yend", "era"))
npgo.coefs <- npgo.coefs |> rename(ocean_region_lab=Ocean.Region, varnam=Coefficient) |>
  mutate(varnam = stringr::str_replace_all(varnam, "Comp", "Competitors"))
all.coefs.npgo <- bind_rows(df.era.npgo, npgo.coefs, .id="Model") # bind together
all.coefs.npgo[is.na(all.coefs.npgo$Era),"Era"] <- "All"
all.coefs.npgo$varnam <- gsub("SST x Competitors", "SST x \nComp.", all.coefs.npgo$varnam) # Have to change this back to fit on plot



# Combine era and stationary model coefficients -- base model
df.era <- df.era |> select(!c("ystart", "yend", "var")) |>
  rename(region=Ocean.Region2, mean=reg_mean, lower_2.5=lower,
                 upper_97.5=upper, var=varnam) |>
  mutate(Era = factor(era, levels=c("Early", "Middle", "Late", "All"))) |>
  select(!c("era"))
all.coefs.base <- bind_rows(df.era, stat.df, .id="Model") # bind together
all.coefs.base[is.na(all.coefs.base$Era), "Era"] <- "All"
names(all.coefs.base)[names(all.coefs.base)=="var"] <- "varnam"
all.coefs.base <- ocean_region_lab(all.coefs.base, var="region")
all.coefs.base$varnam <- gsub("SST x Competitors", "SST x \nComp.", all.coefs.base$varnam) # Have to change this label to fit on plot

# Make stat coefs labelled as "All" era

# Make master plots, one per species:
pdf(file=here(sens.fig.dir, "npgo-era-dot.pdf"), onefile = T)

for(spp in c("Sockeye", "Pink", "Chum")) {
  g <- all.coefs.npgo |> filter(Species==spp) |>
    ggplot() +
    geom_hline(yintercept = 0, color = "grey50", linetype = 2, linewidth = 0.25) +
    geom_point(data=filter(all.coefs.base, Species==spp),
               aes(y=mean, x=Era, shape=Model), col="grey50", size=3, position=position_nudge(x=.1),
               alpha=0.9) +
    geom_segment(data=filter(all.coefs.base, Species==spp),
                 aes(y=lower_2.5, yend=upper_97.5, x=Era, xend=Era),
                 col="grey50", linewidth=.8, position=position_nudge(x=.1), alpha=0.9) + # base models
    geom_point(aes(y=Mean, x=Era, col=ocean_region_lab, shape=Model), size=3, alpha=0.9) +
    geom_segment(aes(y=Lower_95_CI, yend=Upper_95_CI, x=Era, xend=Era, col=ocean_region_lab),
                 linewidth=.8, alpha=0.9) + # npgo models
    facet_grid(cols=vars(factor(varnam, levels=c("SST", "Competitors", "SST x \nComp.", "NPGO"))), rows=vars(ocean_region_lab), scales="free", space="free_x") +
    scale_shape_manual(values=c(18,15), guide="none") +
    scale_colour_manual(values=col.region, guide="none") +
    labs(x="Time Period", y="Covariate effect", title=spp) +  theme_sleek()

  print(g)
}
dev.off()




