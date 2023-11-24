## Bayesian model inference

## Set Species to most recently used -----
speciesFlag = tolower(unique(data_master$Species))


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
names(shp.reg) <- unique(info_master$ocean_label2)


## Two-covariate models ---------------------------------------------------

### --- Eras model: Data

# Stock-specific dataframe
df.era.st.2c <- era_hb_param_df(era.2c, par=c("gamma", "kappa")) # seems to work


# Summarized dataframe (regional-level)
df.era.reg.2c <- era_hb_param_df(era.2c, par=c("gamma", "kappa"), mu = TRUE) # seems to work


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
  geom_point(aes(x = mu, y = Stock, color = interaction(Ocean.Region2, era), shape = Ocean.Region2)) +
  geom_segment(data = df.era.reg.2c, aes(y = ystart, yend = yend, x = reg_mean, xend=reg_mean,
                                           color = interaction(Ocean.Region2, era)), linewidth = 0.25) +
  geom_rect(data = df.era.reg.2c, aes(xmin = lower_10, xmax = upper_90, ymin = ystart,
                                        ymax = yend, fill = interaction(Ocean.Region2, era)), alpha=0) +
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
  geom_point(aes(x = mu, y = Stock, color = Ocean.Region2, shape = Ocean.Region2)) +
  geom_segment(aes(y = Stock, yend = Stock, x = lower_10, xend = upper_90,
                   color = Ocean.Region2), linewidth = 0.25) +
  geom_segment(data = df.era.reg.2c, aes(y = ystart, yend = yend, x = reg_mean, xend=reg_mean,
                                           color = Ocean.Region2), linewidth = 0.25) +
  geom_rect(data = df.era.reg.2c, aes(xmin = lower_10, xmax = upper_90, ymin = ystart,
                                        ymax = yend, fill = Ocean.Region2), alpha=0.2) +
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
  geom_path(aes(x=x, y=dens, group=stock, col=Ocean.Region2), alpha=0.2) + 
  geom_path(data=dens.df.reg.2c, aes(x=x, y=dens, col=Ocean.Region2), alpha=0.85, linewidth=1) +
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

# Summarized dataframe (regional-level)
# gamma/kappa are series-specific; no mu output. Summarize stocks instead
df.dyn.reg.2c <- dplyr::summarize(df.dyn.st.2c, 
                                  reg_mean=mean(mu, na.rm=T), 
                                  n_stk=n_distinct(Stock),
                                  lower_10=quantile(mu, 0.1), 
                                  upper_90=quantile(mu, 0.9), 
                                  .by=c(Ocean.Region2, BY, varnam))
df.dyn.reg.2c <- ddply(df.dyn.reg.2c, .(Ocean.Region2), dplyr::filter, n_stk >= max(n_stk)*0.1) # remove years with less than 10% of stocks observed


### --- Dynamic model: Figures 

# Stacked gamma timeseries / sparkline: 2-covar
g <- ggplot(df.dyn.st.2c) + 
  geom_line(aes(x=BY, y=mu, col=Ocean.Region2)) + 
  facet_grid(rows=vars(Stock), cols=vars(varnam), switch ="y", scales="free_y", as.table=F) + 
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
pdf(here(fig.dir, "dyn_stack_2c.pdf"), width=8, height=10)
print(g)
dev.off()

# Grouped gamma timeseries : 2-covar
g <- ggplot(df.dyn.reg.2c) +
  geom_line(data= df.dyn.st.2c, aes(x=BY, y=mu, group=Stock, col=Ocean.Region2), alpha=0.2) +
  geom_line(aes(x=BY, y=reg_mean, col=Ocean.Region2), linewidth=1) +
  geom_ribbon(aes(x=BY, y=reg_mean, ymin=lower_10, ymax=upper_90, fill=Ocean.Region2), alpha=0.2) +
  facet_grid(rows=vars(Ocean.Region2), cols=vars(varnam)) + 
  ylim(c(-1,1)) + 
  scale_colour_manual(values=col.region, aesthetics=c("colour", "fill")) +
  theme_minimal()
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
  geom_line(data= df.dyn.st.2c, aes(x=BY, y=mu, group=Stock, col=Ocean.Region2), alpha=0.2) +
  geom_line(aes(x=BY, y=reg_mean, col=Ocean.Region2), linewidth=1) +
  geom_ribbon(aes(x=BY, y=reg_mean, ymin=lower_10, ymax=upper_90, fill=Ocean.Region2), alpha=0.2) +
  geom_segment(data=df.era.reg.2c, aes(x=xstart, xend=xend, y=reg_mean, yend=reg_mean), col="gray40", linetype="dashed") +
  facet_grid(rows=vars(Ocean.Region2), cols=vars(varnam)) + 
  ylim(c(-1,1)) + 
  scale_colour_manual(values=rev(col.region), aesthetics=c("colour", "fill")) +
  scale_y_continuous(limits=c(-.75,.75), breaks=c(-0.5, 0, 0.5), oob=scales::squish) +
  theme_sleek() +
  theme(legend.position = "none") + labs(x="Brood Year", y="Mean covariate effects")
pdf(here(fig.dir, "dyn_era_2c.pdf"))  
print(g)
dev.off()




## Single covariate models ----------------------------------------------------


### ----- Eras model: Data  

lst.era <- list(era.sst, era.comp)
probs <- c(0.025, 0.05, 0.10, 0.50, 0.90, 0.95, 0.975)


# Stock-specific dataframes
dfl.era.st <- lapply(seq_along(lst.era), function(i){
  summ <- rstan::summary(lst.era[[i]], pars=c(paste0("gamma", 1:3)), probs=probs)[[1]]
  df <- data.frame(Stock = rep(info_master$Stock, 3),
                   Ocean.Region2 = rep(info_master$Ocean.Region2, 3),
                   ocean_region_lab = ocean_region_lab(info_master$Ocean.Region2),
                   mu = summ[, "mean"],
                   se = summ[, "se_mean"],
                   lower_10 = summ[, "10%"],
                   upper_90 = summ[ , "90%"],
                   era = case_when(str_extract(rownames(summ), "\\d") == "1" ~ "Early",
                                   str_extract(rownames(summ), "\\d") == "2" ~ "Middle",
                                   str_extract(rownames(summ), "\\d") == "3" ~ "Late",
                                   .ptype=factor( levels=c("Early", "Middle", "Late"))) )
  return(df) }
)
names(dfl.era.st) <- c("SST", "Competitors")
dfl.era.st <- lapply(dfl.era.st, ocean_region_lab)

# Summarized dataframes (regional-level)
reg_start <- info_master$Stock[match(unique(info_master$Ocean.Region2), info_master$Ocean.Region2)]
reg_end <- c(info_master$Stock[match(unique(info_master$Ocean.Region2), info_master$Ocean.Region2)-1], 
             info_master$Stock[nrow(info_master)])

dfl.era.reg <- lapply(seq_along(lst.era), function(i){
  summ <- rstan::summary(lst.era[[i]], pars=c(paste0("mu_gamma", 1:3)), probs=probs)[[1]]
  df <- data.frame(Ocean.Region2 = rep(unique(info_master$Ocean.Region2), 3),
                   ystart = rep(reg_start, 3),
                   yend = rep(reg_end, 3),
                   reg_mean = summ[, "mean"],
                   reg_se = summ[, "se_mean"],
                   lower_10 = summ[, "10%"],
                   upper_90 = summ[ , "90%"],
                   era = case_when(str_extract(rownames(summ), "\\d") == "1" ~ "Early",
                                   str_extract(rownames(summ), "\\d") == "2" ~ "Middle",
                                   str_extract(rownames(summ), "\\d") == "3" ~ "Late",
                                   .ptype=factor( levels=c("Early", "Middle", "Late"))) )
  return(df) }
)
names(dfl.era.reg) <- c("SST", "Competitors")
dfl.era.reg <- lapply(dfl.era.reg, ocean_region_lab)


# Density dataframe - stock-level
dens.dfl.st <- list()
for(i in seq_along(lst.era)) {
  post <- rstan::extract(lst.era[[i]], pars=c(paste0("gamma", c(1:3))))
  dens.l <- lapply(post, function(x){
    dens.out <- col_density(x, plot.it=F)
    dens.list <- lapply(dens.out, adply, .margins=c(1,2))
    dens.df.sst <- join(dens.list$x, dens.list$y, by=c("X1", "X2"))
    names(dens.df.sst) <- c("n", "stock", "gamma.x", "gamma.dens")
    dens.df.sst <- mutate(dens.df.sst, stock=levels(data_master$Stock)[as.numeric(stock)] ) 
    return(dens.df.sst)
  } )
  summ.dens <- bind_rows(dens.l, .id="par")
  dens.df.st <- data.frame(summ.dens,
                           Ocean.Region2 = info_master$Ocean.Region2[
                             match(summ.dens$stock, info_master$Stock)],
                           era = case_when(
                             str_extract(summ.dens$par, "\\d") == "1" ~ "Early",
                             str_extract(summ.dens$par, "\\d") == "2" ~ "Middle",
                             str_extract(summ.dens$par, "\\d") == "3" ~ "Late",
                             .ptype=factor( levels=c("Early", "Middle", "Late")))
                            )
  dens.dfl.st[[i]] <- dens.df.st 
}
names(dens.dfl.st) <- c("SST", "Competitors")
dens.dfl.st <- lapply(dens.dfl.st, ocean_region_lab)


# mean (regional) density dataframe 
dens.dfl.reg <- list()
for(i in seq_along(lst.era)) {
  post <- rstan::extract(lst.era[[i]], pars=c(paste0("mu_gamma", c(1:3))))
  dens.l <- lapply(post, function(x){
    dens.out <- col_density(x, plot.it=F)
    dens.list <- lapply(dens.out, adply, .margins=c(1,2))
    dens.df <- join(dens.list$x, dens.list$y, by=c("X1", "X2"))
    names(dens.df) <- c("n", "region", "x", "dens")
    dens.df$Ocean.Region2 <- case_when(dens.df$region == 1 ~ "WC",
                                dens.df$region == 2 ~ "SEAK",
                                dens.df$region == 3 ~ "GOA",
                                dens.df$region == 4 ~ "BS")
    return(dens.df)
  } )
  summ.dens <- bind_rows(dens.l, .id="par")
  dens.df.reg <- data.frame(summ.dens,
                           era = case_when(
                             str_extract(summ.dens$par, "\\d") == "1" ~ "Early",
                             str_extract(summ.dens$par, "\\d") == "2" ~ "Middle",
                             str_extract(summ.dens$par, "\\d") == "3" ~ "Late",
                             .ptype=factor( levels=c("Early", "Middle", "Late")))
  )
  dens.dfl.reg[[i]] <- dens.df.reg
}
names(dens.dfl.reg) <- c("SST", "Competitors")
dens.dfl.reg <- lapply(dens.dfl.reg, ocean_region_lab)



### --- Eras model : Figures

## SST 

# Caterpillar plot (combined) - sst
g <- ggplot(dfl.era.st$SST) +
  geom_vline(xintercept = 0, color = "grey50", linetype = 2, linewidth = 0.25) +
  geom_point(aes(x = mu, y = Stock, color = interaction(ocean_region_lab, era), shape = ocean_region_lab)) +
  geom_segment(data = dfl.era.reg$SST, aes(y = ystart, yend = yend, x = reg_mean, xend=reg_mean,
                                 color = interaction(ocean_region_lab, era)), linewidth = 0.25) +
  geom_rect(data = dfl.era.reg$SST, aes(xmin = lower_10, xmax = upper_90, ymin = ystart,
                         ymax = yend, fill = interaction(ocean_region_lab, era)), alpha=0) +
  scale_color_manual(values = col.eras) +
  scale_shape_manual(values = c(15:18), guide = "legend") +
  scale_fill_manual(values = col.eras, guide="legend") +
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
pdf(here(fig.dir, "eras_coef_dot_sst_comb.pdf"), width=5, height=8)
print(g)
dev.off()

# Caterpillar plot (panel) - sst
g <- ggplot(dfl.era.st$SST) +
  geom_vline(xintercept = 0, color = "grey50", linetype = 2, linewidth = 0.25) +
  geom_point(aes(x = mu, y = Stock, color = ocean_region_lab, shape = ocean_region_lab)) +
  geom_segment(aes(y = Stock, yend = Stock, x = lower_10, xend = upper_90,
                   color = ocean_region_lab), linewidth = 0.25) +
  geom_segment(data = dfl.era.reg$SST, aes(y = ystart, yend = yend, x = reg_mean, xend=reg_mean,
                                           color = ocean_region_lab), linewidth = 0.25) +
  geom_rect(data = dfl.era.reg$SST, aes(xmin = lower_10, xmax = upper_90, ymin = ystart,
                                        ymax = yend, fill = ocean_region_lab), alpha=0.2) +
  facet_wrap(vars(era), ncol=3) +
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
pdf(here(fig.dir, "eras_coef_dot_sst_panel.pdf"), width=5, height=8)
print(g)
dev.off()

# Posterior density plot - sst
g <- ggplot(dens.dfl.st$SST) + 
  geom_vline(xintercept = 0, color = "grey50", linetype = 2, linewidth = 0.25) +
  geom_path(aes(x=gamma.x, y=gamma.dens, group=stock, col=ocean_region_lab), alpha=0.2) + 
  geom_path(data=dens.dfl.reg$SST, aes(x=x, y=dens, col=ocean_region_lab), alpha=0.85, linewidth=1) +
  facet_wrap(~era, nrow=3) +
  scale_colour_manual(values=col.region) +
  coord_cartesian(xlim=c(-1, 1)) +
  theme_minimal() + labs(x="SST effect", y="", col="Ocean Region") +
  theme(axis.text.y=element_blank())
pdf(here(fig.dir, "eras_sst_dens.pdf"))
print(g)
dev.off()


## Competitors

# Caterpillar plot (combined) - comp
g <- ggplot(dfl.era.st$Competitors) +
  geom_vline(xintercept = 0, color = "grey50", linetype = 2, linewidth = 0.25) +
  geom_point(aes(x = mu, y = Stock, color = interaction(ocean_region_lab, era), shape = ocean_region_lab)) +
  #geom_segment(aes(y = Stock, yend = Stock, x = `2.5%`, xend = `97.5%`,
  #                 color = ocean_region_lab), linewidth = 0.25) +
  geom_segment(data = dfl.era.reg$Competitors, aes(y = ystart, yend = yend, x = reg_mean, xend=reg_mean,
                                           color = interaction(ocean_region_lab, era)), linewidth = 0.25) +
  geom_rect(data = dfl.era.reg$Competitors, aes(xmin = lower_10, xmax = upper_90, ymin = ystart,
                                        ymax = yend, fill = interaction(ocean_region_lab, era)), alpha=0) +
  scale_color_manual(values = col.eras) +
  scale_shape_manual(values = c(15:18), guide = "legend") +
  scale_fill_manual(values = col.eras, guide="legend") +
  labs(x = "Competitors Coefficient",
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
pdf(here(fig.dir, "eras_coef_dot_comp_comb.pdf"), width=5, height=8)
print(g)
dev.off()


# Caterpillar plot (panel) - Competitors
g <- ggplot(dfl.era.st$Competitors) +
  geom_vline(xintercept = 0, color = "grey50", linetype = 2, linewidth = 0.25) +
  geom_point(aes(x = mu, y = Stock, color = ocean_region_lab, shape = ocean_region_lab)) +
  geom_segment(aes(y = Stock, yend = Stock, x = lower_10, xend = upper_90,
                   color = ocean_region_lab), linewidth = 0.25) +
  geom_segment(data = dfl.era.reg$Competitors, aes(y = ystart, yend = yend, x = reg_mean, xend=reg_mean,
                                           color = ocean_region_lab), linewidth = 0.25) +
  geom_rect(data = dfl.era.reg$Competitors, aes(xmin = lower_10, xmax = upper_90, ymin = ystart,
                                        ymax = yend, fill = ocean_region_lab), alpha=0.2) +
  facet_wrap(vars(era), ncol=3) +
  scale_color_manual(values = col.region) +
  scale_shape_manual(values = c(15:18), guide = "legend") +
  scale_fill_manual(values = col.region, guide="legend") +
  labs(x = "Competitor Coefficient",
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
pdf(here(fig.dir, "eras_coef_dot_comp_panel.pdf"), width=5, height=8)
print(g)
dev.off()


# Posterior density plot - comp
g <- ggplot(dens.dfl.st$Competitors) + 
  geom_vline(xintercept = 0, color = "grey50", linetype = 2, linewidth = 0.25) +
  geom_path(aes(x=gamma.x, y=gamma.dens, group=stock, col=ocean_region_lab), alpha=0.2) + 
  geom_path(data=dens.dfl.reg$Competitors, aes(x=x, y=dens, col=ocean_region_lab), alpha=0.85, linewidth=1) +
  facet_wrap(~era, nrow=3) +
  scale_colour_manual(values=col.region) +
  coord_cartesian(xlim=c(-1, 1)) +
  theme_minimal() + labs(x="Competitor effect", y="", col="Ocean Region") +
  theme(axis.text.y=element_blank())
pdf(here(fig.dir, "eras_comp_dens.pdf"))
print(g)
dev.off()



### --- Dynamic model - Data 

lst.dyn <- list(SST  = dyn.sst,
                Comp = dyn.comp)


# Stock-specific dataframe
dfl.dyn.st <- lapply(seq_along(lst.dyn), function(i) {
  probs <- c(0.025, 0.05, 0.10, 0.50, 0.90, 0.95, 0.975)
  summ <- rstan::summary(lst.dyn[[i]], pars = "gamma", probs = probs)[[1]]
  df <- data.frame(Stock = data_master$Stock,
                   Ocean.Region2 = data_master$Ocean.Region2,
                   BY = data_master$BY,
                   mu_gamma = summ[ , "mean"],
                   se_gamma = summ[ , "se_mean"],
                   lower_10 = summ[ , "10%"],
                   upper_90 = summ[ , "90%"],
                   # lower = summ[ , "5%"],
                   # upper = summ[ , "95%"],
                   var = names(lst.dyn)[i])
  return(df)
})

names(dfl.dyn.st) <- c("SST", "Comp")
dfl.dyn.st <- lapply(dfl.dyn.st, ocean_region_lab)

# Regional summary dataframe
dfl.dyn.reg <- lapply(seq_along(dfl.dyn.st), function(i){
  df <- dfl.dyn.st[[i]]
  reg_df <- dplyr::summarize(df, 
                   reg_mean_gamma=mean(mu_gamma, na.rm=T), 
                   gamma_10=quantile(mu_gamma, 0.1), 
                   gamma_90=quantile(mu_gamma, 0.9), 
                   n_stk = n_distinct(Stock),
                   .by=c(Ocean.Region2, BY))
  reg_df$even.odd <- ifelse(gtools::odd(reg_df$BY), "odd", "even")
  return(reg_df)
                })
names(dfl.dyn.reg) <- c("SST", "Comp")
dfl.dyn.reg <- lapply(dfl.dyn.reg, ocean_region_lab)


### --- Dynamic model - Figures

## SST

# Stacked gamma timeseries / sparkline: SST
g <- ggplot(dfl.dyn.st[["SST"]]) + 
  geom_line(aes(x=BY, y=mu_gamma, col=Ocean.Region2)) + 
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
pdf(here(fig.dir, "dyn_stack_sst.pdf"), width=5, height=10)
print(g)
dev.off()

# Grouped gamma timeseries : SST
g <- ggplot(dfl.dyn.reg[["SST"]]) +
  geom_line(data= dfl.dyn.st[["SST"]], aes(x=BY, y=mu_gamma, group=Stock, col=ocean_region_lab), alpha=0.2) +
  geom_line(aes(x=BY, y=reg_mean_gamma, col=ocean_region_lab, group=even.odd), linewidth=1) +
  #geom_ribbon(aes(x=BY, y=reg_mean_gamma, ymin=gamma_10, ymax=gamma_90, fill=ocean_region_lab), alpha=0.2) +
  facet_wrap(vars(ocean_region_lab), nrow=2) + ylim(c(-1,1)) + 
  scale_colour_manual(values=col.region, aesthetics=c("colour", "fill")) +
  theme_minimal() + labs(x= "Brood Year", y="SST effect") +
  theme(legend.position="none")
pdf(here(fig.dir, "dyn_sst_grouped-evenodd.pdf"))
print(g)
dev.off()
  

## Competitors

# Stacked gamma timeseries / sparkline: Comp
g <- ggplot(dfl.dyn.st[["Comp"]]) + 
  geom_line(aes(x=BY, y=mu_gamma, col=Ocean.Region2)) + 
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
pdf(here(fig.dir, "dyn_stack_comp.pdf"), width=5, height=10)
print(g)
dev.off()

# Grouped gamma timeseries : Comp
g <- ggplot(dfl.dyn.reg[["Comp"]]) +
  geom_line(data= dfl.dyn.st[["Comp"]], aes(x=BY, y=mu_gamma, group=Stock, col=Ocean.Region2), alpha=0.2) +
  #geom_line(aes(x=BY, y=reg_mean_gamma, col=Ocean.Region2), linewidth=1) +
  #geom_ribbon(aes(x=BY, y=reg_mean_gamma, ymin=gamma_10, ymax=gamma_90, fill=Ocean.Region2), alpha=0.2) +
  facet_wrap(vars(Ocean.Region2), nrow=2) + ylim(c(-1,1)) + 
  scale_colour_manual(values=col.region, aesthetics=c("colour", "fill")) +
  theme_minimal() + labs(x= "Brood Year", y="Competitor effect") +
  theme(legend.position="none")
pdf(here(fig.dir, "dyn_comp_grouped.pdf")) 
print(g)
dev.off()



## Additional plots
## Dot-plots of alpha, beta, sigma 
plot_hbm_dot(era.sst, pdf.file = here(fig.dir, "era_sst_dot.pdf"))
plot_hbm_dot(dyn.sst, pdf.file = here(fig.dir, "dyn_sst_dot.pdf"))

## Posterior density of alpha
plot_hbm_dens(era.sst, pdf.file = here(fig.dir, "era_sst_dens.pdf"))
plot_hbm_dens(dyn.sst, pdf.file = here(fig.dir, "dyn_sst_dens.pdf"))


