## Bayesian model inference

if(!dir.exists("./figures/dyn/hbm_inf/"))
    dir.create("./figures/dyn/hbm_inf/")

# Load fits
for(i in list.files(path = "./output/models/dyn/", pattern = "*.RData$")) {
  load(paste("./output/models/dyn/", i, sep = ""))
}


# Set colours 
col.region <- chroma::qpal(7, luminance = 40)[c(1, 3, 5, 7)]
col.eras <- c("#4db8ff", "#b3a100", "#ff80d7", "#00b39e", 
              "#0070BDFF", "#6D6200FF", "#BC007FFF", "#008070",
              "#002e4d", "#4d4500", "#4d0034",  "#00332d")



## Two-covariate models ---------------------------------------------------

### --- Eras model: Data

# Stock-specific dataframe
probs <- c(0.025, 0.05, 0.10, 0.50, 0.90, 0.95, 0.975)
summ <- rstan::summary(era.2c, pars = c(paste0("gamma", 1:3), paste0("kappa", 1:3)), probs = probs)[[1]]
df.era.st.2c <- data.frame(Stock = rep(sock.info$Stock, 3), 
                           Ocean.Region2 = rep(sock.info$Ocean.Region2, 3),
                           mu = summ[, "mean"],
                           se = summ[, "se_mean"],
                           lower_10 = summ[, "10%"],
                           upper_90 = summ[ , "90%"],
                           var = stringr::str_extract(rownames(summ), "[a-z]+"),
                           era = rep(c("early", "middle", "late"), each=56) #hard coding because I don't understand regular expressions
)

# Summarized dataframe (regional-level)
df.era.reg.2c <- dplyr::summarize(df.era.st.2c, 
                                  reg_mean=mean(mu), 
                                  lower_10=quantile(mu, 0.1), 
                                  upper_90=quantile(mu, 0.9), 
                                  .by=c(Ocean.Region2, var, era))
df.era.reg.2c$ystart <- rep(c("Early Stuart", "Alastair", "Copper", "Nelson"), 3)
df.era.reg.2c$yend <- rep(c("Atnarko", "Chilkoot", "Chignik Lake", "Goodnews"), 3)


# Density dataframe 
post <- rstan::extract(era.2c, pars=c(paste0("gamma", c(1:3)),
                                       paste0("kappa", c(1:3))))
dens.l <- lapply(post, function(x){
  dens.out <- col_density(x, plot.it=F)
  dens.list <- lapply(dens.out, adply, .margins=c(1,2))
  dens.df <- join(dens.list$x, dens.list$y, by=c("X1", "X2"))
  names(dens.df) <- c("n", "stock", "x", "dens")
  dens.df <- mutate(dens.df, stock=levels(sock$Stock)[as.numeric(stock)] ) 
  return(dens.df)
} )
dens.df.st.2c <- bind_rows(dens.l, .id="par") 
dens.df.st.2c$era <- case_when(dens.df.st.2c$par %in% c("gamma1", "kappa1") ~ "early",
                             dens.df.st.2c$par %in% c("gamma2", "kappa2") ~ "middle",
                             dens.df.st.2c$par %in% c("gamma3", "kappa3") ~ "late") 
dens.df.st.2c$era <- factor(dens.df.st.2c$era, levels=c("early", "middle", "late"))
dens.df.st.2c$var <- stringr::str_extract(dens.df.st.2c$par, "[a-z]+")
dens.df.st.2c <- ocean_region_col(dens.df.st.2c, stock.col = "stock")

# Mean (regional) density dataframe 
post <- rstan::extract(era.2c, pars=c(paste0("mu_gamma", c(1:3)),
                                                         paste0("mu_kappa", c(1:3))))
dens.l <- lapply(post, function(x){
  dens.out <- col_density(x, plot.it=F)
  dens.list <- lapply(dens.out, adply, .margins=c(1,2))
  dens.df <- join(dens.list$x, dens.list$y, by=c("X1", "X2"))
  names(dens.df) <- c("n", "region", "x", "dens")
  dens.df$region <- case_when(dens.df$region == 1 ~ "WC",
                              dens.df$region == 2 ~ "SEAK",
                              dens.df$region == 3 ~ "GOA",
                              dens.df$region == 4 ~ "BS")
  return(dens.df)
} )
dens.df.reg.2c <- bind_rows(dens.l, .id="par") 
dens.df.reg.2c$era <- case_when(dens.df.reg.2c$par %in% c("mu_gamma1", "mu_kappa1") ~ "early",
                            dens.df.reg.2c$par %in% c("mu_gamma2", "mu_kappa2") ~ "middle",
                            dens.df.reg.2c$par %in% c("mu_gamma3", "mu_kappa3") ~ "late") 
dens.df.reg.2c$era <- factor(dens.df.reg.2c$era, levels=c("early", "middle", "late"))
dens.df.reg.2c$var <- case_when(dens.df.reg.2c$par %in% c("mu_gamma1", "mu_gamma2","mu_gamma3") ~ "gamma",
                                dens.df.reg.2c$par %in% c("mu_kappa1", "mu_kappa2",  "mu_kappa3") ~ "kappa")


### --- Eras model: Figures

# Caterpillar plot - 2-covar
g <- ggplot(df.era.st.2c) +
  geom_vline(xintercept = 0, color = "grey50", linetype = 2, linewidth = 0.25) +
  geom_point(aes(x = mu, y = Stock, color = interaction(Ocean.Region2, era), shape = Ocean.Region2)) +
  geom_segment(data = df.era.reg.2c, aes(y = ystart, yend = yend, x = reg_mean, xend=reg_mean,
                                           color = interaction(Ocean.Region2, era)), linewidth = 0.25) +
  geom_rect(data = df.era.reg.2c, aes(xmin = lower_10, xmax = upper_90, ymin = ystart,
                                        ymax = yend, fill = interaction(Ocean.Region2, era)), alpha=0) +
  facet_wrap(vars(var), ncol=2) +
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
pdf("./figures/dyn/hbm_inf/eras_coef_dot_c2_comb.pdf", width=5, height=8)
print(g)
dev.off()


# Posterior density plot - 2-covar
g <- ggplot(dens.df.st.2c) + 
  geom_vline(xintercept = 0, color = "grey50", linetype = 2, linewidth = 0.25) +
  geom_path(aes(x=x, y=dens, group=stock, col=region), alpha=0.2) + 
  geom_path(data=dens.df.reg.2c, aes(x=x, y=dens, col=region), alpha=0.85, linewidth=1) +
  scale_colour_manual(values=col.region) +
  facet_grid(rows=vars(era), cols=vars(var)) + 
  coord_cartesian(xlim=c(-1, 1)) +
  theme_minimal() + labs(x="covariate effect", y="", col="Ocean Region") +
  theme(axis.text.y=element_blank())
pdf("./figures/dyn/hbm_inf/eras_2c_dens.pdf")
print(g)
dev.off()




### --- Dynamic model: Data

# Stock-specific dataframe
probs <- c(0.025, 0.05, 0.10, 0.50, 0.90, 0.95, 0.975)
summ <- rstan::summary(dyn.2c, pars = c("gamma", "kappa"), probs = probs)[[1]]
df.dyn.st.2c <- data.frame(Stock = sock$Stock,
                   Ocean.Region2 = sock$Ocean.Region2,
                   BY = sock$BY,
                   mu = summ[, "mean"],
                   se = summ[, "se_mean"],
                   lower_10 = summ[, "10%"],
                   upper_90 = summ[ , "90%"],
                   var = stringr::str_extract(rownames(summ), "[a-z]+")
                   )

# Summarized dataframe (regional-level)
df.dyn.reg.2c <- dplyr::summarize(df.dyn.st.2c, 
                           reg_mean=mean(mu), 
                           lower_10=quantile(mu, 0.1), 
                           upper_90=quantile(mu, 0.9), 
                           .by=c(Ocean.Region2, BY, var))

### --- Dynamic model: Figures 

# Stacked gamma timeseries / sparkline: 2-covar
g <- ggplot(df.dyn.st.2c) + 
  geom_line(aes(x=BY, y=mu, col=Ocean.Region2)) + 
  facet_grid(rows=vars(Stock), cols=vars(var), switch ="y", scales="free_y", as.table=F) + 
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
pdf("./figures/dyn/hbm_inf/dyn_stack_2c.pdf", width=8, height=10)
print(g)
dev.off()

# Grouped gamma timeseries : 2-covar
g <- ggplot(df.dyn.reg.2c) +
  geom_line(data= df.dyn.st.2c, aes(x=BY, y=mu, group=Stock, col=Ocean.Region2), alpha=0.2) +
  geom_line(aes(x=BY, y=reg_mean, col=Ocean.Region2), linewidth=1) +
  geom_ribbon(aes(x=BY, y=reg_mean, ymin=lower_10, ymax=upper_90, fill=Ocean.Region2), alpha=0.2) +
  facet_grid(rows=vars(Ocean.Region2), cols=vars(var)) + 
  ylim(c(-1,1)) + 
  scale_colour_manual(values=col.region, aesthetics=c("colour", "fill")) +
  theme_minimal()
pdf("./figures/dyn/hbm_inf/dyn_2c_grouped.pdf")  
print(g)
dev.off()




## Single covariate models ----------------------------------------------------


### ----- Eras model: Data  
#TODO : make data wrangling below more transparent for me & consistent with dyn models

mcmc.era.sst  <- As.mcmc.list(era.sst, pars = pars_era_1c)
mcmc.era.comp  <- As.mcmc.list(era.comp, pars = pars_era_1c)
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


# Summarized dataframe (regional-level)
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


# Density dataframe - sst
post <- rstan::extract(era.sst, pars=c(paste0("gamma", c(1:3))))
dens.l <- lapply(post, function(x){
  dens.out <- col_density(x, plot.it=F)
  dens.list <- lapply(dens.out, adply, .margins=c(1,2))
  dens.df.sst <- join(dens.list$x, dens.list$y, by=c("X1", "X2"))
  names(dens.df.sst) <- c("n", "stock", "gamma.x", "gamma.dens")
  dens.df.sst <- mutate(dens.df.sst, stock=levels(sock$Stock)[as.numeric(stock)] ) 
  return(dens.df.sst)
} )
dens.df.sst <- bind_rows(dens.l, .id="par") # rename - this is confusing with the above naming
dens.df.sst$era <- case_when(dens.df.sst$par == "gamma1" ~ "early",
                              dens.df.sst$par == "gamma2" ~ "middle",
                              dens.df.sst$par == "gamma3" ~ "late")
dens.df.sst$era <- factor(dens.df.sst$era, levels=c("early", "middle", "late"))
dens.df.sst <- ocean_region_col(dens.df.sst, stock.col = "stock")

# mean (regional) density dataframe 
post <- rstan::extract(era.sst, pars=c(paste0("mu_gamma", c(1:3))))
dens.l <- lapply(post, function(x){
  dens.out <- col_density(x, plot.it=F)
  dens.list <- lapply(dens.out, adply, .margins=c(1,2))
  dens.df <- join(dens.list$x, dens.list$y, by=c("X1", "X2"))
  names(dens.df) <- c("n", "region", "x", "dens")
  dens.df$region <- case_when(dens.df$region == 1 ~ "WC",
                              dens.df$region == 2 ~ "SEAK",
                              dens.df$region == 3 ~ "GOA",
                              dens.df$region == 4 ~ "BS")
  return(dens.df)
} )
dens.df.reg.sst <- bind_rows(dens.l, .id="par") 
dens.df.reg.sst$era <- case_when(dens.df.reg.sst$par == "mu_gamma1" ~ "early",
                                  dens.df.reg.sst$par == "mu_gamma2" ~ "middle",
                                  dens.df.reg.sst$par == "mu_gamma3" ~ "late") 
dens.df.reg.sst$era <- factor(dens.df.reg.sst$era, levels=c("early", "middle", "late"))


# Density dataframe - comp
post <- rstan::extract(era.comp, pars=c(paste0("gamma", c(1:3))))
dens.l <- lapply(post, function(x){
  dens.out <- col_density(x, plot.it=F)
  dens.list <- lapply(dens.out, adply, .margins=c(1,2))
  dens.df.comp <- join(dens.list$x, dens.list$y, by=c("X1", "X2"))
  names(dens.df.comp) <- c("n", "stock", "gamma.x", "gamma.dens")
  dens.df.comp <- mutate(dens.df.comp, stock=levels(sock$Stock)[as.numeric(stock)] ) 
  return(dens.df.comp)
} )
dens.df.comp <- bind_rows(dens.l, .id="par") # rename - this is confusing with the above naming
dens.df.comp$era <- case_when(dens.df.comp$par == "gamma1" ~ "early",
                         dens.df.comp$par == "gamma2" ~ "middle",
                         dens.df.comp$par == "gamma3" ~ "late")
dens.df.comp$era <- factor(dens.df.comp$era, levels=c("early", "middle", "late"))
dens.df.comp <- ocean_region_col(dens.df.comp, stock.col = "stock")

# mean (regional) density dataframe - comp
post <- rstan::extract(era.comp, pars=c(paste0("mu_gamma", c(1:3))))
dens.l <- lapply(post, function(x){
  dens.out <- col_density(x, plot.it=F)
  dens.list <- lapply(dens.out, adply, .margins=c(1,2))
  dens.df <- join(dens.list$x, dens.list$y, by=c("X1", "X2"))
  names(dens.df) <- c("n", "region", "x", "dens")
  dens.df$region <- case_when(dens.df$region == 1 ~ "WC",
                              dens.df$region == 2 ~ "SEAK",
                              dens.df$region == 3 ~ "GOA",
                              dens.df$region == 4 ~ "BS")
  return(dens.df)
} )
dens.df.reg.comp <- bind_rows(dens.l, .id="par") 
dens.df.reg.comp$era <- case_when(dens.df.reg.comp$par == "mu_gamma1" ~ "early",
                                dens.df.reg.comp$par == "mu_gamma2" ~ "middle",
                                dens.df.reg.comp$par == "mu_gamma3" ~ "late") 
dens.df.reg.comp$era <- factor(dens.df.reg.comp$era, levels=c("early", "middle", "late"))


### --- Eras model : Figures

## SST 

# Caterpillar plot - sst
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
pdf("./figures/dyn/hbm_inf/eras_coef_dot_sst_comb.pdf", width=5, height=8)
print(g)
dev.off()


# Posterior density plot - sst
g <- ggplot(dens.df.sst) + 
  geom_vline(xintercept = 0, color = "grey50", linetype = 2, linewidth = 0.25) +
  geom_path(aes(x=gamma.x, y=gamma.dens, group=stock, col=region), alpha=0.2) + 
  geom_path(data=dens.df.reg.sst, aes(x=x, y=dens, col=region), alpha=0.85, linewidth=1) +
  facet_wrap(~era, nrow=3) +
  scale_colour_manual(values=col.region) +
  coord_cartesian(xlim=c(-1, 1)) +
  theme_minimal() + labs(x="covariate effect", y="", col="Ocean Region") +
  theme(axis.text.y=element_blank())
pdf("./figures/dyn/hbm_inf/eras_sst_dens.pdf")
print(g)
dev.off()


## Competitors

# Caterpillar plot - comp
g <- ggplot(dfl.era.st$Comp) +
  geom_vline(xintercept = 0, color = "grey50", linetype = 2, linewidth = 0.25) +
  geom_point(aes(x = gamma_mu, y = Stock, color = interaction(ocean_region_lab, era), shape = ocean_region_lab)) +
  #geom_segment(aes(y = Stock, yend = Stock, x = `2.5%`, xend = `97.5%`,
  #                 color = ocean_region_lab), linewidth = 0.25) +
  geom_segment(data = dfl.era.reg$Comp, aes(y = ystart, yend = yend, x = gamma_mu, xend=gamma_mu,
                                           color = interaction(ocean_region, era)), linewidth = 0.25) +
  geom_rect(data = dfl.era.reg$Comp, aes(xmin = gamma_10, xmax = gamma_90, ymin = ystart,
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
pdf("./figures/dyn/hbm_inf/eras_coef_dot_comp_comb.pdf", width=5, height=8)
print(g)
dev.off()


# Posterior density plot - comp
g <- ggplot(dens.df.comp) + 
  geom_vline(xintercept = 0, color = "grey50", linetype = 2, linewidth = 0.25) +
  geom_path(aes(x=gamma.x, y=gamma.dens, group=stock, col=region), alpha=0.2) + 
  geom_path(data=dens.df.reg.comp, aes(x=x, y=dens, col=region), alpha=0.85, linewidth=1) +
  facet_wrap(~era, nrow=3) +
  scale_colour_manual(values=col.region) +
  coord_cartesian(xlim=c(-1, 1)) +
  theme_minimal() + labs(x="covariate effect", y="", col="Ocean Region") +
  theme(axis.text.y=element_blank())
pdf("./figures/dyn/hbm_inf/eras_comp_dens.pdf")
print(g)
dev.off()



### --- Dynamic model - Data 

lst.dyn <- list(SST  = dyn.sst,
                Comp = dyn.comp)


# Stock-specific dataframe
dfl.dyn.st <- lapply(seq_along(lst.dyn), function(i) {
  probs <- c(0.025, 0.05, 0.10, 0.50, 0.90, 0.95, 0.975)
  summ <- rstan::summary(lst.dyn[[i]], pars = "gamma", probs = probs)[[1]]
  df <- data.frame(Stock = sock$Stock,
                   Ocean.Region2 = sock$Ocean.Region2,
                   BY = sock$BY,
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


# Regional summary dataframe
dfl.dyn.reg <- lapply(seq_along(dfl.dyn.st), function(i){
  df <- dfl.dyn.st[[i]]
  reg_df <- dplyr::summarize(df, 
                   reg_mean_gamma=mean(mu_gamma), 
                   gamma_10=quantile(mu_gamma, 0.1), 
                   gamma_90=quantile(mu_gamma, 0.9), 
                   .by=c(Ocean.Region2, BY))
  return(reg_df)
                })
names(dfl.dyn.reg) <- c("SST", "Comp")


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
pdf("./figures/dyn/hbm_inf/dyn_stack_sst.pdf", width=5, height=10)
print(g)
dev.off()

# Grouped gamma timeseries : SST
g <- ggplot(dfl.dyn.reg[["SST"]]) +
  geom_line(data= dfl.dyn.st[["SST"]], aes(x=BY, y=mu_gamma, group=Stock, col=Ocean.Region2), alpha=0.2) +
  geom_line(aes(x=BY, y=reg_mean_gamma, col=Ocean.Region2), linewidth=1) +
  geom_ribbon(aes(x=BY, y=reg_mean_gamma, ymin=gamma_10, ymax=gamma_90, fill=Ocean.Region2), alpha=0.2) +
  facet_wrap(vars(Ocean.Region2), nrow=2) + ylim(c(-1,1)) + 
  scale_colour_manual(values=col.region, aesthetics=c("colour", "fill")) +
  theme_minimal()
pdf("./figures/dyn/hbm_inf/dyn_sst_grouped.pdf")  
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
pdf("./figures/dyn/hbm_inf/dyn_stack_comp.pdf", width=5, height=10)
print(g)
dev.off()

# Grouped gamma timeseries : Comp
g <- ggplot(dfl.dyn.reg[["Comp"]]) +
  geom_line(data= dfl.dyn.st[["Comp"]], aes(x=BY, y=mu_gamma, group=Stock, col=Ocean.Region2), alpha=0.2) +
  geom_line(aes(x=BY, y=reg_mean_gamma, col=Ocean.Region2), linewidth=1) +
  geom_ribbon(aes(x=BY, y=reg_mean_gamma, ymin=gamma_10, ymax=gamma_90, fill=Ocean.Region2), alpha=0.2) +
  facet_wrap(vars(Ocean.Region2), nrow=2) + ylim(c(-1,1)) + 
  scale_colour_manual(values=col.region, aesthetics=c("colour", "fill")) +
  theme_minimal()
pdf("./figures/dyn/hbm_inf/dyn_comp_grouped.pdf")  
print(g)
dev.off()






## Additional plots
## Dot-plots of alpha, beta, sigma 
plot_hbm_dot(era.sst, pdf.file = "./figures/dyn/hbm_fit/era_sst_dot.pdf")
plot_hbm_dot(dyn.sst, pdf.file = "./figures/dyn/hbm_fit/dyn_sst_dot.pdf")

## Posterior density of alpha
plot_hbm_dens(era.sst, pdf.file = "./figures/dyn/hbm_fit/era_sst_dens.pdf")
plot_hbm_dens(dyn.sst, pdf.file = "./figures/dyn/hbm_fit/dyn_sst_dens.pdf")


