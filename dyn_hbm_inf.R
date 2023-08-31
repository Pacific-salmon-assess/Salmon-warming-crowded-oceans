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
                           var = stringr::str_extract(rownames(summ), "\\D+"),
                           varnam = case_when(grepl("^gamma", rownames(summ)) ~ "SST",
                                              grepl("^kappa", rownames(summ)) ~ "Competitors"),
                           era = case_when(str_extract(rownames(summ), "\\d") == "1" ~ "Early",
                                           str_extract(rownames(summ), "\\d") == "2" ~ "Middle",
                                           str_extract(rownames(summ), "\\d") == "3" ~ "Late",
                                           .ptype=factor( levels=c("Early", "Middle", "Late")))
                             )

# Summarized dataframe (regional-level)
reg_start <- sock.info$Stock[match(unique(sock.info$Ocean.Region2), sock.info$Ocean.Region2)]
reg_end <- c(sock.info$Stock[match(unique(sock.info$Ocean.Region2), sock.info$Ocean.Region2)-1], 
             sock.info$Stock[nrow(sock.info)])
summ <- rstan::summary(era.2c, pars = c(paste0("mu_gamma", 1:3), paste0("mu_kappa", 1:3)), probs = probs)[[1]]
df.era.reg.2c <- data.frame(Ocean.Region2 = rep(unique(sock.info$Ocean.Region2), 3),
                            ystart = rep(reg_start, 3),
                            yend = rep(reg_end, 3),
                            reg_mean = summ[, "mean"], 
                            reg_se = summ[ ,"se_mean"],
                            lower_10 = summ[ , "10%"],
                            upper_90 = summ[ , "90%"], 
                            var = stringr::str_extract(rownames(summ), "\\D+"),
                            varnam = case_when(grepl("^mu_gamma", rownames(summ)) ~ "SST",
                                               grepl("^mu_kappa", rownames(summ)) ~ "Competitors"),
                            era = case_when(str_extract(rownames(summ), "\\d") == "1" ~ "Early",
                                            str_extract(rownames(summ), "\\d") == "2" ~ "Middle",
                                            str_extract(rownames(summ), "\\d") == "3" ~ "Late",
                                            .ptype=factor( levels=c("Early", "Middle", "Late"))))

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
summ.dens <- bind_rows(dens.l, .id="par") 
dens.df.st.2c <- data.frame(summ.dens,
                            Ocean.Region2 = sock.info$Ocean.Region2[match(summ.dens$stock, sock.info$Stock)],
                            era = case_when(
                              str_extract(summ.dens$par, "\\d") == "1" ~ "Early",
                              str_extract(summ.dens$par, "\\d") == "2" ~ "Middle",
                              str_extract(summ.dens$par, "\\d") == "3" ~ "Late",
                              .ptype=factor( levels=c("Early", "Middle", "Late"))),
                            var = stringr::str_extract(summ.dens$par, "\\D+"),
                            varnam = case_when(
                              str_extract(summ.dens$par, "\\D+") == "gamma" ~ "SST",
                              str_extract(summ.dens$par, "\\D+") == "kappa" ~ "Competitors"))

# Mean (regional) density dataframe 
post <- rstan::extract(era.2c, pars=c(paste0("mu_gamma", c(1:3)),
                                      paste0("mu_kappa", c(1:3))))
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
dens.df.reg.2c <- data.frame(summ.dens,
                             era = case_when(
                               str_extract(summ.dens$par, "\\d") == "1" ~ "Early",
                               str_extract(summ.dens$par, "\\d") == "2" ~ "Middle",
                               str_extract(summ.dens$par, "\\d") == "3" ~ "Late",
                              .ptype=factor( levels=c("Early", "Middle", "Late"))),
                             var = stringr::str_extract(summ.dens$par, "\\D+"),
                             varnam = case_when(
                               str_extract(summ.dens$par, "\\D+") == "mu_gamma" ~ "SST",
                               str_extract(summ.dens$par, "\\D+") == "mu_kappa" ~ "Competitors")
)

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
pdf("./figures/dyn/hbm_inf/eras_coef_dot_c2_comb.pdf", width=5, height=8)
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
pdf("./figures/dyn/hbm_inf/eras_coef_dot_2c_panel.pdf", width=8, height=8)
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
                   var = stringr::str_extract(rownames(summ), "[a-z]+"),
                   varnam = case_when(grepl("^gamma", rownames(summ)) ~ "SST",
                                      grepl("^kappa", rownames(summ)) ~ "Competitors")
                   )

# Summarized dataframe (regional-level)
# gamma/kappa are series-specific; no mu output. Summarize stocks instead
df.dyn.reg.2c <- dplyr::summarize(df.dyn.st.2c, 
                           reg_mean=mean(mu, na.rm=T), 
                           lower_10=quantile(mu, 0.1), 
                           upper_90=quantile(mu, 0.9), 
                           .by=c(Ocean.Region2, BY, varnam))
df.dyn.reg.2c[which(df.dyn.reg.2c$Ocean.Region2=="BS" & 
                      df.dyn.reg.2c$BY<1965), c("reg_mean", "lower_10", "upper_90")] <- NA # manually remove some sparse years

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
pdf("./figures/dyn/hbm_inf/dyn_stack_2c.pdf", width=8, height=10)
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
pdf("./figures/dyn/hbm_inf/dyn_2c_grouped.pdf")  
print(g)
dev.off()




## Single covariate models ----------------------------------------------------


### ----- Eras model: Data  

lst.era <- list(era.sst, era.comp)
probs <- c(0.025, 0.05, 0.10, 0.50, 0.90, 0.95, 0.975)


# Stock-specific dataframes
dfl.era.st <- lapply(seq_along(lst.era), function(i){
  summ <- rstan::summary(lst.era[[i]], pars=c(paste0("gamma", 1:3)), probs=probs)[[1]]
  df <- data.frame(Stock = rep(sock.info$Stock, 3),
                   Ocean.Region2 = rep(sock.info$Ocean.Region2, 3),
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


# Summarized dataframes (regional-level)
reg_start <- sock.info$Stock[match(unique(sock.info$Ocean.Region2), sock.info$Ocean.Region2)]
reg_end <- c(sock.info$Stock[match(unique(sock.info$Ocean.Region2), sock.info$Ocean.Region2)-1], 
             sock.info$Stock[nrow(sock.info)])

dfl.era.reg <- lapply(seq_along(lst.era), function(i){
  summ <- rstan::summary(lst.era[[i]], pars=c(paste0("mu_gamma", 1:3)), probs=probs)[[1]]
  df <- data.frame(Ocean.Region2 = rep(unique(sock.info$Ocean.Region2), 3),
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


# Density dataframe - stock-level
dens.dfl.st <- list()
for(i in seq_along(lst.era)) {
  post <- rstan::extract(lst.era[[i]], pars=c(paste0("gamma", c(1:3))))
  dens.l <- lapply(post, function(x){
    dens.out <- col_density(x, plot.it=F)
    dens.list <- lapply(dens.out, adply, .margins=c(1,2))
    dens.df.sst <- join(dens.list$x, dens.list$y, by=c("X1", "X2"))
    names(dens.df.sst) <- c("n", "stock", "gamma.x", "gamma.dens")
    dens.df.sst <- mutate(dens.df.sst, stock=levels(sock$Stock)[as.numeric(stock)] ) 
    return(dens.df.sst)
  } )
  summ.dens <- bind_rows(dens.l, .id="par")
  dens.df.st <- data.frame(summ.dens,
                           Ocean.Region2 = sock.info$Ocean.Region2[
                             match(summ.dens$stock, sock.info$Stock)],
                           era = case_when(
                             str_extract(summ.dens$par, "\\d") == "1" ~ "Early",
                             str_extract(summ.dens$par, "\\d") == "2" ~ "Middle",
                             str_extract(summ.dens$par, "\\d") == "3" ~ "Late",
                             .ptype=factor( levels=c("Early", "Middle", "Late")))
                            )
  dens.dfl.st[[i]] <- dens.df.st 
}
names(dens.dfl.st) <- c("SST", "Competitors")

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



### --- Eras model : Figures

## SST 

# Caterpillar plot (combined) - sst
g <- ggplot(dfl.era.st$SST) +
  geom_vline(xintercept = 0, color = "grey50", linetype = 2, linewidth = 0.25) +
  geom_point(aes(x = mu, y = Stock, color = interaction(Ocean.Region2, era), shape = Ocean.Region2)) +
  geom_segment(data = dfl.era.reg$SST, aes(y = ystart, yend = yend, x = reg_mean, xend=reg_mean,
                                 color = interaction(Ocean.Region2, era)), linewidth = 0.25) +
  geom_rect(data = dfl.era.reg$SST, aes(xmin = lower_10, xmax = upper_90, ymin = ystart,
                         ymax = yend, fill = interaction(Ocean.Region2, era)), alpha=0) +
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
pdf("./figures/dyn/hbm_inf/eras_coef_dot_sst_comb.pdf", width=5, height=8)
print(g)
dev.off()

# Caterpillar plot (panel) - sst
g <- ggplot(dfl.era.st$SST) +
  geom_vline(xintercept = 0, color = "grey50", linetype = 2, linewidth = 0.25) +
  geom_point(aes(x = mu, y = Stock, color = Ocean.Region2, shape = Ocean.Region2)) +
  geom_segment(aes(y = Stock, yend = Stock, x = lower_10, xend = upper_90,
                   color = Ocean.Region2), linewidth = 0.25) +
  geom_segment(data = dfl.era.reg$SST, aes(y = ystart, yend = yend, x = reg_mean, xend=reg_mean,
                                           color = Ocean.Region2), linewidth = 0.25) +
  geom_rect(data = dfl.era.reg$SST, aes(xmin = lower_10, xmax = upper_90, ymin = ystart,
                                        ymax = yend, fill = Ocean.Region2), alpha=0.2) +
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
pdf("./figures/dyn/hbm_inf/eras_coef_dot_sst_panel.pdf", width=5, height=8)
print(g)
dev.off()

# Posterior density plot - sst
g <- ggplot(dens.dfl.st$SST) + 
  geom_vline(xintercept = 0, color = "grey50", linetype = 2, linewidth = 0.25) +
  geom_path(aes(x=gamma.x, y=gamma.dens, group=stock, col=Ocean.Region2), alpha=0.2) + 
  geom_path(data=dens.dfl.reg$SST, aes(x=x, y=dens, col=Ocean.Region2), alpha=0.85, linewidth=1) +
  facet_wrap(~era, nrow=3) +
  scale_colour_manual(values=col.region) +
  coord_cartesian(xlim=c(-1, 1)) +
  theme_minimal() + labs(x="SST effect", y="", col="Ocean Region") +
  theme(axis.text.y=element_blank())
pdf("./figures/dyn/hbm_inf/eras_sst_dens.pdf")
print(g)
dev.off()


## Competitors

# Caterpillar plot (combined) - comp
g <- ggplot(dfl.era.st$Competitors) +
  geom_vline(xintercept = 0, color = "grey50", linetype = 2, linewidth = 0.25) +
  geom_point(aes(x = mu, y = Stock, color = interaction(Ocean.Region2, era), shape = Ocean.Region2)) +
  #geom_segment(aes(y = Stock, yend = Stock, x = `2.5%`, xend = `97.5%`,
  #                 color = Ocean.Region2), linewidth = 0.25) +
  geom_segment(data = dfl.era.reg$Competitors, aes(y = ystart, yend = yend, x = reg_mean, xend=reg_mean,
                                           color = interaction(Ocean.Region2, era)), linewidth = 0.25) +
  geom_rect(data = dfl.era.reg$Competitors, aes(xmin = lower_10, xmax = upper_90, ymin = ystart,
                                        ymax = yend, fill = interaction(Ocean.Region2, era)), alpha=0) +
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
pdf("./figures/dyn/hbm_inf/eras_coef_dot_comp_comb.pdf", width=5, height=8)
print(g)
dev.off()


# Caterpillar plot (panel) - Competitors
g <- ggplot(dfl.era.st$Competitors) +
  geom_vline(xintercept = 0, color = "grey50", linetype = 2, linewidth = 0.25) +
  geom_point(aes(x = mu, y = Stock, color = Ocean.Region2, shape = Ocean.Region2)) +
  geom_segment(aes(y = Stock, yend = Stock, x = lower_10, xend = upper_90,
                   color = Ocean.Region2), linewidth = 0.25) +
  geom_segment(data = dfl.era.reg$Competitors, aes(y = ystart, yend = yend, x = reg_mean, xend=reg_mean,
                                           color = Ocean.Region2), linewidth = 0.25) +
  geom_rect(data = dfl.era.reg$Competitors, aes(xmin = lower_10, xmax = upper_90, ymin = ystart,
                                        ymax = yend, fill = Ocean.Region2), alpha=0.2) +
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
pdf("./figures/dyn/hbm_inf/eras_coef_dot_comp_panel.pdf", width=5, height=8)
print(g)
dev.off()


# Posterior density plot - comp
g <- ggplot(dens.dfl.st$Competitors) + 
  geom_vline(xintercept = 0, color = "grey50", linetype = 2, linewidth = 0.25) +
  geom_path(aes(x=gamma.x, y=gamma.dens, group=stock, col=Ocean.Region2), alpha=0.2) + 
  geom_path(data=dens.dfl.reg$Competitors, aes(x=x, y=dens, col=Ocean.Region2), alpha=0.85, linewidth=1) +
  facet_wrap(~era, nrow=3) +
  scale_colour_manual(values=col.region) +
  coord_cartesian(xlim=c(-1, 1)) +
  theme_minimal() + labs(x="Competitor effect", y="", col="Ocean Region") +
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
                   reg_mean_gamma=mean(mu_gamma, na.rm=T), 
                   gamma_10=quantile(mu_gamma, 0.1), 
                   gamma_90=quantile(mu_gamma, 0.9), 
                   .by=c(Ocean.Region2, BY))
  reg_df[which(reg_df$Ocean.Region2=="BS" & 
                        reg_df$BY<1965), c("reg_mean_gamma", "gamma_10", "gamma_90")] <- NA # manually remove some sparse years
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
  theme_minimal() + labs(x= "Brood Year", y="SST effect") +
  theme(legend.position="none")
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
  theme_minimal() + labs(x= "Brood Year", y="Competitor effect") +
  theme(legend.position="none")
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


