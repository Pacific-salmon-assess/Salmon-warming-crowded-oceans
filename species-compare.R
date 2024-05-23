## comparison plots across species##

# 1) Boxplots of effects

# Load era model fits

#sockeye
load(here('output', 'models', 'dyn', 'sockeye', 'hbm_era_2c_3sub.RData'), verbose=T) # eras
sock.box <- era_density_df(era.2c.3sub, par=c("gamma", "kappa"), mu=T, info=sock.info)
load(here('output', 'models', 'stat', 'sockeye', 'stat_sub.RData'), verbose=T) # stationary 
s.dens <- hb05_density_df(stat_sub, ocean.regions = 4, info_master=sock.info.sub, data_master=sock.sub)$region
s.dens <- s.dens %>% filter(var %in% c("SST", "Comp"))
s.df <- data.frame(n=factor(s.dens$n),
                   Ocean.Region2=rep(c("WC", "SEAK", "GOA", "BS"), each=max(s.dens$n)),
                   x=s.dens$x,
                   dens=s.dens$y,
                   era="All",
                   varnam=gsub("Comp", "Competitors", s.dens$var),
                   sp.id="sockeye")
#pink
load(here('output', 'models', 'dyn', 'pink', 'hbm_era_2c.RData'), verbose=T) # eras
pink.box <- era_density_df(era.2c, par=c("gamma", "kappa"), mu=T, info=pink.info)
pink.box <- dplyr::filter(pink.box, !(Ocean.Region2 %in% c("BS", "GOA") & era=="Late"))
load(here('output', 'models', 'stat', 'pink', 'stat_a.RData'), verbose=T) # stationary 
p.dens <- hb05_density_df(stat_a, ocean.regions = 4, info_master=pink.info, data_master=pink)$region
p.dens <- p.dens %>% filter(var %in% c("SST", "Comp"))
p.df <- data.frame(n=factor(p.dens$n),
                   Ocean.Region2=rep(c("WC", "SEAK", "GOA", "BS"), each=max(p.dens$n)),
                   x=p.dens$x,
                   dens=p.dens$y,
                   era="All",
                   varnam=gsub("Comp", "Competitors", p.dens$var),
                   sp.id="pink")

#pink-odd
load(here('output', 'models', 'dyn', 'pink', 'hbm_era_2c_odd.RData'), verbose=T) # eras
pink.box <- era_density_df(era.2c.odd, par=c("gamma", "kappa"), mu=T, info=pink.info[grep("-Odd", pink.info$Stock),])
pink.box.odd <- dplyr::filter(pink.box, !(Ocean.Region2 %in% c("BS", "GOA") & era=="Late"))
load(here('output', 'models', 'stat', 'pink', 'stat_a_odd.RData'), verbose=T) # stationary 
p.dens <- hb05_density_df(stat_a_odd, ocean.regions = 4, data_master=data_master, info_master=info_master)$region
p.dens <- p.dens %>% filter(var %in% c("SST", "Comp"))
p.odd.df <- data.frame(n=factor(p.dens$n),
                   Ocean.Region2=rep(c("WC", "SEAK", "GOA", "BS"), each=max(p.dens$n)),
                   x=p.dens$x,
                   dens=p.dens$y,
                   era="All",
                   varnam=gsub("Comp", "Competitors", p.dens$var),
                   sp.id="pink-odd")

#pink-even
load(here('output', 'models', 'dyn', 'pink', 'hbm_era_2c_even.RData'), verbose=T) # eras
pink.box <- era_density_df(era.2c.even, par=c("gamma", "kappa"), mu=T, info=pink.info[grep("-Even", pink.info$Stock),])
pink.box.even <- dplyr::filter(pink.box, !(Ocean.Region2 %in% c("BS", "GOA") & era=="Late"))
load(here('output', 'models', 'stat', 'pink', 'stat_a_even.RData'), verbose=T) # stationary 
p.dens <- hb05_density_df(stat_a_even, ocean.regions = 4, data_master=data_master, info_master=info_master)$region
p.dens <- p.dens %>% filter(var %in% c("SST", "Comp"))
p.even.df <- data.frame(n=factor(p.dens$n),
                   Ocean.Region2=rep(c("WC", "SEAK", "GOA", "BS"), each=max(p.dens$n)),
                   x=p.dens$x,
                   dens=p.dens$y,
                   era="All",
                   varnam=gsub("Comp", "Competitors", p.dens$var),
                   sp.id="pink-even")


#chum
load(here('output', 'models', 'dyn', 'chum', 'hbm_era_2c.RData'), verbose=T) # eras
chum.info$Ocean.Region <- gsub("SEAK", "WC", chum.info$Ocean.Region2)
chum.box <- era_density_df(era.2c, par=c("gamma", "kappa"), mu=T, info=chum.info, region.var="Ocean.Region")
chum.box <- dplyr::filter(chum.box, !(Ocean.Region %in% c("BS", "GOA") & era=="Late"))
chum.box <- dplyr::rename(chum.box, Ocean.Region2=Ocean.Region)
load(here('output', 'models', 'stat', 'chum', 'stat_a.RData'), verbose=T) # stationary 
c.dens <- hb05_density_df(stat_a, ocean.regions = 3, info_master=chum.info, data_master=chum)$region
c.dens <- c.dens %>% filter(var %in% c("SST", "Comp"), !is.na(region))
c.df <- data.frame(n=factor(c.dens$n),
                   Ocean.Region2=rep(c("WC", "GOA", "BS"), each=max(c.dens$n)),
                   x=c.dens$x,
                   dens=c.dens$y,
                   era="All",
                   varnam=gsub("Comp", "Competitors", c.dens$var),
                   sp.id="chum")

box.lst <- list(sock.box, pink.box.odd, pink.box.even, chum.box)
names(box.lst) <- c("sockeye", "pink-odd", "pink-even", "chum")
box.eras <- bind_rows(box.lst, .id="sp.id")

box.stat <- bind_rows(s.df, p.odd.df, p.even.df, c.df)
box.stat$x <- log((box.stat$x/100)+1)
box <- bind_rows(box.eras, box.stat)



# boxplot with both covariates and stationary estimates
col.box <- c(RColorBrewer::brewer.pal(n=4, "BuGn")[2:4], "#B0B0B0")
full.box <-  ggplot(box) +
  geom_hline(aes(yintercept=0), linetype="dashed", colour="gray50") +
  geom_boxplot(aes(x=sp.id, y=x, fill=factor(era, levels=c("Early", "Middle", "Late", "All"))), position=position_dodge2(preserve="single")) +
  scale_y_continuous(limits=c(-1, 1)) + 
  facet_grid(rows=vars(Ocean.Region2), cols=vars(factor(varnam, levels=c("SST","Competitors")))) +
  labs(y="Covariate Effect", x="", fill="") +
  scale_fill_manual(values=col.box) +
  theme_sleek()
png(here('figures', 'spp-explore', 'boxplot-2c-compare-evenodd.png'), width=650, height=600)
print(full.box)
dev.off()


# Boxplots with overlaid RW lines - works but looks shitty

rw.mov$sp.id <- rw.mov$spp
box %>% filter(varnam=="SST") %>% mutate(BY = case_when(era=="Early" ~ 1980,
                                                        era=="Middle" ~ 2000,
                                                        era=="Late" ~ 2015,
                                                        era=="All" ~ 2030)) %>%
  ggplot() +
  geom_hline(aes(yintercept=0), linetype="dashed", colour="gray50") +
  geom_line(data=filter(rw.mov, varnam=="SST"), aes(x=BY, y=mov_avg), linewidth=1, col="grey70", alpha=0.8) +
  geom_boxplot(aes(x=BY, y=x, fill=factor(era, levels=c("Early", "Middle", "Late", "All"))), position=position_dodge2(preserve="single"), alpha=0.7) +
  scale_y_continuous(limits=c(-1, 1)) + 
  facet_grid(rows=vars(Ocean.Region2), cols=vars(sp.id)) +
  labs(y="Covariate Effect", x="", fill="") +
  scale_fill_manual(values=col.box) +
  theme_sleek()




# Figures
sst.box <- box %>% filter(par %in% paste0("mu_gamma", 1:3)) %>%
                  ggplot() + geom_boxplot(aes(x=sp.id, y=x, fill=era)) + 
                  coord_cartesian(ylim=c(-1, 1)) + 
                  geom_hline(aes(yintercept=0), linetype="dashed", colour="gray50") +
                  facet_wrap(vars(Ocean.Region2)) +
                  labs(y="SST Effect") +
                  theme_minimal()
png(here('figures', 'spp-explore', 'sst-compare.png'))
print(sst.box)
dev.off()


comp.box <- box %>% filter(par %in% paste0("mu_kappa", 1:3)) %>%
  ggplot() + geom_boxplot(aes(x=sp.id, y=x, fill=era)) + 
  coord_cartesian(ylim=c(-1, 1)) + 
  geom_hline(aes(yintercept=0), linetype="dashed", colour="gray50") +
  scale_fill_manual(values=col.box) +
  facet_wrap(vars(Ocean.Region2)) +
  labs(y="Competition Effect") +
  theme_sleek()
png(here('figures', 'spp-explore', 'comp-compare.png'))
print(comp.box)
dev.off()



# 2) Random walk smoothed - all species


# a) sockeye

stk.sub <- sock.info$Stock.ID[sock.info$yr_start < 1985 & sock.info$yr_end >= 2014]
sock.sub <- sock %>% filter(Stock.ID %in% stk.sub)
sock.info.sub <- sock.info %>% filter(Stock.ID %in% stk.sub)


# load sockeye RW
load("./output/models/dyn/sockeye/hbm_dyn_2c_sub.RData", verbose=T)
summ <- rstan::summary(dyn.2c.sub, pars = c("gamma", "kappa"), probs = NULL)$summary
rw.stk <- data.frame(Stock = sock.sub$Stock,
                           Ocean.Region2 = sock.sub$Ocean.Region2,
                           BY = sock.sub$BY,
                           mu = summ[, "mean"],
                           se = summ[, "se_mean"],
                           var = str_extract(rownames(summ), "[a-z]+"),
                           varnam = case_when(grepl("^gamma", rownames(summ)) ~ "SST",
                                              grepl("^kappa", rownames(summ)) ~ "Competitors")
)



# Summarize at regional lvl
rw.reg <- dplyr::summarize(rw.stk, 
                            reg_mean=mean(mu, na.rm=T), 
                            n_stk=n_distinct(Stock),
                            .by=c(Ocean.Region2, BY, varnam))
rw.reg <- ddply(rw.reg, .(Ocean.Region2), dplyr::filter, n_stk >= 4)

# apply moving average
rw.reg.avg.sock <- moving_average_df(rw.reg, "reg_mean", lag=2)


# b) pink

load("./output/models/dyn/pink/hbm_dyn_2c.RData")

summ <- rstan::summary(dyn.2c, pars = c("gamma", "kappa"), probs = NULL)$summary
rw.stk.pink <- data.frame(Stock = pink$Stock,
                           Ocean.Region2 = pink$Ocean.Region2,
                           BY = pink$BY,
                           mu = summ[, "mean"],
                           se = summ[, "se_mean"],
                           var = str_extract(rownames(summ), "[a-z]+"),
                           varnam = case_when(grepl("^gamma", rownames(summ)) ~ "SST",
                                              grepl("^kappa", rownames(summ)) ~ "Competitors"),
                           even_odd = ifelse(grepl("Even$", pink$Stock), "Even", "Odd")
)

# Summarize at regional lvl
rw.reg.pink <- dplyr::summarize(rw.stk.pink,
                                  reg_mean=mean(mu, na.rm=T), 
                                  n_stk=n_distinct(Stock),
                                  .by=c(Ocean.Region2, BY, varnam, even_odd))
rw.reg.pink <- ddply(rw.reg.pink, .(Ocean.Region2), dplyr::filter, n_stk >= 4)

rw.reg.avg.pink.even <- moving_average_df(subset(rw.reg.pink, even_odd=="Even"), "reg_mean", lag=2)
rw.reg.avg.pink.odd <- moving_average_df(subset(rw.reg.pink, even_odd=="Odd"), "reg_mean", lag=2)

# c) chum 

load("./output/models/dyn/chum/hbm_dyn_2c.RData")

summ <- rstan::summary(dyn.2c, pars = c("gamma", "kappa"), probs = NULL)$summary
rw.stk.chum <- data.frame(Stock = chum$Stock,
                          Ocean.Region2 = chum$Ocean.Region2,
                          BY = chum$BY,
                          mu = summ[, "mean"],
                          se = summ[, "se_mean"],
                          var = str_extract(rownames(summ), "[a-z]+"),
                          varnam = case_when(grepl("^gamma", rownames(summ)) ~ "SST",
                                             grepl("^kappa", rownames(summ)) ~ "Competitors")
)

# Summarize at regional lvl
rw.reg.chum <- dplyr::summarize(rw.stk.chum,
                                reg_mean=mean(mu, na.rm=T), 
                                n_stk=n_distinct(Stock),
                                .by=c(Ocean.Region2, BY, varnam))
rw.reg.chum <- ddply(rw.reg.chum, .(Ocean.Region2), dplyr::filter, n_stk >= 4)

rw.reg.avg.chum <- moving_average_df(rw.reg.chum, "reg_mean", lag=2)

# bind all the moving average dfs
rw.mov <- bind_rows(rw.reg.avg.sock, rw.reg.avg.pink.even, rw.reg.avg.pink.odd, rw.reg.avg.chum, .id="spp") %>% 
  mutate(spp=case_when(spp == 1 ~ "sockeye",
                       spp == 2 ~ "pink-even",
                       spp == 3 ~ "pink-odd",
                       spp == 4 ~ "chum")) %>% 
  tidyr::replace_na(list(even_odd = "Both"))

# era coefficients to overlay # not in use
eras.all.spp <- box %>% summarize(med_x = median(x), .by=c("sp.id", "par", "Ocean.Region2", "era", "varnam")) %>%
                                    mutate(xstart = case_when(era=="Early" ~ 1950,
                                                             era=="Middle" ~ 1989,
                                                             era=="Late" ~ 2011),
                                          xend = case_when(era=="Early" ~ 1988,
                                                           era=="Middle" ~ 2010,
                                                           era=="Late" ~ 2019),
                                          xmid = case_when(era=="Early" ~ 1969,
                                                           era=="Middle" ~ 1999,
                                                           era=="Late" ~ 2013)) %>%
                        rename(spp=sp.id)


# Figure
lty.brd <- c("chum"="solid", "pink-even"="dotdash", "pink-odd"="dashed", "sockeye"="solid")
lwdth.brd <- c("chum"=0.8, "pink-even"=0.65, "pink-odd"=0.9, "sockeye"=0.8)
sp.col.rw <- c("seagreen4", "palevioletred3", "palevioletred3", "orangered4") # species colours 
names(sp.col.rw) <- c("Chum", "Pink-Even", "Pink-Odd", "Sockeye") # name
dyn_compare <-  
  ggplot(rw.mov) +
  geom_vline(xintercept=c(1989,2011), color = "grey80", linetype = 1, linewidth = 0.25, alpha=0.8) +
  geom_hline(yintercept = 0, color = "grey50", linetype = 2, linewidth = 0.25) +
  geom_line(aes(x=BY, y=mov_avg, col=spp, linewidth=spp), alpha=0.8) +
  #geom_point(data=eras.all.spp, aes(x=xmid, y=med_x, col=Ocean.Region2), alpha=0.6) +
  #geom_segment(data=eras.all.spp, aes(x=xstart, xend=xend, y=med_x, yend=med_x, col=Ocean.Region2), alpha=0.2) + 
  scale_colour_manual(values=unname(sp.col.rw)) +
  #scale_linetype_manual(values=lty.brd) +
  scale_linewidth_manual(values=lwdth.brd) +
  facet_grid(cols=vars(factor(varnam, levels=c("SST", "Competitors"))), rows=vars(Ocean.Region2)) + 
  scale_y_continuous(limits=c(-0.5, 0.75), breaks=seq(-1,1,0.5)) + 
  theme_sleek() + theme(aspect.ratio = 0.85,
                        legend.position="right") +
  labs(x="Brood Year", y="Mean covariate effects \n (moving average)", col="Species", linewidth="Species")
  
png(here('figures', 'dyn', 'dyn-compare-v2.png'), width=900, height=1100, res=150)
print(dyn_compare)
dev.off()


overlay2 <- rw.mov %>% filter(spp=="sockeye", ocean_region_lab=="West Coast")  %>% ggplot() +
  #geom_hline(yintercept = 0, color = "grey50", linetype = 2, linewidth = 0.25) +
  geom_line(aes(x=BY, y=mov_avg, col=Ocean.Region2, linewidth=spp), alpha=0.8) +
  #geom_point(data=eras.all.spp, aes(x=xmid, y=med_x, col=Ocean.Region2), alpha=0.6) +
  geom_segment(data=filter(eras.all.spp, spp=="sockeye", Ocean.Region2=="WC"), aes(x=xstart, xend=xend, y=med_x, yend=med_x), col="grey30", alpha=0.8) + 
  scale_colour_manual(values=unname(col.region)) +
  #scale_linetype_manual(values=lty.brd) +
  scale_linewidth_manual(values=lwdth.brd) +
  facet_grid(cols=vars(factor(varnam, levels=c("SST", "Competitors"))), rows=vars(Ocean.Region2)) + 
  #scale_y_continuous(limits=c(-0.5, 0.75), breaks=seq(-1,1,0.5)) + 
  theme_sleek() + theme(aspect.ratio = 0.85,
                        legend.position="none") +
  labs(x="Brood Year", y="Mean covariate effects", col="Species", linewidth="Species")

png(here('figures', 'spp-explore', 'sock_rw_eras_trim.png'), width=548*3, height=304*3, res=72*3)
print(overlay2)
dev.off()


# Comparison table of ERA results

# load sockeye
load(here('output', 'models', 'dyn', 'sockeye', 'hbm_era_2c_3sub.RData'), verbose=T) 
stbl <- rstan::summary(era.2c.3sub, par=c("mu_alpha", paste0("mu_gamma", 1:3), paste0("mu_kappa", 1:3)))$summary
stbl <- as.data.frame(stbl)
stbl$species <- "sockeye"
stbl$id <- rownames(stbl)

# got halfway through this df and left it

# load pink-odd
load(here('output', 'models', 'dyn', 'pink', 'hbm_era_2c_odd.RData'), verbose=T) 
potbl <- rstan::summary(era.2c.odd, par=c("mu_alpha", paste0("mu_gamma", 1:3), paste0("mu_kappa", 1:3)))$summary
potbl <- as.data.frame(potbl)
potbl$species <- "pink-odd"
potbl$id <- rownames(potbl)

# load pink-even
load(here('output', 'models', 'dyn', 'pink', 'hbm_era_2c_even.RData'), verbose=T) 
petbl <- rstan::summary(era.2c.even, par=c("mu_alpha", paste0("mu_gamma", 1:3), paste0("mu_kappa", 1:3)))$summary
petbl <- as.data.frame(petbl)
petbl$species <- "pink-even"
petbl$id <- rownames(petbl)

# load chum
load(here('output', 'models', 'dyn', 'chum', 'hbm_era_2c.RData'), verbose=T) 
ctbl <- rstan::summary(era.2c, par=c("mu_alpha", paste0("mu_gamma", 1:3), paste0("mu_kappa", 1:3)))$summary
ctbl <- as.data.frame(ctbl)
ctbl$species <- "chum"
ctbl$id <- rownames(ctbl)

full_tbl <- bind_rows(stbl, potbl, petbl, ctbl)
save(full_tbl,  file="./output/eras_full_compare.RData")
View(full_tbl)

# Get stationary estimates to compare to eras
statmeans <- list(NA)
load(here('output', 'models', 'stat', 'sockeye', 'stat_sub.RData'), verbose=T)
statmeans[1] <- rstan::summary(stat_sub, par=c("mu_gamma"))$summary[1,1]
load(here('output', 'models', 'stat', 'chum', 'stat_a.RData'), verbose=T)
statmeans[2] <- rstan::summary(stat_a, par=c("mu_gamma"))$summary[1,1]
load(here('output', 'models', 'stat', 'pink', 'stat_a_even.RData'), verbose=T)
statmeans[3] <- rstan::summary(stat_a_even, par=c("mu_gamma"))$summary[1,1]
load(here('output', 'models', 'stat', 'pink', 'stat_a_odd.RData'), verbose=T)
statmeans[4] <- rstan::summary(stat_a_odd, par=c("mu_gamma"))$summary[1,1]

filter(full_tbl, id=="mu_gamma3[1]")[c("mean", "species")]
