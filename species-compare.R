## comparison plots across species##

#### 1) Boxplots of effects ####

# Load era  and stationary model fits

#sockeye
load(here('output', 'models', 'dyn', 'sockeye', 'hbm_era_2c.RData'), verbose=T) # eras
sock.box <- era_density_df(era.2c, par=c("gamma", "kappa"), mu=T, info=sock.info, percent.change = T)
load(here('output', 'models', 'stat', 'sockeye', 'stat_a.RData'), verbose=T) # stationary
s.dens <- hb05_density_df(stat_a, ocean.regions = 4, info=sock.info, data=sock)$region
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
pink.box <- era_density_df(era.2c, par=c("gamma", "kappa"), mu=T, info=pink.info, percent.change = T)
load(here('output', 'models', 'stat', 'pink', 'stat_a.RData'), verbose=T) # stationary
p.dens <- hb05_density_df(stat_a, ocean.regions = 4, info=pink.info, data=pink)$region
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
pink.box <- era_density_df(era.2c.odd, par=c("gamma", "kappa"), mu=T, info=pink.info[grep("-Odd", pink.info$Stock),], percent.change = T)
pink.box.odd <- dplyr::filter(pink.box, !(Ocean.Region2 %in% c("BS", "GOA") & era=="Late"))
load(here('output', 'models', 'stat', 'pink', 'stat_a_odd.RData'), verbose=T) # stationary
p.dens <- hb05_density_df(stat_a_odd, ocean.regions = 4, data=pink[grep("-Odd", pink$Stock),], info=pink.info[grep("-Odd", pink.info$Stock),])$region
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
pink.box <- era_density_df(era.2c.even, par=c("gamma", "kappa"), mu=T, info=pink.info[grep("-Even", pink.info$Stock),], percent.change = T)
pink.box.even <- dplyr::filter(pink.box, !(Ocean.Region2 %in% c("BS", "GOA") & era=="Late"))
load(here('output', 'models', 'stat', 'pink', 'stat_a_even.RData'), verbose=T) # stationary
p.dens <- hb05_density_df(stat_a_even, ocean.regions = 4, data=pink[grep("-Even", pink$Stock),], info=pink.info[grep("-Even", pink.info$Stock),])$region
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
#chum.info$Ocean.Region <- gsub("SEAK", "WC", chum.info$Ocean.Region2)
chum.box <- era_density_df(era.2c, par=c("gamma", "kappa"), mu=T, info=chum.info, region.var="Ocean.Region2", percent.change = T)
load(here('output', 'models', 'stat', 'chum', 'stat_a.RData'), verbose=T) # stationary
c.dens <- hb05_density_df(stat_a, ocean.regions = 3, info=chum.info, data=chum)$region
c.dens <- c.dens %>% filter(var %in% c("SST", "Comp"), !is.na(region))
c.df <- data.frame(n=factor(c.dens$n),
                   Ocean.Region2=rep(c("WC", "GOA", "BS"), each=max(c.dens$n)),
                   x=c.dens$x,
                   dens=c.dens$y,
                   era="All",
                   varnam=gsub("Comp", "Competitors", c.dens$var),
                   sp.id="chum")

box.lst <- list(sock.box, pink.box, chum.box)
names(box.lst) <- c("sockeye", "pink", "chum")
box.eras <- bind_rows(box.lst, .id="sp.id")

box.stat <- bind_rows(s.df, p.df, p.df, c.df)
box <- bind_rows(box.eras, box.stat)
box <- ocean_region_lab(box)


# boxplot with both covariates and stationary estimates
col.box <- c(RColorBrewer::brewer.pal(n=4, "BuGn")[2:4], "#B0B0B0")
full.box <-  ggplot(box) +
  geom_hline(aes(yintercept=0), linetype="dashed", colour="gray50") +
  geom_boxplot(aes(x=sp.id, y=x, fill=factor(era, levels=c("Early", "Middle", "Late", "All"))), position=position_dodge2(preserve="single")) +
  #scale_y_continuous(limits=c(-1, 1)) +
  facet_grid(rows=vars(Ocean.Region2), cols=vars(factor(varnam, levels=c("SST","Competitors"))), scales="free") +
  labs(y="Covariate Effect", x="", fill="") +
  scale_fill_manual(values=col.box) +
  theme_sleek()
png(here('figures', 'spp-explore', 'boxplot-2c-compare.png'), width=650, height=600)
print(full.box)
dev.off()


#### 2) Boxplots with overlaid RW lines ####

#rw.mov$sp.id <- rw.mov$spp
rw.mov <- rw.mov %>% mutate(sp.id = case_when(spp=="sockeye" ~ "sockeye",
                                              spp=="pink-even" ~ "pink",
                                              spp=="pink-odd" ~ "pink",
                                              spp=="chum" ~ "chum"))
rw.mov <- ocean_region_lab(rw.mov)

# use raw RW instead?
rw.reg.all <- ocean_region_lab(rw.reg.all) # rw.reg.all is created lower in this script

# transform both to % change
#box$x <-  (exp(box$x) - 1) * 100 # added code to do so above
rw.reg.all$reg_mean <- (exp(rw.reg.all$reg_mean) - 1) * 100
rw.reg.all$lower_10 <- (exp(rw.reg.all$lower_10) - 1) * 100
rw.reg.all$upper_90 <- (exp(rw.reg.all$upper_90) - 1) * 100


# SST plot - have to plot species separately

# pink
pink.rw.box.fig <- box %>%
  filter(varnam=="SST", era != "All", sp.id=="pink") %>% mutate(BY = case_when(era=="Early" ~ 1975,
                                                        era=="Middle" ~ 2000,
                                                        era=="Late" ~ 2015)) %>%
  mutate(regxera = paste0(ocean_region_lab, sep=".", era)) %>%
  ggplot() +
  geom_hline(aes(yintercept=0), linewidth=1, colour="grey95") +
  geom_line(data=filter(rw.reg.all, varnam=="SST", sp.id=="pink"), aes(x=BY, y=reg_mean, group=even_odd, col=ocean_region_lab), linewidth=1, alpha=0.4) +
  geom_ribbon(data=filter(rw.reg.all, varnam=="SST", sp.id=="pink"), aes(x=BY, ymin=lower_10, ymax=upper_90, group=even_odd, fill=ocean_region_lab), alpha=0.2) +
  #geom_boxplot(aes(x=BY, y=x, fill=factor(era, levels=c("Early", "Middle", "Late", "All"))), position=position_dodge2(preserve="single"), alpha=0.7) +
  geom_boxplot(aes(x=BY, y=x, fill=regxera), position=position_dodge2(preserve="single"), alpha=0.7, outlier.shape = NA) +
  #coord_cartesian(ylim=c(-100,250)) +
  #scale_y_continuous(limits=c(-100, 500)) +
  facet_grid(rows=vars(Ocean.Region2)) +
  labs(y="% Change in R/S (per \u00B0C)", x="", fill="", title="Pink") +
  scale_fill_manual(values=col.eras) +
  scale_colour_manual(values=col.region, guide="none") +
  scale_x_continuous(breaks=c(1970, 1990, 2010)) +
  theme_sleek() +
  theme(#axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.text=element_text(size=12),
        axis.title=element_text(size=15),
        strip.background = element_rect(fill="transparent", colour="transparent"),
        strip.text.y = element_blank(),
        panel.spacing.x = unit(0, unit="pt"),
        panel.spacing.y = unit(30, unit="pt"),
        panel.background = element_rect(fill="white"),
        panel.grid.major.y = element_line(colour="grey95", size=0.5),
        panel.grid.major.x = element_blank(),
        legend.position = "none",
        legend.key = element_blank(),
        plot.title = element_text(hjust=0.5, size=15, colour="grey20"))

# chum
chum.rw.box.fig <- box %>%
  filter(varnam=="SST", era != "All", sp.id=="chum") %>% mutate(BY = case_when(era=="Early" ~ 1975,
                                                                               era=="Middle" ~ 2000,
                                                                               era=="Late" ~ 2015)) %>%
  mutate(regxera = paste0(ocean_region_lab, sep=".", era)) %>%
  ggplot() +
  geom_hline(aes(yintercept=0), linewidth=1, colour="grey95") +
  geom_line(data=filter(rw.reg.all, varnam=="SST", sp.id=="chum"), aes(x=BY, y=reg_mean, group=even_odd, col=ocean_region_lab), linewidth=1, alpha=0.4) +
  geom_ribbon(data=filter(rw.reg.all, varnam=="SST", sp.id=="chum"), aes(x=BY, ymin=lower_10, ymax=upper_90, group=even_odd, fill=ocean_region_lab), alpha=0.2) +
  #geom_boxplot(aes(x=BY, y=x, fill=factor(era, levels=c("Early", "Middle", "Late", "All"))), position=position_dodge2(preserve="single"), alpha=0.7) +
  geom_boxplot(aes(x=BY, y=x, fill=regxera), position=position_dodge2(preserve="single"), alpha=0.7, outlier.shape = NA) +
  #coord_cartesian(ylim=c(-100,100)) +
  #scale_y_continuous(limits=c(-100, 500)) +
  facet_grid(rows=vars(Ocean.Region2)) +
  labs(y="% Change in R/S (per \u00B0C)", x="", fill="", title="Chum") +
  scale_fill_manual(values=col.eras) +
  scale_colour_manual(values=col.region, guide="none") +
  scale_x_continuous(breaks=c(1970, 1990, 2010)) + theme_sleek() +
  theme(#axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.text=element_text(size=12),
        axis.title=element_text(size=15),
        strip.background = element_rect(fill="transparent", colour="transparent"),
        strip.text.y = element_blank(),
        panel.spacing.x = unit(0, unit="pt"),
        panel.spacing.y = unit(30, unit="pt"),
        panel.background = element_rect(fill="white"),
        panel.grid.major.y = element_line(colour="grey95", size=0.5),
        panel.grid.major.x = element_blank(),
        legend.position = "none",
        legend.key = element_blank(),
        plot.title = element_text(hjust=0.5, size=15, colour="grey20"))

# Sockeye
sock.rw.box.fig <- box %>%
  filter(varnam=="SST", era != "All", sp.id=="sockeye") %>% mutate(BY = case_when(era=="Early" ~ 1975,
                                                                               era=="Middle" ~ 2000,
                                                                               era=="Late" ~ 2015)) %>%
  mutate(regxera = paste0(ocean_region_lab, sep=".", era)) %>%
  ggplot() +
  geom_hline(aes(yintercept=0), linewidth=1, colour="grey95") +
  geom_line(data=filter(rw.reg.all, varnam=="SST", sp.id=="sockeye"), aes(x=BY, y=reg_mean, group=even_odd, col=ocean_region_lab), linewidth=1, alpha=0.4) +
  geom_ribbon(data=filter(rw.reg.all, varnam=="SST", sp.id=="sockeye"), aes(x=BY, ymin=lower_10, ymax=upper_90, group=even_odd, fill=ocean_region_lab), alpha=0.2) +
  #geom_boxplot(aes(x=BY, y=x, fill=factor(era, levels=c("Early", "Middle", "Late", "All"))), position=position_dodge2(preserve="single"), alpha=0.7) +
  geom_boxplot(aes(x=BY, y=x, fill=regxera), position=position_dodge2(preserve="single"), alpha=0.7, outlier.shape = NA) +
  #coord_cartesian(ylim=c(-100,100)) +
  #scale_y_continuous(limits=c(-100, 500)) +
  facet_grid(rows=vars(Ocean.Region2)) +
  labs(y="% Change in R/S (per \u00B0C)", x="", fill="", title="Sockeye") +
  scale_fill_manual(values=col.eras) +
  scale_colour_manual(values=col.region, guide="none") +
  scale_x_continuous(breaks=c(1970, 1990, 2010)) + theme_sleek() +
  theme(#axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.text=element_text(size=12),
        axis.title=element_text(size=15),
        strip.background = element_rect(fill="transparent", colour="transparent"),
        strip.text.y = element_blank(),
        panel.spacing.x = unit(0, unit="pt"),
        panel.spacing.y = unit(30, unit="pt"),
        panel.background = element_rect(fill="white"),
        panel.grid.major.y = element_line(colour="grey95", size=0.5),
        panel.grid.major.x = element_blank(),
        legend.position = "none",
        legend.key = element_blank(),
        plot.title = element_text(hjust=0.5, size=15, colour="grey20"))

# Manually make legend (this is repeated below somewhere)
lgnd.dat <- data.frame(name=names(col.eras),
                       region=str_extract(names(col.eras), "[^.]+"),
                       era=str_extract(names(col.eras), "\\w*$"))

lgnd1 <- lgnd.dat %>% mutate(era_lab = paste(era, case_when(era=="Early" ~ "\n <1988",
                                                            era=="Middle" ~ "\n    1989-2010",
                                                            era=="Late" ~ "\n        2011-2019"))) %>%
  ggplot() + geom_tile(aes(x=region, y=1,
                           fill=name), col="white") +
  scale_fill_manual(values=col.eras, guide="none") + labs(x="", y="") +
  facet_grid(rows=vars(factor(era_lab, levels=unique(era_lab))), switch="x") +
  labs(title="Era") +
  theme(panel.background=element_rect(fill="white"),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        strip.background = element_rect(fill="white"),
        strip.text.y.right = element_text(angle=0, margin=margin(l=-15)),
        plot.title = element_text(hjust=0.5)
  )
blank <- ggplot() + theme(panel.background = element_rect(fill="white")) # blank grob to put legend on
legend <- cowplot::ggdraw(blank) + cowplot::draw_plot(lgnd1, 0,0.4,1,0.2)

# Plot each spp with the legend

# pink
p <- cowplot::plot_grid(pink.rw.box.fig, legend, ncol=2, rel_widths = c(1,0.6))
png(here('figures', 'spp-explore', 'sst-box-rw-pink.png'), height=700*2, width=400*2, res=72*2) # save pink
print(p)
dev.off()

#chum
p <- cowplot::plot_grid(chum.rw.box.fig, legend, ncol=2, rel_widths = c(1,0.6))
png(here('figures', 'spp-explore', 'sst-box-rw-chum.png'), height=700*2, width=400*2, res=72*2) # save chum
print(p)
dev.off()

#sockeye
p <- cowplot::plot_grid(sock.rw.box.fig, legend, ncol=2, rel_widths = c(1,0.6))
png(here('figures', 'spp-explore', 'sst-box-rw-sock.png'), height=700*2, width=400*2, res=72*2) # save sock
print(p)
dev.off()


## Plot them together
p <- cowplot::plot_grid(sock.rw.box.fig, chum.rw.box.fig + theme(axis.title.y=element_blank()), pink.rw.box.fig + theme(axis.title.y=element_blank()), legend, ncol=4, rel_widths = c(1,1,1,0.7))
png(here('figures', 'spp-explore', 'sst-box-rw-2026.png'), height=700*2, width=900*2, res=72*2) # save sock
print(p)
dev.off()


# Boxplot + RW for competitors

comp_rw_box <- box %>%
  filter(varnam=="Competitors", era != "All") %>% mutate(BY = case_when(era=="Early" ~ 1975,
                                                                era=="Middle" ~ 2000,
                                                                era=="Late" ~ 2015)) %>%
  mutate(regxera = paste0(ocean_region_lab, sep=".", era)) %>%
  ggplot() +
  geom_hline(aes(yintercept=0), linewidth=1, colour="grey95") +
  geom_line(data=filter(rw.mov, varnam=="Competitors"), aes(x=BY, y=mov_avg, group=spp, col=ocean_region_lab), linewidth=1, alpha=0.4) +
  #geom_boxplot(aes(x=BY, y=x, fill=factor(era, levels=c("Early", "Middle", "Late", "All"))), position=position_dodge2(preserve="single"), alpha=0.7) +
  geom_boxplot(aes(x=BY, y=x, fill=regxera), position=position_dodge2(preserve="single"), alpha=0.7) +
  scale_y_continuous(limits=c(-1, 1)) +
  facet_grid(rows=vars(factor(ocean_region_lab, levels=rev(unique(box$ocean_region_lab)))), cols=vars(sp.id), switch="x") +
  labs(y="", x="", fill="", title="Competitors") +
  scale_fill_manual(values=col.eras) +
  scale_colour_manual(values=col.region, guide="none") +
  scale_x_continuous(breaks=c(1970, 1990, 2010)) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_rect(fill="transparent", colour="transparent"),
        strip.text.y = element_text(margin=margin(l=10)),
        strip.text.x = element_text(margin=margin(b=10,t=5)),
        panel.spacing.x = unit(0, unit="pt"),
        panel.spacing.y = unit(30, unit="pt"),
        panel.background = element_rect(fill="white"),
        panel.grid.major.y = element_line(colour="grey95", size=0.5),
        panel.grid.major.x = element_blank(),
        legend.position = "none",
        legend.key = element_blank(),
        plot.margin = margin(l=0),
        plot.title = element_text(hjust=0.5, size=10, colour="grey20", margin=margin(t=5,b=5)))
# still want to get the lines to change colour. easiest way would be to attach a colour scale based on year to the data frame rw.mov

# Put SST and Comp together
main_box <- cowplot::plot_grid(sst_rw_box, comp_rw_box)


# Make a legend manually
lgnd.dat <- data.frame(name=names(col.eras),
                     region=str_extract(names(col.eras), "[^.]+"),
                     era=str_extract(names(col.eras), "\\w*$"))

lgnd1 <- lgnd.dat %>% mutate(era_lab = paste(era, case_when(era=="Early" ~ "\n <1988",
                                                            era=="Middle" ~ "\n    1989-2010",
                                                            era=="Late" ~ "\n        2011-2019"))) %>%
  ggplot() + geom_tile(aes(x=region, y=1,
                               fill=name), col="white") +
  scale_fill_manual(values=col.eras, guide="none") + labs(x="", y="") +
  facet_grid(rows=vars(factor(era_lab, levels=unique(era_lab))), switch="x") +
  labs(title="Era") +
  theme(panel.background=element_rect(fill="white"),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        strip.background = element_rect(fill="white"),
        strip.text.y.right = element_text(angle=0, margin=margin(l=-15)),
        plot.title = element_text(hjust=0.5)
        )
blank <- ggplot() + theme(panel.background = element_rect(fill="white")) # blank grob to put legend on
legend <- cowplot::ggdraw(blank) + cowplot::draw_plot(lgnd1, 0,0.4,1,0.2)


p <- cowplot::plot_grid(main_box, legend, ncol=2, rel_widths = c(1,0.19))
png(here('figures', 'spp-explore', 'main-fig-12-24.png'), height=700*2, width=900*2, res=72*2)
print(p)
dev.off()

# Save a version with SST only
p <- cowplot::plot_grid(sst_rw_box, legend, ncol=2, rel_widths = c(1,0.19))
png(here('figures', 'spp-explore', 'sst-box-rw.png'), height=700*2, width=900*2, res=72*2)
print(p)
dev.off()


# Plain boxplot
box %>% filter(era != "All") %>%
  mutate(regxera = factor(paste0(ocean_region_lab, sep=".", era),
                             levels=names(col.eras)[matrix(1:12, nrow=3, byrow=T)])
  ) %>%
  ggplot() +
  geom_hline(aes(yintercept=0), linetype="dashed", colour="gray50") +
  geom_boxplot(aes(x=sp.id, y=x, fill=regxera, col=ocean_region_lab), position=position_dodge2(preserve="single"), alpha=0.7) +
  scale_y_continuous(limits=c(-1, 1)) +
  facet_grid(rows=vars(Ocean.Region2), cols=vars(varnam)) +
  labs(y="Covariate Effect", x="", fill="") +
  scale_fill_manual(values=col.eras) +
  scale_colour_manual(values=unname(col.eras[9:12]), guide="none") +
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



#### 3) Random walk smoothed - all species #####


# a) sockeye

#stk.sub <- sock.info$Stock.ID[sock.info$yr_start < 1985 & sock.info$yr_end >= 2014]
#sock.sub <- sock %>% filter(Stock.ID %in% stk.sub)
#sock.info.sub <- sock.info %>% filter(Stock.ID %in% stk.sub)


# load sockeye RW
load("./output/models/dyn/sockeye/hbm_dyn_2c.RData", verbose=T)
summ <- rstan::summary(dyn.2c, pars = c("gamma", "kappa"), probs = NULL)$summary
rw.stk <- data.frame(Stock = sock$Stock,
                           Ocean.Region2 = sock$Ocean.Region2,
                           BY = sock$BY,
                           mu = summ[, "mean"],
                           se = summ[, "se_mean"],
                           var = str_extract(rownames(summ), "[a-z]+"),
                           varnam = case_when(grepl("^gamma", rownames(summ)) ~ "SST",
                                              grepl("^kappa", rownames(summ)) ~ "Competitors")
)



# Summarize at regional lvl
rw.reg <- dplyr::summarize(rw.stk,
                            reg_mean=mean(mu, na.rm=T),
                            lower_10=quantile(mu, 0.1, na.rm=T),
                            upper_90 = quantile(mu, 0.9, na.rm=T),
                            n_stk=n_distinct(Stock),
                            .by=c(Ocean.Region2, BY, varnam))
#rw.reg <- ddply(rw.reg, .(Ocean.Region2), dplyr::filter, n_stk >= 4)

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
                                lower_10=quantile(mu, 0.1, na.rm=T),
                                upper_90 = quantile(mu, 0.9, na.rm=T),
                                  n_stk=n_distinct(Stock),
                                  .by=c(Ocean.Region2, BY, varnam, even_odd))
#rw.reg.pink <- ddply(rw.reg.pink, .(Ocean.Region2), dplyr::filter, n_stk >= 4) # have to remove for now since 2025 data update

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
                                lower_10=quantile(mu, 0.1, na.rm=T),
                                upper_90 = quantile(mu, 0.9, na.rm=T),
                                n_stk=n_distinct(Stock),
                                .by=c(Ocean.Region2, BY, varnam))
#rw.reg.chum <- ddply(rw.reg.chum, .(Ocean.Region2), dplyr::filter, n_stk >= 4)

rw.reg.avg.chum <- moving_average_df(rw.reg.chum, "reg_mean", lag=2)

# bind all the moving average dfs
rw.mov <- bind_rows(rw.reg.avg.sock, rw.reg.avg.pink.even, rw.reg.avg.pink.odd, rw.reg.avg.chum, .id="sp.id") %>%
  mutate(spp=case_when(sp.id == 1 ~ "sockeye",
                       sp.id == 2 ~ "pink-even",
                       sp.id == 3 ~ "pink-odd",
                       sp.id == 4 ~ "chum")) %>%
  tidyr::replace_na(list(even_odd = "Both"))

# bind all the (raw) RW dfs
rw.reg.all <- bind_rows(rw.reg, rw.reg.chum, rw.reg.pink, .id="sp.id")
rw.reg.all <- rw.reg.all %>% mutate(sp.id=case_when(sp.id == 1 ~ "sockeye",
                                                  sp.id == 2 ~ "chum",
                                                  sp.id == 3 ~ "pink"))


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

#### 4) Map + covars re-design ####

## Fig: Map + covars (Hannah's version) --------------------

# colours
col.dk <- rev(chroma::qpal(7, luminance = 20)[c(1, 3, 5, 7)])
names(col.dk) <- unique(info_master$ocean_region_lab)


# Downlad map and convert to sp
na_map <- rnaturalearth::ne_states(country = c("United States of America", "Canada"), returnclass="sf")

axes <- list( xlims=c(-165.5, -121),
              ylims=c(47, 61),
              xbreaks=seq(-160,-120,10),
              xlabels=as.character(seq(-160,-120,10)),
              seq(45, 65, 5),
              ybreaks=seq(50, 60, 5),
              ylabels=as.character(seq(50,60,5)))

#make map data
map.info <- info_master %>% select(Stock, lon, lat, ocean_region_lab) %>%
  mutate(stock.no = 1:nrow(info_master)) %>%
  dplyr::summarize(n.stk = n_distinct(Stock),
                   first.stk = first(stock.no),
                   last.stk = last(stock.no),
                   .by= c("lat", "lon", "ocean_region_lab")) #%>%
#filter(!(first.stk %in% c(54:56))) # remove some overlapping points
map.info[map.info$first.stk %in% c(52:56, 41:45), c("first.stk", "last.stk")] <- NA
map.info <- map.info %>% mutate(num = ifelse(first.stk == last.stk, first.stk, paste(first.stk, last.stk, sep="-")))
misc.lab <- data.frame(x=c(-157.9, -152.5), y=c(59.5, 57.25), label=c("52-56", "41-45"))

map <- ggplot(map.info) +
  geom_sf(data=na_map, color="grey30", fill="white", linewidth=0.1, ) +
  ggspatial::geom_spatial_point(aes(x=lon, y=lat, col=ocean_region_lab, shape=ocean_region_lab, fill=ocean_region_lab),
                                crs=4326, size=1, stroke=1.75, alpha=0.7) +
  #geom_text(aes(x=lon, y=lat, label=num), vjust=1.4, col="gray20", size=3) +
  #ggrepel::geom_text_repel(aes(x=lon, y=lat, label=num), col="gray20", size=2.5, min.segment.length = 0.25, box.padding=0.1) +
  #geom_text(data=misc.lab, aes(x, y, label=label), col="gray20", size=2.5) +
  coord_sf(xlim=axes$xlims, ylim=axes$ylims) +
  scale_x_continuous(breaks=axes$xbreaks, labels=axes$xlabels) +
  scale_y_continuous(breaks=axes$ybreaks, labels=axes$ylabels) +
  scale_colour_manual(values=col.region, name="Ocean Region Grouping") +
  scale_fill_manual(values=col.dk, name="Ocean Region Grouping") +
  scale_shape_manual(values=c(22, 24, 21, 23), name="Ocean Region Grouping") +
  labs(x="Longitude (°E)", y="Latitude (°N)") +
  theme_sleek() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust=0.5),
        legend.position = "none"
  )

## Make dataframe of covariate info
covar.dat.st <- data_master %>% select(Stock, BY, Ocean.Region2, early_sst, np_pinks_sec) %>% tidyr::pivot_longer(cols=c(early_sst, np_pinks_sec), names_to = "covar") %>% mutate(covar_nam = ifelse(covar=="early_sst", "SST Index", "Competitor Index"))
covar.dat.st <- ocean_region_lab(covar.dat.st)
covar.dat.reg <- dplyr::summarize(.data=covar.dat.st, mean_covar = mean(value), .by = c("BY", "ocean_region_lab", "covar_nam"))

# plot comp
comp.covar <- covar.dat.st %>% filter(covar_nam == "Competitor Index") %>%
  ggplot() +
  geom_vline(xintercept=c(1976,1988), color = "grey50", linetype = 2, linewidth = 0.25) +
  geom_line(aes(x=BY, y=value, group=Stock), linewidth=0.5, alpha=0.4, col="darkred") +
  #geom_line(data=covar.dat.reg, aes(x=BY, y=mean_covar, col=ocean_region_lab), linewidth=0.5) +
  scale_colour_manual(values=col.region) +
  scale_y_continuous(n.breaks=4) +
  theme_sleek() +
  theme(legend.position="none",
        axis.title.x = element_text(size=8),
        axis.text = element_text(size=5),
        axis.line.x.top = element_blank(),
        plot.margin = margin(0,0,0,0)) +
  labs(x= "Brood Year", y="")

# plot SST
sst.covar <- covar.dat.st %>% filter(covar_nam == "SST Index") %>% ggplot() +
  geom_vline(xintercept=c(1976,1988), color = "grey50", linetype = 2, linewidth = 0.25) +
  geom_line(aes(x=BY, y=value, group=Stock, colour=ocean_region_lab), linewidth=0.5, alpha=0.4) +
  #geom_line(data=covar.dat.reg, aes(x=BY, y=mean_covar, col=ocean_region_lab), linewidth=0.5) +
  facet_grid(rows=vars(factor(ocean_region_lab, levels=rev(unique(ocean_region_lab)))), scales="free_y") +
  scale_colour_manual(values=col.region) +
  scale_y_continuous(n.breaks=4) +
  theme_sleek() + theme(legend.position="none") +
  labs(x= "Brood Year", y="")

main <- cowplot::ggdraw(map) + cowplot::draw_plot(comp.covar, x=.2,y=.2,width=.4,height=.4)

full <- cowplot::plot_grid(main, sst.covar, ncol=2, nrow=1, rel_widths = c(2, 1), rel_heights=c(1.3,1))


#### Avg productivity (lnRS) / time ####

# Also try productivity residuals?

sock %>% summarize(mean_RS = mean(RS_stnd), .by=c("Ocean.Region2", "BY")) %>%
  ggplot() +
  geom_point(aes(x=BY, y=Ocean.Region2, size=abs(mean_RS), col=mean_RS), alpha=0.5) +
  scale_color_distiller(palette="RdBu") +
  scale_size_continuous(guide="none") +
  theme_sleek() + theme() +
  labs(x="Brood Year", y="Region", col="Average stnd R/S")


#### Comparison table of ERA results ####

# load sockeye
load(here('output', 'models', 'dyn', 'sockeye', 'hbm_era_2c.RData'), verbose=T)
stbl <- rstan::summary(era.2c, par=c("mu_alpha", paste0("mu_gamma", 1:3), paste0("mu_kappa", 1:3)))$summary
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
load(here('output', 'models', 'stat', 'sockeye', 'stat_a.RData'), verbose=T)
statmeans[1] <- rstan::summary(stat_a, par=c("mu_gamma"))$summary[1,1]
load(here('output', 'models', 'stat', 'chum', 'stat_a.RData'), verbose=T)
statmeans[2] <- rstan::summary(stat_a, par=c("mu_gamma"))$summary[1,1]
load(here('output', 'models', 'stat', 'pink', 'stat_a_even.RData'), verbose=T)
statmeans[3] <- rstan::summary(stat_a_even, par=c("mu_gamma"))$summary[1,1]
load(here('output', 'models', 'stat', 'pink', 'stat_a_odd.RData'), verbose=T)
statmeans[4] <- rstan::summary(stat_a_odd, par=c("mu_gamma"))$summary[1,1]

filter(full_tbl, id=="mu_gamma3[1]")[c("mean", "species")]
