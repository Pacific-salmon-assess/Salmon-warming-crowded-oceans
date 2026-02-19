## Era model inference - remake boxplot fig from raw samples

## Pink

load(here('output', 'models', 'dyn', 'pink', 'hbm_era_2c.RData'), verbose=T) # eras

# Region-lvl dataframe
summ_avg <- rstan::extract(era.2c, pars = c(paste0("mu_gamma", 1:3), paste0("mu_kappa", 1:3)))
summ_avg <- lapply(summ_avg, as.data.frame)
summ_avg <- bind_rows(summ_avg, .id="par")
names(summ_avg) <- c("par", "West Coast", "Southeast Alaska", "Gulf of Alaska", "Bering Sea")

pink.df.avg <- summ_avg |> tidyr::pivot_longer(cols=-c("par"), names_to="ocean_region_lab") |>
  mutate(varnam = case_when(grepl("^mu_gamma", par) ~ "SST",
                            grepl("^mu_kappa", par) ~ "Competitors"),
        era = case_when(str_extract(par, "\\d")==1 ~ "Early",
                          str_extract(par, "\\d")==2 ~ "Middle",
                          str_extract(par, "\\d")==3 ~ "Late"),
        species="pink",
        pc = (exp(value)-1)*100)


## Chum

load(here('output', 'models', 'dyn', 'chum', 'hbm_era_2c.RData'), verbose=T) # eras

# Region-lvl dataframe
summ_avg <- rstan::extract(era.2c, pars = c(paste0("mu_gamma", 1:3), paste0("mu_kappa", 1:3)))
summ_avg <- lapply(summ_avg, as.data.frame)
summ_avg <- bind_rows(summ_avg, .id="par")
names(summ_avg) <- c("par", "West Coast", "Gulf of Alaska", "Bering Sea")

chum.df.avg <- summ_avg |> tidyr::pivot_longer(cols=-c("par"), names_to="ocean_region_lab") |>
  mutate(varnam = case_when(grepl("^mu_gamma", par) ~ "SST",
                            grepl("^mu_kappa", par) ~ "Competitors"),
         era = case_when(str_extract(par, "\\d")==1 ~ "Early",
                         str_extract(par, "\\d")==2 ~ "Middle",
                         str_extract(par, "\\d")==3 ~ "Late"),
         species="chum",
         pc = (exp(value)-1)*100)


## Sockeye

load(here('output', 'models', 'dyn', 'sockeye', 'hbm_era_2c.RData'), verbose=T) # eras

# Region-lvl dataframe
summ_avg <- rstan::extract(era.2c, pars = c(paste0("mu_gamma", 1:3), paste0("mu_kappa", 1:3)))
summ_avg <- lapply(summ_avg, as.data.frame)
summ_avg <- bind_rows(summ_avg, .id="par")
names(summ_avg) <- c("par", "West Coast", "Southeast Alaska", "Gulf of Alaska", "Bering Sea")

sock.df.avg <- summ_avg |> tidyr::pivot_longer(cols=-c("par"), names_to="ocean_region_lab") |>
  mutate(varnam = case_when(grepl("^mu_gamma", par) ~ "SST",
                            grepl("^mu_kappa", par) ~ "Competitors"),
         era = case_when(str_extract(par, "\\d")==1 ~ "Early",
                         str_extract(par, "\\d")==2 ~ "Middle",
                         str_extract(par, "\\d")==3 ~ "Late"),
         species="sockeye",
         pc = (exp(value)-1)*100)


# Bind all species together
df.avg <- rbind(sock.df.avg, chum.df.avg, pink.df.avg)


# Make simple boxplots for sst

# get colours
col.region <- rev(chroma::qpal(7, luminance = 40)[c(1, 3, 5, 7)])
col.dk <- rev(chroma::qpal(7, luminance = 20)[c(1, 3, 5, 7)])
names(col.region) <- names(col.dk) <- c("West Coast", "Southeast Alaska", "Gulf of Alaska", "Bering Sea")

# plot
df.avg |> filter(varnam=="SST") |> ggplot() +
  geom_boxplot(aes(x=factor(era, levels=c("Early", "Middle", "Late")), y=pc, fill=ocean_region_lab),
               outliers=FALSE) + # whiskers are not 95% exactly, consider
  facet_grid(rows=vars(ocean_region_lab), cols=vars(species), scales="free") +
  scale_fill_manual(values=col.region) +
  theme_sleek() +
  labs(x="Era", y="% Change in R/S (per degree C)")




## Load Random Walk --

# sockeye
sock.rw <- read.csv(here('figures/dyn/sockeye/hbm_inf/reg_coefficients_rw_sockeye.csv'))
obs_yr <-  sock |> summarize(n_obs=n_distinct(Stock), .by=c(Ocean.Region2, BY)) |> # prop stks observed in each year, region
   left_join(summarize(sock, n_tot = n_distinct(Stock), .by=c(Ocean.Region2))) |>
   mutate(pr_obs = n_obs/n_tot)
sock.rw <- sock.rw |> left_join(obs_yr[,c(1:3,5)], by=c("Ocean.Region2", "BY")) # bind this metric for use in plotting
sock.rw[is.na(sock.rw$n_obs),c("n_obs", "pr_obs")] <- 0 # replace nas with 0s for proportion metric
med_start <- sock.info |> group_by(ocean_region_lab) |> summarize(med_start = round(median(yr_start),0)) # alternatively, get median start date for each region
# convert values from effect size to % change
sock.rw <- sock.rw |> mutate(mu_pc = (exp(mu)-1)*100,
                             lower_10_pc = (exp(lower_10)-1)*100,
                             upper_90_pc = (exp(upper_90)-1)*100) |>
                      ocean_region_lab() |>
                      left_join(med_start, by="ocean_region_lab") |> # bind median start dates
                      filter(BY >= med_start)

# Chum
chum.rw <- read.csv(here('figures/dyn/chum/hbm_inf/reg_coefficients_rw_chum.csv'))
obs_yr <-  chum |> summarize(n_obs=n_distinct(Stock), .by=c(Ocean.Region2, BY)) |> # prop stks observed in each year, region
  left_join(summarize(chum, n_tot = n_distinct(Stock), .by=c(Ocean.Region2))) |>
  mutate(pr_obs = n_obs/n_tot)
chum.rw <- chum.rw |> left_join(obs_yr[,c(1:3,5)], by=c("Ocean.Region2", "BY")) # bind this metric for use in plotting
chum.rw[is.na(chum.rw$n_obs),c("n_obs", "pr_obs")] <- 0 # replace nas with 0s for proportion metric
med_start <- chum.info |> group_by(ocean_region_lab) |> summarize(med_start = round(median(yr_start),0)) # alternatively, get median start date for each region
# convert values from effect size to % change
chum.rw <- chum.rw |> mutate(mu_pc = (exp(mu)-1)*100,
                             lower_10_pc = (exp(lower_10)-1)*100,
                             upper_90_pc = (exp(upper_90)-1)*100) |>
                      ocean_region_lab() |>
                      left_join(med_start, by="ocean_region_lab") |> # bind median start dates
                      filter(BY >= med_start)

# Pink
pink.rw <- read.csv(here('figures/dyn/pink/hbm_inf/reg_coefficients_rw_pink.csv'))
obs_yr <-  pink |> summarize(n_obs=n_distinct(Stock), .by=c(Ocean.Region2, BY)) |> # prop stks observed in each year, region
  left_join(summarize(pink, n_tot = n_distinct(Stock), .by=c(Ocean.Region2))) |>
  mutate(pr_obs = n_obs/n_tot)
pink.rw <- pink.rw |> left_join(obs_yr[,c(1:3,5)], by=c("Ocean.Region2", "BY")) # bind this metric for use in plotting
pink.rw[is.na(pink.rw$n_obs),c("n_obs", "pr_obs")] <- 0 # replace nas with 0s for proportion metric
med_start <- pink.info |> group_by(ocean_region_lab) |> summarize(med_start = round(median(yr_start),0)) # alternatively, get median start date for each region
# convert values from effect size to % change
pink.rw <- pink.rw |> mutate(mu_pc = (exp(mu)-1)*100,
                          lower_10_pc = (exp(lower_10)-1)*100,
                          upper_90_pc = (exp(upper_90)-1)*100) |>
                      ocean_region_lab() |>
                      left_join(med_start, by="ocean_region_lab") |> # bind median start dates
                      filter(BY >= med_start)


# Bind all together & convert to % change
rw.mod.fits <- bind_rows(sock.rw, chum.rw, pink.rw, .id="species")




## Make Era + RW plots -- SST

# Sockeye
sock.box.rw <- sock.df.avg |>
      filter(varnam == "SST") |>
      mutate(BY = case_when(era=="Early" ~ 1975,
                                         era=="Middle" ~ 2000,
                                         era=="Late" ~ 2015)) |>
      ggplot() +
      geom_vline(xintercept = c(1989,2010), lty="dashed", col="gray65", alpha=0.5) +
      geom_ribbon(data=filter(sock.rw, varnam=="SST"),
                  aes(x=BY, ymin=lower_10_pc, ymax=upper_90_pc, fill=ocean_region_lab), alpha=0.2) +
      geom_line(data=filter(sock.rw, varnam=="SST"), aes(x=BY, y=mu_pc, col=ocean_region_lab)) +
      #geom_boxplot(aes(x=BY, y=pc, group=era, fill=ocean_region_lab), outliers=F) +
      geom_violin(aes(x=BY, y=pc, group=era, fill=ocean_region_lab, col=ocean_region_lab), alpha=0.75) +
      scale_y_continuous(breaks=c(-50,0,50)) +
      scale_x_continuous(limits=c(1960,2022), breaks=c(1970,1990,2010)) +
      facet_grid(rows=vars(ocean_region_lab)) +
      scale_fill_manual(values=col.region, guide=NULL) +
      scale_colour_manual(values=col.dk) +
      labs(x = "", y = "% Change in R/S (per \u00B0C)", title="Sockeye") +
      theme_sleek() +
      theme(legend.position="none",
            strip.text.y.right = element_blank(),
            plot.title = element_text(hjust=0.5))
png(here('figures/spp-explore', paste0("violin_rw_sockeye.png")))
print(sock.box.rw)
dev.off()

# Chum
chum.box.rw <- chum.df.avg |>
  filter(varnam == "SST") |>
  mutate(BY = case_when(era=="Early" ~ 1975,
                        era=="Middle" ~ 2000,
                        era=="Late" ~ 2015)) |>
  bind_rows(c(par=NA, ocean_region_lab="Southeast Alaska")) |> # add empty row for SEAK
  ggplot() +
  geom_vline(xintercept = c(1989,2010), lty="dashed", col="gray65", alpha=0.5) +
  geom_ribbon(data=filter(chum.rw, varnam=="SST"),
              aes(x=BY, ymin=lower_10_pc, ymax=upper_90_pc, fill=ocean_region_lab),
              alpha=0.2) +
  geom_line(data=filter(chum.rw, varnam=="SST"), aes(x=BY, y=mu_pc, col=ocean_region_lab)) +
  #geom_boxplot(aes(x=BY, y=pc, group=era, fill=ocean_region_lab), outliers=F) +
  geom_violin(aes(x=BY, y=pc, group=era, fill=ocean_region_lab, col=ocean_region_lab), alpha=0.75) +
  scale_y_continuous(breaks=c(-25,0,25)) +
  scale_x_continuous(limits=c(1960,2022), breaks=c(1970,1990,2010)) +
  facet_grid(rows=vars(ocean_region_lab)) +
  scale_fill_manual(values=col.region, guide=NULL) +
  scale_colour_manual(values=col.dk) +
  labs(x = "Brood Year", y = "", title="Chum") +
  theme_sleek() +
  theme(legend.position="none",
        strip.text.y.right = element_blank(),
        plot.title = element_text(hjust=0.5))

png(here('figures/spp-explore', paste0("violin_rw_chum.png")))
print(chum.box.rw)
dev.off()


# Pink
pink.box.rw <- pink.df.avg |>
  filter(varnam == "SST") |>
  mutate(BY = case_when(era=="Early" ~ 1975,
                        era=="Middle" ~ 2000,
                        era=="Late" ~ 2015)) |>
  ggplot() +
  geom_vline(xintercept = c(1989,2010), lty="dashed", col="gray65", alpha=0.5) +
  geom_ribbon(data=filter(pink.rw, varnam=="SST"),
              aes(x=BY, ymin=lower_10_pc, ymax=upper_90_pc, fill=ocean_region_lab),
              alpha=0.2) +
  geom_line(data=filter(pink.rw, varnam=="SST"), aes(x=BY, y=mu_pc, col=ocean_region_lab)) +
  #geom_boxplot(aes(x=BY, y=pc, group=era, fill=ocean_region_lab), outliers=F) +
  geom_violin(aes(x=BY, y=pc, group=era, fill=ocean_region_lab, col=ocean_region_lab), alpha=0.75) +
  scale_y_continuous(n.breaks=4) +
  scale_x_continuous(limits=c(1960,2022), breaks=c(1970,1990,2010)) +
  facet_grid(rows=vars(ocean_region_lab), scales="free_y") +
  scale_fill_manual(values=col.region, guide=NULL) +
  scale_colour_manual(values=col.dk) +
  labs(x = "", y = "", col="Region", title="Pink") +
  theme_sleek() +
  theme(legend.position="none",
        plot.title = element_text(hjust=0.5))

png(here('figures/spp-explore', paste0("violin_rw_pink.png")))
print(pink.box.rw)
dev.off()

full.violin <- cowplot::plot_grid(sock.box.rw, chum.box.rw + theme(axis.title.y.left = element_blank()), pink.box.rw + theme(axis.title.y.left = element_blank()), ncol=3, rel_widths=c(1,1,1.1))
png(here('figures/spp-explore/violin_full.png'), width=900*3, height=600*3, res=72*5)
print(full.violin)
dev.off()



## Make Era + RW plots -- Competitors

# Sockeye
sock.box.rw <- sock.df.avg |>
  filter(varnam == "Competitors") |>
  mutate(BY = case_when(era=="Early" ~ 1975,
                        era=="Middle" ~ 2000,
                        era=="Late" ~ 2015)) |>
  ggplot() +
  geom_vline(xintercept = c(1989,2010), lty="dashed", col="gray65", alpha=0.5) +
  geom_ribbon(data=filter(sock.rw, varnam=="Competitors"),
              aes(x=BY, ymin=lower_10_pc, ymax=upper_90_pc, fill=ocean_region_lab), alpha=0.2) +
  geom_line(data=filter(sock.rw, varnam=="Competitors"), aes(x=BY, y=mu_pc, col=ocean_region_lab)) +
  #geom_boxplot(aes(x=BY, y=pc, group=era, fill=ocean_region_lab), outliers=F) +
  geom_violin(aes(x=BY, y=pc, group=era, fill=ocean_region_lab, col=ocean_region_lab), alpha=0.75) +
  scale_y_continuous(breaks=c(-50,0,50)) +
  scale_x_continuous(limits=c(1960,2022), breaks=c(1970,1990,2010)) +
  facet_grid(rows=vars(ocean_region_lab)) +
  scale_fill_manual(values=col.region, guide=NULL) +
  scale_colour_manual(values=col.dk) +
  labs(x = "", y = "% Change in R/S (per 155 million pinks)", title="Sockeye") +
  theme_sleek() +
  theme(legend.position="none",
        strip.text.y.right = element_blank(),
        plot.title = element_text(hjust=0.5))
png(here('figures/spp-explore', paste0("violin_rw_sockeye_comp.png")))
print(sock.box.rw)
dev.off()

# Chum
chum.box.rw <- chum.df.avg |>
  filter(varnam == "Competitors") |>
  mutate(BY = case_when(era=="Early" ~ 1975,
                        era=="Middle" ~ 2000,
                        era=="Late" ~ 2015)) |>
  bind_rows(c(par=NA, ocean_region_lab="Southeast Alaska")) |> # add empty row for SEAK
  ggplot() +
  geom_vline(xintercept = c(1989,2010), lty="dashed", col="gray65", alpha=0.5) +
  geom_ribbon(data=filter(chum.rw, varnam=="Competitors"),
              aes(x=BY, ymin=lower_10_pc, ymax=upper_90_pc, fill=ocean_region_lab),
              alpha=0.2) +
  geom_line(data=filter(chum.rw, varnam=="Competitors"), aes(x=BY, y=mu_pc, col=ocean_region_lab)) +
  #geom_boxplot(aes(x=BY, y=pc, group=era, fill=ocean_region_lab), outliers=F) +
  geom_violin(aes(x=BY, y=pc, group=era, fill=ocean_region_lab, col=ocean_region_lab), alpha=0.75) +
  scale_y_continuous(n.breaks=4) +
  scale_x_continuous(limits=c(1960,2022), breaks=c(1970,1990,2010)) +
  facet_grid(rows=vars(ocean_region_lab)) +
  scale_fill_manual(values=col.region, guide=NULL) +
  scale_colour_manual(values=col.dk) +
  labs(x = "Brood Year", y = "", title="Chum") +
  theme_sleek() +
  theme(legend.position="none",
        strip.text.y.right = element_blank(),
        plot.title = element_text(hjust=0.5))

png(here('figures/spp-explore', paste0("violin_rw_chum_comp.png")))
print(chum.box.rw)
dev.off()


# Pink
pink.box.rw <- pink.df.avg |>
  filter(varnam == "Competitors") |>
  mutate(BY = case_when(era=="Early" ~ 1975,
                        era=="Middle" ~ 2000,
                        era=="Late" ~ 2015)) |>
  ggplot() +
  geom_vline(xintercept = c(1989,2010), lty="dashed", col="gray65", alpha=0.5) +
  geom_ribbon(data=filter(pink.rw, varnam=="Competitors"),
              aes(x=BY, ymin=lower_10_pc, ymax=upper_90_pc, fill=ocean_region_lab),
              alpha=0.2) +
  geom_line(data=filter(pink.rw, varnam=="Competitors"), aes(x=BY, y=mu_pc, col=ocean_region_lab)) +
  #geom_boxplot(aes(x=BY, y=pc, group=era, fill=ocean_region_lab), outliers=F) +
  geom_violin(aes(x=BY, y=pc, group=era, fill=ocean_region_lab, col=ocean_region_lab), alpha=0.75) +
  scale_y_continuous(n.breaks=4) +
  scale_x_continuous(limits=c(1960,2022), breaks=c(1970,1990,2010)) +
  facet_grid(rows=vars(ocean_region_lab), scales="free_y") +
  scale_fill_manual(values=col.region, guide=NULL) +
  scale_colour_manual(values=col.dk) +
  labs(x = "", y = "", col="Region", title="Pink") +
  theme_sleek() +
  theme(legend.position="none",
        plot.title = element_text(hjust=0.5))

png(here('figures/spp-explore', paste0("violin_rw_pink_comp.png")))
print(pink.box.rw)
dev.off()

full.violin.comp <- cowplot::plot_grid(sock.box.rw, chum.box.rw + theme(axis.title.y.left = element_blank()), pink.box.rw + theme(axis.title.y.left = element_blank()), ncol=3, rel_widths=c(1,1,1.1))
png(here('figures/spp-explore/violin_full_comp.png'), width=900*3, height=600*3, res=72*5)
print(full.violin.comp)
dev.off()
