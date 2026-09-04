library(ggsidekick) #for theme_sleek() - doesn't work with some vs of R, hence the comment


## All spp. interaction density plot -------------------------------------------
s.df <- m.df <- NULL
for(sp in c("sockeye", "pink", "chum")){
  print(sp)
  # Species
  if(sp=="pink") {
    data_tmp <- pink
    info_tmp <- pink.info} else if (sp=="chum") {
      data_tmp <- chum
      info_tmp <- chum.info } else if(sp=="sockeye"){
        data_tmp <- sock
        info_tmp <- sock.info }
  load(here("output", "models", "stat", sp, "stat_inter.RData"), verbose=T)
  lst <- hb07_density_df(stat_inter, data=data_tmp, info = info_tmp,
                         ocean.regions = ifelse(sp=="chum", 3, 4))
  lst$stock$species <- sp
  lst$region$species <- sp
  s.df <- rbind(s.df, lst$stock)
  m.df <- rbind(m.df, lst$region)
  m.df$region <- factor(m.df$region, levels = c("West Coast", "Gulf of Alaska", "Southeast Alaska", "Bering Sea"))
}

## Covariate and species labels
m.df <- filter(m.df, var != "SST + Comp")
s.df <- filter(s.df, var != "SST + Comp")
vars <- data.frame(var = unique(m.df$var))
vars$lab <- paste0("(", letters[1:nrow(vars)], ") ", vars$var)
vars$var <- factor(vars$var, levels = c("SST", "Comp", "SST x Comp"))
m.df$species <- factor(str_to_sentence(m.df$species), levels=c("Sockeye","Chum","Pink"))
s.df$species <- factor(str_to_sentence(s.df$species), levels=c("Sockeye","Chum","Pink"))

m.df.plot <- m.df
s.df.plot <- s.df
# sneaky trim data to display with different x-axis ranges without huge tails

# m.df.plot <- filter(m.df[m.df$species == "Pink",], x<=100, x>=-50) # trim pinks to -50,100
# m.df.plot <- rbind(m.df.plot, filter(m.df[m.df$species == "Chum",], x<=25, x>=-25)) # trim chum to -25,25
# m.df.plot <- rbind(m.df.plot, filter(m.df[m.df$species == "Sockeye",], x<=50, x>=-50)) # trim sockeye to -50,50
#
# s.df.plot <- filter(s.df[s.df$species == "Pink",], x<=100, x>=-50) # trim pinks to 100
# s.df.plot <- rbind(s.df.plot, filter(s.df[s.df$species == "Chum",], x<=25, x>=-25)) # trim chum to -25,25
# s.df.plot <- rbind(s.df.plot, filter(s.df[s.df$species == "Sockeye",], x<=50, x>=-50)) # trim sockeye to -50,50


g<-ggplot(m.df.plot) +
  geom_vline(xintercept = 0, color = "grey50", linetype = 2, linewidth = 0.25) +
  geom_path(data=s.df.plot, aes(x=x, y=y, group=stock, color = region), alpha=0.3, na.rm=T) +
  geom_path(aes(x = x, y = y, color = region), linewidth = 1, alpha=1,
            na.rm = TRUE) +
  scale_colour_manual(values=col.region) +
  labs(x = "Percent change in R/S",
       y = "Posterior density",
       color = "") +
  scale_y_continuous(breaks=NULL) +
  facet_grid(rows=vars(var), cols=vars(species)) +
  xlim(-50,50) +
  theme_sleek(base_size = 9) +
  theme(legend.justification = c(0, 0),
        legend.position = c(0.83, 0.76),
        legend.key.size = unit(10, "pt"),
        legend.background = element_blank(),
        legend.text = element_text(size = 6),
        panel.spacing.y = unit(-0.5, "pt"),
        strip.background = element_blank(),
        strip.text.x.top = element_text(size=9),
        legend.key.spacing = unit(-0.5, "pt"))


png(here('figures', 'manuscript/main-text', "posterior-percent-change.png"), width = 800*2, height = 500*2, res = 72*4)
print(g)
dev.off()



## High/Low SST, Comp vs Recruits-per-Spawner -----------------------------

# sockeye ----
load(here("output", "models", "stat", "sockeye", "stat_inter.RData"), verbose=T)
sockeye_post <- as.data.frame(stat_inter)

sst_index <- c(-2,2)
comp_index <- seq(0,3,length.out=100)

# calculate what raw all species abundance is at 0-3 SDUs
raw.comp$all_spp_anomaly <- scale(raw.comp$all_spp_numbers_np)
comp_index_0 <- which.min(abs(raw.comp$all_spp_anomaly - 0))
comp_index_1 <- which.min(abs(raw.comp$all_spp_anomaly - 1))
comp_index_2 <- which.min(abs(raw.comp$all_spp_anomaly - 2))
comp_index_3 <- which.min(abs(raw.comp$all_spp_anomaly - 3))
comp_label <- round(c(raw.comp$all_spp_numbers_np[comp_index_0], raw.comp$all_spp_numbers_np[comp_index_1], raw.comp$all_spp_numbers_np[comp_index_2],
                raw.comp$all_spp_numbers_np[comp_index_3]))


pred.all.df <- data.frame(
  prod = numeric(0),
  comp = numeric(0),
  sst = numeric(0),
  region = character(0))


  for (i in 1:1000){ # west coast
    samp <- sample(4000,1)
    alpha <- sockeye_post$`mu_alpha[2]`[samp]
    sst <- sockeye_post$`mu_gamma[1]`[samp] # SST effect
    comp <- sockeye_post$`mu_kappa[1]`[samp] # comp effect
    sstXcomp <- sockeye_post$`mu_chi[1]`[samp] # interaction

    pred.cool <- alpha +(sst*sst_index[1]) + (comp*comp_index) + (sstXcomp*(comp_index*sst_index[1])) # cool years
    pred.warm <- alpha +(sst*sst_index[2]) + (comp*comp_index) + (sstXcomp*(comp_index*sst_index[2])) # cool years

    pred.df <- cbind(c(pred.cool,pred.warm),
                   c(comp_index,comp_index),
                   c(rep("cool",100),rep("warm",100)),
                   c(rep("West Coast",200)))

    pred.df <- as.data.frame(pred.df)
    pred.all.df <- rbind(pred.all.df,pred.df)
  }



for (i in 1:1000){ # SEAK
  samp <- sample(4000,1)
  alpha <- sockeye_post$`mu_alpha[2]`[samp]
  sst <- sockeye_post$`mu_gamma[2]`[samp] # SST effect
  comp <- sockeye_post$`mu_kappa[2]`[samp] # comp effect
  sstXcomp <- sockeye_post$`mu_chi[2]`[samp] # interaction

  pred.cool <- alpha +(sst*sst_index[1]) + (comp*comp_index) + (sstXcomp*(comp_index*sst_index[1])) # cool years
  pred.warm <- alpha +(sst*sst_index[2]) + (comp*comp_index) + (sstXcomp*(comp_index*sst_index[2])) # cool years

  pred.df <- cbind(c(pred.cool,pred.warm),
                   c(comp_index,comp_index),
                   c(rep("cool",100),rep("warm",100)),
                   c(rep("SEAK",200)))

  pred.df <- as.data.frame(pred.df)
  pred.all.df <- rbind(pred.all.df,pred.df)
}

for (i in 1:1000){ # GoA
  samp <- sample(4000,1)
  alpha <- sockeye_post$`mu_alpha[2]`[samp]
  sst <- sockeye_post$`mu_gamma[3]`[samp] # SST effect
  comp <- sockeye_post$`mu_kappa[3]`[samp] # comp effect
  sstXcomp <- sockeye_post$`mu_chi[3]`[samp] # interaction

  pred.cool <- alpha +(sst*sst_index[1]) + (comp*comp_index) + (sstXcomp*(comp_index*sst_index[1])) # cool years
  pred.warm <- alpha +(sst*sst_index[2]) + (comp*comp_index) + (sstXcomp*(comp_index*sst_index[2])) # cool years

  pred.df <- cbind(c(pred.cool,pred.warm),
                   c(comp_index,comp_index),
                   c(rep("cool",100),rep("warm",100)),
                   c(rep("GoA",200)))

  pred.df <- as.data.frame(pred.df)
  pred.all.df <- rbind(pred.all.df,pred.df)
}

for (i in 1:1000){ # Bering Sea
  samp <- sample(4000,1)
  alpha <- sockeye_post$`mu_alpha[2]`[samp]
  sst <- sockeye_post$`mu_gamma[4]`[samp] # SST effect
  comp <- sockeye_post$`mu_kappa[4]`[samp] # comp effect
  sstXcomp <- sockeye_post$`mu_chi[4]`[samp] # interaction

  pred.cool <- alpha +(sst*sst_index[1]) + (comp*comp_index) + (sstXcomp*(comp_index*sst_index[1])) # cool years
  pred.warm <- alpha +(sst*sst_index[2]) + (comp*comp_index) + (sstXcomp*(comp_index*sst_index[2])) # cool years

  pred.df <- cbind(c(pred.cool,pred.warm),
                   c(comp_index,comp_index),
                   c(rep("cool",100),rep("warm",100)),
                   c(rep("Bering Sea",200)))

  pred.df <- as.data.frame(pred.df)
  pred.all.df <- rbind(pred.all.df,pred.df)
}

colnames(pred.all.df) <- c("prod","comp","sst","region")

summ.pred.all.df <- pred.all.df |>
  mutate(prod = as.numeric(prod),
         comp = as.numeric(comp)) |>
  group_by(comp, sst, region) |>
  summarise(med.pred = median(prod),
            lwr.pred = quantile(prod, probs = 0.2, na.rm = TRUE),
            upr.pred = quantile(prod, probs = 0.8, na.rm = TRUE))

custom_colors <- c("cool" = "blue", "warm" = "red")

g <- ggplot(data = summ.pred.all.df, aes(x = comp, y = exp(med.pred), col=sst)) +
  geom_ribbon(aes(ymin=exp(lwr.pred), ymax=exp(upr.pred), col=sst, fill=sst), alpha=0.2) +
  geom_line(aes(col=sst),lwd=1.05) +
  theme_sleek() +
  xlab("North Pacific competitor abundance (millions)") +
  ylab("Recruits-per-spawner") +
  scale_color_manual(values = custom_colors)+
  scale_fill_manual(values = custom_colors) +
  facet_wrap(~region, ncol=2) +
  labs(fill = "SST") +
  scale_x_continuous(breaks=c(0,1,2,3),labels=comp_label)+
  guides(color="none") +
  theme(legend.justification = c(0, 0),
        legend.position = c(0.05, 0.3),
        legend.key.size = unit(11, "pt"),
        legend.background = element_blank(),
        legend.title = element_text( size = 8))

png(here('figures', 'spp-explore', "sockeye-sst-inter.png"), units = "in", width = 6, height = 5, res = 500)
print(g)
dev.off()

# quick stats: % reduction in R/S for each region in warm and cool years
# west coast
wc_hatch_effect <- summ.pred.all.df |>
  filter(region == "West Coast",
         sst == "warm",
         comp %in% c(2,3)) |>
  mutate(rawRS = exp(med.pred))
(1-(wc_hatch_effect[2,7]/wc_hatch_effect[1,7]))*100

wc_hatch_effect <- summ.pred.all.df |>
  filter(region == "West Coast",
         sst == "cool",
         comp %in% c(2,3)) |>
  mutate(rawRS = exp(med.pred))
(1-(wc_hatch_effect[2,7]/wc_hatch_effect[1,7]))*100

# SEAK
SEAK_hatch_effect <- summ.pred.all.df |>
  filter(region == "SEAK",
         sst == "warm",
         comp %in% c(2,3))|>
  mutate(rawRS = exp(med.pred))
(1-(SEAK_hatch_effect[2,7]/SEAK_hatch_effect[1,7]))*100

SEAK_hatch_effect <- summ.pred.all.df |>
  filter(region == "SEAK",
         sst == "cool",
         comp %in% c(2,3))|>
  mutate(rawRS = exp(med.pred))
(1-(SEAK_hatch_effect[2,7]/SEAK_hatch_effect[1,7]))*100

# GoA
GoA_hatch_effect <- summ.pred.all.df |>
  filter(region == "GoA",
         sst == "warm",
         comp %in% c(2,3)) |>
  mutate(rawRS = exp(med.pred))
(1-(GoA_hatch_effect[2,7]/GoA_hatch_effect[1,7]))*100

GoA_hatch_effect <- summ.pred.all.df |>
  filter(region == "GoA",
         sst == "cool",
         comp %in% c(2,3)) |>
  mutate(rawRS = exp(med.pred))
(1-(GoA_hatch_effect[2,7]/GoA_hatch_effect[1,7]))*100

# Bering Sea
BS_hatch_effect <- summ.pred.all.df |>
  filter(region == "Bering Sea",
         sst == "warm",
         comp %in% c(2,3)) |>
  mutate(rawRS = exp(med.pred))
(1-(BS_hatch_effect[2,7]/BS_hatch_effect[1,7]))*100

BS_hatch_effect <- summ.pred.all.df |>
  filter(region == "Bering Sea",
         sst == "cool",
         comp %in% c(2,3))  |>
  mutate(rawRS = exp(med.pred))
(1-(BS_hatch_effect[2,7]/BS_hatch_effect[1,7]))*100

# chum ----
load(here("output", "models", "stat", "chum", "stat_inter.RData"), verbose=T)
chum_post <- as.data.frame(stat_inter)

sst_index <- c(-2,2)
comp_index <- seq(0,3,length.out=100)

pred.all.df <- data.frame(
  prod = numeric(0),
  comp = numeric(0),
  sst = numeric(0),
  region = character(0))


for (i in 1:1000){ # west coast
  samp <- sample(4000,1)
  alpha <- chum_post$`mu_alpha[1]`[samp]
  sst <- chum_post$`mu_gamma[1]`[samp] # SST effect
  comp <- chum_post$`mu_kappa[1]`[samp] # comp effect
  sstXcomp <- chum_post$`mu_chi[1]`[samp] # interaction

  pred.cool <- alpha +(sst*sst_index[1]) + (comp*comp_index) + (sstXcomp*(comp_index*sst_index[1])) # cool years
  pred.warm <- alpha +(sst*sst_index[2]) + (comp*comp_index) + (sstXcomp*(comp_index*sst_index[2])) # cool years

  pred.df <- cbind(c(pred.cool,pred.warm),
                   c(comp_index,comp_index),
                   c(rep("cool",100),rep("warm",100)),
                   c(rep("West Coast",200)))

  pred.df <- as.data.frame(pred.df)
  pred.all.df <- rbind(pred.all.df,pred.df)
}



for (i in 1:1000){ # GoA
  samp <- sample(4000,1)
  alpha <- chum_post$`mu_alpha[1]`[samp]
  sst <- chum_post$`mu_gamma[2]`[samp] # SST effect
  comp <- chum_post$`mu_kappa[2]`[samp] # comp effect
  sstXcomp <- chum_post$`mu_chi[2]`[samp] # interaction

  pred.cool <- alpha +(sst*sst_index[1]) + (comp*comp_index) + (sstXcomp*(comp_index*sst_index[1])) # cool years
  pred.warm <- alpha +(sst*sst_index[2]) + (comp*comp_index) + (sstXcomp*(comp_index*sst_index[2])) # cool years

  pred.df <- cbind(c(pred.cool,pred.warm),
                   c(comp_index,comp_index),
                   c(rep("cool",100),rep("warm",100)),
                   c(rep("GoA",200)))

  pred.df <- as.data.frame(pred.df)
  pred.all.df <- rbind(pred.all.df,pred.df)
}

for (i in 1:1000){ # Bering Sea
  samp <- sample(4000,1)
  alpha <- chum_post$`mu_alpha[1]`[samp]
  sst <- chum_post$`mu_gamma[3]`[samp] # SST effect
  comp <- chum_post$`mu_kappa[3]`[samp] # comp effect
  sstXcomp <- chum_post$`mu_chi[3]`[samp] # interaction

  pred.cool <- alpha +(sst*sst_index[1]) + (comp*comp_index) + (sstXcomp*(comp_index*sst_index[1])) # cool years
  pred.warm <- alpha +(sst*sst_index[2]) + (comp*comp_index) + (sstXcomp*(comp_index*sst_index[2])) # cool years

  pred.df <- cbind(c(pred.cool,pred.warm),
                   c(comp_index,comp_index),
                   c(rep("cool",100),rep("warm",100)),
                   c(rep("Bering Sea",200)))

  pred.df <- as.data.frame(pred.df)
  pred.all.df <- rbind(pred.all.df,pred.df)
}


colnames(pred.all.df) <- c("prod","comp","sst","region")

summ.pred.all.df <- pred.all.df |>
  mutate(prod = as.numeric(prod),
         comp = as.numeric(comp)) |>
  group_by(comp, sst, region) |>
  summarise(med.pred = median(prod),
            lwr.pred = quantile(prod, probs = 0.2, na.rm = TRUE),
            upr.pred = quantile(prod, probs = 0.8, na.rm = TRUE))

custom_colors <- c("cool" = "blue", "warm" = "red")

g <- ggplot(data = summ.pred.all.df, aes(x = comp, y = exp(med.pred), col=sst)) +
  geom_ribbon(aes(ymin=exp(lwr.pred), ymax=exp(upr.pred), col=sst, fill=sst), alpha=0.2) +
  geom_line(aes(col=sst),lwd=1.05) +
  theme_sleek() +
  xlab("North Pacific competitor abundance (millions)") +
  ylab("Recruits-per-spawner") +
  scale_color_manual(values = custom_colors)+
  scale_fill_manual(values = custom_colors) +
  facet_wrap(~region,ncol=2) +
  labs(fill = "SST") +
  scale_x_continuous(breaks=c(0,1,2,3),labels=comp_label)+
  guides(color="none")+
  theme(legend.justification = c(0, 0),
        legend.position = c(0.05, 0.35),
        legend.key.size = unit(11, "pt"),
        legend.background = element_blank(),
        legend.title = element_text( size = 8))

png(here('figures', 'spp-explore', "chum-sst-inter.png"), units = "in", width = 6, height = 5, res = 500)
print(g)
dev.off()


# pink ----
load(here("output", "models", "stat", "pink", "stat_inter.RData"), verbose=T)
pink_post <- as.data.frame(stat_inter)

sst_index <- c(-2,2)
comp_index <- seq(0,3,length.out=100)

pred.all.df <- data.frame(
  prod = numeric(0),
  comp = numeric(0),
  sst = numeric(0),
  region = character(0))


for (i in 1:1000){ # west coast
  samp <- sample(4000,1)
  alpha <- pink_post$`mu_alpha[1]`[samp]
  sst <- pink_post$`mu_gamma[1]`[samp] # SST effect
  comp <- pink_post$`mu_kappa[1]`[samp] # comp effect
  sstXcomp <- pink_post$`mu_chi[1]`[samp] # interaction

  pred.cool <- alpha +(sst*sst_index[1]) + (comp*comp_index) + (sstXcomp*(comp_index*sst_index[1])) # cool years
  pred.warm <- alpha +(sst*sst_index[2]) + (comp*comp_index) + (sstXcomp*(comp_index*sst_index[2])) # cool years

  pred.df <- cbind(c(pred.cool,pred.warm),
                   c(comp_index,comp_index),
                   c(rep("cool",100),rep("warm",100)),
                   c(rep("West Coast",200)))

  pred.df <- as.data.frame(pred.df)
  pred.all.df <- rbind(pred.all.df,pred.df)
}



for (i in 1:1000){ # SEAK
  samp <- sample(4000,1)
  alpha <- pink_post$`mu_alpha[1]`[samp]
  sst <- pink_post$`mu_gamma[2]`[samp] # SST effect
  comp <- pink_post$`mu_kappa[2]`[samp] # comp effect
  sstXcomp <- pink_post$`mu_chi[2]`[samp] # interaction

  pred.cool <- alpha +(sst*sst_index[1]) + (comp*comp_index) + (sstXcomp*(comp_index*sst_index[1])) # cool years
  pred.warm <- alpha +(sst*sst_index[2]) + (comp*comp_index) + (sstXcomp*(comp_index*sst_index[2])) # cool years

  pred.df <- cbind(c(pred.cool,pred.warm),
                   c(comp_index,comp_index),
                   c(rep("cool",100),rep("warm",100)),
                   c(rep("SEAK",200)))

  pred.df <- as.data.frame(pred.df)
  pred.all.df <- rbind(pred.all.df,pred.df)
}

for (i in 1:1000){ # GoA
  samp <- sample(4000,1)
  alpha <- pink_post$`mu_alpha[1]`[samp]
  sst <- pink_post$`mu_gamma[3]`[samp] # SST effect
  comp <- pink_post$`mu_kappa[3]`[samp] # comp effect
  sstXcomp <- pink_post$`mu_chi[3]`[samp] # interaction

  pred.cool <- alpha +(sst*sst_index[1]) + (comp*comp_index) + (sstXcomp*(comp_index*sst_index[1])) # cool years
  pred.warm <- alpha +(sst*sst_index[2]) + (comp*comp_index) + (sstXcomp*(comp_index*sst_index[2])) # cool years

  pred.df <- cbind(c(pred.cool,pred.warm),
                   c(comp_index,comp_index),
                   c(rep("cool",100),rep("warm",100)),
                   c(rep("GoA",200)))

  pred.df <- as.data.frame(pred.df)
  pred.all.df <- rbind(pred.all.df,pred.df)
}

for (i in 1:1000){ # Bering Sea
  samp <- sample(4000,1)
  alpha <- pink_post$`mu_alpha[1]`[samp]
  sst <- pink_post$`mu_gamma[4]`[samp] # SST effect
  comp <- pink_post$`mu_kappa[4]`[samp] # comp effect
  sstXcomp <- pink_post$`mu_chi[4]`[samp] # interaction

  pred.cool <- alpha +(sst*sst_index[1]) + (comp*comp_index) + (sstXcomp*(comp_index*sst_index[1])) # cool years
  pred.warm <- alpha +(sst*sst_index[2]) + (comp*comp_index) + (sstXcomp*(comp_index*sst_index[2])) # cool years

  pred.df <- cbind(c(pred.cool,pred.warm),
                   c(comp_index,comp_index),
                   c(rep("cool",100),rep("warm",100)),
                   c(rep("Bering Sea",200)))

  pred.df <- as.data.frame(pred.df)
  pred.all.df <- rbind(pred.all.df,pred.df)
}

colnames(pred.all.df) <- c("prod","comp","sst","region")

summ.pred.all.df <- pred.all.df |>
  mutate(prod = as.numeric(prod),
         comp = as.numeric(comp)) |>
  group_by(comp, sst, region) |>
  summarise(med.pred = median(prod),
            lwr.pred = quantile(prod, probs = 0.2, na.rm = TRUE),
            upr.pred = quantile(prod, probs = 0.8, na.rm = TRUE))

custom_colors <- c("cool" = "blue", "warm" = "red")

g <- ggplot(data = summ.pred.all.df, aes(x = comp, y = exp(med.pred), col=sst)) +
  geom_ribbon(aes(ymin=exp(lwr.pred), ymax=exp(upr.pred), col=sst, fill=sst), alpha=0.2) +
  geom_line(aes(col=sst),lwd=1.05) +
  theme_sleek() +
  xlab("North Pacific competitor abundance (millions)") +
  ylab("Recruits-per-spawner") +
  scale_color_manual(values = custom_colors)+
  scale_fill_manual(values = custom_colors) +
  facet_wrap(~region, ncol=2) +
  labs(fill = "SST") +
  scale_x_continuous(breaks=c(0,1,2,3),labels=comp_label)+
  guides(color="none")+
  theme(legend.justification = c(0, 0),
        legend.position = c(0.05, 0.35),
        legend.key.size = unit(11, "pt"),
        legend.background = element_blank(),
        legend.title = element_text( size = 8))

png(here('figures', 'spp-explore', "pink-sst-inter.png"), units = "in", width = 6, height = 5, res = 500)
print(g)
dev.off()

# all spp interaction plot ----
# sockeye ----
load(here("output", "models", "stat", "sockeye", "stat_inter.RData"), verbose=T)
sockeye_post <- as.data.frame(stat_inter)

sst_index <- c(-2,2)
comp_index <- seq(0,3,length.out=100)

pred.all.df <- data.frame(
  prod = numeric(0),
  comp = numeric(0),
  sst = numeric(0),
  region = character(0))


for (i in 1:1000){ # west coast
  samp <- sample(4000,1)
  alpha <- sockeye_post$`mu_alpha[1]`[samp]
  sst <- sockeye_post$`mu_gamma[1]`[samp] # SST effect
  comp <- sockeye_post$`mu_kappa[1]`[samp] # comp effect
  sstXcomp <- sockeye_post$`mu_chi[1]`[samp] # interaction

  pred.cool <- alpha +(sst*sst_index[1]) + (comp*comp_index) + (sstXcomp*(comp_index*sst_index[1])) # cool years
  pred.warm <- alpha +(sst*sst_index[2]) + (comp*comp_index) + (sstXcomp*(comp_index*sst_index[2])) # cool years

  pred.df <- cbind(c(pred.cool,pred.warm),
                   c(comp_index,comp_index),
                   c(rep("cool",100),rep("warm",100)),
                   c(rep("West Coast",200)))

  pred.df <- as.data.frame(pred.df)
  pred.all.df <- rbind(pred.all.df,pred.df)
}



for (i in 1:1000){ # SEAK
  samp <- sample(4000,1)
  alpha <- sockeye_post$`mu_alpha[2]`[samp]
  sst <- sockeye_post$`mu_gamma[2]`[samp] # SST effect
  comp <- sockeye_post$`mu_kappa[2]`[samp] # comp effect
  sstXcomp <- sockeye_post$`mu_chi[2]`[samp] # interaction

  pred.cool <- alpha +(sst*sst_index[1]) + (comp*comp_index) + (sstXcomp*(comp_index*sst_index[1])) # cool years
  pred.warm <- alpha +(sst*sst_index[2]) + (comp*comp_index) + (sstXcomp*(comp_index*sst_index[2])) # cool years

  pred.df <- cbind(c(pred.cool,pred.warm),
                   c(comp_index,comp_index),
                   c(rep("cool",100),rep("warm",100)),
                   c(rep("SEAK",200)))

  pred.df <- as.data.frame(pred.df)
  pred.all.df <- rbind(pred.all.df,pred.df)
}

for (i in 1:1000){ # GoA
  samp <- sample(4000,1)
  alpha <- sockeye_post$`mu_alpha[2]`[samp]
  sst <- sockeye_post$`mu_gamma[3]`[samp] # SST effect
  comp <- sockeye_post$`mu_kappa[3]`[samp] # comp effect
  sstXcomp <- sockeye_post$`mu_chi[3]`[samp] # interaction

  pred.cool <- alpha +(sst*sst_index[1]) + (comp*comp_index) + (sstXcomp*(comp_index*sst_index[1])) # cool years
  pred.warm <- alpha +(sst*sst_index[2]) + (comp*comp_index) + (sstXcomp*(comp_index*sst_index[2])) # cool years

  pred.df <- cbind(c(pred.cool,pred.warm),
                   c(comp_index,comp_index),
                   c(rep("cool",100),rep("warm",100)),
                   c(rep("GoA",200)))

  pred.df <- as.data.frame(pred.df)
  pred.all.df <- rbind(pred.all.df,pred.df)
}

for (i in 1:1000){ # Bering Sea
  samp <- sample(4000,1)
  alpha <- sockeye_post$`mu_alpha[2]`[samp]
  sst <- sockeye_post$`mu_gamma[4]`[samp] # SST effect
  comp <- sockeye_post$`mu_kappa[4]`[samp] # comp effect
  sstXcomp <- sockeye_post$`mu_chi[4]`[samp] # interaction

  pred.cool <- alpha +(sst*sst_index[1]) + (comp*comp_index) + (sstXcomp*(comp_index*sst_index[1])) # cool years
  pred.warm <- alpha +(sst*sst_index[2]) + (comp*comp_index) + (sstXcomp*(comp_index*sst_index[2])) # cool years

  pred.df <- cbind(c(pred.cool,pred.warm),
                   c(comp_index,comp_index),
                   c(rep("cool",100),rep("warm",100)),
                   c(rep("Bering Sea",200)))

  pred.df <- as.data.frame(pred.df)
  pred.all.df <- rbind(pred.all.df,pred.df)
}

colnames(pred.all.df) <- c("prod","comp","sst","region")

summ.pred.all.sockeye.df <- pred.all.df |>
  mutate(prod = as.numeric(prod),
         comp = as.numeric(comp),
         spp = "Sockeye") |>
  group_by(comp, sst, region, spp) |>
  summarise(med.pred = median(prod),
            lwr.pred = quantile(prod, probs = 0.2, na.rm = TRUE),
            upr.pred = quantile(prod, probs = 0.8, na.rm = TRUE))

# chum ----
load(here("output", "models", "stat", "chum", "stat_inter.RData"), verbose=T)
chum_post <- as.data.frame(stat_inter)

sst_index <- c(-2,2)
comp_index <- seq(0,3,length.out=100)

pred.all.df <- data.frame(
  prod = numeric(0),
  comp = numeric(0),
  sst = numeric(0),
  region = character(0))


for (i in 1:1000){ # west coast
  samp <- sample(4000,1)
  alpha <- chum_post$`mu_alpha[1]`[samp]
  sst <- chum_post$`mu_gamma[1]`[samp] # SST effect
  comp <- chum_post$`mu_kappa[1]`[samp] # comp effect
  sstXcomp <- chum_post$`mu_chi[1]`[samp] # interaction

  pred.cool <- alpha +(sst*sst_index[1]) + (comp*comp_index) + (sstXcomp*(comp_index*sst_index[1])) # cool years
  pred.warm <- alpha +(sst*sst_index[2]) + (comp*comp_index) + (sstXcomp*(comp_index*sst_index[2])) # cool years

  pred.df <- cbind(c(pred.cool,pred.warm),
                   c(comp_index,comp_index),
                   c(rep("cool",100),rep("warm",100)),
                   c(rep("West Coast",200)))

  pred.df <- as.data.frame(pred.df)
  pred.all.df <- rbind(pred.all.df,pred.df)
}



for (i in 1:1000){ # GoA
  samp <- sample(4000,1)
  alpha <- chum_post$`mu_alpha[1]`[samp]
  sst <- chum_post$`mu_gamma[2]`[samp] # SST effect
  comp <- chum_post$`mu_kappa[2]`[samp] # comp effect
  sstXcomp <- chum_post$`mu_chi[2]`[samp] # interaction

  pred.cool <- alpha +(sst*sst_index[1]) + (comp*comp_index) + (sstXcomp*(comp_index*sst_index[1])) # cool years
  pred.warm <- alpha +(sst*sst_index[2]) + (comp*comp_index) + (sstXcomp*(comp_index*sst_index[2])) # cool years

  pred.df <- cbind(c(pred.cool,pred.warm),
                   c(comp_index,comp_index),
                   c(rep("cool",100),rep("warm",100)),
                   c(rep("GoA",200)))

  pred.df <- as.data.frame(pred.df)
  pred.all.df <- rbind(pred.all.df,pred.df)
}

for (i in 1:1000){ # Bering Sea
  samp <- sample(4000,1)
  alpha <- chum_post$`mu_alpha[1]`[samp]
  sst <- chum_post$`mu_gamma[3]`[samp] # SST effect
  comp <- chum_post$`mu_kappa[3]`[samp] # comp effect
  sstXcomp <- chum_post$`mu_chi[3]`[samp] # interaction

  pred.cool <- alpha +(sst*sst_index[1]) + (comp*comp_index) + (sstXcomp*(comp_index*sst_index[1])) # cool years
  pred.warm <- alpha +(sst*sst_index[2]) + (comp*comp_index) + (sstXcomp*(comp_index*sst_index[2])) # cool years

  pred.df <- cbind(c(pred.cool,pred.warm),
                   c(comp_index,comp_index),
                   c(rep("cool",100),rep("warm",100)),
                   c(rep("Bering Sea",200)))

  pred.df <- as.data.frame(pred.df)
  pred.all.df <- rbind(pred.all.df,pred.df)
}


colnames(pred.all.df) <- c("prod","comp","sst","region")

summ.pred.all.chum.df <- pred.all.df |>
  mutate(prod = as.numeric(prod),
         comp = as.numeric(comp),
         spp = "Chum") |>
  group_by(comp, sst, region, spp) |>
  summarise(med.pred = median(prod),
            lwr.pred = quantile(prod, probs = 0.2, na.rm = TRUE),
            upr.pred = quantile(prod, probs = 0.8, na.rm = TRUE))

# pink ----
load(here("output", "models", "stat", "pink", "stat_inter.RData"), verbose=T)
pink_post <- as.data.frame(stat_inter)

sst_index <- c(-2,2)
comp_index <- seq(0,3,length.out=100)

pred.all.df <- data.frame(
  prod = numeric(0),
  comp = numeric(0),
  sst = numeric(0),
  region = character(0))


for (i in 1:1000){ # west coast
  samp <- sample(4000,1)
  alpha <- pink_post$`mu_alpha[1]`[samp]
  sst <- pink_post$`mu_gamma[1]`[samp] # SST effect
  comp <- pink_post$`mu_kappa[1]`[samp] # comp effect
  sstXcomp <- pink_post$`mu_chi[1]`[samp] # interaction

  pred.cool <- alpha +(sst*sst_index[1]) + (comp*comp_index) + (sstXcomp*(comp_index*sst_index[1])) # cool years
  pred.warm <- alpha +(sst*sst_index[2]) + (comp*comp_index) + (sstXcomp*(comp_index*sst_index[2])) # cool years

  pred.df <- cbind(c(pred.cool,pred.warm),
                   c(comp_index,comp_index),
                   c(rep("cool",100),rep("warm",100)),
                   c(rep("West Coast",200)))

  pred.df <- as.data.frame(pred.df)
  pred.all.df <- rbind(pred.all.df,pred.df)
}



for (i in 1:1000){ # SEAK
  samp <- sample(4000,1)
  alpha <- pink_post$`mu_alpha[1]`[samp]
  sst <- pink_post$`mu_gamma[2]`[samp] # SST effect
  comp <- pink_post$`mu_kappa[2]`[samp] # comp effect
  sstXcomp <- pink_post$`mu_chi[2]`[samp] # interaction

  pred.cool <- alpha +(sst*sst_index[1]) + (comp*comp_index) + (sstXcomp*(comp_index*sst_index[1])) # cool years
  pred.warm <- alpha +(sst*sst_index[2]) + (comp*comp_index) + (sstXcomp*(comp_index*sst_index[2])) # cool years

  pred.df <- cbind(c(pred.cool,pred.warm),
                   c(comp_index,comp_index),
                   c(rep("cool",100),rep("warm",100)),
                   c(rep("SEAK",200)))

  pred.df <- as.data.frame(pred.df)
  pred.all.df <- rbind(pred.all.df,pred.df)
}

for (i in 1:1000){ # GoA
  samp <- sample(4000,1)
  alpha <- pink_post$`mu_alpha[1]`[samp]
  sst <- pink_post$`mu_gamma[3]`[samp] # SST effect
  comp <- pink_post$`mu_kappa[3]`[samp] # comp effect
  sstXcomp <- pink_post$`mu_chi[3]`[samp] # interaction

  pred.cool <- alpha +(sst*sst_index[1]) + (comp*comp_index) + (sstXcomp*(comp_index*sst_index[1])) # cool years
  pred.warm <- alpha +(sst*sst_index[2]) + (comp*comp_index) + (sstXcomp*(comp_index*sst_index[2])) # cool years

  pred.df <- cbind(c(pred.cool,pred.warm),
                   c(comp_index,comp_index),
                   c(rep("cool",100),rep("warm",100)),
                   c(rep("GoA",200)))

  pred.df <- as.data.frame(pred.df)
  pred.all.df <- rbind(pred.all.df,pred.df)
}

for (i in 1:1000){ # Bering Sea
  samp <- sample(4000,1)
  alpha <- pink_post$`mu_alpha[1]`[samp]
  sst <- pink_post$`mu_gamma[4]`[samp] # SST effect
  comp <- pink_post$`mu_kappa[4]`[samp] # comp effect
  sstXcomp <- pink_post$`mu_chi[4]`[samp] # interaction

  pred.cool <- alpha +(sst*sst_index[1]) + (comp*comp_index) + (sstXcomp*(comp_index*sst_index[1])) # cool years
  pred.warm <- alpha +(sst*sst_index[2]) + (comp*comp_index) + (sstXcomp*(comp_index*sst_index[2])) # cool years

  pred.df <- cbind(c(pred.cool,pred.warm),
                   c(comp_index,comp_index),
                   c(rep("cool",100),rep("warm",100)),
                   c(rep("Bering Sea",200)))

  pred.df <- as.data.frame(pred.df)
  pred.all.df <- rbind(pred.all.df,pred.df)
}

colnames(pred.all.df) <- c("prod","comp","sst","region")

summ.pred.all.pink.df <- pred.all.df |>
  mutate(prod = as.numeric(prod),
         comp = as.numeric(comp),
         spp = "Pink") |>
  group_by(comp, sst, region, spp) |>
  summarise(med.pred = median(prod),
            lwr.pred = quantile(prod, probs = 0.2, na.rm = TRUE),
            upr.pred = quantile(prod, probs = 0.8, na.rm = TRUE))



summ.pred.all.spp.df <- rbind(summ.pred.all.pink.df,summ.pred.all.chum.df,summ.pred.all.sockeye.df)

custom_colors <- c("cool" = "blue", "warm" = "red")


# Sockeye
sox.box<-summ.pred.all.spp.df|>
  filter(spp=="Sockeye")


s <- ggplot(data=sox.box, aes(x = comp, y = exp(med.pred), col=sst)) +
  geom_ribbon(aes(ymin=exp(lwr.pred), ymax=exp(upr.pred), col=sst, fill=sst), alpha=0.125) +
  geom_line(aes(col=sst),lwd=1) +
  theme_sleek() +
  xlab("") +
  ylab("Recruits-per-spawner") +
  scale_color_manual(values = custom_colors)+
  scale_fill_manual(values = custom_colors) +
  facet_grid(rows=vars(region), scales = "free") +
  labs(fill = "SST", title="Sockeye") +
  scale_x_continuous(breaks=c(0,1,2,3),labels=comp_label)+
  scale_y_continuous(n.breaks=4) +
  guides(color="none")+
  theme(legend.position="none",
        strip.text.y.right = element_blank(),
        plot.title = element_text(hjust=0.5),
        axis.title.y = element_text(size = 10))


# Chum
chum.box<-summ.pred.all.spp.df|>
  filter(spp=="Chum") |>
  ungroup() |>
  add_row(region = "SEAK")  # add empty row for SEAK

c <- ggplot(data=chum.box, aes(x = comp, y = exp(med.pred), col=sst)) +
  geom_ribbon(aes(ymin=exp(lwr.pred), ymax=exp(upr.pred), col=sst, fill=sst), alpha=0.125) +
  geom_line(aes(col=sst),lwd=1) +
  theme_sleek() +
  xlab("North Pacific competitor abundance (m)") +
  ylab("Recruits-per-spawner") +
  scale_color_manual(values = custom_colors)+
  scale_fill_manual(values = custom_colors) +
  facet_grid(rows=vars(region), scales = "free") +
  labs(fill = "SST", title="Chum") +
  scale_x_continuous(breaks=c(0,1,2,3),labels=comp_label)+
  scale_y_continuous(n.breaks=4) +
  guides(color="none")+
  theme(legend.position="none",
        strip.text.y.right = element_blank(),
        plot.title = element_text(hjust=0.5),
        axis.title.x = element_text(size = 10))


# Pink
pink.box<-summ.pred.all.spp.df|>
  filter(spp=="Pink")


p<- ggplot(data=pink.box, aes(x = comp, y = exp(med.pred), col=sst)) +
  geom_ribbon(aes(ymin=exp(lwr.pred), ymax=exp(upr.pred), col=sst, fill=sst), alpha=0.125) +
  geom_line(aes(col=sst),lwd=1) +
  theme_sleek() +
  xlab("") +
  ylab("Recruits-per-spawner") +
  scale_color_manual(values = custom_colors)+
  scale_fill_manual(values = custom_colors) +
  facet_grid(rows=vars(region), scales = "free") +
  labs(fill = "SST",  title="Pink") +
  scale_x_continuous(breaks=c(0,1,2,3),labels=comp_label)+
  scale_y_continuous(n.breaks=4) +
  guides(color="none")+
  theme(legend.justification = c(0, 0),
        legend.position = c(0.05, 0.35),
        legend.key.size = unit(11, "pt"),
        legend.background = element_blank(),
        legend.title = element_text( size = 8),
        plot.title = element_text(hjust=0.5))

full.interaction <- cowplot::plot_grid(s, c + theme(axis.title.y.left = element_blank()), p + theme(axis.title.y.left = element_blank()), ncol=3, rel_widths=c(1.1,1,1.1))+
  cowplot::draw_grob(
    grid::rectGrob(
      x = 0.5, y = 0.40, width = 0.3, height = 0.225,
      gp = grid::gpar(fill = "white", col = NA)))


png(here('figures/manuscript/main-text/sst-comp-interaction.png'), width=900*3, height=600*3, res=72*5)
print(full.interaction)
dev.off()
