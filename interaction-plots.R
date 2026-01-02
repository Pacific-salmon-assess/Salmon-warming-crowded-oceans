library(ggsidekick) #for theme_sleek() - doesn't work with some vs of R, hence the comment

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

summ.pred.all.df <- pred.all.df |>
  mutate(prod = as.numeric(prod),
         comp = as.numeric(comp)) |>
  group_by(comp, sst, region) |>
  summarise(med.pred = median(prod),
            lwr.pred = quantile(prod, probs = 0.2, na.rm = TRUE),
            upr.pred = quantile(prod, probs = 0.8, na.rm = TRUE))

custom_colors <- c("cool" = "blue", "warm" = "red")
comp_label <- c(400,500,600,700)

g <- ggplot(data = summ.pred.all.df, aes(x = comp, y = exp(med.pred), col=sst)) +
  geom_ribbon(aes(ymin=exp(lwr.pred), ymax=exp(upr.pred), col=sst, fill=sst), alpha=0.2) +
  geom_line(aes(col=sst),lwd=1.05) +
  theme_sleek() +
  xlab("North Pacific pink abundance (millions)") +
  ylab("Recruits-per-spawner") +
  scale_color_manual(values = custom_colors)+
  scale_fill_manual(values = custom_colors) +
  facet_wrap(~region, ncol=2, scales="free_y") +
  labs(fill = "SST") +
  scale_x_continuous(breaks=c(0,1,2,3),labels=comp_label)+
  guides(color="none") +
  theme(legend.justification = c(0, 0),
        legend.position = c(0.05, 0.35),
        legend.key.size = unit(11, "pt"),
        legend.background = element_blank())

png(here('figures', 'spp-explore', "sockeye-sst-inter.png"), units = "in", width = 6, height = 6, res = 500)
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
comp_label <- c(400,500,600,700)

g <- ggplot(data = summ.pred.all.df, aes(x = comp, y = exp(med.pred), col=sst)) +
  geom_ribbon(aes(ymin=exp(lwr.pred), ymax=exp(upr.pred), col=sst, fill=sst), alpha=0.2) +
  geom_line(aes(col=sst),lwd=1.05) +
  theme_sleek() +
  xlab("North Pacific pink abundance (millions)") +
  ylab("Recruits-per-spawner") +
  scale_color_manual(values = custom_colors)+
  scale_fill_manual(values = custom_colors) +
  facet_wrap(~region, scales = "free_y",ncol=2) +
  labs(fill = "SST") +
  scale_x_continuous(breaks=c(0,1,2,3),labels=comp_label)+
  guides(color="none")+
  theme(legend.justification = c(0, 0),
        legend.position = c(0.05, 0.35),
        legend.key.size = unit(11, "pt"),
        legend.background = element_blank())

png(here('figures', 'spp-explore', "chum-sst-inter.png"), units = "in", width = 6, height = 6, res = 500)
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
comp_label <- c(400,500,600,700)

g <- ggplot(data = summ.pred.all.df, aes(x = comp, y = exp(med.pred), col=sst)) +
  geom_ribbon(aes(ymin=exp(lwr.pred), ymax=exp(upr.pred), col=sst, fill=sst), alpha=0.2) +
  geom_line(aes(col=sst),lwd=1.05) +
  theme_sleek() +
  xlab("North Pacific pink abundance (millions)") +
  ylab("Recruits-per-spawner") +
  scale_color_manual(values = custom_colors)+
  scale_fill_manual(values = custom_colors) +
  facet_wrap(~region, scales = "free_y", ncol=2) +
  labs(fill = "SST") +
  scale_x_continuous(breaks=c(0,1,2,3),labels=comp_label)+
  guides(color="none")+
  theme(legend.justification = c(0, 0),
        legend.position = c(0.05, 0.35),
        legend.key.size = unit(11, "pt"),
        legend.background = element_blank())

png(here('figures', 'spp-explore', "pink-sst-inter.png"), units = "in", width = 6, height = 6, res = 500)
print(g)
dev.off()
