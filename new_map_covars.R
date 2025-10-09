# Create NEW map /covars figure
library(sf)
library(cowplot)

# Downlad map and convert to sp
na_map <- rnaturalearth::ne_countries(country = c("United States of America", "Canada"), scale='medium', returnclass="sf")


axes <- list( xlims=c(-167, -121),
              ylims=c(46, 67),
              xbreaks=seq(-160,-120,10),
              xlabels=as.character(seq(-160,-120,10)),
              seq(45, 65, 5),
              ybreaks=seq(45, 65, 5),
              ylabels=as.character(seq(45,65,5)))

#make map data
sock.info$Species <- "Sockeye"
chum.info$Species <- "Chum"
pink.info$Species <- "Pink"
map.info <- bind_rows(sock.info, chum.info, pink.info)
map.info <- ocean_region_lab(map.info)
# colour
sp.col <- c("seagreen4", "palevioletred3", "orangered4")
names(sp.col) <- c("Chum", "Pink", "Sockeye")

col.dk <- rev(chroma::qpal(7, luminance = 20)[c(1, 3, 5, 7)])
names(col.dk) <- unique(map.info$ocean_region_lab)


# Map with Albers projection
na_map_albers <- sf::st_transform(na_map, crs=3579) # transform map
#transform points:
points_albers <- sf::st_as_sf(map.info, coords=c("lon", "lat"), crs="EPSG:4326")
points.transformed <- sf::st_transform(points_albers,
                                       crs="EPSG:3579")
st_crs(points.transformed)
# transform axis limits:
lims <- data.frame(x=axes$xlims, y=axes$ylims)
lims_proj <- st_as_sf(lims, coords=c("x", "y"), crs="EPSG:4326")
lims_alb <- st_transform(lims_proj, crs="EPSG:3579")

legend_labels <- c("West Coast (WC)", "Southeast Alaska (SEAK)",
                   "Gulf of Alaska (GoA)", "Bering Sea (BS)")

map_albers <-
  ggplot(points.transformed) +
  geom_sf(data=na_map_albers, fill="grey85", colour="grey70",
          linewidth=0.2) +
  geom_sf(aes(col=ocean_region_lab, fill=ocean_region_lab, shape=Species),size=1.2, alpha=0.6) +
  coord_sf(xlim=axes$xlims, ylim=axes$ylims, default_crs = st_crs(4326)) +
  scale_x_continuous(breaks=axes$xbreaks, labels=axes$xlabels) +
  scale_y_continuous(breaks=axes$ybreaks, labels=axes$ylabels) +
  scale_colour_manual(values=col.dk,
                      guide=guide_legend(reverse=T,
                                         override.aes =
                                           list(size=3,
                                           alpha=.8,
                                           shape="square")),
                      labels = legend_labels) +
  scale_fill_manual(values=col.region, guide="none") +
  scale_shape_manual(values=c(22,21,24), guide="none") +
  labs(x="Longitude (°E)", y="Latitude (°N)", col="") +
  theme_sleek() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.3,0.3),
        legend.title = element_blank(),
        legend.box.background = element_rect(fill="grey95",
                                             colour="transparent"),
        aspect.ratio = 0.7)

# SST anomalies - multipanel
# Load data
unique.oc.entry <- distinct(map.info, lat, lon, .keep_all=TRUE)
sst.anom <- read.csv("data/sst_raw_anomalies_extend.csv")
## Calculate average SST anomaly within area where stock spends first few months of marine life
sst_anom <- sst.averager(unique.oc.entry, sst.anom, distance = 400)
sst_anom <- left_join(sst_anom, unique.oc.entry[,c("Stock.ID", "Ocean.Region2")], by=c("stock.id" = "Stock.ID"))
sst_anom <- bind_rows(data.frame(Ocean.Region2 = rep(unique(sst_anom$Ocean.Region2), each=2),
           year = 1990,
           lims = c(-2,2,-5,5,-5,5,-5,5)), sst_anom) # add points to hack scales

sst.fig <-
  ggplot(ocean_region_lab(sst_anom)) +
  geom_line(aes(x=year, y=sst.anom, col=ocean_region_lab, group=stock.id), alpha=0.2) +
  #geom_point(aes(x=year, y=lims), col="transparent") +
  geom_hline(aes(yintercept=0), linetype=1, col="gray50") +
  geom_vline(xintercept=c(1989,2011), color = "grey80", linetype = 2, linewidth = 0.25, alpha=0.8) +
  facet_grid(rows=vars(as.character(ocean_region_lab)),
             scales="free_y", labeller =
               as_labeller(c("West Coast" = "WC",
                             "Southeast Alaska" = "SEAK",
                             "Gulf of Alaska" = "GoA",
                             "Bering Sea" = "BS"))) +
  scale_colour_manual(values=col.region) +
  labs(x="Year", y="SST anomalies") +
  coord_cartesian(xlim=c(1950,2020)) +
  scale_y_continuous(breaks=seq(-4,4,1), labels=c("-4", "", "-2", "", " 0", " ", " 2", "", " 4"), position="left") +
  theme_sleek() +
  theme(legend.position="none",
        strip.text.y.right = element_text(angle=0, hjust = 0),
        plot.margin = margin(r=-2, l=10),
        axis.text.y.left = element_text(margin=margin(l=-2)))

# Pink salmon timeseries
comp.fig <- ggplot(raw.comp) +
  geom_line(aes(x=Year, y=pink_numbers_np), col="grey20") +
  geom_vline(xintercept=c(1989,2011), color = "grey80", linetype = 2, linewidth = 0.25, alpha=0.8) +
  labs(x="Year", y="Pink salmon\n(millions)") +
  scale_y_continuous(limits=c(0,800), breaks = seq(0,750,250),
                     position="left") +
  scale_x_continuous(breaks=seq(1960,2020,20)) +
  theme_sleek() +
  theme(axis.text=element_text(size=7),
        axis.title = element_text(size=12),
        plot.margin = margin(l=25, r=5, unit = "pt"))


# Raw SST boxplot horizontal
sst.anom <- read.csv("data/sst_raw_anomalies_extend.csv")
sst_full <- sst.averager(map.info, sst.anom, distance=400)
sst_full <- left_join(sst_full, map.info[,c("Stock.ID", "Ocean.Region2")], by=c("stock.id" = "Stock.ID"))

sst_mean <- sst_full %>% dplyr::summarize(ts_mean=mean(sst), .by=c(stock.id, Ocean.Region2))
sst_raw_inset_hor <- ggplot(ocean_region_lab(sst_mean)) +
  geom_boxplot(aes(x=ts_mean, y=ocean_region_lab, fill=ocean_region_lab), alpha=0.8) +
  scale_fill_manual(values=col.region, guide="none") +
  scale_y_discrete(position="right",
                   labels=c("West \n Coast",
                            "Southeast \n Alaska",
                            "Gulf of \n Alaska",
                            "Bering \n Sea")) +
  theme_sleek() +
  labs(x="Mean SST (°C)", y="") +
  theme(plot.background=element_rect(fill="transparent"),
        axis.title=element_text(size=10),
        axis.title.y = element_blank(),
        plot.margin = margin(l=18, t=10))

# Raw SST dot-and-whiskers
sst_raw_dot <- sst_full |>
  group_by(Ocean.Region2) |>
  summarize(q2.5 = quantile(sst, 0.025),
            q25 = quantile(sst, 0.25),
            q50 = median(sst, 0.5),
            q75 = quantile(sst, 0.75),
            q97.5 = quantile(sst, 0.975)) |>
  ocean_region_lab() |>
  mutate(ocean_region_lab =
           factor(ocean_region_lab,
             levels=rev(levels(ocean_region_lab)))) |>
  ggplot() +
  geom_segment(aes(x=q2.5, xend=q97.5,
                   y=ocean_region_lab, yend=ocean_region_lab,
                   col=ocean_region_lab), alpha=0.5) +
  geom_segment(aes(x=q25, xend=q75,
                   y=ocean_region_lab, yend=ocean_region_lab,
                   col=ocean_region_lab), linewidth=1.5, alpha=0.5) +
  geom_point(aes(x=q50, y=ocean_region_lab, col=ocean_region_lab),
             size=2) +
  scale_x_continuous(breaks=seq(6, 12, 2)) +
  scale_y_discrete(labels=c("WC", "SEAK", "GoA", "BS"), position = "right") +
  scale_colour_manual(values=col.region) +
  theme_sleek() + labs(x="SST (°C)", y="") +
  theme(legend.position = "none",
        plot.margin = margin(l=30, t=10, r=-5, unit="pt"),
        axis.text.y.right = element_text(margin = margin(l=4, r=-7, unit="pt")),
        axis.title = element_text(size=10))

# Draw general legend for the plot
lgnd.dat <- data.frame(region = unique(map.info$ocean_region_lab),
                       name = c("West \n Coast",
                                "Southeast \n Alaska",
                                "Gulf of \n Alaska",
                                "Bering \n Sea"))

lgnd2 <- lgnd.dat |> ggplot() +
  geom_tile(aes(x=.5, y=.5, fill=region, col=region),
            linewidth=4, linejoin="round", alpha=0.5) +
  facet_grid(rows=vars(name)) +
  scale_fill_manual(values=col.region, guide="none") +
  scale_colour_manual(values=col.dk, guide="none") +
  labs(x="", y="") +
  theme(panel.background=element_rect(fill="white"),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        strip.background = element_rect(fill="white"),
        strip.text.y.right = element_text(angle=0, size=7,
                                          hjust = 0,
                                          margin=margin(l=10, unit="pt")),
        plot.title = element_text(hjust=0.5),
        aspect.ratio = 1)


blank <- ggplot() + theme(panel.background = element_rect(fill="white")) # blank grob to put legend on
legend <- cowplot::ggdraw(blank) + cowplot::draw_plot(lgnd2)



# arrangement 1
right <- cowplot::plot_grid(sst.fig, sst_raw_dot, rel_heights=c(2.5,1), nrow=2, labels=c("c", "d"), hjust=-1.5) # horizonal justification doesn't work either as
left <- cowplot::plot_grid(map_albers, (comp.fig + theme_sleek()), ncol=1, rel_heights = c(2.5, 1), labels="auto", hjust=-1.7, vjust=c(1.5, .7))
intro.plot.2025 <- cowplot::plot_grid(left, right, ncol=2, rel_widths = c(1.7, 1))
png(here("figures", "spp-explore", "multi-intro-2025-al.png"), height=721*1.5, width=1000*1.5, res=72*3)
print(intro.plot.2025)
dev.off()


# Alternative arrangement 2
right <- cowplot::plot_grid(sst.fig, comp.fig, rel_heights=c(2.5,1), nrow=2, labels=c("c", "d"))

test <- cowplot::plot_grid(left, right, ncol=2, rel_widths = c(1.8, 1))

png(here("figures", "spp-explore", "multi-intro-2025-al2.png"), height=721*1.5, width=1000*1.5, res=72*3)
print(test)
dev.off()

# Arrangement 3
top <- plot_grid(map_albers, sst.fig, ncol=2, labels="auto", rel_widths = c(1.8, 1))
bottom <- plot_grid(comp.fig, (sst_raw_inset_hor + theme(axis.text = element_blank())), legend, nrow=1, labels = c("c", "d", ""), rel_widths=c(1.5,1.5,1))

test3 <- plot_grid(top, bottom, nrow=2, rel_heights=c(1,0.4))

png(here("figures", "spp-explore", "multi-intro-2025-al3.png"), height=721*1.5, width=1000*1.5, res=72*3)
print(test3)
dev.off()
