## ---- Making plots for project overview/background ------ ##

# Clear environment
#rm(list=ls())

# Dependencies 
library(pacman)
p_load(ggplot2, ggspatial, sf, dplyr, chroma)

## Read in stock.info
stock.info <- read.csv("data/all_stocks_info_may2023.csv")
sockeye.info <- read.csv("./data-downloaded/raw_stock_info_2023_05_19.csv")

# Clean data 
stock.info$lat[ which(stock.info$stock.id==144) ] <- 47.16 # fix a typo
stock.info$lon <- -abs(stock.info$lon) # add neg to Lon 

# Combine pink-odd and pink-even
stock.info$species <- gsub("-.*", "", stock.info$species)
stock.info$stock.name <- gsub("-Pink-.*", "-Pink", stock.info$stock.name)

# Make summary table that has # stocks at each ocean entry 
stock.map <- distinct(stock.info[,c("species", "stock.name", "lat", "lon")])
stock.map <- stock.map %>% group_by(species, lat, lon) %>% summarize(n.stk = n())

## Define colors
col.lt <- rev(chroma::qpal(7)[c(1, 4, 6)])
col.dk <- rev(chroma::qpal(7, luminance = 30)[c(1, 4, 6)])
col.sp <- rev(chroma::qpal(5))
shp.region <- c("WC" = 21, "GOA" = 22, "BS" = 24)

## Labels and axes
axes <- list( xlims=c(-165.5, -121), 
              ylims=c(44.5, 66.5),
              xbreaks=seq(-160,-120,10), 
              xlabels=as.character(seq(-160,-120,10)),
              seq(45, 65, 5), 
              ylabels=as.character(seq(45,65,5)))
reg.labs <- c("BS" = "Bering Sea", "GOA" = "Gulf of Alaska", "WC" = "West Coast")

#species index 
sp.ind <- unique(stock.info$species)


## Load maps 
cl <- rnaturalearth::ne_states(country = c("United States of America", "Canada"))
na_map <- st_as_sf(cl)




## ---- MAP with no regions, point size = # stocks at ocean entry ---- ##

breaks <- list( sock = c(1, 10, 20),
                chum = c(1,3),
                pink = c(1),
                chin = c(1, 15, 30),
                coho = c(1, 5) )
pt.sz <- list( sock = c(2.5, 8),
               chum = c(2.5, 5),
               pink = 2.5, 
               chin = c(2.5, 7),
               coho = c(2.5, 5.5) )
  

for(i in 1:length(sp.ind)) {
  
  n.tot <- sum(stock.map$n.stk[which(stock.map$species==sp.ind[i])])
  g <- ggplot() + 
    geom_sf(data=na_map, color="grey50", fill="white", linewidth=0.3) + 
    geom_spatial_point(data = filter(stock.map, species==sp.ind[i]),
                       aes(x=lon, y=lat, size=n.stk), 
                       crs=4326, alpha=0.8, colour=col.sp[i]) +
    coord_sf(xlim=axes$xlims, ylim=axes$ylims) +
    scale_x_continuous(breaks=axes$xbreaks, labels=axes$xlabels) +
    scale_y_continuous(breaks=axes$ybreaks, labels=axes$ylabels) +
    scale_size_continuous(range=pt.sz[[i]], name="No. of Stocks", breaks=breaks[[i]]) +
    geom_text(data=NULL, aes(x=-153, y=45, label=paste("Total N = ", n.tot))) +
    labs(x="Longitude (°E)", y="Latitude (°N)") +
    theme_bw() + 
    theme( panel.grid = element_blank(),
           plot.title = element_text(hjust=0.5),
           legend.position = c(0.3, 0.25)
           )

  ggsave(plot=g, paste0("./figures/background-presentation/",sp.ind[i], "_map.png"), width=5, height=5, dpi=600)
  
}



## ---- SST Timeseries --- ##

# load data - SST
sst.dat <- read.csv("./data/sst_yr_1_stock_anomalies.csv")
ind.stks <- as.character(c( 102, 160, 142))
sst.regs <- c("102"="Southern BC", "160"="Bering Sea", "142"="Gulf of Alaska")
sst.plot.dat <- filter(sst.dat, Stock.ID %in% as.integer(ind.stks))

# multipanel plot:
g <- ggplot(sst.plot.dat, aes(x=Year, y=sst_anomaly, col=as.factor(Stock.ID))) + 
  labs(x="Year", y="SST Anomaly") +
  scale_colour_manual(values=col.dk, labels=sst.regs) +
  geom_line(linewidth=1) + facet_wrap(vars(as.factor(Stock.ID)), nrow=1, labeller=as_labeller(sst.regs)) + 
  theme_bw() + theme(panel.grid = element_blank(),
                     legend.position = "none", 
                     strip.background = element_rect(fill="white", colour="white"),
                     strip.text = element_text(size=12, face="bold", colour="gray20"))
ggsave(plot=g, "./figures/background-presentation/SST_timeseries.png", width=6, height=2.5, dpi=500)




## -- NOT USED -- map with ocean groupings

#for (i in 1: length(sp.ind)) {
g <- ggplot() + geom_sf(data=na_map, color="grey50", fill="white", linewidth=0.3) + 
  geom_spatial_point(data = filter(stock.info, species=="Chum"), #sp.ind[i]), 
                     aes(x=lon, y=lat, colour=ocean.basin, fill=ocean.basin, shape=ocean.basin), 
                     crs=4326, size=3, stroke=1.75) +
  coord_sf(xlim=axes$xlims, ylim=axes$ylims) +
  scale_x_continuous(breaks=axes$xbreaks, labels=axes$xlabels) +
  scale_y_continuous(breaks=axes$ybreaks, labels=axes$ylabels) +
  scale_fill_manual(values=col.lt, 
                    labels=reg.labs) + 
  scale_colour_manual(values=col.dk, 
                      labels=reg.labs) + 
  scale_shape_manual(values=shp.region,
                     labels=reg.labs) + 
  labs(x="Longitude (°E)", y="Latitude (°N)", title="Chum") + 
  theme_bw() + 
  theme( panel.grid = element_blank(),
         legend.position = c(0.25, 0.25),
         legend.title = element_blank(), 
         plot.title = element_text(hjust=0.5)
  )
# }
print(g)
dev.off()

