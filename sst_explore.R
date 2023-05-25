## (4) R script to explore the SST anomalies ##
## ----------------------------------------- ##


if(!dir.exists("./figures/sst-explore/"))
    dir.create("./figures/sst-explore/")

sst.dat <- read.csv("./data/sst_raw_anomalies.csv")
head(sst.dat)
tail(sst.dat)
sapply(sst.dat, class)
summary(sst.dat)


## Plot centers of SST grid cells
dat1 <- sst.dat[sst.dat$year == 2015 & sst.dat$month == 3, ]
pdf("./figures/sst-explore/sst_grid_lon.pdf", width = 8, height = 6)
plot(dat1$lon, dat1$lat, pch = 16)
maps::map('world2', add = TRUE, col = "grey50")
dev.off()


## Plot monthly SST anomalies
## These files are large and are ignored by git.
cols <- chroma::dpal(500, hue = c(240, 0), chroma = 70, power = 1.0)
at <- seq(-5.5, 4.5, 0.5)
dir  <- "./figures/sst-explore/sst-anomaly-maps/"
sst.map(data = sst.dat,
        plot.dir = dir,
        col.regions = cols,
        at = at)
