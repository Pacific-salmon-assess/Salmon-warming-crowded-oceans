## Exploratory graphics and summaries for salmon and covariate data 

# Species
if(speciesFlag=="pink") {
  data_master <- pink
  info_master <- pink.info} else if (speciesFlag=="pinkeven"){
    data_master <- pinkeven
    info_master <- pinkeven.info } else if (speciesFlag=="pinkodd"){
      data_master <- pinkodd
      info_master <- pinkodd.info } else if (speciesFlag=="chum") {
        data_master <- chum
        info_master <- chum.info } else if(speciesFlag=="sockeye"){
          data_master <- sock
          info_master <- sock.info }


# Set paths to output locations - dependent on species 
fig.dir <- here("figures", "spp-explore", speciesFlag) # place to store figures generated in this script

# Make it if it doesn't exist
if(!dir.exists(fig.dir))  dir.create(fig.dir, recursive = T)


## Data summaries ------------------------------------------
head(data_master)
tail(data_master)
levels(data_master$Stock)
paste("brood tbl ordered geographically? -", all(levels(data_master$Stock) == unique(data_master$Stock))) 
paste("info tbl ordered geographically? -", all(levels(info_master$Stock) == unique(info_master$Stock))) 
paste("levels in brood tbl and info tbl the same? -", all(levels(info_master$Stock) == levels(data_master$Stock)))
summary(data_master)

length(unique(data_master$Stock))
min(data_master$BY)
max(data_master$BY)

data_fill <- fill.time.series(data_master)
stk.summary <- plyr::ddply(data_fill, .(Stock), summarize,
                stock.id = unique(Stock.ID),           
                min.yr = min(BY),
                max.yr = max(BY),
                n.total = length(BY),
                n.na = sum(is.na(R)),
                n.data  = sum(!is.na(R)))
stk.summary


## Avg prod x covar correlations ---------------------------
## Plot the average correlation of between productivity and the covariates
## grouping by ocean region.

cor.stock <- plyr::ddply(data_master, .(Stock.ID), plyr::summarize,
                         Ocean.Region2 = unique(Ocean.Region2),
                         early_sst = cor(lnRS, early_sst, use = "pairwise.complete.obs"),
                         np_pinks_sec = cor(lnRS, np_pinks_sec, use = "pairwise.complete.obs"))

cor.stock$Stock.ID <- NULL
cor.stock <- reshape2::melt(cor.stock, id.vars = "Ocean.Region2")

cor.ocean <- plyr::ddply(cor.stock, .(Ocean.Region2, variable), summarize,
                         cor.avg = mean(value))

pdf(here(fig.dir, "cor_prod_covar_bar.pdf"), width = 12, height = 9)
g <- barchart(cor.avg ~ variable, data = cor.ocean, groups = Ocean.Region2,
              origin = 0, par.settings = theme.mjm(),
              auto.key = list(space = "right"),
              ylab = "Average correlation",
              main = "Across stock average correlation b/w lnRS and covariate")
print(g)
dev.off()



## Avg covar correlations ----------------------------------
## Plot the average cross correlation of the specified covariate indices. The
## correlations are computed for each stock and then averaged across all stocks.

## specify covar names in data
covars <- c("early_sst",
            "np_pinks_sec")

## create empty 3D array
array.cor <- array(NA, dim = c(length(covars),
                                length(covars),
                                nrow(info_master)))


## calculate stock specific covar correlations
for(i in seq_along(info_master$Stock.ID)) {
    stk.i <- data_master[data_master$Stock.ID == info_master$Stock.ID[i], ]
    covar.i <- subset(stk.i, select = covars)
    cor.i <- cor(covar.i, use = "pairwise.complete.obs")
    array.cor[ , , i] <- cor.i
}

## Average across stocks
cor.covars <- apply(array.cor, c(1, 2), mean)
row.names(cor.covars) <- covars
colnames(cor.covars) <- covars

pdf(here(fig.dir, "cor_covar_matrix.pdf"), width = 12, height = 9)
cols <- chroma::dpal(500, hue = c(240, 0), chroma = 70, power = 1.0)
at.index <- seq(-1, 1, 0.1)
g <- levelplot(cor.covars, xlab = "", ylab = "",
               col.regions = cols,
               at = at.index,
               main = "Average stock-specific covariate correlations",
               panel = function(x, y, ...) {
                   panel.levelplot(x, y, ...)
                   panel.abline(v = seq(1.5, max(x) - 0.5), col = "grey50")
                   panel.abline(h = seq(1.5, max(y) - 0.5), col = "grey50")
               })
print(g)
dev.off()



## Map ocean entry -----------------------------------------
pdf(here(fig.dir, "map_ocean_entry.pdf"), width = 7, height = 6)
plot(info_master$lon, info_master$lat,
     xlim = c(-170, -120),
     ylim = c(45, 70),
     main = "Unique ocean entry locations",
     type = "n",
     las = 1,
     xlab = "Longitude",
     ylab = "Latitude")
plot(countriesLow, add = TRUE, border = "grey50")
points(info_master$lon, info_master$lat, pch = 16, col = as.factor(info_master$Ocean.Region2))
dev.off()


## Productivity time series --------------------------------
# Make sure Stock is ordered by latitude
pdf(here(fig.dir, "productivity_by_stock.pdf"), width = 19, height = 9)
g <- xyplot(R/S ~ BY | Stock, data = data_master,
            type = "l",
            scales = list(y = list(relation = "free")),
            ylab = "Productivity (R/S)",
            xlab = "Brood year",
            main = "Productivity by stock",
            par.settings = theme.mjm(),
            panel = function(x, y, ...) {
                panel.xyplot(x, y, ...)
            })
print(g)
dev.off()

pdf(here(fig.dir, "productivity_by_stock_stnd.pdf"), width = 19, height = 9)
g <- xyplot(RS_stnd ~ BY | Stock, data = data_master,
            type = "l",
            ylab = "Standardized productivity",
            xlab = "Brood year",
            main = "Standardized productivity by stock",
            par.settings = theme.mjm(),
            panel = function(x, y, ...) {
                panel.abline(h = 0, col = "grey50", lty = 2)
                panel.xyplot(x, y, ...)
            })
print(g)
dev.off()

pdf(here(fig.dir, "spawner_recruit.pdf"), width = 19, height = 9)
g <- xyplot(R ~ S | Stock, data = data_master,
            type = "p",
            pch = 16, cex = 0.5,
            scales = list(relation = "free"),
            ylab = "Recruits",
            xlab = "Spawners",
            main = "Recruits vs. spawners w/ fitted Ricker curve",
            par.settings = theme.mjm(),
            panel = function(x, y, ...) {
                panel.xyplot(x, y, ...)
                rk <- lm(log(y/x) ~ x)
                a  <- exp(coef(rk)[1])
                b  <- coef(rk)[2]
                panel.curve(a*x*exp(b*x),
                            from = min(x, na.rm = TRUE),
                            to = max(x, na.rm = TRUE))
            })
print(g)
dev.off()

pdf(here(fig.dir, "spawner_recruit_ln.pdf"), width = 19, height = 9)
g <- xyplot(lnRS ~ S | Stock, data = data_master,
            type = "p",
            pch = 16, cex = 0.5,
            scales = list(relation = "free"),
            ylab = "ln(R/S)",
            xlab = "Spawners",
            main = "ln(R/S) vs. spawners w/ fitted Ricker line",
            par.settings = theme.mjm(),
            panel = function(x, y, ...) {
                panel.xyplot(x, y, ...)
                panel.abline(lm(y ~ x))
            })
print(g)
dev.off()


## Detailed R/S plots for checking data quality ------------------

require(ggforce)

# make high R/S rule, i.e. points >50 get highlighted
data_fill <- data_fill %>% group_by(Stock) %>% 
  mutate(outlier = ifelse(RS>50, RS, NA))
# get right number of pages to put plots on
pg <- ceiling(nlevels(data_master$Stock)/18)

pdf(here(fig.dir, "productivity_detail.pdf"))
for(i in 1:pg){
g<- ggplot(data_fill) + geom_line(aes(x=BY, y=RS)) + 
    geom_point(aes(x=BY, y=outlier), col="red") +
    ggforce::facet_wrap_paginate(vars(Stock), scales="free_y", nrow=6, ncol=3, page=i) + theme_sleek()
print(g)
}
dev.off()



## SST indices ---------------------------------------------
pdf(here(fig.dir, "sst_early_index.pdf"), width = 19, height = 9)
g <- xyplot(early_sst ~ BY | Stock, data = data_master,
            type = "l",
            par.settings = theme.mjm(),
            xlab = "Brood year",
            ylab = "Early SST index",
            main = "Early SST index by stock",
            panel = function(x, y, ...) {
                panel.abline(h = 0, col = "grey50", lty = 2)
                panel.loess(x, y, col = "red2", lwd = 2)
                panel.xyplot(x, y, ...)
            })
print(g)
dev.off()




## Pink abundance indices ----------------------------------

pdf(here(fig.dir, "pink_index_second_year.pdf"), width = 19, height = 9)
g <- xyplot(np_pinks_sec ~ BY | Stock, data = data_master,
            type = "l",
            xlab = "Brood year",
            ylab = "Second year pink abundance index",
            main = "Second year pink abundance index by stock",
            par.settings = theme.mjm(),
            panel = function(x, y, ...) {
                panel.xyplot(x, y, ...)
                panel.loess(x, y, col = "red2", lwd = 2)
            })
print(g)
dev.off()



## Comparisons of covars -----------------------------------


pdf(here(fig.dir, "comparison_pink_second_sst_early.pdf"), width = 19, height = 9)
g <- xyplot(np_pinks_sec_stnd ~ early_sst_stnd | Stock, data = data_master,
            type = "p", col = 1, cex = 0.5,
            xlim = c(-3.25, 3.25),
            ylim = c(-3.25, 3.25),
            par.settings = theme.mjm(),
            xlab = "Early SST index",
            ylab = "Second year pink abundance index",
            main = "Comparison of pink second year and early SST indices",
            panel = function(x, y, ...) {
                panel.xyplot(x, y, ...)
                panel.loess(x, y, col = "red2", lwd = 2)
                panel.text(-2.75, 2.75,
                           paste("r =", round(cor(x, y, use = "pairwise.complete.obs"), 2)),
                           cex = 0.8, adj = 0, col = "black")
            })
print(g)
dev.off()



## Age structure -------------------------------------------
if(speciesFlag == "sockeye") {
  pdf(here(fig.dir, "age_ocean_entry_proportions.pdf"), width = 19, height = 9)
  oc <- subset(data_master, select = c("Stock", "BY", "ocean_0", "ocean_1", "ocean_2", "ocean_3", "ocean_4"))
  oc <- reshape2::melt(oc, id.vars = c("Stock", "BY"))
  g <- xyplot(value ~ BY | Stock, data = oc,
              groups = variable,
              auto.key = list(space = "right"),
              type = "l",
              par.settings = theme.mjm(),
              xlab = "Brood year",
              ylab = "Proportion of fish entering ocean",
              main = "Proportion of fish entering ocean in BY+",
              panel = function(x, y, ...) {
                  panel.xyplot(x, y, ...)
              })
  print(g)
  dev.off()
  
  
  pdf(here(fig.dir, "age_class_proportions.pdf"), width = 19, height = 14)
  oc <- subset(data_master, select = c("Stock", "BY", grep("^R[0-9]", names(data_master), value = TRUE)))
  oc <- reshape2::melt(oc, id.vars = c("Stock", "BY"))
  g <- levelplot(value ~ BY*variable | Stock, data = oc,
              type = "l",
              par.settings = theme.mjm(),
              xlab = "Brood year",
              ylab = "Age class",
              main = "Age class proportions",
              panel = function(x, y, ...) {
                  panel.grid(h = -18, v = -1)
                  panel.levelplot(x, y, ...)
              })
  print(g)
  dev.off()
}


## Bivariate -----------------------------------------------
sst.pos.neg <- ifelse(data_master$early_sst_stnd >= 0, "sst_pos", "sst_neg")
bi.panel <- function(x, y, ...) {
    panel.xyplot(x, y, ...)
    panel.abline(lm(y ~ x), ...)
}
pdf(here(fig.dir, "productivity_pinks_conditional_sst.pdf"), width = 19, height = 14)
g <- xyplot(lnRS ~ np_pinks_sec_stnd | Stock, data = data_master,
            groups = sst.pos.neg,
            type = "p",
            par.settings = theme.mjm(),
            xlab = "Second year pink index",
            ylab = "ln(R/S)",
            main = "Productivity vs second year pinks conditional on early SST",
            auto.key = list(space = "right"),
            panel = function(x, y, ...) {
                panel.superpose(x, y, panel.groups = bi.panel, ...)
            })
print(g)
dev.off()

##  for SST -----------------------------------------------
pink.pos.neg <- ifelse(data_master$np_pinks_sec_stnd >= 0, "pink_pos", "pink_neg")
bi.panel <- function(x, y, ...) {
  panel.xyplot(x, y, ...)
  panel.abline(lm(y ~ x), ...)
}
pdf(here(fig.dir, "productivity_sst_conditional_pink.pdf"), width = 19, height = 14)
g <- xyplot(lnRS ~ early_sst_stnd | Stock, data = data_master,
            groups = pink.pos.neg,
            type = "p",
            par.settings = theme.mjm(),
            xlab = "early SST anomaly",
            ylab = "ln(R/S)",
            main = "Productivity vs early sst conditional on second year pinks",
            auto.key = list(space = "right"),
            panel = function(x, y, ...) {
              panel.superpose(x, y, panel.groups = bi.panel, ...)
            })
print(g)
dev.off()


## Fig: Map ALL species (Hannah's version) --------------------


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
# colour
sp.col <- c("seagreen4", "palevioletred3", "orangered4")
names(sp.col) <- c("Chum", "Pink", "Sockeye")

col.dk <- rev(chroma::qpal(7, luminance = 20)[c(1, 3, 5, 7)])
names(col.dk) <- unique(info_master$ocean_region_lab)


map <- ggplot(map.info) + 
  geom_sf(data=na_map, color="grey40", fill="white", linewidth=0.1) + 
  ggspatial::geom_spatial_point(aes(x=lon, y=lat, col=Species, shape=Species), 
                                crs=4326, size=1.5, alpha=0.8, position=position_jitter(w=0.2, h=0.2)) +
  coord_sf(xlim=axes$xlims, ylim=axes$ylims) +
  scale_x_continuous(breaks=axes$xbreaks, labels=axes$xlabels) +
  scale_y_continuous(breaks=axes$ybreaks, labels=axes$ylabels) +
  #scale_colour_brewer(palette="Dark2") + 
  scale_colour_manual(values=sp.col, labels=c("Chum (n=44)", "Pink (n=70)", "Sockeye (n=52)")) +
  scale_shape_manual(values=c(15,16,17), labels=c("Chum (n=44)", "Pink (n=70)", "Sockeye (n=52)")) +
  labs(x="Longitude (°E)", y="Latitude (°N)") +
  theme_sleek() + 
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust=0.5),
        legend.position = c(0.85,0.75),
        legend.background = element_rect(colour="grey75"),
        aspect.ratio = 0.55) 

png(here('figures', 'spp-explore', "all_sp_map_pres.png"), res=288, height=1216, width=2192)
print(map)
dev.off()

ggsave(here('figures', 'spp-explore', "all_sp_map_pres.png"), map)
#pdf(here("figures", "spp-explore", "all_sp_map.pdf"), height=4, width=7)
#print(map)
#dev.off()


map.alt <- 
  ggplot(map.info) + 
  geom_sf(data=na_map, color="grey40", fill="white", linewidth=0.2) + 
  ggspatial::geom_spatial_point(aes(x=lon, y=lat, col=ocean_region_lab, fill=ocean_region_lab, shape=Species), 
                                crs=4326, size=1.2, alpha=0.6, position=position_jitter(w=0.5, h=0.2)) +
  coord_sf(xlim=axes$xlims, ylim=axes$ylims) +
  scale_x_continuous(breaks=axes$xbreaks, labels=axes$xlabels) +
  scale_y_continuous(breaks=axes$ybreaks, labels=axes$ylabels) +
  #scale_colour_brewer(palette="Dark2") + 
  scale_colour_manual(values=col.dk) +
  scale_fill_manual(values=col.region) +
  scale_shape_manual(values=c(22,21,24)) +
  labs(x="Longitude (°E)", y="Latitude (°N)") +
  theme_sleek() + 
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust=0.5),
        legend.position = "none",
        legend.background = element_rect(colour="grey75"),
        aspect.ratio = 0.7)

pdf(here("figures", "spp-explore", "all_sp_map_alt.pdf"), height=4, width=7)
print(map.alt)
dev.off()



## Fig: Time series of covariates

# 1) Competitor time series with raw data

comp.fig <- ggplot(raw.comp) + 
        geom_line(aes(x=Year, y=pink_numbers_np), col="darkred") +
        geom_vline(xintercept=c(1989,2011), color = "grey80", linetype = 1, linewidth = 0.25, alpha=0.8) +
        labs(x="Year", y="Pink abundance (millions)") +
        scale_y_continuous(limits=c(0,800), breaks = seq(0,750,250)) +
        theme_sleek() +
        theme(aspect.ratio=0.35)
# for presentation
comp.fig.pres <- comp.fig + theme(aspect.ratio=0.55)
png(here("figures", "spp-explore", "comp-pres.png"), width=550*2, height=300*2, res=72*4)
print(comp.fig.pres)
dev.off()
  

# 2) SST with raw data

# Load data
unique.oc.entry <- distinct(map.info, lat, lon, .keep_all=TRUE)
sst.anom <- read.csv("data/sst_raw_anomalies_extend.csv")
## Calculate average SST anomaly within area where stock spends first few months of marine life 
sst_anom <- sst.averager(unique.oc.entry, sst.anom, distance = 400)
sst_anom <- left_join(sst_anom, unique.oc.entry[,c("Stock.ID", "Ocean.Region2")], by=c("stock.id" = "Stock.ID"))

sst.fig <- 
  ggplot(sst_anom) +
  geom_line(aes(x=year, y=sst.anom, col=ocean_region_lab, group=stock.id), alpha=0.2) +
  geom_hline(aes(yintercept=0), linetype="dashed", col="gray50") +
  geom_vline(xintercept=c(1989,2011), color = "grey80", linetype = 1, linewidth = 0.25, alpha=0.8) +
  facet_grid(rows=vars(as.character(ocean_region_lab)), switch="y") +
  scale_colour_manual(values=col.region) +
  labs(x="Year", y="SST anomalies") +
  coord_cartesian(xlim=c(1950,2020)) +
  scale_y_continuous(breaks=seq(-3, 3, 3), labels=c("-3", " 0", " 3"), position="right") +
  theme_sleek() +
  theme(legend.position="none")



left <- cowplot::plot_grid(map.alt, comp.fig, ncol=1, rel_heights = c(1.7, 1), labels="auto", label_x = 0.075, label_y=1.03)
full.intro <- cowplot::plot_grid(left, sst.fig, ncol=2, rel_widths = c(1, .7), labels=list("", "c"))

png(here("figures", "spp-explore", "multi-intro.png"), height=721*1.5, width=1000*1.5, res=72*3)
print(full.intro)
dev.off()


