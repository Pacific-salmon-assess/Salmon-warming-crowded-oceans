## (6) Exploratory graphics and summaries for sockeye and covariate data ## 
## --------------------------------------------------------------------- ##


if(!dir.exists("./figures/sock-explore/"))
    dir.create("./figures/sock-explore/")


## Data summaries ------------------------------------------
head(sock)
tail(sock)
levels(sock$Stock)
summary(sock)

length(unique(sock$Stock))
min(sock$BY)
max(sock$BY)

plyr::ddply(sock, .(Stock), summarize,
            min.yr = min(BY),
            max.yr = max(BY),
            n.total = length(BY),
            n.na = sum(is.na(R)),
            n.data  = sum(!is.na(R)))



## Avg prod x covar correlations ---------------------------
## Plot the average correlation of between productivity and the covariates
## grouping by ocean region.

cor.stock <- plyr::ddply(sock, .(Stock.ID), summarize,
                         Ocean.Region = unique(Ocean.Region),
                         early_sst = cor(lnRS, early_sst, use = "pairwise.complete.obs"),
                         np_pinks_sec = cor(lnRS, np_pinks_sec, use = "pairwise.complete.obs"))

cor.stock$Stock.ID <- NULL
cor.stock <- reshape2::melt(cor.stock, id.vars = "Ocean.Region")

cor.ocean <- plyr::ddply(cor.stock, .(Ocean.Region, variable), summarize,
                         cor.avg = mean(value))

pdf("./figures/sock-explore/cor_prod_covar_bar.pdf", width = 12, height = 9)
g <- barchart(cor.avg ~ variable, data = cor.ocean, groups = Ocean.Region,
              origin = 0, par.settings = theme.mjm(),
              auto.key = list(space = "right"),
              ylab = "Average correlation",
              main = "Across stock average correlation b/w lnRS and covariate")
print(g)
dev.off()



## Avg covar correlations ----------------------------------
## Plot the average cross correlation of the specified covariate indices. The
## correlations are computed for each stock and then averaged across all stocks.

## specify covar names in `sock`
covars <- c("early_sst",
            "np_pinks_sec")

## create empty 3D array
array.cor <- array(NA, dim = c(length(covars),
                                length(covars),
                                nrow(sock.info)))


## calculate stock specific covar correlations
for(i in seq_along(sock.info$Stock.ID)) {
    sock.i <- sock[sock$Stock.ID == sock.info$Stock.ID[i], ]
    covar.i <- subset(sock.i, select = covars)
    cor.i <- cor(covar.i, use = "pairwise.complete.obs")
    array.cor[ , , i] <- cor.i
}

## Average across stocks
cor.covars <- apply(array.cor, c(1, 2), mean)
row.names(cor.covars) <- covars
colnames(cor.covars) <- covars

pdf("./figures/sock-explore/cor_covar_matrix.pdf", width = 12, height = 9)
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
pdf("./figures/sock-explore/map_ocean_entry.pdf", width = 7, height = 6)
plot(sock.info$lon, sock.info$lat,
     xlim = c(-170, -120),
     ylim = c(45, 70),
     main = "Unique ocean entry locations",
     type = "n",
     las = 1,
     xlab = "Longitude",
     ylab = "Latitude")
plot(countriesLow, add = TRUE, border = "grey50")
points(sock.info$lon, sock.info$lat, pch = 16, col = "red2")
dev.off()


## Productivity time series --------------------------------
pdf("./figures/sock-explore/productivity_by_stock.pdf", width = 19, height = 9)
g <- xyplot(R/S ~ BY | Stock, data = sock,
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

pdf("./figures/sock-explore/productivity_by_stock_stnd.pdf", width = 19, height = 9)
g <- xyplot(RS_stnd ~ BY | Stock, data = sock,
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

pdf("./figures/sock-explore/spawner_recruit.pdf", width = 19, height = 9)
g <- xyplot(R ~ S | Stock, data = sock,
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

pdf("./figures/sock-explore/spawner_recruit_ln.pdf", width = 19, height = 9)
g <- xyplot(lnRS ~ S | Stock, data = sock,
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



## SST indices ---------------------------------------------
pdf("./figures/sock-explore/sst_early_index.pdf", width = 19, height = 9)
g <- xyplot(early_sst ~ BY | Stock, data = sock,
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

pdf("./figures/sock-explore/pink_index_second_year.pdf", width = 19, height = 9)
g <- xyplot(np_pinks_sec ~ BY | Stock, data = sock,
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


pdf("./figures/sock-explore/comparison_pink_second_sst_early.pdf", width = 19, height = 9)
g <- xyplot(np_pinks_sec_stnd ~ early_sst_stnd | Stock, data = sock,
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
pdf("./figures/sock-explore/age_ocean_entry_proportions.pdf", width = 19, height = 9)
oc <- subset(sock, select = c("Stock", "BY", "ocean_0", "ocean_1", "ocean_2", "ocean_3", "ocean_4"))
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

pdf("./figures/sock-explore/age_class_proportions.pdf", width = 19, height = 14)
oc <- subset(sock, select = c("Stock", "BY", grep("^R[0-9]", names(sock), value = TRUE)))
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



## Bivariate -----------------------------------------------
sst.pos.neg <- ifelse(sock$early_sst_stnd >= 0, "sst_pos", "sst_neg")
bi.panel <- function(x, y, ...) {
    panel.xyplot(x, y, ...)
    panel.abline(lm(y ~ x), ...)
}
pdf("./figures/sock-explore/productivity_pinks_conditional_sst.pdf", width = 19, height = 14)
g <- xyplot(lnRS ~ np_pinks_sec_stnd | Stock, data = sock,
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
pink.pos.neg <- ifelse(sock$np_pinks_sec_stnd >= 0, "pink_pos", "pink_neg")
bi.panel <- function(x, y, ...) {
  panel.xyplot(x, y, ...)
  panel.abline(lm(y ~ x), ...)
}
pdf("./figures/sock-explore/productivity_sst_conditional_pink.pdf", width = 19, height = 14)
g <- xyplot(lnRS ~ early_sst_stnd | Stock, data = sock,
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
