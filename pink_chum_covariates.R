## Pink and Chum covariates

##   Load necessary data

## 1st year climate var 
raw.clim <- read.csv(file="data/sst_yr_1_stock_anomalies.csv", header=TRUE)
head(raw.clim)

## Pink competition
raw.comp <- read.csv(file="data-downloaded/competitor_indices_2023_05_16.csv", header = TRUE)
head(raw.comp)

# req'd columns # Obsolete now I think?
#bt.sock <- read.csv("./data/sockeye/master_brood_table.csv", header=T)
#r.cols <- grep("^R[[:digit:]]\\.[[:digit:]]", names(bt.sock), value = TRUE)

### PINK

## Master brood table
raw.clim.pink <- filter(raw.clim, Species=="Pink")
bt.raw.pink <- read.csv("./data/pink/master_pink_brood_table.csv", header=T)
bt.complete.pink <- bt.raw.pink[complete.cases(bt.raw.pink),]


## Climate index: SST at ocean entry point in 1st yr marine life
## SST during early marine life : no age weighting required
raw.clim.pink$BY <- raw.clim.pink$Year - 2 
early.sst.pink <- dplyr::left_join(bt.complete.pink, raw.clim.pink, by=c("Stock.ID", "BY"))

## competitor index: Pinks in 1st year marine life for pinks (2nd for all other sp)
## competitors in second year of marine life
raw.comp$BY <- raw.comp$Year - 2
np.pink.sec.pink <- dplyr::left_join(bt.complete.pink[,c("Stock.ID", "BY")], raw.comp[, c("pink_numbers_np", "BY")], by= "BY")


## Merge datasets 
master.pink <- dplyr::left_join(early.sst.pink, np.pink.sec.pink, by=c("BY","Stock.ID"))
master.pink.bt_w_cov1 <- master.pink
master.pink.bt_w_cov1 <- geographic.order(master.pink.bt_w_cov1) # Ordered factor
head(master.pink.bt_w_cov1)
summary(master.pink.bt_w_cov1)
unique(master.pink.bt_w_cov1$Stock.ID)


## Add derived columns
master.pink.bt_w_cov2 <- ddply(master.pink.bt_w_cov1, .(Stock), transform,
                          RS = R/S,
                          RS_stnd = scale(R/S)[ , 1],
                          lnRS = log(R/S),
                          S_stnd = scale(S)[ , 1],
                          early_sst_stnd = scale(sst_anomaly)[ , 1], 
                          np_pinks_sec_stnd = scale(pink_numbers_np)[ , 1])


## Fill in missing years that fall w/in min and max BY for each stock
master.pink.bt_w_cov3 <- fill.time.series(master.pink.bt_w_cov2) # do we want to use this for pinks?

# Export to output
pink <- master.pink.bt_w_cov3
write.csv(pink, "data/pink/master_pink_brood_table_covar.csv", row.names=FALSE)


### CHUM

#Data
raw.clim.chum <- filter(raw.clim, Species=="Chum")
bt.raw.chum <- read.csv("./data/chum/master_chum_brood_table.csv", header=T)
bt.complete.chum <- bt.raw.chum[complete.cases(bt.raw.chum),]

## Climate index: SST at ocean entry point in 1st yr marine life
## SST during early marine life : no age weighting required
raw.clim.chum$BY <- raw.clim.chum$Year - 2 
early.sst.chum <- dplyr::left_join(bt.complete.chum, raw.clim.chum, by=c("Stock.ID", "BY", "Species"))

## competitor index: pinks in 2nd year of marine life
## competitors in second year of marine life
raw.comp$BY <- raw.comp$Year - 3
np.pink.sec.chum <- dplyr::left_join(bt.complete.chum[,c("Stock.ID", "BY")], raw.comp[, c("pink_numbers_np", "BY")], by= "BY") # this isn't quite right


## Merge datasets 
master.chum <- dplyr::left_join(early.sst.chum, np.pink.sec.chum, by=c("BY","Stock.ID"))
## Stopped here because of duplication issue
master.chum.bt_w_cov1 <- master.chum
master.chum.bt_w_cov1 <- geographic.order(master.chum.bt_w_cov1) # Ordered factor
head(master.chum.bt_w_cov1)
summary(master.chum.bt_w_cov1)
unique(master.chum.bt_w_cov1$Stock.ID)


## Add derived columns
master.chum.bt_w_cov2 <- ddply(master.chum.bt_w_cov1, .(Stock), transform,
                               RS = R/S,
                               RS_stnd = scale(R/S)[ , 1],
                               lnRS = log(R/S),
                               S_stnd = scale(S)[ , 1],
                               early_sst_stnd = scale(sst_anomaly)[ , 1], 
                               np_pinks_sec_stnd = scale(pink_numbers_np)[ , 1])


## Fill in missing years that fall w/in min and max BY for each stock
master.chum.bt_w_cov3 <- fill.time.series(master.chum.bt_w_cov2) 


# Export to output
chum <- master.chum.bt_w_cov3
write.csv(chum, "data/chum/master_chum_brood_table_covar.csv", row.names=FALSE)

