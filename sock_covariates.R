## (5) R scripts to create the stock specific covariates used in the analysis ##
## -------------------------------------------------------------------------- ##

# A simplified version of covariates script, including only:
# first year SST, competitors in 2nd year

##   Load necessary data

## Master brood table
bt.raw <- read.csv("./data/master_brood_table.csv", header=T)
bt.complete <- bt.raw[complete.cases(bt.raw),]
head(bt.complete)
tail(bt.complete)
sapply(bt.complete, class)

## 1st year climate var 
raw.clim <- read.csv(file="data/sst_yr_1_stock_anomalies.csv",header=TRUE)
head(raw.clim)

## Pink competition
raw.comp <- read.csv(file="data-downloaded/competitor_indices_2023_05_16.csv", header = TRUE)
head(raw.comp)


## Age weighted climate index: SST at ocean entry point in 1st yr marine life

## SST during early marine life 
early.sst <- clim.wgt.avg(brood.table = bt.complete,
                          env.data = raw.clim,
                          env.covar = "sst_anomaly",
                          type = "first_year",
                          out.covar = "early_sst")


## Age weighted competitor index: Pinks in 2nd yr marine life

## competitors in second year of marine life
np.pink.sec <- pink.wgt.avg(brood.table = bt.complete,
                            pink.data = raw.comp,
                            pink.covar = "pink_numbers_np",
                            type = "second_year",
                            out.covar = "np_pinks_sec")


## Merge datasets 
master <- dplyr::left_join(bt.complete, early.sst, by=c("BY","Stock.ID"))
master <- dplyr::left_join(master, np.pink.sec, by=c("BY","Stock.ID"))
master.bt_w_cov1 <- master
master.bt_w_cov1$Stock <- geographic.order(master.bt_w_cov1) # Ordered factor
head(master.bt_w_cov1)
tail(master.bt_w_cov1)
summary(master.bt_w_cov1)
sapply(master.bt_w_cov1, class)
unique(master.bt_w_cov1$Stock.ID)


## Add derived columns
master.bt_w_cov2 <- ddply(master.bt_w_cov1, .(Stock), transform,
                                RS = R/S,
                                RS_stnd = scale(R/S)[ , 1],
                                lnRS = log(R/S),
                                S_stnd = scale(S)[ , 1],
                                early_sst_stnd = scale(early_sst)[ , 1], 
                                np_pinks_sec_stnd = scale(np_pinks_sec)[ , 1])


## Fill in missing years that fall w/in min and max BY for each stock
master.bt_w_cov3 <- fill.time.series(master.bt_w_cov2)

# Export to output
sock <- master.bt_w_cov3
write.csv(sock, "data/master_brood_table_covar.csv", row.names=FALSE)

