## (5) R scripts to create the stock specific covariates used in the analysis ##
## -------------------------------------------------------------------------- ##

# A simplified version of covariates script, including only:
# first year SST, competitors in 2nd year, 

##   Load necessary data

## Master brood table
bt.raw <- read.csv("data/master_brood_table.csv", header=T)
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
master <- merge(bt.complete, early.sst, by=c("BY","Stock.ID"),all.x=T)
master <- merge(master, np.pink.sec, by=c("BY","Stock.ID"),all.x=T)
master.bt_w_cov1 <- master[order(master$Stock.ID),]
head(master.bt_w_cov1)
tail(master.bt_w_cov1)
summary(master.bt_w_cov1)
sapply(master.bt_w_cov1, class)
unique(master.bt_w_cov1$Stock.ID)


## Add derived columns
master.bt_w_cov2 <- plyr::ddply(master.bt_w_cov1, .(Stock.ID), transform,
                                RS = R/S,
                                RS_stnd = scale(R/S)[ , 1],
                                lnRS = log(R/S),
                                S_stnd = scale(S)[ , 1],
                                early_sst_stnd = scale(early_sst)[ , 1], 
                                np_pinks_sec_stnd = scale(np_pinks_sec)[ , 1])


## Fill in missing years that fall w/in min and max BY for each stock
master.bt_w_cov3 <- fill.time.series(master.bt_w_cov2)

## Export to output
master.bt_w_cov <- master.bt_w_cov3
write.csv(master.bt_w_cov, "data/master_brood_table_covar.csv", row.names=FALSE)


sock <- master.bt_w_cov
sock$Stock <- geographic.order(sock) # make Stock an ordered factor for plotting

