## (3) Reads in SST, generates anomalies and then calculates average 
## SST over a specified period and region ## 
## -------------------------------------- ##


## Load necessary data 

sst.raw <- read.csv("./data-downloaded/climate-data/sst_raw.csv")
head(sst.raw)
tail(sst.raw)
sapply(sst.raw, class)
summary(sst.raw)


##   Calculate SST anomalies and average across specified period and region 

## Calculate SST anomalies
sst.anom <- sst.anomaly(sst.raw, ref.years = 1950:2022)
head(sst.anom)
tail(sst.anom)
summary(sst.anom)
sapply(sst.anom, class)

## Convert longitude to match salmon data
sst.anom$lon2 <- ifelse(sst.anom$lon > 180, sst.anom$lon - 360, sst.anom$lon)

write.csv(sst.anom, "data/sst_raw_anomalies.csv", row.names=F)

## Read in population summary table
summary_table <- read.csv("data/master_stock_info.csv", header=T)
head(summary_table)


## Calculate average SST anomaly within area where stock spends first few months of marine life 
sst_yr_1_stock_anomalies <- sst.averager(summary_table, sst.anom, distance = 400)
colnames(sst_yr_1_stock_anomalies) <- c("Year","sst_raw","sst_anomaly","Stock.ID")

write.csv(sst_yr_1_stock_anomalies, "data/sst_yr_1_stock_anomalies.csv", row.names=F)
