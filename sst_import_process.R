## (3) Reads in SST, generates anomalies and then calculates average 
## SST over a specified period and region ## 
## -------------------------------------- ##


### Regular time period: 1950-2022 ------------------------------------

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

# (1) Sockeye
## Read in population summary table 
sock_summary_table <- read.csv("data/sockeye/master_sockeye_stock_info.csv", header=T)
head(sock_summary_table)

## Calculate average SST anomaly within area where stock spends first few months of marine life 
sock_sst_yr_1_stock_anomalies <- sst.averager(sock_summary_table, sst.anom, distance = 400)
colnames(sock_sst_yr_1_stock_anomalies) <- c("Year","sst_raw","sst_anomaly","Stock.ID")
sock_sst_yr_1_stock_anomalies$Species <- "Sockeye"

# (2) Pink
## Read in population summary table 
pink_summary_table <- read.csv("data/pink/master_pink_stock_info.csv", header=T)
head(pink_summary_table)

## Calculate average SST anomaly within area where stock spends first few months of marine life 
pink_sst_yr_1_stock_anomalies <- sst.averager(pink_summary_table, sst.anom, distance = 400)
colnames(pink_sst_yr_1_stock_anomalies) <- c("Year","sst_raw","sst_anomaly","Stock.ID")
pink_sst_yr_1_stock_anomalies$Species <- "Pink"

# (3) Chum
chum_summary_table <- read.csv("data/chum/master_chum_stock_info.csv", header=T)
head(chum_summary_table)

## Calculate average SST anomaly within area where stock spends first few months of marine life 
chum_sst_yr_1_stock_anomalies <- sst.averager(chum_summary_table, sst.anom, distance = 400)
colnames(chum_sst_yr_1_stock_anomalies) <- c("Year","sst_raw","sst_anomaly","Stock.ID")
chum_sst_yr_1_stock_anomalies$Species <- "Chum"

# bind them
sst_yr_1_stock_anomalies <- rbind(sock_sst_yr_1_stock_anomalies, pink_sst_yr_1_stock_anomalies, chum_sst_yr_1_stock_anomalies)

write.csv(sst_yr_1_stock_anomalies, "./data/sst_yr_1_stock_anomalies.csv", row.names=F)


### Extended time period: 1855-2022 ------------------------------------

## Load data for extended SST time series

sst.raw <- read.csv("./data-downloaded/climate-data/sst_raw_extend.csv")
head(sst.raw)
tail(sst.raw)
sapply(sst.raw, class)
summary(sst.raw)


##   Calculate SST anomalies and average across specified period and region 

## Calculate SST anomalies
sst.anom <- sst.anomaly(sst.raw, ref.years = 1855:2022)
head(sst.anom)
tail(sst.anom)
summary(sst.anom)
sapply(sst.anom, class)

## Convert longitude to match salmon data
sst.anom$lon2 <- ifelse(sst.anom$lon > 180, sst.anom$lon - 360, sst.anom$lon)

write.csv(sst.anom, "data/sst_raw_anomalies_extend.csv", row.names=F)

# (1) Sockeye
## Read in population summary table 
sock_summary_table <- read.csv("data/sockeye/master_sockeye_stock_info.csv", header=T)
head(sock_summary_table)

## Calculate average SST anomaly within area where stock spends first few months of marine life 
sock_sst_yr_1_stock_anomalies <- sst.averager(sock_summary_table, sst.anom, distance = 400)
colnames(sock_sst_yr_1_stock_anomalies) <- c("Year","sst_raw","sst_anomaly","Stock.ID")
sock_sst_yr_1_stock_anomalies$Species <- "Sockeye"

# (2) Pink
## Read in population summary table 
pink_summary_table <- read.csv("data/pink/master_pink_stock_info.csv", header=T)
head(pink_summary_table)

## Calculate average SST anomaly within area where stock spends first few months of marine life 
pink_sst_yr_1_stock_anomalies <- sst.averager(pink_summary_table, sst.anom, distance = 400)
colnames(pink_sst_yr_1_stock_anomalies) <- c("Year","sst_raw","sst_anomaly","Stock.ID")
pink_sst_yr_1_stock_anomalies$Species <- "Pink"

# (3) Chum
chum_summary_table <- read.csv("data/chum/master_chum_stock_info.csv", header=T)
head(chum_summary_table)

## Calculate average SST anomaly within area where stock spends first few months of marine life 
chum_sst_yr_1_stock_anomalies <- sst.averager(chum_summary_table, sst.anom, distance = 400)
colnames(chum_sst_yr_1_stock_anomalies) <- c("Year","sst_raw","sst_anomaly","Stock.ID")
chum_sst_yr_1_stock_anomalies$Species <- "Chum"

# bind them
sst_yr_1_stock_anomalies_extend <- rbind(sock_sst_yr_1_stock_anomalies, pink_sst_yr_1_stock_anomalies, chum_sst_yr_1_stock_anomalies)

write.csv(sst_yr_1_stock_anomalies_extend, "./data/sst_yr_1_stock_anomalies_extend.csv", row.names=F)
