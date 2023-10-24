## Pink & Chum data cleaning

## This script cleans/processes the raw compiled Pink and Chum data. The output of
## this script is a master brood table for the analysis and a summary info table

## PINKS!!

## Read in downloaded data
p.brood <- read.csv("./data-downloaded/pink/raw_pink_brood_table2023-09-20.csv",
                      sep = ",",skip = 0,
                      stringsAsFactors = FALSE, header = TRUE)

p.info <- read.csv("./data-downloaded/pink/pink_info2023-10-19.csv", sep = ",",
                     skip = 0, stringsAsFactors = FALSE, header = TRUE)


## Create lat/lon table ---------------------------------
p.info$lat[p.info$lat==447.16] <- 47.16
p.info$lon <- -abs(p.info$lon)
coord.lookup <- distinct(p.info[,c("stock", "lat", "lon")])
names(coord.lookup) <- str_to_title(names(coord.lookup))


## Clean-up master brood table -----------------------------
head(p.brood)
tail(p.brood)
sapply(p.brood, class)
nrow(p.brood)

## Subset usable data points
p.brood <- p.brood[p.brood$use == 1, ]
head(p.brood)
tail(p.brood)
nrow(p.brood)
summary(p.brood)


bt <- data.frame(Stock.ID = p.brood$stock.id,
                 Species = "Pink",
                 Stock = p.brood$stock,
                 Region = p.brood$region,
                 Ocean.Region2 = NA,
                 Sub.Region = p.brood$sub.region,
                 BY = p.brood$broodyear,
                 R0.1 = p.brood$recruits,
                 R = p.brood$recruits,
                 S = p.brood$spawners )
                
bt.out.1 <- bt[!is.na(bt$S),]  # drop years with missing data
summary(bt.out.1)

bt.out.2 <- subset(bt.out.1,BY <= 2019) # currently have pink-NP data up to 2021 (+2 yr)

## check against stock.info table
all(bt.out.2$Stock %in% p.info$stock)
all(bt.out.2$Stock.ID %in% p.info$stock.id)

# Join Lat/Lon to bt
bt.out.3 <- left_join(bt.out.2, coord.lookup, by="Stock")
head(bt.out.3)

## Add ocean.region groupings
bt.out.4 <- bt.out.3
bt.out.4$Ocean.Region2 <- p.info$ocean.region[match(bt.out.3$Stock, p.info$stock)]
bt.out.4$Ocean.Region2[bt.out.4$Sub.Region %in% c("SEAK", "BC North", "Yakutat")] <- "SEAK"

# In sockeye there is a step here to trim and add NAs for missing years

## Order stocks geographically to make plotting easier
bt.out.5 <- geographic.order(bt.out.4)
bt.out.6 <- dplyr::arrange(bt.out.5, factor(Stock, levels=levels(bt.out.5$Stock)))

# Summary
bt.out <- bt.out.6
head(bt.out)
tail(bt.out)
sapply(bt.out, class)
summary(bt.out)

write.csv(bt.out, "./data/pink/master_pink_brood_table.csv", row.names = FALSE)


## Make stock info table ------------------------------------
p.info.brood <- ddply(bt.out.6, .(Stock.ID), plyr::summarize,
                      Stock = unique(Stock),
                      Region = unique(Region),
                      Ocean.Region2 = unique(Ocean.Region2),
                      Sub.Region = unique(Sub.Region),
                      lat = unique(Lat),
                      lon = unique(Lon),
                      na_count = sum(is.na(R)),
                      n_years = sum(!is.na(R)),
                      yr_start = min(BY),
                      yr_end = max(BY))

p.info.brood <- ocean_region_lab(p.info.brood)
p.info.brood <- geographic.order(p.info.brood) # Order stocks geographically to make comparison easier
p.info.brood <- dplyr::arrange(p.info.brood, factor(Stock, levels=levels(p.info.brood$Stock)))


write.csv(p.info.brood, "./data/pink/master_pink_stock_info.csv", row.names = FALSE)

pink.info <- p.info.brood


## CHUM

## Read in downloaded data
c.brood <- read.csv("./data-downloaded/chum/raw_chum_brood_table2023-10-10.csv",
                      sep = ",",skip = 0,
                      na.string = c("NA", "ND", "<NA>"),
                      stringsAsFactors = FALSE, header = TRUE, row.names=NULL)

c.info <- read.csv("./data-downloaded/chum/chum_info2023-10-10.csv", sep = ",",
                     skip = 0, stringsAsFactors = FALSE, header = TRUE, row.names=NULL)


# Create Lat/Lon lookup table
c.info$lon <- -abs(c.info$lon)
coord.lookup <- distinct(c.info[,c("stock", "lat", "lon")])
names(coord.lookup) <- str_to_title(names(coord.lookup))


## Clean-up master brood table -----------------------------
head(c.brood)
tail(c.brood)
sapply(c.brood, class)

# NA values need to be replaced with 0s for years with no recruits
# This replacement is only done for the "recruit" columns
r.cols <- grep(".\\.[[:digit:]]", names(c.brood), value = TRUE)
c.brood[ , r.cols][is.na(c.brood[ , r.cols])] <- 0
head(c.brood)
tail(c.brood)
sapply(c.brood, class)

#Rename some columns
names(c.brood) <- str_to_title(names(c.brood))
names(c.brood)[names(c.brood) %in% c("Stock.id", "Sub.region", "Broodyear", "Use")] <- c("Stock.ID", "Sub.Region", "BY", "UseFlag")
names(c.brood)[names(c.brood) %in% c("Spawners", "Recruits")] <- c("S", "R")
names(c.brood)[grepl("\\d$", names(c.brood))] <- paste0("R0.", as.numeric(str_extract(names(c.brood)[grepl("\\d$", names(c.brood))], "\\d$"))-1)

## Subset usable data points
c.brood.use <- c.brood[c.brood$UseFlag == 1, ]
c.brood.use <- c.brood.use[ , names(c.brood.use) != "Age"]
c.brood.use <- dplyr::filter(c.brood.use, Stock.ID != 366, Stock.ID != 363) # Remove 2 stocks with no recruitment detail for now
c.brood.use$Species <- "Chum" 

head(c.brood.use)
tail(c.brood.use)
nrow(c.brood.use)
summary(c.brood.use)


r.cols.new <- grep("^R\\d", names(c.brood.use), value = TRUE)
bt = c.brood.use
bt[,r.cols.new] <- lapply(bt[,r.cols.new], function(x) x/bt$R) # Make "R" columns proportions

bt.out.1 <- bt[complete.cases(bt),]                # drop years with missing data
bt.out.2 <- subset(bt.out.1, BY <= 2019) # currently have pink-NP data up to 2021 (+2 yr)

# join lat/lon
all(bt.out.2$Stock %in% c.info$stock)
all(bt.out.2$Stock.ID %in% c.info$stock.id)
all(bt.out.2$Stock %in% coord.lookup$Stock)
bt.out.3 <- left_join(bt.out.2, coord.lookup, by="Stock")
head(bt.out.3)

# Add ocean regions
bt.out.4 <- bt.out.3
bt.out.4$Ocean.Region2 <- c.info$ocean.region[match(bt.out.3$Stock, c.info$stock)]
bt.out.4$Ocean.Region2[bt.out.4$Sub.Region %in% c("SEAK", "BC North", "Yakutat")] <- "SEAK"

# In sockeye there is a step here to trim and add NAs for missing years

## Order stocks geographically to make plotting easier
bt.out.5 <- geographic.order(bt.out.4)
bt.out.6 <- dplyr::arrange(bt.out.5, factor(Stock, levels=levels(bt.out.5$Stock)))

# Summary
bt.out <- bt.out.6
head(bt.out)
tail(bt.out)
sapply(bt.out, class)
summary(bt.out)

write.csv(bt.out, "./data/chum/master_chum_brood_table.csv", row.names = FALSE)


## Create stock info table ---------------------------------
c.info.brood <- ddply(bt.out, .(Stock.ID), summarize,
                      Stock = unique(Stock),
                      Region = unique(Region),
                      Ocean.Region2 = unique(Ocean.Region2),
                      Sub.Region = unique(Sub.Region),
                      lat = unique(Lat),
                      lon = unique(Lon),
                      na_count = sum(is.na(R)),
                      n_years = sum(!is.na(R)),
                      yr_start = min(BY),
                      yr_end = max(BY))

c.info.brood <- ocean_region_lab(c.info.brood)

## Order stocks geographically to make comparison easier
c.info.brood <- geographic.order(c.info.brood)
c.info.brood <- dplyr::arrange(c.info.brood, factor(Stock, levels=levels(c.info.brood$Stock)))


write.csv(c.info.brood, "./data/chum/master_stock_info.csv", row.names = FALSE)

chum.info <- c.info.brood


