## (2) Clean-up the downloaded sockeye data ##
## ----------------------------------------  ##


## This script cleans/processes the raw compiled sockeye data. The output of
## this script is a master brood table for the analysis and a summary info table

## Read in downloaded data
data_full <- read.csv("./data-downloaded/salmon_productivity_compilation2023-12-28.csv", row.names=1) 
info_full <- read.csv("./data-downloaded/stock_info2023-12-28.csv", row.names=1)
# Data source: https://github.com/Pacific-salmon-assess/dfo_salmon_compilation

## Filter for Sockeye
s.brood <- data_full %>% dplyr::filter(species == "Sockeye")
s.info <- info_full %>% dplyr::filter(species == "Sockeye")


## Clean-up master brood table -----------------------------
head(s.brood)
tail(s.brood)
sapply(s.brood, class)
summary(s.brood)

# NAs in 'RX.X' columns need to be replaced with 0s
# Ragged ends have been dealt with in data compilation
r.cols <- grep("^r[[:digit:]]\\.[[:digit:]]", names(s.brood), value = TRUE)
s.brood[ , r.cols][is.na(s.brood[ , r.cols])] <- 0

# Remove 'RX' columns for Sockeye
r.cols.drop <- grep("^r[[:digit:]]{1}$", names(s.brood))
s.brood <- s.brood[, -r.cols.drop]

head(s.brood)
tail(s.brood)
sapply(s.brood, class)

#Rename some columns
names(s.brood) <- str_to_title(names(s.brood))
s.brood <- dplyr::rename(s.brood, "BY" = "Broodyear", "Stock.ID"="Stock.id", 
                         "R" = "Recruits", "S" = "Spawners")
head(s.brood)


## Lat/Lon lookup table
s.info$lon <- -abs(s.info$lon)
coord.lookup <- distinct(s.info[,c("stock.name", "lat", "lon")])
names(coord.lookup) <- str_to_title(names(coord.lookup))

# join lat/lon
all(s.brood$Stock %in% s.info$stock.name)
all(s.brood$Stock.ID %in% s.info$stock.id)
all(s.brood$Stock %in% coord.lookup$Stock)
s.brood <- left_join(s.brood, coord.lookup, by=c("Stock"="Stock.name"))
head(s.brood)


# Remove species from Stock names
s.brood$Stock <- stringr::str_remove(s.brood$Stock, "-Sockeye")
s.info$stock.name <- stringr::str_remove(s.info$stock.name, "-Sockeye")


## Create data frame with estimates of recruits, spawners and proportion of
## recruits that entered ocean at age 0, 1 and 2 and from each age class
 # Use detailed recruit info if it exists, and average proportions if not
r.cols <- grep("^R[[:digit:]]\\.[[:digit:]]", names(s.brood), value = TRUE)

# add flag indicating if there is recruitment detail (i.e. age structured)
s.brood <- s.brood %>% dplyr::mutate(detailFlag = if_else(rowSums(.[r.cols])==0, 0, 1))

# make Recruits the sum of RX.X cols (if populated)
s.brood.2 <- s.brood
s.brood.2[s.brood.2$detailFlag==1, "R"] <- rowSums(s.brood.2[s.brood.2$detailFlag==1, r.cols]) 
summary(s.brood$R - s.brood.2$R)

# Call function to calculate proportion entering ocean at age 0, 1, ...4
# These are calculated differently for stocks with recruitment detail (detailFlag==1) and without.
bt <- get_ocean_entry_prop(s.brood.2, r.cols=r.cols)

# Make all "R" columns proportions
bt[,r.cols] <- lapply(bt[,r.cols], function(x) x/bt$R)
## Check that all Rx.x columns were included
r1 <- grep("^R[[:digit:]]\\.[[:digit:]]", names(s.brood), value = TRUE)
r2 <- grep("^R[[:digit:]]\\.[[:digit:]]", names(bt), value = TRUE)
all.equal(sort(r1), sort(r2))


# More cleaning steps 
bt.out.1 <- bt[complete.cases(bt),]                # drop years with missing data
(nrow(bt) - nrow(bt.out.1)) # how many rows dropped
bt.out.2 <- subset(bt.out.1, !(Stock %in% c("Osoyoos", "Frazer", "Washington", "Sustut"))) # remove stocks:
    # Frazer -hatchery influence, Oosoyoos -unk ocean entry coords, Washington - ?, Sustut-duplicated
bt.out.3 <- subset(bt.out.2,BY > 1949)        # do this because pre 1950 data is very sparse
bt.out.4 <- subset(bt.out.3, BY<=2015) #- currently have pink-NP data up to 2021 (+6 yr)


# Add ocean regions
bt.out.4$Ocean.Region2 <- s.info$ocean.basin[match(bt.out.4$Stock, s.info$stock.name)]
# Add SEAK grouping
bt.out.4$Ocean.Region2[bt.out.4$Ocean.Region2=="WC" & bt.out.4$Lat >= 54.09] <- "SEAK"

# A step for infilling then re-trimming missing values was removed here. If timeseries gaps are an issue, can be added back.

## Order stocks geographically to make plotting easier
bt.out.5 <- geographic.order(bt.out.4)
bt.out.6 <- dplyr::arrange(bt.out.5, factor(Stock, levels=levels(bt.out.5$Stock)))

# Summary
bt.out <- bt.out.6
head(bt.out)
tail(bt.out)
sapply(bt.out, class)
summary(bt.out)

write.csv(bt.out, "./data/sockeye/master_sockeye_brood_table.csv", row.names = FALSE)



## Create stock info table ---------------------------------
s.info.brood <- ddply(bt.out, .(Stock.ID), summarize,
                            Stock = unique(Stock),
                            Ocean.Region2 = unique(Ocean.Region2),
                            lat = unique(Lat),
                            lon = unique(Lon),
                            na_count = sum(is.na(R)),
                            n_years = sum(!is.na(R)),
                            yr_start = min(BY),
                            yr_end = max(BY))

s.info.brood <- ocean_region_lab(s.info.brood)

## Order stocks geographically to make comparison easier
s.info.brood <- geographic.order(s.info.brood)
s.info.brood <- dplyr::arrange(s.info.brood, factor(Stock, levels=levels(s.info.brood$Stock)))


write.csv(s.info.brood, "./data/sockeye/master_sockeye_stock_info.csv", row.names = FALSE)

sock.info <- s.info.brood

