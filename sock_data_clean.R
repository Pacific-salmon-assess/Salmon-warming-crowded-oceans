## (2) Clean-up the downloaded sockeye data ##
## ----------------------------------------  ##


## This script cleans/processes the raw compiled sockeye data. The output of
## this script is a master brood table for the analysis and a summary info table

## Read in downloaded data
data_full <- read.csv("./data-downloaded/salmon_productivity_compilation2023-11-20.csv", row.names=1) 
info_full <- read.csv("./data-downloaded/stock_info2023-11-20.csv", row.names=1)

## Filter for Sockeye
s.brood <- data_full %>% dplyr::filter(species == "Sockeye")
s.info <- info_full %>% dplyr::filter(species == "Sockeye")


## Clean-up master brood table -----------------------------
head(s.brood)
tail(s.brood)
sapply(s.brood, class)
summary(s.brood)

# NA values need to be replaced with 0s for years with no recruits
# This replacement is only done for the "recruit" columns
r.cols.old <- grep("^r[[:digit:]]\\.[[:digit:]]", names(s.brood), value = TRUE)
s.brood[ , r.cols.old][is.na(s.brood[ , r.cols.old])] <- 0
if(!"R2.5" %in% names(s.brood)) {
  s.brood$r2.5 <- 0 }
head(s.brood)
tail(s.brood)
sapply(s.brood, class)

#Rename some columns
names(s.brood) <- str_to_title(names(s.brood))
names(s.brood)[names(s.brood) %in% c("Broodyear", "Stock.id")] <- c("BY", "Stock.ID")
#names(s.brood)[names(s.brood) %in% c("Spawners", "Recruits")] <- c("S", "R")
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


# Shorten Stock names
s.brood$Stock <- stringr::str_remove(s.brood$Stock, "-Sockeye")
# temporarily remove stock with duplicated BYs:
s.brood <- s.brood %>% filter(Stock!="Sustut")


## Create data frame with estimates of recruits, spawners and proportion of
## recruits that entered ocean at age 0, 1 and 2 and from each age class
 # Use detailed recruit info if it exists, and average proportions if not
r.cols <- grep("^R[[:digit:]]\\.[[:digit:]]", names(s.brood), value = TRUE)

source("recruitment_means.R")

bt <- ddply(s.brood, c("Stock.ID", "BY"),function(x) {
  DetailFlag <- ifelse(all(x[r.cols]==0), 0, 1)
  if(nrow(x) > 1) print("Duplicates present!", x) 
  cond <- DetailFlag == 1
	R <- if_else(cond, sum(x[ , grep("^R[[:digit:]]\\.[[:digit:]]", names(x))], na.rm = TRUE), 
	             x$Recruits) # keep total recruits as is if no detail, otherwise sum recruits
	S <- x$Spawners # spawners
	ocean_0 <- ifelse(cond, sum(x[,grep("^R0\\.", names(x))],na.rm=T)/R, 
	                  lifehist$prop_age0[match(x$Stock, lifehist$Stock)])
	ocean_1 <- ifelse(cond, sum(x[,grep("^R1\\.", names(x))],na.rm=T)/R, 
	                  lifehist$prop_age1[match(x$Stock, lifehist$Stock)])
	ocean_2 <- ifelse(cond, sum(x[,grep("^R2\\.", names(x))],na.rm=T)/R, 
	                  lifehist$prop_age2[match(x$Stock, lifehist$Stock)])
	ocean_3 <- ifelse(cond, sum(x[,grep("^R3\\.", names(x))],na.rm=T)/R, 
	                  lifehist$prop_age3[match(x$Stock, lifehist$Stock)])
	ocean_4 <- ifelse(cond, sum(x[,grep("^R4\\.", names(x))],na.rm=T)/R, 
	                  lifehist$prop_age4[match(x$Stock, lifehist$Stock)])
    data.frame(Species="Sockeye",
               Stock = x$Stock,
               Lat = x$Lat,
               Lon = x$Lon,
               R, S, ocean_0, ocean_1, ocean_2, ocean_3, ocean_4, x$R0.1, x$R0.2, x$R0.3, x$R0.4, x$R0.5, x$R1.1, x$R1.2, x$R1.3, x$R1.4, x$R1.5, x$R2.1, x$R2.2, x$R2.3, x$R2.4, x$R2.5, x$R3.1, x$R3.2, x$R3.3, x$R3.4, x$R4.1, x$R4.2, x$R4.3, stringsAsFactors = FALSE)
	}) 

## Need to figure out what to do with new sockeye data in the above

bt[,r.cols] <- lapply(bt[,r.cols], function(x) x/bt$R) # Make "R" columns proportions



## Check that all Rx.x columns were included
r1 <- grep("^R[[:digit:]]\\.[[:digit:]]", names(s.brood.use), value = TRUE)
r2 <- grep("^R[[:digit:]]\\.[[:digit:]]", names(bt), value = TRUE)
all.equal(sort(r1), sort(r2))

bt.out.1 <- bt[complete.cases(bt),]                # drop years with missing data
bt.out.2 <- subset(bt.out.1, !(Stock %in% c("Osoyoos", "Frazer", "Washington"))) # post hoc removal of stocks:
    # Frazer -hatchery influence, Oosoyoos -unk ocean entry coords, Washington - ?
bt.out.2 <- subset(bt.out.2, !(Stock == "Nushagak" & BY < 1985)) # Nushagak unreliable <1985
bt.out.3 <- subset(bt.out.2,BY <= 2015) # currently have pink-NP data up to 2021 (+6 yr)
bt.out.4 <- subset(bt.out.3,BY > 1949)        # do this because pre 1950 data is very sparse

## Add alternate (now MAIN) ocean.region grouping
bt.out.4$Ocean.Region2 <- bt.out.4$Ocean.Region
bt.out.4$Ocean.Region2 <- ifelse(bt.out.4$Region %in% c("SEAK", "BC North"), "SEAK", bt.out.4$Ocean.Region2)

## Fill in missing years that fall w/in min and max BY for each stock
bt.out.5 <- fill.time.series(bt.out.4) # this adds NAs and is supposed to! 


## Trim time series of NA values by selecting the longest
## string of consecutive brood years for stocks with NA values.
bt.out.6 <- ddply(bt.out.5, .(Stock.ID), function(i) {
                        BY <- i$BY
                        R  <- i$R
                        BY.na <- ifelse(is.na(R), NA, BY)
                        ind <- consec_years(BY.na)
                        return(i[i$BY %in% ind, ])
                       })

## Order stocks geographically to make plotting easier
bt.out.6 <- geographic.order(bt.out.6)
bt.out.7 <- dplyr::arrange(bt.out.6, factor(Stock, levels=levels(bt.out.6$Stock)))

# Summary
bt.out <- bt.out.7
head(bt.out)
tail(bt.out)
sapply(bt.out, class)
summary(bt.out)

write.csv(bt.out, "./data/sockeye/master_brood_table.csv", row.names = FALSE)



## Create stock info table ---------------------------------
s.info.brood <- ddply(bt.out, .(Stock.ID), summarize,
                            Stock = unique(Stock),
                            Region = unique(Region),
                            Ocean.Region = unique(Ocean.Region),
                            Ocean.Region2 = unique(Ocean.Region2),
                            Sub.Region = unique(Sub.Region),
                            lat = unique(Lat),
                            lon = unique(Lon),
                            na_count = sum(is.na(R)),
                            n_years = sum(!is.na(R)),
                            yr_start = min(BY),
                            yr_end = max(BY))

s.info.brood$alpha_region <- ifelse(s.info.brood$Region == "Fraser River", "FR", "AK")
s.info.brood <- ocean_region_lab(s.info.brood)

## Order stocks geographically to make comparison easier
s.info.brood <- geographic.order(s.info.brood)
s.info.brood <- dplyr::arrange(s.info.brood, factor(Stock, levels=levels(s.info.brood$Stock)))


write.csv(s.info.brood, "./data/sockeye/master_stock_info.csv", row.names = FALSE)

sock.info <- s.info.brood

