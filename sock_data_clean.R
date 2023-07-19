## (2) Clean-up the downloaded sockeye data ##
## ----------------------------------------  ##


## This script cleans/processes the raw downloaded sockeye data. The output of
## this script is a master brood table for the anlaysis and a summary info table


## Read in downloaded data
s.brood <- read.table("./data-downloaded/raw_brood_table_2023_06_18.csv",
                      sep = ",",skip = 0,
                      na.string = c("NA", "ND", "Canadians:  Please check the G-R ages.  "),
                      stringsAsFactors = FALSE, header = TRUE)

s.info <- read.table("./data-downloaded/raw_stock_info_2023_06_18.csv", sep = ",",
                     skip = 0, stringsAsFactors = FALSE, header = TRUE, quote="", fill=TRUE)


## Clean-up master brood table -----------------------------
head(s.brood)
tail(s.brood)
sapply(s.brood, class)

# NA values need to be replaced with 0s for years with no recruits
# This replacement is only done for the "recruit" columns
r.cols <- grep("^R[[:digit:]]\\.[[:digit:]]", names(s.brood), value = TRUE)
s.brood[ , r.cols][is.na(s.brood[ , r.cols])] <- 0
head(s.brood)
tail(s.brood)
sapply(s.brood, class)

## Don't use Nushagak prior to 1985
s.brood[s.brood$Stock == "Nushagak", ]
s.brood$UseFlag[s.brood$Stock == "Nushagak" & s.brood$BroodYear < 1985] <- 0

## Don't use L. Washington
s.brood$UseFlag[s.brood$Stock == "Washington"] <- 0

## Subset usable data points
s.brood.use <- s.brood[s.brood$UseFlag == 1, ]
head(s.brood.use)
tail(s.brood.use)
nrow(s.brood.use)
summary(s.brood.use)

names(s.brood.use)[names(s.brood.use) == "BroodYear"] <- "BY"
if(!"R2.5" %in% names(s.brood)) {
    s.brood.use$R2.5 <- 0
}


## Create data frame with estimates of recruits, spawners and proportion of
## recruits that entered ocean at age 0, 1 and 2 and from each age class
 # Use detailed recruit info if it exists, and average proportions if not
source("recruitment_means.R")

bt <- ddply(s.brood.use, c("Stock.ID", "BY"),function(x) {
  cond <- x$DetailFlag==1
	R <- if_else(cond, sum(x[ , grep("^R[[:digit:]]\\.[[:digit:]]", names(x))], na.rm = TRUE), 
	             x$TotalRecruits) # keep total recruits as is if no detail, otherwise sum recruits
	S <- x$TotalEscapement # spawners
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
	DetailFlag <- x$DetailFlag # keep this flag in for now
	R0.1 <- x[,"R0.1"]/R
	R0.2 <- x[,"R0.2"]/R
	R0.3 <- x[,"R0.3"]/R
	R0.4 <- x[,"R0.4"]/R
	R0.5 <- x[,"R0.5"]/R
	R1.1 <- x[,"R1.1"]/R
	R1.2 <- x[,"R1.2"]/R
	R1.3 <- x[,"R1.3"]/R
	R1.4 <- x[,"R1.4"]/R
	R1.5 <- x[,"R1.5"]/R
	R2.1 <- x[,"R2.1"]/R
	R2.2 <- x[,"R2.2"]/R
	R2.3 <- x[,"R2.3"]/R
	R2.4 <- x[,"R2.4"]/R
	R2.5 <- x[,"R2.5"]/R
	R3.1 <- x[,"R3.1"]/R
	R3.2 <- x[,"R3.2"]/R
	R3.3 <- x[,"R3.3"]/R
	R3.4 <- x[,"R3.4"]/R
	R4.1 <- x[,"R4.1"]/R
	R4.2 <- x[,"R4.2"]/R
	R4.3 <- x[,"R4.3"]/R
    data.frame(Stock = x$Stock,
               Region = x$Region,
               Ocean.Region = x$Ocean.Region,
               Sub.Region = x$Sub.Region,
               Lat = x$Lat,
               Lon = x$Lon,
               R, S, ocean_0, ocean_1, ocean_2, ocean_3, ocean_4, DetailFlag, R0.1, R0.2, R0.3, R0.4,
               R0.5, R1.1, R1.2, R1.3, R1.4, R1.5, R2.1, R2.2, R2.3, R2.4, R2.5, R3.1,
               R3.2, R3.3, R3.4, R4.1, R4.2, R4.3, stringsAsFactors = FALSE)
	})


## Check that all Rx.x columns were included
r1 <- grep("^R[[:digit:]]\\.[[:digit:]]", names(s.brood.use), value = TRUE)
r2 <- grep("^R[[:digit:]]\\.[[:digit:]]", names(bt), value = TRUE)
all.equal(sort(r1), sort(r2))

bt.out.1 <- bt[complete.cases(bt),]                # drop years with missing data
bt.out.2 <- subset(bt.out.1, Stock != "Osoyoos") # drop Osoyooos for now, do not have ocean entry coordinates
bt.out.2 <- subset(bt.out.2, Stock != "Frazer") # drop Frazer, hatchery influence
bt.out.3 <- subset(bt.out.2,BY <= 2015) # currently have pink-NP data up to 2021 (+6 yr)
bt.out.4 <- subset(bt.out.3,BY > 1949)        # do this because pre 1950 data is very sparse

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


bt.out <- bt.out.6
head(bt.out)
tail(bt.out)
sapply(bt.out, class)
summary(bt.out)

## Add alternate ocean.region grouping
bt.out$Ocean.Region2 <- bt.out$Ocean.Region
bt.out$Ocean.Region2 <- ifelse(bt.out$Region %in% c("SEAK", "BC North"), "SEAK", bt.out$Ocean.Region2)


write.csv(bt.out, "./data/master_brood_table.csv", row.names = FALSE)



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
s.info.brood$ocean_label <- case_when(s.info.brood$Ocean.Region == "WC" ~ "West Coast",
                                      s.info.brood$Ocean.Region == "GOA" ~ "Gulf of Alaska",
                                      s.info.brood$Ocean.Region == "BS" ~ "Bering Sea")
s.info.brood$ocean_label2 <- case_when(s.info.brood$Ocean.Region2 == "WC" ~ "West Coast", 
                                       s.info.brood$Ocean.Region2 == "GOA" ~ "Gulf of Alaska", 
                                       s.info.brood$Ocean.Region2 == "BS" ~ "Bering Sea",
                                       s.info.brood$Ocean.Region2 == "SEAK" ~ "Southeast Alaska")

write.csv(s.info.brood, "./data/master_stock_info.csv", row.names = FALSE)

sock.info <- s.info.brood
sock.info$Stock <- factor(sock.info$Stock, levels = unique(sock.info$Stock))

