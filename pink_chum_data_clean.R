## Pink & Chum data cleaning

## This script cleans/processes the raw compiled Pink and Chum data
## The output is a master brood table for the analysis and a summary info table


## Should R be sum of RX.X columns?!

data_full <- read.csv("./data-downloaded/salmon_productivity_compilation2023-12-28.csv") 

info_full <- read.csv("./data-downloaded/stock_info2023-12-28.csv")

## PINKS!!

## Read in downloaded data
p.brood <- data_full %>% dplyr::filter(species %in% c("Pink-Odd", "Pink-Even"))
p.info <- info_full %>% dplyr::filter(species %in% c("Pink-Odd", "Pink-Even"))

## Create lat/lon table ---------------------------------
p.info$lat[p.info$lat==447.16] <- 47.16 ## report this issue
p.info$lon <- -abs(p.info$lon)
coord.lookup <- distinct(p.info[,c("stock.name", "lat", "lon")])
names(coord.lookup) <- str_to_title(names(coord.lookup))


## Clean-up master brood table -----------------------------
head(p.brood)
tail(p.brood)
sapply(p.brood, class)
nrow(p.brood)


bt <- data.frame(Stock.ID = p.brood$stock.id,
                 Species = "Pink",
                 Stock = p.brood$stock,
                 Ocean.Region2 = NA,
                 BY = p.brood$broodyear,
                 R0.1 = p.brood$recruits,
                 R = p.brood$recruits,
                 S = p.brood$spawners )
                
bt.out.1 <- bt[!is.na(bt$S),]  # drop years with missing data
summary(bt.out.1)

bt.out.2 <- subset(bt.out.1, BY <= 2020) # currently have pink-NP data up to 2021 (+1 yr)

## check against stock.info table
all(bt.out.2$Stock %in% p.info$stock.name)
all(bt.out.2$Stock.ID %in% p.info$stock.id)

# Join Lat/Lon to bt
bt.out.3 <- left_join(bt.out.2, coord.lookup, by=c("Stock"="Stock.name"))
head(bt.out.3)

## Add ocean.region groupings
bt.out.4 <- bt.out.3
bt.out.4$Ocean.Region2 <- p.info$ocean.basin[match(bt.out.3$Stock, p.info$stock.name)]
# Add SEAK grouping
bt.out.4$Ocean.Region2[bt.out.4$Ocean.Region2=="WC" & bt.out.4$Lat >= 54.09] <- "SEAK"
#bt.out.4$Ocean.Region2[bt.out.4$Sub.Region %in% c("SEAK", "BC North", "Yakutat")] <- "SEAK"

# In sockeye there is a step here to trim and add NAs for missing years

## Order stocks geographically to make plotting easier
bt.out.5 <- geographic.order(bt.out.4)
bt.out.6 <- dplyr::arrange(bt.out.5, factor(Stock, levels=levels(bt.out.5$Stock)))

# Simplify stock names
bt.out.7 <- bt.out.6
#bt.out.7$Even.Odd <- stringr::str_extract(bt.out.7$Stock, "\\w{3,4}$")
bt.out.7$Stock <- stringr::str_remove(bt.out.7$Stock, ".\\(\\w{3,4}\\)")
bt.out.7$Stock <- stringr::str_remove(bt.out.7$Stock, "-Pink")

# Drop stocks (temporarily?!):
bt.out.7[bt.out.7$R==0, "R"] # check for zero recruits
bt.out.7[bt.out.7$S==0, "S"] # check for zero spawners
bt.out.7 <- bt.out.7 %>% filter(!Stock.ID %in% c(148, 149)) # remove Nisqually and S. South Misc temporarily
bt.out.7 <- bt.out.7 %>% filter(!Stock.ID %in% c(202, 203, 214, 204, 209)) # remove Haida Gwaii stocks & Homathko
# Summary
bt.out <- bt.out.7
head(bt.out)
tail(bt.out)
sapply(bt.out, class)
summary(bt.out)

write.csv(bt.out, "./data/pink/master_pink_brood_table.csv", row.names = FALSE)


## Make stock info table ------------------------------------
p.info.brood <- ddply(bt.out.7, .(Stock.ID), plyr::summarize,
                      Stock = unique(Stock),
                      Ocean.Region2 = unique(Ocean.Region2),
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


## CHUM -----------------------------------------------

## Read in downloaded data
c.brood <- data_full %>% dplyr::filter(species == "Chum")
c.info <- info_full %>% dplyr::filter(species  == "Chum")

# Create Lat/Lon lookup table
c.info$lon <- -abs(c.info$lon)
coord.lookup <- distinct(c.info[,c("stock.name", "lat", "lon")])
names(coord.lookup) <- str_to_title(names(coord.lookup))


## Clean-up master brood table -----------------------------
c.brood <- c.brood[,names(c.brood)!="X"]
head(c.brood)
tail(c.brood)
sapply(c.brood, class)


# Remove empty "RX.X columns and rename "RX" columns to Eur naming convention
names(c.brood)[grepl("\\w{1}\\d{1}$", names(c.brood))] <- paste0("R0.", as.numeric(str_extract(names(c.brood)[grepl("\\w{1}\\d{1}$", names(c.brood))], "\\d$"))-1)
r.cols.old <- names(c.brood)[grepl("^r\\d", names(c.brood))]
c.brood <- c.brood[,!names(c.brood) %in% r.cols.old]

#Rename some columns
names(c.brood) <- str_to_title(names(c.brood))
names(c.brood)[names(c.brood) %in% c("Broodyear", "Stock.id")] <- c("BY", "Stock.ID")
names(c.brood)[names(c.brood) %in% c("Spawners", "Recruits")] <- c("S", "R")


head(c.brood)
tail(c.brood)
nrow(c.brood)
summary(c.brood)


bt <- c.brood
r.cols <- names(bt)[grepl("R0.\\d$", names(bt))]
# Replace NAs with 0 in recruitment cols
bt[ , r.cols][is.na(c.brood[,r.cols])] <- 0
rowSums(bt[,r.cols]) == bt$R
bt[,r.cols] <- lapply(bt[,r.cols], function(x) x/bt$R) # Make "R" columns proportions

bt.out.1 <- bt[complete.cases(bt),]                # drop years with missing data
bt.out.2 <- subset(bt.out.1, BY <= 2019) # currently have pink-NP data up to 2021 (+2 yr)

head(bt.out.2)
summary(bt.out.2)

# join lat/lon
all(bt.out.2$Stock %in% c.info$stock.name)
all(bt.out.2$Stock.ID %in% c.info$stock.id)
all(bt.out.2$Stock %in% coord.lookup$Stock)
bt.out.3 <- left_join(bt.out.2, coord.lookup, by=c("Stock"="Stock.name"))
head(bt.out.3)

# Add ocean regions
bt.out.4 <- bt.out.3
bt.out.4$Ocean.Region2 <- c.info$ocean.basin[match(bt.out.4$Stock, c.info$stock.name)]
# Add SEAK grouping
bt.out.4$Ocean.Region2[bt.out.4$Ocean.Region2=="WC" & bt.out.4$Lat >= 54.09] <- "SEAK"


# In sockeye there is a step here to trim and add NAs for missing years

## Order stocks geographically to make plotting easier
bt.out.5 <- geographic.order(bt.out.4)
bt.out.6 <- dplyr::arrange(bt.out.5, factor(Stock, levels=levels(bt.out.5$Stock)))

# Simplify stock names
bt.out.7 <- bt.out.6
bt.out.7$Stock <- stringr::str_remove(bt.out.7$Stock, "-Chum")


# Summary
bt.out <- bt.out.7
head(bt.out)
tail(bt.out)
sapply(bt.out, class)
summary(bt.out)

write.csv(bt.out, "./data/chum/master_chum_brood_table.csv", row.names = FALSE)


## Create stock info table ---------------------------------
c.info.brood <- ddply(bt.out, .(Stock.ID), summarize,
                      Stock = unique(Stock),
                      Ocean.Region2 = unique(Ocean.Region2),
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


