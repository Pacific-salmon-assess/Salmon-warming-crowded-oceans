## Exploratory things - HIDDEN ##

## -- Quick look at filtering/infilling in skeena dataset --## 
library(dplyr)

skeena <- read.csv("C:/Users/hunterh/Documents/sockeye-climate-competition/nonstationary_dynamics/data/raw data/sockeye/Skeena_Nass_sockeye.csv")

main <- filter(skeena, Label=="Main")
infill <- filter(skeena, Label=="Main_Infill")
filter45 <- filter(skeena, Label=="Filter45")

# See what was infilled
summary(filter(infill, Infill==TRUE)) # 36 total infills
unique(infill$Stock[which(infill$Infill==TRUE)]) # across 13 stocks #all but four are stocks we removed from analysis due to ts length


# See what was filtered
summary(filter(filter45, Filter==TRUE)) # 7 obs filtered out
unique(filter45$Stock[which(filter45$Filter==TRUE)]) # 5 stocks #all removed from analysis due to ts length
