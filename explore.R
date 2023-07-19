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


## -- Some unused plots from early HMM exploration (single-stock) --##

# Covariate effects for each state by stock
## Density plots of beta1 estimates # not that useful
hmm_tidy %>% ggplot(aes(x=beta1, fill=state_K)) + 
  geom_histogram() + theme_minimal()
# Dot plot of the covar effects by stock
asc.beta <- hmm_tidy$stock[match(sort(unique(hmm_tidy$beta1[which(hmm_tidy$state_K == "State 1")])), hmm_tidy$beta1)]
hmm_tidy %>% ggplot(aes(x=beta1, y=factor(stock, levels=asc.beta), colour=state_K)) +
  geom_point() # this actually hurts my eyes
