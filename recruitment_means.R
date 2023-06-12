##-- "BYPASS" script for estimating age @ ocean entry proportions 
##    for stocks without recruitment detail --##

## Skeena proportions are Table S3 of Price et al (2021) -- "Modern"
## Bristol Bay proportions are average prop. from last 10 years of complete data

library(dplyr) 

##Take Bristol Bay averages

#s.brood.use = brood table filtered for UseFlag=1 and a few other things
bbay_lh <- s.brood.use %>% filter(Region == "Bristol Bay", DetailFlag ==1) #Only keep BBay years with detail

bbay_lh2 <- bbay_lh %>%      #new columns with age@ entry props  
  mutate( prop_age0 = rowSums(select(., starts_with("R0")), na.rm=T)/TotalRecruits, 
         prop_age1 = rowSums(select(., starts_with("R1")), na.rm=T)/TotalRecruits,
         prop_age2 = rowSums(select(., starts_with("R2")), na.rm=T)/TotalRecruits,
         prop_age3 = rowSums(select(., starts_with("R3")), na.rm=T)/TotalRecruits,
         prop_age4 = rowSums(select(., starts_with("R4")), na.rm=T)/TotalRecruits)


bbay_lh3 <- bbay_lh2 %>% 
  dplyr::arrange("Stock.ID", "BY") %>% 
  group_by(Stock) %>% 
  slice_tail(n=10) %>%    # Keep last 10 years of data for each stk
  dplyr::summarize_at(vars(starts_with("prop_")), mean)  # average


## Read in skeena lifehistory
skeena_lh <- read.csv("./data-downloaded/skeena_lifehistory_prop.csv")

## Bind skeena and bristol bay stocks
lifehist <- bind_rows(bbay_lh3, skeena_lh) %>% replace(is.na(.), 0)

## Double check that proportions sum to 1
if(any(round(rowSums(lifehist[,2:6])) != 1)) print("Problem calculating ocean entry proportions")
