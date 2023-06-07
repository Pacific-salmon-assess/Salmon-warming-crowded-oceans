##-- "BYPASS" script for estimating age @ ocean entry proportions 
##    for stocks without recruitment detail --##

## Skeena proportions are Table S3 of Price et al (2021) -- "Modern"
## Bristol Bay proportions are average prop. from last 10 years of complete data

library(dplyr) 

##Take Bristol Bay averages

#s.brood.use = brood table filtered for UseFlag=1 and a few other things
bbay_lh <- filter(s.brood.use, Region == "Bristol Bay") %>% 
           drop_na(num_range("R", 0.1:3.4)) #drop recent years with NA Rx.x cols

bbay_lh2 <- bbay_lh %>% rowwise() %>%       #new columns with age@ entry props  
  mutate(prop_age0 = sum(across(starts_with("R0")), na.rm=T)/TotalRecruits,
         prop_age1 = sum(across(starts_with("R1")), na.rm=T)/TotalRecruits,
         prop_age2 = sum(across(starts_with("R2")), na.rm=T)/TotalRecruits,
         prop_age3 = sum(across(starts_with("R3")), na.rm=T)/TotalRecruits,
         prop_age4 = sum(across(starts_with("R4")), na.rm=T)/TotalRecruits) %>%
  ungroup()     # rowwise() groups by rows, remove these groupings

bbay_lh3 <- bbay_lh2 %>% 
  arrange("Stock.ID", "BY") %>% 
  group_by(Stock) %>% 
  slice_tail(n=10) %>%    # Keep last 10 years of data for each stk
  summarize(across(starts_with("prop_"), ~mean(., na.rm=T)))  # average


## Read in skeena lifehistory
skeena_lh <- read.csv("./data-downloaded/skeena_lifehistory_prop.csv")

## Bind skeena and bristol bay stocks
lifehist <- bind_rows(bbay_lh3, skeena_lh) %>% replace(is.na(.), 0)
