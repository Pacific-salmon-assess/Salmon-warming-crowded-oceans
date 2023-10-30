# code to wrangle density measurements and calculate effective density
library(tidyverse)
library(here)

#read competitors, wrangle to long, make more meningful cols. 
competitors <- read.csv(here("data-downloaded/competitor_indices_2023_05_16.csv")) |>
  pivot_longer(-Year, names_to = "type") |>
  mutate(species = ifelse(grepl("pink", type), "pink", "all"),
         metric = ifelse(grepl("numbers", type), "n", "biomass"),
         area = case_when(grepl("np_bs", type) ~ "np-bs",
                          grepl("na_bs", type) ~ "na-bs",
                          grepl("np_goa", type) ~ "np-goa",
                          grepl("na_goa", type) ~ "na-goa",
                          grepl("np_wc", type) ~ "np-wc",
                          grepl("na_wc", type) ~ "na-wc",
                          grepl("np_hatch", type) ~ "na-goa", 
                          type %in% c("pink_numbers_np", "all_spp_biomass_np",
                                      "all_spp_numbers_np") ~ "np",
                          type == "pink_numbers_na" ~ "na")) |>
  select(-type) |>
  filter(!is.na(value))

#break df into biomass only then overwrite w/ effective density
eff_dens <- filter(competitors, metric == "biomass") |>
  mutate(metric = "eff-dens",
         value = ((value)^(1/3))^2) #eff dens calculation

#make it into one object
competitors_long <- bind_rows(competitors, eff_dens)

#write it to repo
write.csv(competitors_long, here("data/competitor_density_long.csv"))