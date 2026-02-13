## Generate an abridged stock info table for including in publication.

info.all <- bind_rows(sock.info, chum.info, pink.info)
info.all <- info.all |> rename(Region = ocean_region_lab,
                              `First Brood Year` = yr_start,
                              `Last Brood Year` = yr_end,
                              Latitude = lat,
                              Longitude = lon) |>
  select(Species, Stock, Region, `First Brood Year`, `Last Brood Year`,
         Latitude, Longitude)
write.csv(info.all, file=here("data", "all_species_info_pub.csv"), row.names=F)
