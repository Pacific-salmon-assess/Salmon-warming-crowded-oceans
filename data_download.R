## Download data for project ## 
## ONLY if data have changed) ##
## -------------------------- ##

## Locations of data:
## Master brood table: "./data-downloaded/raw_brood_table_2023_05_19"
## Pink abundance: "./data-downloaded/competitor_indices_2023_05_16.csv"


## Climate indices download

## SST raw -------------------------------------------------
ersst::sst_download(years = 1950:2022,
                    months = 1:12,
                    save.dir = "./data-downloaded/climate-data/sst_raw/",
                    version = 5)

sst.raw.full <- ersst::sst_load(years = 1950:2022,
                                months = 1:12,
                                read.dir = "./data-downloaded/climate-data/sst_raw/",
                                version = 5)

sst.raw.np <- ersst::sst_subset_space(sst.raw.full,
                                      lat.min = 36,
                                      lat.max = 80,
                                      lon.min = 170,
                                      lon.max = 250)

sst.raw.df <- ersst::sst_dataframe(sst.raw.np)

write.csv(sst.raw.df, "./data-downloaded/climate-data/sst_raw.csv", row.names = FALSE)



## PDO + NPGO ----------------------------------------------
years <- 1950:2022
pdo <- rsoi::download_pdo(years) 
npgo <- rsoi::download_npgo(years)

write.csv(pdo, "./data-downloaded/climate-data/pdo.csv", row.names = FALSE)
write.csv(npgo, "./data-downloaded/climate-data/npgo.csv", row.names = FALSE)


## SST raw extended (1855-present) -------------------------------------------------
ersst::sst_download(years = 1855:2022,
                    months = 1:12,
                    save.dir = "./data-downloaded/climate-data/sst_raw_extend/",
                    version = 5)

sst.raw.full <- ersst::sst_load(years = 1855:2022,
                                months = 1:12,
                                read.dir = "./data-downloaded/climate-data/sst_raw_extend/",
                                version = 5)

sst.raw.np <- ersst::sst_subset_space(sst.raw.full,
                                      lat.min = 36,
                                      lat.max = 80,
                                      lon.min = 170,
                                      lon.max = 250)

sst.raw.df <- ersst::sst_dataframe(sst.raw.np)

write.csv(sst.raw.df, "./data-downloaded/climate-data/sst_raw_extend.csv", row.names = FALSE)


