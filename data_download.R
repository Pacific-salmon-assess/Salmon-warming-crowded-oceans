## Download data for project ## 
## ONLY if data have changed) ##
## -------------------------- ##

## Locations of data:
## Master brood table: "./data-downloaded/raw_brood_table_2023_05_19"
## Pink abundance: "./data/competitor_indices_2023_05_16.csv"


## Climate indices download

## SST raw -------------------------------------------------
ersst::sst_download(years = 1950:2019,
                    months = 1:12,
                    save.dir = "./data-downloaded/sst_raw/",
                    version = 5)

sst.raw.full <- ersst::sst_load(years = 1950:2019,
                                months = 1:12,
                                read.dir = "./data-downloaded/sst_raw/",
                                version = 5)

sst.raw.np <- ersst::sst_subset_space(sst.raw.full,
                                      lat.min = 36,
                                      lat.max = 80,
                                      lon.min = 170,
                                      lon.max = 250)

sst.raw.df <- ersst::sst_dataframe(sst.raw.np)

write.csv(sst.raw.df, "./data/sst_raw.csv", row.names = FALSE)



## PDO + NPGO ----------------------------------------------
years <- 1950:2019
pdo <- rsoi::download_pdo(years) 
npgo <- rsoi::download_npgo(years)

write.csv(pdo, "./data/pdo.csv", row.names = FALSE)
write.csv(npgo, "./data/npgo.csv", row.names = FALSE)

