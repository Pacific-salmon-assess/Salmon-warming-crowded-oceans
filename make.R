## Source this file to reproduce the (simplified) project ##
## Fits only one model to 2 covariates ##
## ----------------------------------- ##

##These scripts are modified from: https://github.com/brendanmichaelconnors/sockeye-climate-competition 

rm(list = ls())

## remove reproducible directories
unlink("./figures", recursive = TRUE)
unlink("./output", recursive = TRUE)
unlink("./pub", recursive = TRUE)


time.start <- proc.time()

#source("data_download.R") # ONLY if data have changed
source ("functions.R")
source("load.R") #(1)
# source("data_download.R")  ## only needs to be run if data have changed
source("sock_data_clean.R") #(2)
source("sst_import_process.R") #(3)
source("sst_explore.R") #(4)
source("sock_covariates.R") #(5)
source("sock_explore.R") #(6)
source("single_stock.R") # (7)
source("hbm_fit.R") #(8)
source("pub.R") #(9)
source("dynamic_hbm_fit.R") #(10)
source("dynamic_hbm_inf.R") #(11)
source("hmm_test_fit.R") # (12)


time.run <- proc.time() - time.start
round(time.run[3] / 60, 4)
