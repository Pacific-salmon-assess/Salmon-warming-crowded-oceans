## Source this file to reproduce the project ##
## ----------------------------------- ##

##These scripts are modified from: https://github.com/brendanmichaelconnors/sockeye-climate-competition 

rm(list = ls())

## remove reproducible directories if starting fresh
#unlink("./figures", recursive = TRUE)
#unlink("./output", recursive = TRUE)

time.start <- proc.time()

source ("functions.R")
suppressWarnings(source("load.R")) 
# source("data_download.R")  ## ONLY if climate data have changed - very long run time
source("sock_data_clean.R") 
source("pink_chum_data_clean.R")
source("sst_import_process.R") 
source("sst_explore.R") 
source("sock_covariates.R")
source("pink_chum_covariates.R")
source("sock_explore.R")  
source("stat_single_stock_fit.R") 
source("stat_hbm_fit.R") 
source("stat_hbm_inf.R")  
source("dyn_hbm_fit.R") 
source("dyn_hbm_inf.R") 
source("hmm_single_stock_fit.R") 
source("hmm_single_stock_inf.R") 


time.run <- proc.time() - time.start
round(time.run[3] / 60, 4)
