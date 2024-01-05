## Source this file to reproduce the project ##
## ----------------------------------- ##

##These scripts are modified from: https://github.com/brendanmichaelconnors/sockeye-climate-competition 

rm(list = ls())

## remove reproducible directories if starting fresh
#unlink("./figures", recursive = TRUE)
#unlink("./output", recursive = TRUE)

time.start <- proc.time()

# 1. Precursors
source ("functions.R")
suppressWarnings(source("load.R")) 
# source("data_download.R")  ## ONLY if climate data have changed - very long run time

# 2. Data cleaning & processing
source("sock_data_clean.R") 
source("pink_chum_data_clean.R")
source("sst_import_process.R") 
source("sst_explore.R") 
source("sock_covariates.R")
source("pink_chum_covariates.R")

# Toggle one species 'ON' for #3-6
#speciesFlag = "pink"
#speciesFlag = "chum"
speciesFlag = "sockeye"

# 3. Exploratory plots and models
source("data_explore.R")  
source("stat_single_stock_fit.R") 

# 4. Stationary Hierarchical Bayesian Models & figures
source("stat_hbm_fit.R") 
source("stat_hbm_inf.R")  

# 5. Dynamic (time-varying) Hierarchical Bayesian Models & figures
source("dyn_hbm_fit.R") 
source("dyn_hbm_inf.R") 

# 6. Hidden Markov Models (optional)
source("hmm_single_stock_fit.R") 
source("hmm_single_stock_inf.R") 


time.run <- proc.time() - time.start
round(time.run[3] / 60, 4)
