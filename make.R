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
#source("sst_import_process.R")
#source("sst_explore.R")
source("sock_covariates.R")
source("pink_chum_covariates.R")

# Colour scheme for all plots
col.region <- rev(chroma::qpal(7, luminance = 40)[c(1, 3, 5, 7)])
names(col.region) <- c("West Coast", "Southeast Alaska", "Gulf of Alaska", "Bering Sea")
col.scale.reg <- scale_colour_manual(name = "Ocean Region", values=col.region)

# Steps 3-6 run one species at a time
species = c("sockeye", "pink", "chum")

for(sp in species) {

  speciesFlag = sp # set species

  # 3. Exploratory plots and models
  source("data_explore.R")
  source("stat_single_stock_fit.R")

  # 4. Stationary Hierarchical Bayesian Models & figures
  source("stat_hbm_fit.R")
  source("stat_hbm_inf.R")
  source("interaction-fits.R") # Run stationary models with interaction term

  # 5. Dynamic (time-varying) Hierarchical Bayesian Models & figures
  source("era_hbm_fit.R")
  source("era_hbm_inf.R")
  source("rw_model_fit.R") # run NEW random walk models (have done for all spp, but not below)

  # 6. Sensitivity analyses
  source("sensitivity-analyses/comp-sens.R") # Alternative competitor indices; long run time
  source("sensitivity-analyses/pdo-npgo-sens.R") # NPGO as additional covariate
}

# 7. All species figures
source("interaction-plots.R") # Make sst x comp interaction figures
source("era_rw_fig.R") # Make main era + RW results figures
source("sensitivity-analyses/pdo-npgo-figs.R")


time.run <- proc.time() - time.start
round(time.run[3] / 60, 4)
