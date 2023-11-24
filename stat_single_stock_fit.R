## Fit single-stock generalized Ricker models
##
## This script fits and explores single-stock versions of our proposed
## hierarchical models. The parameters in these single-stock models are
## estimated only using data from a single stock, rather than being informed
## from data from multiple stocks in the hierarchical models.
##
## The reasons for exploring these single-stock models include:
##  1. Have a base set of results in which to compare the hierarchical models
##  2. Determine the most appropriate hierarchical structure for the parameters
##  3. Explore model assumptions, e.g., normality, co-linearity, etc.

## Set Species -----
speciesFlag = "pink"
#speciesFlag = "chum"
#speciesFlag = "sockeye"

if(speciesFlag=="pink") 
  data_master <- pink else if(speciesFlag=="chum") 
    data_master <- chum else if(speciesFlag=="sockeye")
      data_master <- sock

## Make directories to save plots
fig.dir <- here("figures", "stat", speciesFlag, "single-stock")
fit.dir <- here("output", "models", "stat", speciesFlag, "single-stock")

if(!dir.exists(fig.dir)) dir.create(fig.dir, recursive = T)
if(!dir.exists(fit.dir)) dir.create(fit.dir, recursive = T)


## Specify model formulas ----------------------------------

m1a.formula <- lnRS ~ S | Stock
m2a.formula <- lnRS ~ S + early_sst_stnd | Stock
m3a.formula <- lnRS ~ S + np_pinks_sec_stnd | Stock
m4a.formula <- lnRS ~ S + early_sst_stnd + np_pinks_sec_stnd | Stock


mod.list <- list(model1a = m1a.formula,
                 model2a = m2a.formula,
                 model3a = m3a.formula,
                 model4a = m4a.formula)

## Fit single-stock models ---------------------------------

## All years of data
ss.all.yrs <- single.stock.fit(mod.list,
                               years = seq(min(sock$BY), max(sock$BY)),
                               plot.path = fig.dir)

# Save
save(ss.all.yrs, file = here(fit.dir, "single_stock_lms.Rdata"))



