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
                               plot.path = "./figures/stat/single-stock/")

if(!dir.exists("./output/models/stat/single-stock")) dir.create("./output/models/stat/single-stock")

save(ss.all.yrs, file = "./output/models/stat/single-stock/single_stock_lms.Rdata")



