## Fit single-stock generalized Ricker models 
## Non-hierarchical equivalents of HBM for comparison, testing model assumptions

# Species
if(speciesFlag=="pink") {
  data_master <- pink
  info_master <- pink.info} else if (speciesFlag=="pinkeven"){
    data_master <- pinkeven
    info_master <- pinkeven.info } else if (speciesFlag=="pinkodd"){
      data_master <- pinkodd
      info_master <- pinkodd.info } else if (speciesFlag=="chum") {
        data_master <- chum
        info_master <- chum.info } else if(speciesFlag=="sockeye"){
          data_master <- sock
          info_master <- sock.info }

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
                               years = seq(min(data_master$BY), max(data_master$BY)),
                               plot.path = fig.dir)

# Save
save(ss.all.yrs, file = here(fit.dir, "single_stock_lms.Rdata"))



