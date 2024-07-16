## Competitor Index sensitivity analyses ##
# Sockeye only for now, using era models


# Load competitor data
raw.comp <- read.csv(file="data-downloaded/competitor_indices_2024.csv", header = TRUE)


comp.sec <- data.frame(BY=NA, Stock.ID=NA)
comp.cols <- names(raw.comp)[-1]
comp.cols.stnd <- paste0(comp.cols, "_stnd")

for( ind in comp.cols) { # This loop creates and binds stock-specific comp indices
  temp <- pink.wgt.avg(brood.table = bt.complete, # source sock_covariates.R to make sure bt.complete is right
                       pink.data = raw.comp,
                       pink.covar = ind,
                       type = "second_year",
                       out.covar = ind)
  comp.sec <- full_join(comp.sec, temp, by=c("BY", "Stock.ID"))
}

sock_comp <- dplyr::left_join(sock, comp.sec, by=c("BY","Stock.ID")) # bind to brood table

# Scale all indices
sock_comp[,comp.cols.stnd] <- NA
sock_comp[,comp.cols.stnd] <- lapply(sock_comp[comp.cols], scale)


## Get correlation between indices and with response


## create empty 3D array
array.cor <- array(NA, dim = c(length(comp.cols)+1,
                               length(comp.cols)+1,
                               nrow(info_master)))


## calculate stock specific covar correlations
for(i in seq_along(info_master$Stock.ID)) { 
  stk.i <- sock_comp[sock_comp$Stock.ID == info_master$Stock.ID[i], ]
  covar.i <- subset(stk.i, select = c(comp.cols, "lnRS"))
  cor.i <- cor(covar.i, use = "pairwise.complete.obs")
  array.cor[ , , i] <- cor.i
}

## Average across stocks
cor.covars <- apply(array.cor, c(1, 2), mean)
row.names(cor.covars) <- c(comp.cols, "lnRS")
colnames(cor.covars) <- c(comp.cols, "lnRS")

cols <- chroma::dpal(500, hue = c(240, 0), chroma = 70, power = 1.0)
levelplot(cor.covars, xlab = "", ylab = "",
               col.regions = cols,
               at = seq(-0.3, 1, 0.1),
               main = "Average stock-specific covariate correlations",
               scales = list(x=list(rot=45)),
               panel = function(x, y, ...) {
                 panel.levelplot(x, y, ...)
                 panel.abline(v = seq(1.5, max(x) - 0.5), col = "grey50")
                 panel.abline(h = seq(1.5, max(y) - 0.5), col = "grey50")
               })
# better to split the above into 2 figs due to different cor values


# Run senstivity analyses

# Set up model params
pars_era_2c <- c("alpha", "beta", "sigma", "phi", "mu_alpha", "sigma_alpha",
                 "gamma1", "gamma2", "gamma3", 
                 "mu_gamma1", "mu_gamma2", "mu_gamma3", "sigma_gamma", 
                 "kappa1", "kappa2", "kappa3", 
                 "mu_kappa1", "mu_kappa2", "mu_kappa3", "sigma_kappa" )

pars.gen.quant <- c("log_lik", "yrep", "yhat", "yresid") ## Generated quantities to monitor

# Get model data in a list
  
  ## INDEX:
  #  1. North American Pinks (pink_numbers_na_stnd)
  # 2. North American Pink, Chum, Sockeye (all_spp_numbers_na_stnd)
  # 3. North Pacific Pink, Chum, Sockeye (all_spp_numbers_np_stnd)
  # 4. North American Pink Biomass (pink_biomass_na_stnd)
  # 5. North American Pink, Chum, Sockeye Biomass (all_spp_biomass_na_stnd)
  # 6. North Pacific Pink Biomass (pink_biomass_np_stnd)
  # 7. (BASE) North Pacific Pinks (pink_numbers_np_stnd)
  # 8. North Pacific Pink, Chum, Sockeye Biomass (all_spp_biomass_np_stnd)

data_list <- vector("list", length=length(comp.cols))
names(data_list) <- comp.cols

for(i in 1:length(comp.cols)){
  name_temp <- as.character(comp.cols[i])
  data_list[[i]] <- stan_data_dyn(sock_comp, 
                                var.x2 = "early_sst_stnd",
                                var.x3 = name_temp,
                                breakpoint1 = 1989,
                                breakpoint2 = 2011,
                                var.region="Ocean.Region2", 
                                scale.x1 = TRUE,
                                alpha.group = TRUE)
  }


# Run era models in a loop and extract key parameters
fit_list <- vector("list", length = length(comp.cols))
summary_list <- vector("list", length = length(comp.cols))
names(fit_list) <- comp.cols
for(i in 1:length(comp.cols)){
  fit <- rstan::stan(file = "./stan/hbm_era_2c.stan",
                               data = data_list[[i]],
                               pars = c(pars_era_2c, pars.gen.quant),
                               warmup = 1000,
                               iter = 2000,
                               cores = 4,
                               chains = 4,
                               seed = 123,
                               control = list(adapt_delta = 0.9,
                                              max_treedepth = 20))
  fit.list[[i]] <- fit  
  save(fit.list[[i]], file=here('sensitivity-analyses', 'fits', paste0(comp_, comp.cols[i], .RData)))
  
  summary_list[[i]] <- rstan::summary(fit, pars=c("alpha", "beta", "sigma",
                                                  "mu_gamma1", "mu_gamma2", "mu_gamma3", "sigma_gamma", 
                                                  "mu_kappa1", "mu_kappa2", "mu_kappa3", "sigma_kappa", "log_lik"), probs=c(0.025, 0.5, 97.5))$summary
  
}

# Run model selection
# Make master comparison table


