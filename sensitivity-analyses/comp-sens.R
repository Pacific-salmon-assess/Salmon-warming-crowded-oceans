## Competitor Index sensitivity analyses ##

# Load competitor data
raw.comp <- read.csv(file="data-downloaded/competitor_indices_2024.csv", header = TRUE)


comp.sec <- data.frame(BY=NA, Stock.ID=NA)
comp.cols <- names(raw.comp)[-1]

for( ind in comp.cols) { # This loop creates and binds stock-specific comp indices
  temp <- pink.wgt.avg(brood.table = bt.complete,
                       pink.data = raw.comp,
                       pink.covar = ind,
                       type = "second_year",
                       out.covar = ind)
  comp.sec <- full_join(comp.sec, temp, by=c("BY", "Stock.ID"))
}

sock_comp <- dplyr::left_join(sock, comp.sec, by=c("BY","Stock.ID")) # bind to brood table

# Scale all indices
sock_comp[comp.cols] <- lapply(sock_comp[comp.cols], scale)
names(sock_comp)[names(sock_comp) %in% comp.cols] <- paste0(comp.cols, "_stnd")



## most of the below not necessary with this setup 


# Sockeye only for now

# load updated data
comp.dat <- read.csv(here('data-downloaded', 'competitor_indices_2024.csv'), na.strings = "")

# 1. run the models for one species with each alternative index (start with era)

source("sock_data_clean.R")
#source("sst_import_process.R")
sock_comp <- bt.out


pars_era_2c <- c("alpha", "beta", "sigma", "phi", "mu_alpha", "sigma_alpha",
                 "gamma1", "gamma2", "gamma3", 
                 "mu_gamma1", "mu_gamma2", "mu_gamma3", "sigma_gamma", 
                 "kappa1", "kappa2", "kappa3", 
                 "mu_kappa1", "mu_kappa2", "mu_kappa3", "sigma_kappa" )

pars.gen.quant <- c("log_lik", "yrep", "yhat", "yresid") ## Generated quantities to monitor


  # 1)'BASE' model - NP pink abundance
  
  # 2) NP all-species abundance [pink, chum, sockeye]
  
all.sp.num.wgt <- pink.wgt.avg(brood.table = sock_comp,
                            pink.data = comp.dat,
                            pink.covar = "all_sp_numbers_np",
                            type = "second_year",
                            out.covar = "all_sp_numbers_np")

sock_comp_ext <- dplyr::left_join(sock_comp, all.sp.num.wgt, by=c("BY","Stock.ID"))

stan.dat.2c <- stan_data_dyn(sock_comp_ext, 
                             var.x2 = "early_sst_stnd",
                             var.x3 = "all_sp_numbers_np",
                             breakpoint1 = 1989,
                             breakpoint2 = 2011,
                             var.region="Ocean.Region2", 
                             scale.x1 = TRUE,
                             alpha.group = TRUE) # make FALSE for pink/chum

era.comp.np.all <- rstan::stan(file = "./stan/hbm_era_2c.stan",
                                data = stan.dat.2c,
                                pars = c(pars_era_2c, pars.gen.quant),
                                warmup = 1000,
                                iter = 2000,
                                cores = 4,
                                chains = 4,
                                seed = 123,
                                control = list(adapt_delta = 0.99,
                                               max_treedepth = 20))

  # 3) NA pink abundance
  
na.pink.num.wgt <- pink.wgt.avg(brood.table = sock_comp,
                                pink.data = comp.dat,
                                pink.covar = "pink_numbers_na",
                                type = "second_year",
                                out.covar = "pink_numbers_na")

  # 4) NP Hatchery pink abundance [not avail.]
  
  # 5) NP Hatchery all-species abundance [not avail.]
  
  # 6) NP all-species biomass
all.sp.bio.wgt <- pink.wgt.avg(brood.table = sock_comp,
                               pink.data = comp.dat,
                               pink.covar = "all_sp_biomass_np",
                               type = "second_year",
                               out.covar = "all_sp_biomass_np")
  
  
# 2. store coefficients in table
# 3. run model selection
# 4. store model selection results
# 5. make master comparison table


