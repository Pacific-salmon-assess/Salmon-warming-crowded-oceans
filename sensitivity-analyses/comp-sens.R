## Competitor Index sensitivity analyses ##

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

# Load competitor data
raw.comp <- read.csv(file="data-downloaded/competitor_indices_2024.csv", header = TRUE)
# Choose alt. competitor indices to run
comp.sec <- data.frame(BY=NA, Stock.ID=NA)
comp.cols <- names(raw.comp)[c(5,8,9)]#names(raw.comp)[-1]
comp.cols.stnd <- paste0(comp.cols, "_stnd")

# Load correct brood table
brood <- if(speciesFlag=="sockeye"){
  read.csv("./data/sockeye/master_sockeye_brood_table.csv", header=T)} else if(speciesFlag=="pink") {
    read.csv("./data/pink/master_pink_brood_table.csv", header=T)
  } else if (speciesFlag=="chum"){
    read.csv("./data/chum/master_chum_brood_table.csv", header=T)}
brood <- brood[complete.cases(brood),]

# Create and bind stock-specific (for sockeye) comp indices in a loop
if(speciesFlag == "sockeye"){
  for( ind in comp.cols) {  # stock-specific weighted average
    temp <- pink.wgt.avg(brood.table = brood,
                       pink.data = raw.comp,
                       pink.covar = ind,
                       type = "second_year",
                       out.covar = ind)
  comp.sec <- full_join(comp.sec, temp, by=c("BY", "Stock.ID"))
  }
brood_comp <- dplyr::left_join(data_master, comp.sec, by=c("BY","Stock.ID")) # bind new comp indices to full brood table
} else {
  raw.comp$BY <- raw.comp$Year - 2
  brood_comp <- left_join(data_master, raw.comp[, c(comp.cols, "BY")], by= "BY")
}


# Scale all indices
brood_comp[,comp.cols.stnd] <- NA
brood_comp[,comp.cols.stnd] <- lapply(brood_comp[comp.cols], function(x){scale(x)[,1]})

# Order geographically
brood_comp <- geographic.order(brood_comp)


## Get correlation between indices and with response


## create empty 3D array
array.cor <- array(NA, dim = c(length(comp.cols)+1,
                               length(comp.cols)+1,
                               nrow(info_master)))


## calculate stock specific covar correlations
for(i in seq_along(info_master$Stock.ID)) {
  stk.i <- brood_comp[brood_comp$Stock.ID == info_master$Stock.ID[i], ]
  covar.i <- subset(stk.i, select = c(comp.cols, "lnRS"))
  cor.i <- cor(covar.i, use = "pairwise.complete.obs")
  array.cor[ , , i] <- cor.i
}

## Average across stocks
cor.covars <- apply(array.cor, c(1, 2), mean)
row.names(cor.covars) <- c(comp.cols, "lnRS")
colnames(cor.covars) <- c(comp.cols, "lnRS")

cols <- chroma::dpal(500, hue = c(240, 0), chroma = 70, power = 1.0)
# Plot correlation between Comp covariates
g <- levelplot(cor.covars[comp.cols, comp.cols], xlab = "", ylab = "",
               col.regions = cols,
               #at = seq(-0.3, 1, 0.1),
               main = "Average stock-specific covariate correlations",
               scales = list(x=list(rot=45)),
               panel = function(x, y, ...) {
                 panel.levelplot(x, y, ...)
                 panel.abline(v = seq(1.5, max(x) - 0.5), col = "grey50")
                 panel.abline(h = seq(1.5, max(y) - 0.5), col = "grey50")
               })
pdf(here('sensitivity-analyses', 'corr-comp-covars.pdf'))
print(g)
dev.off()

# Plot correlation between each and lnRS
g <- ggplot(data.frame(covar=comp.cols,
                       cor=cor.covars[comp.cols,"lnRS"])) +
  geom_col(aes(x=covar, y=cor, fill=covar), position=position_dodge(), alpha=0.7) +
  scale_fill_brewer(palette="Paired") +
  theme_sleek() + theme(axis.text.x=element_text(angle=90),
                        legend.position="none") +
  labs(x="", y="Average Correlation", fill="Covariate",
                       title="All Years")
pdf(here('sensitivity-analyses', 'corr-comp-lnRS.pdf'))
print(g)
dev.off()

# Run sensitivity analyses

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

era_data_list <- vector("list", length=length(comp.cols))
stat_data_list <- vector("list", length=length(comp.cols))
names(era_data_list) <- names(stat_data_list) <- comp.cols

for(i in 1:length(comp.cols.stnd)){
  name_temp <- as.character(comp.cols.stnd[i])
  era_data_list[[i]] <- stan_data_dyn(brood_comp,
                                var.x2 = "early_sst_stnd",
                                var.x3 = name_temp,
                                breakpoint1 = 1989,
                                breakpoint2 = 2011,
                                var.region="Ocean.Region2",
                                scale.x1 = TRUE,
                                alpha.group = ifelse(speciesFlag=="sockeye", TRUE, FALSE))

  stat_data_list[[i]] <- stan_data_stat(brood_comp,
                                        scale.x1 = TRUE,
                                        var.x2 = "early_sst_stnd",
                                        var.x3 = name_temp,
                                        var.region = "Ocean.Region2",
                                        alpha.group = ifelse(speciesFlag=="sockeye", TRUE, FALSE))
  }


# Run models in a loop and extract key parameters
era_fit_list <- vector("list", length = length(comp.cols))
stat_fit_list <- vector("list", length = length(comp.cols))
era_summary_list <- vector("list", length = length(comp.cols))
stat_summary_list <- vector("list", length = length(comp.cols))
names(era_fit_list) <- names(stat_fit_list) <- comp.cols


for(i in 1:length(comp.cols.stnd)){ # Loop to run stan models:

  # era model
  era.fit <- rstan::stan(file = "./stan/hbm_era_2c.stan",
                               data = era_data_list[[i]],
                               pars = c(pars_era_2c, pars.gen.quant),
                               warmup = 1000,
                               iter = 3000,
                               cores = 4,
                               chains = 4,
                               seed = 123,
                               control = list(adapt_delta = 0.999,
                                              max_treedepth = 20))

    save(era.fit, file=here('sensitivity-analyses', 'fits', speciesFlag, paste0("era_comp_", comp.cols[i], ".RData")))


  era_summary_list[[i]] <- rstan::summary(era.fit, pars=c("alpha", "beta", "sigma",
                                                  "mu_gamma1", "mu_gamma2", "mu_gamma3",
                                                  "mu_kappa1", "mu_kappa2", "mu_kappa3",
                                                  "log_lik"),
                                          probs=c(0.025, 0.1, 0.5, 0.9, 0.975))$summary


  stat.fit <- rstan::stan(file = "./stan/hbm_stat_inter.stan",
                         data = stat_data_list[[i]],
                         warmup = 1000,
                         iter = 3000,
                         cores = 4,
                         chains = 4,
                         seed = 123,
                         control = list(adapt_delta = 0.999,
                                        max_treedepth = 20))


  save(stat.fit, file = here('sensitivity-analyses', 'fits', speciesFlag, paste0("stat_comp_", comp.cols[i], ".RData")))



  # store summary
  stat_summary_list[[i]] <- rstan::summary(stat.fit, pars=c("alpha", "beta", "sigma", "mu_gamma",
                                                            "mu_kappa", "mu_chi", "log_lik"),
                                           probs=c(0.025, 0.1, 0.5, 0.9, 0.975))$summary


}

# Make master comparison table

comp_tbl_list_era <- list()
comp_tbl_list_era <- lapply(era_summary_list, function(x){
        summ <- x[c(grepl("^mu_", rownames(x))),
                          c("mean", "2.5%", "10%", "50%", "90%", "97.5%")]
        df <- data.frame(
                   Parameter = case_when(grepl("^mu_gamma", rownames(summ)) ~ "SST effect",
                                         grepl("^mu_kappa", rownames(summ)) ~ "Competitor effect"),
                   Era = case_when(str_sub(rownames(summ), start=-4, end=-4) == "1" ~ "Early",
                                   str_sub(rownames(summ), start=-4, end=-4) == "2" ~ "Middle",
                                   str_sub(rownames(summ), start=-4, end=-4) == "3" ~ "Late"),
                   Region = case_when(str_sub(rownames(summ), start=-2,end=-2) == "1" ~ "West Coast",
                                      str_sub(rownames(summ), start=-2,end=-2) == "2" ~ "Southeast Alaska",
                                      str_sub(rownames(summ), start=-2,end=-2) == "3" ~ "Gulf of Alaska",
                                      str_sub(rownames(summ), start=-2,end=-2) == "4" ~ "Bering Sea",),
                   Posterior_Mean = summ[,"mean"],
                   Lower_2.5_CI = summ[,"2.5%"],
                   Lower_10_CI = summ[,"10%"],
                   Posterior_Median = summ[,"50%"],
                   Upper_90_CI = summ[,"90%"],
                   Upper_97.5_CI = summ[,"97.5%"],
                   row.names=NULL)
        return(df)
      }
    )

comp_tbl_list_stat <- lapply(stat_summary_list, function(x){
  summ <- x[c(grepl("^mu_", rownames(x))),
            c("mean", "2.5%", "10%", "50%", "90%", "97.5%")]
  df <- data.frame(
    Parameter = case_when(grepl("^mu_gamma", rownames(summ)) ~ "SST effect",
                          grepl("^mu_kappa", rownames(summ)) ~ "Competitor effect",
                          grepl("^mu_chi", rownames(summ)) ~ "SST x Competitor effect"),
    Era = "All",
    Region = case_when(str_sub(rownames(summ), start=-2,end=-2) == "1" ~ "West Coast",
                       str_sub(rownames(summ), start=-2,end=-2) == "2" ~ "Southeast Alaska",
                       str_sub(rownames(summ), start=-2,end=-2) == "3" ~ "Gulf of Alaska",
                       str_sub(rownames(summ), start=-2,end=-2) == "4" ~ "Bering Sea",),
    Posterior_Mean = summ[,"mean"],
    Lower_2.5_CI = summ[,"2.5%"],
    Lower_10_CI = summ[,"10%"],
    Posterior_Median = summ[,"50%"],
    Upper_90_CI = summ[,"90%"],
    Upper_97.5_CI = summ[,"97.5%"],
    row.names=NULL)
  return(df)
}
)

names(comp_tbl_list_era) <- names(comp_tbl_list_stat) <- comp.cols
comp_tbl_era <- bind_rows(comp_tbl_list_era, .id="Competitor index")
comp_tbl_stat <- bind_rows(comp_tbl_list_stat, .id="Competitor index")

comp_tbl <- bind_rows(comp_tbl_era, comp_tbl_stat, .id="Model")

write.csv(comp_tbl, file=here('sensitivity-analyses', 'alt-comp', paste0(speciesFlag, '-comp-index-fit-compare.csv')))

