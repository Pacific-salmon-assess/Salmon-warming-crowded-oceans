# Create table output of era results


## Sockeye ##

load(here('output', 'models', 'dyn', 'sockeye', 'hbm_era_2c.RData'), verbose=T) # eras

# Stock-specific dataframe
probs <- c(0.025, 0.50, 0.975)
summ <- rstan::summary(era.2c, pars = c(paste0("gamma", 1:3), paste0("kappa", 1:3)), probs = probs)[[1]]
sock.df <- data.frame( Species = "Sockeye",
                       Stock = rep(sock.info$Stock, times=6),
                       Stock.ID = rep(sock.info$Stock.ID, times=6),
                       mean = summ[, "mean"],
                       lower_2.5 = summ[, "2.5%"],
                       upper_97.5 = summ[ , "97.5%"],
                       varnam = case_when(grepl("^gamma", rownames(summ)) ~ "SST",
                                          grepl("^kappa", rownames(summ)) ~ "Competitors"),
                       era = case_when(str_extract(rownames(summ), "\\d")==1 ~ "Early",
                                       str_extract(rownames(summ), "\\d")==2 ~ "Middle",
                                       str_extract(rownames(summ), "\\d")==3 ~ "Late"),
                       region = rep(sock.info$Ocean.Region2, times=6))
# Region-lvl dataframe
summ_avg <- rstan::summary(era.2c, pars = c(paste0("mu_gamma", 1:3), paste0("mu_kappa", 1:3)), probs = probs)[[1]]
sock.df.avg <- data.frame( Species = "Sockeye",
                           region = rep(unique(sock.info$Ocean.Region2), times=6),
                           mean = summ_avg[, "mean"],
                           lower_2.5 = summ_avg[, "2.5%"],
                           upper_97.5 = summ_avg[ , "97.5%"],
                           varnam = case_when(grepl("^mu_gamma", rownames(summ_avg)) ~ "SST",
                                              grepl("^mu_kappa", rownames(summ_avg)) ~ "Competitors"),
                           era = case_when(str_extract(rownames(summ_avg), "\\d")==1 ~ "Early",
                                           str_extract(rownames(summ_avg), "\\d")==2 ~ "Middle",
                                           str_extract(rownames(summ_avg), "\\d")==3 ~ "Late")
                           )

# Output SST values - anomalies and raw
raw.sst.sock <- clim.wgt.avg(brood.table = sock,
                              env.data = raw.clim.sock,
                              env.covar = "sst_raw",
                              type = "first_year",
                              out.covar = "early_sst")

cov.dat.sock <- full_join(raw.sst.sock, sock[,c("Stock.ID", "Stock", "Ocean.Region2", "BY", "early_sst", "early_sst_stnd", "lnRS")] , by=c("Stock.ID", "BY"))
cov.dat.sock <- cov.dat.sock %>% rename(early_sst_raw = early_sst.x, 
                                        early_sst_anom = early_sst.y,
                                        early_sst_anom_stnd = early_sst_stnd)


## Pink ##

load(here('output', 'models', 'dyn', 'pink', 'hbm_era_2c.RData'), verbose=T) # eras

# Stock-specific dataframe
probs <- c(0.025, 0.50, 0.975)
summ <- rstan::summary(era.2c, pars = c(paste0("gamma", 1:3), paste0("kappa", 1:3)), probs = probs)[[1]]
pink.df <- data.frame( Species = "Pink",
                       Stock = rep(pink.info$Stock, times=6),
                       mean = summ[, "mean"],
                       lower_2.5 = summ[, "2.5%"],
                       upper_97.5 = summ[ , "97.5%"],
                       varnam = case_when(grepl("^gamma", rownames(summ)) ~ "SST",
                                          grepl("^kappa", rownames(summ)) ~ "Competitors"),
                       era = case_when(str_extract(rownames(summ), "\\d")==1 ~ "Early",
                                       str_extract(rownames(summ), "\\d")==2 ~ "Middle",
                                       str_extract(rownames(summ), "\\d")==3 ~ "Late"),
                       region = rep(pink.info$Ocean.Region2, times=6))
# Region-lvl dataframe
summ_avg <- rstan::summary(era.2c, pars = c(paste0("mu_gamma", 1:3), paste0("mu_kappa", 1:3)), probs = probs)[[1]]
pink.df.avg <- data.frame( Species = "pink",
                           region = rep(unique(pink.info$Ocean.Region2), times=6),
                           mean = summ_avg[, "mean"],
                           lower_2.5 = summ_avg[, "2.5%"],
                           upper_97.5 = summ_avg[ , "97.5%"],
                           varnam = case_when(grepl("^mu_gamma", rownames(summ_avg)) ~ "SST",
                                              grepl("^mu_kappa", rownames(summ_avg)) ~ "Competitors"),
                           era = case_when(str_extract(rownames(summ_avg), "\\d")==1 ~ "Early",
                                           str_extract(rownames(summ_avg), "\\d")==2 ~ "Middle",
                                           str_extract(rownames(summ_avg), "\\d")==3 ~ "Late")
)

# Output SST values - anomalies and raw
raw.sst.pink <- raw.clim.pink[, c("BY","Stock.ID", "sst_raw")]

cov.dat.pink <- right_join(raw.sst.pink, pink[,c("Stock.ID", "Stock", "Ocean.Region2", "BY", "early_sst", "early_sst_stnd", "lnRS")] , by=c("Stock.ID", "BY"))
cov.dat.pink <- cov.dat.pink %>% rename(early_sst_raw = sst_raw, 
                                        early_sst_anom = early_sst,
                                        early_sst_anom_stnd = early_sst_stnd)
  


## Chum ##

load(here('output', 'models', 'dyn', 'chum', 'hbm_era_2c.RData'), verbose=T) # eras

# Stock-specific dataframe
probs <- c(0.025, 0.50, 0.975)
summ <- rstan::summary(era.2c, pars = c(paste0("gamma", 1:3), paste0("kappa", 1:3)), probs = probs)[[1]]
chum.df <- data.frame( Species = "Chum",
                       Stock = rep(chum.info$Stock, times=6),
                       mean = summ[, "mean"],
                       lower_2.5 = summ[, "2.5%"],
                       upper_97.5 = summ[ , "97.5%"],
                       varnam = case_when(grepl("^gamma", rownames(summ)) ~ "SST",
                                          grepl("^kappa", rownames(summ)) ~ "Competitors"),
                       era = case_when(str_extract(rownames(summ), "\\d")==1 ~ "Early",
                                       str_extract(rownames(summ), "\\d")==2 ~ "Middle",
                                       str_extract(rownames(summ), "\\d")==3 ~ "Late"),
                       region = rep(chum.info$Ocean.Region2, times=6))
# Region-lvl dataframe
summ_avg <- rstan::summary(era.2c, pars = c(paste0("mu_gamma", 1:3), paste0("mu_kappa", 1:3)), probs = probs)[[1]]
chum.df.avg <- data.frame( Species = "chum",
                           region = rep(unique(chum.info$Ocean.Region2), times=6),
                           mean = summ_avg[, "mean"],
                           lower_2.5 = summ_avg[, "2.5%"],
                           upper_97.5 = summ_avg[ , "97.5%"],
                           varnam = case_when(grepl("^mu_gamma", rownames(summ_avg)) ~ "SST",
                                              grepl("^mu_kappa", rownames(summ_avg)) ~ "Competitors"),
                           era = case_when(str_extract(rownames(summ_avg), "\\d")==1 ~ "Early",
                                           str_extract(rownames(summ_avg), "\\d")==2 ~ "Middle",
                                           str_extract(rownames(summ_avg), "\\d")==3 ~ "Late")
)

# Output SST values - anomalies and raw
raw.sst.chum <- raw.clim.chum[, c("BY","Stock.ID", "sst_raw")]

cov.dat.chum <- right_join(raw.sst.chum, chum[,c("Stock.ID", "Stock", "Ocean.Region2", "BY", "early_sst", "early_sst_stnd", "lnRS")] , by=c("Stock.ID", "BY"))
cov.dat.chum <- cov.dat.chum %>% rename(early_sst_raw = sst_raw, 
                                        early_sst_anom = early_sst,
                                        early_sst_anom_stnd = early_sst_stnd)




# Bind dataframes
eras_out_df <- rbind(sock.df, pink.df, chum.df)
rownames(eras_out_df) <- NULL

eras_avg_out_df <- rbind(sock.df.avg, pink.df.avg, chum.df.avg)
rownames(eras_avg_out_df) <- NULL

raw_sst_out <- bind_rows(cov.dat.sock, cov.dat.pink, cov.dat.chum, .id="Species")
raw_sst_out <- raw_sst_out %>% mutate(Species = 
                                      case_when(Species == 1 ~ "Sockeye",
                                                Species == 2 ~ "Pink",
                                                Species == 3 ~ "Chum")) %>%
  relocate(early_sst_raw, .after=Ocean.Region2)


write.csv(eras_out_df, file="./output/models/era_stk_coefs.csv", row.names = F)
write.csv(eras_avg_out_df, file="./output/models/era_avg_coefs.csv", row.names = F)
write.csv(raw_sst_out, file = "./data/raw_sst_formatted.csv", row.names = F)

