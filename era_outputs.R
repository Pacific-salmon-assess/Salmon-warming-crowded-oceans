# Create table output of era results


## Sockeye ##

load(here('output', 'models', 'dyn', 'sockeye', 'hbm_era_2c.RData'), verbose=T) # eras

# Stock-specific dataframe
probs <- c(0.025, 0.50, 0.975)
summ <- rstan::summary(era.2c, pars = c(paste0("gamma", 1:3), paste0("kappa", 1:3)), probs = probs)[[1]]
sock.df <- data.frame( Species = "Sockeye",
                       Stock = rep(sock.info$Stock, times=6),
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

# Bind dataframes
eras_out_df <- rbind(sock.df, pink.df, chum.df)
rownames(eras_out_df) <- NULL

eras_avg_out_df <- rbind(sock.df.avg, pink.df.avg, chum.df.avg)
rownames(eras_avg_out_df) <- NULL

write.csv(eras_out_df, file="./output/models/era_stk_coefs.csv", row.names = F)
write.csv(eras_avg_out_df, file="./output/models/era_avg_coefs.csv", row.names = F)


