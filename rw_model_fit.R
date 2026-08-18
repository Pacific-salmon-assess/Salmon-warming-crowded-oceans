
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


# Set paths to output locations - dependent on species
fit.dir <- here("output", "models", "dyn", speciesFlag)
fig.dir <- here("figures", "dyn", speciesFlag, "hbm_inf")

#funcs
make_design_matrix <- function(x,grp) {
  x2 <- matrix(nrow = length(x), ncol = length(unique(grp)))
  for(i in 1:length(unique(grp))){
    x2[,i] <- ifelse(grp == levels(factor(grp))[i], 1, 0)*x
  }
  return(x2)
}

# Run and extract with rstan

# Make dataframe
smax <- data_master %>% group_by(Stock.ID) %>% summarize(max.S = max(S))
S.mat <- make_design_matrix(data_master$S, grp = data_master$Stock.ID)

stk.year <- expand.grid(levels(factor(data_master$Stock.ID)), seq(min(data_master$BY), max(data_master$BY)))
stk.year[,3] <- paste(stk.year[,1], stk.year[,2], sep='.')
stk.year <- stk.year[order(stk.year[,1], stk.year[,2]),]
stk.year.id <- match(paste(data_master$Stock.ID, data_master$BY, sep='.'), stk.year[,3])
if(speciesFlag=="sockeye"){a_grp <- ifelse(info_master$lat == 49.120, 1, 2)} else { # alpha groups - fraser sockeye
  a_grp <- rep(1, nrow(info_master))}

stan.dat.rw <- list(N = nrow(data_master),
           J = length(unique(data_master$Stock.ID)),
           L = max(data_master$BY) - min(data_master$BY) + 1,
           R = length(unique(info_master$Ocean.Region2)),
           #Na = ifelse(speciesFlag=="sockeye", 2, 1),
           #a_grp = a_grp,
           J_i = as.numeric(factor(data_master$Stock.ID)),
           J_ii = stk.year.id,
           J_or = as.numeric(factor(info_master$Ocean.Region2,
                              levels=unique(info_master$Ocean.Region2))),
           R_S = data_master$lnRS,
           S = S.mat,
           X1 = data_master$early_sst_stnd,
           X2 = data_master$np_all_spp_sec_stnd,
           pSmax_mean = 0.5*smax$max.S,
           pSmax_sig = 2*smax$max.S)


# Fit
rw.fit <- rstan::stan(file = "./stan/ind_tvalpha_stat_kappa_ricker.stan", # CHANGED: point at the stationary-kappa Stan model
                      data = stan.dat.rw,
                      warmup = 1000,
                      iter = 2000,
                      cores = 4,
                      chains = 4,
                      seed = 123,
                      control = list(adapt_delta = 0.99,
                                     max_treedepth = 20))
save(rw.fit, file = here("output/models/dyn", speciesFlag, "dyn_new_2025_stat_kappa.RData")) # CHANGED: renamed output so it doesn't overwrite the fully time-varying rw model fit

# Load fit
load(here(fit.dir, "dyn_new_2025_stat_kappa.RData"), verbose=T) # CHANGED: match renamed file above

# Extract

# Stock-specific dataframe
probs <- c(0.025, 0.05, 0.10, 0.50, 0.90, 0.95, 0.975)

# CHANGED: gamma (g_t) is still time-varying [L,J], but kappa is now stationary [J] with no
# year dimension, so it can no longer share one rstan::summary() call / one data.frame
# construction with g_t (their row counts no longer match). Gamma and kappa are now summarized
# and assembled separately, then combined with rbind. Kappa's single per-stock value is
# repeated across all L years so it still plots as a (flat) line alongside gamma's time-varying
# line, using the same Stock/BY layout as before.

summ.g <- rstan::summary(rw.fit, pars = "g_t", probs = probs)[[1]] # CHANGED: gamma summarized on its own

df.dyn.st.2c.g <- data.frame(Stock = rep(info_master$Stock, times=stan.dat.rw$L),
                             Ocean.Region2 = rep(info_master$Ocean.Region2, times=stan.dat.rw$L),
                             BY = rep(min(info_master$yr_start):(min(info_master$yr_start) +
                                                                  stan.dat.rw$L -1),
                                      each=nrow(info_master)),
                             mu = summ.g[, "mean"],
                             se = summ.g[, "se_mean"],
                             lower_2.5 = summ.g[,"2.5%"],
                             lower_10 = summ.g[, "10%"],
                             med = summ.g[,"50%"],
                             upper_90 = summ.g[ , "90%"],
                             upper_97.5 = summ.g[, "97.5%"],
                             var = "gamma",
                             varnam = "SST")

summ.k <- rstan::summary(rw.fit, pars = "kappa", probs = probs)[[1]] # CHANGED: kappa summarized on its own (stationary, no L dimension)

df.dyn.st.2c.k <- data.frame(Stock = rep(info_master$Stock, times=stan.dat.rw$L),                          # CHANGED: kappa block built separately, its single value repeated across years
                             Ocean.Region2 = rep(info_master$Ocean.Region2, times=stan.dat.rw$L),
                             BY = rep(min(info_master$yr_start):(min(info_master$yr_start) +
                                                                  stan.dat.rw$L -1),
                                      each=nrow(info_master)),
                             mu = rep(summ.k[, "mean"], times=stan.dat.rw$L),                                # CHANGED: replicate stationary kappa value across years
                             se = rep(summ.k[, "se_mean"], times=stan.dat.rw$L),
                             lower_2.5 = rep(summ.k[,"2.5%"], times=stan.dat.rw$L),
                             lower_10 = rep(summ.k[, "10%"], times=stan.dat.rw$L),
                             med = rep(summ.k[,"50%"], times=stan.dat.rw$L),
                             upper_90 = rep(summ.k[ , "90%"], times=stan.dat.rw$L),
                             upper_97.5 = rep(summ.k[, "97.5%"], times=stan.dat.rw$L),
                             var = "kappa",
                             varnam = "Competitors")

df.dyn.st.2c <- rbind(df.dyn.st.2c.g, df.dyn.st.2c.k) # CHANGED: combine time-varying gamma rows with the (replicated) stationary kappa rows

if(speciesFlag == "pink"){
  df.dyn.st.2c$even_odd <- ifelse(gtools::odd(df.dyn.st.2c$BY), "odd", "even")
}

# Save outputs
write.csv(df.dyn.st.2c, here(fig.dir, paste0('stk_coefficients_rw_', speciesFlag, '.csv')), row.names=F)


# Summarized dataframe (regional-level)
# CHANGED: same split as the stock-level block above -- mu_g_r stays [L,R], but mu_k_r is now
# a stationary vector[R], so gamma and kappa are summarized/assembled separately and combined.

summ.r.g <- rstan::summary(rw.fit, pars = "mu_g_r", probs = probs)[[1]] # CHANGED: gamma summarized on its own

df.dyn.reg.2c.g <- data.frame(Ocean.Region2 = rep(unique(info_master$Ocean.Region2), times=stan.dat.rw$L),
                              BY = rep(min(info_master$yr_start):(min(info_master$yr_start) +
                                                                    stan.dat.rw$L -1),
                                       each=stan.dat.rw$R),
                              mu = summ.r.g[, "mean"],
                              se = summ.r.g[, "se_mean"],
                              lower_2.5 = summ.r.g[,"2.5%"],
                              lower_10 = summ.r.g[, "10%"],
                              med = summ.r.g[,"50%"],
                              upper_90 = summ.r.g[ , "90%"],
                              upper_97.5 = summ.r.g[, "97.5%"],
                              var = "gamma",
                              varnam = "SST")

summ.r.k <- rstan::summary(rw.fit, pars = "mu_k_r", probs = probs)[[1]] # CHANGED: mu_k_r is now vector[R] (stationary), no year dimension

df.dyn.reg.2c.k <- data.frame(Ocean.Region2 = rep(unique(info_master$Ocean.Region2), times=stan.dat.rw$L), # CHANGED: replicate kappa's single regional value across years for plotting
                              BY = rep(min(info_master$yr_start):(min(info_master$yr_start) +
                                                                    stan.dat.rw$L -1),
                                       each=stan.dat.rw$R),
                              mu = rep(summ.r.k[, "mean"], times=stan.dat.rw$L),
                              se = rep(summ.r.k[, "se_mean"], times=stan.dat.rw$L),
                              lower_2.5 = rep(summ.r.k[,"2.5%"], times=stan.dat.rw$L),
                              lower_10 = rep(summ.r.k[, "10%"], times=stan.dat.rw$L),
                              med = rep(summ.r.k[,"50%"], times=stan.dat.rw$L),
                              upper_90 = rep(summ.r.k[ , "90%"], times=stan.dat.rw$L),
                              upper_97.5 = rep(summ.r.k[, "97.5%"], times=stan.dat.rw$L),
                              var = "kappa",
                              varnam = "Competitors")

df.dyn.reg.2c <- rbind(df.dyn.reg.2c.g, df.dyn.reg.2c.k) # CHANGED: combine time-varying gamma with replicated stationary kappa

if(speciesFlag == "pink"){
  df.dyn.reg.2c$even_odd <- ifelse(gtools::odd(df.dyn.reg.2c$BY), "odd", "even")
}

# Save outputs
write.csv(df.dyn.reg.2c, here(fig.dir, paste0('reg_coefficients_rw_', speciesFlag, '.csv')), row.names=F)


# Trim early year estimates that are not based on data - by average start year of region
reg_start_yr <- info_master %>% group_by(Ocean.Region2) %>% summarize(avg_start=round(mean(yr_start), 0))
df.dyn.reg.2c <- df.dyn.reg.2c %>% left_join(reg_start_yr,  by="Ocean.Region2") %>% filter(BY >= avg_start)
df.dyn.reg.2c <- ocean_region_lab(df.dyn.reg.2c)

df.dyn.st.2c <- ocean_region_lab(df.dyn.st.2c)

# Plot
g <- ggplot(df.dyn.st.2c) +
  geom_hline(aes(yintercept=0), linetype="dashed", col="grey70") +
  geom_line(aes(x=BY, y=mu, group=Stock, col=ocean_region_lab), alpha=0.3) +
  geom_line(data=df.dyn.reg.2c, aes(x=BY, y=mu, col=ocean_region_lab), linewidth=1) +
  geom_ribbon(data=df.dyn.reg.2c, aes(x=BY, ymin=lower_10, ymax=upper_90, fill=ocean_region_lab),
              alpha=0.3) +
  scale_colour_manual(values=col.region, aesthetics = c("colour", "fill"), guide="none") +
  facet_grid(rows=vars(Ocean.Region2), cols=vars(varnam)) +
  ylim(-1,1) + labs(x="Brood Year", y="Covariate effect") +
  theme_sleek()

png(file=here(fig.dir, paste0("rw_2025", speciesFlag, ".png")))
print(g)
dev.off()


