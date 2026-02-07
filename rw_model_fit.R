
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
           X2 = data_master$np_pinks_sec_stnd,
           pSmax_mean = 0.5*smax$max.S,
           pSmax_sig = 2*smax$max.S)


# Fit
rw.fit <- rstan::stan(file = "./stan/ind_tvalpha_ricker.stan", # test running "-Copy"
                      data = stan.dat.rw,
                      warmup = 1000,
                      iter = 3000,
                      cores = 4,
                      chains = 4,
                      seed = 123,
                      control = list(adapt_delta = 0.999,
                                     max_treedepth = 20))
save(rw.fit, file = here("output/models/dyn", speciesFlag, "dyn_new_2025.RData"))

# Load fit
load(here("output/models/dyn", speciesFlag, "dyn_new_2025.RData"), verbose=T)

# Extract

# Stock-specific dataframe
probs <- c(0.025, 0.05, 0.10, 0.50, 0.90, 0.95, 0.975)
summ <- rstan::summary(rw.fit, pars = c("g_t", "k_t"), probs = probs)[[1]]

df.dyn.st.2c <- data.frame(Stock = rep(info_master$Stock, times=stan.dat.rw$L*2),
                           Ocean.Region2 = rep(info_master$Ocean.Region2, times=stan.dat.rw$L*2),
                           BY = rep(rep(min(info_master$yr_start):(min(info_master$yr_start) +
                                                                   stan.dat.rw$L -1),
                                        each=nrow(info_master)), times=2),
                           mu = summ[, "mean"],
                           se = summ[, "se_mean"],
                           lower_10 = summ[, "10%"],
                           upper_90 = summ[ , "90%"],
                           var = case_when(str_extract(rownames(summ), "[a-z]+") == "g" ~ "gamma",
                                        str_extract(rownames(summ), "[a-z]+") == "k" ~ "kappa"),
                           varnam = case_when(grepl("^g", rownames(summ)) ~ "SST",
                                              grepl("^k", rownames(summ)) ~ "Competitors"))
df.dyn.st.2c <- ocean_region_lab(df.dyn.st.2c)


# Summarized dataframe (regional-level)
summ.r <- rstan::summary(rw.fit, pars = c("mu_g_r", "mu_k_r"), probs = probs)[[1]]
df.dyn.reg.2c <- data.frame(Ocean.Region2 = rep(unique(info_master$Ocean.Region2),
                                                  times=stan.dat.rw$L*2),
                              BY = rep(rep(min(info_master$yr_start):(min(info_master$yr_start) +
                                                                      stan.dat.rw$L -1),
                                           each=stan.dat.rw$R), times=2),
                              mu = summ.r[, "mean"],
                              se = summ.r[, "se_mean"],
                              lower_10 = summ.r[, "10%"],
                              upper_90 = summ.r[ , "90%"],
                              var = case_when(grepl("^mu_g", rownames(summ.r)) ~ "gamma",
                                              grepl("^mu_k", rownames(summ.r)) ~ "kappa"),
                              varnam = case_when(grepl("^mu_g", rownames(summ.r)) ~ "SST",
                                                 grepl("^mu_k", rownames(summ.r)) ~ "Competitors")
                              )
if(speciesFlag == "pink"){
  df.dyn.reg.2c$even_odd <- ifelse(gtools::odd(df.dyn.reg.2c$BY), "odd", "even")
}

# Save outputs
write.csv(df.dyn.reg.2c, here('output', paste0('rw_2025_mod_outputs_', speciesFlag, '.csv')), row.names=F)



# Trim early year estimates that are not based on data - by average start year of region (NOT an ideal way, more nuance needed...)
reg_start_yr <- info_master %>% group_by(Ocean.Region2) %>% summarize(avg_start=round(mean(yr_start), 0))
df.dyn.reg.2c <- df.dyn.reg.2c %>% left_join(reg_start_yr,  by="Ocean.Region2") %>% filter(BY >= avg_start)
df.dyn.reg.2c <- ocean_region_lab(df.dyn.reg.2c)


# Plot
ggplot(df.dyn.st.2c) +
  geom_line(aes(x=BY, y=mu, group=Stock), col="grey50") +
  geom_line(data=df.dyn.reg.2c, aes(x=BY, y=mu, col=Ocean.Region2), linewidth=1) +
  facet_grid(rows=vars(Ocean.Region2), cols=vars(varnam)) +
  ylim(-1,1)


# Compare to "old" RW fits

# load and process old fits
load(here("output/models/dyn", speciesFlag, "hbm_dyn_2c.RData"), verbose=T)

# Stock-specific dataframe
probs <- c(0.025, 0.05, 0.10, 0.50, 0.90, 0.95, 0.975)
summ <- rstan::summary(dyn.2c, pars = c("gamma", "kappa"), probs = probs)[[1]]
df.dyn.st.2c.old <- data.frame(Stock = data_master$Stock,
                           Ocean.Region2 = data_master$Ocean.Region2,
                           BY = data_master$BY,
                           mu = summ[, "mean"],
                           se = summ[, "se_mean"],
                           lower_10 = summ[, "10%"],
                           upper_90 = summ[ , "90%"],
                           var = str_extract(rownames(summ), "[a-z]+"),
                           varnam = case_when(grepl("^gamma", rownames(summ)) ~ "SST",
                                              grepl("^kappa", rownames(summ)) ~ "Competitors")
)
df.dyn.st.2c.old <- ocean_region_lab(df.dyn.st.2c.old)

# Summarized dataframe (regional-level)
# gamma/kappa are series-specific; no mu output. Summarize stocks instead
if(speciesFlag=="pink"){
  df.dyn.reg.2c.old <- df.dyn.st.2c.old %>%
    mutate(even_odd = ifelse(grepl("Even$", Stock), "Even", "Odd")) %>%
    dplyr::summarize(reg_mean=mean(mu, na.rm=T),
                     n_stk=n_distinct(Stock),
                     lower_10=quantile(mu, 0.1),
                     upper_90=quantile(mu, 0.9),
                     .by=c(Ocean.Region2, BY, varnam, even_odd))
} else {
  df.dyn.reg.2c.old <- dplyr::summarize(df.dyn.st.2c.old,
                                    reg_mean=mean(mu, na.rm=T),
                                    n_stk=n_distinct(Stock),
                                    lower_10=quantile(mu, 0.1),
                                    upper_90=quantile(mu, 0.9),
                                    .by=c(Ocean.Region2, BY, varnam))
}
df.dyn.reg.2c.old <- ddply(df.dyn.reg.2c.old, .(Ocean.Region2), dplyr::filter, n_stk >= max(n_stk)*0.1) # remove years with less than 10% of stocks observed
df.dyn.reg.2c.old <- ocean_region_lab(df.dyn.reg.2c.old)


#Compare new vs old fits on one plot

g <- ggplot(df.dyn.st.2c) + # plot stock trajectories
  geom_line(aes(x=BY, y=mu, group=Stock), col="darkred") +
  geom_line(data=df.dyn.st.2c.old, aes(x=BY, y=mu, group=Stock), col="grey50") +
  facet_grid(rows=vars(Ocean.Region2), cols=vars(varnam)) +
  theme_sleek() +
  ylim(-1,1)
png(filename=here('figures/dyn', paste0('dyn-stk-compare-2026-', speciesFlag, '.png')))
print(g)
dev.off()

g <- ggplot() + # plot region trajectories
  geom_line(data=df.dyn.reg.2c, aes(x=BY, y=mu, col=Ocean.Region2), linewidth=1, linetype="solid") +
  geom_ribbon(data=df.dyn.reg.2c, aes(x=BY, ymin=lower_10, ymax=upper_90, fill=Ocean.Region2), alpha=0.2) +
  facet_grid(rows=vars(Ocean.Region2), cols=vars(varnam)) +
  theme_sleek() +
  ylim(-1,1)
if(speciesFlag=="pink"){
  g <- g + geom_line(data=df.dyn.reg.2c.old, aes(x=BY, y=reg_mean, linetype=even_odd), col="grey25", linewidth=1, alpha=0.5) +
    geom_ribbon(data=df.dyn.reg.2c.old, aes(x=BY, ymin=lower_10, ymax=upper_90, group=even_odd), fill="grey25", alpha=0.2)
} else {
  g <- g + geom_line(data=df.dyn.reg.2c.old, aes(x=BY, y=reg_mean), col="grey25", linewidth=1, alpha=0.5) +
    geom_ribbon(data=df.dyn.reg.2c.old, aes(x=BY, ymin=lower_10, ymax=upper_90), fill="grey25", alpha=0.2)
}
png(filename=here('figures/dyn', paste0('dyn-reg-compare-2026-', speciesFlag, '.png')))
print(g)
dev.off()




