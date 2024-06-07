## Fit Random Walk time-varying covariate models ## 

# Species
if(speciesFlag=="pink") 
  data_master <- pink else if(speciesFlag=="chum") 
    data_master <- chum else if(speciesFlag=="sockeye")
      data_master <- sock
    
    # Set paths to output locations - dependent on species 
    fig.dir <- here("figures", "dyn", speciesFlag, "hbm_fit") # place to store figures generated in this script
    fit.dir <- here("output", "models", "dyn", speciesFlag) # place to store model fits
    diag.dir <- here("output", "diagnostics", "dyn", speciesFlag) # place to store diagnostics
    
    
    # Make them if they don't exist
    if(!dir.exists(fig.dir))
      dir.create(fig.dir, recursive = T)
    
    if(!dir.exists(fit.dir))
      dir.create(fit.dir, recursive = T)
    
    if(!dir.exists(diag.dir))
      dir.create(diag.dir, recursive = T)

    
## -- 1. SHARED RW model with 2 covariates -------------------------------------
 
# In this model, all data are pooled w/in regions and a single RW series is calculated for each region.

# filter out short timeseries
stk_20 <- info_master$Stock[info_master$n_years >= 20]
data_dyn <- dplyr::filter(data_master, Stock %in% stk_20)
    
stan.dat.2c <- stan_data_dyn(data_dyn, 
                                 var.x2 = "early_sst_stnd",
                                 var.x3 = "np_pinks_sec_stnd",
                                 breakpoint1 = 1989,
                                 breakpoint2 = 2011,
                                 var.region="Ocean.Region2", 
                                 scale.x1 = TRUE,
                                 alpha.group = ifelse(speciesFlag=="sockeye", TRUE, FALSE))
    
# set pars to monitor
pars_dyn_2c <- c("alpha", "beta", "sigma", "phi", "mu_alpha", "sigma_alpha",
                 "gamma", "sigma_gamma", "signal_noise_g",
                 "kappa", "sigma_kappa", "signal_noise_k")

pars.gen.quant <- c("log_lik", "yrep") ## Generated quantities to monitor


# Run MCMC

dyn.grp.2c <- rstan::stan(file = "./stan/hbm_dyn_2c_shared.stan",
                      data = stan.dat.2c,
                      pars = c(pars_dyn_2c, pars.gen.quant),
                      warmup = 1000,
                      iter = 4000, # try more iterations to fix Rhat & Neff
                      cores = 4,
                      chains = 4,
                      seed = 123,
                      control = list(adapt_delta = 0.99,
                                     max_treedepth = 20))
save(dyn.grp.2c, file = here(fit.dir, "hbm_dyn_grp_2c.Rdata"))

# Diagnostics
rstan::check_hmc_diagnostics(dyn.grp.2c)
neff_lowest(dyn.grp.2c, pars = pars_dyn_2c)
rhat_highest(dyn.grp.2c, pars = pars_dyn_2c)

pdf(here(diag.dir, "dyn_grp_2c_diag.pdf"), width = 7, height = 5)
coda_neff(get_neff(dyn.grp.2c, pars = pars_dyn_2c), total_draws(dyn.grp.2c))
coda_rhat(get_rhat(dyn.grp.2c, pars = pars_dyn_2c))
coda_diag(As.mcmc.list(dyn.grp.2c, pars = pars_dyn_2c))
dev.off()

## Posterior predictive checks
plot_post_pc(dyn.grp.2c, stan.dat.2c$y, pdf.path = here(diag.dir, "dyn_grp_2c_yrep.pdf"))

## LOOIC
loo.dyn.grp.2c <- rstan::loo(dyn.grp.2c, cores = 4)
save(loo.dyn.grp.2c, file = here(diag.dir, "loo_dyn_grp_2c.RData"))
sum(pareto_k_values(loo.dyn.grp.2c) > 0.7)
pdf(here(diag.dir, "dyn_grp_2c_loo.pdf"), width = 7, height = 5)
plot(loo.dyn.grp.2c, label_points = TRUE)
dev.off()


### --- Dynamic model: Data

# Stock-specific dataframe
probs <- c(0.025, 0.05, 0.10, 0.50, 0.90, 0.95, 0.975)
summ <- rstan::summary(dyn.grp.2c, pars = c("gamma", "kappa"), probs = probs)[[1]]
df.dyn.reg.2c <- data.frame(Ocean.Region2 = unique((data_dyn[,c("Ocean.Region2", "BY")]))$Ocean.Region2,
                           BY = unique((data_dyn[,c("Ocean.Region2", "BY")]))$BY,
                           mu = summ[, "mean"],
                           se = summ[, "se_mean"],
                           lower_10 = summ[, "10%"],
                           upper_90 = summ[ , "90%"],
                           var = str_extract(rownames(summ), "[a-z]+"),
                           varnam = case_when(grepl("^gamma", rownames(summ)) ~ "SST",
                                              grepl("^kappa", rownames(summ)) ~ "Competitors")
)
df.dyn.reg.2c <- ocean_region_lab(df.dyn.reg.2c)

### --- Dynamic model: Figures 

# Grouped gamma timeseries : 2-covar
g <- ggplot(df.dyn.reg.2c) +
  geom_line(aes(x=BY, y=mu, col=ocean_region_lab), linewidth=1) + 
  geom_ribbon(aes(x=BY, y=mu, ymin=lower_10, ymax=upper_90, fill=ocean_region_lab), alpha=0.2) +
  hline_0() + 
  facet_grid(rows=vars(ocean_region_lab), cols=vars(varnam)) + 
  ylim(c(-1,1)) + 
  scale_colour_manual(values=col.region, aesthetics=c("colour", "fill")) +
  theme_minimal()

pdf(here(fig.dir, "dyn_2c_grp_grouped.pdf"))
print(g)
dev.off()

## Overlay on stock-specific RW estimates
load(here('output', 'models', 'dyn', speciesFlag, "hbm_dyn_2c.RData"), verbose=T)
# Make the summary dataframe
summ.stk <- rstan::summary(dyn.2c, pars = c("gamma", "kappa"), probs = probs)[[1]]
df.dyn.st.2c <- data.frame(Stock = data_master$Stock,
                           Ocean.Region2 = data_master$Ocean.Region2,
                           BY = data_master$BY,
                           mu = summ.stk[, "mean"],
                           se = summ.stk[, "se_mean"],
                           lower_10 = summ.stk[, "10%"],
                           upper_90 = summ.stk[ , "90%"],
                           var = str_extract(rownames(summ.stk), "[a-z]+"),
                           varnam = case_when(grepl("^gamma", rownames(summ.stk)) ~ "SST",
                                              grepl("^kappa", rownames(summ.stk)) ~ "Competitors")
)
df.dyn.st.2c <- ocean_region_lab(df.dyn.st.2c)

pdf(here(fig.dir, "dyn_2c_grp_stkoverlay.pdf"))
g<- ggplot(df.dyn.reg.2c) +
  geom_line(data= df.dyn.st.2c, aes(x=BY, y=mu, group=Stock, col=ocean_region_lab), alpha=0.2) +
  geom_line(aes(x=BY, y=mu, col=ocean_region_lab), linewidth=1) +
  hline_0(colour="grey30", linetype="dashed") + 
  #geom_ribbon(aes(x=BY, y=mu, ymin=lower_10, ymax=upper_90, fill=ocean_region_lab), alpha=0.2) +
  facet_grid(rows=vars(ocean_region_lab), cols=vars(varnam)) + 
  ylim(c(-1,1)) + 
  scale_colour_manual(values=col.region, aesthetics=c("colour", "fill")) +
  theme_minimal()
print(g)
dev.off()


## -- Run the SHARED RW model with 1 covariate to test -------------------------------------

# filter out short timeseries
stk_20 <- info_master$Stock[info_master$n_years >= 20]
data_dyn <- dplyr::filter(data_master, Stock %in% stk_20)

stan.dat.1c <- stan_data_dyn(data_dyn, 
                             var.x2 = "early_sst_stnd",
                             var.x3 = "np_pinks_sec_stnd",
                             breakpoint1 = 1989,
                             breakpoint2 = 2011,
                             var.region="Ocean.Region2", 
                             scale.x1 = TRUE,
                             alpha.group = ifelse(speciesFlag=="sockeye", TRUE, FALSE))

# set pars to monitor
pars_dyn_1c <- c("alpha", "beta", "sigma", "phi", "mu_alpha", "sigma_alpha",
                 "gamma", "sigma_gamma", "signal_noise")

pars.gen.quant <- c("log_lik", "yrep") ## Generated quantities to monitor


# Run MCMC

dyn.grp.1c <- rstan::stan(file = "./stan/hbm_dyn_1c_shared.stan",
                          data = stan.dat.1c,
                          pars = c(pars_dyn_1c, pars.gen.quant),
                          warmup = 1000,
                          iter = 4000, 
                          cores = 4,
                          chains = 4,
                          seed = 123,
                          control = list(adapt_delta = 0.8))
save(dyn.grp.1c, file = here(fit.dir, "hbm_dyn_grp_1c.Rdata"))

# Prior PC
stan.dat.1c.prior <- stan.dat.1c
stan.dat.1c.prior$priors_only <- 1
dyn.grp.1c.priors <- rstan::stan(file = "./stan/hbm_dyn_1c_shared.stan",
                          data = stan.dat.1c.prior,
                          pars = c(pars_dyn_1c, pars.gen.quant),
                          warmup = 1000,
                          iter = 4000, # try more iterations to fix Rhat & Neff
                          cores = 4,
                          chains = 4,
                          seed = 123,
                          control = list(adapt_delta = 0.8))


# Diagnostics
rstan::check_hmc_diagnostics(dyn.grp.1c)
neff_lowest(dyn.grp.1c, pars = pars_dyn_1c)
rhat_highest(dyn.grp.1c, pars = pars_dyn_1c)

pdf(here(diag.dir, "dyn_grp_1c_diag.pdf"), width = 7, height = 5)
coda_neff(get_neff(dyn.grp.1c, pars = pars_dyn_1c), total_draws(dyn.grp.1c))
coda_rhat(get_rhat(dyn.grp.1c, pars = pars_dyn_1c))
coda_diag(As.mcmc.list(dyn.grp.1c, pars = pars_dyn_1c))
dev.off()

## Posterior predictive checks
plot_post_pc(dyn.grp.1c, stan.dat.1c$y, pdf.path = here(diag.dir, "dyn_grp_1c_yrep.pdf"))


## LOOIC
loo.dyn.grp.1c <- rstan::loo(dyn.grp.1c, cores = 4)
save(loo.dyn.grp.1c, file = here(diag.dir, "loo_dyn_grp_1c.RData"))
sum(pareto_k_values(loo.dyn.grp.1c) > 0.7)
pdf(here(diag.dir, "dyn_grp_1c_loo.pdf"), width = 7, height = 5)
plot(loo.dyn.grp.1c, label_points = TRUE)
dev.off()

### --- Dynamic model: Data

# Stock-specific dataframe
probs <- c(0.025, 0.05, 0.10, 0.50, 0.90, 0.95, 0.975)
summ <- rstan::summary(dyn.grp.1c, pars = c("gamma"), probs = probs)[[1]]
df.dyn.reg.1c <- data.frame(Ocean.Region2 = unique((data_dyn[,c("Ocean.Region2", "BY")]))$Ocean.Region2,
                            BY = unique((data_dyn[,c("Ocean.Region2", "BY")]))$BY,
                            mu = summ[, "mean"],
                            se = summ[, "se_mean"],
                            lower_10 = summ[, "10%"],
                            upper_90 = summ[ , "90%"],
                            var = str_extract(rownames(summ), "[a-z]+"),
                            varnam = case_when(grepl("^gamma", rownames(summ)) ~ "SST")
)
df.dyn.reg.1c <- ocean_region_lab(df.dyn.reg.1c)

### --- Dynamic model: Figures 

# Grouped gamma timeseries : 1-covar
g <- ggplot(df.dyn.reg.1c) +
  geom_line(aes(x=BY, y=mu, col=ocean_region_lab), linewidth=1) + 
  geom_ribbon(aes(x=BY, y=mu, ymin=lower_10, ymax=upper_90, fill=ocean_region_lab), alpha=0.2) +
  hline_0() + 
  facet_grid(rows=vars(ocean_region_lab), cols=vars(varnam)) + 
  ylim(c(-3,2)) + 
  scale_colour_manual(values=col.region, aesthetics=c("colour", "fill")) +
  theme_minimal()

pdf(here(fig.dir, "dyn_1c_grp_grouped.pdf"))
print(g)
dev.off()


## -- Run the Shared-sigma RW model with 1 covariate to test -------------------------------------

# filter out short timeseries
stk_20 <- info_master$Stock[info_master$n_years >= 20]
data_dyn <- dplyr::filter(data_master, Stock %in% stk_20)

stan.dat.1c <- stan_data_dyn(data_dyn, 
                             var.x2 = "early_sst_stnd",
                             var.x3 = "np_pinks_sec_stnd",
                             breakpoint1 = 1989,
                             breakpoint2 = 2011,
                             var.region="Ocean.Region2", 
                             scale.x1 = TRUE,
                             alpha.group = ifelse(speciesFlag=="sockeye", TRUE, FALSE))

# set pars to monitor
pars_dyn_1c <- c("alpha", "beta", "sigma", "phi", "mu_alpha", "sigma_alpha",
                 "gamma", "sigma_gamma", "signal_noise")

pars.gen.quant <- c("log_lik", "yrep") ## Generated quantities to monitor


# Run MCMC

dyn.sigma.1c <- rstan::stan(file = "./stan/hbm_dyn_1c_sigma_shared.stan",
                          data = stan.dat.1c,
                          pars = c(pars_dyn_1c, pars.gen.quant),
                          warmup = 1000,
                          iter = 4000, 
                          cores = 4,
                          chains = 4,
                          seed = 123,
                          control = list(adapt_delta = 0.8))
save(dyn.sigma.1c, file = here(fit.dir, "hbm_dyn_sigma_shared_1c.Rdata"))


# Diagnostics
rstan::check_hmc_diagnostics(dyn.sigma.1c)
neff_lowest(dyn.sigma.1c, pars = pars_dyn_1c)
rhat_highest(dyn.sigma.1c, pars = pars_dyn_1c)

pdf(here(diag.dir, "dyn_sigma_1c_diag.pdf"), width = 7, height = 5)
coda_neff(get_neff(dyn.sigma.1c, pars = pars_dyn_1c), total_draws(dyn.sigma.1c))
coda_rhat(get_rhat(dyn.sigma.1c, pars = pars_dyn_1c))
coda_diag(As.mcmc.list(dyn.sigma.1c, pars = pars_dyn_1c))
dev.off()

## Posterior predictive checks
plot_post_pc(dyn.sigma.1c, stan.dat.1c$y, pdf.path = here(diag.dir, "dyn_grp_1c_yrep.pdf"))


## LOOIC
loo.dyn.sigma.1c <- rstan::loo(dyn.sigma.1c, cores = 4)
save(loo.dyn.sigma.1c, file = here(diag.dir, "loo_dyn_sigma_1c.RData"))
sum(pareto_k_values(loo.dyn.sigma.1c) > 0.7)
pdf(here(diag.dir, "dyn_sigma_1c_loo.pdf"), width = 7, height = 5)
plot(loo.dyn.sigma.1c, label_points = TRUE)
dev.off()

### --- Dynamic model: Data

# Stock-specific dataframe
probs <- c(0.025, 0.10, 0.50, 0.90, 0.975)
summ <- rstan::summary(dyn.sigma.1c, pars = c("gamma"), probs = probs)[[1]]
df.dyn.reg.1c <- data.frame(Stock = data_dyn$Stock,
                            Ocean.Region2 = data_dyn$Ocean.Region2,
                            BY = data_dyn$BY,
                            mu = summ[, "mean"],
                            se = summ[, "se_mean"],
                            lower_10 = summ[, "10%"],
                            upper_90 = summ[ , "90%"],
                            var = str_extract(rownames(summ), "[a-z]+"),
                            varnam = case_when(grepl("^gamma", rownames(summ)) ~ "SST")
)
df.dyn.reg.1c <- ocean_region_lab(df.dyn.reg.1c)

### --- Dynamic model: Figures 

# Grouped gamma timeseries : 1-covar
g <- ggplot(df.dyn.reg.1c) +
  geom_line(aes(x=BY, y=mu, col=ocean_region_lab, group=Stock), linewidth=0.5, alpha=0.3) + 
  #geom_ribbon(aes(x=BY, y=mu, ymin=lower_10, ymax=upper_90, fill=ocean_region_lab), alpha=0.2) +
  hline_0() + 
  facet_wrap(vars(ocean_region_lab)) + 
  scale_colour_manual(values=col.region, aesthetics=c("colour", "fill")) +
  theme_minimal()

pdf(here(fig.dir, "dyn_1c_sigma_grouped.pdf"))
print(g)
dev.off()


## -- Run the Shared-sigma RW model with 2 covariates to test -------------------------------------

# filter out short timeseries
stk_20 <- info_master$Stock[info_master$n_years >= 20]
data_dyn <- dplyr::filter(data_master, Stock %in% stk_20)

stan.dat.2c <- stan_data_dyn(data_dyn, 
                             var.x2 = "early_sst_stnd",
                             var.x3 = "np_pinks_sec_stnd",
                             breakpoint1 = 1989,
                             breakpoint2 = 2011,
                             var.region="Ocean.Region2", 
                             scale.x1 = TRUE,
                             alpha.group = ifelse(speciesFlag=="sockeye", TRUE, FALSE))

# set pars to monitor
pars_dyn_2c <- c("alpha", "beta", "sigma", "phi", "mu_alpha", "sigma_alpha",
                 "gamma", "kappa", "sigma_gamma", "sigma_kappa", "signal_noise_g", "signal_noise_k")

pars.gen.quant <- c("log_lik", "yrep") ## Generated quantities to monitor


# Run MCMC

dyn.sigma.2c <- rstan::stan(file = "./stan/hbm_dyn_2c_sigma_shared.stan",
                            data = stan.dat.2c,
                            pars = c(pars_dyn_2c, pars.gen.quant),
                            warmup = 1000,
                            iter = 3000, 
                            cores = 4,
                            chains = 4,
                            seed = 123,
                            control = list(adapt_delta = 0.8))
save(dyn.sigma.2c, file = here(fit.dir, "hbm_dyn_sigma_shared_2c.Rdata"))


# Diagnostics
rstan::check_hmc_diagnostics(dyn.sigma.2c)
neff_lowest(dyn.sigma.2c, pars = pars_dyn_2c)
rhat_highest(dyn.sigma.2c, pars = pars_dyn_2c)

pdf(here(diag.dir, "dyn_sigma_2c_diag.pdf"), width = 7, height = 5)
coda_neff(get_neff(dyn.sigma.2c, pars = pars_dyn_2c), total_draws(dyn.sigma.2c))
coda_rhat(get_rhat(dyn.sigma.2c, pars = pars_dyn_2c))
coda_diag(As.mcmc.list(dyn.sigma.2c, pars = pars_dyn_2c))
dev.off()

## Posterior predictive checks
plot_post_pc(dyn.sigma.2c, stan.dat.2c$y, pdf.path = here(diag.dir, "dyn_grp_2c_yrep.pdf"))


## LOOIC
loo.dyn.sigma.2c <- rstan::loo(dyn.sigma.2c, cores = 4)
save(loo.dyn.sigma.2c, file = here(diag.dir, "loo_dyn_sigma_2c.RData"))
sum(pareto_k_values(loo.dyn.sigma.2c) > 0.7)
pdf(here(diag.dir, "dyn_sigma_2c_loo.pdf"), width = 7, height = 5)
plot(loo.dyn.sigma.2c, label_points = TRUE)
dev.off()

### --- Dynamic model: Data

# Stock-specific dataframe
probs <- c(0.025, 0.10, 0.50, 0.90, 0.975)
summ <- rstan::summary(dyn.sigma.2c, pars = c("gamma","kappa"), probs = probs)[[1]]
df.dyn.st.2c <- data.frame(Stock = data_dyn$Stock,
                           Ocean.Region2 = data_dyn$Ocean.Region2,
                           BY = data_dyn$BY,
                           mu = summ[, "mean"],
                           se = summ[, "se_mean"],
                           lower_10 = summ[, "10%"],
                           upper_90 = summ[ , "90%"],
                           var = str_extract(rownames(summ), "[a-z]+"),
                           varnam = case_when(grepl("^gamma", rownames(summ)) ~ "SST",
                                              grepl("^kappa", rownames(summ)) ~ "Competitors")
)
df.dyn.st.2c <- ocean_region_lab(df.dyn.st.2c)





## -- DG's RW model ------------------------------------

# subset of stocks to test with
mvrw.dat <- filter(data_dyn, Stock.ID %in% 10:28)

#stock-year index for model
stock_year=expand.grid(unique(rw.dg.dat$Stock),unique(rw.dg.dat$BY))
stock_year= stock_year[order(stock_year[,1]),]
stock_year[,3]=paste(stock_year[,1],stock_year[,2],sep="_")
mvrw.dat$stock_yr=match(paste(rw.dg.dat$Stock,rw.dg.dat$BY,sep='_'),stock_year[,3])
# Make design matrices
X_s=make_design_matrix(mvrw.dat$S, grp=mvrw.dat$Stock)
X_x=make_design_matrix(mvrw.dat$early_sst_stnd, grp=mvrw.dat$Stock)
# get smax prior 
smax_prior=mvrw.dat%>%group_by(Stock) %>%summarize(m=max(S))

# list of data
stan.dat <- list(N=nrow(mvrw.dat), #number of observations
     L=max(mvrw.dat$BY)-min(mvrw.dat$BY)+1, #total years
     J_i=as.numeric(factor(mvrw.dat$Stock)), #stock index
     J_ii=mvrw.dat$stock_yr, #stock-year index
     J=length(unique(mvrw.dat$Stock)),
     R_S=mvrw.dat$lnRS,
     S=X_s,
     X=mvrw.dat$early_sst_stnd,
     #X1=mvrw.dat$early_sst_stnd,
     #X2=mvrw.dat$np_pinks_sec_stnd,
     pSmax_mean=smax_prior$m*0.5,
     pSmax_sig=smax_prior$m*0.5)
  
# Run
mvrw.fr.c1.simp <- rstan::stan(file = "./stan/rwa_mv_HH_dag.stan",
                     data = stan.dat,
                     #pars = c(pars_dyn_2c, pars.gen.quant),
                     warmup = 1000,
                     iter = 3000, 
                     cores = 4,
                     chains = 4,
                     seed = 123,
                     control = list(adapt_delta = 0.9,
                                    max_treedepth = 20))

rstan::check_hmc_diagnostics(mvrw.fr)

# actually ran Fraser stks not Bering Sea
save(mvrw.fr, file=here(fit.dir, 'mvrw_fraser.RData'))

# Plot

summ <- rstan::summary(mvrw.fr.c1.simp, pars=c("g_t", "gamma_t_m"))$summary #"k_t", "kappa_t_m"))$summary
mvrw.out <- data.frame(mu_gt = summ[grepl("g_t*", rownames(summ)),"mean"],
                        #mu_kt = summ[grepl("k_t*", rownames(summ)),"mean"],
                        mu_gmean = rep(summ[grepl("gamma_t_m*", rownames(summ)), "mean"], each=stan.dat$J),
                        #mu_kmean = rep(summ[grepl("kappa_t_m*", rownames(summ)), "mean"], each=stan.dat$J),
                        Stock = rep(unique(mvrw.dat$Stock), times=stan.dat$L),
                        BY = rep(sort(unique(mvrw.dat$BY)), each=n_distinct(mvrw.dat$Stock)),
                        BY_start = rep(dplyr::summarize(mvrw.dat, min(BY), .by = "Stock")[,2], times=stan.dat$L),
                        BY_end = rep(dplyr::summarize(mvrw.dat, max(BY), .by = "Stock")[,2], times=stan.dat$L)
                        )
mvrw.out <- mutate(mvrw.out, from_data = as.factor(if_else(dplyr::between(BY, BY_start, BY_end), 1, 0)))

ggplot(data=mvrw.out) + geom_line(aes(BY, mu_gt, group=Stock, colour=from_data)) + geom_line(aes(BY, mu_gmean), col="darkblue", linewidth=1)
ggplot(data=mvrw.out) + geom_line(aes(BY, mu_kt, group=Stock, colour=from_data)) + geom_line(aes(BY, mu_kmean), col="darkblue", linewidth=1)
  
# re-run with prior on that one parameter, see if that resolves most divergences - yes
# get SD or CI around mean
# plot this subgroup again with individual and mean trajectories
# compare to stock-specific trajectories

# Altered Dan's code for pulling out the correlation matrix
corr_lst <- rstan::summary(mvrw.fr, pars="Cor_J")$summary[,"50%"]

create_corr_mat=function(x,J){
  #X - should be a draws matrix from cmdstanr
  #J - number of groups
  mat=matrix(nrow=J,ncol=J)
  for(j in 1:J){
    x_j=x[grepl(paste(',',j,']',sep=''),names(x))]
    mat[,j]=x_j
  }
  return(mat)  
}

corr_ff <- create_corr_mat(corr_lst,J=19)

corrplot::corrplot(corr_ff, diag=FALSE, is.corr=FALSE, method='shade', col = corrplot::COL2('RdYlBu', 10), addCoef.col = 'black', tl.srt=0, tl.cex = 0.9, tl.pos='d')



## ----------- COMPARE all the RW models ------------------ ##

# 1. Independent random walk
load(here(fit.dir, "hbm_dyn_2c_sub.Rdata"), verbose=T)
rw1.sum <- rstan::summary(dyn.2c.sub, pars = c("gamma", "kappa"), probs = c())[[1]]
rw1.df <- data.frame(Stock = data_dyn$Stock,
                           Ocean.Region2 = data_dyn$Ocean.Region2,
                           BY = data_dyn$BY,
                           mu = rw1.sum[, "mean"],
                           varnam = case_when(grepl("^gamma", rownames(rw1.sum)) ~ "SST",
                                              grepl("^kappa", rownames(rw1.sum)) ~ "Competitors"))
rw1.df <- ocean_region_lab(rw1.df)
rw1.df.mu <- rw1.df %>% summarize(reg_mean=mean(mu, na.rm=T), 
                                  .by=c(ocean_region_lab, BY, varnam))
rw1.df %>% filter(ocean_region_lab=="West Coast") %>% ggplot() + geom_line(aes(x=BY, y=mu, group=Stock, col=ocean_region_lab), alpha=0.2) + 
  geom_line(data=subset(rw1.df.mu, ocean_region_lab=="West Coast"), aes(x=BY, y=reg_mean, col=ocean_region_lab), linewidth=1) +
  hline_0(col="grey75") +
  facet_grid(rows=vars(ocean_region_lab), cols=vars(varnam)) +
  scale_colour_manual(values=col.region) + 
  scale_y_continuous(limits=c(-1,.75)) +
  theme_sleek() + theme(legend.position="none") +
  labs(x="Brood Year", y="Effect size")


# 2. Grouped random walk
load(here(fit.dir, "hbm_dyn_grp_2c.Rdata"), verbose=T)
probs <- c(0.025, 0.975)
rw2.sum <- rstan::summary(dyn.grp.2c, pars = c("gamma", "kappa"), probs = probs)[[1]]
rw2.df <- data.frame(Ocean.Region2 = unique((data_master[,c("Ocean.Region2", "BY")]))$Ocean.Region2,
                            BY = unique((data_master[,c("Ocean.Region2", "BY")]))$BY,
                            mu = rw2.sum[, "mean"],
                            lower_2.5 = rw2.sum[, "2.5%"],
                            upper_97.5 = rw2.sum[ , "97.5%"],
                            varnam = case_when(grepl("^gamma", rownames(rw2.sum)) ~ "SST",
                                               grepl("^kappa", rownames(rw2.sum)) ~ "Competitors"))
rw2.df <- ocean_region_lab(rw2.df)
rw2.df %>% filter(ocean_region_lab == "West Coast") %>%  ggplot() + geom_line(aes(x=BY, y=mu, col=ocean_region_lab), linewidth=1) + 
  geom_ribbon(aes(x=BY, ymin=lower_2.5, ymax=upper_97.5, fill=ocean_region_lab), alpha=0.2) +
  hline_0(col="grey75") +
  facet_grid(rows=vars(ocean_region_lab), cols=vars(varnam)) +
  scale_colour_manual(values=col.region, aesthetics = c("colour", "fill")) + 
  coord_cartesian(ylim=c(-1.5,1)) +
  theme_sleek() + theme(legend.position="none") +
  labs(x="Brood Year", y="Effect size")


# 3. Shared Sigma random walk
load(here(fit.dir, "hbm_dyn_sigma_shared_2c.Rdata"), verbose=T)
rw3.sum <- rstan::summary(dyn.sigma.2c, pars = c("gamma","kappa"), probs = c())[[1]]
rw3.df <- data.frame(Stock = data_dyn$Stock,
                           Ocean.Region2 = data_dyn$Ocean.Region2,
                           BY = data_dyn$BY,
                           mu = rw3.sum[, "mean"],
                           varnam = case_when(grepl("^gamma", rownames(rw3.sum)) ~ "SST",
                                              grepl("^kappa", rownames(rw3.sum)) ~ "Competitors"))
rw3.df <- ocean_region_lab(rw3.df)
rw3.df.mu <- rw3.df %>% summarize(reg_mean=mean(mu, na.rm=T), 
                                  .by=c(ocean_region_lab, BY, varnam))
rw3.df %>% filter(ocean_region_lab=="West Coast") %>% ggplot() + geom_line(aes(x=BY, y=mu, group=Stock, col=ocean_region_lab), alpha=0.2) + 
  geom_line(data=subset(rw3.df.mu, ocean_region_lab=="West Coast"), aes(x=BY, y=reg_mean, col=ocean_region_lab), linewidth=1) +
  hline_0(col="grey75") +
  facet_grid(rows=vars(ocean_region_lab), cols=vars(varnam)) +
  scale_colour_manual(values=col.region) + 
  scale_y_continuous(limits=c(-1,.75)) +
  theme_sleek() + theme(legend.position="none") +
  labs(x="Brood Year", y="Effect size")

# 4. Multivariate random walk
load(here(fit.dir, "mvrw_fraser.RData"), verbose=T)
rw4.sum <- rstan::summary(mvrw.fr, pars=c("g_t", "gamma_t_m","k_t", "kappa_t_m"))$summary
rw4.df <- data.frame(mu_gt = rw4.sum[grepl("g_t*", rownames(rw4.sum)),"mean"],
           mu_kt = rw4.sum[grepl("k_t*", rownames(rw4.sum)),"mean"],
           mu_gmean = rep(rw4.sum[grepl("gamma_t_m*", rownames(rw4.sum)), "mean"], each=stan.dat$J),
           mu_kmean = rep(rw4.sum[grepl("kappa_t_m*", rownames(rw4.sum)), "mean"], each=stan.dat$J),
           Stock = rep(unique(mvrw.dat$Stock), times=stan.dat$L),
           BY = rep(sort(unique(mvrw.dat$BY)), each=n_distinct(mvrw.dat$Stock)),
           BY_start = rep(dplyr::summarize(mvrw.dat, min(BY), .by = "Stock")[,2], times=stan.dat$L),
           BY_end = rep(dplyr::summarize(mvrw.dat, max(BY), .by = "Stock")[,2], times=stan.dat$L))

rw4.df <- mutate(rw4.df, from_data = as.factor(if_else(dplyr::between(BY, BY_start, BY_end), 1, 0)))


ggplot(rw4.df) + geom_line(aes(BY, mu_gt, group=Stock, colour=from_data), alpha=0.3) + 
  scale_colour_manual(values=c("black", unname(col.region)[1])) +
  hline_0(col="gray75") + 
  geom_line(aes(BY, mu_gmean), col=unname(col.region)[1], linewidth=1) +
  theme_sleek() + labs(x="Brood Year", y="Effect size", title="SST") +
  theme(legend.position="none",
        plot.title = element_text(hjust=0.5))
ggplot(rw4.df) + geom_line(aes(BY, mu_kt, group=Stock, colour=from_data), alpha=0.3) + 
  scale_colour_manual(values=c("black", unname(col.region)[1])) +
  hline_0(col="gray75") + 
  geom_line(aes(BY, mu_kmean), col=unname(col.region)[1], linewidth=1) +
  theme_sleek() + labs(x="Brood Year", y="Effect size", title="Competitors") +
  theme(legend.position="none",
        plot.title = element_text(hjust=0.5))
