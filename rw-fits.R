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
    
    
# Get stan data    

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


## -- Run the model with 1 covariate to test -------------------------------------

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

## -- DG's RW model ------------------------------------

# subset of stocks to test with
rw.dg.dat <- filter(data_dyn, Ocean.Region2 == "BS")

#stock-year index for model
stock_year=expand.grid(unique(rw.dg.dat$Stock),unique(rw.dg.dat$BY))
stock_year= stock_year[order(stock_year[,1]),]
stock_year[,3]=paste(stock_year[,1],stock_year[,2],sep="_")
rw.dg.dat$stock_yr=match(paste(rw.dg.dat$Stock,rw.dg.dat$BY,sep='_'),stock_year[,3])


#util function for spawners design matrix
make_design_matrix=function(x,grp){
  x2=matrix(nrow=length(x),ncol=length(unique(grp)))
  for(i in 1:length(unique(grp))){
    x2[,i]=ifelse(grp==levels(factor(grp))[i],1,0)*x
  }
  return(x2)
}

# Make design matrices
X_s=make_design_matrix(rw.dg.dat$S, grp=rw.dg.dat$Stock)
X_x=make_design_matrix(rw.dg.dat$early_sst_stnd, grp=rw.dg.dat$Stock)


# get smax prior 
smax_prior=rw.dg.dat%>%group_by(Stock) %>%summarize(m=max(S))

# list of data
stan.dat.dg <- list(N=nrow(rw.dg.dat), #number of observations
     L=max(rw.dg.dat$BY)-min(rw.dg.dat$BY)+1, #total years
     J_i=as.numeric(factor(rw.dg.dat$Stock)), #stock index
     J_ii=rw.dg.dat$stock_yr, #stock-year index
     #L_i=rw.dg.dat$broodyear-min(rw.dg.dat$broodyear)+1, #year index
     #Y_i=as.numeric(factor(rw.dg.dat$BY)),
     #start_y=start_n,
     #end_y=end_n,
     #start_t=s_y$s,
     #end_t=s_y$e,
     J=length(unique(rw.dg.dat$Stock)),
     R_S=rw.dg.dat$lnRS,
     #S_X=X_s,
     S=X_s,
     X=rw.dg.dat$early_sst_stnd,
     pSmax_mean=smax_prior$m*0.5,
     pSmax_sig=smax_prior$m*0.5)
  
# Run
rw.dg.bs <- rstan::stan(file = "./stan/rwa_mv_HH_dag.stan",
                     data = stan.dat.dg,
                     #pars = c(pars_dyn_2c, pars.gen.quant),
                     warmup = 1000,
                     iter = 3000, 
                     cores = 4,
                     chains = 4,
                     seed = 123,
                     control = list(adapt_delta = 0.995,
                                    max_treedepth = 20))

rstan::check_hmc_diagnostics(rw.dg.bs)

# actually ran Fraser stks not Bering Sea
save(rw.dg.fraser, file=here(fit.dir, 'rw_dg_fraser.RData'))

# Plot
summ <- rstan::summary(rw.dg.bs, pars=c("g_t"))$summary
dg.rw.dat <- data.frame(mu = summ[,"mean"],
                        lower_2.5 = summ[,"2.5%"],
                        upper_97.5 = summ[, "97.5%"],
                        Stock = rep(unique(rw.dg.dat$Stock), times=stan.dat.dg$L),
                        BY = rep(sort(unique(rw.dg.dat$BY)), each=n_distinct(rw.dg.dat$Stock)),
                        BY_start = rep(dplyr::summarize(rw.dg.dat, min(BY), .by = "Stock")[,2], times=stan.dat.dg$L),
                        BY_end = rep(dplyr::summarize(rw.dg.dat, max(BY), .by = "Stock")[,2], times=stan.dat.dg$L)
                        )
dg.rw.dat <- mutate(dg.rw.dat, from_data = as.factor(if_else(dplyr::between(BY, BY_start, BY_end), 1, 0)))

ggplot() + geom_line(data=dg.rw.dat, aes(BY, mu, group=Stock, colour=from_data))
  


