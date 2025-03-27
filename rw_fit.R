## Fit Random Walk 

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


## Get data for Stan ---------------------------------------

rw.dat <- sock %>% group_by(Ocean.Region2) %>% arrange(Stock, BY, by_group=T) # this line may not be necessary

#stock-year index for model
stk_yr <- ddply(rw.dat, .(factor(Ocean.Region2, levels=unique(rw.dat$Ocean.Region2))), function(x){
            stock_year = expand.grid(unique(x$Stock),unique(x$BY))
            stock_year = stock_year[order(stock_year[,2]),]
            stock_year= stock_year[order(stock_year[,1]),] 
            stock_year[,3] = paste(stock_year[,1],stock_year[,2],sep="_")
            return(stock_year)
          })
rw.dat$stock_yr=match(paste(rw.dat$Stock,rw.dat$BY,sep='_'),stk_yr[,4])


# Make design matrices
X_s=make_design_matrix(rw.dat$S, grp=rw.dat$Stock)
# smax priors
smax_prior=rw.dat%>%group_by(Stock) %>%dplyr::summarize(m=max(S))
# summaries by region
rw.sum <- rw.dat %>% group_by(Ocean.Region2) %>% summarize(L = n_distinct(BY), J = n_distinct(Stock)) %>% arrange(factor(Ocean.Region2, levels=c("WC", "SEAK", "GOA", "BS")))

j_l <- vector()
j_init <- vector()
Lstart <- vector()
Lend <- vector()
pos=1
for(r in 1:4){
 mat <- matrix(pos:(pos+(rw.sum$L*rw.sum$J)[r]-1), nrow=rw.sum$J[r], ncol=rw.sum$L[r], byrow=T)  
 j_l <- c(j_l, as.vector(mat))
 j_init <- c(j_init, as.vector(mat[,1]))
 Lstart <- c(Lstart, seq(from=pos, by=rw.sum$J[r], length.out=rw.sum$L[r]))
 Lend <- c(Lend, seq(from=(pos-1+rw.sum$J[r]), by=rw.sum$J[r], length.out=rw.sum$L[r]))
 pos = max(mat)+1
}

# list of data
stan.dat <- list(N=nrow(rw.dat), #number of observations
                 L= rw.sum$L,
                 J = rw.sum$J,
                 R = 4,
                 L_tot = sum(rw.sum$L),
                 J_tot = sum(rw.sum$J),
                 L_c = cumsum(rw.sum$L),
                 Ng = sum(rw.sum$L*rw.sum$J),
                 J_l = j_l,
                 j_init = j_init,
                 Lstart = Lstart,
                 Lend = Lend,
                 J_i=as.numeric(factor(rw.dat$Stock)), #stock index
                 J_ii=rw.dat$stock_yr, #stock-year index
                 R_S=rw.dat$lnRS,
                 S=X_s,
                 X1=rw.dat$early_sst_stnd,
                 X2=rw.dat$np_pinks_sec_stnd,
                 pSmax_mean=smax_prior$m*0.5,
                 pSmax_sig=smax_prior$m*0.5)



## Run MCMC  -----------------------------------------
mvrw <- rstan::stan(file = "./stan/rwa_mv_2c_sep_Lcorr_Reg.stan",
                    data = stan.dat,
                    warmup = 10,
                    iter = 50, 
                    cores = 4,
                    init = "random",
                    chains = 1,
                    seed = 123,
                    control = list(adapt_delta = 0.9,
                                   max_treedepth = 20))

save(mvrw, file=here::here(fit.dir, paste0('mvrw_2corr_', Sys.Date(), '.RData')))

cmdstanr::cmdstan_model("./stan/rwa_mv_2c_sep_Lcorr_Reg.stan")

## Diagnostics -----------------------------------------

## Check pathology 

rstan::check_hmc_diagnostics(mvrw)

neff_lowest(mvrw, pars = pars_dyn_2c)

rhat_highest(mvrw, pars = pars_dyn_2c)

rstan::get_elapsed_time(mvrw)

## MCMC diagnostics

pdf(here(diag.dir, "dyn_2c_diag.pdf"), width = 7, height = 5)
coda_neff(get_neff(dyn.2c, pars = pars_dyn_2c), total_draws(dyn.2c))
coda_rhat(get_rhat(dyn.2c, pars = pars_dyn_2c))
coda_diag(As.mcmc.list(dyn.2c, pars = pars_dyn_2c))
dev.off()


## Posterior predictive checks
plot_post_pc(mvrw, stan.dat.2c$y, pdf.path = here(diag.dir, "dyn_2c_yrep.pdf"))

## LOOIC
loo.dyn.2c <- rstan::loo(mvrw, cores = 4)

save(loo.dyn.2c, file = here(diag.dir, "loo_dyn_2c.RData"))

sum(pareto_k_values(loo.dyn.2c) > 0.7)

pdf(here(diag.dir, "dyn_2c_loo.pdf"), width = 7, height = 5)
plot(loo.dyn.2c, label_points = TRUE)
dev.off()
