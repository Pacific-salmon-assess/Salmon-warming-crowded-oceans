## Fit to "in-development" random walk models using a small subset of data ##
# Intended to be reproducible outside rest of repo.


# -- Load packages
library(ggplot2)
library(dplyr)
library(rstan)
library(here)

rw.compare.dir <- here::here('output', 'models', 'dyn', 'sockeye', 'rw-model-compare')

# -- Get (sockeye) data 
all.sock <- read.csv(here::here('data', 'sockeye', 'master_sockeye_brood_table_covar.csv'))
sample.ids <- c(67:68,72:78) #c(58:66, 69, 79, 80) # fit to BS stks #10:28 # fit to 19 fraser stocks 
rw.dat <- all.sock %>% filter(Stock.ID %in% sample.ids) # make rw.dat dataframe 
n_stk_yr <- rw.dat %>% group_by(BY) %>% dplyr::summarize(n=n())


# Define data wrangling functions: ####

 levels_start_end <- function(factor) {
      ## Find start and end point of levels in a factor
      ## factor = factor (vector)
      n <- as.numeric(factor)
      internal.end <- which(diff(n) != 0)
      internal.start <- internal.end + 1
      end <- c(internal.end, length(n))
      start <- c(1, internal.start)
      list(start = start, end = end)
    }
    if(FALSE) {
      f <- as.factor(c("a", "a", "a", "b", "b", "c", "c", "c"))
      levels_start_end(f)
      f <- as.factor(c("a", "a", "a", "b", "c", "c", "c"))
      levels_start_end(f)
      levels_start_end(sock$Stock)
    }

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
 

#### Model '2a' : grouped RW with stationary coefficients ####

# Get stan data
start.end <- levels_start_end(as.factor(rw.dat$Stock))
start.end.grp.lst <- lapply(split(rw.dat, rw.dat[["Ocean.Region2"]]),
                            function(x) min(x$BY):max(x$BY))
start.end.grp.vec <- unlist(lapply(seq_along(start.end.grp.lst),
                                   function(x) rep(x, length(start.end.grp.lst[[x]]))))
start.end.grp <- levels_start_end(as.factor(start.end.grp.vec))

smax_priors <- rw.dat %>% group_by(Stock.ID) %>% dplyr::summarize(stock=unique(Stock),maxS=max(S))

stan.dat.rw2a <- list(N = nrow(rw.dat), # total observations
                      n_series = n_distinct(rw.dat$Stock.ID), # no. of series
                      n_years = max(rw.dat$BY) - min(rw.dat$BY) +1, # max span of years
                      Ng_groups = n_distinct(rw.dat$Ocean.Region2), # just one gamma group in training
                      Na_groups = n_distinct(rw.dat$Ocean.Region2), # fraser v non fraser
                      a_group = rep(1:n_distinct(rw.dat$Ocean.Region2), 
                                    n_distinct(rw.dat$Stock.ID)), # all fraser
                      g_group = rep(1:n_distinct(rw.dat$Ocean.Region2), 
                                    n_distinct(rw.dat$Stock.ID)), # all same region
                      year = as.numeric(as.factor(rw.dat$BY)), # year index
                      y_start = start.end$start, # start point of each stk
                      y_end = start.end$end, # end point of each stk
                      g_start = as.array(start.end.grp$start), # start point of group #needs to be an array if 1 
                      g_end = as.array(start.end.grp$end), # end point of group
                      Ng = as.integer(max(diff(start.end$start))), # no. of gammas to be estimated (=max ts length)
                      x1 = rw.dat$S,
                      x2 = rw.dat$early_sst_stnd,
                      x3 = rw.dat$np_pinks_sec_stnd,
                      y = rw.dat$lnRS,
                      priors_only = 0, 
                      pSmax_mean = 0.5*smax_priors$maxS,
                      pSmax_sig = 0.5*smax_priors$maxS)

# Alternative way to get data 
#stan.dat.rw2a <- stan_data_dyn(sock, var.x2 = "early_sst_stnd", var.x3 = "np_pinks_sec_stnd",breakpoint1 = 1989,breakpoint2 = 2011, var.region="Ocean.Region2", scale.x1 = TRUE, alpha.group = FALSE)

# Run
rw_2a <- rstan::stan(file = "./stan/hbm_dyn_2c_shared_dag.stan",
            data = stan.dat.rw2a,
            warmup = 1000,
            iter = 3000, 
            cores = 4,
            chains = 4,
            seed = 123,
            control = list(adapt_delta = 0.9, 
                           max_treedepth = 20))
save(rw_2a, file = here::here('output', 'models', 'dyn', 'sockeye', 'rw-model-compare', paste0('rw_2a_', Sys.Date(), '.RData')))

# Diagnostics
rstan::check_hmc_diagnostics(rw_2a)

# Extract gammas/kappas
probs <- c(0.025, 0.975)
summ.g <- rstan::summary(rw_2a, pars = c("gamma", "kappa"), probs = probs)[[1]]
summ.g.i <- rstan::summary(rw_2a, pars=c("gamma_i", "kappa_i"), probs=probs)[[1]]
rw.2a.dat.r <- data.frame(BY = sort(unique((rw.dat[,c("BY")]))),
                        mu = summ.g[, "mean"],
                        lower = summ.g[, "2.5%"],
                        upper = summ.g[ , "97.5%"],
                        var = sub("\\[.*", "", rownames(summ.g)),
                        stock = NA)
rw.2a.dat <- rbind(rw.2a.dat.r, data.frame(BY = NA,
                                         mu = summ.g.i[, "mean"],
                                         lower = summ.g.i[, "2.5%"],
                                         upper = summ.g.i[ , "97.5%"],
                                         var = sub("_.*", "", rownames(summ.g.i)),
                                         stock = rep(unique(rw.dat$Stock), 2)
                                         ))

# Plot gamma & kappa time series
ggplot(rw.2a.dat) +
  geom_line(aes(x=BY, y=mu), col = "navyblue", linewidth=1) + 
  geom_ribbon(aes(x=BY, y=mu, ymin=lower, ymax=upper), fill="navyblue", alpha=0.2) +
  geom_hline(data=rw.2a.dat[!is.na(rw.2a.dat$stock),], aes(yintercept=mu), col="grey40") + 
  facet_wrap(vars(var)) + 
  ylim(-2,1.2) + 
  theme_minimal() + labs(x="brood year", y="covar coefficient est.")



#### Model '5a' - grouped RW with stock-specific deviances ####

# get stan data
stan.dat.rw5a <- stan.dat.rw2a
stan.dat.rw5a$priors_only = 0

# run
rw_5a <- rstan::stan(file="./stan/hbm_dyn_2c_shared_tvraneff_dag.stan",
                     data = stan.dat.rw5a,
                     warmup = 1000,
                     iter = 3000,
                     cores = 4,
                     chains = 4,
                     seed = 123,
                     control = list(adapt_delta = 0.99,
                                    max_treedepth = 20))
save(rw_5a, file = here::here('output', 'models', 'dyn', 'sockeye', 'rw-model-compare', paste0('rw_5a_', Sys.Date(), '.RData')))


# Diagnostics
rstan::check_hmc_diagnostics(rw_5a)
#ppc.y <- stan.dat.rw5a$y
#ppc.yrep <- rstan::extract(rw_5a)$yrep
#shinystan::launch_shinystan(rw_5a)
#mcmc_pairs(rw_5a, pars=c(paste0("d_gamma_it[1,", 1:4, "]"), paste0("sigma_gamma_it[", 10:14, "]")))
#ppc_hist(ppc.y, ppc.yrep[1:10,]) + xlim(-100,100)
#ppc_dens_overlay(ppc.y, ppc.yrep[1:100,]) + xlim(-100,100)


# Extract gammas/kappas
probs <- c(0.025, 0.975)
summ.g <- rstan::summary(rw_5a, pars = c("gamma", "kappa"), probs = probs)[[1]]
summ.g.i <- rstan::summary(rw_5a, pars=c("gamma_it", "kappa_it"), probs=probs)[[1]]
rw.5a.dat.r <- data.frame(BY = sort(unique((rw.dat[,c("BY")]))),
                        mu = summ.g[, "mean"],
                        lower = summ.g[, "2.5%"],
                        upper = summ.g[ , "97.5%"],
                        var = sub("\\[.*", "", rownames(summ.g)),
                        stock = NA)
rw.5a.dat <- rbind(rw.5a.dat.r, data.frame(BY = rep(min(rw.dat$BY):max(rw.dat$BY),stan.dat.rw5a$n_series*2),
                                         mu = summ.g.i[, "mean"],
                                         lower = summ.g.i[, "2.5%"],
                                         upper = summ.g.i[ , "97.5%"],
                                         var = sub("_.*", "", rownames(summ.g.i)) ,
                                         stock = rep(unique(rw.dat$Stock), each=stan.dat.rw5a$n_years)
))

# Plot gamma & kappa time series
ggplot(rw.5a.dat[is.na(rw.5a.dat$stock),]) +
  geom_line(aes(x=BY, y=mu), col = "navyblue", linewidth=1) + 
  #geom_line(data=rw.5a.dat[!is.na(rw.5a.dat$stock),], aes(x=BY, y=mu, group=stock), col="grey40") + 
  geom_ribbon(aes(x=BY, y=mu, ymin=lower, ymax=upper), fill="navyblue", alpha=0.2) +
  facet_wrap(vars(var)) + 
  ylim(-2,1.2) + 
  theme_minimal() + labs(x="brood year", y="covar coefficient est.")


#### Multivariate model ####

# Get data for stan

#stock-year index for model
stock_year=expand.grid(unique(rw.dat$Stock),unique(rw.dat$BY))
stock_year = stock_year[order(stock_year[,2]),] # added this to maybe resolve stuff?
stock_year= stock_year[order(stock_year[,1]),]
stock_year[,3]=paste(stock_year[,1],stock_year[,2],sep="_")
rw.dat$stock_yr=match(paste(rw.dat$Stock,rw.dat$BY,sep='_'),stock_year[,3])
# Make design matrices
X_s=make_design_matrix(rw.dat$S, grp=rw.dat$Stock)
X_x=make_design_matrix(rw.dat$early_sst_stnd, grp=rw.dat$Stock)
# smax priors
smax_prior=rw.dat%>%group_by(Stock) %>%dplyr::summarize(m=max(S))

# list of data
stan.dat <- list(N=nrow(rw.dat), #number of observations
                 L=max(rw.dat$BY)-min(rw.dat$BY)+1, #total years
                 J_i=as.numeric(factor(rw.dat$Stock)), #stock index
                 J_ii=rw.dat$stock_yr, #stock-year index
                 J=length(unique(rw.dat$Stock)),
                 R_S=rw.dat$lnRS,
                 S=X_s,
                 X1=rw.dat$early_sst_stnd,
                 X2=rw.dat$np_pinks_sec_stnd,
                 pSmax_mean=smax_prior$m*0.5,
                 pSmax_sig=smax_prior$m*0.5)

# RUn multivariate model
mvrw <- rstan::stan(file = "./stan/rwa_mv_2c_sep_Lcorr.stan",
                               data = stan.dat,
                               warmup = 1000,
                               iter = 3000, 
                               cores = 4,
                               chains = 4,
                               seed = 123,
                               control = list(adapt_delta = 0.9,
                                              max_treedepth = 20))

save(mvrw, file=here::here(rw.compare.dir, paste0('mvrw_2corr_', Sys.Date(), '.RData')))

# OR load an old version:
#load(here::here('output', 'models', 'dyn', 'sockeye', 'rw-model-compare', 'mvrw_2corr_2024-07-31.RData'), verbose=T)

# Get output
probs <- c(0.025, 0.975)
summ <- rstan::summary(mvrw, pars=c("gamma_t_m", "kappa_t_m"), probs=probs)$summary
mvrw.out <- data.frame(BY = sort(unique(rw.dat$BY)),
                       mu = summ[, "mean"],
                       lower = summ[, "2.5%"],
                       upper = summ[, "97.5%"],
                       var = sub("_.*", "", rownames(summ)),
                       stock = NA
                       )
summ.i <- rstan::summary(mvrw, pars=c("g_t", "k_t"), probs=probs)$summary
mvrw.out.i <- data.frame(mu = summ.i[,"mean"],
                         Stock = rep(unique(rw.dat$Stock), times=stan.dat$L),
                         BY = rep(sort(unique(rw.dat$BY)), each=n_distinct(rw.dat$Stock)),
                         var = rep(c("gamma", "kappa"), each=stan.dat$L*stan.dat$J),
                         BY_start = rep(dplyr::summarize(rw.dat, min(BY), .by = "Stock")[,2], times=stan.dat$L),
                         BY_end = rep(dplyr::summarize(rw.dat, max(BY), .by = "Stock")[,2], times=stan.dat$L))


#Plot results 
ggplot() + geom_line(data= mvrw.out, aes(x=BY, y=mu), linewidth=1, col="navyblue") +
  hline_0(col="grey50") +
  geom_ribbon(data=mvrw.out, aes(x=BY, ymin=lower, ymax=upper, y=mu), fill="navyblue", alpha=0.2) +
  geom_line(data=mvrw.out.i, aes(x=BY, y=mu, group=Stock), col="grey20") +
  theme_sleek() +
  facet_wrap(vars(var))

#altered Dan's code for pulling out the correlation matrix
corr_lst <- list(gamma = rstan::summary(mvrw, pars="Cor_J_g")$summary[,"50%"],
                 kappa = rstan::summary(mvrw, pars="Cor_J_k")$summary[,"50%"])
corr_ff <- list(gamma = create_corr_mat(corr_lst$gamma, J=stan.dat$J),
                kappa = create_corr_mat(corr_lst$kappa, J=stan.dat$J))
corr_ff <- lapply(corr_ff, function(cor){for(i in 1:nrow(cor)){cor[i,i] <- 0}
                              return(cor)})
png(here::here(rw.compare.dir, paste0('Lcorr_', Sys.Date(), '.png')), width=1300, height=739, res=195)
par(mfrow=c(1,2))
corrplot::corrplot(corr_ff$gamma, diag=FALSE, method='shade', is.corr=F, 
                   col.lim = c(-0.2,0.2), title="gamma",
                   addCoef.col = 'black', tl.srt=0, tl.cex = 0.9, type='upper')
corrplot::corrplot(corr_ff$kappa, diag=FALSE, method='shade', is.corr=F, 
                   col.lim = c(-0.2,0.2), title="kappa",
                   addCoef.col = 'black', tl.srt=0, tl.cex = 0.9, type='upper')
dev.off()


#### Independent ('original') model ####
load(here::here('output', 'models', 'dyn', 'sockeye', 'hbm_dyn_2c_sub.RData'), verbose=T)

# Get data
probs <- c(0.025, 0.975)
summ <- rstan::summary(dyn.2c.sub, pars = c("gamma", "kappa"), probs = probs)[[1]]
ind.rw.out <- data.frame(BY = sock$BY,
                           mu = summ[, "mean"],
                           lower = summ[, "2.5%"],
                           upper = summ[ , "97.5%"],
                           var = sub("\\[.*", "", rownames(summ)),
                           stock = sock$Stock
                           )
ind.rw.out <- ind.rw.out %>% filter(stock %in% rw.dat$Stock, BY %in% rw.dat$BY) %>% dplyr::summarize(mu = mean(mu, na.rm=T), lower=mean(lower, na.rm=T), upper=mean(upper, na.rm = T), .by=c(BY, var))


#### Assess model fits

  # residual plots over year
  # posterior pc - see plot_post_pc in functions - not all models have the right outputs
  # would need to calculate ypred manually for dan's models

#### Compare models ------------------------------------------ #### 
df_all <- bind_rows(rw.2a.dat.r, rw.5a.dat.r, mvrw.out, ind.rw.out, .id="model_id")
rownames(df_all) <- NULL
df_all <- mutate(df_all, model=case_when(model_id==1 ~ "shared rw + stat. dev. (2a)",
                                         model_id==2 ~ "shared rw + time-var dev. (5a)",
                                         model_id==3 ~ "multivariate (4)",
                                         model_id==4 ~ "independent rw avg (1)"))
save(df_all, file=here::here(rw.compare.dir, 'goa-model-outputs.RData'))

png(here::here(rw.compare.dir, paste0('rw_compare_', Sys.Date(), '.png')), width=1300, height=739, res=195)
ggplot(df_all) + geom_line(aes(x=BY, y=mu, col=model), linewidth=.75, alpha=0.7) +
  #geom_line(data=n_stk_yr, aes(x=BY, y=0.6, linewidth=n), col="grey", alpha=0.5) +
  scale_linewidth(range=c(0.5,2)) +
  geom_ribbon(aes(x=BY, ymin=lower, ymax=upper, fill=model), alpha=0.1) +
  geom_hline(yintercept=0, col="grey70") + lims(x=c(1950,2020)) +
  facet_wrap(vars(var), labeller=as_labeller(c("gamma"="SST", "kappa"="Competitors"))) + theme_sleek() + labs(x="Brood Year", y="Mean Effect", col="Model", fill="Model")
dev.off()

# Make comparison table

comp.tbl <- data.frame(delta_gamma = NA,
                       delta_kappa = NA,
                       max_gamma = NA,
                       min_gamma = NA,
                       max_kappa = NA,
                       min_kappa = NA,
                       runtime = NA
                       ) 
comp.tbl[1:4,] <- NA
rownames(comp.tbl) <- unique(df_all$model)


# runtime 
rt <- list(rw_2a = get_elapsed_time(rw_2a),
           rw_5a = get_elapsed_time(rw_5a),
           mvrw = get_elapsed_time(mvrw),
           indep = get_elapsed_time(dyn.2c.sub)
)
comp.tbl[,"runtime"] <- c((max(rt$rw_2a[,"warmup"])+max(rt$rw_2a[,"sample"]))/60,
                          (max(rt$rw_5a[,"warmup"])+max(rt$rw_5a[,"sample"]))/60,
                          (max(rt$mvrw[,"warmup"])+max(rt$mvrw[,"sample"]))/60, # this isn't accurate for fraser run
                          (max(rt$indep[,"warmup"])+max(rt$indep[,"sample"]))/(60*4) # guess at ind runtime with 19 stks
                          )


# delta gamma/kappa
mean.start <- df_all %>% group_by(var, model_id) %>% slice_head(n=5) %>% dplyr::summarize(mean_start = mean(mu))
mean.end <- df_all %>% group_by(var, model_id) %>% slice_tail(n=5) %>% dplyr::summarize(mean_end = mean(mu))
comp.tbl[,"delta_gamma"] <- unname(mean.end[mean.end$var=="gamma","mean_end"] - mean.start[mean.start$var=="gamma","mean_start"])
comp.tbl[,"delta_kappa"] <- unname(mean.end[mean.end$var=="kappa","mean_end"] - mean.start[mean.start$var=="kappa","mean_start"])

# max and min values
max.min <- df_all %>% group_by(var, model_id) %>% dplyr::summarize(max = max(mu), min = min(mu))
comp.tbl[,c("max_gamma", "min_gamma")] <- max.min[max.min$var=="gamma",c("max", "min")]
comp.tbl[,c("max_kappa", "min_kappa")] <- max.min[max.min$var=="kappa",c("max", "min")]

# add notes
# nature of global trends, nature of stock-specific trends, alpha parameterization, gamma parameterization...
comp.tbl$stk_est <- c("stationary deviations from regional mean", "time-varying deviations from regional mean", "estimated", "estimated")
comp.tbl$reg_est <- c("estimated",  "estimated", "arithmetic mean", "arithmetic mean")

# transpose the table
comp.tbl.t <- t(comp.tbl)
write.csv(comp.tbl.t, file=here::here(rw.compare.dir, 'goa-run-comp.csv'))

