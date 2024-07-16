## Fit to "in-development" random walk models using a small subset of data ##
# Intended to be reproducible outside rest of repo.


# -- Load packages
library(ggplot2)
library(dplyr)
library(rstan)
library(here)

# Define a data wrangling function:
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

 
# -- Get (sockeye) data 
all.sock <- read.csv(here::here('data', 'sockeye', 'master_sockeye_brood_table_covar.csv'))
sample.ids <- 10:28 # fit to 19 fraser stocks 
rw.dat <- filter(all.sock, Stock.ID %in% sample.ids) # make rw.dat dataframe
 
 

# -- Model '2a' : grouped RW with stationary coefficients

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
            init=0,
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
rw.2a.dat.r <- data.frame(BY = unique((rw.dat[,c("Ocean.Region2", "BY")]))$BY,
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



## -- Model '5a' - grouped RW with stock-specific deviances

# get stan data
stan.dat.rw5a <- stan.dat.rw2a
stan.dat.rw5a$priors_only = 0

# run
rw_5a <- rstan::stan(file="./stan/hbm_dyn_2c_shared_tvraneff_dag.stan",
                     data = stan.dat.rw5a,
                     warmup = 100,
                     iter = 300,
                     cores = 1,
                     chains = 1,
                     seed = 123,
                     control = list(adapt_delta = 0.9,
                                    max_treedepth = 20))
#save(rw_5a, file = here::here('output', 'models', 'dyn', 'sockeye', 'rw-model-compare', paste0('rw_5a_', Sys.Date(), '.RData')))


# Diagnostics
rstan::check_hmc_diagnostics(rw_5a)
ppc.y <- stan.dat.rw5a$y
ppc.yrep <- rstan::extract(rw_5a)$yrep
#shinystan::launch_shinystan(rw_5a)
mcmc_pairs(rw_5a, pars=c(paste0("d_gamma_it[1,", 1:4, "]"), paste0("sigma_gamma_it[", 10:14, "]")))
ppc_hist(ppc.y, ppc.yrep[1:10,]) + xlim(-100,100)
ppc_dens_overlay(ppc.y, ppc.yrep[1:100,]) + xlim(-100,100)


# Extract gammas/kappas
probs <- c(0.025, 0.975)
summ.g <- rstan::summary(rw_5a, pars = c("gamma", "kappa", "mu_alpha"), probs = probs)[[1]]
summ.g.i <- rstan::summary(rw_5a, pars=c("gamma_it", "kappa_it"), probs=probs)[[1]]
rw.5a.dat.r <- data.frame(BY = unique((rw.dat[,c("Ocean.Region2", "BY")]))$BY,
                        mu = summ.g[, "mean"],
                        lower = summ.g[, "2.5%"],
                        upper = summ.g[ , "97.5%"],
                        var = sub("\\[.*", "", rownames(summ.g)),
                        stock = NA)
rw.5a.dat <- rbind(rw.5a.dat.r, data.frame(BY = rep(min(rw.dat$BY):max(rw.dat$BY),19*2),
                                         mu = summ.g.i[, "mean"],
                                         lower = summ.g.i[, "2.5%"],
                                         upper = summ.g.i[ , "97.5%"],
                                         var = sub("_.*", "", rownames(summ.g.i)) ,
                                         stock = rep(unique(rw.dat$Stock), each=66)
))

# Plot gamma & kappa time series
ggplot(rw.5a.dat[is.na(rw.5a.dat$stock),]) +
  geom_line(aes(x=BY, y=mu), col = "navyblue", linewidth=1) + 
  geom_line(data=rw.5a.dat[!is.na(rw.5a.dat$stock),], aes(x=BY, y=mu, group=stock), col="grey40") + 
  geom_ribbon(aes(x=BY, y=mu, ymin=lower, ymax=upper), fill="navyblue", alpha=0.2) +
  facet_wrap(vars(var)) + 
  ylim(-2,1.2) + 
  theme_minimal() + labs(x="brood year", y="covar coefficient est.")


## -- Multivariate model 

load(here::here('output', 'models', 'dyn', 'sockeye', 'mvrw_fraser.RData'), verbose=T)

# Get data
probs <- c(0.025, 0.975)
summ <- rstan::summary(mvrw.fr, pars=c("gamma_t_m", "kappa_t_m"), probs=probs)$summary
mvrw.out <- data.frame(BY = unique(mvrw.dat$BY),
                       mu = summ[, "mean"],
                       lower = summ[, "2.5%"],
                       upper = summ[, "97.5%"],
                       var = sub("_.*", "", rownames(summ)),
                       stock = NA)



## -- Independent ('original') model
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
ind.rw.out <- ind.rw.out %>% filter(stock %in% rw.dat$Stock) %>% dplyr::summarize(mu = mean(mu, na.rm=T), lower=mean(lower, na.rm=T), upper=mean(upper, na.rm = T), .by=c(BY, var))


## Compare models ------------------------------------------ ## 
df_all <- bind_rows(rw.2a.dat.r, rw.5a.dat.r, mvrw.out, ind.rw.out, .id="model") 
rownames(df_all) <- NULL
df_all <- mutate(df_all, model=case_when(model==1 ~ "shared rw + stat. dev. (2a)",
                                         model==2 ~ "shared rw + time-var dev. (5a)",
                                         model==3 ~ "multivariate",
                                         model==4 ~ "independent rw avg"))

ggplot(df_all) + geom_line(aes(x=BY, y=mu, col=model), linewidth=.75, alpha=0.7) +
  geom_ribbon(aes(x=BY, ymin=lower, ymax=upper, fill=model), alpha=0.1) +
  geom_hline(yintercept=0, col="grey70") +
  facet_wrap(vars(var)) + theme_sleek() 

# Make comparison table

# add alpha
# add 

