## Fit to "in-development" random walk models using a small subset of data ##
# Intended to be reproducible outside rest of repo.


# -- Load packages
library(dplyr)
library(rstan)
library(here)
library(ggplot2)

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
all.sock <- read.csv(here('data', 'sockeye', 'master_sockeye_brood_table_covar.csv'))
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

smax_priors=rw.dat %>% group_by(Stock.ID) %>%summarise(stock=Stock,maxS=max(S))
smax_priors=distinct(smax_priors,Stock.ID,.keep_all=T)


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
            warmup = 200,
            init=0,
            iter = 800, 
            cores = 4,
            chains = 4,
            seed = 123,
            control = list(adapt_delta = 0.9, 
                           max_treedepth = 20)) # increasing treedepth didn't work?
save(rw_2a, file = here('output', 'models', 'dyn', 'sockeye', 'rw-model-compare', paste0('rw_2a_', Sys.Date(), '.RData')))

# Diagnostics
rstan::check_hmc_diagnostics(rw_2a)

# Extract gammas/kappas
probs <- c(0.025, 0.975)
summ.g <- rstan::summary(rw_2a, pars = c("gamma", "kappa"), probs = probs)[[1]]
summ.g.i <- rstan::summary(rw_2a, pars=c("gamma_i", "kappa_i"), probs=probs)[[1]]
rw.2a.dat <- data.frame(BY = unique((rw.dat[,c("Ocean.Region2", "BY")]))$BY,
                        mu = summ.g[, "mean"],
                        lower = summ.g[, "2.5%"],
                        upper = summ.g[ , "97.5%"],
                        var = sub("\\[.*", "", rownames(summ.g)),
                        stock = NA)
rw.2a.dat <- rbind(rw.2a.dat, data.frame(BY = NA,
                                         mu = summ.g.i[, "mean"],
                                         lower = summ.g.i[, "2.5%"],
                                         upper = summ.g.i[ , "97.5%"],
                                         var = sub("_.*", "", rownames(summ.g.i)),
                                         stock = rep(unique(rw.dat$Stock), 2)
                                         ))

# Plot gamma & kappa time series
ggplot(rw.2a.dat) +
  geom_line(aes(x=BY, y=mu), col = "navyblue", linewidth=1) + 
  geom_ribbon(aes(x=BY, y=mu, ymin=lower, ymax=upper), fill="lightblue4", alpha=0.2) +
  #geom_hline(data=rw.2a.dat[!is.na(rw.2a.dat$Stock)], aes(yintercept=mu)) + 
  facet_wrap(vars(var)) + 
  #ylim(-5,5) + 
  theme_minimal() + labs(x="brood year", y="covar coefficient est.")

