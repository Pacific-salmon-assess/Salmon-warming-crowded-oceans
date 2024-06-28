## This script compares the remaining Random Walk model options using a small subset of data ##

# -- Get data:
sample.ids <- 10:28 # fit to 19 fraser stocks for now
rw.dat <- sock %>% filter(Stock.ID %in% sample.ids)


# -- Model '2a' : grouped RW with stationary coefficients

# Data
stan.dat.rw2a <- list(N = nrow(rw.dat), # total observations
                      n_series = n_distinct(rw.dat$Stock.ID), # no. of series
                      Ng_groups = as.integer(1), # just one gamma group
                      Na_groups = as.integer(1), # fraser v non fraser
                      a_group = rep(1, n_distinct(rw.dat$Stock.ID)), # all fraser
                      g_group = rep(1, n_distinct(rw.dat$Stock.ID)), # all same region
                      year = as.numeric(as.factor(rw.dat$BY)), # year index
                      y_start = levels_start_end(rw.dat$Stock)$start, # start point of each stk
                      y_end = levels_start_end(rw.dat$Stock)$end, # end point of each stk
                      g_start = as.array(1), # start point of group #needs to be an array to keep dimensions
                      g_end = as.array(max(diff(levels_start_end(rw.dat$Stock)$start))), # end point of group
                      Ng = as.integer(max(diff(levels_start_end(rw.dat$Stock)$start))), # no. of gammas to be estimated (=max ts length)
                      x1 = rw.dat$S,
                      x2 = rw.dat$early_sst_stnd,
                      x3 = rw.dat$np_pinks_sec_stnd,
                      y = rw.dat$lnRS,
                      priors_only = 1) ## change to 0!

# Run
rw_2a <- rstan::stan(file = "./stan/hbm_dyn_2c_shared_dag.stan",
            data = stan.dat.rw2a,
            warmup = 1000,
            iter = 3000, 
            cores = 4,
            chains = 4,
            seed = 123,
            control = list(adapt_delta = 0.9))
save(rw_2a, file = here('output', 'models', 'dyn', 'sockeye', 'rw-model-compare', 'rw_2a.RData'))

prior<-ppc_dens_overlay(y=stan.dat.rw2a$y, yrep=extract(rw_2a, "yrep")[[1]][1:25,])
pdf(here('hannahland', 'model2a_priorpc.pdf'))
print(prior)
dev.off()
ycheck <- stan.dat.rw2a$y


