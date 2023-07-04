## HMM Fitting ## 

# 'hmm_test.stan' Accepts a single stock for now

# Choose one stock
hmm.dat <- dplyr::filter(sock, Stock == "Gates") 

# Define data
hmm.stan.dat <- list( N = nrow(hmm.dat),
     R_S = hmm.dat$lnRS,
     S = hmm.dat$S,
     K = 2,
     X1 = hmm.dat$early_sst_stnd, # covariate = SST
     alpha_dirichlet = matrix(c(2,1,1,2), nrow=2) )


# Kept all the specs from hmm.fit script ??
hmm_test <- rstan::stan(file = "./stan/hmm_test.stan",
                     data = hmm.stan.dat,
                     warmup = 1000,
                     iter = 4000,
                     cores = 4,
                     chains = 4,
                     seed = 123,
                     control = list(adapt_delta = 0.90,
                                    max_treedepth = 10))
save(hmm_test, file = "./output/models/hmm_test.Rdata")

post=rstan::extract(hmm_test)
