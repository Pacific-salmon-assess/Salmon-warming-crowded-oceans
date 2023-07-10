## HMM Fitting ## 

## 'hmm_test.stan' Accepts a single stock for now
## Fitting to each stock individually

params_out <- c("beta1", "gamma", "alpha", "beta", "log_lik")
nlvl <- 2
# list to store outputs
hmm_out <- list()

#loop over stocks
for ( i in 1:nlevels(sock$Stock) ) {

  
# Filter to one stock
hmm.dat <- dplyr::filter(sock, Stock == levels(sock$Stock)[i]) 

#output array 
hmm_outputs <- array(NA, dim=c(length(params_out), nlvl, nrow(hmm.dat)))
dimnames(hmm_outputs) <- list( params_out, c("State 1", "State 2"), hmm.dat$BY)


# Define data
hmm.stan.dat <- list( 
     N = nrow(hmm.dat),
     R_S = hmm.dat$lnRS,
     S = hmm.dat$S,
     K = nlvl,
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
save(hmm_test, file = paste0("./output/models/hmm_test_", levels(sock$Stock)[i], ".Rdata"))

## Extract parameters of interest
post=rstan::extract(hmm_test)


# dim 6000, 2 (K?) # beta1 is first covar
hmm_outputs[1,,] <- rep(apply(post$beta1, 2, median), times=hmm.stan.dat$N) 
# gamma : dim (6000, 48=N years, 2=K)
hmm_outputs[2,,] <- t( apply(post$gamma, c(2,3), median) )
# alpha : dim (6000, 48=N years, 2=K)
hmm_outputs[3,,] <- t( apply(post$alpha, c(2,3), median) )
# beta : dim (6000, 48=N years, 2=K) # beta is ricker b 
hmm_outputs[4,,] <- t( apply(post$beta, c(2,3), median) )
# log_lik : dim(6000, 48) 
hmm_outputs[5,,] <- rep( apply(post$log_lik, c(2), median), each=nlvl )


hmm_out[[i]] <- hmm_outputs


}

## Data managing and plots 

# Wrangle the list of arrays into a dataframe 
names(hmm_out) <- levels(sock$Stock)
hmm_out_1 <- lapply(hmm_out, plyr::adply, .margins=c(2,3))
hmm_out_2 <- dplyr::bind_rows(hmm_out_1, .id="stock")
names(hmm_out_2)[2:3] <- c("state_K", "BY")
hmm_out_2$BY <- as.numeric(as.character(hmm_out_2$BY)) # this is so dumb

hmm_tidy <- hmm_out_2
## Plot results
# Plot gamma over time (prob of being in each state?)
  # horiz line at 0.5
# Plot (??) the linearized ricker with alpha, beta, and beta1 (stock-specific)
# density plot of covariate terms (beta1) across stocks?

## Gamma (State 1) over time 
pdf("./figures/hmm_test/gamma_timeseries.pdf")
for(i in 1:nlevels(sock$Stock)){
g <- hmm_tidy %>% filter(state_K == "State 1", stock==levels(sock$Stock)[i]) %>%
  ggplot() +
  geom_abline(intercept=0.5, slope=0, colour="grey70", lty="dashed") +
  geom_point(aes(x=BY, y=gamma)) +
  geom_line(aes(x=BY, y=gamma)) + 
  labs(title=levels(sock$Stock)[i]) + theme_minimal() 
print(g)
}
dev.off()

## Density plots of beta1 estimates
hmm_tidy %>% ggplot(aes(x=beta1, fill=state_K)) + 
  geom_histogram() + theme_minimal()


# Plot the linearized ricker with alpha, beta, and beta1 (stock-specific) - look at single_stock function code