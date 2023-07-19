## HMM Fitting ## 
## Fitting to each stock individually


if(!dir.exists("./output/models/hmm-ss/")) dir.create("./figures/stat/pub/")
if(!dir.exists("./figures/hmm-ss/")) dir.create("./figures/stat/pub/")


# Parameters of interest
params_out <- c("beta1", "gamma", "alpha", "beta", "log_lik")



## hmm_ss1 [Single Stock version 1] : 2 states, original transition matrix (2,1,1,2)
## --------------------------------------------------------------------------------

# list to store outputs
hmm_out_ss1 <- list()

#loop over all stocks
for ( i in 1:nlevels(sock$Stock) ) {
  
  # Number of states
  nlvl = 2 
  
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
  hmm_ss1 <- rstan::stan(file = "./stan/hmm_test.stan",
                       data = hmm.stan.dat,
                       warmup = 1000,
                       iter = 4000,
                       cores = 4,
                       chains = 4,
                       seed = 123,
                       control = list(adapt_delta = 0.90,
                                      max_treedepth = 10))
  save(hmm_ss1, file = paste0("./output/models/hmm-ss/hmm_ss1_", levels(sock$Stock)[i], ".Rdata"))
  
  ## Extract parameters of interest
  post=rstan::extract(hmm_ss1)
  
  
  # beta1 : dim (6000, 2=K)
  hmm_outputs[1,,] <- rep(apply(post$beta1, 2, median), times=hmm.stan.dat$N) 
  # gamma : dim (6000, N years, 2=K)
  hmm_outputs[2,,] <- t( apply(post$gamma, c(2,3), median) )
  # alpha : dim (6000, N years, 2=K)
  hmm_outputs[3,,] <- t( apply(post$alpha, c(2,3), median) )
  # beta : dim (6000, N years, 2=K) # beta is ricker b 
  hmm_outputs[4,,] <- t( apply(post$beta, c(2,3), median) )
  # log_lik : dim(6000, N years) 
  hmm_outputs[5,,] <- rep( apply(post$log_lik, c(2), median), each=nlvl )
  
  
  hmm_out_ss1[[i]] <- hmm_outputs


}

## Visualizations

# Wrangle the list of arrays into a dataframe 
names(hmm_out_ss1) <- levels(sock$Stock)
hmm_out <- lapply(hmm_out_ss1, plyr::adply, .margins=c(2,3))
hmm_out <- dplyr::bind_rows(hmm_out, .id="stock")
names(hmm_out)[2:3] <- c("state_K", "BY")
hmm_out$BY <- as.numeric(as.character(hmm_out$BY)) # this is so dumb
hmm_tidy_ss1 <- hmm_out

## Plot gamma (state 1) over time
pdf("./figures/hmm-ss/ss1_gamma_timeseries.pdf")
for(i in 1:nlevels(sock$Stock)){
g <- hmm_tidy_ss1 %>% filter(state_K == "State 1", stock==levels(sock$Stock)[i]) %>%
  ggplot() +
  geom_abline(intercept=0.5, slope=0, colour="grey70", lty="dashed") +
  geom_point(aes(x=BY, y=gamma)) +
  geom_line(aes(x=BY, y=gamma)) + 
  labs(title=levels(sock$Stock)[i]) + theme_minimal() 
print(g)
}
dev.off()


## Plot the expected beta1 (covariate term)
# Load dynamic model outputs
load("./output/models/dyn/hbm5_sst.Rdata", verbose=T)
rw_dat <- rstan::extract(hbm5.sst)
rw_covar <- apply(rw_dat$gamma, 2, median )
# Need to run dyn model with same dataset to continue

pdf("./figures/hmm-ss/ss1_expb1_timeseries.pdf")
for(i in 1:nlevels(sock$Stock)){
  covar_ts <- hmm_tidy_ss1 %>% filter(stock==levels(sock$Stock)[i]) %>% 
    mutate(gamma_prod = beta1*gamma) %>% dplyr::summarize(covar = sum(gamma_prod), .by=BY)
  
  g <- ggplot(covar_ts) +
    geom_point(aes(x=BY, y=covar)) +
    geom_line(aes(x=BY, y=covar)) + 
    labs(title=levels(sock$Stock)[i]) + theme_minimal() 
  print(g)
}
dev.off()


## hmm_ss2 [Single Stock version 2] : 3 states, 2,1,1 transition matrix
## --------------------------------------------------------------------------------

# list to store outputs
hmm_out_ss2 <- list()

#loop over all stocks
for ( i in 1:nlevels(sock$Stock) ) {
  
  # Number of states
  nlvl = 3 
  
  # Filter to one stock
  hmm.dat <- dplyr::filter(sock, Stock == levels(sock$Stock)[i]) 
  
  #output array 
  hmm_outputs <- array(NA, dim=c(length(params_out), nlvl, nrow(hmm.dat)))
  dimnames(hmm_outputs) <- list( params_out, c("State 1", "State 2", "State 3"), hmm.dat$BY)
  
  
  # Define data
  hmm.stan.dat <- list( 
    N = nrow(hmm.dat),
    R_S = hmm.dat$lnRS,
    S = hmm.dat$S,
    K = nlvl,
    X1 = hmm.dat$early_sst_stnd, # covariate = SST
    alpha_dirichlet = matrix(c(2,1,1,1,2,1,1,1,2), nrow=3) )
  
  
  # Kept all the specs from hmm.fit script ??
  hmm_ss2 <- rstan::stan(file = "./stan/hmm_test.stan",
                          data = hmm.stan.dat,
                          warmup = 1000,
                          iter = 4000,
                          cores = 4,
                          chains = 4,
                          seed = 123,
                          control = list(adapt_delta = 0.90,
                                         max_treedepth = 10))
  save(hmm_ss2, file = paste0("./output/models/hmm-ss/hmm_ss2_", levels(sock$Stock)[i], ".Rdata"))
  
  ## Extract parameters of interest
  post=rstan::extract(hmm_ss2)
  
  
  # beta1 : dim (6000, 2=K)
  hmm_outputs[1,,] <- rep(apply(post$beta1, 2, median), times=hmm.stan.dat$N) 
  # gamma : dim (6000, N years, 2=K)
  hmm_outputs[2,,] <- t( apply(post$gamma, c(2,3), median) )
  # alpha : dim (6000, N years, 2=K)
  hmm_outputs[3,,] <- t( apply(post$alpha, c(2,3), median) )
  # beta : dim (6000, N years, 2=K) # beta is ricker b 
  hmm_outputs[4,,] <- t( apply(post$beta, c(2,3), median) )
  # log_lik : dim(6000, N years) 
  hmm_outputs[5,,] <- rep( apply(post$log_lik, c(2), median), each=nlvl )
  
  
  hmm_out_ss2[[i]] <- hmm_outputs
  
  
}

## Visualizations

# Wrangle the list of arrays into a dataframe 
names(hmm_out_ss2) <- levels(sock$Stock)
hmm_out <- lapply(hmm_out_ss2, plyr::adply, .margins=c(2,3))
hmm_out <- dplyr::bind_rows(hmm_out, .id="stock")
names(hmm_out)[2:3] <- c("state_K", "BY")
hmm_out$BY <- as.numeric(as.character(hmm_out$BY)) # this is so dumb
hmm_tidy_ss2 <- hmm_out

## Plot gamma (state 1) over time
pdf("./figures/hmm-ss/ss2_gamma_timeseries.pdf")
for(i in 1:nlevels(sock$Stock)){
  dat_temp <- filter(hmm_tidy_ss2, stock==levels(sock$Stock)[i])
  g <- ggplot() +
    geom_abline(intercept=0.5, slope=0, colour="grey70", lty="dashed") +
    geom_point(data=dat_temp[which(dat_temp$state_K == "State 1"),], aes(x=BY, y=gamma), col="black") +
    geom_point(data=dat_temp[which(dat_temp$state_K == "State 2"),], aes(x=BY, y=gamma), col="blue") +
    geom_point(data=dat_temp[which(dat_temp$state_K == "State 3"),], aes(x=BY, y=gamma), col="lightblue") +
    geom_line(data=dat_temp[which(dat_temp$state_K == "State 1"),], aes(x=BY, y=gamma), col="black") +
    geom_line(data=dat_temp[which(dat_temp$state_K == "State 2"),], aes(x=BY, y=gamma), col="blue") +
    geom_line(data=dat_temp[which(dat_temp$state_K == "State 3"),], aes(x=BY, y=gamma), col="lightblue") +
    labs(title=levels(sock$Stock)[i]) + theme_minimal() 
  print(g)
}
dev.off()



## hmm_ss3 [Single Stock version 3] : 2 states, conservative transition matrix (4,1,1,4)
## --------------------------------------------------------------------------------

# list to store outputs
hmm_out_ss3 <- list()

#loop over all stocks
for ( i in 1:nlevels(sock$Stock) ) {
  
  # Number of states
  nlvl = 2 
  
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
    alpha_dirichlet = matrix(c(4,1,1,4), nrow=2) )
  
  
  # Kept all the specs from hmm.fit script ??
  hmm_ss3 <- rstan::stan(file = "./stan/hmm_test.stan",
                         data = hmm.stan.dat,
                         warmup = 1000,
                         iter = 4000,
                         cores = 4,
                         chains = 4,
                         seed = 123,
                         control = list(adapt_delta = 0.90,
                                        max_treedepth = 10))
  save(hmm_ss3, file = paste0("./output/models/hmm-ss/hmm_ss3_", levels(sock$Stock)[i], ".Rdata"))
  
  ## Extract parameters of interest
  post=rstan::extract(hmm_ss3)
  
  
  # beta1 : dim (6000, 2=K)
  hmm_outputs[1,,] <- rep(apply(post$beta1, 2, median), times=hmm.stan.dat$N) 
  # gamma : dim (6000, N years, 2=K)
  hmm_outputs[2,,] <- t( apply(post$gamma, c(2,3), median) )
  # alpha : dim (6000, N years, 2=K)
  hmm_outputs[3,,] <- t( apply(post$alpha, c(2,3), median) )
  # beta : dim (6000, N years, 2=K) # beta is ricker b 
  hmm_outputs[4,,] <- t( apply(post$beta, c(2,3), median) )
  # log_lik : dim(6000, N years) 
  hmm_outputs[5,,] <- rep( apply(post$log_lik, c(2), median), each=nlvl )
  
  
  hmm_out_ss3[[i]] <- hmm_outputs
  
  
}


## Visualizations

# Wrangle the list of arrays into a dataframe 
names(hmm_out_ss3) <- levels(sock$Stock)
hmm_out <- lapply(hmm_out_ss3, plyr::adply, .margins=c(2,3))
hmm_out <- dplyr::bind_rows(hmm_out, .id="stock")
names(hmm_out)[2:3] <- c("state_K", "BY")
hmm_out$BY <- as.numeric(as.character(hmm_out$BY)) # this is so dumb
hmm_tidy_ss3 <- hmm_out

## Plot gamma (state 1) over time
pdf("./figures/hmm-ss/ss3_gamma_timeseries.pdf")
for(i in 1:nlevels(sock$Stock)){
  g <- hmm_tidy_ss3 %>% filter(state_K == "State 1", stock==levels(sock$Stock)[i]) %>%
    ggplot() +
    geom_abline(intercept=0.5, slope=0, colour="grey70", lty="dashed") +
    geom_point(aes(x=BY, y=gamma)) +
    geom_line(aes(x=BY, y=gamma)) + 
    labs(title=levels(sock$Stock)[i]) + theme_minimal() 
  print(g)
}
dev.off()


## Diagnostics
elpd_ss1 <- list()
elpd_ss2 <- list()
elpd_ss3 <- list()
comp_ss1_ss2 <- list()
comp_ss1_ss3 <- list()

for(i in 1:nlevels(sock$Stock)){
  ## Load fits for version 1 
  load(paste0("./output/models/hmm-ss/hmm_ss1_", levels(sock$Stock)[i], ".Rdata"), verbose=T)  
  post_ss1 <- rstan::extract(hmm_ss1)
  # Save elpd
  elpd_ss1[[i]] <- loo::loo(posterior::as_draws_matrix(post_ss1$log_lik))

  ## Load fits for version 2
  load(paste0("./output/models/hmm-ss/hmm_ss2_", levels(sock$Stock)[i], ".Rdata"), verbose=T)
  post_ss2 <- rstan::extract(hmm_ss2)
  # Save elpd
  elpd_ss2[[i]] <- loo::loo(posterior::as_draws_matrix(post_ss2$log_lik))
  
  ## Load fits for version 2
  load(paste0("./output/models/hmm-ss/hmm_ss3_", levels(sock$Stock)[i], ".Rdata"), verbose=T)
  post_ss3 <- rstan::extract(hmm_ss3)
  # Save elpd
  elpd_ss3[[i]] <- loo::loo(posterior::as_draws_matrix(post_ss3$log_lik))
  
  ## -- Comparisons -- ##
  # 2-state versus 3-state model
  comp_ss1_ss2[[i]] <- loo::loo_compare(elpd_ss1[[i]], elpd_ss2[[i]])
  
  # 2,1 versus 4,1 transition matrix
  comp_ss1_ss3[[i]] <- loo::loo_compare(elpd_ss1[[i]], elpd_ss3[[i]])
}


# LOO example
# elpd1=loo::loo(posterior::as_draws_matrix(ll1_df)) # insert the log-likelihood
# loo::loo_compare(elpd1,elpd2)



## IOHMM ----- ##

hmm.dat <- dplyr::filter(sock, Stock == levels(sock$Stock)[1]) 

hmm.stan.dat <- list( 
  N = nrow(hmm.dat),
  R_S = hmm.dat$lnRS,
  S = hmm.dat$S,
  K = nlvl,
  X1 = hmm.dat$early_sst_stnd, # covariate = SST
  alpha_dirichlet = matrix(c(2,1,1,2), nrow=2) )


# Kept all the specs from hmm.fit script ??
iohmm <- rstan::stan(file = "./stan/iohmm_test.stan",
                       data = hmm.stan.dat,
                       warmup = 1000,
                       iter = 4000,
                       cores = 4,
                       chains = 4,
                       seed = 123,
                       control = list(adapt_delta = 0.90,
                                      max_treedepth = 10))
save(hmm_ss1, file = paste0("./output/models/hmm-ss/hmm_ss1_", levels(sock$Stock)[i], ".Rdata"))

