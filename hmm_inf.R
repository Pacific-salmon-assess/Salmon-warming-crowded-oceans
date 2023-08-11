## Inference from HMM models ##
## ----------------------------------------------------------------------------



# Load saved outputs (if necessary)


# Set colours
col.stock  <- chroma::qpal(7, alpha = 0.4)[c(1, 3, 5, 7)]
col.region <- chroma::qpal(7, luminance = 40)[c(1, 3, 5, 7)]
col.lt <- chroma::qpal(7)[c(1, 3, 5, 7)]
col.dk <- chroma::qpal(7, luminance = 30)[c(1, 3, 5, 7)]

## hmm_ss1 [Single Stock version 1] : 2 states, original transition matrix (2,1,1,2)
## --------------------------------------------------------------------------------

# Wrangle the list of arrays into a dataframe 
names(hmm_out_ss1) <- levels(sock$Stock)
hmm_out <- lapply(hmm_out_ss1, plyr::adply, .margins=c(2,3))
hmm_out <- dplyr::bind_rows(hmm_out, .id="stock")
names(hmm_out)[2:3] <- c("state_K", "BY")
hmm_out$BY <- as.numeric(as.character(hmm_out$BY)) # this is so dumb
hmm_tidy_ss1 <- hmm_out


## Plot gamma (state 1) over time
  # (1) Page for each stock
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

  # (1a) Posterior estimates of covariate states - beta
    # Load and format data
for(i in 1:nlevels(sock$Stock)){
  load(paste0("./output/models/hmm-ss/hmm_ss1_", levels(sock$Stock)[i], ".Rdata"))
  post <- rstan::extract(hmm_ss1)
  post_beta <- lapply(col_density(post$beta1, plot.it=F), plyr::adply, .margins=c(1,2))
  beta_df <- join(post_beta$x, post_beta$y, by=c("X1","X2"))
  beta_df$Stock = levels(sock$Stock)[i]
  names(beta_df) <- c("n", "state", "x", "y", "Stock")
  if(i==1) 
    beta_df_master <- beta_df else
    beta_df_master <- rbind(beta_df_master, beta_df)
}
beta_df_master <- left_join(beta_df_master, sock.info[,c("Stock", "Ocean.Region2")], by='Stock')
beta_df_master$Stock <- factor(beta_df_master$Stock, levels=levels(sock$Stock))

    # Plot  
g <- ggplot(beta_df_master) + geom_path(aes(x=x, y=y, col=Ocean.Region2)) + 
  facet_grid(rows=vars(Stock), as.table=F, switch="y") + 
  scale_colour_manual(values=col.region) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.y.left = element_text(angle=0),
        strip.background = element_rect(fill="transparent", colour="transparent"),
        strip.text = element_text(size=7),
        panel.spacing.y = unit(0, unit="cm"),
        panel.background = element_rect(fill="white"),
        legend.position = "none") +
  labs(y="", x="Covariate effect")
pdf("./figures/hmm-ss/beta1_post.pdf", width=5, height=10)
print(g)
dev.off()


  # (2) Multi-panel grid
hmm_tidy_ss1 <- left_join(hmm_tidy_ss1, sock.info[,c("Stock", "Ocean.Region2")], by=c('stock'='Stock'))
hmm_tidy_ss1$stock <- factor(hmm_tidy_ss1$stock, levels=levels(sock$Stock))
g <- hmm_tidy_ss1 %>% filter(state_K == "State 1") %>%
      ggplot() + geom_line(aes(x=BY, y=gamma, col=Ocean.Region2)) +
      geom_abline(intercept=0.5, slope=0, colour="grey70", lty="dashed") +
      scale_colour_manual(values=col.region) +
      facet_wrap(vars(stock), as.table=F) +
      theme_minimal() +
      labs(x="Brood Year", y="Gamma / Pr(state 1)", col=NULL)
pdf("figures/hmm-ss/ss1_gamma_grid.pdf")
print(g)
dev.off()


  # (3) Grouped by stock
ss1_avg <- hmm_tidy_ss1 %>% group_by(state_K, BY, Ocean.Region2) %>% dplyr::summarize(mean_gamma = mean(gamma, na.rm=TRUE), sd = sd(gamma, na.rm=T))
ss1_med <- hmm_tidy_ss1 %>% group_by(state_K, BY, Ocean.Region2) %>% dplyr::summarize(med_gamma = median(gamma), g25 = quantile(gamma, 0.25), g75 = quantile(gamma, 0.75))
g <- ggplot(hmm_tidy_ss1[hmm_tidy_ss1$state_K=="State 1",]) + 
  geom_line(aes(x=BY, y=gamma, group=stock, col=Ocean.Region2), alpha=0.2) +
  geom_line(data=ss1_avg[ss1_avg$state_K == "State 1",], aes(x=BY, y=mean_gamma, col=Ocean.Region2), linewidth=1) +
  geom_ribbon(data=ss1_avg[ss1_avg$state_K == "State 1",], aes(x=BY, y=mean_gamma, ymin=mean_gamma-sd, ymax=mean_gamma+sd, fill=Ocean.Region2), alpha=0.2) +
  facet_wrap(vars(Ocean.Region2), nrow=2) +
  scale_colour_manual(values=col.region, aesthetics = c("colour", "fill")) +
  theme_minimal()
pdf("figures/hmm-ss/ss1_gamma_grouped.pdf")
print(g)
dev.off()


## Plot ss1 covariates (beta1) 

  # (1) SST effects of each state by stock
g <- ggplot(hmm_tidy_ss1) +
      geom_point(aes(x=beta1, y=stock, col=Ocean.Region2, shape=Ocean.Region2)) + 
      geom_line(aes(x=beta1, y=stock, col=Ocean.Region2), alpha=0.5) +
      scale_colour_manual(values=col.region, name="Region") +
      scale_shape_manual(values=c(22, 21, 24, 23), name="Region") +
      theme_minimal() +
      labs(x="State SST effects", y="")
pdf("./figures/hmm-ss/b1_dot_ss1.pdf")
print(g)
dev.off()

  # (2) Dot plot of timeseries means of posterior mean covariate effect & 95% quant 

b1_ss1_q <- hmm_tidy_ss1 %>% 
  mutate(gamma_prod = beta1*gamma) %>% group_by(stock, BY) %>% 
  dplyr::summarize(covar = sum(gamma_prod)) %>% ungroup() %>%
  dplyr::summarize(mean_covar = mean(covar), lower=quantile(covar, 0.025), upper=quantile(covar, 0.975), .by=stock) %>%
  left_join(sock.info[,c("Stock", "Ocean.Region2")], by=c('stock'='Stock'))

# This one - but think about if we want some kind of regional grouping
g <- ggplot(b1_ss1_q) + geom_point(aes(x=mean_covar, y=stock, col=Ocean.Region2, shape=Ocean.Region2)) + 
      geom_segment(aes(y=stock, yend=stock, x=lower, xend=upper, col=Ocean.Region2)) +
      scale_colour_manual(values=col.region) + theme_minimal() +
      scale_shape_manual(values=c(15:18)) + 
      labs(y="", x="SST effect timeseries means", col="", shape="")
pdf("./figures/hmm-ss/ss1_dot_covar.pdf")
print(g)
dev.off()


  # (3) Stacked time series of posterior means
b1_ss1 <- hmm_tidy_ss1 %>% 
  mutate(gamma_prod = beta1*gamma) %>% group_by(stock, BY) %>% 
  dplyr::summarize(covar = sum(gamma_prod)) %>%
  left_join(sock.info[,c("Stock", "Ocean.Region2")], by=c('stock'='Stock'))


g <- ggplot(b1_ss1) + geom_line(aes(x=BY, y=covar, col=Ocean.Region2)) + 
  facet_grid(rows=vars(stock), switch ="y", scales="free_y", as.table=F) + 
  scale_colour_manual(values=col.region) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.y.left = element_text(angle=0),
        strip.background = element_rect(fill="transparent", colour="transparent"),
        strip.text = element_text(size=7, ),
        panel.spacing.y = unit(0, unit="cm"),
        panel.background = element_rect(fill="white"),
        legend.position = "none") +
  labs(y="", x="Brood Year")

pdf("./figures/hmm-ss/ss1_stacked_ts.pdf", width=5, height=10)
print(g)
dev.off()




## hmm_ss2 [Single Stock version 2] : 3 states, 2,1,1 transition matrix
## --------------------------------------------------------------------------------

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



## Input-output HMM (IOHMM) : varying transition matrix
## ------------------------------------------------------------ ##

# Wrangle the list of arrays into a dataframe 
names(iohmm_out_ss1) <- levels(sock$Stock)
hmm_out <- lapply(iohmm_out_ss1, plyr::adply, .margins=c(2,3))
hmm_out <- dplyr::bind_rows(hmm_out, .id="stock")
names(hmm_out)[2:3] <- c("state_K", "BY")
hmm_out$BY <- as.numeric(as.character(hmm_out$BY)) # this is so dumb
iohmm_tidy_ss1 <- hmm_out

## Plot gamma (state 1) over time
pdf("./figures/hmm-ss/iohmm_ss1_gamma_timeseries.pdf")
for(i in 1:nlevels(sock$Stock)){
  g <- iohmm_tidy_ss1 %>% filter(state_K == "State 1", stock==levels(sock$Stock)[i]) %>%
    ggplot() +
    geom_abline(intercept=0.5, slope=0, colour="grey70", lty="dashed") +
    geom_point(aes(x=BY, y=gamma)) +
    geom_line(aes(x=BY, y=gamma)) + 
    labs(title=levels(sock$Stock)[i]) + theme_minimal() 
  print(g)
}
dev.off()



## HMM with Autocorrelation test 
## ------------------------------------------------------------ ##
  #SST
names(hmm_ac_out_ss1) <- levels(sock$Stock)[1:length(hmm_ac_out_ss1)]
hmm_out <- lapply(hmm_ac_out_ss1, plyr::adply, .margins=c(2,3))
hmm_out <- dplyr::bind_rows(hmm_out, .id="stock")
names(hmm_out)[2:3] <- c("state_K", "BY")
hmm_out$BY <- as.numeric(as.character(hmm_out$BY)) # this is so dumb
hmm_ac_tidy_ss1 <- hmm_out

## Plot gamma (state 1) over time
pdf("./figures/hmm-ss/hmm_ac_ss1_gamma_timeseries.pdf")
for(i in 1:length(hmm_ac_out_ss1)){
  g <- hmm_ac_tidy_ss1 %>% filter(state_K == "State 1", stock==levels(sock$Stock)[i]) %>%
    ggplot() +
    geom_abline(intercept=0.5, slope=0, colour="grey70", lty="dashed") +
    geom_point(aes(x=BY, y=gamma)) +
    geom_line(aes(x=BY, y=gamma)) + 
    labs(title=levels(sock$Stock)[i]) + theme_minimal() 
  print(g)
}
dev.off()

  # Comp
names(hmm_ac_out_comp) <- levels(sock$Stock)[1:length(hmm_ac_out_comp)]
hmm_out <- lapply(hmm_ac_out_comp, plyr::adply, .margins=c(2,3))
hmm_out <- dplyr::bind_rows(hmm_out, .id="stock")
names(hmm_out)[2:3] <- c("state_K", "BY")
hmm_out$BY <- as.numeric(as.character(hmm_out$BY)) # this is so dumb
hmm_ac_tidy_comp <- hmm_out

## Plot gamma (state 1) over time
pdf("./figures/hmm-ss/hmm_ac_comp_gamma_timeseries.pdf")
for(i in 1:length(hmm_ac_out_comp)){
  g <- hmm_ac_tidy_comp %>% filter(state_K == "State 1", stock==levels(sock$Stock)[i]) %>%
    ggplot() +
    geom_abline(intercept=0.5, slope=0, colour="grey70", lty="dashed") +
    geom_point(aes(x=BY, y=gamma)) +
    geom_line(aes(x=BY, y=gamma)) + 
    labs(title=levels(sock$Stock)[i]) + theme_minimal() 
  print(g)
}
dev.off()
    


## Model comparison 
## --------------------------------------------------------------------------

# Realized covariate (beta1) over time - dynamic (RW) + 2x HMM comparison

# Load dynamic model outputs
load("./output/hbm_gamma_diff.RData", verbose=T)

pdf("./figures/hmm-ss/model-compare_timeseries.pdf")
for(i in 1:nlevels(sock$Stock)){
  rw_covar_ts <- filter(hbm.gamma.diff, Stock == levels(sock$Stock)[i])
  
  b1_ss1 <- hmm_tidy_ss1 %>% filter(stock==levels(sock$Stock)[i]) %>% 
    mutate(gamma_prod = beta1*gamma) %>% dplyr::summarize(covar = sum(gamma_prod), .by=BY)
  b1_ss2 <- hmm_tidy_ss2 %>% filter(stock==levels(sock$Stock)[i]) %>% 
    mutate(gamma_prod = beta1*gamma) %>% dplyr::summarize(covar = sum(gamma_prod), .by=BY)
  b1_iohmm <- iohmm_tidy_ss1 %>% filter(stock==levels(sock$Stock)[i]) %>% 
    mutate(gamma_prod = beta1*gamma) %>% dplyr::summarize(covar = sum(gamma_prod), .by=BY)
  
  g <- ggplot() +
    geom_line(data = b1_ss1, aes(x=BY, y=covar, col='Single-Stock HMM')) + # HMM ss1
    geom_line(data = rw_covar_ts, aes(x=BY, y=gamma, col = 'Hierarchical RW')) + # RW
    geom_line(data = b1_iohmm, aes(x=BY, y=covar, col='Single-Stock \n Input-Output HMM')) + # IOHMM
    #geom_line(data = b1_ss2, aes(x=BY, y=covar, col='Three-state HMM')) +
    labs(title=levels(sock$Stock)[i], y = "SST effect", x = "Year", col = NULL) + 
    scale_colour_manual(values = c('Single-Stock HMM'='mediumblue', 'Hierarchical RW'='maroon', 'Single-Stock \n Input-Output HMM'='orange', 'Three-state HMM'='darkgreen')) +
    theme_minimal()
  print(g)
}
dev.off()




# Realized covariate (beta1) over time - Compare HMM with Autocorrelation vs without

rho <- vector(length=length(hmm_ac_out_ss1))

pdf("./figures/hmm-ss/AC-model-compare_timeseries.pdf")
for(i in 1:length(hmm_ac_out_ss1)){
  rw_covar_ts <- filter(hbm.gamma.diff, Stock == levels(sock$Stock)[i])
  
  b1_ss1 <- hmm_tidy_ss1 %>% filter(stock==levels(sock$Stock)[i]) %>% 
    mutate(gamma_prod = beta1*gamma) %>% dplyr::summarize(covar = sum(gamma_prod), .by=BY)
  b1_ac <- hmm_ac_tidy_ss1 %>% filter(stock==levels(sock$Stock)[i]) %>% 
    mutate(gamma_prod = beta1*gamma) %>% dplyr::summarize(covar = sum(gamma_prod), .by=BY)
  rho[i] <- filter(hmm_ac_tidy_ss1, stock==levels(sock$Stock)[i])[1,"rho"]
  
  g <- ggplot() +
    geom_line(data = b1_ss1, aes(x=BY, y=covar, col='HMM (no AC)')) +
    geom_line(data = b1_ac, aes(x=BY, y=covar, col='HMM with Autocorrelation')) + 
    labs(title=levels(sock$Stock)[i], subtitle = paste("Rho =", round(rho[i], 2)), y = "SST effect", x = "Year", col = NULL) + 
    scale_colour_manual(values = c('HMM (no AC)'='mediumblue', 'HMM with Autocorrelation'='orange')) +
    theme_minimal()
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






