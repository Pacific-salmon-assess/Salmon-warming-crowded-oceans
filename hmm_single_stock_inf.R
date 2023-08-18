## Inference from Single-stock HMM models ##
## ----------------------------------------------------------------------------

if(!dir.exists("./figures/hmm/single-stock/")) dir.create("./figures/hmm/single-stock/")


# Set colours
col.stock  <- chroma::qpal(7, alpha = 0.4)[c(1, 3, 5, 7)]
col.region <- chroma::qpal(7, luminance = 40)[c(1, 3, 5, 7)]
col.lt <- chroma::qpal(7)[c(1, 3, 5, 7)]
col.dk <- chroma::qpal(7, luminance = 30)[c(1, 3, 5, 7)]


## HMM with Autocorrelation : SST
## ------------------------------------------------------------ ##
# Wrangle list into output dataframe
names(hmm_ac_out_sst) <- levels(sock$Stock)
post.df <- plyr::ldply(hmm_ac_out_sst, .id = "stock", .fun=function(x) { 
  data.frame(beta1_mu = x[grep('^beta1', rownames(x)), "mean"],
             beta1_2.5 = x[grep('^beta1', rownames(x)), "2.5%"],
             beta1_20 = x[grep('^beta1', rownames(x)), "20%"],
             beta1_80 = x[grep('^beta1', rownames(x)), "80%"],
             beta1_97.5 = x[grep('^beta1', rownames(x)), "97.5%"],
             gamma_mu = x[grep('^gamma', rownames(x)), "mean"],
             gamma_2.5 = x[grep('^gamma', rownames(x)), "2.5%"],
             gamma_97.5 = x[grep('^gamma', rownames(x)), "97.5%"],
             state = c(1,2),
             row.names=NULL) } )
post.df <- plyr::ddply(post.df, .variables="stock", 
                          .fun=function(x){
                            BY <- rep(sock$BY[sock$Stock == x$stock[1]], each=2)
                            region <- rep(sock$Ocean.Region2[sock$Stock == x$stock[1]], each=2)
                            cbind(x, BY, region)})


## (1) Plot gamma1 (state 1 probability) by stock
pdf("./figures/hmm/single-stock/hmm_ac_sst_gamma.pdf")
for(i in 1:length(hmm_ac_out_sst)){
  g <- post.df %>% filter(state == 1, stock==levels(sock$Stock)[i]) %>%
    ggplot() +
    geom_abline(intercept=0.5, slope=0, colour="grey70", lty="dashed") +
    geom_point(aes(x=BY, y=gamma_mu)) +
    geom_line(aes(x=BY, y=gamma_mu)) + 
    labs(title=levels(sock$Stock)[i], x="Brood Year", y="State 1 probability") +
    theme_minimal() 
  print(g)
}
dev.off()

## (2) Multi-panel grid of realized covariate "realb"
stk.summ.df <- post.df %>% mutate(prod_mu = beta1_mu*gamma_mu,
                              prod_2.5 = beta1_2.5*gamma_2.5,
                              prod_97.5 = beta1_97.5*gamma_97.5) %>% 
  dplyr::summarize(realb_mu = sum(prod_mu),
                   realb_2.5 = sum(prod_2.5),
                   realb_97.5 = sum(prod_97.5),
                   .by=c("stock", "BY", "region") )
g <- ggplot(stk.summ.df) + 
  geom_line(aes(x=BY, y=realb_mu, col=region)) +
  scale_colour_manual(values=col.region) +
  facet_wrap(vars(stock), as.table=F, nrow = 10, ncol=6) +
  theme_minimal() + theme(strip.text = element_text(size=7),
                          axis.text.x = element_text(angle=45)) +
  labs(x="Brood Year", y="SST coefficient", col=NULL) + ylim(c(-.5, .5))
pdf("figures/hmm/single-stock/hmm_ac_sst_grid.pdf", height=10, width=6)
print(g)
dev.off()

## (3) Grouped realized covariate "realb"
reg.summ.df <- stk.summ.df %>% dplyr::summarize(realb_reg_mu = mean(realb_mu, na.rm=T),
                                                realb_reg_2.5_avg = quantile(realb_mu, 0.2),
                                                realb_reg_97.5_avg = quantile(realb_mu, 0.8),
                                         .by=c("BY", "region"))
g <- ggplot(reg.summ.df) + 
  geom_line(data=stk.summ.df, aes(x=BY, y=realb_mu, group=stock, col=region), alpha=0.2) +
  geom_line(aes(x=BY, y=realb_reg_mu, col=region), linewidth=1) +
  geom_ribbon(aes(x=BY, ymin=realb_reg_2.5_avg, ymax=realb_reg_97.5_avg, fill=region), alpha=0.2) +
  facet_wrap(vars(region), nrow=2) +
  scale_colour_manual(values=col.region, aesthetics = c("colour", "fill")) +
  theme_minimal() + theme(legend.position="none") +
  ylim(c(-0.6, 0.6)) + labs(x="Brood Year", y="SST coefficient")
pdf("figures/hmm/single-stock/hmm_ac_sst_grouped.pdf")
print(g)
dev.off()

## (4) Stacked timeseries of realized covariate
g <- ggplot(stk.summ.df) + 
  geom_line(aes(x=BY, y=realb_mu, col=region)) + 
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
pdf("./figures/hmm/single-stock/hmm_ac_sst_stack.pdf", width=5, height=10)
print(g)
dev.off()

## (5) Posterior states -- density
# Posterior density data ## prefer not to do this loop here
for(i in 1:nlevels(sock$Stock)){
  load(paste0("./output/models/hmm-ss/hmm_ac_sst_", levels(sock$Stock)[i], ".Rdata"))
  post <- rstan::extract(hmm_ac)
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
g <- ggplot(beta_df_master) + 
  geom_path(aes(x=x, y=y, group=state, col=Ocean.Region2), linewidth=0.3) + 
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
  labs(y="", x="Covariate effect") + xlim(c(-1, 1))
pdf("./figures/hmm/single-stock/hmm_ac_sst_statedens.pdf", width=5, height=10)
print(g)
dev.off()

## (6) Posterior states -- dot and whiskers
# fitx alpha issue
g <- ggplot(post.df) + 
  geom_point(aes(x=beta1_mu, y=stock, col=region)) +
  geom_segment(aes(x=beta1_2.5, xend=beta1_97.5, y=stock, yend=stock, col=region)) +
  geom_segment(aes(x=beta1_20, xend=beta1_80, y=stock, yend=stock, col=region), linewidth=3.5, alpha=0.01) +
  scale_colour_manual(values=col.region) + theme_minimal() + labs()
pdf("./figures/hmm/single-stock/hmm_ac_sst_statedot.pdf", width=5, height=10) 
print(g)
dev.off()


## HMM with autocorrelation: Competitors --------------------
## ---------------------------------------------------------

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
    

## 2-covariate HMM
## ------------------------------------------------ ## 
names(hmm_out_2covar) <- levels(sock$Stock)
hmm_out <- lapply(hmm_out_2covar, plyr::adply, .margins=c(2,3))
hmm_out <- dplyr::bind_rows(hmm_out, .id="stock")
names(hmm_out)[2:3] <- c("state_K", "BY")
hmm_out$BY <- as.numeric(as.character(hmm_out$BY)) # this is so dumb
hmm_tidy_2covar <- hmm_out
hmm_tidy_2covar <- left_join(hmm_tidy_2covar, sock.info[,c("Stock", "Ocean.Region2")], by=c('stock'='Stock'))
hmm_tidy_2covar$stock <- factor(hmm_tidy_2covar$stock, levels=levels(sock$Stock))

## Plot gamma (state 1) over time
# (1) Page for each stock
pdf("./figures/hmm-ss/2covar_gamma_timeseries.pdf")
for(i in 1:nlevels(sock$Stock)){
  g <- hmm_tidy_2covar %>% filter(state_K == 1, stock==levels(sock$Stock)[i]) %>%
    ggplot() +
    geom_abline(intercept=0.5, slope=0, colour="grey70", lty="dashed") +
    geom_point(aes(x=BY, y=gamma)) +
    geom_line(aes(x=BY, y=gamma)) + 
    labs(title=levels(sock$Stock)[i]) + theme_minimal() 
  print(g)
}
dev.off()

## Plot beta1 and beta2 over time
g <- ggplot(hmm_tidy_2covar) +
  geom_point(aes(x=beta1, y=stock, col=Ocean.Region2, shape=Ocean.Region2)) + 
  geom_line(aes(x=beta1, y=stock, col=Ocean.Region2), alpha=0.5) +
  scale_colour_manual(values=col.region, name="Region") +
  scale_shape_manual(values=c(22, 21, 24, 23), name="Region") +
  theme_minimal() + theme(legend.position="none") +
  labs(x="State SST effects", y="")
g2 <- ggplot(hmm_tidy_2covar) +
  geom_point(aes(x=beta2, y=stock, col=Ocean.Region2, shape=Ocean.Region2)) + 
  geom_line(aes(x=beta2, y=stock, col=Ocean.Region2), alpha=0.5) +
  scale_colour_manual(values=col.region, name="Region") +
  scale_shape_manual(values=c(22, 21, 24, 23), name="Region") +
  theme_minimal() + 
  labs(x="State Comp effects", y="")

pdf("./figures/hmm-ss/b1_dot_2covar.pdf", width = 8, height=6)
cowplot::plot_grid(g, g2, ncol=2, rel_widths = c(1, 1.3))
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






