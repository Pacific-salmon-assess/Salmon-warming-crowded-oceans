## Inference from Single-stock HMM models ##
## ----------------------------------------------------------------------------

if(!dir.exists("./figures/hmm/")) dir.create("./figures/hmm/")
if(!dir.exists("./figures/hmm/single-stock/")) dir.create("./figures/hmm/single-stock/")


# Set colours
col.stock  <- chroma::qpal(7, alpha = 0.4)[c(1, 3, 5, 7)]
col.region <- chroma::qpal(7, luminance = 40)[c(1, 3, 5, 7)]
col.lt <- chroma::qpal(7)[c(1, 3, 5, 7)]
col.dk <- chroma::qpal(7, luminance = 30)[c(1, 3, 5, 7)]


## HMM with Autocorrelation : SST
## ------------------------------------------------------------ ##

# Load fits 
if(!exists("hmm_ac_out_sst")){
  hmm_ac_out_sst <- list() 
  params_out <- c("beta1", "gamma", "alpha", "beta", "log_lik", "rho") # Parameters of interest
    for(i in 1:nlevels(sock$Stock)){
      load(file=paste0("./output/models/hmm-ss/hmm_ac_sst_", levels(sock$Stock)[i], ".Rdata"))
      hmm_ac_out_sst[[i]] <- rstan::summary(hmm_ac, pars = params_out, probs=c(0.025, 0.2, 0.5, 0.8, 0.975))$summary
      print(i)
    }
  save(hmm_ac_out_sst, file="./output/hmm_ac_out_sst.Rdata")
}

## Data ------------------------------
# Wrangle output list into a dataframe
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
                            BY <- rep(sock$BY[factor(sock$Stock, ordered=F) == x$stock[1]], each=2)
                            region <- rep(sock$Ocean.Region2[factor(sock$Stock, ordered=F) == x$stock[1]], each=2)
                            cbind(x, BY, region)})

# Make stock lvl summary dataframe
stk.summ.df <- post.df %>% mutate(prod_mu = beta1_mu*gamma_mu,
                                  prod_2.5 = beta1_2.5*gamma_2.5,
                                  prod_97.5 = beta1_97.5*gamma_97.5) %>% 
  dplyr::summarize(realb_mu = sum(prod_mu),
                   realb_2.5 = sum(prod_2.5),
                   realb_97.5 = sum(prod_97.5),
                   .by=c("stock", "BY", "region") )

# Make regional lvl summary dataframe
reg.summ.df <- stk.summ.df %>% dplyr::summarize(realb_reg_mu = mean(realb_mu, na.rm=T),
                                                realb_reg_2.5_avg = quantile(realb_mu, 0.2),
                                                realb_reg_97.5_avg = quantile(realb_mu, 0.8),
                                                .by=c("BY", "region"))

# Get posterior density data 
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


## Figures ----------------------------------
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


## (5) Posterior states -- dot and whiskers
g <- ggplot(unique(select(.data=post.df, stock, region, contains("beta1")))) + 
  geom_point(aes(x=beta1_mu, y=stock, col=region)) +
  geom_segment(aes(x=beta1_2.5, xend=beta1_97.5, y=stock, yend=stock, col=region)) +
  geom_segment(aes(x=beta1_20, xend=beta1_80, y=stock, yend=stock, col=region), linewidth=3.5, alpha=0.2) +
  scale_colour_manual(values=col.region) + theme_minimal() + 
  labs(col="Ocean Region", y="") # + coord_cartesian(xlim=c(-1.5, 1.5)) # option to zoom in
pdf("./figures/hmm/single-stock/hmm_ac_sst_statedot.pdf", width=5, height=10) 
print(g)
dev.off()


## (6) Posterior states -- density
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


## HMM with autocorrelation: Competitors --------------------
## ---------------------------------------------------------

# Load fits 
if(!exists("hmm_ac_out_comp")){
  hmm_ac_out_comp <- list() 
  params_out <- c("beta1", "gamma", "alpha", "beta", "log_lik", "rho") # Parameters of interest
  for(i in 1:nlevels(sock$Stock)){
    load(file=paste0("./output/models/hmm-ss/hmm_ac_comp_", levels(sock$Stock)[i], ".Rdata"))
    hmm_ac_out_comp[[i]] <- rstan::summary(hmm_ac, pars = params_out, probs=c(0.025, 0.2, 0.5, 0.8, 0.975))$summary
    print(i)
  }
  save(hmm_ac_out_comp, file="./output/hmm_ac_out_comp.Rdata")
}

## Data ------------------------------
# Wrangle output list into a dataframe
names(hmm_ac_out_comp) <- levels(sock$Stock)
post.df <- plyr::ldply(hmm_ac_out_comp, .id = "stock", .fun=function(x) { 
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
                         BY <- rep(sock$BY[factor(sock$Stock, ordered=F) == x$stock[1]], each=2)
                         region <- rep(sock$Ocean.Region2[factor(sock$Stock, ordered=F) == x$stock[1]], each=2)
                         cbind(x, BY, region)})

# Make stock lvl summary dataframe
stk.summ.df <- post.df %>% mutate(prod_mu = beta1_mu*gamma_mu,
                                  prod_2.5 = beta1_2.5*gamma_2.5,
                                  prod_97.5 = beta1_97.5*gamma_97.5) %>% 
  dplyr::summarize(realb_mu = sum(prod_mu),
                   realb_2.5 = sum(prod_2.5),
                   realb_97.5 = sum(prod_97.5),
                   .by=c("stock", "BY", "region") )

# Make regional lvl summary dataframe
reg.summ.df <- stk.summ.df %>% dplyr::summarize(realb_reg_mu = mean(realb_mu, na.rm=T),
                                                realb_reg_2.5_avg = quantile(realb_mu, 0.2),
                                                realb_reg_97.5_avg = quantile(realb_mu, 0.8),
                                                .by=c("BY", "region"))

# Get posterior density data 
for(i in 1:nlevels(sock$Stock)){
  load(paste0("./output/models/hmm-ss/hmm_ac_comp_", levels(sock$Stock)[i], ".Rdata"))
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


## Figures ----------------------------------
## (1) Plot gamma1 (state 1 probability) by stock
pdf("./figures/hmm/single-stock/hmm_ac_comp_gamma.pdf")
for(i in 1:length(hmm_ac_out_comp)){
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
g <- ggplot(stk.summ.df) + 
  geom_line(aes(x=BY, y=realb_mu, col=region)) +
  scale_colour_manual(values=col.region) +
  facet_wrap(vars(stock), as.table=F, nrow = 10, ncol=6) +
  theme_minimal() + theme(strip.text = element_text(size=7),
                          axis.text.x = element_text(angle=45)) +
  labs(x="Brood Year", y="Comp coefficient", col=NULL) + ylim(c(-.5, .5))
pdf("figures/hmm/single-stock/hmm_ac_comp_grid.pdf", height=10, width=6)
print(g)
dev.off()

## (3) Grouped realized covariate "realb"
g <- ggplot(reg.summ.df) + 
  geom_line(data=stk.summ.df, aes(x=BY, y=realb_mu, group=stock, col=region), alpha=0.2) +
  geom_line(aes(x=BY, y=realb_reg_mu, col=region), linewidth=1) +
  geom_ribbon(aes(x=BY, ymin=realb_reg_2.5_avg, ymax=realb_reg_97.5_avg, fill=region), alpha=0.2) +
  facet_wrap(vars(region), nrow=2) +
  scale_colour_manual(values=col.region, aesthetics = c("colour", "fill")) +
  theme_minimal() + theme(legend.position="none") +
  ylim(c(-1, 0.6)) + labs(x="Brood Year", y="Comp coefficient")
pdf("figures/hmm/single-stock/hmm_ac_comp_grouped.pdf")
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
pdf("./figures/hmm/single-stock/hmm_ac_comp_stack.pdf", width=5, height=10)
print(g)
dev.off()


## (5) Posterior states -- dot and whiskers
g <- ggplot(unique(select(.data=post.df, stock, region, contains("beta1")))) + 
  geom_point(aes(x=beta1_mu, y=stock, col=region)) +
  geom_segment(aes(x=beta1_2.5, xend=beta1_97.5, y=stock, yend=stock, col=region)) +
  geom_segment(aes(x=beta1_20, xend=beta1_80, y=stock, yend=stock, col=region), linewidth=3.5, alpha=0.2) +
  scale_colour_manual(values=col.region) + theme_minimal() + 
  labs(col="Ocean Region", y="") # + coord_cartesian(xlim=c(-1.5, 1.5)) # option to zoom in
pdf("./figures/hmm/single-stock/hmm_ac_comp_statedot.pdf", width=5, height=10) 
print(g)
dev.off()


## (6) Posterior states -- density
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
pdf("./figures/hmm/single-stock/hmm_ac_comp_statedens.pdf", width=5, height=10)
print(g)
dev.off()

    

## 2-covariate HMM --------------------------------
## ------------------------------------------------  

# Load fits 
if(!exists("hmm_ac_out_2c")){
  hmm_ac_out_2c <- list() 
  params_out <- c("beta1", "beta2", "gamma", "alpha", "beta", "log_lik")
  for(i in 1:nlevels(sock$Stock)){
    load(file=paste0("./output/models/hmm-ss/hmm_ac_2c_", levels(sock$Stock)[i], ".Rdata"))
    hmm_ac_out_2c[[i]] <- rstan::summary(hmm_ac_2c, pars = params_out, probs=c(0.025, 0.1, 0.5, 0.9, 0.975))$summary
  }
  save(hmm_ac_out_2c, file="./output/hmm_ac_out_2c.Rdata")
}

## Data
post.df <- hmm_param_df(hmm_ac_out_2c, summary="none")

stk.summ.df <- hmm_param_df(hmm_ac_out_2c, summary="stock")

reg.summ.df <- hmm_param_df(hmm_ac_out_2c, summary="reg")

## Figures --------------------------------------

## (1) Plot gamma1 (state 1 probability) by stock
pdf("./figures/hmm/single-stock/hmm_ac_2c_gamma.pdf")
for(i in 1:nlevels(sock$Stock)){
  g <- post.df %>% filter(state == 1, stock==levels(sock$Stock)[i]) %>%
    ggplot() +
    geom_abline(intercept=0.5, slope=0, colour="grey70", lty="dashed") +
    geom_point(aes(x=BY, y=gamma_mu) ) +
    geom_line(aes(x=BY, y=gamma_mu)) +
    labs(title=levels(sock$Stock)[i], x="Brood Year", y="State 1 probability") + theme_minimal() 
  print(g)
}
dev.off()

## (2) Multi-panel grid of realized covariate "realb"
# skip

## (3) Grouped realized covariate "realb" 
g_b1 <- ggplot(reg.summ.df) + 
  geom_line(data=stk.summ.df, aes(x=BY, y=realb1_mu, group=stock, col=region), alpha=0.2) +
  geom_line(aes(x=BY, y=realb1_reg_mu, col=region), linewidth=1) +
  geom_ribbon(aes(x=BY, ymin=realb1_reg_2.5_avg, ymax=realb1_reg_97.5_avg, fill=region), alpha=0.2) +
  facet_wrap(vars(region), nrow=1) +
  scale_colour_manual(values=col.region, aesthetics = c("colour", "fill")) +
  theme_minimal() + theme(legend.position="none", axis.text.x = element_blank()) +
  ylim(c(min(stk.summ.df$realb1_mu)-0.025, max(stk.summ.df$realb1_mu)+0.025)) + 
  labs(x="", y="SST coefficient")
g_b2 <- ggplot(reg.summ.df) + 
  geom_line(data=stk.summ.df, aes(x=BY, y=realb2_mu, group=stock, col=region), alpha=0.2) +
  geom_line(aes(x=BY, y=realb2_reg_mu, col=region), linewidth=1) +
  geom_ribbon(aes(x=BY, ymin=realb2_reg_2.5_avg, ymax=realb2_reg_97.5_avg, fill=region), alpha=0.2) +
  facet_wrap(vars(region), nrow=1) +
  scale_colour_manual(values=col.region, aesthetics = c("colour", "fill")) +
  theme_minimal() + theme(legend.position="none", strip.text = element_blank()) +
  ylim(c(min(stk.summ.df$realb2_mu)-0.025, max(stk.summ.df$realb2_mu)+0.025)) + 
  labs(x="Brood Year", y="Comp coefficient")
g <- cowplot::plot_grid(g_b1, g_b2, nrow=2)
pdf("figures/hmm/single-stock/hmm_ac_2c_grouped.pdf")
print(g)
dev.off()


## (4) Stacked timeseries of realized covariate  
g_b1 <- ggplot(stk.summ.df) + 
  geom_line(aes(x=BY, y=realb1_mu, col=region)) + 
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
g_b2 <- ggplot(stk.summ.df) + 
  geom_line(aes(x=BY, y=realb2_mu, col=region)) + 
  facet_grid(rows=vars(stock), switch ="y", scales="free_y", as.table=F) + 
  scale_colour_manual(values=col.region) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.spacing.y = unit(0, unit="cm"),
        panel.background = element_rect(fill="white"),
        legend.position = "none") +
  labs(y="", x="Brood Year")
g <- cowplot::plot_grid(g_b1, g_b2, nrow=1, rel_widths = c(1.25, 1))
pdf("./figures/hmm/single-stock/hmm_ac_2c_stack.pdf", width=8, height=10)
print(g)
dev.off()


## (5) Posterior states -- dot and whiskers 
g_b1 <- ggplot(unique(select(.data=post.df, stock, region, contains("beta")))) + 
  geom_point(aes(x=beta1_mu, y=stock, col=region)) +
  geom_segment(aes(x=beta1_2.5, xend=beta1_97.5, y=stock, yend=stock, col=region)) +
  geom_segment(aes(x=beta1_20, xend=beta1_80, y=stock, yend=stock, col=region), linewidth=3.5, alpha=0.2) +
  scale_colour_manual(values=col.region) + theme_minimal() + theme(legend.position="none") +
  labs(x="SST Coefficient", col="Ocean Region", y="") # + coord_cartesian(xlim=c(-1.5, 1.5)) # option to zoom in

g_b2 <- ggplot(unique(select(.data=post.df, stock, region, contains("beta")))) + 
  geom_point(aes(x=beta2_mu, y=stock, col=region)) +
  geom_segment(aes(x=beta2_2.5, xend=beta2_97.5, y=stock, yend=stock, col=region)) +
  geom_segment(aes(x=beta2_20, xend=beta2_80, y=stock, yend=stock, col=region), linewidth=3.5, alpha=0.2) +
  scale_colour_manual(values=col.region) + theme_minimal() + 
  theme(axis.text.y = element_blank(),
        legend.position="none") +
  labs(x="Comp Coefficient", col="Ocean Region", y="") # + coord_cartesian(xlim=c(-1.5, 1.5)) # option to zoom in

g <- cowplot::plot_grid(g_b1, g_b2, nrow=1, rel_widths = c(1.25, 1))
pdf("./figures/hmm/single-stock/hmm_ac_2c_statedot.pdf", width=8, height=10) 
print(g)
dev.off()



