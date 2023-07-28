## Exploratory things - HIDDEN ##

## -- Quick look at filtering/infilling in skeena dataset --## 
library(dplyr)

skeena <- read.csv("C:/Users/hunterh/Documents/sockeye-climate-competition/nonstationary_dynamics/data/raw data/sockeye/Skeena_Nass_sockeye.csv")

main <- filter(skeena, Label=="Main")
infill <- filter(skeena, Label=="Main_Infill")
filter45 <- filter(skeena, Label=="Filter45")

# See what was infilled
summary(filter(infill, Infill==TRUE)) # 36 total infills
unique(infill$Stock[which(infill$Infill==TRUE)]) # across 13 stocks #all but four are stocks we removed from analysis due to ts length


# See what was filtered
summary(filter(filter45, Filter==TRUE)) # 7 obs filtered out
unique(filter45$Stock[which(filter45$Filter==TRUE)]) # 5 stocks #all removed from analysis due to ts length


## -- Some unused plots from early HMM exploration (single-stock) --##

# Covariate effects for each state by stock
## Density plots of beta1 estimates # not that useful
hmm_tidy %>% ggplot(aes(x=beta1, fill=state_K)) + 
  geom_histogram() + theme_minimal()
# Dot plot of the covar effects by stock
asc.beta <- hmm_tidy$stock[match(sort(unique(hmm_tidy$beta1[which(hmm_tidy$state_K == "State 1")])), hmm_tidy$beta1)]
hmm_tidy %>% ggplot(aes(x=beta1, y=factor(stock, levels=asc.beta), colour=state_K)) +
  geom_point() # this actually hurts my eyes

## Many iterations of trying to plot timeseries means of posterior means 

  # Boxplots
b1_ss1_q <- hmm_tidy_ss1 %>% 
  mutate(gamma_prod = beta1*gamma) %>% group_by(stock, BY) %>% 
  dplyr::summarize(covar = sum(gamma_prod)) %>% ungroup() %>%
  dplyr::summarize(mean_covar = mean(covar), lower=quantile(covar, 0.05), upper=quantile(covar, 0.95), .by=stock) %>%
  left_join(sock.info[,c("Stock", "Ocean.Region2")], by=c('stock'='Stock'))

ggplot(b1_ss1_q) + geom_boxplot(aes(y=stock, xmin=lower, 
                                    xlower=lower, xmiddle=mean_covar, 
                                    xmax=upper, xupper=upper, 
                                    fill=Ocean.Region2, col=Ocean.Region2), alpha=0.5, stat="identity") +
  geom_point(data=hmm_tidy_ss1, aes(x=beta1, y=stock, col=Ocean.Region2), alpha=0.5, size=0.75) +
  scale_colour_manual(values=col.region, aesthetics=c("colour", "fill")) + theme_minimal()

  # Violin plots
  # Nice but not really enough space to show all the info
b1_ss1 <- hmm_tidy_ss1 %>% 
  mutate(gamma_prod = beta1*gamma) %>% group_by(stock, BY) %>% 
  dplyr::summarize(covar = sum(gamma_prod)) %>%
  left_join(sock.info[,c("Stock", "Ocean.Region2")], by=c('stock'='Stock'))

g <- ggplot(b1_ss1) + geom_violin(aes(x=covar, y=stock, col=Ocean.Region2, fill=Ocean.Region2)) + 
  geom_point(data=hmm_tidy_ss1, aes(x=beta1, y=stock, col=Ocean.Region2), alpha=0.5, size=0.75) +
  scale_colour_manual(values=col.region, aesthetics=c("colour", "fill")) +
  facet_wrap(vars(Ocean.Region2), scales="free_y") +
  theme_minimal() + theme(legend.position="none")

pdf("./figures/hmm-ss/ss1_coef_violin.pdf")
print(g)
dev.off()


  # Dot plots with both states and mean realized covar
b1_ss1 <- hmm_tidy_ss1 %>% 
  mutate(gamma_prod = beta1*gamma) %>% group_by(stock, BY) %>% 
  dplyr::summarize(covar = sum(gamma_prod)) %>%
  ungroup() %>% dplyr::summarize(mean_covar = mean(covar), .by=stock) %>% 
  left_join(sock.info[,c("Stock", "Ocean.Region2")], by=c('stock'='Stock'))

g <- ggplot(hmm_tidy_ss1) +
  geom_point(aes(x=beta1, y=stock, col=Ocean.Region2, shape=Ocean.Region2)) + 
  geom_line(aes(x=beta1, y=stock, col=Ocean.Region2)) +
  scale_colour_manual(values=col.region, name="Region") +
  geom_point(data=b1_ss1, aes(x=mean_covar, y=stock, col=Ocean.Region2, shape=Ocean.Region2), size=2) +
  scale_shape_manual(values=c(15,16,17,18), name="Region") +
  theme_minimal() +
  labs(x="Covariate effects", y="")



## At the end of IOHMM fitting - doing some exploration of the parameter outputs 
  
  ## TEST STOCKS ##
  test_nam <- c("Meziadin", "Bowron", "Egegik")
  test_ind <- which(levels(sock$Stock) %in% test_nam)
  
  ## compare w1 with gamma for iohmm
  load(paste0("./output/models/hmm-ss/iohmm_ss1_", test_nam[2], ".Rdata"), verbose=T)
  fit_summ <- rstan::summary(iohmm)$summary
  sst_ts <- sock$early_sst_stnd[sock$Stock == test_nam[2]]
  unA1 <- fit_summ["w[1]","mean"] * sst_ts 
  unA2 <- fit_summ["w[2]","mean"] * sst_ts
  A1 <- fit_summ[paste0("A[", 1:66, ",1]"),"mean"]
  softmax1 <- exp(unA1) / (exp(unA1) + exp(unA2))
  softmax2 <- exp(unA2) / (exp(unA1) + exp(unA2))
  A1 == softmax
  ggplot(aes(x=)) + geom_line()
  
  temp_iohmm <- iohmm_tidy_ss1 %>% filter(state_K == "State 1", stock==test_nam[2])
  temp_hmm <- hmm_tidy_ss1 %>% filter(state_K == "State 1", stock==test_nam[2])
  g <- ggplot() +
    geom_abline(intercept=0.5, slope=0, colour="grey70", lty="dashed") +
    geom_line(data=temp_iohmm, aes(x=BY, y=gamma, col='smoothed gamma[1]')) + 
    geom_line(data=NULL, aes(x=temp_iohmm$BY, y=softmax, col='softmax w[1]')) + 
    scale_colour_manual(values=c('smoothed gamma[1]'='orange', 'softmax w[1]'='purple')) +
    labs(title=test_nam[2], y=NULL, color=NULL) + theme_minimal() 
  png("./figures/hmm-ss/w1_gamma1_Meziadin.png")
  print(g)
  dev.off()
  
