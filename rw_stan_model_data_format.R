#library(here);library(dplyr)

#funcs
make_design_matrix <- function(x,grp) {
  x2 <- matrix(nrow = length(x), ncol = length(unique(grp)))
  for(i in 1:length(unique(grp))){
    x2[,i] <- ifelse(grp == levels(factor(grp))[i], 1, 0)*x
  }
  return(x2)
}


cmdstanr::set_cmdstan_path(path='')
library(cmdstanr)
options(mc.cores=4)

# Data
sock = read.csv(here('data','sockeye', 'master_sockeye_brood_table_covar.csv'))
sock.info = read.csv(here('data','sockeye', 'master_sockeye_stock_info.csv'))


#file.tv=file.path(cmdstanr::cmdstan_path(),'sr models', "ind_tvalpha_ricker.stan")
file.tv=file.path(here('stan', 'ind_tvalpha_ricker.stan'))
m.rw=cmdstanr::cmdstan_model(file.tv) #COMPILE


#random walk data setup
scy.info=sock.info #subset(stk.info,species=='Sockeye')
scy.dat=sock #stk.dat[stk.dat$stock.id%in%scy.info$stock.id,]

smax.scy=scy.dat%>%group_by(Stock.ID)%>%summarize(max.S=max(S))

S.mat.scy=make_design_matrix(scy.dat$S,grp=scy.dat$Stock.ID)

stk.year=expand.grid(levels(factor(scy.dat$Stock.ID)),seq(min(scy.dat$BY),max(scy.dat$BY)))
stk.year[,3]=paste(stk.year[,1],stk.year[,2],sep='.')
stk.year=stk.year[order(stk.year[,1],stk.year[,2]),]
scy.dat$stk.year=match(paste(scy.dat$Stock.ID,scy.dat$BY,sep='.'),stk.year[,3])

dl.rw=list(N=nrow(scy.dat),
           J=length(unique(scy.dat$Stock.ID)),
           L=max(scy.dat$BY)-min(scy.dat$BY)+1,
           R=length(unique(scy.info$Ocean.Region2)),
           J_i=as.numeric(factor(scy.dat$Stock.ID)),
           J_ii=scy.dat$stk.year,
           J_or=as.numeric(factor(scy.info$Ocean.Region2)),
           R_S=scy.dat$lnRS,
           S=S.mat.scy,
           X1=scy.dat$early_sst_stnd,
           X2=scy.dat$np_pinks_sec_stnd,
           pSmax_mean=0.5*smax.scy$max.S,
           pSmax_sig=2*smax.scy$max.S)

init_fun <- function() list(
  sigma_t = rep(0.5, dl.rw$J),
  p_sig = replicate(dl.rw$J, c(0.4, 0.3, 0.3), simplify = FALSE)
)

tv.scy_rw=     m.rw$sample(data=dl.rw,
                            chains = 1,
                            iter_warmup = 500,
                            iter_sampling =1000,
                            #refresh = 250,
                            #adapt_delta = 0.999,
                            #max_treedepth = 20,
                           init = init_fun) # SAMPLE


gamma.t.scy=tv.scy_rw$draws(variable='mu_g_r',format='draws_matrix')
m.gamma.t.scy=apply(gamma.t.scy,2,median)
l90.gamma.t.scy=apply(gamma.t.scy,2,quantile,0.025)
u90.gamma.t.scy=apply(gamma.t.scy,2,quantile,0.975)
gamma.jt.scy=tv.scy_rw$draws(variable='g_t',format='draws_matrix')

plot(m.gamma.t.scy[199:(199+66-1)]~seq(min(scy.dat$BY),max(scy.dat$BY)),ylim=c(-1,1),lwd=3,col='darkred',ylab='gamma',xlab='Brood Cohort Year',type='n',bty='l',main='Sockeye')
for(j in 1:dl.rw$J){
  logaj=apply(gamma.jt.scy[,grepl(paste(',',j,']',sep=''),colnames(gamma.jt.scy))],2,median)
  lines(logaj~seq(min(scy.dat$BY),max(scy.dat$BY)),lty=3,col=adjustcolor('darkgray',alpha.f=0.4))
  st=scy.info$yr_start[j]-min(scy.dat$BY)+1
  end=scy.info$yr_end[j]-min(scy.dat$BY)+1
  lines(logaj[st:end]~seq(scy.info$yr_start[j],scy.info$yr_end[j]),col=adjustcolor('darkgray',alpha.f=0.8))
}
lines(m.gamma.t.scy[199:(199+66-1)]~seq(min(scy.dat$BY),max(scy.dat$BY)),lwd=3,col='darkred')
lines(l90.gamma.t.scy[199:(199+66-1)]~seq(min(scy.dat$BY),max(scy.dat$BY)),lty=5,col='darkred')
lines(u90.gamma.t.scy[199:(199+66-1)]~seq(min(scy.dat$BY),max(scy.dat$BY)),lty=5,col='darkred')


# Run and extract with rstan

#Run
smax <- sock %>% group_by(Stock.ID) %>% summarize(max.S = max(S))
S.mat <- make_design_matrix(sock$S, grp = sock$Stock.ID)

stk.year <- expand.grid(levels(factor(sock$Stock.ID)), seq(min(sock$BY), max(sock$BY)))
stk.year[,3] <- paste(stk.year[,1], stk.year[,2], sep='.')
stk.year <- stk.year[order(stk.year[,1], stk.year[,2]),]
sock$stk.year <- match(paste(sock$Stock.ID, sock$BY, sep='.'), stk.year[,3])

stan.dat.rw <- list(N = nrow(sock),
           J = length(unique(sock$Stock.ID)),
           L = max(sock$BY) - min(sock$BY) + 1,
           R = length(unique(sock.info$Ocean.Region2)),
           J_i = as.numeric(factor(sock$Stock.ID)),
           J_ii = sock$stk.year,
           J_or = as.numeric(factor(sock.info$Ocean.Region2,
                              levels=unique(sock.info$Ocean.Region2))),
           R_S = sock$lnRS,
           S = S.mat.scy,
           X1 = sock$early_sst_stnd,
           X2 = sock$np_pinks_sec_stnd,
           pSmax_mean = 0.5*smax.scy$max.S,
           pSmax_sig = 2*smax.scy$max.S)



rw.fit <- rstan::stan(file = "./stan/ind_tvalpha_ricker.stan",
                      data = stan.dat.rw,
                      warmup = 1000,
                      iter = 2000,
                      cores = 4,
                      chains = 4,
                      seed = 123,
                      control = list(adapt_delta = 0.99,
                                     max_treedepth = 20))
save(rw.fit, file = here("output/models/dyn", speciesFlag, "dyn_new_2025.RData"))

# Load fit
load(here("output/models/dyn", speciesFlag, "dyn_new_2025.RData"), verbose=T)

# Extract

# Stock-specific dataframe
probs <- c(0.025, 0.05, 0.10, 0.50, 0.90, 0.95, 0.975)
summ <- rstan::summary(rw.fit, pars = c("g_t", "k_t"), probs = probs)[[1]]

df.dyn.st.2c <- data.frame(Stock = rep(info_master$Stock, times=stan.dat.rw$L*2),
                           Ocean.Region2 = rep(info_master$Ocean.Region2, times=stan.dat.rw$L*2),
                           BY = rep(rep(min(sock.info$yr_start):(min(sock.info$yr_start) +
                                                                   stan.dat.rw$L -1),
                                        each=nrow(info_master)), times=2),
                           mu = summ[, "mean"],
                           se = summ[, "se_mean"],
                           lower_10 = summ[, "10%"],
                           upper_90 = summ[ , "90%"],
                           var = case_when(str_extract(rownames(summ), "[a-z]+") == "g" ~ "gamma",
                                        str_extract(rownames(summ), "[a-z]+") == "k" ~ "kappa"),
                           varnam = case_when(grepl("^g", rownames(summ)) ~ "SST",
                                              grepl("^k", rownames(summ)) ~ "Competitors")
)
df.dyn.st.2c <- ocean_region_lab(df.dyn.st.2c)


# Summarized dataframe (regional-level)
summ.r <- rstan::summary(rw.fit, pars = c("mu_g_r", "mu_k_r"), probs = probs)[[1]]
df.dyn.reg.2c <- data.frame(Ocean.Region2 = rep(rev(unique(info_master$Ocean.Region2)),
                                                  times=stan.dat.rw$L*2), # REMOVE REV if models rerun
                              BY = rep(rep(min(sock.info$yr_start):(min(sock.info$yr_start) +
                                                                      stan.dat.rw$L -1),
                                           each=4), times=2),
                              mu = summ.r[, "mean"],
                              se = summ.r[, "se_mean"],
                              lower_10 = summ.r[, "10%"],
                              upper_90 = summ.r[ , "90%"],
                              var = case_when(grepl("^mu_g", rownames(summ.r)) ~ "gamma",
                                              grepl("^mu_k", rownames(summ.r)) ~ "kappa"),
                              varnam = case_when(grepl("^mu_g", rownames(summ.r)) ~ "SST",
                                                 grepl("^mu_k", rownames(summ.r)) ~ "Competitors")
                              )
if(speciesFlag == "pink"){
  df.dyn.reg.2c$even_odd <- ifelse(gtools::odd(df.dyn.reg.2c$BY), "odd", "even")
}

# Trim early year estimates that are not based on data - by average start year of region (NOT an ideal way, more nuance needed...)
reg_start_yr <- info_master %>% group_by(Ocean.Region2) %>% summarize(avg_start=round(mean(yr_start), 0))
df.dyn.reg.2c <- df.dyn.reg.2c %>% left_join(reg_start_yr,  by="Ocean.Region2") %>% filter(BY >= avg_start)
df.dyn.reg.2c <- ocean_region_lab(df.dyn.reg.2c)

write.csv(df.dyn.reg.2c, here('output', paste0('rw_2025_model_forplot_', speciesFlag, '.csv')))

# Plot
ggplot(df.dyn.st.2c) +
  geom_line(aes(x=BY, y=mu, group=Stock), col="grey50") +
  geom_line(data=df.dyn.reg.2c, aes(x=BY, y=mu, col=Ocean.Region2), linewidth=1) +
  facet_grid(rows=vars(Ocean.Region2), cols=vars(varnam)) +
  ylim(-1,1)


# Compare to "old" RW fits

# load and process old fits
load(here("output/models/dyn", speciesFlag, "hbm_dyn_2c.RData"), verbose=T)

# Stock-specific dataframe
probs <- c(0.025, 0.05, 0.10, 0.50, 0.90, 0.95, 0.975)
summ <- rstan::summary(dyn.2c, pars = c("gamma", "kappa"), probs = probs)[[1]]
df.dyn.st.2c.old <- data.frame(Stock = data_master$Stock,
                           Ocean.Region2 = data_master$Ocean.Region2,
                           BY = data_master$BY,
                           mu = summ[, "mean"],
                           se = summ[, "se_mean"],
                           lower_10 = summ[, "10%"],
                           upper_90 = summ[ , "90%"],
                           var = str_extract(rownames(summ), "[a-z]+"),
                           varnam = case_when(grepl("^gamma", rownames(summ)) ~ "SST",
                                              grepl("^kappa", rownames(summ)) ~ "Competitors")
)
df.dyn.st.2c.old <- ocean_region_lab(df.dyn.st.2c.old)

# Summarized dataframe (regional-level)
# gamma/kappa are series-specific; no mu output. Summarize stocks instead
if(speciesFlag=="pink"){
  df.dyn.reg.2c.old <- df.dyn.st.2c.old %>%
    mutate(even_odd = ifelse(grepl("Even$", Stock), "Even", "Odd")) %>%
    dplyr::summarize(reg_mean=mean(mu, na.rm=T),
                     n_stk=n_distinct(Stock),
                     lower_10=quantile(mu, 0.1),
                     upper_90=quantile(mu, 0.9),
                     .by=c(Ocean.Region2, BY, varnam, even_odd))
} else {
  df.dyn.reg.2c.old <- dplyr::summarize(df.dyn.st.2c.old,
                                    reg_mean=mean(mu, na.rm=T),
                                    n_stk=n_distinct(Stock),
                                    lower_10=quantile(mu, 0.1),
                                    upper_90=quantile(mu, 0.9),
                                    .by=c(Ocean.Region2, BY, varnam))
}
df.dyn.reg.2c.old <- ddply(df.dyn.reg.2c.old, .(Ocean.Region2), dplyr::filter, n_stk >= max(n_stk)*0.1) # remove years with less than 10% of stocks observed
df.dyn.reg.2c.old <- ocean_region_lab(df.dyn.reg.2c.old)


#Compare new vs old fits on one plot

ggplot(df.dyn.st.2c) + # plot stock trajectories
  geom_line(aes(x=BY, y=mu, group=Stock), col="darkred") +
  geom_line(data=df.dyn.st.2c.old, aes(x=BY, y=mu, group=Stock), col="grey50") +
  facet_grid(rows=vars(Ocean.Region2), cols=vars(varnam)) +
  theme_sleek() +
  ylim(-1,1)


ggplot() + # plot region trajectories
  geom_line(data=df.dyn.reg.2c, aes(x=BY, y=mu, col=Ocean.Region2), linewidth=1, linetype="solid") +
  geom_ribbon(data=df.dyn.reg.2c, aes(x=BY, ymin=lower_10, ymax=upper_90, fill=Ocean.Region2), alpha=0.2) +
  geom_line(data=df.dyn.reg.2c.old, aes(x=BY, y=reg_mean), col="grey25", linewidth=1) +
  geom_ribbon(data=df.dyn.reg.2c.old, aes(x=BY, ymin=lower_10, ymax=upper_90), fill="grey25", alpha=0.2) +
  facet_grid(rows=vars(Ocean.Region2), cols=vars(varnam)) +
  theme_sleek() +
  ylim(-1,1)



#GAM setup

file.gam3=file.path(cmdstanr::cmdstan_path(),'sr models', "gam_tvalpha_hier.stan")
m.gam=cmdstanr::cmdstan_model(file.gam3)  # skip over GAM for now

scy.info=subset(stk.info,species=='Sockeye')
scy.dat=stk.dat[stk.dat$stock.id%in%scy.info$stock.id,]

smax.scy=scy.dat%>%group_by(stock.id)%>%summarize(max.S=max(spawners))

#create knots at every 2 years - cubic (df=3) piecewise polynomials
bspline_year=splines::bs(scy.dat$broodyear,knots=seq(min(scy.dat$broodyear),max(scy.dat$broodyear),by=2),degree=3,intercept=F)
predspline_year=splines::bs(seq(min(scy.dat$broodyear),max(scy.dat$broodyear)),knots=seq(min(scy.dat$broodyear),max(scy.dat$broodyear),by=2),degree=3,intercept=F)

sr.scy=dplyr::distinct(scy.info,region,.keep_all=T)
sr.scy=sr.scy[order(sr.scy$region),] #dataframe of ordered subregions

dl2=list(N=nrow(scy.dat),
         J=length(unique(scy.dat$stock.id)),
         L=max(scy.dat$broodyear)-min(scy.dat$broodyear)+1,
         J_i=as.numeric(factor(scy.dat$stock.id)),
         sr_or=as.numeric(factor(sr.scy$ocean.basin)),
         j_or=as.numeric(factor(scy.info$ocean.basin)),
         j_sr=as.numeric(factor(scy.info$region)),
         R=length(unique(scy.info$ocean.basin)),
         SR=length(unique(scy.info$region)),
         R_S=scy.dat$logRS,
         S=scy.dat$spawners,
         Nb=ncol(bspline_year),
         B_year=bspline_year,
         B_pred_year=predspline_year,
         pSmax_mean=0.5*smax.scy$max.S,
         pSmax_sig=2*smax.scy$max.S)

tv.scy_k2=     m.gam$sample(data=dl2,
                             chains = 5,
                             iter_warmup = 500,
                             iter_sampling =500,
                             refresh = 250,
                             adapt_delta = 0.995,
                             max_treedepth = 20)

log.a.t.scy=tv.scy_k2$draws(variable='log_a_t',format='draws_matrix')
m.log.a.t.scy=apply(log.a.t.scy,2,median)
l90.log.a.t.scy=apply(log.a.t.scy,2,quantile,0.025)
u90.log.a.t.scy=apply(log.a.t.scy,2,quantile,0.975)
log.a.jt.scy=tv.scy_k2$draws(variable='log_a_jt',format='draws_matrix')

plot(m.log.a.t.scy~seq(min(scy.dat$broodyear),max(scy.dat$broodyear)),ylim=c(-1,4),lwd=3,col='darkred',ylab='Productivity (alpha) - Recruits/Spawner',xlab='Brood Cohort Year',type='n',bty='l',main='Sockeye')
for(j in 1:dl2$J){
  logaj=apply(log.a.jt.scy[,grepl(paste(',',j,']',sep=''),colnames(log.a.jt.scy))],2,median)
  lines(logaj~seq(min(scy.dat$broodyear),max(scy.dat$broodyear)),lty=3,col=adjustcolor('darkgray',alpha.f=0.4))
  st=scy.info$begin[j]-min(scy.dat$broodyear)+1
  end=scy.info$end[j]-min(scy.dat$broodyear)+1
  lines(logaj[st:end]~seq(scy.info$begin[j],scy.info$end[j]),col=adjustcolor('darkgray',alpha.f=0.8))
}
lines(m.log.a.t.scy~seq(min(scy.dat$broodyear),max(scy.dat$broodyear)),lwd=3,col='darkred')
lines(l90.log.a.t.scy~seq(min(scy.dat$broodyear),max(scy.dat$broodyear)),lty=5,col='darkred')
lines(u90.log.a.t.scy~seq(min(scy.dat$broodyear),max(scy.dat$broodyear)),lty=5,col='darkred')


log.a.r.scy=tv.scy_k2$draws(variable='log_a_rt',format='draws_matrix')
log.a.r.bs=log.a.r.scy[,grepl(',1]',colnames(log.a.r.scy))]
log.a.r.goa=log.a.r.scy[,grepl(',2]',colnames(log.a.r.scy))]
log.a.r.nc=log.a.r.scy[,grepl(',3]',colnames(log.a.r.scy))]
log.a.r.sc=log.a.r.scy[,grepl(',4]',colnames(log.a.r.scy))]

plot(m.log.a.t.scy~seq(min(scy.dat$broodyear),max(scy.dat$broodyear)),ylim=c(-1,3),lwd=3,col='darkred',ylab='Productivity (alpha) - Recruits/Spawner',xlab='Brood Cohort Year',type='n',bty='l',main='Sockeye')
lines(apply(log.a.r.bs,2,median)~seq(min(scy.dat$broodyear),max(scy.dat$broodyear)),col='goldenrod',lwd=2)
lines(apply(log.a.r.goa,2,median)~seq(min(scy.dat$broodyear),max(scy.dat$broodyear)),col='darkgreen',lwd=2)
lines(apply(log.a.r.nc,2,median)~seq(min(scy.dat$broodyear),max(scy.dat$broodyear)),col='navy',lwd=2)
lines(apply(log.a.r.sc,2,median)~seq(min(scy.dat$broodyear),max(scy.dat$broodyear)),col='orange',lwd=2)
