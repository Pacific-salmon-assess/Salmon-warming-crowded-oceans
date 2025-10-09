## Test models with an interaction term ##
# Stationary model only to start #

# Species
if(speciesFlag=="pink") {
  data_master <- pink
  info_master <- pink.info} else if (speciesFlag=="pinkeven"){
    data_master <- pinkeven
    info_master <- pinkeven.info } else if (speciesFlag=="pinkodd"){
      data_master <- pinkodd
      info_master <- pinkodd.info } else if (speciesFlag=="chum") {
        data_master <- chum
        info_master <- chum.info } else if(speciesFlag=="sockeye"){
          data_master <- sock
          info_master <- sock.info }


# Set paths to output locations - dependent on species
fig.dir <- here("sensitivity-analyses", "interaction", speciesFlag) # place to store figures generated in this script
fit.dir <- here("sensitivity-analyses", "fits", speciesFlag) # place to store model fits
diag.dir <- here("sensitivity-analyses", "interaction", speciesFlag) # place to store diagnostics

# Make them if they don't exist
if(!dir.exists(fig.dir))
  dir.create(fig.dir, recursive = TRUE)
if(!dir.exists(fit.dir))
  dir.create(fit.dir, recursive = TRUE)
if(!dir.exists(diag.dir))
  dir.create(diag.dir, recursive = TRUE)

## Monitor params
pars.stat <- c("alpha", "beta", "sigma", "phi", "mu_alpha", "sigma_alpha",
               "gamma", "mu_gamma", "sigma_gamma",
               "kappa", "mu_kappa", "sigma_kappa",
               "chi", "mu_chi", "sigma_chi")
pars.gen.quant <- c("log_lik", "yhat", "yrep", "yresid") ## Generated quantities to monitor


## all yrs of data

## Run MCMC
stan.dat.all <- stan_data_stat(data_master,
                               scale.x1 = TRUE,
                               var.x2 = "early_sst_stnd",
                               var.x3 = "np_pinks_sec_stnd", # comp = pink abundance
                               var.region = "Ocean.Region2",
                               alpha.group = ifelse(speciesFlag=="sockeye", TRUE, FALSE)) # set to TRUE for sockeye
stat_inter <- rstan::stan(file = "./stan/hbm_stat_inter.stan",
                      data = stan.dat.all,
                      pars = c(pars.stat, pars.gen.quant),
                      warmup = 1000,
                      iter = 2000,
                      cores = 4,
                      chains = 4,
                      seed = 123,
                      control = list(adapt_delta = 0.99,
                                     max_treedepth = 20)) # increased treedepth
save(stat_inter, file = here(fit.dir, "stat_inter.RData"))


## Load model
load(here(fit.dir, "stat_inter.RData"), verbose=T)

## Diagnostic plots
pdf(here(fig.dir, "stat_inter_diag.pdf"), width = 7, height = 5)
coda_neff(get_neff(stat_inter, pars = pars.stat), total_draws(stat_inter))
coda_rhat(get_rhat(stat_inter, pars = pars.stat))
coda_diag(As.mcmc.list(stat_inter, pars = pars.stat))
dev.off()

plot_post_pc(stat_inter, stan.dat.all$y, pdf.path = here(fig.dir, "stat_inter_yrep.pdf")) # Working again?

loo.stat_inter <- rstan::loo(stat_inter, cores = 4)
save(loo.stat_inter, file = here(diag.dir, "loo_stat_inter.RData"))
waic.stat_inter <- loo::waic(loo::extract_log_lik(stat_inter, "log_lik"))
save(waic.stat_inter, file = here(diag.dir, "waic_stat_inter.RData"))
pdf(here(fig.dir, "stat_inter_loo.pdf"), width = 7, height = 5)
plot(loo.stat_inter, label_points = TRUE)
dev.off()

pdf( here(fig.dir, "stat_inter_resid.pdf"), width = 8, height = 8)
plot_hbm_resids(stat_inter, data_master)
dev.off()


## Density plot

# Get data
lst <- hb07_density_df(stat_inter, ocean.regions=ifelse(speciesFlag=="chum", 3, 4))
s.df <- lst$stock
m.df <- lst$region
m.df$region <- factor(m.df$region, levels = c("West Coast", "Gulf of Alaska", "Southeast Alaska", "Bering Sea"))

## Covariate labels
vars <- data.frame(var = levels(m.df$var))
vars$lab <- gsub("Comp", "Competitors", vars$var)

# Make density plot
g <- m.df %>% filter(var != "SST + Comp") %>%
  ggplot() +
  geom_vline(xintercept = 0, color = "grey50", linetype = 2, linewidth = 0.25) +
  geom_path(data = s.df[s.df$var != "SST + Comp",],
            aes(x = x, y = y, group = stock, color = region), alpha=0.3,
            na.rm = TRUE) +
  geom_path(aes(x = x, y = y, color = region), linewidth = 1, alpha=1,
            na.rm = TRUE) +
  scale_colour_manual(name = "Ocean Region", values=col.region, guide="legend") +
  labs(x = "Percent change in R/S",
       y = "Posterior density",
       color = "") +
  scale_x_continuous(limits = c(-50, 50), expand = c(0, 0)) +
  scale_y_continuous(breaks=NULL) +
  geom_text(data = vars[1:3,],
            aes(x = -48.1,
                y = max(m.df$y) - 0.008,
                label = lab),
            hjust = 0,
            size = 3.5,
            color = "grey30") +
  facet_wrap( vars(var), ncol = 1) +
  theme_sleek(base_size = 9) +
  theme(legend.justification = c(0, 0),
        legend.position = c(0.75, 0.75),
        legend.key.size = unit(10, "pt"),
        legend.background = element_blank(),
        legend.text = element_text(size = 9),
        panel.spacing.y = unit(-0.5, "pt"),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.title = element_text(size=12),
        axis.text.x = element_text(size=10),
        # make background transparent
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA)
        )

png(file=here('hannahland', 'BC-plots-Jan24', 'stat_inter_dens.png'), width=500*4, height=400*4, res=72*5, bg='transparent')
print(g)
dev.off()

pdf( here(fig.dir, "stat_inter_dens.pdf"))
print(g)
dev.off()


# Boxplots instead of density
m.df %>% filter(var=="Comp") %>%
  ggplot() +
  geom_hline(aes(yintercept=0), linetype="dashed", colour="gray50") +
  geom_vline(aes(xintercept=0), linetype="dashed", colour="gray50") +
  geom_boxplot(aes(x=region, y=x, col=region, fill=region),
               position=position_dodge2(preserve="single"), alpha=0.4, linewidth=0.80) +
  #scale_y_continuous(limits=c(-50,50), labels=seq(-50,50,25)) +
  coord_cartesian(ylim=c(-50,150)) +
  scale_fill_manual(values=sp.col, aesthetics = c("fill", "col")) +
  labs(y="Percent change in R/S", x="", col="Species", fill="Species") +
  theme_sleek() +
  theme(legend.position="right",
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA))

