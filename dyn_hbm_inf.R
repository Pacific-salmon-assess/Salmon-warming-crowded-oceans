## Bayesian model inference

if(!dir.exists("./figures/hbm_inf/"))
    dir.create("./figures/hbm_inf/")

## Load model fits
for(i in list.files(path = "./output/hbm_fit/", pattern = "*.RData$")) {
    load(paste("./output/hbm_fit/", i, sep = ""))
}


## Gamma: era ----------------------------------------------

mcmc.hbm2.sst  <- As.mcmc.list(hbm2.sst,  pars = pars.hbm2)
mcmc.hbm2.npgo <- As.mcmc.list(hbm2.npgo,  pars = pars.hbm2)
mcmc.hbm2.pdo  <- As.mcmc.list(hbm2.pdo,  pars = pars.hbm2)
mcmc.hbm3.sst  <- As.mcmc.list(hbm3.sst,  pars = pars.hbm3)
mcmc.hbm3.npgo <- As.mcmc.list(hbm3.npgo, pars = pars.hbm3)
mcmc.hbm3.pdo  <- As.mcmc.list(hbm3.pdo,  pars = pars.hbm3)

lst.hbm2 <- list(SST  = mcmc.hbm2.sst,
                 NPGO = mcmc.hbm2.npgo,
                 PDO  = mcmc.hbm2.pdo)
lst.hbm3 <- list(SST  = mcmc.hbm3.sst,
                 NPGO = mcmc.hbm3.npgo,
                 PDO  = mcmc.hbm3.pdo)


dfl.hbm2 <- lapply(seq_along(lst.hbm2), function(i) {
    mcmc <- lst.hbm2[[i]]
    df.wide <- coda_df(mcmc, grep("mu_gamma", coda::varnames(mcmc),
                                  value = TRUE))
    df.long <- reshape2::melt(df.wide, id.vars = c("chain", "iter"),
                              variable.name = "par")
    # browser()
    info <- data.frame(par = as.character(unique(df.long$par)),
                       ocean_region = c("WC", "GOA", "BS"),
                       era = rep("all", 3))
    df <- plyr::join(df.long, info, by = "par")
    df$var <- names(lst.hbm2)[i]
    return(df)
})

dfl.hbm3 <- lapply(seq_along(lst.hbm3), function(i) {
    mcmc <- lst.hbm3[[i]]
    df.wide <- coda_df(mcmc, grep("mu_gamma", coda::varnames(mcmc),
                                  value = TRUE))
    df.long <- reshape2::melt(df.wide, id.vars = c("chain", "iter"),
                              variable.name = "par")
    info <- data.frame(par = as.character(unique(df.long$par)),
                       ocean_region = rep(c("WC", "GOA", "BS"), 3),
                       era = c(rep("early", 3),
                               rep("middle", 3),
                               rep("late", 3)))
    df <- plyr::join(df.long, info, by = "par")
    df$var <- names(lst.hbm3)[i]
    return(df)
})

df.hbm2.3 <- rbind(plyr::rbind.fill(dfl.hbm3), plyr::rbind.fill(dfl.hbm2))
df.hbm2.3$era <- factor(df.hbm2.3$era, levels = unique(df.hbm2.3$era))
df.hbm2.3 <- ocean_region_lab(df.hbm2.3)

hbm.gamma.era <- df.hbm2.3
save(hbm.gamma.era, file = "./output/hbm_gamma_era.RData")




## Stock specific era
dfl.hbm2.st <- lapply(seq_along(lst.hbm2), function(i) {
    mcmc <- lst.hbm2[[i]]
    df.wide <- coda_df(mcmc, grep("^gamma", coda::varnames(mcmc),
                                  value = TRUE))
    df.long <- reshape2::melt(df.wide, id.vars = c("chain", "iter"),
                              variable.name = "par")
    info <- sock.info[ , c("stock_no", "stock", "ocean_region")]
    info$par <- paste0("gamma[", info$stock_no, "]")
    info$era <- "all"
    info <- ocean_region_lab(info)
    df <- plyr::join(info, df.long, by = "par")
    df$var <- names(lst.hbm2)[i]
    return(df)
})

dfl.hbm3.st <- lapply(seq_along(lst.hbm3), function(i) {
    mcmc <- lst.hbm3[[i]]
    df.wide <- coda_df(mcmc, grep("^gamma", coda::varnames(mcmc),
                                  value = TRUE))
    df.long <- reshape2::melt(df.wide, id.vars = c("chain", "iter"),
                              variable.name = "par")
    info <- sock.info[ , c("stock_no", "stock", "ocean_region")]
    info1 <- info; info2 <- info; info3 <- info
    info1$par <- paste0("gamma1[", info1$stock_no, "]")
    info2$par <- paste0("gamma2[", info2$stock_no, "]")
    info3$par <- paste0("gamma3[", info3$stock_no, "]")
    info1$era <- "early"
    info2$era <- "middle"
    info3$era <- "late"
    info <- rbind(info1, info2, info3)
    info <- ocean_region_lab(info)
    df <- plyr::join(info, df.long, by = "par")
    df$var <- names(lst.hbm3)[i]
    return(df)
})

df.hbm2.3.st <- rbind(plyr::rbind.fill(dfl.hbm3.st), plyr::rbind.fill(dfl.hbm2.st))
df.hbm2.3.st$era <- factor(df.hbm2.3.st$era, levels = unique(df.hbm2.3.st$era))

hbm.gamma.era.stock <- df.hbm2.3.st
save(hbm.gamma.era.stock, file = "./output/hbm_gamma_era_stock.RData")



## Gamma s2: era -------------------------------------------

mcmc.hbm2.sst.s2  <- As.mcmc.list(hbm2.sst.s2,  pars = pars.hbm2)
mcmc.hbm2.npgo.s2 <- As.mcmc.list(hbm2.npgo.s2,  pars = pars.hbm2)
mcmc.hbm2.pdo.s2  <- As.mcmc.list(hbm2.pdo.s2,  pars = pars.hbm2)
mcmc.hbm3.sst.s2  <- As.mcmc.list(hbm3.sst.s2,  pars = pars.hbm3)
mcmc.hbm3.npgo.s2 <- As.mcmc.list(hbm3.npgo.s2, pars = pars.hbm3)
mcmc.hbm3.pdo.s2  <- As.mcmc.list(hbm3.pdo.s2,  pars = pars.hbm3)

lst.hbm2.s2 <- list(SST  = mcmc.hbm2.sst.s2,
                    NPGO = mcmc.hbm2.npgo.s2,
                    PDO  = mcmc.hbm2.pdo.s2)
lst.hbm3.s2 <- list(SST  = mcmc.hbm3.sst.s2,
                    NPGO = mcmc.hbm3.npgo.s2,
                    PDO  = mcmc.hbm3.pdo.s2)


dfl.hbm2.s2 <- lapply(seq_along(lst.hbm2.s2), function(i) {
    mcmc <- lst.hbm2.s2[[i]]
    df.wide <- coda_df(mcmc, grep("mu_gamma", coda::varnames(mcmc),
                                  value = TRUE))
    df.long <- reshape2::melt(df.wide, id.vars = c("chain", "iter"),
                              variable.name = "par")
    # browser()
    info <- data.frame(par = as.character(unique(df.long$par)),
                       ocean_region = c("WC", "GOA", "BS"),
                       era = rep("all", 3))
    df <- plyr::join(df.long, info, by = "par")
    df$var <- names(lst.hbm2.s2)[i]
    return(df)
})

dfl.hbm3.s2 <- lapply(seq_along(lst.hbm3.s2), function(i) {
    mcmc <- lst.hbm3.s2[[i]]
    df.wide <- coda_df(mcmc, grep("mu_gamma", coda::varnames(mcmc),
                                  value = TRUE))
    df.long <- reshape2::melt(df.wide, id.vars = c("chain", "iter"),
                              variable.name = "par")
    info <- data.frame(par = as.character(unique(df.long$par)),
                       ocean_region = rep(c("WC", "GOA", "BS"), 3),
                       era = c(rep("early", 3),
                               rep("middle", 3),
                               rep("late", 3)))
    df <- plyr::join(df.long, info, by = "par")
    df$var <- names(lst.hbm3.s2)[i]
    return(df)
})

df.hbm2.3.s2 <- rbind(plyr::rbind.fill(dfl.hbm3.s2), plyr::rbind.fill(dfl.hbm2.s2))
df.hbm2.3.s2$era <- factor(df.hbm2.3.s2$era, levels = unique(df.hbm2.3.s2$era))
df.hbm2.3.s2 <- ocean_region_lab(df.hbm2.3.s2)

hbm.gamma.era.s2 <- df.hbm2.3.s2
save(hbm.gamma.era.s2, file = "./output/hbm_gamma_era_s2.RData")



## Stock specific era
dfl.hbm2.st.s2 <- lapply(seq_along(lst.hbm2.s2), function(i) {
    mcmc <- lst.hbm2.s2[[i]]
    df.wide <- coda_df(mcmc, grep("^gamma", coda::varnames(mcmc),
                                  value = TRUE))
    df.long <- reshape2::melt(df.wide, id.vars = c("chain", "iter"),
                              variable.name = "par")
    info <- sock.info[ , c("stock_no", "stock", "ocean_region")]
    info$par <- paste0("gamma[", info$stock_no, "]")
    info$era <- "all"
    info <- ocean_region_lab(info)
    df <- plyr::join(info, df.long, by = "par")
    df$var <- names(lst.hbm2.s2)[i]
    return(df)
})

dfl.hbm3.st.s2 <- lapply(seq_along(lst.hbm3.s2), function(i) {
    mcmc <- lst.hbm3.s2[[i]]
    df.wide <- coda_df(mcmc, grep("^gamma", coda::varnames(mcmc),
                                  value = TRUE))
    df.long <- reshape2::melt(df.wide, id.vars = c("chain", "iter"),
                              variable.name = "par")
    info <- sock.info[ , c("stock_no", "stock", "ocean_region")]
    info1 <- info; info2 <- info; info3 <- info
    info1$par <- paste0("gamma1[", info1$stock_no, "]")
    info2$par <- paste0("gamma2[", info2$stock_no, "]")
    info3$par <- paste0("gamma3[", info3$stock_no, "]")
    info1$era <- "early"
    info2$era <- "middle"
    info3$era <- "late"
    info <- rbind(info1, info2, info3)
    info <- ocean_region_lab(info)
    df <- plyr::join(info, df.long, by = "par")
    df$var <- names(lst.hbm3.s2)[i]
    return(df)
})

df.hbm2.3.st.s2 <- rbind(plyr::rbind.fill(dfl.hbm3.st.s2), plyr::rbind.fill(dfl.hbm2.st.s2))
df.hbm2.3.st.s2$era <- factor(df.hbm2.3.st.s2$era, levels = unique(df.hbm2.3.st.s2$era))

hbm.gamma.era.stock.s2 <- df.hbm2.3.st.s2
save(hbm.gamma.era.stock.s2, file = "./output/hbm_gamma_era_stock_s2.RData")



## Gamma: different ----------------------------------------

lst.hbm5 <- list(SST  = hbm5.sst,
                 NPGO = hbm5.npgo,
                 PDO  = hbm5.pdo)

dfl.hbm5 <- lapply(seq_along(lst.hbm5), function(i) {
    probs <- c(0.025, 0.05, 0.10, 0.50, 0.90, 0.95, 0.975)
    summ <- summary(lst.hbm5[[i]], pars = "gamma", probs = probs)[[1]]
    df <- data.frame(stock = sock.covar$stock,
                     ocean_region = sock.covar$ocean_region,
                     brood_yr = sock.covar$brood_yr,
                     gamma = summ[ , "50%"],
                     lower = summ[ , "10%"],
                     upper = summ[ , "90%"],
                     # lower = summ[ , "5%"],
                     # upper = summ[ , "95%"],
                     var = names(lst.hbm5)[i])
    return(df)
})

df.hbm5 <- plyr::rbind.fill(dfl.hbm5)
df.hbm5 <- ocean_region_lab(df.hbm5)
df.hbm5 <- plyr::ddply(df.hbm5, .(ocean_region_lab, var, stock), transform,
                       sig = (lower < 0 & upper < 0) |
                             (lower > 0 & upper > 0))
df.hbm5$var <- factor(df.hbm5$var, levels = c("NPGO", "PDO", "SST"))

hbm.gamma.diff <- df.hbm5
save(hbm.gamma.diff, file = "./output/hbm_gamma_diff.RData")



## Gamma s1: different -------------------------------------

lst.hbm5.s1 <- list(SST  = hbm5.sst.s1,
                    NPGO = hbm5.npgo.s1,
                    PDO  = hbm5.pdo.s1)

dfl.hbm5.s1 <- lapply(seq_along(lst.hbm5.s1), function(i) {
    probs <- c(0.025, 0.05, 0.10, 0.50, 0.90, 0.95, 0.975)
    summ <- summary(lst.hbm5.s1[[i]], pars = "gamma", probs = probs)[[1]]
    df <- data.frame(stock = sock.covar$stock,
                     ocean_region = sock.covar$ocean_region,
                     brood_yr = sock.covar$brood_yr,
                     gamma = summ[ , "50%"],
                     lower = summ[ , "10%"],
                     upper = summ[ , "90%"],
                     var = names(lst.hbm5.s1)[i])
    return(df)
})

df.hbm5.s1 <- plyr::rbind.fill(dfl.hbm5.s1)
df.hbm5.s1 <- ocean_region_lab(df.hbm5.s1)
df.hbm5.s1 <- plyr::ddply(df.hbm5.s1, .(ocean_region_lab, var, stock), transform,
                          sig = any((lower < 0 & upper < 0) |
                                    (lower > 0 & upper > 0)))
df.hbm5.s1$var <- factor(df.hbm5.s1$var, levels = c("NPGO", "PDO", "SST"))

hbm.gamma.diff.s1 <- df.hbm5.s1
save(hbm.gamma.diff.s1, file = "./output/hbm_gamma_diff_s1.RData")



## Gamma: same ---------------------------------------------

lst.hbm6 <- list(SST  = hbm6.sst,
                 NPGO = hbm6.npgo,
                 PDO  = hbm6.pdo)

dfl.hbm6 <- lapply(seq_along(lst.hbm6), function(i) {
    probs <- c(0.025, 0.05, 0.10, 0.50, 0.90, 0.95, 0.975)
    summ <- summary(lst.hbm6[[i]], pars = "gamma", probs = probs)[[1]]
    dat <- sock.covar
    dat$ocean_region <- factor(dat$ocean_region,
                               levels = unique(dat$ocean_region))
    df <- plyr::ddply(dat, .(ocean_region), function(x) {
                  brood_yr <- min(x$brood_yr):max(x$brood_yr)
                  ocean_region <- rep(unique(x$ocean_region),
                                      length(unique(x$brood_yr)))
                  data.frame(brood_yr = brood_yr,
                             ocean_region = ocean_region) })
    df$gamma <- summ[ , "50%"]
    df$lower <- summ[ , "5%"]
    df$upper <- summ[ , "95%"]
    df$var <- names(lst.hbm6)[i]
    return(df)
})

df.hbm6 <- plyr::rbind.fill(dfl.hbm6)
df.hbm6 <- ocean_region_lab(df.hbm6)
df.hbm6$var <- factor(df.hbm6$var, levels = c("NPGO", "PDO", "SST"))

hbm.gamma.same <- df.hbm6
save(hbm.gamma.same, file = "./output/hbm_gamma_same.RData")



## Graphics ------------------------------------------------

## Base model
plot_hbm_dot(hbm1, pdf.file = "./figures/hbm_inf/hbm1_dot.pdf")
plot_hbm_dens(hbm1, "./figures/hbm_inf/hbm1_sst_dens.pdf")


## SST
plot_hbm_dot(hbm2.sst, pdf.file = "./figures/hbm_inf/hbm2_sst_dot.pdf")
plot_hbm_dot(hbm3.sst, pdf.file = "./figures/hbm_inf/hbm3_sst_dot.pdf")
plot_hbm_dot(hbm4.sst, pdf.file = "./figures/hbm_inf/hbm4_sst_dot.pdf")
plot_hbm_dot(hbm5.sst, pdf.file = "./figures/hbm_inf/hbm5_sst_dot.pdf")
plot_hbm_dot(hbm6.sst, pdf.file = "./figures/hbm_inf/hbm6_sst_dot.pdf")

plot_hbm_dens(hbm2.sst, "./figures/hbm_inf/hbm2_sst_dens.pdf")
plot_hbm_dens(hbm3.sst, "./figures/hbm_inf/hbm3_sst_dens.pdf")
plot_hbm_dens(hbm4.sst, "./figures/hbm_inf/hbm4_sst_dens.pdf")
plot_hbm_dens(hbm5.sst, "./figures/hbm_inf/hbm5_sst_dens.pdf")
plot_hbm_dens(hbm6.sst, "./figures/hbm_inf/hbm6_sst_dens.pdf")


## NPGO
plot_hbm_dot(hbm2.npgo, pdf.file = "./figures/hbm_inf/hbm2_npgo_dot.pdf")
plot_hbm_dot(hbm3.npgo, pdf.file = "./figures/hbm_inf/hbm3_npgo_dot.pdf")
plot_hbm_dot(hbm4.npgo, pdf.file = "./figures/hbm_inf/hbm4_npgo_dot.pdf")
plot_hbm_dot(hbm5.npgo, pdf.file = "./figures/hbm_inf/hbm5_npgo_dot.pdf")
plot_hbm_dot(hbm6.npgo, pdf.file = "./figures/hbm_inf/hbm6_npgo_dot.pdf")

plot_hbm_dens(hbm2.npgo, "./figures/hbm_inf/hbm2_npgo_dens.pdf")
plot_hbm_dens(hbm3.npgo, "./figures/hbm_inf/hbm3_npgo_dens.pdf")
plot_hbm_dens(hbm4.npgo, "./figures/hbm_inf/hbm4_npgo_dens.pdf")
plot_hbm_dens(hbm5.npgo, "./figures/hbm_inf/hbm5_npgo_dens.pdf")
plot_hbm_dens(hbm6.npgo, "./figures/hbm_inf/hbm6_npgo_dens.pdf")


## PDO
plot_hbm_dot(hbm2.pdo, pdf.file = "./figures/hbm_inf/hbm2_pdo_dot.pdf")
plot_hbm_dot(hbm3.pdo, pdf.file = "./figures/hbm_inf/hbm3_pdo_dot.pdf")
plot_hbm_dot(hbm4.pdo, pdf.file = "./figures/hbm_inf/hbm4_pdo_dot.pdf")
plot_hbm_dot(hbm5.pdo, pdf.file = "./figures/hbm_inf/hbm5_pdo_dot.pdf")
plot_hbm_dot(hbm6.pdo, pdf.file = "./figures/hbm_inf/hbm6_pdo_dot.pdf")

plot_hbm_dens(hbm2.pdo, "./figures/hbm_inf/hbm2_pdo_dens.pdf")
plot_hbm_dens(hbm3.pdo, "./figures/hbm_inf/hbm3_pdo_dens.pdf")
plot_hbm_dens(hbm4.pdo, "./figures/hbm_inf/hbm4_pdo_dens.pdf")
plot_hbm_dens(hbm5.pdo, "./figures/hbm_inf/hbm5_pdo_dens.pdf")
plot_hbm_dens(hbm6.pdo, "./figures/hbm_inf/hbm6_pdo_dens.pdf")
