## (8) Recreate the tables and figures in the publication ##
## ------------------------------------------------------ ##


if(!dir.exists("./pub/")) dir.create("./pub/")

## Which model fit ?
fit <- hb07a
fitnam <- "hb07a"


## Define colors
col.stock  <- rev(chroma::qpal(7, alpha = 0.4)[c(1, 4, 6)])
col.region <- rev(chroma::qpal(7, luminance = 40)[c(1, 4, 6)])
col.lt <- rev(chroma::qpal(7)[c(1, 4, 6)])
col.dk <- rev(chroma::qpal(7, luminance = 30)[c(1, 4, 6)])

## Table: coefficients 

gamma <- summary(fit, pars = "mu_gamma")$summary
kappa <- summary(fit, pars = "mu_kappa")$summary
chi   <- summary(fit, pars = "mu_chi")$summary
reg   <- c("West Coast", "Gulf of Alaska", "Bering Sea")

tab.g <- data.frame(reg = reg,
                    coef = "SST",
                    lower = gamma[ , "2.5%"],
                    mean = gamma[ , "mean"],
                    upper = gamma[ , "97.5%"])
tab.k <- data.frame(reg = reg,
                    coef = "Comp",
                    lower = kappa[ , "2.5%"],
                    mean = kappa[ , "mean"],
                    upper = kappa[ , "97.5%"])
tab.c <- data.frame(reg = reg,
                    coef = "SST x Comp",
                    lower = chi[ , "2.5%"],
                    mean = chi[ , "mean"],
                    upper = chi[ , "97.5%"])
tab.coef <- rbind(tab.g, tab.k, tab.c)
tab.coef$perc <- (exp(tab.coef$mean) - 1) * 100
row.names(tab.coef) <- NULL
names(tab.coef) <- c("Ecosystem", "Coefficient", "Lower 95% CI", "Mean",
                    "Upper 95% CI", "Mean % change in R/S")

write.csv(tab.coef, file = paste0("./pub/model_coefficients_", fitnam,  ".csv"))

## Tables for SI ------------------------------------------- 

## Stock-specific params and CIs of main model 
probs <- c(0.025, 0.975)
alpha    <- as.data.frame(summary(fit, pars = "alpha", probs = probs)[[1]])
beta     <- as.data.frame(summary(fit, pars = "beta", probs = probs)[[1]])
sigma    <- as.data.frame(summary(fit, pars = "sigma", probs = probs)[[1]])
sst      <- as.data.frame(summary(fit, pars = "gamma", probs = probs)[[1]])
comp     <- as.data.frame(summary(fit, pars = "kappa", probs = probs)[[1]])
interact <- as.data.frame(summary(fit, pars = "chi", probs = probs)[[1]])

stock.tbl <- data.frame(stock = sock.info$Stock,
                        alpha_mean = alpha[["mean"]],
                        alpha_lowerCI = alpha[["2.5%"]],
                        alpha_upperCI = alpha[["97.5%"]],
                        beta_mean = beta[["mean"]],
                        beta_lowerCI = beta[["2.5%"]],
                        beta_upperCI = beta[["97.5%"]],
                        sst_mean = sst[["mean"]],
                        sst_lowerCI = sst[["2.5%"]],
                        sst_upperCI = sst[["97.5%"]],
                        comp_mean = comp[["mean"]],
                        comp_lowerCI = comp[["2.5%"]],
                        comp_upperCI = comp[["97.5%"]],
                        interact_mean = interact[["mean"]],
                        interact_lowerCI = interact[["2.5%"]],
                        interact_upperCI = interact[["97.5%"]],
                        sigma_mean = sigma[["mean"]],
                        sigma_lowerCI = sigma[["2.5%"]],
                        sigma_upperCI = sigma[["97.5%"]])
write.csv(stock.tbl, file = paste0("./pub/si_main_model_stock_params_", fitnam,  ".csv"))



## Fig: Map + covars ---------------------------------------

pdf(paste0("./pub/fig_map_covars_", fitnam, ".pdf"),
     width = 5, height = 6.0)

## hi-res map
cl <- rnaturalearth::ne_states(country = c("United States of America",
                                           "Mexico",
                                           "Canada"))
namerica.state.sp <- sp::spTransform(cl, sp::CRS("EPSG:4326"))


## Define graphic layout
m <- rbind(c(1, 1, 1, 1, 1, 1),
           c(1, 1, 1, 1, 1, 1),
           c(1, 1, 1, 1, 1, 1),
           c(2, 2, 3, 3, 4, 4),
           c(5, 5, 6, 6, 7, 7),
           c(8, 8, 9, 9, 10, 10))
layout(m)
# print(m)
# layout.show(10)


## Set colors

oc.reg     <- unique(sock.info$Ocean.Region)
rep.stocks <- c("Chilko", "Copper", "Wood")
p.title    <- c("West Coast", "Gulf of Alaska", "Bering Sea")
f.ind      <- 1

par(ps = 8,
    mar = c(0,0,0,0),
    oma = c(4, 4, 0.8, 0.5))


## ~~~~~~~~~~~~~~~~~~~~~~
## Map
## ~~~~~~~~~~~~~~~~~~~~~~

## Setup map labels
labs <- sock.info
labs$no <- 1:nrow(sock.info)
labs <- plyr::ddply(labs, .(lon, lat), summarize,
                    min.no = min(no),
                    max.no = max(no),
                    lab = ifelse(min.no == max.no, min.no,
                                 paste0(min.no, "-", max.no)),
                    x.off = ifelse(min.no == max.no, 1.1, 1.6),
                    y.off = 0.27)
labs$use <- 1
labs$use <- ifelse(labs$min.no == 26, 0, labs$use)
labs$use <- ifelse(labs$min.no == 29, 0, labs$use)
labs$use <- ifelse(labs$min.no == 30, 0, labs$use)
labs$use <- ifelse(labs$min.no == 43, 0, labs$use)
labs$use <- ifelse(labs$min.no == 44, 0, labs$use)
labs$use <- ifelse(labs$min.no == 41, 0, labs$use)
labs <- labs[labs$use == 1, ]
labs$y.off <- ifelse(labs$min.no == 36, -0.1, labs$y.off)

sock.info$stock_no <- 1:nrow(sock.info)
dat.map <- ocean_region_lab(sock.info, factor = FALSE)
labs <- plyr::ddply(dat.map, .(ocean_region_lab, lon, lat), summarize,
                  n = length(stock_no),
                  no_min = min(stock_no),
                  no_max = max(stock_no))
labs$label <- ifelse(labs$no_min == labs$no_max, labs$no_min, paste0(labs$no_min, "-", labs$no_max))
labs$x <- labs$lon + 0.6
labs$y <- labs$lat + 0.5
labs$x[labs$no_min == 20] <- labs$lon[labs$no_min == 20] + 0.95
labs$y[labs$no_min == 20] <- labs$lat[labs$no_min == 20] - 0.1
labs$x[labs$no_min == 26] <- labs$lon[labs$no_min == 26] + 0.95
labs$x[labs$no_min == 27] <- labs$lon[labs$no_min == 27] + 0.95
labs$y[labs$no_min == 26] <- labs$lat[labs$no_min == 26] - 0.1
labs$y[labs$no_min == 27] <- labs$lat[labs$no_min == 27] - 0.1
labs$y[labs$no_min == 28] <- labs$lat[labs$no_min == 28] - 0.65
labs$y[labs$no_min == 30] <- labs$lat[labs$no_min == 30] - 0.3
labs$x[labs$no_min == 30] <- labs$lon[labs$no_min == 30] + 1.7
labs$y[labs$no_min == 35] <- labs$lat[labs$no_min == 35] - 0.6
labs$x[labs$no_min == 37] <- labs$lon[labs$no_min == 37] - 0.8
labs$y[labs$no_min == 37] <- labs$lat[labs$no_min == 37] + 0.3
labs$y[labs$no_min == 39] <- labs$lat[labs$no_min == 39] + 0
labs$x[labs$no_min == 39] <- labs$lon[labs$no_min == 39] + 1
labs$y[labs$no_min == 40] <- labs$lat[labs$no_min == 40] + 0
labs$x[labs$no_min == 40] <- labs$lon[labs$no_min == 40] + 1
labs$x[labs$no_min == 41] <- labs$lon[labs$no_min == 41] + 0.8
## labs$x[labs$no_min == 45] <- labs$lon[labs$no_min == 42] + 0
labs$y[labs$no_min == 45] <- labs$lat[labs$no_min == 45] + 0.8
labs$label[labs$no_min == 29] <- ""
labs$label[labs$no_min == 28] <- "28-29"
labs$label[labs$no_max == 34] <- ""
labs$label[labs$no_max == 32] <- ""
labs$label[labs$no_min == 30] <- "30-34"
labs$label[labs$no_min == 42] <- ""
labs$label[labs$no_min == 43] <- ""
labs$label[labs$no_min == 45] <- ""
labs$label[labs$no_min == 44] <- "42-45"


## Plot map
pch <- c(21, 22, 24)
par(mar = c(5, 0, 0, 0))
plot(sock.info$lon, sock.info$lat,
     type = "n",
     xlab = "",
     ylab = "",
     xlim = c(-165, -120),
     ylim = c(47, 61),
     col = "red",
     axes = FALSE)
box(col = "grey65", lwd = 0.6)
axis(side = 1, lwd = 0, lwd.tick = 1, col = "grey65",
     padj = -1.9, tck = -0.02)
axis(side = 2, lwd = 0, lwd.tick = 1, las = 1, col = "grey65",
     hadj = 0.3, tck = -0.02)
mtext(expression(paste("Longitude (", degree, "E)")),
      side = 1, line = 1.8, cex = 1.0, col = "grey25")
mtext(expression(paste("Latitude (", degree, "N)")),
      side = 2, line = 1.8, cex = 1.0, col = "grey25")
plot(namerica.state.sp, add = TRUE, border = "grey50", lwd = 0.5)
add.label("a", font = 2, cex = 2, xfrac = -0.1, yfrac = 0, xpd = NA)


## Plot stock points
regs <- unique(sock.info$Ocean.Region)
for(i in seq_along(regs)) {
    points(sock.info$lon[sock.info$Ocean.Region == regs[i]],
           sock.info$lat[sock.info$Ocean.Region == regs[i]],
           pch = pch[i], col = col.dk[i], bg = col.lt[i],
           cex = 1.5, lwd = 2)
}

text(x = labs$x, y = labs$y, cex = 0.9, labels = labs$lab)
## text(x = labs$lon + labs$x.off, y = labs$lat + labs$y.off, cex = 0.9,
##      labels = labs$lab)

legend(x = -159, y = 53,
       legend = rev(p.title),
       pch = rev(pch), col = rev(col.dk), pt.bg = rev(col.lt),
       bty = "n", cex = 1.4)

## reset mar
par(mar = c(0, 0, 0, 0))
f.ind <- f.ind + 1


## ~~~~~~~~~~~~~~~~~~~~~~
## Smolt Age-1
## ~~~~~~~~~~~~~~~~~~~~~~
for(i in rev(seq_along(oc.reg))) {

    dat.i <- sock[sock$Ocean.Region == oc.reg[i], ]
    age.mean <- plyr::ddply(dat.i, .(BY), summarize,
                            mean = mean(ocean_1, na.rm = TRUE))

    plot(0, type = "n",
         ylim = c(0, 1),
         xlim = c(1950, 2010),
         ylab = "",
         xlab = "",
         axes = FALSE)
    mtext(p.title[i], col = "grey25")
    abline(v=1975,lty=2,col=" dark grey",lwd=1)
    if(i == 3) {
        mtext("Smolt age-1",
              side = 2, line = 2, cex = 1.0, col = "grey25")
        axis(side = 2, lwd = 0, lwd.tick = 1, col = "grey65", las = 1, hadj = 0.3, tck = -0.04)
        add.label("b", font = 2, cex = 2, xfrac = -0.3, yfrac = -0.1, xpd = NA)
    }
    box(col = "grey65", lwd = 0.6)

    x <- lapply(split(dat.i, dat.i$Stock), function(x) {
                lines(x$BY, x$ocean_1, col = col.stock[i])
         })
    lines(age.mean$BY, age.mean$mean, col = col.dk[i], lwd = 2)
    f.ind <- f.ind + 1
}


## ~~~~~~~~~~~~~~~~~~~~~~
## Covariates
## ~~~~~~~~~~~~~~~~~~~~~~
covars <- c("np_pinks_sec", "early_sst")
scale <- c(100, 1)
at <- list(c(0, 2, 4, 6), -1:1)
ind <- 1
for(k in 1:2 ) {  ## loop over covars

    for(i in rev(seq_along(oc.reg))) {  ## loop over regions

        dat.i <- sock[sock$Ocean.Region == oc.reg[i], ]
        mm <- plyr::ddply(dat.i, .(BY), function(x) {
                          min  <- min(x[ , covars[k]] / scale[k], na.rm = TRUE)
                          max  <- max(x[ , covars[k]] / scale[k], na.rm = TRUE)
                          mean <- mean(x[ , covars[k]] / scale[k], na.rm = TRUE)
                          data.frame(min = min, max = max, mean = mean)
                        })
        if(k == 1)
            ylim <- c(0, max(sock[ , covars[k]] / scale[k], na.rm = TRUE))
        if(k == 2)
            ylim <- range(sock[ , covars[k]] / scale[k], na.rm = TRUE)
        plot(0, type = "n",
             xlim = c(1950, 2010),
             ylim = ylim,
             axes = FALSE)
        abline(v=1975,lty=2,col="dark grey",lwd=1)
        if(i == 3) {
            axis(side = 2, lwd = 0, lwd.tick = 1, las = 1,
                 col = "grey65", at = at[[k]], hadj = 0, tck = -0.04)
        }
        if(k == 2) {
            axis(side = 1, lwd = 0, lwd.tick = 1, las = 1, col = "grey65",
                 at = c(1960, 1980, 2000), padj = -1.9, tck = -0.04)
        }
        if(ind == 5) {
            mtext("Brood year", side = 1, line = 1.8, cex = 1.0, col = "grey25")
        }
        if(ind == 1) {
            mtext("Pink index", side = 2, line = 2, cex = 1.0, col = "grey25")
        }
        if(ind == 4) {
            mtext("SST index", side = 2, line = 2, cex = 1.0, col = "grey25")
        }
        box(col = "grey65", lwd = 0.6)

        x <- lapply(split(dat.i, dat.i$Stock), function(x) {
                    lines(x$BY, x[ , covars[k]] / scale[k], col = col.stock[i])
             })
        if(k == 1) {
            x <- lapply(split(dat.i, dat.i$Stock), function(x) {
                        lines(x$BY, x[ , "np_pinks_sec"] / scale[k], col = col.stock[i], lty = 2)
            })
            ph <- plyr::ddply(dat.i, .(BY), function(x) {
                min  <- min(x[ , "np_pinks_sec"] / scale[k], na.rm = TRUE)
                max  <- max(x[ , "np_pinks_sec"] / scale[k], na.rm = TRUE)
                mean <- mean(x[ , "np_pinks_sec"] / scale[k], na.rm = TRUE)
                data.frame(min = min, max = max, mean = mean)
            })
            lines(ph$BY, ph$mean, col = col.dk[i], lwd = 2, lty = 2)
        }
        lines(mm$BY, mm$mean, col = col.dk[i], lwd = 2)
        ind <- ind + 1
        f.ind <- f.ind + 1
    }
}



## ~~~~~~~~~~~~~~~~~~~~~~
## Inset map
## ~~~~~~~~~~~~~~~~~~~~~~
# c(x_left, x_right, y_bottom, y_top)
par(fig=c(0.78, 0.99, 0.85, 0.99), mar = c(0,0,0,0), new = TRUE)
plot(1, 1,
     xlim = c(-180, -57),
     ylim = c(12.5, 80),
     col  = "white",
     xlab = "",
     ylab = "",
     axes = FALSE)
inc.cont <- c("United States of America",
              "Canada",
              "Mexico",
              "Guatemala",
              "Honduras",
              "Nicaragua",
              "Costa Rica",
              "Panama")
plot(countriesLow[countriesLow$ADMIN %in% inc.cont, ],
     col = "grey55", border = "grey55", add = TRUE)
rect(xleft = -165, ybottom = 47, xright = -120, ytop = 61,
     lwd = 2, border = 1)


dev.off()


## Fig: Posterior percent change density ------------------- 
lst <- hb07_density_df(fit)
s.df <- lst$stock
m.df <- lst$region

## Covariate labels
vars <- data.frame(var = levels(m.df$var))
vars$lab <- paste0("(", letters[1:nrow(vars)], ") ", vars$var)
g <- ggplot(m.df) +
    geom_vline(xintercept = 0, color = "grey50", linetype = 2, linewidth = 0.25) +
    geom_path(data = s.df[s.df$region == "West Coast", ],
              aes(x = x, y = y, group = stock), color = col.stock[1],
              na.rm = TRUE) +
    geom_path(data = s.df[s.df$region == "Gulf of Alaska", ],
              aes(x = x, y = y, group = stock), color = col.stock[2],
              na.rm = TRUE) +
    geom_path(data = s.df[s.df$region == "Bering Sea", ],
              aes(x = x, y = y, group = stock), color = col.stock[3],
              na.rm = TRUE) +
    geom_path(aes(x = x, y = y, color = region), linewidth = 1,
              na.rm = TRUE) +
    scale_color_manual(values = rev(col.region)) +
    labs(x = "Percent change in R/S",
         y = "Posterior density",
         color = "") +
    scale_x_continuous(limits = c(-50, 50), expand = c(0, 0)) +
    scale_y_continuous(breaks=NULL) +
    geom_text(data = vars,
              aes(x = -48.1,
                  y = max(m.df$y) - 0.008,
                  label = lab),
              hjust = 0,
              size = 2.7,
              color = "grey30") +
    facet_wrap( ~ var, ncol = 1) +
    theme_sleek(base_size = 9) +
    theme(legend.justification = c(0, 0),
          legend.position = c(0.7, 0.91),
          legend.key.size = unit(10, "pt"),
          legend.background = element_blank(),
          legend.text = element_text(size = 8),
          panel.spacing.y = unit(-0.5, "pt"),
          strip.background = element_blank(),
          strip.text.x = element_blank())
print(g)

pdf(paste0("./pub/fig_post_density_", fitnam, ".pdf"), width = 4, height = 6)
    print(g)
dev.off()


## Fig: dot + density main --------------------------------- 
gamma.stock <- hb_param_df(fit, "gamma", "Ocean.Region", "SST")
kappa.stock <- hb_param_df(fit, "kappa", "Ocean.Region", "Comp")
chi.stock   <- hb_param_df(fit, "chi", "Ocean.Region", "SST x Comp")
df.dot <- rbind(gamma.stock, kappa.stock, chi.stock)
df.dot <- ocean_region_lab(df.dot, "region", FALSE)
df.dot$var <- factor(df.dot$var, levels = c("SST", "Comp", "SST x Comp"))
df.mu <- plyr::ddply(df.dot, .(region, var), summarize,
                     mu_mean = unique(mu_mean),
                     mu_2.5 = unique(`mu_2.5%`),
                     mu_97.5 = unique(`mu_97.5%`),
                     ocean_region_lab = unique(ocean_region_lab),
                     ystart = Stock[1],
                     yend = Stock[length(Stock)])

g <- ggplot(df.dot) +
    geom_vline(xintercept = 0, color = "grey50", linetype = 2, size = 0.25) +
    geom_point(aes(x = mean, y = Stock, color = ocean_region_lab, shape = ocean_region_lab)) +
    geom_segment(aes(y = Stock, yend = Stock, x = `2.5%`, xend = `97.5%`,
                     color = ocean_region_lab), size = 0.25) +
    geom_segment(data = df.mu, aes(y = ystart, yend = yend, x = mu_mean, xend = mu_mean,
                                   color = ocean_region_lab), size = 0.25) +
    geom_rect(data = df.mu, aes(xmin = mu_2.5, xmax = mu_97.5, ymin = ystart,
                                ymax = yend, fill = ocean_region_lab),
              alpha = 0.2) +
    scale_color_manual(values = rev(col.region)) +
    scale_shape_manual(values = c(15, 16, 17)) +
    scale_fill_manual(values = rev(col.region), guide="none") +
    labs(x = "Coefficient",
         y = "",
         color = "",
         shape = "") +
    facet_wrap( ~ var) +
    scale_x_continuous(breaks=c(-0.25,0,0.25))+
    theme_sleek(base_size = 10) +
    theme(legend.justification = c(0, 0),
          legend.position = c(0.01, 0.87),
          legend.key.size = unit(10, "pt"),
          legend.background = element_blank(),
          legend.text = element_text(size = 8),
          panel.spacing.x = unit(-0.5, "pt"))
print(g)

pdf(paste0("./pub/Sfig_coef_dot_main_", fitnam, ".pdf"), width = 6.5, height = 6.0)
    print(g)
dev.off()
