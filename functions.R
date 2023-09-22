## Functions for analysis ##
## ---------------------- ##

## geographic.order --------------------------------------

geographic.order <- function(x) {
  # Accepts a dataframe with the columns Stock, Ocean.Region, Lat, Lon
  # Returns Stock as an ordered factor appropriate for plotting
      # i.e. 1 is southmost stk on WC, max is northmost in BS, and 
                                # GOA are ordered E->W 

  wc.ind <- bs.ind <- seak.ind <- which(names(x) %in% c("Stock", "Lat", "lat")) # WC and BS stocks organized by latitude
  goa.ind <- which(names(x) %in% c("Stock", "Lon", "lon")) # GOA stocks organized by longitude

  # Get stock names in each region
  WC_stk <- unique( x[ which (x$Ocean.Region2=="WC"), wc.ind  ] )
  SEAK_stk <- unique( x[ which (x$Ocean.Region2=="SEAK"), seak.ind  ] )
  GOA_stk <- unique( x[ which (x$Ocean.Region2=="GOA"), goa.ind ] )
  BS_stk <- unique( x[ which (x$Ocean.Region2=="BS"), bs.ind ] )
  
  # Rank by Lat (WC, BS) or Lon (GOA)
  wc.ind <- which(names(WC_stk) %in% c("Lat", "lat"))
  seak.ind <- which(names(SEAK_stk) %in% c("Lat", "lat"))
  goa.ind <- which(names(GOA_stk) %in% c("Lon", "lon"))
  bs.ind <- which(names(BS_stk) %in% c("Lat", "lat"))
  
  
  WC_stk$geo_id <- data.table::frankv(WC_stk, cols=wc.ind, ties.method = "first")
  SEAK_stk$geo_id <- nrow(WC_stk) + data.table::frankv(SEAK_stk, cols=seak.ind, order=-1L, ties.method = "first")
  GOA_stk$geo_id <- nrow(WC_stk) + nrow(SEAK_stk) + data.table::frankv(GOA_stk, cols=goa.ind, order=-1L, ties.method = "first")
  BS_stk$geo_id <- nrow(WC_stk) + nrow(SEAK_stk) + nrow(GOA_stk) + data.table::frankv(BS_stk, cols=bs.ind, ties.method="first")
  
  #combine them
  geo.id <-  rbind(WC_stk[, c("Stock", "geo_id")], SEAK_stk[, c("Stock", "geo_id")], GOA_stk[, c("Stock", "geo_id")], BS_stk[, c("Stock", "geo_id")])
  geo.id <- geo.id[order(geo.id$geo_id), ] # sort in ascending order
  
  x$Stock <- factor(x$Stock, levels = geo.id$Stock, ordered = T)
  
  return(x$Stock)
}



## hb05_density_df -----------------------------------------
hb05_density_df <- function(stanfit, ocean.regions = 3) {

  # this function has been modified a lot to fit specific purposes but needs work
  # to make it more generally applicable
  
  # Get stanfit name 
  fitnam <- deparse(substitute(stanfit))
  
  ## Get posterior matrices
  lst.f <- hb.posterior.list(stanfit)
  
  ## Density smoothness
  adjust <- 1.5
  
  if(ocean.regions == 3) region_col <- sock.info$Ocean.Region 
  else if (ocean.regions == 4) region_col <- sock.info$Ocean.Region2
  
  ## Define region column indices
  if (fitnam == "stat_ctrl"){ #Use 2020 data for control model
    control_dat <- distinct(ctrl_dat[,c("Stock", "Ocean.Region")])
    ind.wc  <- which(control_dat$Ocean.Region == "WC")
    ind.goa <- which(control_dat$Ocean.Region == "GOA")
    ind.bs  <- which(control_dat$Ocean.Region == "BS")
    
  } else if( fitnam == "stat_tr"){
    ind.wc  <- match(sock.info$Stock[which(region_col == "WC" & sock.info$yr_end >= 1975)], levels(sock$Stock) ) 
    ind.goa <- match(sock.info$Stock[which(region_col == "GOA" & sock.info$yr_end >= 1975)], levels(sock$Stock) )
    ind.bs  <- match(sock.info$Stock[which(region_col == "BS" & sock.info$yr_end >= 1975)], levels(sock$Stock) )
    ind.seak <- match(sock.info$Stock[which(region_col == "SEAK" & sock.info$yr_end >= 1975)], levels(sock$Stock))
    ind.reg <- list(ind.bs, ind.goa, ind.wc, ind.seak)
    
  } else {
    ind.wc  <- match(sock.info$Stock[which(region_col == "WC")], levels(sock$Stock) ) 
    ind.goa <- match(sock.info$Stock[which(region_col == "GOA")], levels(sock$Stock) )
    ind.bs  <- match(sock.info$Stock[which(region_col == "BS")], levels(sock$Stock) )
    ind.seak <- match(sock.info$Stock[which(region_col == "SEAK")], levels(sock$Stock))
    ind.reg <- ifelse(ocean.regions == 3, list(ind.bs, ind.goa, ind.wc), list(ind.bs, ind.goa, ind.wc, ind.seak) )
  }
  
  
  ## Extract stock-specific param matrices
  s.alpha <- lst.f[["alpha"]]
  s.sigma <- lst.f[["sigma"]]
  s.gamma <- lst.f[["gamma"]]
  s.kappa <- lst.f[["kappa"]]
  #s.chi   <- lst.f[["chi"]]
  
  ## Extract regional-level param matrices
  m.gamma <- lst.f[["mu_gamma"]]
  m.kappa <- lst.f[["mu_kappa"]]
  #m.chi   <- lst.f[["mu_chi"]]
  
  ## Get percent change values
  pc.s.gamma <- (exp(s.gamma) - 1) * 100
  pc.s.kappa <- (exp(s.kappa) - 1) * 100
  #pc.s.chi   <- (exp(s.chi) - 1) * 100
  
  pc.m.gamma <- (exp(m.gamma) - 1) * 100
  pc.m.kappa <- (exp(m.kappa) - 1) * 100
  #pc.m.chi   <- (exp(m.chi) - 1) * 100
  
  pc.s.joint <- (exp(s.gamma + s.kappa ) - 1) * 100
  pc.m.joint <- (exp(m.gamma + m.kappa ) - 1) * 100
  

  ## Calculate kernel densities
  s.den.wc.gamma  <- col_density(pc.s.gamma[ , ind.wc], plot.it = FALSE, adjust = adjust)
  s.den.goa.gamma <- col_density(pc.s.gamma[ , ind.goa], plot.it = FALSE, adjust = adjust)
  s.den.bs.gamma  <- col_density(pc.s.gamma[ , ind.bs], plot.it = FALSE, adjust = adjust)
  if(ocean.regions == 4) s.den.seak.gamma <- col_density(pc.s.gamma[ , ind.seak], plot.it = FALSE, adjust = adjust)
  s.den.wc.kappa  <- col_density(pc.s.kappa[ , ind.wc], plot.it = FALSE, adjust = adjust)
  s.den.goa.kappa <- col_density(pc.s.kappa[ , ind.goa], plot.it = FALSE, adjust = adjust)
  s.den.bs.kappa  <- col_density(pc.s.kappa[ , ind.bs], plot.it = FALSE, adjust = adjust)
  if(ocean.regions == 4) s.den.seak.kappa <- col_density(pc.s.kappa[ , ind.seak], plot.it = FALSE, adjust = adjust)
  #s.den.wc.chi    <- col_density(pc.s.chi[ , ind.wc], plot.it = FALSE, adjust = adjust)
  #s.den.goa.chi   <- col_density(pc.s.chi[ , ind.goa], plot.it = FALSE, adjust = adjust)
  #s.den.bs.chi    <- col_density(pc.s.chi[ , ind.bs], plot.it = FALSE, adjust = adjust)
  
  s.den.wc.joint  <- col_density(pc.s.joint[ , ind.wc], plot.it = FALSE, adjust = adjust)
  s.den.goa.joint <- col_density(pc.s.joint[ , ind.goa], plot.it = FALSE, adjust = adjust)
  s.den.bs.joint  <- col_density(pc.s.joint[ , ind.bs], plot.it = FALSE, adjust = adjust)
  
  if(ocean.regions == 4)  s.den.seak.joint  <- col_density(pc.s.joint[ , ind.seak], plot.it = FALSE, adjust = adjust)
  
  m.den.gamma <- col_density(pc.m.gamma, plot.it = FALSE, adjust = adjust)
  m.den.kappa <- col_density(pc.m.kappa, plot.it = FALSE, adjust = adjust)
  #m.den.chi   <- col_density(pc.m.chi, plot.it = FALSE, adjust = adjust)
  m.den.joint <- col_density(pc.m.joint, plot.it = FALSE, adjust = adjust)
  
  ## Density data frame for mu
  m.lst <- list(m.den.gamma,
                m.den.kappa,
                #m.den.chi,
                m.den.joint)
  m.lst.df <- lapply(seq_along(m.lst), function(i) {
    x <- m.lst[[i]]
    x.df <- reshape2::melt(x$x, varnames = c("n", "region"), value.name = "x")
    y.df <- reshape2::melt(x$y, varnames = c("n", "region"), value.name = "y")
    if(ocean.regions == 4) {
      x.df$region <- dplyr::case_when(x.df$region == 1 ~ "West Coast",
                                      x.df$region == 2 ~ "Southeast Alaska",
                                      x.df$region == 3 ~ "Gulf of Alaska",
                                      x.df$region == 4 ~ "Bering Sea" )
    } else {
      x.df$region <- dplyr::case_when( x.df$region == 1 ~ "West Coast",
                                       x.df$region == 2 ~ "Gulf of Alaska",
                                       x.df$region == 3 ~ "Bering Sea" )
    }
    
    y.df$region <- x.df$region
    var  <- ifelse(i == 1, "SST", NA)
    var  <- ifelse(i == 2, "Comp", var)
    var  <- ifelse(i == 3, "SST + Comp", var)
    #var  <- ifelse(i == 4, "SST + Comp + SST x Comp", var)
    x.df$var <- var
    y.df$var <- var
    merge(x.df, y.df, by = c("n", "region", "var"), sort = FALSE)
  })
  m.df <- do.call(rbind, m.lst.df)
  
  
  ## Density data frame for mu
  s.wc.lst <- list(s.den.wc.gamma,
                   s.den.wc.kappa,
                   #s.den.wc.chi,
                   s.den.wc.joint)
  s.wc.lst.df <- lapply(seq_along(s.wc.lst), function(i) {
    x <- s.wc.lst[[i]]
    x.df <- reshape2::melt(x$x, varnames = c("n", "stock"), value.name = "x")
    y.df <- reshape2::melt(x$y, varnames = c("n", "stock"), value.name = "y")
    x.df$region <- "West Coast"
    y.df$region <- "West Coast"
    var  <- ifelse(i == 1, "SST", NA)
    var  <- ifelse(i == 2, "Comp", var)
    var  <- ifelse(i == 3, "SST + Comp", var)
    #var  <- ifelse(i == 4, "SST + Comp + SST x Comp", var)
    x.df$var <- var
    y.df$var <- var
    merge(x.df, y.df, by = c("n", "stock", "region", "var"), sort = FALSE)
  })
  s.wc.df <- do.call(rbind, s.wc.lst.df)
  
  s.goa.lst <- list(s.den.goa.gamma,
                    s.den.goa.kappa,
                    #s.den.goa.chi,
                    s.den.goa.joint)
  s.goa.lst.df <- lapply(seq_along(s.goa.lst), function(i) {
    x <- s.goa.lst[[i]]
    x.df <- reshape2::melt(x$x, varnames = c("n", "stock"), value.name = "x")
    y.df <- reshape2::melt(x$y, varnames = c("n", "stock"), value.name = "y")
    x.df$region <- "Gulf of Alaska"
    y.df$region <- "Gulf of Alaska"
    var  <- ifelse(i == 1, "SST", NA)
    var  <- ifelse(i == 2, "Comp", var)
    var  <- ifelse(i == 3, "SST + Comp", var)
    #var  <- ifelse(i == 4, "SST + Comp + SST x Comp", var)
    x.df$var <- var
    y.df$var <- var
    merge(x.df, y.df, by = c("n", "stock", "region", "var"), sort = FALSE)
  })
  s.goa.df <- do.call(rbind, s.goa.lst.df)
  
  s.bs.lst <- list(s.den.bs.gamma,
                   s.den.bs.kappa,
                   #s.den.bs.chi,
                   s.den.bs.joint)
  s.bs.lst.df <- lapply(seq_along(s.bs.lst), function(i) {
    x <- s.bs.lst[[i]]
    x.df <- reshape2::melt(x$x, varnames = c("n", "stock"), value.name = "x")
    y.df <- reshape2::melt(x$y, varnames = c("n", "stock"), value.name = "y")
    x.df$region <- "Bering Sea"
    y.df$region <- "Bering Sea"
    var  <- ifelse(i == 1, "SST", NA)
    var  <- ifelse(i == 2, "Comp", var)
    var  <- ifelse(i == 3, "SST + Comp", var)
    #var  <- ifelse(i == 4, "SST + Comp + SST x Comp", var)
    x.df$var <- var
    y.df$var <- var
    merge(x.df, y.df, by = c("n", "stock", "region", "var"), sort = FALSE)
  })
  s.bs.df <- do.call(rbind, s.bs.lst.df)
  
  if(ocean.regions == 4) {
    s.seak.lst <- list(s.den.seak.gamma,
                     s.den.seak.kappa,
                     #s.den.seak.chi,
                     s.den.seak.joint)
    s.seak.lst.df <- lapply(seq_along(s.seak.lst), function(i) {
      x <- s.seak.lst[[i]]
      x.df <- reshape2::melt(x$x, varnames = c("n", "stock"), value.name = "x")
      y.df <- reshape2::melt(x$y, varnames = c("n", "stock"), value.name = "y")
      x.df$region <- "Southeast Alaska"
      y.df$region <- "Southeast Alaska"
      var  <- ifelse(i == 1, "SST", NA)
      var  <- ifelse(i == 2, "Comp", var)
      var  <- ifelse(i == 3, "SST + Comp", var)
      #var  <- ifelse(i == 4, "SST + Comp + SST x Comp", var)
      x.df$var <- var
      y.df$var <- var
      merge(x.df, y.df, by = c("n", "stock", "region", "var"), sort = FALSE)
    })
    s.seak.df <- do.call(rbind, s.seak.lst.df)
  }

  ## Combine region stock-specific data frames
  if(ocean.regions == 3) s.df <- rbind(s.wc.df, s.goa.df, s.bs.df) else
  if(ocean.regions == 4) s.df <- rbind(s.wc.df, s.goa.df, s.bs.df, s.seak.df)
  
  ## Set factor levels
  s.df$var <- factor(s.df$var, levels = c("SST", "Comp" , "SST + Comp")) #, "SST + Comp + SST x Comp"))
  m.df$var <- factor(m.df$var, levels = c("SST", "Comp" , "SST + Comp")) #, "SST + Comp + SST x Comp"))
  #m.df$region <- factor(m.df$region, levels=c("West Coast", "Gulf of Alaska", "Bering Sea", "Southeast Alaska"))
  #s.df$region <- factor(s.df$region, levels=c("West Coast", "Gulf of Alaska", "Bering Sea", "Southeast Alaska"))
  
  
  out <- list(stock = s.df,
              region = m.df)
  return(out)
}


## era_density_df -----------------------------------------
era_density_df <- function(stanfit, par, region.var="Ocean.Region2", mu = FALSE, percent.change=FALSE, neras=3, info=sock.info){
  
  if(mu){ ## regional mean df
    mu.parxera <- paste0(rep(paste0("mu_", par), each=neras), 1:neras)
    
    post <- rstan::extract(stanfit, pars=c(mu.parxera))
    
    dens.l <- lapply(post, function(x){
      if(percent.change) {
        x<- (exp(x)- 1) * 100
      }
      dens.out <- col_density(x, plot.it=F)
      dens.list <- lapply(dens.out, adply, .margins=c(1,2))
      dens.df <- join(dens.list$x, dens.list$y, by=c("X1", "X2"))
      names(dens.df) <- c("n", "region", "x", "dens")
      dens.df$Ocean.Region2 <- case_when(dens.df$region == 1 ~ "WC",
                                         dens.df$region == 2 ~ "SEAK",
                                         dens.df$region == 3 ~ "GOA",
                                         dens.df$region == 4 ~ "BS")
      return(dens.df)
    } )
    
    summ.dens <- bind_rows(dens.l, .id="par") 
    dens.df.reg.2c <- data.frame(summ.dens,
                                 era = case_when(
                                   str_extract(summ.dens$par, "\\d") == "1" ~ "Early",
                                   str_extract(summ.dens$par, "\\d") == "2" ~ "Middle",
                                   str_extract(summ.dens$par, "\\d") == "3" ~ "Late",
                                   .ptype=factor( levels=c("Early", "Middle", "Late"))),
                                 var = str_extract(summ.dens$par, "\\D+"),
                                 varnam = case_when(
                                   str_extract(summ.dens$par, "\\D+") == "mu_gamma" ~ "SST",
                                   str_extract(summ.dens$par, "\\D+") == "mu_kappa" ~ "Competitors",
                                   .default = NA)
    )
    
    
    
  } else { ## Stock specific df
    parxera <- paste0(rep(par, each=neras), 1:neras)
    post <- rstan::extract(stanfit, pars=parxera)
  dens.l <- lapply(post, function(x){
    if(percent.change) {
      x<- (exp(x) - 1) * 100
    }
    dens.out <- col_density(x, plot.it=F)
    dens.list <- lapply(dens.out, adply, .margins=c(1,2))
    dens.df <- join(dens.list$x, dens.list$y, by=c("X1", "X2"))
    names(dens.df) <- c("n", "stock", "x", "dens")
    dens.df <- mutate(dens.df, stock=levels(info$Stock)[as.numeric(stock)] ) 
    return(dens.df)
  } )
  summ.dens <- bind_rows(dens.l, .id="par") 
  dens.df.st.2c <- data.frame(summ.dens,
                              Ocean.Region2 = info[[region.var]][match(summ.dens$stock, info$Stock)],
                              era = case_when(
                                str_extract(summ.dens$par, "\\d") == "1" ~ "Early",
                                str_extract(summ.dens$par, "\\d") == "2" ~ "Middle",
                                str_extract(summ.dens$par, "\\d") == "3" ~ "Late",
                                .ptype=factor( levels=c("Early", "Middle", "Late"))),
                              var = str_extract(summ.dens$par, "\\D+"),
                              varnam = case_when(
                                str_extract(summ.dens$par, "\\D+") == "gamma" ~ "SST",
                                str_extract(summ.dens$par, "\\D+") == "kappa" ~ "Competitors",
                                .default = NA))
  
  }
}



## ocean_region_lab ----------------------------------------
ocean_region_lab <- function(data, var = "Ocean.Region2", factor = TRUE) {
    ## Add ocean_region_lab column to data.frame
    ##
    ## data = a data.frame with a "Ocean.Region" column

    lab <- rep(NA, nrow(data))
    lab <- ifelse(data[[var]] == "WC", "West Coast", lab)
    lab <- ifelse(data[[var]] == "GOA", "Gulf of Alaska", lab)
    lab <- ifelse(data[[var]] == "BS", "Bering Sea", lab)
    lab <- ifelse(data[[var]] == "SEAK", "Southeast Alaska", lab)
    if(factor)
        lab <- factor(lab, levels = unique(lab))
    data[["ocean_region_lab"]] <- lab
    return(data)
}


## plot_hbm_resids -----------------------------------------
plot_hbm_resids <- function(stanfit, data,
                            pdf.path = NULL,
                            var.yhat = "yhat",
                            var.resid = "yresid") {
    ## Plot realized residuals
    ##
    ## stanfit = stanfit object with posterior samples
    ## data = data.frame of data that model was fit with
    ## pdf.path = file path to save graphics
    ## var.yrep = name of yhat parameter in stanfit
    ## var.yrep = name of resid parameter in stanfit

    yresid <- summary(stanfit, pars = var.resid)$summary
    yhat   <- summary(stanfit, pars = var.yhat)$summary
    df <- data.frame(stock = data[["Stock"]],
                     ocean_region = data[["Ocean.Region2"]],
                     BY = data[["BY"]],
                     y = data[["lnRS"]],
                     yhat = yhat[ , "mean"],
                     yresid = yresid[ , "mean"])
    df$regime <- ifelse(df$BY < 1989, "early", "late")
    row.names(df) <- NULL

    if(!is.null(pdf.path)) {
        pdf(pdf.path, width = 8, height = 8)
    }

    ## resids vs. brood yr (point)
    g <- ggplot(df) +
        geom_hline(yintercept = 0, color = "grey50", linetype = 2) +
        geom_point(aes(x = BY, y = yresid), alpha = 0.5) +
        geom_smooth(aes(x = BY, y = yresid),
                    method = "gam",
                    formula = y ~ s(x, bs = "cs", k = 5),
                    se = FALSE) +
        labs(x = "Brood year",
             y = "Average realized residual") +
        facet_wrap( ~ ocean_region, nrow = 3) +
        theme_sleek()
    print(g)

    ## resids vs. brood yr (line)
    g <- ggplot(df) +
        geom_hline(yintercept = 0, color = "grey50", linetype = 2) +
        geom_line(aes(x = BY, y = yresid, group = stock), alpha = 0.5) +
        labs(x = "Brood year",
             y = "Average realized residual") +
        facet_wrap( ~ ocean_region, nrow = 3) +
        theme_sleek()
    print(g)

    ## resid density by regime
    g <- ggplot(df) +
        geom_vline(xintercept = 0, color = "grey50", linetype = 2) +
        geom_density(aes(x = yresid, color = regime), adjust = 1.5) +
        labs(x = "Average realized residual",
             y = "Density") +
        scale_color_manual(values = c("tomato", "steelblue")) +
        facet_wrap( ~ ocean_region, nrow = 3) +
        theme_sleek()
    print(g)

    ## resids vs. yhat
    g <- ggplot(df) +
        geom_hline(yintercept = 0, color = "grey50", linetype = 2) +
        geom_point(aes(x = yhat, y = yresid), alpha = 0.5) +
        geom_smooth(aes(x = yhat, y = yresid),
                    method = "gam",
                    formula = y ~ s(x, bs = "cs", k = 5),
                    se = FALSE) +
        labs(x = "Averaged fitted y",
             y = "Average realized residual") +
        facet_wrap( ~ ocean_region, nrow = 3) +
        theme_sleek()
    print(g)

    if(!is.null(pdf.path)) {
        dev.off()
    }
}



## bayes_R2 ------------------------------------------------
bayes_R2 <- function(y, ypred) {
    ## Calculate Bayesian R^2 using Gelman et al. (2018) approach
    ## see: https://avehtari.github.io/bayes_R2/bayes_R2.html
    ##
    ## y = vector of observed data
    ## ypred = matrix of predicted values: dims = [iter x N]

    e <- sweep(ypred, 2, y, FUN = "-") * -1  ## residuals
    var_ypred <- apply(ypred, 1, var)        ## predicted variance
    var_e <- apply(e, 1, var)                ## residual variance
    var_ypred / (var_ypred + var_e)          ## R^2
}

if(FALSE) {
    r2 <- bayes_R2(sock$lnRS, as.matrix(hb01a, pars = "yhat"))
    summary(r2)
    hist(r2)
}



## plot_hbm_dot_mu -----------------------------------------
plot_hbm_dot_mu <- function(stanfit, par, region.var, var = NULL,
                            title = NULL, mu = TRUE, base_size = 11) {

    model <- deparse(substitute(stanfit))
    if(is.null(title)) title <- model

    lst <- vector("list", length(par))
    for(i in 1:length(par)) {
        lst[[i]] <- hb_param_df(stanfit, par[i], region.var, var[i])
    }
    df <- plyr::rbind.fill(lst)
    df$var <- factor(df$var, levels = unique(df$var))

    g <- ggplot(df) +
        geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
        geom_point(aes(y = Stock, x = `50%`, color = region)) +
        geom_segment(aes(y = Stock, yend = Stock, x = `2.5%`, xend = `97.5%`,
                         color = region)) +
        labs(x = "Coefficient",
             y = "",
             fill = "",
             color = "") +
        ggtitle(title) +
        theme_sleek(base_size = base_size) +
        scale_color_manual(values = c("#00A0EDFF", "#DE732DFF", "#00AB06FF")) +
        scale_fill_manual(values = c("#00A0EDFF", "#DE732DFF", "#00AB06FF")) +
        facet_wrap( ~ var)

        if(mu) {
            g <- g + geom_segment(aes(y = as.numeric(Stock) - 0.5,
                             yend = as.numeric(Stock) + 0.5,
                             x = `mu_50%`,
                             xend = `mu_50%`,
                             color = region), size = 1) +
            geom_rect(data = df, aes(fill = region),
                      xmin = df[["mu_2.5%"]],
                      xmax = df[["mu_97.5%"]],
                      ymin = as.numeric(df$Stock) - 0.5,
                      ymax = as.numeric(df$Stock) + 0.5, alpha = 0.25)
        }


    print(g)

}



## plot_hbm_dot --------------------------------------------
plot_hbm_dot <- function(stanfit,
                         gamma = FALSE,
                         kappa = FALSE,
                         chi = FALSE,
                         pdf.file = NULL,
                         width = 6,
                         height = 10) {
    ## HMB posterior dot plots
    ##
    ## CI inner = 80%
    ## CI outer = 95%
    ##
    ## stanfit = stanfit model object
    ## pdf.file = path to pdf file to save
    ## width = pdf file width
    ## height = pdf file height

    model <- deparse(substitute(stanfit))

    p.alpha <- mcmc_intervals(as.matrix(stanfit, pars = c("alpha", "mu_alpha")),
                              prob = 0.8, prob_outer = 0.95) +
    ggtitle(model)

    p.beta <- mcmc_intervals(as.matrix(stanfit, pars = "beta"),
                             prob = 0.8, prob_outer = 0.95) +
    ggtitle(model)

    p.sigma <- mcmc_intervals(as.matrix(stanfit, pars = "sigma"),
                              prob = 0.8, prob_outer = 0.95) +
    ggtitle(model)

    if(gamma) {
        p.gamma <- mcmc_intervals(as.matrix(stanfit, pars = c("gamma", "mu_gamma")),
                                  prob = 0.8, prob_outer = 0.95) +
        ggtitle(model)
    }

    if(kappa) {
        p.kappa <- mcmc_intervals(as.matrix(stanfit, pars = c("kappa", "mu_kappa")),
                                  prob = 0.8, prob_outer = 0.95) +
        ggtitle(model)
    }

    if(chi) {
        p.chi <- mcmc_intervals(as.matrix(stanfit, pars = c("chi", "mu_chi")),
                                  prob = 0.8, prob_outer = 0.95) +
        ggtitle(model)
    }


    if(!is.null(pdf.file)) pdf(pdf.file, width = width, height = height)
        print(p.alpha)
        print(p.beta)
        print(p.sigma)
        if(gamma) print(p.gamma)
        if(kappa) print(p.kappa)
        if(chi) print(p.chi)
    if(!is.null(pdf.file)) dev.off()

}




## plot_hbm_dens -------------------------------------------
plot_hbm_dens <- function(stanfit,
                          gamma = FALSE,
                          kappa = FALSE,
                          chi = FALSE,
                          pdf.file = NULL,
                          width = 10,
                          height = 7) {
    ## HMB posterior density plots
    ##
    ## CI inner = 95%
    ##
    ## stanfit = stanfit model object
    ## pdf.filt = path to pdf file to save
    ## width = pdf file width
    ## height = pdf file height

    model <- deparse(substitute(stanfit))

    mu.alpha <- mcmc_areas(as.matrix(stanfit, pars = "mu_alpha"),
                           prob = 0.95, adjust = 1.5) +
    ggtitle(model)

    phi <- mcmc_areas(as.matrix(stanfit, pars = "phi"),
                      prob = 0.95, adjust = 1.5) +
    ggtitle(model)

    sig.alpha <- mcmc_areas(as.matrix(stanfit, pars = "sigma_alpha"),
                            prob = 0.95, adjust = 1.5) +
    ggtitle(model)


    if(gamma) {
        mu.gamma <- mcmc_areas(as.matrix(stanfit, pars = "mu_gamma"),
                               prob = 0.95, adjust = 1.5) +
        ggtitle(model)
        sig.gamma <- mcmc_areas(as.matrix(stanfit, pars = "sigma_gamma"),
                                prob = 0.95, adjust = 1.5) +
        ggtitle(model)
    }

    if(kappa) {
        mu.kappa <- mcmc_areas(as.matrix(stanfit, pars = "mu_kappa"),
                               prob = 0.95, adjust = 1.5) +
        ggtitle(model)
        sig.kappa <- mcmc_areas(as.matrix(stanfit, pars = "sigma_kappa"),
                                prob = 0.95, adjust = 1.5) +
        ggtitle(model)
    }

    if(chi) {
        mu.chi <- mcmc_areas(as.matrix(stanfit, pars = "mu_chi"),
                               prob = 0.95, adjust = 1.5) +
        ggtitle(model)
        sig.chi <- mcmc_areas(as.matrix(stanfit, pars = "sigma_chi"),
                                prob = 0.95, adjust = 1.5) +
        ggtitle(model)
    }



    if(!is.null(pdf.file)) pdf(pdf.file, width = width, height = height)
        print(mu.alpha)
        print(phi)
        print(sig.alpha)
        if(gamma) print(mu.gamma)
        if(gamma) print(sig.gamma)
        if(kappa) print(mu.kappa)
        if(kappa) print(sig.kappa)
        if(chi) print(mu.chi)
        if(chi) print(sig.chi)
    if(!is.null(pdf.file)) dev.off()

}



## plot_post_pc --------------------------------------------
plot_post_pc <- function(stanfit, y,
                         pdf.path = NULL,
                         var.yrep = "yrep",
                         data = sock) {
    ## Plot posterior predictive distributions
    ##
    ## stanfit = stanfit object sampled using priors only (no likelihood)
    ## y = response variable used to fit model
    ## pdf.path = file path to save graphics
    ## var.yrep = name of yrep parameter in stanfit

    yrep <- rstan::extract(stanfit, pars = var.yrep)

    if(!is.null(pdf.path)) {
        pdf(pdf.path, width = 10, height = 8)
    }
    
    g <- ppc_scatter_avg(y, yrep[[1]]) +
        labs(title="Observed vs. predicted")
    print(g)

    g <- ppc_stat_2d(y, yrep[[1]]) +
        labs(title="Y rep: target distributions")
    print(g)

    g <- ppc_dens_overlay(y = y, yrep = yrep[[1]][1:100, ]) +
        labs(title="Y rep: posterior predictive check")
    print(g)
    
    g <- ppc_dens_overlay_grouped(y = y, yrep = yrep[[1]][1:50, ], 
                                  group=data$Stock)  +
        labs(title="Y rep: posterior predictive check")
    print(g) 

    if(!is.null(pdf.path)) {
        dev.off()
    }
}



## hb_param_df ---------------------------------------------
hb_param_df <- function(stanfit, par, region.var, var = NULL, info = sock.info) {
    ## Wrangle hierarchical parameter summary into a data frame
    ## Output data.frame includes stock-specific and mean params
    ##
    ## stanfit = stanfit object,
    ## par = parameter name
    ## region.var = column name in sock.info of region variable
    ## var = optional variable name, add as column "var" to output data.frame
  
    probs <- c(0.025, 0.05, 0.10, 0.50, 0.90, 0.95, 0.975)
    s.par <- as.data.frame(summary(stanfit, pars = par, probs = probs)[[1]])
    s.par$par <- par
    s.par$Stock <- info$Stock
    s.par$Ocean.Region2 <- info[[region.var]]
    s.par[[region.var]] <- info[[region.var]]

    s.mu <- as.data.frame(summary(stanfit, pars = paste0("mu_", par),
                                  probs = probs)[[1]])
    names(s.mu) <- paste0("mu_", names(s.mu))
    s.mu[[region.var]] <- unique(info[[region.var]])
    s.out <- plyr::join(s.par, s.mu, by = region.var)
    names(s.out)[names(s.out) == region.var] <- "region"
    s.out$var <- var
    return(s.out)
}


era_hb_param_df <- function(stanfit, par, mu = FALSE, region.var = "Ocean.Region2", neras = 3, info = sock.info){
  
  ## Parameter posteriors from Eras models, wrangle into dataframe
  
  if(mu){
    
    ## Regional summary dataframe (mu_gamma/ mu_kappa)
    probs <- c(0.025, 0.05, 0.10, 0.50, 0.90, 0.95, 0.975)
    reg_start <- info$Stock[match(unique(info[[region.var]]), info[[region.var]])]
    reg_end <- c(info$Stock[match(unique(info[[region.var]]), info[[region.var]])-1], 
                 info$Stock[nrow(info)])
    mu.parsxera <- paste0(rep(paste0("mu_", par), each=neras), 1:neras)
    summ <- rstan::summary(stanfit, pars = mu.parsxera, probs = probs)[[1]]
    df.era.reg.2c <- data.frame(Ocean.Region2 = rep(unique(info[[region.var]]), neras),
                                ystart = rep(reg_start, neras),
                                yend = rep(reg_end, neras),
                                reg_mean = summ[, "mean"], 
                                reg_se = summ[ ,"se_mean"],
                                lower_10 = summ[ , "10%"],
                                upper_90 = summ[ , "90%"], 
                                var = str_extract(rownames(summ), "\\D+"),
                                varnam = case_when(grepl("^mu_gamma", rownames(summ)) ~ "SST",
                                                   grepl("^mu_kappa", rownames(summ)) ~ "Competitors"),
                                era = case_when(str_extract(rownames(summ), "\\d") == "1" ~ "Early",
                                                str_extract(rownames(summ), "\\d") == "2" ~ "Middle",
                                                str_extract(rownames(summ), "\\d") == "3" ~ "Late",
                                                .ptype=factor( levels=c("Early", "Middle", "Late"))))
    return(df.era.reg.2c)
    
    
  } else {
    
  ## Stock lvl dataframe 
  probs <- c(0.025, 0.05, 0.10, 0.50, 0.90, 0.95, 0.975)
  parxera <- paste0(rep(par, each=neras), 1:neras)
  summ <- rstan::summary(stanfit, pars = parxera, probs = probs)[[1]]
  df.era.st.2c <- data.frame(Stock = rep(info$Stock, neras), 
                             Ocean.Region2 = rep(info[[region.var]], neras),
                             mu = summ[, "mean"],
                             se = summ[, "se_mean"],
                             lower_10 = summ[, "10%"],
                             upper_90 = summ[ , "90%"],
                             var = str_extract(rownames(summ), "\\D+"),
                             varnam = case_when(grepl("^gamma", rownames(summ)) ~ "SST",
                                                grepl("^kappa", rownames(summ)) ~ "Competitors"),
                             era = case_when(str_extract(rownames(summ), "\\d") == "1" ~ "Early",
                                             str_extract(rownames(summ), "\\d") == "2" ~ "Middle",
                                             str_extract(rownames(summ), "\\d") == "3" ~ "Late",
                                             .ptype=factor( levels=c("Early", "Middle", "Late"))))
  return(df.era.st.2c)

  }
}
  



## hb.posterior.list ---------------------------------------
hb.posterior.list <- function(stanfit) {
    ## Extract parameter posteriors as matrices
    ##
    ## This function outputs a list of matrices with one list element per model
    ## parameter. Each list element is a matrix with each column being a vector
    ## of mcmc samples. The columns represent stock or region level parameters.
    ##
    ## stanfit = stanfit object

    mcmc <- As.mcmc.list(stanfit)
    params.gr <- c("alpha", "sigma\\[", "gamma", "kappa",
                   "chi", "mu_alpha", "mu_gamma", "mu_kappa",
                   "mu_chi", "phi")
    params.nm <- gsub("\\[", "", params.gr, fixed = TRUE)
    lst <- vector("list", length(params.gr))
    names(lst) <- params.nm
    for(i in seq_along(params.gr)) {
        val <- codatools::coda_grep(mcmc, paste0("^", params.gr[i]),
                                    return.matrix = TRUE)
        if(!empty(val)) {
            lst[[i]] <- val
        }
    }
    lst.f <- Filter(Negate(is.null), lst) ## remove empty elements
    return(lst.f)
}

if(FALSE) {

    lst.f <- hb.posterior.list(hb07a)
    names(lst.f)
    lapply(lst.f, is.null)
    lapply(lst.f, dim)
    head(lst.f[[1]])
    head(lst.f[[7]])
    head(lst.f[["mu_gamma"]])

}



## stan utilities ------------------------------------------
rhat_highest <- function(stanfit, k = 4, pars) {
    rhat <- get_rhat(stanfit, pars = pars)
    rhat.max <- rev(sort(rhat)[(length(rhat) - k):length(rhat)])
    return(rhat.max)
}

neff_lowest <- function(stanfit, k = 4, pars) {
    neff <- get_neff(stanfit, pars = pars)
    neff.min <- sort(neff)[1:k]
    return(neff.min)
}

pairs_lowest <- function(stanfit, k = 4, pars) {
    n <- get_neff(stanfit, pars = pars)
    n.min <- names(sort(n))[1:k]
    pairs(stanfit, pars = n.min)
}

get_rhat <- function(stanfit, pars) {
    if(!is(stanfit, "stanfit")) {
        stop("Input not of class stanfit")
    }
    summary(stanfit, pars = pars)$summary[ , "Rhat"]
}

get_neff <- function(stanfit, pars) {
    if(!is(stanfit, "stanfit")) {
        stop("Input not of class stanfit")
    }
    summary(stanfit, pars = pars)$summary[ , "n_eff"]
}

total_draws <- function(stanfit) {
    ## N chains * N draws  -- post warmup
    if(!is(stanfit, "stanfit")) {
        stop("Input not of class stanfit")
    }
    dim(stanfit)[1] * dim(stanfit)[2]
}



## stan_data_stat -----------------------------------------------
stan_data_stat <- function(data,
                      var.x2 = "early_sst_stnd",
                      var.x3 = "np_pinks_stnd",
                      var.region = "Ocean.Region2",
                      scale.x1 = FALSE,
                      priors.only = FALSE) {
    ## Get list of data for input to Stan for stationary models 
    ##
    ## data = data.frame of salmon data
    ## var.x2 = column name in `data` of x2 variable
    ## var.x3 = column name in `data` of x3 variable
    ## scale.x1 = logical, if TRUE spawners variable is scaled
    ## priors.only = logical indicating if a likelihood should be calculated
    ##               TRUE indicates prior predictive distributions will be
    ##               sampled only
  

    sock.stan <- data

    ## Set factor levels for Ocean.Region
    sock.stan[[var.region]] <- factor(sock.stan[[var.region]],
                                     levels = unique(sock.stan[[var.region]]))

#browser()
    ## Get start/end for each stock
    start.end  <- levels_start_end(sock.stan$Stock)

    ## Get grouping factor for series
    sock.stan[["OC_REGION_DUMMY"]] <- sock.stan[[var.region]]
    #browser()
    grp.df <- plyr::ddply(sock.stan, .(Stock), plyr::summarize,
                          group = unique(OC_REGION_DUMMY))
    g.grp <- as.numeric(factor(grp.df$group, levels = unique(grp.df$group)))

    a.grp.df <- plyr::ddply(sock.stan, .(Stock), plyr::summarize,
                            group = unique(Region))
    a.grp <- ifelse(a.grp.df$group == "Fraser River", 1, 2)


    ## Get start/end for group-specific gamma series
    start.end.grp.lst <- lapply(split(sock.stan, sock.stan[[var.region]]),
                                function(x) min(x$BY):max(x$BY))
    start.end.grp.vec <- unlist(lapply(seq_along(start.end.grp.lst),
                                       function(x) rep(x, length(start.end.grp.lst[[x]]))))
    start.end.grp <- levels_start_end(as.factor(start.end.grp.vec))


    ## Get year indices (map gamma -> y)
    year.lst  <- lapply(split(sock.stan, sock.stan[[var.region]]),
                        function(x) as.numeric(as.factor(x$BY)))
    if(length(unique(sock.stan[[var.region]])) == 1) {
        year <- year.lst[[1]]
    } else {
        year <- c(year.lst$WC, year.lst$GOA + max(year.lst$WC))
        year <- c(year, year.lst$BS + max(year))
    }

    if(scale.x1) {
        x1 <- plyr::ddply(data, .(Stock), transform,
                          x1 = scale(S)[ , 1])$x1
    } else {
        x1 = data$S / 1e5
    }

    lst <- list(y = sock.stan$lnRS,
                x1 = x1,
                x2 = sock.stan[[var.x2]],
                x3 = sock.stan[[var.x3]],
                g_group = g.grp,
                a_group = a.grp,
                year = year,
                n_series = length(unique(sock.stan$Stock)),
                Ng_groups = length(unique(sock.stan[[var.region]])),
                Na_groups = length(unique(a.grp)),
                y_start = start.end$start,
                y_end = start.end$end,
                g_start = start.end.grp$start,
                g_end = start.end.grp$end,
                Ng = length(start.end.grp.vec),
                priors_only = ifelse(priors.only, 1, 0),
                N = nrow(sock.stan))
    out <- Filter(Negate(is.null), lst)
    return(out)
}

if(FALSE) {

    stan_data(sock)

}

## -- stan_data_dyn --------------------------------------------
stan_data_dyn <- function(data,
                      var.x2 = "early_sst_stnd",
                      var.x3 = "np_pinks_sec_stnd", 
                      breakpoint1 = 1989,
                      breakpoint2 = NULL,
                      scale.x1 = FALSE,
                      var.region = "Ocean.Region2",
                      prior.sigma.gamma = c(df = 3, mu = 0, sd = 0.5),
                      prior.sigma.kappa = c(df = 3, mu = 0, sd = 0.5),
                      priors.only = FALSE) {
  ## Get list of data for input to Stan 'dynamic' models (i.e. random walk and era)
  ##
  ## data = data.frame of salmon data
  ## var.x2 = column name in `data` of x2 variable
  ## var.x3 = column name in 'data' of x3 variable, if any
  ## breakpoint1 = year for start of second era
  ## breakpoint2 = year for start of third era, if NULL, only two eras are
  ##               returned defined by breapoint1
  ## scale.x1 = logical, should the x1 (scaled) be scaled to N(0, 1)
  ## prior.sigma.gamma = Student-t prior parameters for sigma_gamma
  ## priors.only = logical indicating if a likelihood should be calculated
  ##               TRUE indicates prior predictive distributions will be
  ##               sampled only
  
  sock.stan <- data
  ## Set factor levels for var.region
  sock.stan[[var.region]] <- factor(sock.stan[[var.region]],
                                   levels = unique(sock.stan[[var.region]]))
  
  
  ## Get start/end for each Stock
  start.end  <- levels_start_end(sock.stan$Stock)
  
  ## Get grouping factor for series
  sock.stan[["OC_REGION_DUMMY"]] <- sock.stan[[var.region]]
  
  grp.df <- plyr::ddply(sock.stan, .(Stock), plyr::summarize,
                        group = unique(OC_REGION_DUMMY))
  g.grp <- as.numeric(factor(grp.df$group, levels = unique(grp.df$group)))
  
  a.grp.df <- plyr::ddply(sock.stan, .(Stock), plyr::summarize,
                          group = unique(Region))
  a.grp <- ifelse(a.grp.df$group == "Fraser River", 1, 2)
  
  
  ## Get start/end for group-specific gamma series
  start.end.grp.lst <- lapply(split(sock.stan, sock.stan[[var.region]]),
                              function(x) min(x$BY):max(x$BY))
  start.end.grp.vec <- unlist(lapply(seq_along(start.end.grp.lst),
                                     function(x) rep(x, length(start.end.grp.lst[[x]]))))
  start.end.grp <- levels_start_end(as.factor(start.end.grp.vec))
  
  
  ## Get year indices (map gamma -> y) ## Caution: Regions hard coded here?
   year.lst  <- lapply(split(sock.stan, sock.stan[[var.region]]),
                      function(x) as.numeric(as.factor(x$BY)))
  if(length(unique(sock.stan[[var.region]])) == 1) {
    year <- year.lst[[1]]
  } else {
    year <- c(year.lst$WC, year.lst$GOA + max(year.lst$WC))
    year <- c(year, year.lst$BS + max(year))
  }
  
  if(scale.x1) {
    x1 <- plyr::ddply(data, .(Stock), transform,
                      x1 = scale(S)[ , 1])$x1
  } else {
    x1 = data$S
  }
  
  if(is.null(breakpoint2)) {
    era1 <- ifelse(sock.stan$BY < breakpoint1, 1, 0)
    era2 <- ifelse(sock.stan$BY >= breakpoint1, 1, 0)
    era3 <- NULL
  } else {
    era1 <- ifelse(sock.stan$BY < breakpoint1, 1, 0)
    era2 <- ifelse(sock.stan$BY >= breakpoint1 &
                     sock.stan$BY < breakpoint2, 1, 0)
    era3 <- ifelse(sock.stan$BY >= breakpoint2, 1, 0)
  }
  
  
  lst <- list(y = sock.stan$lnRS,
              x1 = x1,
              x2 = sock.stan[[var.x2]],
              x3 = sock.stan[[var.x3]],
              g_group = g.grp,
              a_group = a.grp,
              year = year,
              era1 = era1,
              era2 = era2,
              era3 = era3,
              n_series = length(unique(sock.stan$Stock)),
              Ng_groups = length(unique(sock.stan[[var.region]])),
              Na_groups = length(unique(a.grp)),
              y_start = start.end$start,
              y_end = start.end$end,
              g_start = start.end.grp$start,
              g_end = start.end.grp$end,
              Ng = length(start.end.grp.vec),
              priors_only = ifelse(priors.only, 1, 0),
              sigma_gamma_df = prior.sigma.gamma[1],
              sigma_gamma_mu = prior.sigma.gamma[2],
              sigma_gamma_sd = prior.sigma.gamma[3],
              sigma_kappa_df = prior.sigma.kappa[1],
              sigma_kappa_mu = prior.sigma.kappa[2],
              sigma_kappa_sd = prior.sigma.kappa[3],
              N = nrow(sock.stan))
  out <- Filter(Negate(is.null), lst)
  return(out)
}

if(FALSE) {
  stan_data(sock.covar, "early_sst_stnd")
}




## level_start_end -----------------------------------------
levels_start_end <- function(factor) {
    ## Find start and end point of levelss in a factor
    ##
    ## factor = factor (vector)

    n <- as.numeric(factor)
    internal.end <- which(diff(n) != 0)
    internal.start <- internal.end + 1

    end <- c(internal.end, length(n))
    start <- c(1, internal.start)

    list(start = start, end = end)
}

if(FALSE) {
    f <- as.factor(c("a", "a", "a", "b", "b", "c", "c", "c"))
    levels_start_end(f)

    f <- as.factor(c("a", "a", "a", "b", "c", "c", "c"))
    levels_start_end(f)

    levels_start_end(sock$Stock)
}




## consec_years --------------------------------------------
consec_years <- function(x) {
    ## Find longest consective years
    ##
    ## This function takes as input a vector of years with NA values and finds
    ## the longest set of consecutive years without an NA. The function returns
    ## the vector of consecutive years.
    ##
    ## x = vector of years

    n <- ifelse(is.na(x), 0, x)
    brk <- c(0, which(diff(n) != 1), length(n))
    vec <- lapply(seq(length(brk) - 1),
                  function(i) n[(brk[i] + 1):brk[i+1]])
    consec <- vec[[which.max(lapply(vec, length))]]
    ret <- n[n %in% consec]

    return(ret)
}

if(FALSE) {

    n <- c(NA, 1950:1955)
    consec_years(n)

    n <- c(1949, NA, 1951:1955)
    consec_years(n)

    n <- c(NA, 1951:1955)
    consec_years(n)

    n <- c(1951:1955, NA)
    consec_years(n)

    n <- c(NA, 1950:1955, NA, 1957:1965, NA, 1967)
    consec_years(n)

    n <- c(1950:1955, 1957:1965, NA, 1967:1970)
    consec_years(n)

}



## load_rdata ----------------------------------------------
load_rdata <- function(path = "./output/", verbose = TRUE) {
    ## Load saved .RData output
    ##
    ## This function loads all *.RData files found in the input directory.
    ##
    ## path = path to directory to read output
    ## verbose = logical, print status

    if(!dir.exists(path))
        stop("Input path doesn't exist")

    fls <- list.files(path = path)
    fls <- fls[fls != "models"]

    for(i in fls) {

        if(verbose)
            cat("Reading file", which(fls == i), "of",
                length(fls), "...", "\n")

        load(paste(path, i, sep = ""), envir = .GlobalEnv)
    }

    if(verbose) {
        if(length(fls) == 0)
            cat("No files to read", "\n")
        else
            cat("Done!", "\n")
    }
}



## clim.wgt.avg --------------------------------------------
clim.wgt.avg <- function(brood.table,
                         env.data,
                         env.covar,
                         type,
                         out.covar = "index") {
    ## Calculate climate index weighted by age-class
    ##
    ## Create data frame with estimates of climate weighted by the proportion of
    ## 	recruits for each age-class. Indices can be created for either the first
    ## 	ocean year or second ocean year.
    ##
    ## brood.table = sockeye salmon brood table (need columns $BY $Stock.ID $Rx.x)
    ## env.data = data.frame of pink data (need column ($Year)
    ## env.covar = column name in `pink.data` for calc covar with
    ## 	type = c("first_year", "second_year")
    ## out.covar = name of column for new covar in output data.frame


    if("year" %in% names(env.data))
        names(env.data)[names(env.data) == "year"] <- "Year"

    env.mat <- matrix(NA, length(unique(brood.table$BY)),
                      length(unique(brood.table$Stock.ID)),
                      dimnames=list(unique(brood.table$BY),unique(brood.table$Stock.ID)))
    env.mat <- env.mat[order(row.names(env.mat)), ]

    for (i in unique(brood.table$Stock.ID)){
        brood <- subset(brood.table, Stock.ID == i)
        if(type == "first_year")
            climate <- subset(env.data, Stock.ID == i)
        for (j in unique(brood$BY)){

            if(type == "first_year") {
                env.mat[as.character(j),as.character(i)] <-
                    brood$ocean_0[brood$BY == j] * climate[climate$Year == j+1, env.covar] +
                    brood$ocean_1[brood$BY == j] * climate[climate$Year == j+2, env.covar] +
                    brood$ocean_2[brood$BY == j] * climate[climate$Year == j+3, env.covar] 
                    #brood$ocean_3[brood$BY == j] * climate[climate$Year == j+4, env.covar] 
                    #brood$ocean_4[brood$BY == j] * climate[climate$Year == j+5, env.covar]
            }
            if(type == "second_year") {
                env.mat[as.character(j),as.character(i)] <-
                    brood$ocean_0[brood$BY == j] * env.data[env.data$Year == j+2, env.covar] +
                    brood$ocean_1[brood$BY == j] * env.data[env.data$Year == j+3, env.covar] +
                    brood$ocean_2[brood$BY == j] * env.data[env.data$Year == j+4, env.covar] +
                    brood$ocean_3[brood$BY == j] * env.data[env.data$Year == j+5, env.covar] +
                    brood$ocean_4[brood$BY == j] * env.data[env.data$Year == j+6, env.covar]
            }
        } # end j loop
    } # end i loop
    long.df <- reshape2::melt(env.mat,
                              measure.vars=c((min(brood.table$BY):max(brood.table$BY)),
                                             unique(brood.table$Stock.ID)))
    colnames(long.df) <- c("BY","Stock.ID",out.covar)
    long.df <- long.df[complete.cases(long.df), ]
    return(long.df)
}

if(FALSE) {
    bt.raw <- read.csv("data/master_brood_table.csv", header=T)
    bt.complete <- bt.raw[complete.cases(bt.raw),]
    raw.clim <- read.csv(file="data/sst_yr_1_stock_anomalies.csv",header=TRUE)
    raw.clim.late <- read.csv(file="data/second_year_ocean_indices.csv",
                              header=TRUE)
    raw.clim.mld.late <- read.csv(file="data-downloaded/rob-second-year-indices-2018-02-16.csv",
                                  header=TRUE)

    sst1 <- clim.wgt.avg(bt.complete, raw.clim, "sst_anomaly",
                         "first_year", "early_sst")
    sst2 <- clim.wgt.avg(bt.complete, raw.clim.late, "sst.winter.djfm",
                         "second_year", "late_sst")


    mld.wi2.wc <- clim.wgt.avg(bt.complete[bt.complete$Ocean.Region2 != "BS", ],
                               raw.clim.mld.late, "GoAmld_winter",
                               "second_year", "late_mld_winter")
    mld.wi2.bs <- clim.wgt.avg(bt.complete[bt.complete$Ocean.Region2 == "BS", ],
                               raw.clim.mld.late, "Amld_winter",
                               "second_year", "late_mld_winter")
    mld.wi2 <- rbind(mld.wi2.wc, mld.wi2.bs)

    mld.su2.wc <- clim.wgt.avg(bt.complete[bt.complete$Ocean.Region2 != "BS", ],
                               raw.clim.mld.late, "GoAmld_summer",
                               "second_year", "late_mld_summer")
    mld.su2.bs <- clim.wgt.avg(bt.complete[bt.complete$Ocean.Region2 == "BS", ],
                               raw.clim.mld.late, "Amld_summer",
                               "second_year", "late_mld_summer")
    mld.su2 <- rbind(mld.su2.wc, mld.su2.bs)

    # all.equal(sst1, na.omit(long.early.sst), check.attributes = FALSE)
    # all.equal(sst2, na.omit(long.late.sst), check.attributes = FALSE)
    # all.equal(mld.wi2, na.omit(long.late.mld.winter), check.attributes = FALSE)
    # all.equal(mld.su2, na.omit(long.late.mld.summer), check.attributes = FALSE)

}



## pink.wgt.avg --------------------------------------------
pink.wgt.avg <- function(brood.table,
                         pink.data,
                         pink.covar,
                         type,
                         include_4.X = FALSE,
                         out.covar = "pink_index") {
    ## Calculate pink abundance index weighted by age-class
    ##
    ## Create data frame with estimates of competitors weighted by
    ## the proportion of recruits for each age-class. Competitor index is either
    ## (1) just second year competitor abundance, (2) a geographically varying
    ## second year competitor index (by+4 for WC and GOA stocks and by+3 for BS
    ## stocks)
    ##
    ## brood.table = sockeye salmon brood table (need columns $BY $Stock.ID $Rx.x)
    ## pink.data = data.frame of pink data (need column ($Year)
    ## pink.covar = column name in `pink.data` to calc covar with, can be a
    ##              single string or vector of length three specifying different
    ##              columns for the Ocean.Regions. In the latter case the names
    ##              of the vector need to be the Ocean.Regions
    ## type = c("second_year", "geographic")
    ## out.covar = name of column for new covar in output data.frame

    np.pink <- matrix(NA,length(unique(brood.table$BY)),
                          length(unique(brood.table$Stock.ID)),
                          dimnames=list(unique(brood.table$BY), unique(brood.table$Stock.ID)))
    np.pink <- np.pink[order(row.names(np.pink)), ]

    if(length(pink.covar) == 1) {
        pink.covar <- rep(pink.covar, 4)
        names(pink.covar) <- unique(brood.table$Ocean.Region2)
    }

    for (i in unique(brood.table$Stock.ID)){
        brood <- subset(brood.table, Stock.ID == i)
        for (j in unique(brood$BY)){

            if(type == "second_year") {
                reg_i <- as.vector(unique(brood.table$Ocean.Region2[brood.table$Stock.ID == i]))

                  if(brood$DetailFlag[brood$BY == j] ==1){ # HH added detail condition to accommodate new stocks that don't have detailed recruit information 
                    np.pink[as.character(j),as.character(i)] <-
                      (brood$R0.1[brood$BY == j] * 0) +
                      (brood$R0.2[brood$BY == j] * pink.data[pink.data$Year == j+3, pink.covar[reg_i]]) +
                      (brood$R0.3[brood$BY == j] * pink.data[pink.data$Year == j+3, pink.covar[reg_i]]) +
                      (brood$R0.4[brood$BY == j] * pink.data[pink.data$Year == j+3, pink.covar[reg_i]]) +
                      (brood$R0.5[brood$BY == j] * pink.data[pink.data$Year == j+3, pink.covar[reg_i]]) +
                      brood$R1.1[brood$BY == j] * 0 +
                      (brood$R1.2[brood$BY == j] * pink.data[pink.data$Year == j+4, pink.covar[reg_i]]) +
                      (brood$R1.3[brood$BY == j] * pink.data[pink.data$Year == j+4, pink.covar[reg_i]]) +
                      (brood$R1.4[brood$BY == j] * pink.data[pink.data$Year == j+4, pink.covar[reg_i]]) +
                      (brood$R1.5[brood$BY == j] * pink.data[pink.data$Year == j+4, pink.covar[reg_i]]) +
                      brood$R2.1[brood$BY == j] * 0 +
                      (brood$R2.2[brood$BY == j] * pink.data[pink.data$Year == j+5, pink.covar[reg_i]])+
                      (brood$R2.3[brood$BY == j] * pink.data[pink.data$Year == j+5, pink.covar[reg_i]])+
                      (brood$R2.4[brood$BY == j] * pink.data[pink.data$Year == j+5, pink.covar[reg_i]])+
                      (brood$R2.5[brood$BY == j] * pink.data[pink.data$Year == j+5, pink.covar[reg_i]])+
                      brood$R3.1[brood$BY == j] * 0 +
                      (brood$R3.2[brood$BY == j] * pink.data[pink.data$Year == j+6, pink.covar[reg_i]])+
                      (brood$R3.3[brood$BY == j] * pink.data[pink.data$Year == j+6, pink.covar[reg_i]])+
                      (brood$R3.4[brood$BY == j] * pink.data[pink.data$Year == j+6, pink.covar[reg_i]])
                    } else {
                    np.pink[as.character(j),as.character(i)] <-
                      (brood$ocean_0[brood$BY == j] * pink.data[pink.data$Year == j+3, pink.covar[reg_i]]) +
                      (brood$ocean_1[brood$BY == j] * pink.data[pink.data$Year == j+4, pink.covar[reg_i]]) +
                      (brood$ocean_2[brood$BY == j] * pink.data[pink.data$Year == j+5, pink.covar[reg_i]]) +  
                      (brood$ocean_3[brood$BY == j] * pink.data[pink.data$Year == j+6, pink.covar[reg_i]])  
                    } # end detail condition
            }

            if(type == "geographic") {
                reg_i <- as.vector(unique(brood.table$Ocean.Region2[brood.table$Stock.ID == i]))
                isBS <- as.vector(unique(brood.table$Ocean.Region2[brood.table$Stock.ID == i])) == "BS"
                if(!isBS){
                    np.pink[as.character(j),as.character(i)] <-
                        (brood$R0.1[brood$BY == j] * 0) +
                        (brood$R0.2[brood$BY == j] * pink.data[pink.data$Year == j+3, pink.covar[reg_i]]) +
                        (brood$R0.3[brood$BY == j] * pink.data[pink.data$Year == j+3, pink.covar[reg_i]]) +
                        (brood$R0.4[brood$BY == j] * pink.data[pink.data$Year == j+3, pink.covar[reg_i]]) +
                        (brood$R0.5[brood$BY == j] * pink.data[pink.data$Year == j+3, pink.covar[reg_i]]) +
                        brood$R1.1[brood$BY == j] * 0 +
                        (brood$R1.2[brood$BY == j] * pink.data[pink.data$Year == j+4, pink.covar[reg_i]]) +
                        (brood$R1.3[brood$BY == j] * pink.data[pink.data$Year == j+4, pink.covar[reg_i]]) +
                        (brood$R1.4[brood$BY == j] * pink.data[pink.data$Year == j+4, pink.covar[reg_i]]) +
                        (brood$R1.5[brood$BY == j] * pink.data[pink.data$Year == j+4, pink.covar[reg_i]]) +
                        brood$R2.1[brood$BY == j] * 0 +
                        (brood$R2.2[brood$BY == j] * pink.data[pink.data$Year == j+5, pink.covar[reg_i]])+
                        (brood$R2.3[brood$BY == j] * pink.data[pink.data$Year == j+5, pink.covar[reg_i]])+
                        (brood$R2.4[brood$BY == j] * pink.data[pink.data$Year == j+5, pink.covar[reg_i]])+
                        (brood$R2.5[brood$BY == j] * pink.data[pink.data$Year == j+5, pink.covar[reg_i]])+
                        brood$R3.1[brood$BY == j] * 0 +
                        (brood$R3.2[brood$BY == j] * pink.data[pink.data$Year == j+6, pink.covar[reg_i]])+
                        (brood$R3.3[brood$BY == j] * pink.data[pink.data$Year == j+6, pink.covar[reg_i]])+
                        (brood$R3.4[brood$BY == j] * pink.data[pink.data$Year == j+6, pink.covar[reg_i]])
                } else {
                    np.pink[as.character(j),as.character(i)] <-
                        (brood$R0.1[brood$BY == j] * 0) +
                        (brood$R0.2[brood$BY == j] * pink.data[pink.data$Year == j+2, pink.covar[reg_i]]) +
                        (brood$R0.3[brood$BY == j] * pink.data[pink.data$Year == j+2, pink.covar[reg_i]]) +
                        (brood$R0.4[brood$BY == j] * pink.data[pink.data$Year == j+2, pink.covar[reg_i]]) +
                        (brood$R0.5[brood$BY == j] * pink.data[pink.data$Year == j+2, pink.covar[reg_i]]) +
                        brood$R1.1[brood$BY == j] * 0 +
                        (brood$R1.2[brood$BY == j] * pink.data[pink.data$Year == j+3, pink.covar[reg_i]]) +
                        (brood$R1.3[brood$BY == j] * pink.data[pink.data$Year == j+3, pink.covar[reg_i]]) +
                        (brood$R1.4[brood$BY == j] * pink.data[pink.data$Year == j+3, pink.covar[reg_i]]) +
                        (brood$R1.5[brood$BY == j] * pink.data[pink.data$Year == j+3, pink.covar[reg_i]]) +
                        brood$R2.1[brood$BY == j] * 0 +
                        (brood$R2.2[brood$BY == j] * pink.data[pink.data$Year == j+4, pink.covar[reg_i]])+
                        (brood$R2.3[brood$BY == j] * pink.data[pink.data$Year == j+4, pink.covar[reg_i]])+
                        (brood$R2.4[brood$BY == j] * pink.data[pink.data$Year == j+4, pink.covar[reg_i]])+
                        (brood$R2.5[brood$BY == j] * pink.data[pink.data$Year == j+4, pink.covar[reg_i]])+
                        brood$R3.1[brood$BY == j] * 0 +
                        (brood$R3.2[brood$BY == j] * pink.data[pink.data$Year == j+5, pink.covar[reg_i]])+
                        (brood$R3.3[brood$BY == j] * pink.data[pink.data$Year == j+5, pink.covar[reg_i]])+
                        (brood$R3.4[brood$BY == j] * pink.data[pink.data$Year == j+5, pink.covar[reg_i]])
                }
            }

        } # end j loop
    } # end i loop

    long.df <- reshape2::melt(np.pink,
                              measure.vars=c((min(brood.table $BY):max(brood.table $BY)),
                                             unique(brood.table$Stock.ID)))
    colnames(long.df) <- c("BY","Stock.ID",out.covar)
    return(long.df)
}

if(FALSE) {
    bt.raw <- read.csv("data/master_brood_table.csv", header=T)
    bt.complete <- bt.raw[complete.cases(bt.raw),]
    raw.comp <- read.csv(file="data-downloaded/pink_abundance_2017_12_08.csv",
                         header=TRUE)

    np <- pink.wgt.avg(bt.complete, raw.comp, "Total",
                       "second_year", "np_pinks")
    go <- pink.wgt.avg(bt.complete, raw.comp, "Total",
                       "geographic", "geo_np_pinks")

    # all.equal(cu, long.cum.np.pink)
    # all.equal(np, long.np.pink)
    # all.equal(go, long.geo.np.pink)
}



## add.label -----------------------------------------------
add.label <- function(label, xfrac = 0.01, yfrac = 0.07, pos = 4, ...) {
    u <- par("usr")
    x <- u[1] + xfrac * (u[2] - u[1])
    y <- u[4] - yfrac * (u[4] - u[3])
    text(x, y, label, pos = pos, ...)
}


## single.stock.fit ----------------------------------------
single.stock.fit <- function(formulas, years, plot.path) {
    ## Fit single-stock Ricker models
    ##
    ## This function takes as input a list of linear model formulas and uses the
    ## nlme::lmList function to fit each formula in the list to all salmon
    ## stocks. The function outputs a list with two elements:
    ##  $coef = list of data.frames of model coeffients with one df for each formula
    ##  $rsq  = list of data.frames of r-squared values with one df for each formula
    ##
    ## The function create several useful graphics of the model fits that are
    ## saved in the directory `plot.path`.
    ##
    ## formulas = a list of linear model formulas to input into lmList()
    ## years = vectors of years to use
    ## plot.path = path to directory to save graphics

    if(!dir.exists(plot.path)) {
        dir.create(plot.path, recursive = TRUE)
    }

    mod.list <- formulas
    mod.coef <- vector("list", length(mod.list))
    mod.coef.oc2 <- vector("list", length(mod.list))
    mod.rsq  <- vector("list", length(mod.list))
    names(mod.coef) <- names(mod.list)
    names(mod.coef.oc2) <- names(mod.list)

    ## get unique ocean.regions and stock names
    reg1 <- subset(sock, select = c("Stock", "Ocean.Region"))
    regions <- reg1[!duplicated(reg1), ]
    reg2 <- subset(sock, select = c("Stock", "Ocean.Region2"))
    regions2 <- reg2[!duplicated(reg2), ]
    

    for(i in seq_along(mod.list)) {
        m.formula <- mod.list[[i]]
        m.name <- names(mod.list)[i]
        g.sub <- paste0(m.name, ": ", gsub("  ", "", paste(deparse(m.formula), collapse = "")))
        i.path <- paste0(plot.path, m.name, "/")

        if(!dir.exists(i.path)) {
            dir.create(i.path, recursive = TRUE)
        }

        i.path.f <- paste0(i.path, m.name)

        ##### Debug
        # cat(m.name, "\n")
        # if(m.name == "model9a") browser()
        #####

        ## Fit lm models
        dat <- sock[sock$BY %in% years, ]
        m.fit <- nlme::lmList(m.formula, data = dat, na.action = na.omit)

        ## Create coefficient data.frame
        m.coef <- coef(m.fit)
        m.coef$Stock <- names(m.fit)
        m.coef <- reshape2::melt(m.coef, id.vars = "Stock")
        m.coef$Stock <- factor(m.coef$Stock, levels = names(m.fit), ordered=is.ordered(sock$geo_id))
        m.coef <- merge(m.coef, regions, by = "Stock", all.x = TRUE, sort = FALSE)
        mod.coef[[i]] <- m.coef
        m.coef.oc2 <- merge(m.coef, regions2, by = "Stock", all.x = TRUE, sort = FALSE)
        mod.coef.oc2[[i]] <- m.coef.oc2

        ## Create r.squared data.frame
        m.rsq <- summary(m.fit)$r.squared
        m.rsq.df <- data.frame(model = m.name, stock = names(m.fit), r.squared = m.rsq)
        mod.rsq[[i]] <- m.rsq.df


        ## Create residual data.frame
        rf <- data.frame(Stock = names(fitted(m.fit)),
                         fitted = fitted(m.fit),
                         residuals = resid(m.fit, type = "pooled.pearson"))
        rf$Stock <- factor(rf$Stock, levels = names(m.fit))

        n.covar <- sum(!unique(as.character(m.coef$variable)) %in% c("(Intercept)", "S"))
        dot.width <- 10
        if(n.covar == 1)
            dot.width <- 6
        if(n.covar == 2)
            dot.width <- 8
        if(n.covar == 3)
            dot.width <- 10
        if(n.covar >= 4)
            dot.width <- 12
        if(n.covar >= 5)
            dot.width <- 14

        pdf(paste0(i.path.f, "_coef_hist.pdf"), width = 10, height = 8)
        g <- histogram( ~ value | variable, data = m.coef,
                       par.settings = theme.mjm(),
                       par.strip.text = list(cex = 0.7),
                       breaks = 7,
                       type = "density",
                       xlab = "Coefficient",
                       ylab = "Density",
                       scales = list(relation = "free"),
                       main = "Histograms of single-stock model coefficients",
                       sub = g.sub,
                       panel = function(x, ...) {
                           panel.histogram(x, ...)
                           panel.rug(x, col = "grey50")
                           panel.densityplot(x, col = "red2", pch = "",
                                             dargs = list(adjust = 1.3))
                       })
        print(g)
        dev.off()

        pdf(paste0(i.path.f, "_coef_hist_region.pdf"), width = 10, height = 8)
        g <- densityplot( ~ value | variable, data = m.coef,
                         groups = Ocean.Region,
                         adjust = 1.3,
                         scales = list(relation = "free"),
                         par.settings = theme.mjm(lwd = 2),
                         par.strip.text = list(cex = 0.7),
                         auto.key = list(space = "right"),
                         xlab = "Coefficient",
                         ylab = "Density",
                         main = "Density of single-stock model coefficients by region",
                         sub = g.sub,
                         panel = function(x, groups = groups, ...) {
                             panel.abline(v = 0, col = "grey60", lty = 2)
                             panel.densityplot(x, groups = groups, ...)
                         })
        print(g)
        dev.off()
        
        # colours = first ocean region grouping (3 groups)
        pdf(paste0(i.path.f, "_coef_dot_all.pdf"), width = dot.width + 2, height = 8)
        g <- xyplot(as.factor(Stock) ~ value | variable, data = m.coef,
                    groups = Ocean.Region,
                    par.settings = theme.mjm(),
                    par.strip.text = list(cex = 0.7),
                    layout = c(NA, 1),
                    scales = list(x = list(relation = "free")),
                    ylab = "",
                    xlab = "Coefficient",
                    main = "Single-stock model coefficients by stock",
                    sub = g.sub,
                    auto.key = list(space = "right"),
                    panel = function(x, y, ...) {
                        panel.abline(v = 0, col = "grey60", lty = 2)
                        panel.xyplot(x, y, ...)
                    })
        print(g)
        dev.off()
        
        # colours = first ocean region grouping (3 groups)
        if(m.name != "model1a") {
          pdf(paste0(i.path.f, "_coef_dot_covars.pdf"), width = dot.width, height = 8)
          
          g <- xyplot(as.factor(Stock) ~ value | variable, data = m.coef,
                      subset = !variable %in% c("(Intercept)", "S"),
                      groups = Ocean.Region,
                      par.settings = theme.mjm(),
                      par.strip.text = list(cex = 0.7),
                      layout = c(NA, 1),
                      scales = list(x = list(relation = "same")),
                      ylab = "",
                      xlab = "Coefficient",
                      main = "Single-stock model coefficients by stock",
                      sub = g.sub,
                      auto.key = list(space = "right"),
                      panel = function(x, y, ...) {
                        panel.abline(v = 0, col = "grey60", lty = 2)
                        panel.xyplot(x, y, ...)
                      })
          print(g)
          
          dev.off()
        }

        # colours = second ocean region grouping (4 groups)
        pdf(paste0(i.path.f, "_coef_dot_all_oc2.pdf"), width = dot.width + 2, height = 8)
        g <- xyplot(as.factor(Stock) ~ value | variable, data = m.coef.oc2,
                    groups = Ocean.Region2,
                    par.settings = theme.mjm(),
                    par.strip.text = list(cex = 0.7),
                    layout = c(NA, 1),
                    scales = list(x = list(relation = "free")),
                    ylab = "",
                    xlab = "Coefficient",
                    main = "Single-stock model coefficients by stock",
                    sub = g.sub,
                    auto.key = list(space = "right"),
                    panel = function(x, y, ...) {
                      panel.abline(v = 0, col = "grey60", lty = 2)
                      panel.xyplot(x, y, ...)
                    })
        print(g)
        dev.off()

        # colours = second ocean region grouping (4 groups)
        if(m.name != "model1a") {
          pdf(paste0(i.path.f, "_coef_dot_covars_oc2.pdf"), width = dot.width, height = 8)
          
          g <- xyplot(as.factor(Stock) ~ value | variable, data = m.coef.oc2,
                      subset = !variable %in% c("(Intercept)", "S"),
                      groups = Ocean.Region2,
                      par.settings = theme.mjm(),
                      par.strip.text = list(cex = 0.7),
                      layout = c(NA, 1),
                      scales = list(x = list(relation = "same")),
                      ylab = "",
                      xlab = "Coefficient",
                      main = "Single-stock model coefficients by stock",
                      sub = g.sub,
                      auto.key = list(space = "right"),
                      panel = function(x, y, ...) {
                        panel.abline(v = 0, col = "grey60", lty = 2)
                        panel.xyplot(x, y, ...)
                      })
          print(g)
          
          dev.off()
        }

        
        pdf(paste0(i.path.f, "_resid_scatter.pdf"), width = 10, height = 8)
        g <- xyplot(residuals ~ fitted, data = rf,
                    par.settings = theme.mjm(),
                    par.strip.text = list(cex = 0.7),
                    ylab = "Standardized residuals",
                    xlab = "Fitted values",
                    main = "Residuals vs. fitted values across all stocks",
                    sub = g.sub,
                    panel = function(x, y, ...) {
                        panel.abline(h = 0, col = "grey60", lty = 1)
                        panel.rug(x, y, col = "grey50", lty = 1)
                        panel.xyplot(x, y, ...)
                    })
        print(g)
        dev.off()

        pdf(paste0(i.path.f, "_resid_scatter_stock.pdf"), width = 14, height = 9)
        g <- xyplot(residuals ~ fitted | Stock, data = rf,
                    par.settings = theme.mjm(),
                    par.strip.text = list(cex = 0.7),
                    scales = list(relation = "free"),
                    ylab = "Standardized residuals",
                    xlab = "Fitted values",
                    main = "Residuals vs. fitted values by stock",
                    sub = g.sub,
                    panel = function(x, y, ...) {
                        panel.abline(h = 0, col = "grey60", lty = 1)
                        panel.xyplot(x, y, ...)
                    })
        print(g)
        dev.off()

        pdf(paste0(i.path.f, "_resid_qq.pdf"), width = 10, height = 8)
        g <- qqmath( ~ residuals, data = rf,
                    par.settings = theme.mjm(),
                    par.strip.text = list(cex = 0.7),
                    distribution = qnorm,
                    ylab = "Empirical quantiles",
                    xlab = "Theoretical quantiles",
                    main = "QQ-plot of standardized residuals",
                    sub = g.sub,
                    panel = function(x, y, ...) {
                        panel.qqmathline(x, col = "grey50")
                        panel.qqmath(x, ...)
                    })
        print(g)
        dev.off()

        pdf(paste0(i.path.f, "_resid_qq_stock.pdf"), width = 14, height = 9)
        g <- qqmath( ~ residuals | Stock, data = rf,
                    par.settings = theme.mjm(),
                    par.strip.text = list(cex = 0.7),
                    distribution = qnorm,
                    ylab = "Empirical quantiles",
                    xlab = "Theoretical quantiles",
                    main = "QQ-plot of standardized residuals by stock",
                    sub = g.sub,
                    panel = function(x, y, ...) {
                        panel.qqmathline(x, col = "grey50")
                        panel.qqmath(x, ...)
                    })
        print(g)
        dev.off()

        ## plot(m.fit, resid(., type = "pooled.pearson") ~ fitted(.),
        ##      par.settings = theme.mjm(),
        ##      main = "Residuals vs. fitted values across all stock",
        ##      sub = g.sub,
        ##      grid = TRUE)
        ##
        ## qqnorm(m.fit, ~ resid(., type = "pooled.pearson"),
        ##        par.settings = theme.mjm(),
        ##        main = "QQ-plot of standardized residuals",
        ##        sub = g.sub,
        ##        abline = c(0, 1))

        ## Variance inflation factor
        if(m.name != "model1a") {
            lst <- lapply(m.fit, function(x) {
                              v <- suppressWarnings(car::vif(x))
                              df <- as.data.frame(matrix(v, ncol = length(v)))
                              names(df) <- names(v)
                              return(df)
                    })
            vf <- do.call("rbind", c(lst, make.row.names = FALSE))
            vf$Stock <- names(m.fit)
            vf$Stock <- factor(vf$Stock, levels = names(m.fit))
            vf <- reshape2::melt(vf, id.vars = "Stock")
            vf$value[is.nan(vf$value)] <- NA

            # if(any(is.na(vf$value))) browser()
            if(max(vf$value, na.rm = TRUE) <= 4.5) {
                x.max <- 4.5
            } else {
                x.max <- NA
            }

            pdf(paste0(i.path.f, "_vif_dot.pdf"), width = 12, height = 8)
            g <- xyplot(as.factor(Stock) ~ value | variable, data = vf,
                        par.settings = theme.mjm(),
                        par.strip.text = list(cex = 0.7),
                        layout = c(NA, 1),
                        xlim = c(NA, x.max),
                        ylab = "",
                        xlab = "VIF",
                        main = "Variance inflation factors",
                        sub = g.sub,
                        panel = function(x, y, ...) {
                            panel.abline(v = 4, col = "red2", lty = 2)
                            panel.xyplot(x, y, ...)
                        })
            print(g)
            dev.off()
        }
    }
    out <- list(coef = mod.coef.oc2, rsq = mod.rsq)
    return(out)
}



## single.stock.compare ------------------------------------
single.stock.compare <- function(data,
                                 model1,
                                 model2,
                                 compare,
                                 levels = NULL,
                                 main = "") {
    ## Compare single stock model coefficients
    ##
    ## data = list of single-stock model coefficients
    ## model1 = name of first model to compare
    ## model2 = name of second model to compare
    ## compare = names of model coefficients to compare
    ## levels = alternative names for 'compare' variables (used for legend)
    ##          needs to be in same order as 'compare'

    if(!is.null(levels) & length(compare) != length(levels))
        stop("Length of compare and levels is not equal", call. = FALSE)

    if(!model1 %in% names(data))
        stop("model1 one not in data", call. = FALSE)

    if(!model2 %in% names(data))
        stop("model2 one not in data", call. = FALSE)

    m1 <- model1
    m2 <- model2

    d1 <- data[[m1]]
    d2 <- data[[m2]]
    d1$model <- m1
    d2$model <- m2
    df <- rbind(d1, d2)
    df <- df[df$variable %in% compare, ]
    df$variable <- as.character(df$variable)

    ## If levels are present replace the var names in data
    if(!is.null(levels)) {
        for(i in seq_along(compare)) {
            df$variable[df$variable == compare[i]] <- levels[i]
        }
    }

    g <- xyplot(as.factor(Stock) ~ value, data = df,
                groups = variable,
                par.settings = theme.mjm(),
                layout = c(NA, 1),
                scales = list(x = list(relation = "free")),
                ylab = "",
                xlab = "Coefficient",
                main = paste("Compare", m1, "and", m2, main),
                auto.key = list(space = "top"),
                panel = function(x, y, ...) {
                    panel.abline(v = 0, col = "grey60", lty = 2)
                    panel.xyplot(x, y, ...)
                })
    print(g)

}



## fill.time.series ----------------------------------------
fill.time.series <- function(data) {

    ## Fill salmon data time series so that all brood years are consecutive,
    ## filling missing years with NA.
    ##
    ## This function takes as input, a brood table with columns `Stock.ID`,
    ## `Stock`, and `BY` and fills in any non-consecutive BY within a stocks
    ## times series with NA values for all other columns present in `data`.
    ##
    ## When filtering the data by the `use` column, data values in the middle of
    ## the time series for a particular salmon stocks also get removed if `use`
    ## is set to 0. This function adds back in those data points, setting them
    ## to NA.

    id <- unique(data$Stock.ID)
    lst <- vector("list", length(id))
    for(i in seq_along(id)) {
        sub <- data[data$Stock.ID == id[i], ]
        BY <- min(sub$BY):max(sub$BY)
        Stock.ID <- unique(sub$Stock.ID)
        Stock <- unique(sub$Stock)
        df <- data.frame(Stock.ID = Stock.ID, Stock = Stock, BY = BY,
                         stringsAsFactors = FALSE)
        lst[[i]] <- merge(df, sub, by = c("Stock.ID", "Stock", "BY"), all.x = TRUE)
    }
    df <- do.call("rbind", c(lst, make.row.names = FALSE))

    ## Don't want NA in these columns
    out <- plyr::ddply(df, .(Stock), transform,
                       Region = unique(na.omit(Region)),
                       Sub.Region = unique(na.omit(Sub.Region)),
                       Ocean.Region2 = unique(na.omit(Ocean.Region2)),
                       Lat = unique(na.omit(Lat)),
                       Lon = unique(na.omit(Lon)))
    return(out)
}



## get.npgo ------------------------------------------------
get.npgo <- function(years) {

    ## This function takes as input a range of years and downloads and processes
    ## the NPGO index. The output of the function is a dataframe in 'long'
    ## format with a column for year, month, and the NPGO index.
    ##
    ## years = vector of years

    if(min(years) < 1950)
        stop("Earliest NPGO year is 1950")

    npgo    <- read.table("http://www.o3d.org/npgo/npgo.php", sep = "\t",
                          strip.white = TRUE)
    n.npgo  <- length(npgo[ , 1])
    npgo    <- npgo[4:n.npgo, ]
    n.npgo  <- length(npgo)
    rm.tail <- n.npgo - 3
    npgo    <- npgo[1:rm.tail]
    npgo    <- as.character(npgo)
    npgo    <- strsplit(npgo, "  ")
    npgo    <- do.call("rbind", npgo)
    npgo    <- data.frame(year = as.numeric(npgo[ , 1]),
                          month = as.numeric(npgo[ , 2]),
                          npgo = as.numeric(npgo[ , 3]))
    npgo    <- npgo[npgo$year >= min(years) & npgo$year <= max(years), ]

    return(npgo)
}

if(FALSE) {

    get.npgo(1950:1950)
    get.npgo(1950:2013)

    npgo <- get.npgo(1950:2016)
    head(npgo)
    tail(npgo)
    sapply(npgo, class)
    summary(npgo)

}



## sst.averager --------------------------------------------
sst.averager <- function(info, sst, distance = 400) {

    ## This function takes as input a data.frame of sst data output from the
    ## sst.anomaly() function and computes sst averages for each stock only
    ## including sst grid cells that are within a specified distance from the
    ## ocean entry location of a particular stock.
    ##
    ## Different months are used to calculate the SST average based on where the
    ## salmon stock enters the ocean:
    ##  WA, BC, SEAK: April-July
    ##  GOA: May-August
    ##  BB, AYK: Jun-September
    ##
    ## The function outputs a data frame with columns:
    ##  $year = one year for each year in input SST data and stock.id
    ##  $stock.id = id number for a salmon stock
    ##  $sst = averaged raw SST values
    ##  $sst.anom = averaged SST anomalies
    ##
    ## Function arguments:
    ##   info = stock.info data.frame w/ stock.id number, lon, and lat, should
    ##          have one row per stock
    ##   sst = sst data output from sst.anomaly()
    ##   distance = distance in km from ocean entry location of stock a grid
    ##              cell can be to be included in the averaging. This distance
    ##              is measured to the center of the SST grid cell

    info$stock.id <- info$Stock.ID
    stock.id <- info$stock.id
    cells    <- unique(subset(sst, select = c(id, lat, lon2)))
    cells    <- cells[order(cells$id), ]
    n.cells  <- length(cells[ , 1])
    row.names(cells) <- NULL

    sst.out <- vector("list", length(stock.id))
    for(i in seq_along(stock.id)) {

        info.i     <- info[i , ]
        stock.id.i <- info.i$stock.id
        lat.i      <- info.i$lat
        lon.i      <- info.i$lon # Brendan removed "* -1" from this line of code; longitude was already converted to negative degrees east

        dist <- rep(NA, n.cells)
        for(j in 1:n.cells)
            dist[j] <- haversine(lat.i, lon.i, cells$lat[j], cells$lon2[j])

        cells.sub <- cells[which(dist <= distance), ]
        sst.sub   <- sst[sst$id %in% cells.sub$id, ]

        if(stock.id.i <= 137)
            months <- 4:7  ## WA, BC, SEAK

        if(stock.id.i > 137 & stock.id.i <= 152)
            months <- 5:8  ## GOA

        if(stock.id.i > 152)
            months <- 6:9  ## BB and AYK

        sst.sub.mnths <- sst.sub[sst.sub$month %in% months, ]

        sst.avg <- ddply(sst.sub.mnths, .(year), plyr::summarize,
                               sst = mean(sst, na.rm = TRUE),
                               sst.anom = mean(sst.anom, na.rm = TRUE))

        sst.avg$stock.id <- stock.id.i

        sst.out[[i]] <- sst.avg
    }
    sst.out <- rbind.fill(sst.out)
    return(sst.out)
}



## haversine -----------------------------------------------
haversine <- function(lat1, lon1, lat2, lon2) {
    ## This function computes the great circle distance between two points given
    ## their latitiude and longitude (in decimal degrees) using the haversine
    ## formula. The output is the distance between the two points in km.
    ##
    ## lat1 = latitude of first point
    ## lon1 = longitude of first point
    ## lat2 = latitude of second point
    ## lon2 = longitude of second point

    # Convert degrees to radians
    lat1 <- lat1 * pi / 180
    lon1 <- lon1 * pi / 180
    lat2 <- lat2 * pi / 180
    lon2 <- lon2 * pi / 180

    R <- 6371 # earth mean radius [km]
    delta.lon <- (lon2 - lon1)
    delta.lat <- (lat2 - lat1)
    a <- sin(delta.lat/2)^2 + cos(lat1) * cos(lat2) * sin(delta.lon/2)^2
    d <- 2 * R * asin(min(1, sqrt(a)))

    return(d) # distance in km
}

if(FALSE) {

    haversine(48, -147, 49, -154)

}



## enviro.avg.months ---------------------------------------
enviro.avg.months <- function(data,
                              first.month,
                              last.month,
                              avg.var,
                              month.var = "month",
                              year.var = "year",
                              grid.id.var = NULL,
                              lat.var = "lat",
                              lon.var = "lon",
                              type = avg.var) {
    ## Average an environmental variable over the specified months within a year
    ##
    ## Input data should be in a "long" format with a column for year, month,
    ## and the environmental variable to average. `first.month` and
    ## `last.month` give the first and last months of a continuous range to
    ## average the environmental variable over.
    ##
    ## If `first.month` is less than `last.month` the environmental variable is
    ## averaged over the months first.month:last.month within each year. For
    ## example, if `first.month` = 3 and `last.month` = 4, the environmental
    ## variable will be averaged over Mar and Apr for each year.
    ##
    ## If `first.month` equals `last.month`, that month is returned with no
    ## averaging.
    ##
    ## If `first.month` is greater than `last.month` the environmental variable
    ## is averaged for year t starting in `first.month` of year t - 1 and ending
    ## in `last.month` of year t. The output year corresponds to the year
    ## January occurs within the average. For example if `first.month` = 12 and
    ## `last.month` = 3, then the average for the environmental variable will
    ## occur over Dec, Jan, Feb, March and the year is specified by the year
    ## for Jan occurs in.
    ##
    ## The function outputs a data.frame with a `year` column and an `index`
    ## column. When `first.month` is greater than `last.month`, the output
    ## data.frame will have one less year than the input data.frame with no
    ## value for the minimum year in the input data frame.
    ##
    ## If `grid.id.var` is non-null, the averaging is done on a per-grid-cell
    ## basis within a year.
    ##
    ## data = a data.frame
    ## first.month = numeric giving the month to start annual average
    ## last.month = numeric giving the month to stop annual average
    ## avg.var = string, column name in `data` of the variable to average
    ## month.var = string, column name in `data` of the month variable
    ## year.var = string, column name in `data` of the year variable
    ## grid.id.var = string, column name in `data` of the grid cell id column
    ## lon.var = string, column name in `data` of for longitude, only used if
    ##           grid.id.var is non-null
    ## lat.var = string, column name in `data` of for latitude, only used if
    ##           grid.id.var is non-null
    ## type = string, value to set a `type` column in the output data.frame.
    ##        Useful if you want to rbind multiple indices together

    if(!is.data.frame(data))
        stop("Input data is not a data.frame", call. = FALSE)

    if(!first.month %in% 1:12 | !last.month %in% 1:12)
        stop("Months not between 1 and 12", call. = FALSE)

    if(!is.numeric(data[ , month.var]) | !is.numeric(data[ , year.var]))
       stop("Month variable must be numeric", call. = FALSE)

    if(first.month < last.month | first.month == last.month) {
        months <- first.month:last.month
        df <- data[data[ , month.var] %in% months, ]
    }

    if(first.month > last.month) {

        ## Remove months prior to `first.month` in first year and
        ## months after `last.month` in last year
        min.yr <- min(data[ , year.var])
        max.yr <- max(data[ , year.var])
        min.rm <- which(data[ , year.var] == min.yr &
                        data[ , month.var] < first.month)
        max.rm <- which(data[ , year.var] == max.yr &
                        data[ , month.var] > last.month)
        sub <- data[-c(min.rm, max.rm), ]

        ## Remove months not being averaged over
        months <- c(first.month:12, 1:last.month)
        sub2 <- sub[sub[ , month.var] %in% months, ]

        ## Create new year index to average over
        sp <- split(sub2, sub2[ , year.var])
        lst <- lapply(sp, function(x) {
                   x$yr.avg <- ifelse(x[ , month.var] %in% first.month:12,
                                      x[ , year.var] + 1, x[ , year.var])
                   return(x)
               })
        df <- do.call("rbind", c(lst, make.row.names = FALSE))
        df[ , year.var] <- df$yr.avg

    }

    ## Calculate averages
    if(is.null(grid.id.var)) {
        sp.avg <- split(df, df[ , year.var])
        lst.avg <- lapply(sp.avg, function(x) {
                    data.frame(year = unique(x[ , year.var]),
                               index = mean(x[ , avg.var]))
                })
        enviro <- do.call("rbind", c(lst.avg, make.row.names = FALSE))
    } else {
        sp.avg <- split(df, list(df[ , grid.id.var], df[ , year.var]))
        lst.avg <- lapply(sp.avg, function(x) {
                    data.frame(year = unique(x[ , year.var]),
                               id = unique(x[ , grid.id.var]),
                               lon = unique(x[ , lon.var]),
                               lat = unique(x[ , lat.var]),
                               index = mean(x[ , avg.var]))
                })
        enviro <- do.call("rbind", c(lst.avg, make.row.names = FALSE))
    }

    enviro$type <- type

    char.months <- c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug",
                     "sep", "oct", "nov", "dec")
    enviro$months <- paste(char.months[first.month],
                           char.months[last.month], sep = "-")
    enviro$months <- as.factor(enviro$months)

    return(enviro)
}

if(FALSE) {

    set.seed(101)
    yr <- c(rep(1950, 24), rep(1951, 24))
    id <- c(rep(1, 12), rep(2, 12), rep(1, 12), rep(2, 12))
    mt <- rep(1:12, 4)
    lt <- c(rep(48, 12), rep(50, 12), rep(48, 12), rep(50, 12))
    ln <- c(rep(120, 12), rep(122, 12), rep(120, 12), rep(122, 12))
    vl <- rnorm(48)
    df <- data.frame(year = yr, month = mt, lat = lt, lon = ln, id = id, value = vl)

    enviro.avg.months(df, 1, 12, "value", grid.id.var = NULL)
    enviro.avg.months(df, 1, 12, "value", grid.id.var = "id")

    enviro.avg.months(df, 1, 3, "value", grid.id.var = NULL)
    enviro.avg.months(df, 1, 3, "value", grid.id.var = "id")

    enviro.avg.months(df, 12, 3, "value", grid.id.var = NULL)
    enviro.avg.months(df, 12, 3, "value", grid.id.var = "id")

    enviro.avg.months(df, 1, 3, "value")
    enviro.avg.months(df, 3, 3, "value")
    enviro.avg.months(df, 12, 3, "value")
    enviro.avg.months(df, 6, 5, "value", group = "test")
}



## get.pdo -------------------------------------------------
get.pdo <- function(years = 1900:2016) {

    ## This function takes as input a range of years and downloads and processes
    ## the PDO index. The output of the function is a dataframe in 'long' format
    ## with a column for year, month, and the PDO index.
    ##
    ## years = vector of years

    if(min(years) < 1900)
        stop("Earliest PDO year is 1900")

    ## Read data
    pdo.r <- readLines("http://jisao.washington.edu/pdo/PDO.latest.txt")

    ## Trim header and footer
    start.line <- grep("^YEAR", pdo.r)
    end.line <- grep(paste0("^", max(years)), pdo.r)
    pdo.s <- pdo.r[start.line:end.line]

    ## Split strings
    pdo.c <- gsub("[/*]", "", pdo.s)
    pdo.c <- gsub("   ", "  ", pdo.c)
    pdo.c <- gsub("    ", "  ", pdo.c)
    pdo.c <- gsub("      ", "  ", pdo.c)
    pdo.c <- strsplit(pdo.c, "  ")
    pdo.c <- do.call("rbind", pdo.c)

    ## Convert to data.frame
    df.names <- trimws(tolower(pdo.c[1, ]), "both")
    pdo.d <- pdo.c[2:nrow(pdo.c), ]
    pdo.d <- apply(pdo.d, 2, function(x) as.numeric(x))
    pdo.d <- as.data.frame(pdo.d)
    names(pdo.d) <- df.names

    ## Reshape into "long" format
    months <- df.names[2:length(df.names)]
    pdo <- reshape(pdo.d, direction = "long",
                   varying = months,
                   v.names = "pdo",
                   times = months,
                   timevar = "month")
    pdo <- pdo[order(pdo$year), ]
    pdo <- pdo[pdo$year >= min(years) & pdo$year <= max(years), ]
    row.names(pdo) <- NULL
    pdo$id <- NULL

    ## Convert months to numerics
    pdo$month <- ifelse(pdo$month == "jan", 1, pdo$month)
    pdo$month <- ifelse(pdo$month == "feb", 2, pdo$month)
    pdo$month <- ifelse(pdo$month == "mar", 3, pdo$month)
    pdo$month <- ifelse(pdo$month == "apr", 4, pdo$month)
    pdo$month <- ifelse(pdo$month == "may", 5, pdo$month)
    pdo$month <- ifelse(pdo$month == "jun", 6, pdo$month)
    pdo$month <- ifelse(pdo$month == "jul", 7, pdo$month)
    pdo$month <- ifelse(pdo$month == "aug", 8, pdo$month)
    pdo$month <- ifelse(pdo$month == "sep", 9, pdo$month)
    pdo$month <- ifelse(pdo$month == "oct", 10, pdo$month)
    pdo$month <- ifelse(pdo$month == "nov", 11, pdo$month)
    pdo$month <- ifelse(pdo$month == "dec", 12, pdo$month)
    pdo$month <- as.numeric(pdo$month)

    return(pdo)
}

if(FALSE) {

    get.pdo(1950:1950)
    get.pdo(1950:2013)

    pdo <- get.pdo(1900:2016)
    head(pdo)
    tail(pdo)
    sapply(pdo, class)
    summary(pdo)

}



## sst.map.index -------------------------------------------
sst.map.index <- function(data, plot.dir = NULL, ...) {
    ## Map SST values using heatmaps
    ##
    ## This functions creates annual heat maps of SST. One file (map) is
    ## created for each year of input SST data.
    ##
    ## data = data.frame of sst data in "long" format
    ## plot.dir = directory to save plots, if null, plot is not saved
    ## ... = passed to levelplot()

    if(!is.null(plot.dir) & !dir.exists(plot.dir))
        dir.create(plot.dir, recursive = TRUE)

    years  <- unique(data$year)

    ## nep <- countriesLow[countriesLow$ADMIN %in%
    ##                     c("United States of America", "Russia", "Canada"), ]

    pb <- txtProgressBar(min = 0, max = length(years), style = 3)
    ind <- 1
    for(i in seq_along(years)) {

        dsub <- data[data$year == years[i], ]
        plot.name <- paste("sst_", years[i], ".jpeg", sep = "")

        if(!is.null(plot.dir))
            jpeg(paste(plot.dir, plot.name, sep = ""), 600, 400)

        g <- levelplot(index ~ lon * lat,
                       data = dsub,
                       cuts = 200,
                       col = "grey60",
                       contour = TRUE,
                       labels = FALSE,
                       main = as.character(years[i]),
                       ylab = "Latitude",
                       xlab = "Longitude",
                       ...,
                       panel = function(...) {
                           panel.fill(col = "grey40")
                           panel.levelplot(...)
                           # sp.polygons(nep, fill = "white", lwd = 0.7)
                           mp <- map('world2', fill = TRUE, plot = FALSE)
                           lpolygon(mp$x, mp$y,
                                    col = "white",
                                    border = "grey25")
                       })
        print(g)

        if(!is.null(plot.dir))
            dev.off()

        setTxtProgressBar(pb, ind)
        ind <- ind + 1
    }
    close(pb)
}



## sst.map -------------------------------------------------
sst.map <- function(data, plot.dir = NULL, ...) {
    ## Map SST values using heatmaps
    ##
    ## This functions creates monthly heat maps of SST. One file (map) is
    ## created for each month and year of input SST data.
    ##
    ## data = data.frame of sst data in "long" format
    ## plot.dir = directory to save plots, if null, plot is not saved
    ## ... = passed to levelplot()

    if(!is.null(plot.dir) & !dir.exists(plot.dir))
        dir.create(plot.dir, recursive = TRUE)

    years  <- unique(data$year)
    months <- unique(data$month)

    ## nep <- countriesLow[countriesLow$ADMIN %in%
    ##                     c("United States of America", "Russia", "Canada"), ]

    pb <- txtProgressBar(min = 0, max = length(years) * length(months),
                         style = 3)
    ind <- 1
    for(i in seq_along(years)) {
        for(j in seq_along(months)) {

            dsub <- data[data$year == years[i] & data$month == months[j], ]
            plot.name <- paste("sst_", years[i], "_", months[j], ".jpeg",
                               sep = "")

            if(!is.null(plot.dir))
                jpeg(paste(plot.dir, plot.name, sep = ""), 600, 500)

            g <- levelplot(sst.anom ~ lon * lat,
                           data = dsub,
                           cuts = 200,
                           col = "grey60",
                           contour = TRUE,
                           labels = FALSE,
                           main = paste("year = ", years[i],
                                        "   month = ", months[j],
                                        sep = ""),
                           ylab = "Latitude",
                           xlab = "Longitude",
                           ...,
                           panel = function(...) {
                               panel.fill(col = "grey40")
                               panel.levelplot(...)
                               # sp.polygons(nep, fill = "white", lwd = 0.7)
                               mp <- map('world2', fill = TRUE, plot = FALSE)
                               lpolygon(mp$x, mp$y,
                                        col = "white",
                                        border = "grey25")
                           })
            print(g)

            if(!is.null(plot.dir))
                dev.off()

            setTxtProgressBar(pb, ind)
            ind <- ind + 1
        }
    }
    close(pb)
}



## sst.anomaly ---------------------------------------------
sst.anomaly <- function(data, ref.years) {
    ## Calculate monthly, per grid cell, anomalies of SST
    ##
    ## Anomalies are calculated as the difference between a grid cell specific
    ## SST value for a given year/month and the long-term monthly mean (defined
    ## by ref.years) for that grid cell. This follows the methods outlined in
    ## Mueter et al. 2002, CJFAS (https://doi.org/10.1139/f02-020).
    ##
    ## data = data.frame of SST data
    ## ref.years = reference years in which to calculate the long-term mean,
    ##             should be a continuous sequence, e.g., 1950:2016

    ## make sure ref.years is a continuous sequence
    if(length(ref.years) != length(min(ref.years):max(ref.years)))
        stop("years vector is not a sequence ascending by 1")

    ## make sure ref.years are available in input data
    if(sum(ref.years %in% data$year) != length(ref.years))
        stop("ref.years not contained in input data")

    ## subset reference years from data
    ref.sst <- data[data$year >= min(ref.years) & data$year <= max(ref.years), ]

    ## calculate monthly long-term mean for each grid cell
    ## NA's are removed in the calculation of long-term mean
    mnth.avg <- aggregate(sst ~ month + id, data = ref.sst,
                          function(x) mean(x, na.rm = TRUE),
                          na.action = na.pass)
    names(mnth.avg)[names(mnth.avg) == 'sst'] <- 'long.avg'
    sst.merge <- merge(data, mnth.avg)
    sst.merge <- sst.merge[order(sst.merge$year,
                                 sst.merge$month,
                                 sst.merge$lat,
                                 sst.merge$lon), ]

    ## calculate sst anomaly
    sst.merge$anom <- sst.merge$sst - sst.merge$long.avg
    row.names(sst.merge) <- NULL
    sst <- data.frame(year = sst.merge$year,
                      month = sst.merge$month,
                      lon = sst.merge$lon,
                      lat = sst.merge$lat,
                      id = sst.merge$id,
                      sst = sst.merge$sst,
                      sst.anom = sst.merge$anom)
    return(sst)
}



## ggplot: theme_sleek -------------------------------------
# see: https://github.com/seananderson/ggsidekick
theme_sleek <- function(base_size = 11, base_family = "") {
  half_line <- base_size/2
  theme_light(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks.length = unit(half_line / 2.2, "pt"),
      strip.background = element_rect(fill = NA, colour = NA),
      strip.text.x = element_text(colour = "grey30"),
      strip.text.y = element_text(colour = "grey30"),
      axis.text = element_text(colour = "grey30"),
      axis.title = element_text(colour = "grey30"),
      legend.title = element_text(colour = "grey30", size = rel(0.9)),
      panel.border = element_rect(fill = NA, colour = "grey70", size = rel(1.1)),
      legend.key.size = unit(0.9, "lines"),
      legend.text = element_text(size = rel(0.7), colour = "grey30"),
      legend.key = element_rect(colour = NA, fill = NA),
      legend.background = element_rect(colour = NA, fill = NA),
      plot.title = element_text(colour = "grey30", size = rel(1)),
      plot.subtitle = element_text(colour = "grey30", size = rel(.85))
    )
}



## lattice: theme.mjm --------------------------------------
theme.mjm <- function(fontsize = 11, ...) {
    ## see latticeExtra::custom.theme
    ## trellis.par.get()

    col.symbol <- "#00A0EDFF"

    col.superpose     <- c("#00A0EDFF", "#DE732DFF", "#00AB06FF", "#BE6CF5FF",
                           "#A09500FF", "#F352B1FF", "#00B3A3FF")
    col.superpose.pol <- c("#00B9FFFF", "#FA8F59FF", "#0DC649FF", "#D58CFFFF",
                           "#BBAF00FF", "#FF75CAFF", "#00CCBDFF")

    col.regions <- c("#004B9FFF", "#004AA2FF", "#0049A4FF", "#0047A6FF",
                     "#0046A8FF", "#0045A9FF", "#0043ABFF", "#0042ADFF",
                     "#0041AEFF", "#003FB0FF", "#003EB1FF", "#003CB2FF",
                     "#003AB4FF", "#0039B5FF", "#0B37B6FF", "#2936B7FF",
                     "#3934B7FF", "#4633B8FF", "#5031B9FF", "#592FBAFF",
                     "#622EBAFF", "#692DBBFF", "#702BBBFF", "#772ABBFF",
                     "#7D29BCFF", "#8328BCFF", "#8827BCFF", "#8D26BCFF",
                     "#9325BCFF", "#9725BCFF", "#9C25BCFF", "#A125BBFF",
                     "#A525BBFF", "#A925BBFF", "#AD26BAFF", "#B127BAFF",
                     "#B528B9FF", "#B929B9FF", "#BC2BB8FF", "#C02DB7FF",
                     "#C32EB6FF", "#C730B5FF", "#CA32B4FF", "#CD35B3FF",
                     "#D037B2FF", "#D339B1FF", "#D63CB0FF", "#D93EAFFF",
                     "#DB41ADFF", "#DE44ACFF", "#E146ABFF", "#E349A9FF",
                     "#E54CA8FF", "#E84FA6FF", "#EA52A4FF", "#EC54A3FF",
                     "#EE57A1FF", "#F05A9FFF", "#F25D9DFF", "#F4609BFF",
                     "#F66399FF", "#F86697FF", "#FA6994FF", "#FB6D92FF",
                     "#FD7090FF", "#FE738DFF", "#FF768BFF", "#FF7988FF",
                     "#FF7C86FF", "#FF7F83FF", "#FF8280FF", "#FF867EFF",
                     "#FF897BFF", "#FF8C78FF", "#FF8F75FF", "#FF9271FF",
                     "#FF966EFF", "#FF996BFF", "#FF9C68FF", "#FF9F64FF",
                     "#FFA261FF", "#FFA65DFF", "#FFA959FF", "#FFAC56FF",
                     "#FFAF52FF", "#FFB34EFF", "#FFB64AFF", "#FFB946FF",
                     "#FFBC42FF", "#FFC03EFF", "#FFC339FF", "#FFC635FF",
                     "#FFC931FF", "#FFCD2DFF", "#FFD029FF", "#FFD325FF",
                     "#FFD722FF", "#FFDA1FFF", "#FFDD1CFF", "#FFE01BFF")
    lwd <- 0.8

    theme <- list(
        fontsize          = list(text = fontsize),
        par.main.text     = list(font = 1),
        strip.background  = list(col  = c("grey93", "grey93")),
        strip.shingle     = list(col  = c("grey65", "grey65")),
        strip.border      = list(col  = "grey50", lwd = rep(lwd, 7)),
        axis.components   = list(right  = list(tck = 0.5),
                                 top    = list(tck = 0.5),
                                 left   = list(tck = 0.5),
                                 bottom = list(tck = 0.5)),
        axis.line         = list(col = "grey50", lwd = lwd),
        plot.symbol       = list(pch = 16, col = col.symbol),
        plot.line         = list(pch = 16, col = "grey10"),
        add.line          = list(col = "grey10"),
        superpose.symbol  = list(col = col.superpose, pch = 16),
        superpose.line    = list(col = col.superpose),
        regions           = list(col = col.regions),
        dot.symbol        = list(col = col.symbol),
        plot.polygon      = list(col = "grey50", border = "white"),
        superpose.polygon = list(col = col.superpose.pol, border = "white"),
        box.rectangle     = list(col = "grey40"),
        box.dot           = list(col = "grey10"),
        box.umbrella      = list(col = "grey40", lty = 1),
        box.3d            = list(col = "grey70", lwd = lwd))
    modifyList(modifyList(standard.theme("pdf"), theme), simpleTheme(...))
}


## moving average_df
moving_average_df <- function(x, value, lag = 2) {  # Create user-defined function
  ma=NA
  ma_sd=NA
  for(t in c(lag+1):c(nrow(x))){
    x1=x[c(t-lag):c(t+lag), value] 
    ma[t]=mean(as.vector(x1))
    ma_sd[t]=sd(as.vector(x1))
  }
  return(cbind(x, ma, ma_sd))
}


## ocean_region_col

ocean_region_col <- function(x, stock.col="Stock", region="Ocean.Region2"){
  # Add an 'Ocean.Region' column to a dataframe that lists stock names

  stk <- x[[stock.col]]
  lookup <- data.frame(Stock = sock.info$Stock, region = sock.info[[region]])
  stk <- data.frame(Stock = stk)
  
  df <- left_join(stk, lookup, by="Stock")
  x$region <- df$region
  return(x)
}
  
