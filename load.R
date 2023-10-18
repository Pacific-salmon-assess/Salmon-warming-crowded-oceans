## (1) Load functions and libraries for analysis ##
## --------------------------------------------- ##

library(lattice)
library(chroma) # see: https://github.com/michaelmalick/r-chroma
library(nlme)
library(car) # needed for vif() function
library(maps)
library(plyr)
library(dplyr)
library(reshape2)
library(rworldmap)
library(ersst)     # see: https://github.com/michaelmalick/r-ersst
library(rsoi)
library(here)
require(fields)
require(plotrix)
require(data.table)
data(countriesLow)

## bayesian model packages
library(parallel)
library(coda)
library(rjags)
#library(jagstools)
library(codatools) # see: https://github.com/MichaelMalick/r-codatools/
library(loo)
library(rstan)
library(bayesplot)
library(ggplot2)
library(stringr)


## Install packages from github
#library(devtools)
#install_github(repo = "michaelmalick/r-chroma")
#install_github(repo = "michaelmalick/r-jagstools")
#install_github(repo = "michaelmalick/r-codatools")
#install_github(repo = "michaelmalick/r-ersst")


## Create directories
if(!dir.exists("./figures/"))
    dir.create("./figures/")
if(!dir.exists("./output/"))
    dir.create("./output/")


## Load saved HBM outputs
#if(dir.exists("./output/"))
#    load_rdata(path = "./output/", verbose = FALSE)
#if(dir.exists("./output/models/"))
#    load_rdata(path = "./output/models/", verbose = TRUE)


## Data
sock <- read.csv("./data/master_brood_table_covar.csv",
                 stringsAsFactors = FALSE)
sock$Stock <- factor(sock$Stock, levels = unique(sock$Stock))

sock.info <- read.csv("./data/master_stock_info.csv",
                      stringsAsFactors = FALSE)
sock.info$Stock <- factor(sock.info$Stock, levels = unique(sock.info$Stock))


## Set bayesplot theme
bayesplot::bayesplot_theme_set(new = theme_sleek())
