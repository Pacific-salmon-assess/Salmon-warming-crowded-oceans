Work in progress repository with code and data to undertake a comparative analysis of spatial and temporal variation in Pacific salmon responses to a warming and more crowded North Pacific Ocean. 

This repository currently just has data and models fit to sockeye spawner-recruitment data, and largely mirrors the analyses origionally described in [Connors et al. 2020](https://cdnsciencepub.com/doi/10.1139/cjfas-2019-0422#:~:text=In%20the%20south%2C%20a%20warm,effects%20of%20competition%20at%20sea.).
These analyses will be expanded to include additional species and models that allow relationships to vary over time and space and across species and life histories (e.g., as in [Malick 2020](https://onlinelibrary.wiley.com/doi/abs/10.1111/fog.12469)).

## Files
- `make.R`: source this file to reproduce the project.

- `load.R`: load packages and scripts necessary for analysis. 

- `functions.R`: all functions written for the analysis are placed in this file.
  
- `data_download.R`: download data needed for project and write to CSV.

- `sock_data_clean.R`: clean/process the raw downloaded sockeye data.
  
- `sst_import_process.R`: read in SST data, generate anomalies and calculate average SST over a specified period and region.

- `sst_explore.R`: explore the SST anomalies. 

- `sock_covariates.R`:  create stock specific covariates used in analysis.

- `sock_explore.R`:  exploratory graphics and summaries of sockeye and covariate data.

- `single_stock.R`:  exploratory single-stock generalized Ricker models.

- `hbm_fit.R`:  fit stationary hierarchical bayesian models.

- `pub.R`:  generate all tables and figures.

- `background_plots.R`:  background figures.

- `recruitment_means.R`: script for estimating age @ ocean entry proportions for stocks without detailed age composition data
