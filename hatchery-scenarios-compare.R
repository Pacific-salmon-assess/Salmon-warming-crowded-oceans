#predicted % change in R/S for sockeye and chum for recent decade average pink salmon hatchery production across range of SSTs - 1) long-term average and 2) largest magnitude positive observed


# need: long-term average sst, recent (heatwave?) sst, recent pinks
# SST-prod relationships (coefficients)
# SST coefs (% change per degree C) x degrees above average + Comp coefs x million pinks

# 1) comp coef * SD unit increase in pinks + sst coef
# 2) comp coef * SD unit increase in pinks + sst coef * degree C increase

comp.raw <- read.csv(here("data-downloaded/competitor_indices_2024.csv"))

pink.increase <- comp.raw |> select(Year, pink_numbers_np) |>
  mutate(std = (pink_numbers_np-mean(pink_numbers_np))/sd(pink_numbers_np)) |>
  filter(Year > 2006) |>
  summarize(recent_pinks_std = mean(std))

sst.raw <- read.csv(here("data/sst_yr_1_stock_anomalies.csv"))
sst.increase <- sst.raw |> filter(Species != "Pink") |>
  left_join(map.info[,1:3], by="Stock.ID") |>
  group_by(Stock.ID) |>
  mutate(std = sst_anomaly / sd(sst_anomaly)) |>
  filter(Year > 2006) |>
  ungroup() |>
  summarize(recent_sst_stnd = mean(std), .by=Ocean.Region2)

box |>
  filter(par %in% c("mu_gamma3", "mu_kappa3"),
         sp.id == "sockeye") |>
  left_join(sst.increase, by="Ocean.Region2") |>
  left_join(pink.increase, by="Ocean.Region2") |>
  mutate(p_change_tot = case_when(
    varnam == "SST" ~ x*recent_sst_stnd,
    varnam == "Competitors" ~ x*mean_comp)) |>
  ggplot() +
  geom_boxplot(aes(x=sp.id, y=p_change_tot, col=Ocean.Region2))


