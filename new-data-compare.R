## Compare old and new AK datasets

new.info <- read.csv(here("data-downloaded/stock_info2025-09-25.csv")) |>
  mutate(stock.name = paste(stock.name, "1"))
old.info <- read.csv(here("data-downloaded/stock_info2024-01-12.csv")) |>
  mutate(stock.name = paste(stock.name, "2"))


new.dat <- read.csv(here('data-downloaded/salmon_productivity_compilation2025-09-25.csv')) |>
  left_join(new.info[,c("stock.id", "ocean.basin")], by="stock.id") |>
  mutate(stock.name = paste(stock, "1"))
old.dat <- read.csv(here('data-downloaded/salmon_productivity_compilation2024-01-12.csv')) |>
  left_join(old.info[,c("stock.id", "ocean.basin")], by="stock.id") |>
  mutate(stock.name = paste(stock, "2"))


## Pink - Odd

pink.odd.stk.new <- new.info |>
  filter(species == "Pink-Odd",
         state == "AK") |>
  select(stock.id, stock.name, species, ocean.basin)

pink.odd.stk.old <- old.info |>
  filter(species == "Pink-Odd",
         state == "AK") |>
  select(stock.id, stock.name, species, ocean.basin)


pink.odd <- bind_rows(filter(new.dat,
                             species=="Pink-Odd",
                             stock.id %in% pink.odd.stk.new$stock.id),
                      filter(old.dat,
                             species=="Pink-Odd",
                             stock.id %in% pink.odd.stk.old$stock.id),
                      .id = "new_dat1")

pink.odd$ocean.basin <- gsub("WC", "NC", pink.odd$ocean.basin)
pink.odd$ocean.basin <- gsub("GOA", "GoA", pink.odd$ocean.basin)



pink.odd |>
  #filter(new_dat1==1) |>
  ggplot() +
  geom_line(aes(x=broodyear, y=logRS, group=stock, col=new_dat1), alpha=0.5) +
  facet_wrap(vars(ocean.basin), nrow=3)


## Pink - even

pink.even.stk.new <- new.info |>
  filter(species == "Pink-Even",
         state == "AK") |>
  select(stock.id, stock.name, species, ocean.basin)

pink.even.stk.old <- old.info |>
  filter(species == "Pink-Even",
         state == "AK") |>
  select(stock.id, stock.name, species, ocean.basin)


pink.even <- bind_rows(filter(new.dat,
                             species=="Pink-Even",
                             stock.id %in% pink.even.stk.new$stock.id),
                      filter(old.dat,
                             species=="Pink-Even",
                             stock.id %in% pink.even.stk.old$stock.id),
                      .id = "new_dat1")

pink.even$ocean.basin <- gsub("WC", "NC", pink.even$ocean.basin)
pink.even$ocean.basin <- gsub("GOA", "GoA", pink.even$ocean.basin)



pink.even |>
  #filter(new_dat1==1) |>
  ggplot() +
  geom_line(aes(x=broodyear, y=logRS, group=stock, col=new_dat1), alpha=0.5) +
  facet_wrap(vars(ocean.basin), nrow=3)



## Chum

chum.stk.new <- new.info |>
  filter(species == "Chum") |>
  select(stock.id, stock.name, species, ocean.basin)

chum.stk.old <- old.info |>
  filter(species == "Chum") |>
  select(stock.id, stock.name, species, ocean.basin)


chum <- bind_rows(filter(new.dat,
                              species=="Chum",
                              stock.id %in% chum.stk.new$stock.id),
                       filter(old.dat,
                              species=="Chum",
                              stock.id %in% chum.stk.old$stock.id),
                       .id = "new_dat1")

#chum$ocean.basin <- gsub("NC", "WC", chum$ocean.basin)
chum$ocean.basin <- gsub("SC", "WC", chum$ocean.basin)
chum$ocean.basin <- gsub("GOA", "GoA", chum$ocean.basin)



chum |>
  #filter(new_dat1==1) |>
  ggplot() +
  geom_line(aes(x=broodyear, y=logRS, group=stock, col=new_dat1), alpha=0.5) +
  facet_wrap(vars(ocean.basin), nrow=4)


