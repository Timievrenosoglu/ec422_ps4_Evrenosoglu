library(tidyverse)
library(fpp3)
library(BVAR)

liquor_raw <- read_csv(
  "data/MRTSSM4453USN.csv",
  na = ".",
  col_types = cols(
    observation_date = col_date(format = "%Y-%m-%d"),
    MRTSSM4453USN = col_double()
  )
)

#2.1
liquor <- liquor_raw |>
  transmute(Month = yearmonth(observation_date), sales = MRTSSM4453USN) |>
  as_tsibble(index = Month) |>
  filter(Month <= yearmonth("2019 Dec"))

has_gaps(liquor)

# plot sales & log(sales)
liquor |>
  autoplot(sales) +
  labs(title = "Monthly sales",
       x = "Month",
       y = "Sales")

liquor |>
  mutate(log_sales = log(sales)) |>
  gg_season(log_sales, labels = "both") +
  labs(title = "Seasonal plot of log sales",
       x = "Month",
       y = "log(Sales)")

#It looks like Decembers are always the peaks, and Januarys are always the lows. It also looks like the seasonal swings get larger as time goes on, suggesting that we need to log the sales data. Using logs standardizes the swings so they are the same size at all dates, this can be observed in out gg_seasonal plot where the swings are of similar magnitude.

#2.2a
liquor_fit <- liquor |>
  model(arima = ARIMA(log(sales)),
        ets   = ETS(log(sales)))

liquor_fit |> select(arima) |> report()
liquor_fit |> select(ets) |> report()

#Chosen Models:ARIMA(1,0,2)(2,1,1)[12] w/ drift & ETS(A,A,A)
# d= 0 means there was no ordinary difference, D = 1 means there was one seasonal difference made
# Error = A (additive), Trend = A (additive), Seasonal = A (additive)

#2.2b
liquor_fit |> select(arima) |> gg_tsresiduals()

liquor_fit |>
  select(arima) |>
  augment() |>
  features(.innov, ljung_box, lag = 24, dof = 6)

# No it doesn't look like white noise, our p value is extremely small, meaning that we can reject the null that the residuals are white noise

#2.3a
liquor_cv <- liquor |>
  stretch_tsibble(.init = 120, .step = 1) |>
  filter(.id < max(.id))

liquor_cv_fit <- liquor_cv |>
  model(snaive = SNAIVE(sales),
        ets    = ETS(log(sales)),
        arima  = ARIMA(log(sales)))

saveRDS(liquor_cv_fit, "results/liquor_cv_fit.rds")

# This TS-CV has 216 forecast origins

# 2.3b
liquor_cv_fc <- liquor_cv_fit |>
  forecast(h = 1)

liquor_cv_fc |>
  accuracy(liquor, measures = list(RMSE = RMSE, MAE = MAE))

# First Date:
# Last Date:  
#### CHECK INSTRUCTIONS

#2.4
liquor_errs <- liquor_cv_fc |>
  as_tibble() |>
  left_join(as_tibble(liquor) |> select(Month, actual = sales),
            by = "Month") |>
  mutate(e = actual - .mean) |>
  select(Month, .model, e) |>
  pivot_wider(names_from = .model, values_from = e)

forecast::dm.test(liquor_errs$arima, liquor_errs$ets, h = 1, power = 2)

#


