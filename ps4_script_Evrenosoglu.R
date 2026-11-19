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

# First Date: Jan 2002
# Last Date: Dec 2019
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

# Test 2: ARIMA vs seasonal random walk
forecast::dm.test(liquor_errs$arima, liquor_errs$snaive, h = 1, power = 2)

# Test 1 --> DM: -2.9062, p-value: 0.004042
# Test 2 --> DM: -10.755, p-value< 0.2.2e-16

# Both reject the null hypothesis of equal accuracy, both DM test signs are negative, so ARIMA has a  smallers squared one-step forecast error than ETS 
# As the liquor store chain, I would use ARIMA. It beat both alternatives on 216 rolling forecast origins, and its advantage over ETS is statistically significant

#3a
liquor_combo_fc <- liquor_cv_fit |>
  mutate(combo = (arima + ets) / 2) |>
  forecast(h = 1)

liquor_combo_fc |>
  accuracy(liquor, measures = list(RMSE = RMSE, MAE = MAE)) |>
  select(.model, RMSE, MAE) |>
  arrange(RMSE)

# It looks like combo is in between the arima and the ets. So it improve on the ets, but drags down the arima model.

cor(liquor_errs$arima, liquor_errs$ets)

# We got a correlation of 0.878 which means the ARIMA and ETS forecast errors are very similar.
#The combination did not beat the ARIMA model, because the ARIMA and ETS forecasts are highly correlated, meaning that if one overforecasts, the other one will as well.
# The DM test in 2.4 showed that the ETS model was less accurate. Giving ETS half the weight pulls the combination toward the weaker model without adding anything that the ARIMA doesn't have.

#4.1a
qd_raw <- read_csv("data/fredqd_2026-09.csv", show_col_types = FALSE) |>
  slice(-(1:2))

growth <- function(x) 100 * difference(log(x))

qd <- qd_raw |>
  mutate(date = as.Date(sasdate, "%m/%d/%Y")) |>
  transmute(date,
            gdp    = growth(GDPC1),     # real GDP
            cons   = growth(PCECC96),   # real consumption
            inv    = growth(GPDIC1),    # real investment
            emp    = growth(PAYEMS),    # payroll employment
            infl   = growth(GDPCTPI),   # GDP price index
            unrate = difference(UNRATE),     # unemployment rate
            ffr    = difference(FEDFUNDS),   # federal funds rate
            gs10   = difference(GS10)) |>    # 10-year Treasury yield
  filter(date >= as.Date("1960-03-01"), date <= as.Date("2019-12-01"))

y <- qd |> select(-date) |> as.matrix()

nrow(y)
head(y)

# 4.1c
N <- nrow(y)
origins <- 120:(N - 1)
actual <- y[121:N, ]

# We get 120 forecasts, we start with Q1 1990 and end with Q4 2019

#4.2
