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

#3b
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
var_ols_fc <- function(z, p) {
  e <- embed(z, p + 1)            # columns: y_t, then its lags 1 to p
  M <- ncol(z)                    # number of variables
  B <- coef(lm(e[, 1:M] ~ e[, -(1:M)]))
  x <- c(1, t(z[nrow(z):(nrow(z) - p + 1), ]))   # constant, latest p lags
  x %*% B
}

fc_var1 <- NULL
fc_var4 <- NULL
for (T in origins) {
  fc_var1 <- rbind(fc_var1, var_ols_fc(y[1:T, ], 1))
  fc_var4 <- rbind(fc_var4, var_ols_fc(y[1:T, ], 4))
}

rmse_by_var <- function(f) sqrt(colMeans((actual - f)^2))

rmse_by_var(fc_var1)
rmse_by_var(fc_var4)


# The VAR(1) forecasts better because it has a lower RMSE for 7/8 variables
# Even though VAR(4) contains VAR(1), VAR(4) has to estimate its extra coefficients from data. Estimating them with error adds noise to its forecasts.

# 4.3a
bvar_cv <- function(lambda) {
  pr <- bv_priors(hyper = "lambda",
                  mn = bv_mn(lambda = bv_lambda(mode = lambda, sd = 0.4,
                                                min = 0.999 * lambda,
                                                max = 1.001 * lambda),
                             b = 0))
  fc <- NULL
  for (T in origins) {
    fit_T <- bvar(y[1:T, ], lags = 4, priors = pr,
                  n_draw = 2000, n_burn = 1000, verbose = FALSE)
    x <- c(1, t(y[T:(T - 3), ]))         # constant and 4 lags of each variable
    fc <- rbind(fc, x %*% coef(fit_T))   # forecast of y_{T+1}
  }
  fc
}

set.seed(422)
fc_bvar_03 <- bvar_cv(0.3)
saveRDS(fc_bvar_03, "results/fc_bvar_03.rds")

#4.3b
set.seed(422)
fc_bvar_005 <- bvar_cv(0.05)

var_results <- tibble(
  Variable = colnames(y),
  `VAR(1)` = rmse_by_var(fc_var1),
  `VAR(4)` = rmse_by_var(fc_var4),
  `BVAR(4), lambda = 0.3` = rmse_by_var(fc_bvar_03),
  `BVAR(4), lambda = 0.05` = rmse_by_var(fc_bvar_005)
)

saveRDS(var_results, "results/var_results.rds")

var_results |> knitr::kable(digits = 3)

# BVAR(4) w/ lambda 0.3 has a lower RMSE for every variable compared to VAR(4) 
# Because the BVAR(4) and OLS VAR(4) share the same variables and lags, any difference in accuracy comes entirely from the Minnesota prior, which shrinks the coefficients toward a prior mean. In this case it made the BVAR more accurate
# OLS Var(1) has 3/8 variables that have a lower RMSE than BVAR(4) w/ lambda 0.3.
# With lambda 0.05 the RMSE is higher for every variable except for one when compared to bvar with lambda 0.3. This tells us that we maybe overshrunk
# No, the winner's RMSE is not a fair estimate of how accurately it will forecast in the future b/c we are the ones choosing the lambda which minimizes the RMSE, this is the model that worked the best on this specifc data (1990-2019)

#4.3c
e_var4 <- actual[, "gdp"] - fc_var4[, 1]
e_bvar <- actual[, "gdp"] - fc_bvar_03[, 1]

forecast::dm.test(e_var4, e_bvar, h = 1, power = 2)

#Repeat
e_var4 <- actual[, "infl"] - fc_var4[, 1]
e_bvar <- actual[, "infl"] - fc_bvar_03[, 1]

forecast::dm.test(e_var4, e_bvar, h = 1, power = 2)

# We can see that for inflation our DM stat is 2.13 and p value = 0.034
# We can see that for inflation our DM stat is 3.02 and p value = 0.003
# In both tests DM is positive, so the OLS VAR(4) had the larger squared errors and the BVAR with lambda 0.3 did better on both variables.
# In relation to 4.3b, the BVAR’s lower RMSE on every variable doesn’t mean it’s meaningfully better on every variable.

#5.1a
md_raw <- read_csv("data/fredmd_2026-09.csv", show_col_types = FALSE)
tcode <- unlist(md_raw[1, -1])
md <- md_raw[-1, ] |> mutate(date = as.Date(sasdate, "%m/%d/%Y"))

transform_series <- function(x, code) {
  d <- function(v) c(NA, diff(v))
  switch(as.character(code),
         "1" = x,                      # level
         "2" = d(x),                   # first difference
         "3" = d(d(x)),                # second difference
         "4" = log(x),                 # log
         "5" = 100 * d(log(x)),        # growth rate
         "6" = 100 * d(d(log(x))),     # change in the growth rate
         "7" = d(x / dplyr::lag(x) - 1))
}

#5.1b
fm <- map2_dfc(md[names(tcode)], tcode, transform_series) |>
  mutate(date = md$date, .before = 1) |>
  filter(date >= as.Date("1960-03-01"), date <= as.Date("2019-12-01"))

fm <- fm[, colSums(is.na(fm)) == 0]

dim(fm)

#5.1c
tcode[c("INDPRO", "PAYEMS", "CPIAUCSL", "FEDFUNDS")]

#INDOPRO is transformed with growth rate, which is just taking the long & the first difference
#PAYEMS is transformed with growth rate, which is just taking the long & the first difference
#CPIAUCSL is transformed with change in growth rate
#FEDFUNDS is transformed with first difference
#FEDFUNDS and INDOPRO are different because INDPRO is an index of physical output that grows by roughly a constant percentage over time.

#5.2a
Z <- fm |> select(-date) |> as.matrix()

pc <- prcomp(scale(Z))

summary(pc)$importance[, 1:5]

#The first three PC factor for around 29% of all variance

#5.2b
F1 <- pc$x[, 1] * sign(cor(pc$x[, 1], fm$INDPRO))

tibble(date = fm$date, F1 = F1) |>
  ggplot(aes(date, F1)) +
  geom_line() +
  labs(title = "First Factor from FRED-MD", x = NULL, y = NULL) +
  theme_minimal()

#5.2c
sort(abs(cor(pc$x[, 1], Z))[1, ], decreasing = TRUE)[1:8]
sort(abs(cor(pc$x[, 2], Z))[1, ], decreasing = TRUE)[1:8]

# The most correlated series are employment (PAYEMS, USGOOD, MANEMP) and industrial production and capacity utilization (INDPRO, IPMANSICS, CUMFNS). This measures the strength of the business cycle, and in the plot it moves sharply during recessions. 
# The second factor is a price factor. This factor loads on CPI and PCE price indexes, especially goods, commodity and transportation prices, so it captures changes in inflation driven largely by energy and goods prices, with spikes around oil price shocks.

#5.3
ip <- fm$INDPRO
n <- nrow(Z)
init <- which(fm$date == as.Date("1989-12-01"))

e_ar <- NULL
e_fa <- NULL
for (T in init:(n - 1)) {
  F_T <- prcomp(scale(Z[1:T, ]))$x[, 1:3]         # factors using data through T
  d <- data.frame(y_next = ip[2:T], y = ip[1:(T - 1)], F_T[1:(T - 1), ])
  new <- data.frame(y = ip[T], t(F_T[T, ]))
  e_ar <- c(e_ar, ip[T + 1] - predict(lm(y_next ~ y, d), new))
  e_fa <- c(e_fa, ip[T + 1] - predict(lm(y_next ~ ., d), new))
}

c(RMSE_AR = sqrt(mean(e_ar^2)), RMSE_FA = sqrt(mean(e_fa^2)))

# I can confirm that the RMSE_AR is 0.617 & RMSE_FA is 0.589
# It would be cheating because the factors would contain information from the future, resulting in information leakage

#5.4
direct_race <- function(series, h = 12, r = 3) {
  lv <- md |>
    mutate(log_x = log(.data[[series]]),
           g = 100 * (log_x - lag(log_x))) |>
    filter(date >= as.Date("1960-03-01"), date <= as.Date("2019-12-01"))
  yh <- 100 * (lead(lv$log_x, h) - lv$log_x)      # growth over the next h months
  g  <- lv$g                                      # growth over the latest month

  e_ar <- NULL
  e_fa <- NULL
  for (T in init:(n - h)) {
    F_T <- prcomp(scale(Z[1:T, ]))$x[, 1:r, drop = FALSE]
    i <- 1:(T - h)
    d <- data.frame(yh = yh[i], g = g[i], F_T[i, , drop = FALSE])
    new <- data.frame(g = g[T], F_T[T, , drop = FALSE])
    e_ar <- c(e_ar, yh[T] - predict(lm(yh ~ g, d), new))
    e_fa <- c(e_fa, yh[T] - predict(lm(yh ~ ., d), new))
  }
  tibble(series, forecasts = length(e_ar),
         RMSE_AR = sqrt(mean(e_ar^2)), RMSE_FA = sqrt(mean(e_fa^2)),
         MAE_AR = mean(abs(e_ar)), MAE_FA = mean(abs(e_fa)),
         DM_pvalue = forecast::dm.test(e_ar, e_fa, h = h, power = 2,
                                     varestimator = "bartlett")$p.value)
}

direct_race("CPIAUCSL")

#The RMSE shrunk around 10.4%. Furthermore, the DM test, adjusted for the autocorrelation in overlapping 12 month forecast errors, gives p = 0.0004, so the improvement is statistically significant.
