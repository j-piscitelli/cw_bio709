#' DESCRIPTION:
#' Script for time-series

# in-class ----------------------------------------------------------------
pacman::p_load(tidyverse,
               forecast,
               lterdatasampler,
               daymetr,
               glarma)

url <- "https://raw.githubusercontent.com/aterui/biostats/master/data_raw/data_ts_anormaly.csv"
(df_ts <- read_csv(url))

# Plot time-series anomalies
df_ts %>% 
  ggplot(aes(x = year,          # Map 'year' to the x-axis
             y = anormaly)) +  # Map 'anormaly' to the y-axis
  geom_line() +                  # Add a line connecting the points
  geom_point() +                 # Add points at each observation
  theme_bw() +                   # Use a clean black-and-white theme
  labs(
    x = "Year",                     # Label x-axis
    y = "Anomaly"                   # Label y-axis
  )

# Air quality anomaly appears to increase over time...
m_lm <- lm(anormaly ~ year,
           data = df_ts)
summary(m_lm)
# ... and a linear model appears to support this conclusion,
# but the data points aren't independent. If each year varies up or down from the
# previous year by a random amount, it will trend one way or another by chance.
# Because each year is correlated with its temporal neighbors, the assumption of
# independence is false.

# Random walk
y <- NULL
y[1] <- 0

for (i in 1:99) {
  y[i + 1] <- y[i] + rnorm(1, mean = 0, sd = 1)
}

tibble(y = y,
       x = 1:length(y)) %>% 
  ggplot(aes(x = x,
             y = y)) +
  geom_point() +
  geom_line()

# Temporal autocorrelation isn't just a problem if the X axis is time. With any
# time series data, the effect of temporal confounding must be accounted for.


## Autoregressive (AR) models
# AR models assume stationarity: their expected values and variance are constant
# over time. 

# Plot time series of Lake Huron water level
df_huron <- tibble(year = time(LakeHuron),
                   water_level = as.numeric(LakeHuron)) %>% 
  arrange(year)

df_huron %>% 
  ggplot(aes(x = year,
             y = water_level)) +
  geom_point(alpha = 0.25) +
  geom_line(linetype = "dotted") +
  geom_smooth(method = "lm",
              color = "black",
              linewidth = 0.5) +
  theme_bw()

# Autoregressive model using Arima() function
# This assumes that the data is presented in temporal order.
# If it is not, it must be arranged by time first.
m_ar1 <- Arima(
  df_huron$water_level,
  order = c(
    1, # = order, the number of time steps that influence is expect to last
    0,
    0)
)

df_huron_ar1 <- df_huron %>% 
  mutate(fit = fitted(m_ar1) %>% 
           as.numeric())

# Plot observed and fitted values
df_huron_ar1 %>% 
  ggplot() +
  geom_point(aes(x = year, 
                 y = water_level),
             alpha = 0.25) +        # Plot observed water levels
  geom_line(aes(x = year, 
                y = fit),           # Plot AR(1) fitted values
            color = "steelblue") +
  theme_bw()                        


# The coefficient for ar1 is 0.8375: this shows considerable influence of
# the past time.

## Moving average (MA) models
# Instead of working from past values of the response variable, MA models
# consider the error from the past several (q) time steps.
m_ma1 <- Arima(
  df_huron$water_level,
  order = c(
    0,
    0,
    1) # use the third argument for the MA order
)

## ARMA models
# You can simply include an AR and an MA together and get an ARMA

m_arma <- Arima(
  df_huron$water_level,
  order = c(
    1, # autoregressive order
    0,
    1) # MA order
)

## ARIMA models
# ARIMA models look for autocorrelation in the *difference* between time steps

m_arma <- Arima(
  df_huron$water_level,
  order = c(
    1, # autoregressive order
    1, # order of difference
    0) 
)

## Model selection
# There are many possible models, since each of the three order values can be
# anything from 0 to the number of time steps in the data.
# There is a function automate model selection based on AICs
auto.arima(
  y = df_huron$water_level,
  stepwise = FALSE, # Search all possible models
  ic = "aic" # Use AIC for model selection
)


## ARIMAX
# ARIMAX adds external predictors to ARIMA.
# Instead of looking for a single mean that the measurements move around, we
# often want to see the trend related to some other variable after accounting for
# temporal autocorrelation.

# Length of ice cover on Lake Mendota
data("ntl_icecover")

df_ice <- ntl_icecover %>% 
  as_tibble() %>% 
  filter(between(year, 1980, 2014),
         lakeid == "Lake Mendota") %>% 
  arrange(year) # Make sure data is in order

df_ice %>% 
  ggplot(aes(x = year,
             y =  ice_duration)) +
  geom_line(linetype = "dotted") +
  geom_point(alpha = 0.25) +
  theme_bw()

# We want to determine the effect of temperature on ice duration

# We can download climate data from the internet! Amazing!
list_mendota <- download_daymet(
  site = "Mendota",
  lat = 43.1,
  lon = -89.4,
  start = 1980,
  end = 2024,
  internal = TRUE
)

# Get annual minimum temperature in Celsius and add to df_ice
df_temp <- list_mendota$data %>% 
  as_tibble() %>% 
  janitor::clean_names() %>% 
  mutate(
    date = as.Date(
      paste(year, yday, sep = "-"),
      format = "%Y-%j"
    ),
    month = month(date)
  ) %>% 
  arrange(year, yday) %>% 
  group_by(year) %>% 
  summarize(temp_min = round(mean(tmin_deg_c), 2 ))

df_ice <- df_ice %>% 
  left_join(df_temp,
            by = "year")

# We can't just regress ice duration by minimum temperature, due to the possible
# temporal autocorrelation. Instead:
obj_arima <- auto.arima(
  df_ice$ice_duration, # Response variable
  xreg = df_ice$temp_min, # Variable whose effect we're interested in
  stepwise = FALSE
)

confint(obj_arima)

# After accounting for autocorrelation, we have about 8 days fewer of ice coverage
# for every increase of 1 degree Celsius. The standard error is 2.67, so the 95%
# confidence interval is 5.34 days on either side of. So we can be pretty sure
# there is a negative effect of minimum temperate on ice duration.


# lab ---------------------------------------------------------------------

# ============================================================
# EXERCISE: Bison Body Mass, Climate, and Time-Series Analysis
# ============================================================

library(lterdatasampler)

# The "knz_bison" dataset contains long-term monitoring data
# on bison captured at Konza Prairie Biological Station.
#
# ------------------------------------------------------------
# Key columns may include:
# rec_year      : Year of capture
# animal_sex    : Sex of the individual (e.g., female, male)
# animal_weight : Body mass of bison
# ------------------------------------------------------------
#
# In this exercise, you will explore long-term trends in bison
# body mass and evaluate how climate variability may influence
# weight dynamics over time.

# 1. Explore the structure of the knz_bison dataset.
#    - Inspect variable types and missing values.
#    - Reformat variables as needed for analysis.

df_bison <- knz_bison %>%
  drop_na() %>% 
  mutate(date = as.Date(
    paste(rec_year, rec_month, rec_day,
          sep = "-")
  ))

# 2. Subset the data to include observations from 1994–2012.
df_bison <- df_bison %>% 
  filter(between(rec_year, 1994, 2012))

# 3. Calculate the average body mass for female and male bison
#    for each year in the selected time period.
df_bison_mean_wt <- df_bison %>% 
  group_by(rec_year, animal_sex) %>% 
  mutate(mean_wt = round(mean(animal_weight), 2)) %>% 
  ungroup() %>% 
  select(c(rec_year, animal_sex, mean_wt)) %>% 
  distinct()

# 4. Obtain climate data from the daymetr dataset.
#    - Identify relevant climate variables (e.g., temperature,
#      precipitation).
#    - Associate climate data with knz_bison by year.
#    - Coordinates: Lat 39.09300	Lon -96.57500

list_knz <- download_daymet(site = "Konza",
                            lat = 39.093,
                            lon = -96.575,
                            start = 1994,
                            end = 2012,
                            internal = TRUE)

df_precip <- list_knz$data %>% 
  as_tibble() %>% 
  janitor::clean_names() %>% 
  mutate(
    date = as.Date(
      paste(year, yday, sep = "-"),
      format = "%Y-%j"
    ),
    month = month(date)
  ) %>% 
  arrange(year, yday) %>% 
  group_by(year) %>% 
  summarize(ann_precip = (sum(prcp_mm_day))) %>%
  rename("rec_year" = year)


df_bison_mean_wt <- left_join(df_bison_mean_wt,
                              df_precip,
                              by = "rec_year")

# 5. Perform a time-series analysis to examine whether selected
#    climate variables influence annual bison body mass.
#    - Consider temporal autocorrelation and lag effects.
#    - Model males and females separately


df_bison_mean_wt_f <- df_bison_mean_wt %>% 
  filter(animal_sex == "F")

df_bison_mean_wt_m <- df_bison_mean_wt %>% 
  filter(animal_sex == "M")

obj_arima_f <- auto.arima(df_bison_mean_wt_f$mean_wt,
                          xreg = df_bison_mean_wt_f$ann_precip,
                          stepwise = FALSE,
                          d = 0)

obj_arima_m <- auto.arima(df_bison_mean_wt_m$mean_wt,
                          xreg = df_bison_mean_wt_m$ann_precip,
                          stepwise = FALSE,
                          d = 0)

# 6. Using your fitted model, compare observed bison body mass
#    with predicted values for the period 2014–2020.
#    - Evaluate model performance and discuss sources of uncertainty.
