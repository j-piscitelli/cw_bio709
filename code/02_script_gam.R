#' DESCRIPTION:
#' Script for GAMs

# in-class ----------------------------------------------------------------
pacman::p_load(tidyverse,
               ggeffects,
               mgcv)

# lab ---------------------------------------------------------------------

## water temperature data
link <- "https://raw.githubusercontent.com/aterui/biostats/master/data_raw/data_water_temp.csv"

(df_wt_raw <- read_csv(link))

# Check column classes
sapply(df_wt_raw, class)

# Extract dates from date_time column and filter to March-October 2022
df_wt <- df_wt_raw %>% 
  mutate(date = as.Date(date_time,
                        format = "%m/%d/%Y"),
         year = year(date),
         month = month(date)) %>% 
  filter(year == 2022,
         between(month, 3, 10))

# Aggregate data by day of observation
df_wt_daily <- df_wt %>% 
  group_by(date, site) %>% 
  summarize(temp = mean(temp, na.rm = TRUE) %>% 
              round(3),
            .groups = "drop")

# Visualize daily water temperature at each site
df_wt_daily %>% 
  ggplot(aes( x = date,
              y = temp,
              color = site)) +
  geom_point(alpha = 0.25) +
  theme_bw() +
  labs( x = "Date",
        y = "Water Temperature",
        color = "Wetland Type")

# Add Julian date column (and change "site" class to factor)
df_wt_daily <- df_wt_daily %>% 
  mutate(j_date = yday(date),
         site = factor(site))

# Fit a GLM with Gaussian family
m_glm <- glm(temp ~ j_date + site,
             data = df_wt_daily,
             family = "gaussian")

summary(m_glm) # The model suggests that water temperature increases with Julian
# date, when it is actually a non-linear relationship, due to the limitations of
# the GLM structure. It also shows no difference between wooded and open sites.

# Predict water temperature using the GLM
df_pred <- ggpredict(m_glm,
                     terms = c("j_date [all]",
                               "site [all]")) %>% 
  rename(site = group,
         j_date = x)

# Add prediction to plot
df_wt_daily %>% 
  ggplot(aes(x = j_date,
             y = temp,
             color = site)) +
  geom_point(alpha = 0.25) +
  geom_line(data = df_pred,
            aes(y = predicted)) +
  theme_bw() +
  labs(x = "Julian Date",
       y =  "Water Temperature",
       color = "Wetland Type")
# The GLM prediction is badly wrong. It can't capture the seasonal pattern.

## Generalized additive model with smooth term
m_gam <- gam(temp ~ site + s(j_date), # s() is the smoothing function
             data = df_wt_daily,
             family = "gaussian")

summary(m_gam) # With the better model, site type has a significant effect.
# edf (effective degrees of freedom) indicates the degree of departure from
# a linear function. There is no coefficient for j_date because there is no
# linear representation of it, so no constant coefficient is possible.

# Visualize with GAM prediction
df_pred_gam <- ggpredict(m_gam,
                         terms = c(
                           "j_date [all]", 
                           "site [all]")) %>% 
  rename(site = group,
         j_date = x)

df_wt_daily %>% 
  ggplot(aes(x = j_date,
             y = temp, 
             color = site)) +
  geom_point(alpha = 0.25) +
  geom_line(data = df_pred_gam,
            aes(y = predicted)) +
  theme_bw() +
  labs(x = "Julian Date",
       y = "Water Temperature",
       color = "Wetland Type")

# The prediction now follows the true shape of the curve


# 1. Read directly from the raw GitHub URL
url <- "https://raw.githubusercontent.com/aterui/public-proj_restore-aqua-complex/v.1.0/data_raw/data_bat.csv"

# Try reading normally
df_bat <- read_csv(url, show_col_types = FALSE)

# ============================================================
# DATA GUIDE: Bat Detector Data
# ============================================================

# ----------------------------
# Raw data columns
# ----------------------------

# Site
#   Location where bat detectors are deployed.
#   Levels:
#     "RECCON"  = prairie site without wetland
#     "RECWET"  = prairie site with constructed wetland
#     "WOODCON" = woody site without wetland
#     "WOODWET" = woody site with constructed wetland

# DATE
#   Calendar date of each bat pass record.
#   Expected format: YYYY-MM-DD (verify and standardize).

# TIME
#   Time of bat pass detection.
#   Expected format: HH:MM:SS (verify and standardize).

# AUTO ID*
#   Automatically identified bat species.
#   Species IDs may contain misclassifications or unknown labels
#   that should be carefully reviewed during data cleaning.

# ============================================================
# GOAL 1: Clean data
# ============================================================

# 1. Format column names
#   - Convert column names to a clean format

df_bat <- df_bat %>% 
  janitor::clean_names()

# 2. Examine each column carefully
#   - Check for missing values, inconsistent formats, and typos
unique(df_bat$site) # SITE: No missing values, no unexpected values/typos

unique(df_bat$auto_id) # SPECIES: Some NAs, no apparent typos.

df_bat %>% 
  ggplot() +
  geom_histogram(aes(x = time)) # The time pattern makes sense for bats, and it
# is clear that there has not been confusion between 12-hour and 24-hour formats,
# since no activity is recorded in the middle of the day.

#   - Confirm DATE and TIME are properly parsed as date/time objects
df_bat <- df_bat %>% 
  mutate(date = as.Date(date,
                        format = "%m/%d/%Y"))

df_bat %>% 
  ggplot() +
  geom_histogram(aes(x = date)) # No dates look especially out of place, and the
# overall seasonal pattern makes sense. If there is a format error, it seems
# unlikely to be found.

#   - Inspect AUTO ID values for NA
df_bat %>% 
  filter(is.na(auto_id)) # some NAs, as distinct from NoID.

#   - Remove or correct invalid or unusable records as needed
# Removing NA values but not "No ID"
df_bat <- df_bat %>% 
  drop_na(auto_id)

# New derived columns to create:
# Site-level categories:
#   Prairie sites: "RECCON", "RECWET"
#   Woody sites:   "WOODCON", "WOODWET"

# 3. habitat_type
#   Broad site classification:
#     "prairie" = RECCON, RECWET
#     "woody"   = WOODCON, WOODWET

# 4. wetland_status
#   Presence/absence of wetland:
#     "no_wetland" = RECCON, WOODCON
#     "wetland"    = RECWET, WOODWET

df_bat <- df_bat %>% 
  mutate(habitat_type = if_else(site %in% c("RECCON", "RECWET"),
                                "prairie",
                                "woody"),
         wetland_status = if_else(site %in% c("RECCON", "WOODCON"),
                                  "no_wetland",
                                  "wetland"))
  
# ============================================================
# GOAL 2: Visualize daily bat activity
# ============================================================

# Objective:
#   Quantify and visualize bat activity as the number of bat passes per day.

# Steps:
#   - Aggregate data to calculate daily bat passes
df_bat_daily <- df_bat %>% 
  group_by(date, site) %>% 
  summarize(passes = n()) %>% 
  mutate(habitat_type = if_else(site %in% c("RECCON", "RECWET"),
                                "prairie",
                                "woody"),
         wetland_status = if_else(site %in% c("RECCON", "WOODCON"),
                                  "no_wetland",
                                  "wetland"))

#   - Convert DATE to Julian date
df_bat_daily <- df_bat_daily %>% 
  mutate(j_date = yday(date))

#   - Plot number of bat passes as a function of Julian date
#   - Optionally:
#       * Color or facet plots by site
#       * Smooth trends to visualize seasonal patterns
df_bat_daily %>% 
  ggplot(aes(x = j_date,
             y = passes,
             color = wetland_status)) +
  geom_smooth(alpha = 0.1) +
  geom_point(alpha = 0.5) +
  facet_wrap(facets =~ habitat_type) +
  theme_bw()

# ============================================================
# GOAL 3: Model differences among sites
# ============================================================

# Objective:
#   Test whether bat activity differs among the four detector sites.
#   Does the presence of wetland affect bat activity?
#   Is the effect of wetland is site-dependent?

# Modeling considerations:
#   - Response variable: daily bat passes
#   - Predictors may include:
#       * habitat_type
#       * wetland_status
#       * site (four-level factor)
#       * Julian date (to account for seasonality)
#   - Consider appropriate count models

m_bat <- gam(passes ~ wetland_status * habitat_type + s(j_date),
             data = df_bat_daily,
             family = "nb")

summary(m_bat)