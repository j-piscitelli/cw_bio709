#' DESCRIPTION:
#' Script for piecewise SEM

# in-class ----------------------------------------------------------------
## The Piecewise SEM method relaxes the normality assumption of a regular SEM.
## However, latent variables can't be included.
pacman::p_load(tidyverse,
               GGally,
               piecewiseSEM,
               glmmTMB)

data("keeley")

(df_keeley <- keeley %>% 
    as_tibble())

# psem() requires each constituent dependent variable to have its own model.
# Here we don't change the assumption of normality.
m1 <- lm(abiotic ~ distance, data = df_keeley)
m2 <- lm(hetero ~ distance, data = df_keeley)
m3 <- lm(firesev ~ age, data = df_keeley)
m4 <- lm(cover ~ firesev, data = df_keeley)
m5 <- lm(rich ~ cover + abiotic + hetero, data = df_keeley)

sem_model <- psem(m1, m2, m3, m4, m5)
summary(sem_model)

## Now without that assumption:
m1 <- lm(abiotic ~ distance, data = df_keeley)
m2 <- lm(hetero ~ distance, data = df_keeley)
m3 <- lm(firesev ~ age, data = df_keeley)
m4 <- lm(cover ~ firesev + hetero, data = df_keeley)
# Richness is modeled with a negative binomial distribution
m5 <- MASS::glm.nb(rich ~ cover + abiotic + hetero + distance,
                   data = df_keeley)

sem_model <- psem(m1, m2, m3, m4, m5)
summary(sem_model)
# Fischer's C compares this model against a complete graph and adds -2 * the sum
# of p-values of the links not included in your model (high C means that missing
# links were actually significant). Fischer's C is itself tested to get a p-value:
# this is significant if your model is missing important links.
# The "Tests of directed separation" section shows the paths you didn't include
# and their significance.

plot(sem_model)


## Random effects (they can be included in a piecewise SEM)
data("shipley")

df_shipley <- shipley %>% 
  as_tibble() %>% 
  janitor::clean_names() %>% 
  drop_na(growth)
# Here individual trees have been measured repeatedly over several years, with
# five trees at each of 20 sites. We might want to account for unknown effects of
# each site and each tree.

df_shipley %>% 
  ggpairs(columns = c("dd","date","growth","live")) +
  theme_bw()
# Our model is that dd affects date, date affects growth, and growth affects survival

m1 <- glmmTMB(date ~ dd + (1 | site) +  (1 | tree),
              data = df_shipley,
              family = "gaussian")
m2 <- glmmTMB(growth ~ date + (1 | site) + (1 | tree),
              data = df_shipley,
              family = "gaussian")
m3 <- glmmTMB(live ~ growth + (1 | site) + (1 | tree),
              data = df_shipley,
              family = "binomial")

sem_glmm <- psem(m1, m2, m3)
summary(sem_glmm)
# Note "Marginal" and "Conditional" R-squareds.
# Marginal is with only fixed effects.
# Conditional is both fixed and random effects.
# Here, growth was mostly explained by the random effects (site and tree)

# lab ---------------------------------------------------------------------

library(piecewiseSEM)
data("meadows")

# =========================================
# EXERCISE: Piecewise SEM with Meadows Data
# =========================================
#
# ------------------------------------------------------------
# Dataset: meadows (from piecewiseSEM package)
# Variables:
#   grazed - 0 = ungrazed, 1 = grazed
#   mass   - plant biomass (g/m²)
#   elev   - plot elevation above sea level
#   rich   - plant species richness per m²
# ------------------------------------------------------------
#
# 1. Explore the dataset (structure, summary, plots).
df_meadows <- as_tibble(meadows)
summary(df_meadows)
ggpairs(data = df_meadows)

# 2. Develop a conceptual model: decide which variables influence others.
#    - Consider direct and indirect effects.
#    - Think about grazing as a disturbance factor.

# There are correlations between all pairs of variables. Richness is likely to be
# directly affected by elevation and grazing, but higher altitudes are apparently
# more grazed. I looked up the paper this data set comes from and it seems that
# ungrazed meadows were selected to be as similar as possible to grazed meadows,
# so it's strange that it should be related to elevation (grazing doesn't change
# elevation, I assume!). I expect that grazing reduces species richness directly
# by removing the species preferred by the grazers (possibly this would increase
# richness as a second-order effect, but the fact that the overall trend is negative
# seems to suggest that that isn't happening). Grazing presumably decreases plant
# biomass directly, too. The negative relationship between biomass and richness
# is odd. If the explanation is that, say, large bunchgrasses are shading out
# competition, and grazing removes them but still manages to lower richness,
# and elevation primarily affects which species can grow in a given meadow, I
# would expect a graph like this:

# 3. Fit component models (e.g., lm) for each hypothesized relationship.
m1 <- MASS::glm.nb(rich ~ elev + grazed + mass,
         data = df_meadows)
m2 <- glm(mass ~ grazed,
         data = df_meadows)

# 4. Combine models into a piecewise SEM using psem().
sem_meadows <- psem(m1, m2,
     data = df_meadows)

plot(sem_meadows)

summary(sem_meadows)

# 5. Evaluate the SEM: path coefficients, significance, variance explained.
# I haven't come up with a particularly good model. My Fischer's C p-values have
# all been near 0, and the links I left out were always significant.
# The model explains only a small amount of the total variance.

# 6. Optional: try alternative models if your model deviates from the expectation.
#
# Deliverables:
# - Code for component models and combined SEM
# - Conceptual SEM diagram
# - Short reasoning about your SEM results
