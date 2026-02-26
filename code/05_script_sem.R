#' DESCRIPTION:
#' Script for SEM
pacman::p_load(tidyverse,
               GGally,
               vegan,
               lavaan,
               lavaanPlot)

# in-class ----------------------------------------------------------------
url <- "https://raw.githubusercontent.com/aterui/biostats/master/data_raw/data_foodweb.csv"

(df_fw <- read_csv(url))

### Path analysis

# We can model the causal structure of our hypothesis with a directed graph. If
# it is acyclic, we can model it with Path Analysis.

# Visualize the relationships between variables to help consider which ones are
# reasonable to include in the directed acyclic graph, taking into account
# plausibility.
df_fw %>% 
  select(-plot_id) %>% 
  ggpairs() +
  theme_bw()

# Formula for a SEM model (first row is first causal layer, second row is second)
m1 <- '
  mass_herbiv ~ mass_plant + cv_h_plant
  mass_pred ~ mass_herbiv
  '

fit1 <- sem(model = m1,
            data = df_fw)
# The p-value here has an unusual interpretation: it is the significance of the
# difference between the data and the model's prediction of the data. Therefore,
# a significant p-value means that the model is inaccurate.

# Degrees of freedom is the number of arrows dropped from the full possible DAG
# in your hypothesized model. The full graph won't include arrows from one
# non-response variable to another, however.

summary(fit1,
        standardize = TRUE) # See the Std.all column for mutually comparable figures

# Plot with strengths of influences
lavaanPlot(model = fit1,
           coefs = TRUE,
           stand = TRUE)

## Model comparison
# We could make a more complex model with an additional causal pathway.
m2 <- '
mass_herbiv ~ mass_plant + cv_h_plant
mass_pred ~ mass_herbiv + cv_h_plant
'
fit2 <- sem(model = m2,
            data = df_fw)

# We can compare our two models with an ANOVA
anova(fit1, fit2)
# The first model has a lower AIC, indicating greater predictive value. Also,
# the  p-value shows that the two models aren't significantly different. Thus,
# the additional arrow isn't adding explanatory power, and so isn't likely to
# represent a real direct causal pathway.


### Structural Equation Modeling
# SEM seeks "latent variables" that are unmeasured. Path analysis is SEM without
# latent variables.

url <- "https://raw.githubusercontent.com/aterui/biostats/master/data_raw/data_herbivory.csv"

(df_herbv <- read_csv(url))

# Visualize relationships between variables:
# We might want sla, cn_ratio, and per_lignin to be separate intermediate causal
# steps between soil N and herbivory, but if they are correlated with each other
# then there may really only be one causal pathway through a latent variable that
# influences all three.
# (SEM requires _linear_ correlation, like a PCA.)
df_herbv %>% 
  ggpairs(columns = c("soil_n",
                      "sla",
                      "cn_ratio",
                      "per_lignin")) +
  theme_bw()

# To specify that our hypothesis includes a latent variable, we use =~ in sem()
m_sem <- '
# latent variable
palatability =~ sla + cn_ratio + per_lignin

# regressions
palatability ~ soil_n
herbivory ~ palatability
'

fit_sem <- sem(m_sem,
               data = df_herbv)

summary(fit_sem)

# lab ---------------------------------------------------------------------

# ============================================================
# EXERCISE: Path Analysis and Covariance Visualization
# ============================================================

library(piecewiseSEM)
data("keeley")

# The "keeley" dataset contains fire-related vegetation data
# collected from shrublands in California.
#
# ------------------------------------------------------------
# Column descriptions:
# elev  : Elevation of the site
# abiotic : Overall suitability of abiotic conditions for plants
# hetero: Environmental heterogeneity
# heat  : Heat load index (a function of slope and aspect)
# firesev: Fire severity
# age   : Time since last fire
# cover : Vegetation cover
# rich  : Plant species richness
# ------------------------------------------------------------
#
# In this exercise, you will explore relationships among variables
# using covariance and path analysis. You will replicate a published
# path model and propose an alternative.

# 1. For the variables depicted in Figure 22.1, draw a figure
#    showing the covariance between variables.

ggpairs(keeley,
        columns = c("rich",
                    "abiotic",
                    "hetero",
                    "cover",
                    "firesev",
                    "age",
                    "distance"))

# 2. Following Figure 22.1, develop a path model using the
#    same variables and relationships. Examine if this model
#    captures the data structure using a Chi-Square test.

m_k <- '
rich ~ abiotic + hetero + cover
cover ~ firesev
firesev ~ age
age ~ distance
hetero ~ distance
abiotic ~ distance
'

fit_k <- sem(m_k,
             data = keeley)

summary(fit_k)
# The p-value derived from the Chi-square test is very low, indicating that the
# model prediction is significantly different from the data. Not so good!


# 3. Develop an alternative path model that you consider more
#    appropriate based on theory or observed data patterns.
m_k2 <- '
rich ~ hetero + cover + abiotic + distance
cover ~ firesev + hetero
hetero ~ distance
firesev ~ age
abiotic ~ distance
'

fit_k2 <- sem(m_k2,
              data = keeley)

summary(fit_k)
# 4. Compare the performance of the published model (Figure 22.1)
#    and your alternative model.
#    - Consider fit indices, path coefficients, and interpretability.

# I haven't found anything much better than the original. Many things probably
# relate to each other--for example, poor abiotic conditions for plants might
# increase age at fire by making the area harder to burn.