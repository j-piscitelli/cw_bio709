# NOTE:
# When instructed to "test XXX", you must report the outcome as comments
# that clearly summarize the relevant statistical results (e.g., effect size,
# direction, significance, and interpretation).
# Providing code alone without documenting and interpreting the results
# in comments will result in point deductions.

pacman::p_load(tidyverse,
               mgcv,
               GGally,
               vegan,
               lavaan,
               lavaanPlot,
               glarma,
               piecewiseSEM,
               lterdatasampler,
               forecast,
               daymetr)


# dataset 1 ---------------------------------------------------------------

link1 <- "https://raw.githubusercontent.com/aterui/biostats/master/data_raw/data_insect_emergence.rds"
df_emg <- readRDS(url(link1, "rb"))

# This dataset ('df_emg') contains daily measurements of aquatic insect emergence
# from two wetland sites over a full calendar year (Jan 1–Dec 31).

# Data structure:
# t           : Day of the year (integer), where 1 = January 1 and 365 = December 31
# site        : Site identifier (factor), with "s1" and "s2" representing the two wetlands
# emergence   : Emergence flux of aquatic insects (g/day)

# Q1. Visualize seasonal patterns in emergence flux at both sites
#     (e.g., plot emergence vs. day of year, with separate lines or colors for each site).
#     [1 point]

ggplot() +
  geom_point(data = df_emg,
             aes(x = t,
                 y = emergence,
                 color = site)) +
  theme_bw()

# Q2. Test whether emergence flux differs significantly between the two sites,
#     while appropriately accounting for seasonal variation
#     [4 points]

m_gam <- gam(emergence ~ site + s(t),
             data = df_emg,
             family = "gaussian")

summary(m_gam)

## After accounting for seasonal variation with the smoothing function, there
# is a significant difference between the two sites. The effect of seasonal 
# variation is, unsurprisingly, much stronger than the difference between sites.

## I notice looking at the graph that the main difference between sites appears
# to be that Site 1 is slightly ahead of Site 2 phenologically, rather than that
# one site is overall higher than the other (though Site 1 is a bit higher in the
# summer emergence period). I am not sure how to tell this information from the
# GAM output, though.

# dataset 2 ---------------------------------------------------------------

link2 <- "https://raw.githubusercontent.com/aterui/cw_bio709/master/data_fmt/data_lake_invert.rds"
df_inv <- readRDS(url(link2, "rb"))

# This dataset 'df_inv' contains 100 observations from 10 lakes.
# Within each lake, 10 plots were established, spaced ~500 m apart.
# At each plot, the following variables were measured:

# s          : Species richness of invertebrates associated with aquatic plants at each plot
# hb         : Standing biomass of invertebrates associated with aquatic plants at each plot
# prod       : Production rate of aquatic plants (macrophytes), measured as g/month
# substrate  : Median diameter of substrate materials (mm)
# cond       : Water electrical conductivity (µS/cm);
#              a proxy for ionized nutrient levels (higher values may indicate eutrophication)
# lake       : lake ID

# Researcher's hypothesis was that: 
# (a) conductivity influences the productivity of macrophyes.
# (b) macrophyte's production rate ('prod') dictates invertebrate biomass ('hb') through bottom-up effects
# (c) macrophyte's production rate ('prod') dictates invertebrate richness ('s') through bottom-up effects 

# Q1. Create a scatter plot of macrophyte production ('prod', y-axis)
#     versus water conductivity ('cond', x-axis), with points colored by lake identity.
#     [1 point]

ggplot() +
  geom_point(data = df_inv,
             aes(x = cond,
                 y = prod,
                 color = lake))

# Q2. Create a scatter plot of raw invertebrate biomass ('hb', y-axis)
#     versus macrophyte production ('prod', x-axis), with points colored by lake identity.
#     [1 point]

ggplot() +
  geom_point(data = df_inv,
             aes(x = prod,
                 y = hb,
                 color = lake))

# Q3. Create a scatter plot of "log-transformed" invertebrate biomass ('hb', y-axis)
#     versus macrophyte production ('prod', x-axis), with points colored by lake identity.
#     [1 point]

ggplot() +
  geom_point(data = df_inv,
             aes(x = prod,
                 y = log(hb),
                 color = lake))

# Q4. Test hypothesis (a) by modeling macrophyte production while
#     statistically controlling for potential confounding variables ('substrate', 'lake').
#     [3 points]

df_inv <- df_inv %>% 
  mutate(log_hb = log(hb),
         log_prod = log(prod),
         log_s = log(s))

m_a <- glmmTMB::glmmTMB(prod ~ cond + substrate + (1 | lake),
                        data = df_inv,
                        family = "gaussian")

summary(m_a)

## This GLMM model finds a statistically significant effect of conductivity on
# productivity, after controlling for substrate size and lake identity.

# Q5. Test hypotheses (a–c) simultaneously using a unified modeling framework.
#     Based on the resulting statistical tests, determine whether the overarching
#     hypothesis (a–c, combined) is supported or rejected.
#     - Use appropriate probability distributions.
#     - Use variable transformation if appropriate given the data.
#     [4 points]

ggpairs(df_inv %>% 
          select(c(s, hb, prod, cond)))

m_b <- glmmTMB::glmmTMB(hb ~ prod + substrate + (1 | lake),
               data = df_inv,
               family = "gaussian")

m_c <- glmmTMB::glmmTMB(s ~ prod + substrate + (1 | lake),
               data = df_inv,
               family = "poisson")

psem_abc <- psem(m_a, m_b, m_c)

summary(psem_abc)

## The three hypotheses together are supported by the data. Not only are the
# component hypotheses each individually significant (with positive relationships
# between cond and prod, prod and s, and prod and hb), but the alternative causal
# pathways (such as s ~ cond) were not significant, suggesting that the proposed
# causal structure is largely accurate. There is a near-significant connection
# between s and hb, suggesting that this causal pathway may also play a role.

## I do not see a strong reason to transform the variables. I tried log(hb) as
# a response variable and got qualitatively similar results. Although log(hb) 
# produces slightly stronger correlations, I don't see a theoretical reason to
# prefer the transformed version in the context of the effect of macrophyte
# productivity--I would expect the relationship to be more or less linear. The
# difference in performance between log(hb) and hb is not great enough for me
# to feel confident that the (seemingly, to me) less theoretically-supported
# measure is really better, rather than that it did better here by chance.

# dataset 3 ---------------------------------------------------------------

link3 <- "https://raw.githubusercontent.com/aterui/cw_bio709/master/data_fmt/nutrient.rds"
nutrient <- readRDS(url(link3, "rb"))

print(trees)

# This dataset ('trees') contains measurements of 31 felled black cherry trees.
# The three variables represent tree diameter, height, and timber volume.
# Note: the variable 'Girth' is actually the diameter measured at 4 ft 6 in above ground.

# Data structure:
# Girth   : Numeric, tree diameter in inches (mislabelled as girth)
# Height  : Numeric, tree height in feet
# Volume  : Numeric, timber volume in cubic feet

# Q1. Visualize relationships among tree diameter ('Girth'), height ('Height'),
#     and timber volume ('Volume') (e.g., using scatterplot matrix or pairwise scatter plots).
#     [1 point]

ggpairs(trees)
## There are apparent positive relationships between all three variables.

# Q2. Perform an appropriate ordination or dimension reduction method to 
#     summarize these three variables into fewer composite axes.
#     Then, identify and retain axes that explain meaningful variation in the original variables
#     [3 points]

obj_pca_trees <- prcomp(x = trees,
                        center = TRUE,
                        scale = TRUE)

summary(obj_pca_trees)

## Most variance (80%) is explained by PC1, which is positively related to all
# three variables.

df_trees_pca <- bind_cols(trees,
                          as_tibble(obj_pca_trees$x)) %>% 
  select(- c(PC2, PC3))

# Q3. If justified, test whether the retained axis (or axes) is significantly 
#     related to "nutrient"; 
#     skip regression if the ordination does not support meaningful interpretation.
#     [1 point]

df_trees_pca <- bind_cols(df_trees_pca,
                          nutrient = nutrient)

lm(PC1 ~ nutrient,
   data = df_trees_pca) %>% 
  summary()

## There is a strong and significant positive correlation of nutrient with PC1.
# Given that our three variables all represent aspects of tree size, it is not
# surprising that PC1, which was positively related to all three variables, is
# also positively related to nutrient availability.

# dataset 4 ---------------------------------------------------------------

df_nile <- dplyr::tibble(
  year = time(Nile), # observation year
  discharge = as.numeric(Nile) # discharge
)

df_sunspot <- dplyr::tibble(
  year = time(sunspot.year), # observation year
  sunspots = as.numeric(sunspot.year) # the number of sunspots
)

# These datasets contain:
# - df_nile    : Annual discharge of the Nile River (Nile dataset)
# - df_sunspot : Annual sunspot counts (sunspot.year dataset)

# Q1. Create a combined data frame aligning the observation years
#     (i.e., only include years present in both datasets)
#     [1 point]

df_ns <- df_sunspot %>% 
  filter(year %in% df_nile$year) %>% 
  left_join(df_nile %>% 
              filter(year %in% df_sunspot$year),
            by = "year") %>% 
  arrange(year)

# Q2. Test whether the number of sunspots is significantly related to Nile's discharge
#     [4 points]

obj_ns_arima <- auto.arima(df_ns$discharge,
                           xreg = df_ns$sunspots,
                           stepwise = FALSE)

summary(obj_ns_arima)
confint(obj_ns_arima, level = .95)
confint(obj_ns_arima, level = .5)

## The confidence interval for xreg (i.e. number of sunspots) in this ARIMAX
# analysis includes zero. In fact, it still includes 0 (barely) if we lower our
# threshold to 50%! This suggests that it is very unlikely that there is a true
# relationship between sunspot count and Nile discharge.
