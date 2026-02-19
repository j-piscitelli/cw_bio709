#' DESCRIPTION:
#' Script for Constrained Ordination

# in-class ----------------------------------------------------------------
pacman::p_load(tidyverse,
               GGally,
               vegan)

data("varespec", "varechem")

m_y <- varespec

colnames(m_y) <- str_to_lower(colnames(varespec))

df_env <- as_tibble(varechem) %>% 
  janitor::clean_names()

m_y %>% 
  ggpairs(columns = 1:3,
          aes(alpha = 0.5)) +
  theme_bw()

## RDA (Redundancy Analysis)
# Constrained version of PCA; thus, assumes linearity and normality
(obj_rda <- rda(m_y ~ n + p + ca,
                data = df_env))
# The proportion of variance explained by the predictors is 0.26.
# The first two RDA axes explain most of the constrained variance;
# however, most variance is from the unconstrained portion. 

# Permutation test on the RDA object
# by = "margin" tests the marginal effect of each predictor, after accounting
# for other predictors (a Type II ANOVA). The default in anova.cca() is Type I.
anova.cca(obj_rda,
          by = "margin",
          permutations = 999)

# Visualization
df_rda <- scores(obj_rda,
                 display = "site", # show values at each point
                 scaling = 2) %>% 
  bind_cols(df_env) %>% 
  janitor::clean_names()

df_bp <- scores(obj_rda,
                display = "bp", # show relationship to constraining variables
                scaling = 2) %>% 
  as_tibble(rownames = "variable") %>% 
  janitor::clean_names()

df_rda %>% 
  ggplot(aes(x = rda1,
             y = rda2)) +
  geom_point(aes(color = n)) +
  geom_segment(data = df_bp,
               aes(x = 0, xend = rda1 * 10,
                   y = 0, yend = rda2 * 10),
               arrow = arrow(length = unit(0.2, "cm"))) +
  geom_text(data = df_bp,
            aes(x = rda1 * 10.5,
                y = rda2 * 10.5,
                label = variable),
            size = 4) +
  theme_bw() +
  labs(x = "RDA1",
       y = "RDA2",
       color = "Nitrogen") +
  scale_color_viridis_c()

# The "horseshoe" distribution of points is a sign that there is not a linear
# relationship between the variables. The RDA results won't be reliable.

## dbRDA (distance-based RDA)
# Useful in the absence of linear correlation or normal distribution
(obj_db <- dbrda(m_y ~ n + p + ca,
                 data = df_env,
                 distance = "bray")) # Bray-Curtis distance
# Summary results appear more or less similar to the RDA, so far...

anova.cca(obj_db,
          by = "margin",
          permutations = 999)
# ... however, we no longer see the seemingly-significant result for N; it was an
# artifact of the RDA's linearity assumption applied to inappropriate data.

# Visualization
df_db <- scores(obj_db, 
                display = "sites",
                scaling = 2) %>% 
  as_tibble() %>%              
  bind_cols(df_env) %>%        
  janitor::clean_names()       

df_bp <- scores(obj_db, 
                display = "bp", 
                scaling = 2) %>% 
  as_tibble(rownames = "variable") %>% 
  janitor::clean_names()

df_db %>% 
  ggplot(aes(x = db_rda1,
             y = db_rda2)) +        # color sites by nitrogen level
  geom_point(aes(color = n)) +
  geom_segment(data = df_bp,
               aes(x = 0, xend = db_rda1,
                   y = 0, yend = db_rda2),
               arrow = arrow(length = unit(0.2, "cm"))
  ) +
  geom_text(data = df_bp,
            aes(x = db_rda1 * 1.1,    # slightly beyond arrow tip
                y = db_rda2 * 1.1,
                label = variable),  # or use a variable column
            size = 4) +
  theme_bw() +
  labs(x = "dbRDA1",
       y = "dbRDA2",
       color = "Nitrogen") +
  scale_color_viridis_c()

# lab ---------------------------------------------------------------------

# ============================================================
# EXERCISE: Community Ordination and Environmental Gradients
# ============================================================

library(vegan)
data("mite", "mite.env")

# The mite datasets contain information on Oribatid mite communities
# sampled from a small peatland area (2.5 m × 10 m).
#
# There are linked datasets:
# ------------------------------------------------------------
# mite     : Species abundance data (35 mite species × 70 sites)
# mite.env : Environmental variables measured at the same sites
# ------------------------------------------------------------
#
# Environmental variable descriptions (mite.env):
# ------------------------------------------------------------
# SubsDens : Substrate density (g/L)
# WatrCont : Water content of the substrate (g/L)
# Substrate: Substrate type (factor with multiple levels)
# Shrub    : Shrub density (ordered factor: low → high)
# Topo     : Microtopography (Blanket vs Hummock)
# ------------------------------------------------------------

mite <- mite %>% 
  janitor::clean_names() %>% 
  vegan::wisconsin() # Wisconsin standardization

mite.env <- mite.env %>% 
  janitor::clean_names()

# 1. Explore and visualize interrelationships among species abundances.
#    - Examine patterns of co-occurrence.
#    - Assess whether relationships among species appear linear or nonlinear.
mite %>% 
  ggpairs(columns = 1:10,
          aes(alpha = 0.5)) +
  theme_bw()
# A few pairs are show linear relationships, but most do not.


# 2. Fit a redundancy analysis (RDA) model using environmental variables of your choice.
#    - Visualize the ordination results.
#    - Examine gradients and species–environment relationships.
#    - Evaluate whether the assumptions of RDA are appropriate for these data.

# The mite data includes many 0-abundances, and the relationships between species
# are not generally linear. RDA is a bad fit for this situation.

obj_mrda <- rda(mite ~ subs_dens + watr_cont,
                data = mite.env)

df_mrda <- scores(obj_mrda,
                  display = "site",
                  scaling = 2) %>% 
  bind_cols(mite.env) %>% 
  janitor::clean_names()

df_mbp <- scores(obj_mrda,
                 display = "bp",
                 scaling = 2) %>% 
  as_tibble(rownames = "variable") %>% 
  janitor::clean_names()

df_mrda %>% 
  ggplot(aes(x = rda1,
             y = rda2)) +
  geom_point() +
  geom_segment(data = df_mbp,
               aes(x = 0, xend = rda1,
                   y = 0, yend = rda2),
               arrow = arrow(length = unit(0.2, "cm"))) +
  geom_text(data = df_mbp,
            aes(x = rda1 * 1.05,
                y = rda2 * 1.05,
                label = variable),
            size = 4) +
  theme_bw()

# 3. Apply alternative ordination methods.
#    - Canonical correspondence analysis (CCA; see ?cca()).
obj_mcca <- cca(mite ~ subs_dens + watr_cont,
                data = mite.env)


df_mcca <- scores(obj_mcca,
                  display = "site",
                  scaling = 2) %>% 
  bind_cols(mite.env) %>% 
  janitor::clean_names()

df_mbpc <- scores(obj_mcca,
                  display = "bp",
                  scaling = 2) %>% 
  as_tibble(rownames = "variable") %>% 
  janitor::clean_names()

df_mcca %>% 
  ggplot(aes(x = cca1,
             y = cca2)) +
  geom_point() +
  geom_segment(data = df_mbpc,
               aes(x = 0, xend = cca1 * 5,
                   y = 0, yend = cca2 * 5),
               arrow = arrow(length = unit(0.2, "cm"))) +
  geom_text(data = df_mbpc,
            aes(x = cca1 * 5.25,
                y = cca2 * 5.25,
                label = variable),
            size = 4) +
  theme_bw()

#    - Distance-based RDA (dbRDA).
obj_mdb <- dbrda(mite ~ subs_dens + watr_cont,
                data = mite.env,
                distance = "bray")

df_mdb <- scores(obj_mdb,
                  display = "site",
                  scaling = 2) %>% 
  bind_cols(mite.env) %>% 
  janitor::clean_names()

df_mbpd <- scores(obj_mdb,
                  display = "bp",
                  scaling = 2) %>% 
  as_tibble(rownames = "variable") %>% 
  janitor::clean_names()

df_mdb %>% 
  ggplot(aes(x = db_rda1,
             y = db_rda2)) +
  geom_point() +
  geom_segment(data = df_mbpd,
               aes(x = 0, xend = db_rda1,
                   y = 0, yend = db_rda2),
               arrow = arrow(length = unit(0.2, "cm"))) +
  geom_text(data = df_mbpd,
            aes(x = db_rda1 * 1.05,
                y = db_rda2 * 1.05,
                label = variable),
            size = 4) +
  theme_bw()

# 4. Compare RDA, CCA, and dbRDA.
#    - Perform permutation analysis to examine the significance of predictor variables
#    - Discuss which method is most appropriate for these data and why.

anova.cca(obj_mrda,
          by = "margin",
          permutations = 999)

anova.cca(obj_mcca,
          by = "margin",
          permutations = 999)

anova.cca(obj_mdb,
          by = "margin",
          permutations = 999)

# After Wisconsin standardization, all three methods find a significant effect
# for both substrate density and water content. Before standardization, the RDA
# did not show a significant result for substrate density. The lack of consistent
# linear relationships between species makes RDA unreliable with this data; CCA
# or dbRDA would be preferable.