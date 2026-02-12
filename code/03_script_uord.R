#' DESCRIPTION:
#' Script for Unconstrained Ordination

# in-class ----------------------------------------------------------------

pacman::p_load(tidyverse,
               GGally,
               vegan)

# PCA ---------------------------------------------------------------------


df_iris <- iris %>% 
  as_tibble() %>% 
  janitor::clean_names()

# Check correlations (are variables linearly related?)
df_iris %>%
  ggpairs(
    progress = FALSE,
    columns = c("sepal_length",
                "sepal_width",
                "petal_length",
                "petal_width"),
    aes(color = species,
        alpha = 0.5)
  ) +
  theme_bw()

## PCA using petal traits
# Reduce to columns that we want to use
df_petal <- df_iris %>% 
  select(starts_with("petal_"))

# PCA using prcomp()
obj_pca <- prcomp(
  x = df_petal,
  center = TRUE, # center at 0 before analysis
  scale = TRUE # scale by standard deviation
)

summary(obj_pca)
#  - Standard deviation shows how much the data varies on that axis
#  - Proportion of Variance shows how much information is in that axis
#  - Cumulative proportion shows the amount of information found cumulatively
# in that axis and those more important than it. We want at least 70% explained
# in the group axes we'll keep

# Add PCA coordinates to original data
# obj_pca$x shows the coordinates of each row along each PC axis
df_pca <- df_iris %>% 
  bind_cols(obj_pca$x)

# Here PC1 is very dominant; we can use it to summarize both petal variables as one
ggplot(df_pca,
       aes(x = species,
           y = PC1,
           fill = species)) +
  geom_boxplot() +
  theme_bw()


# NMDS --------------------------------------------------------------------

# NMDS doesn't require normality in the original data, but retains less information

# Visualize a sample of dune plant community composition data
data(dune)

dune %>% 
  as_tibble() %>% 
  select(1:3) %>% 
  ggpairs() +
  theme_bw()

# Compute pairwise dissimilarities between sites using Bray-Curtis distance
m_bray <- vegdist(dune,
                  method = "bray")

# Run NMDS
obj_nmds <- metaMDS(comm = m_bray,
                    k=2) # number of dimensions to reduce to

# Load environmental data at sites
data(dune.env)

df_nmds <- dune.env %>% 
  as_tibble() %>% 
  bind_cols(obj_nmds$points) %>% 
  janitor::clean_names()

ggplot(df_nmds,
       aes(x = mds1,
           y = mds2,
           color = use)) +
  geom_point(size = 3) +
  stat_ellipse(level = .95,
               linetype = 2) +
  theme_bw() +
  labs(color = "Land-use intensity",
       x = "NMDS1",
       y = "NMDS2")
# Ellipses overlap considerably--land use probably doesn't matter much

# Test significance with PERMANOVA
adonis2(m_bray ~ use,
        data = df_nmds)

# lab ---------------------------------------------------------------------

# ============================================================
# EXERCISE: PCA using the iris dataset
# ============================================================

# In this exercise, you will perform a Principal Component
# Analysis (PCA) using all morphological measurements in the
# iris dataset and visualize multivariate trait patterns
# among species.

# 1. Using all four morphological variables
#    (Sepal.Length, Sepal.Width, Petal.Length, Petal.Width),
#    perform a PCA.
df_morph <- df_iris %>% 
  select(-species)

pca_morph <- prcomp(x = df_morph,
                    center = TRUE,
                    scale = TRUE)

# 2. Visualize flower morphology in PC axes whose cumulative
#    contribution exceeds 90%; color points by species.

summary(pca_morph) # First two components cover ~96% of variance

df_pca_morph <- df_iris %>% 
  bind_cols(pca_morph$x)

ggplot(df_pca_morph,
       aes(x = PC1,
           y = PC2,
           color = species)) +
  geom_point() +
  theme_bw()

## Species differ from each other along PC1, and have considerable within-species
## variation along PC2.

# 3. Which morphological traits contribute most strongly to
#    the first and second principal components? How?

## PC1 is roughly equal parts sepal_length, petal_length, and petal_width, with
## a smaller negative contribution from sepal_width.

## PC2 is primarily sepal_width, with a smaller contribution from sepal_length.
## Both contributions are negative, so a higher PC2 value ~means a smaller sepal.


# ============================================================
# EXERCISE: NMDS using the BCI dataset
# ============================================================

# In this exercise, you will perform a Non-metric Multidimensional
# Scaling (NMDS) using the BCI tree community dataset and explore
# patterns in species composition among sites.

data("BCI", "BCI.env")

# 1. Using the BCI dataset, calculate a dissimilarity matrix
#    (e.g., Bray-Curtis) and perform NMDS.

m_BCI_bray <- vegdist(BCI,
                      method = "bray")

nmds_BCI <- metaMDS(m_BCI_bray,
                    k =2)

# 2. Visualize sites in NMDS space.
#    - How are sites positioned relative to each other?
#    - Color or shape points by environmental groups or site
#      characteristics of your choice.

df_nmds_BCI <- BCI.env %>% 
  as_tibble() %>% 
  bind_cols(nmds_BCI$points)

ggplot(df_nmds_BCI,
       aes(x = MDS1,
           y = MDS2,
           color = Habitat)) +
  geom_point() +
  stat_ellipse(level = .95,
               linetype = 2) +
  theme_bw()
## The different habitat types appear fairly distinct in community composition,
# though "OldLow" overlaps the others somewhat.

# 3. Perform PERMANOVA to examine if communities are grouped
#    by the environmental variable you selected.
adonis2(m_BCI_bray ~ Habitat,
        data = df_nmds_BCI)
## The PERMANOVA confirms that there is a significant difference in community
# composition between habitats
