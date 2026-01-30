#' DESCRIPTION:
#' Script for introductory work

library(tidyverse)

# in-class ----------------------------------------------------------------
# Create fictional data ---------------------------------------------------
# Tibble of imaginary data, group x generally higher than group y
x <- c(3.2, 5.0, 10.0, 100, 50)
y <- c(1.0, 3.0, 2.0, 2.1, 1.2)

df_xy <- tibble(group = c(rep("x", length (x)),
                          rep("y", length(y))),
                value = c(x, y))

# Non-parametric versions of common parametric tests ----------------------
# Non-parametric tests
## Non-parametric tests generally convert data to rank

# t-test (parametric)
t.test(x, y)

# non-parametric version of t-test
## U-test/Wilcoxon test
wilcox.test(x, y)
# calculates W


# ANOVA (parametric) - for more than 2 groups
aov(weight ~ group,
    data = PlantGrowth)

# non-parametric version of ANOVA: Kruskal-Wallis test
kruskal.test(weight ~ group,
             data = PlantGrowth)



# Confidence intervals -----------------------------------------------------
# Confidence intervals are another way of seeing p-values
# If the 95% confidence interval doesn't overlap 0, generally p < .05

m <- lm(Petal.Length ~ Petal.Width,
        data = iris)

summary(m)

# confint() function gives bounds of 95% confidence interval
confint(m)


# Covariance and correlation --------------------------------------------------------------
# Covariance between x and y is based on how much the a data point's deviation
# from the mean of x is related to its deviation from the mean of y. The average
# product of these values (technically, the summed product over *one less than*
# the number of data points) is the covariance.
# Covariance is also the correlation coefficient (Pearson's rho) times the stds
# of x and y.
# Correlation is thus covariance standardized from -1 to 1.

x <- rnorm(100, mean = 0, sd = 1)
y <- rnorm(100, mean = 0.8 * x, sd = 1)

plot(x ~ y)

# Parametric test of correlation: cor.test() by default uses Pearson's rho,
# assuming a normal distribution
cor.test(x, y)

# Non-parametric version: cor.test() with method = "spearman
cor.test(x, y,
         method = "spearman")

# Covariance calculation: cov()
cov(x, y)

# Dividing by product of standard deviations gets Pearson correlation
cov(x, y)/(sd(x) * sd(y))
