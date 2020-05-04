
<!-- README.md is generated from README.Rmd. Please edit that file -->

# asp20boost

  - Supervisor: Thomas Kneib
  - Students: Johannes, Levin, Sebastián
  - Tasks: Parameter estimation and variable selection using the
    boosting ensemble learning algorithm

## The Model

``` r
n <- 500
x <- runif(n)
y <- x + rnorm(n, sd = exp(-3 + 2 * x))
plot(x, y)
abline(0, 1, lwd = 2)
curve(x + 1.96 * exp(-3 + 2 * x), -0.1, 1.1, add = TRUE)
curve(x - 1.96 * exp(-3 + 2 * x), -0.1, 1.1, add = TRUE)
```

<img src="man/figures/README-data-1.png" width="100%" />

## Literature

  - Section **4.3** in “Regression: Models, methods, and applications”
    by Ludwig Fahrmeir, Thomas Kneib, Stefan Lang, and Brian Marx,
    especially Section 2.9.1 on “Regression models for location, scale,
    and shape”

<!--
more datasets:
- gamlss.data::TODO
- MASS::Boston
- https://archive.ics.uci.edu/ml/datasets.php
-->
