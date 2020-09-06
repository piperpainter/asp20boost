
# asp20boost

Welcome to *asp20boost*. We implement both simple boosting and
component-wise boosting ensemble algorithms for location and scale
settings.

  - Supervisor: Thomas Kneib
  - Students: Johannes, Levin, Sebastián

We have successfully implemented both boosting variants. Currently, we
work on the publication of the package i.e. we perform unit tests and
expand the documentation.

### Literature

  - Section **4.3** on “Boosting Linear Regression Models” in
    “Regression: Models, methods, and applications” by Ludwig
    Fahrmeir, Thomas Kneib, Stefan Lang, and Brian Marx
  - Section **2.9.1** on “Regression models for location, scale, and
    shape” in ibd.

<!--
more datasets:
- gamlss.data::TODO
- MASS::Boston
- https://archive.ics.uci.edu/ml/datasets.php
-->

This is a visual representation of a location scale Regression Model
(taken from asp20model)

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

To install and work with out package you need first to install the
`asp20model` & `asp20boost` package from GitLab and load it. The
installing of devtools requires attention.

``` r
## #install.packages("devtools")
## #devtools::install_gitlab("asp20/asp20model", host = "gitlab.gwdg.de")
#devtools::install_gitlab("asp20/asp20boost", host = "gitlab.gwdg.de")
```
