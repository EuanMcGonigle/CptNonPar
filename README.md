
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CptNonPar

<!-- badges: start -->

[![R-CMD-check](https://github.com/EuanMcGonigle/CptNonPar/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/EuanMcGonigle/CptNonPar/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/EuanMcGonigle/CptNonPar/branch/main/graph/badge.svg)](https://app.codecov.io/gh/EuanMcGonigle/CptNonPar?branch=main)
<!-- badges: end -->

Nonparametric change point detection for multivariate time series.
Implements the NP-MOJO methodology proposed in [McGonigle and Cho
(2025)](https://doi.org/doi:10.1093/biomet/asaf024).

## Installation

You can install the released version of `CptNonPar` from
[CRAN](https://CRAN.R-project.org) with:

    install.packages("CptNonPar")

You can install the development version of `CptNonPar` from
[GitHub](https://github.com/) with:

    devtools::install_github("https://github.com/EuanMcGonigle/CptNonPar")

## Usage

For further examples, see the help files within the package. We can
generate an example for change point detection as follows.

We generate a univariate time series of length 1000, with a mean change
at time 300, and an autocovariance (but not marginal) change at time
650. Then, we perform the multi-lag NP-MOJO algorithm with lags 0 and 1,
and print the estimated change points and the associated clusters:

``` r
library(CptNonPar)

n <- 1000
set.seed(123)

noise1 <- stats::arima.sim(model = list(ar = -0.5), n = n, sd = sqrt(1 - 0.5^2))
noise2 <- stats::arima.sim(model = list(ar = 0.5), n = n, sd = sqrt(1 - 0.5^2))

noise <- c(noise1[1:650], noise2[651:n])

signal <- c(rep(0, 300), rep(0.7, 700))

x <- signal + noise

x.c <- np.mojo.multilag(x, G = 166, lags = c(0, 1))

x.c$cpts
#>      cpt lag score
#> [1,] 295   0  1.00
#> [2,] 648   1  0.99

x.c$cpt.clusters
#> [[1]]
#>      cpt lag score
#> [1,] 295   0     1
#> [2,] 296   1     1
#> 
#> [[2]]
#>      cpt lag score
#> [1,] 648   1  0.99
```
