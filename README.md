
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CptNonPar

<!-- badges: start -->

[![R-CMD-check](https://github.com/EuanMcGonigle/CptNonPar/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/EuanMcGonigle/CptNonPar/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/EuanMcGonigle/CptNonPar/branch/main/graph/badge.svg)](https://app.codecov.io/gh/EuanMcGonigle/CptNonPar?branch=main)
<!-- badges: end -->

Nonparametric change point detection for multivariate time series.
Implements the NP-MOJO methodology proposed in

> McGonigle, E. T., Cho, H. (2023). Nonparametric data segmentation in
> multivariate time series via joint characteristic functions. arXiv
> preprint [arXiv:2305.07581](https://arxiv.org/abs/2305.07581).

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
#>       cp lag p.val
#> [1,] 295   0 0.000
#> [2,] 648   1 0.005

x.c$cpt.clusters
#> [[1]]
#>       cp lag p.val
#> [1,] 295   0     0
#> [2,] 296   1     0
#> 
#> [[2]]
#>       cp lag p.val
#> [1,] 648   1 0.005
```
