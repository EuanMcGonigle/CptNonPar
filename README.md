<!-- badges: start -->
[![R-CMD-check](https://github.com/EuanMcGonigle/CptNonPar/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/EuanMcGonigle/CptNonPar/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->
# CptNonPar
Nonparametric change point detection for multivariate time series. Implements the NP-MOJO methodology proposed in

> McGonigle, E. T., Cho, H. (2023). Nonparametric data segmentation in multivariate time series via joint characteristic functions. 
> [(link to paper here)](https://arxiv.org/abs/2305.07581).

## Installation

To install `CptNonPar` from CRAN:

```
install.packages("CptNonPar")
```


To install `CptNonPar` from GitHub:

```
devtools::install_github("https://github.com/EuanMcGonigle/CptNonPar")
```

## Usage

For further examples, see the help files within the package. We can generate an example for change point detection as follows.

Generate a univariate time series of length 1000, with a mean change at time 300, and an autocovariance (but not marginal) change at time 650:

```
n <- 1000
set.seed(123)

noise1 <- stats::arima.sim(model = list(ar = -0.5), n = n, sd = sqrt(1-0.5^2))
noise2 <- stats::arima.sim(model = list(ar = 0.5), n = n, sd = sqrt(1-0.5^2))

noise <- c(noise1[1:650],noise2[651:n])

signal <- c(rep(0,300),rep(0.7,700))

x <- signal + noise

```
Perform the multilag NP-MOJO algorithm with lags 0 and 1:

```
x.c <- np.mojo.multilag(x,G=166, lags = c(0,1))
```

Print the estimated change points and the associated clusters:

```
x.c$merged.cpts$cpts

x.c$merged.cpts$cpt.clusters
```


