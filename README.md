<!-- badges: start -->
[![R-CMD-check](https://github.com/EuanMcGonigle/CptNonPar/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/EuanMcGonigle/CptNonPar/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->
# CptNonPar
Nonparametric change point detection for multivariate time series. Implements the NP-MOJO methdology proposed in

> McGonigle, E. T., Cho, H. (2023). Nonparametric data segmentation in multivariate time series via joint characteristic functions. 
> [(link to paper here)](https://arxiv.org/abs/2305.07581).

## Installation

To install `CptNonPar` from GitHub:

```
devtools::install_github("https://github.com/EuanMcGonigle/CptNonPar")
```

## Usage

For detailed examples, see the help files within the package. We can generate a small example for change point detection as follows:

```
set.seed(123)
n <- 1000
noise <- c(rep(1,300),rep(0.4,700))*stats::arima.sim(model = list(ar =0.3), n = 1000)
signal <- c(rep(0,700),rep(0.5,300))
x <- signal + noise
x.c <- np.mojo(x, G = 166, lag = 0)
x.c$cpts



