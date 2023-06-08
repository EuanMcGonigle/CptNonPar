#Tests for np.mojo function

set.seed(1)
x <- stats::rnorm(500)+c(rep(0,200),rep(2,300))
G <- 80

kernels <- c('quad.exp', 'gauss', 'laplace', 'sine', 'euclidean')
criteria <- c('eta', 'epsilon', 'eta.and.epsilon')

for(i in seq_len(length(kernels))){
  test_that(paste0('np.mojo with Kernel #',i, " executes with data.driven.par FALSE"), {
    skip_on_cran()
    expect_equal(class(np.mojo(x,G, kernel.f = kernels[i], data.driven.kern.par = FALSE)), "list")
  })
}

for(i in (seq_len(length(kernels)-1))){
  test_that(paste0('np.mojo with Kernel #',i, " executes with data.driven.par TRUE"), {
    skip_on_cran()
    expect_equal(class(np.mojo(x,G, kernel.f = kernels[i], data.driven.kern.par = TRUE)), "list")
  })
}

for(i in seq_len(length(criteria))){
  test_that(paste0('np.mojo with criteria #',i, " executes"), {
    skip_on_cran()
    expect_equal(class(np.mojo(x,G, criterion = criteria[i])), "list")
  })
}

test_that("np.mojo executes with param 'boot.method' = 'no.mean.subtract'", {
  skip_on_cran()
  expect_equal(class(np.mojo(x,G, boot.method = "no.mean.subtract")), "list")
})


test_that("np.mojo executes with default params", {
  skip_on_cran()
  expect_equal(class(np.mojo(x,G)), "list")
})

test_that("np.mojo executes with parallel", {
   skip_on_cran()
   expect_equal(class(np.mojo(x,G, parallel = TRUE)), "list")
 })

test_that("np.mojo executes with use.mean", {
  skip_on_cran()
  expect_equal(class(np.mojo(x,G, use.mean = TRUE)), "list")
})

test_that("np.mojo executes with matrix data", {
  skip_on_cran()
  X <- matrix(stats::rnorm(1000), ncol = 2, nrow = 500)
  expect_equal(class(np.mojo(X,G)), "list")
})


test_that("Criterion choice is recognised",{
  expect_error(np.mojo(x = x, G = G, criterion = "eat"),
               "Error: change point detection criterion must be one of 'eta', 'epsilon', or 'eta.and.epsilon'.")
})

test_that("Kernel.f choice is recognised",{
  expect_error(np.mojo(x = x, G = G, kernel.f = "normal"),
               "The kernel.f function must be either 'quad.exp', 'gauss', 'laplace', 'sine', or 'euclidean'")
})

test_that("boot.method choice is recognised",{
  expect_error(np.mojo(x = x, G = G, boot.method = "multiplier"),
               "Parameter 'boot.method' must be set to be either 'mean.subtract' or 'no.mean.subtract'. Highly
         recommended to set boot.method = 'mean.subtract'.")
})


test_that("parallel argument is logical",{
  expect_error(np.mojo(x = x, G = G, parallel = "true"),
               "Error: 'parallel' argument must be logical variable.")
})

test_that("use.mean argument is logical",{
  expect_error(np.mojo(x = x, G = G, use.mean = "true"),
               "Error: 'use.mean' argument must be logical variable.")
})

test_that("data.driven.kern.par argument is logical",{
  expect_error(np.mojo(x = x, G = G, data.driven.kern.par = "true"),
               "Error: 'data.driven.kern.par' argument must be logical variable.")
})


test_that("boot.dep argument is positive",{
  expect_error(np.mojo(x = x, G = G, boot.dep = -2),
               "The bootstrap dependence parameter 'boot.dep' must be a positive value.")
})

test_that("kern.par argument is positive",{
  expect_error(np.mojo(x = x, G = G, kern.par = -2),
               "The kernel parameter must be a positive value.")
})


test_that("reps argument is positive",{
  expect_error(np.mojo(x = x, G = G, reps = -2),
               "Number of bootstrap replications should be a single positive integer.")
})

test_that("reps argument is numeric",{
  expect_error(np.mojo(x = x, G = G, reps = "2"),
               "Number of bootstrap replications should be a single positive integer.")
})

test_that("lag argument is positive",{
  expect_error(np.mojo(x = x, G = G, lag = -1),
               "The lag parameter should be a single positive integer.")
})

test_that("lag argument is numeric",{
  expect_error(np.mojo(x = x, G = G, lag = "1"),
               "The lag parameter should be a single positive integer.")
})

test_that("G is numeric",{
  expect_error(np.mojo(x = x, G =  "100"),
               "Bandwidth parameter G must be numeric positive integer.")
})

test_that("G is larger than lag",{
  expect_error(np.mojo(x = x, G =  20, lag = 25),
               "Bandwidth parameter G must be larger than the chosen lag.")
})

test_that("x has no NAs",{
  expect_error(np.mojo(x = c(x[1:499],NA), G =  G),
               "Missing values in data: NA is not allowed in the data.")
})

test_that("x is numeric",{
  expect_error(np.mojo(x = rep("1",500), G =  G),
               "Data must be numeric.")
})


test_that("threshold is recognised",{
  expect_error(np.mojo(x = x, G =  G, threshold = "asymptotic"),
               "Threshold type parameter 'threshold' must be set to be either 'bootstrap' or 'manual'.")
})


test_that("threshold type 'manual' has threshold.val set",{
  expect_error(np.mojo(x = x, G =  G, threshold = "manual"),
               "Threshold type has been set to 'manual', but threshold.val has not been set.")
})

test_that("threshold type 'manual' has positive numeric threshold.val",{
  expect_error(np.mojo(x = x, G =  G, threshold = "manual", threshold.val = -0.05),
               "Parameter threshold.val must be a nonnegative number.")
})

test_that("reps only used with threshold = 'bootstrap",{
  skip_on_cran()
  expect_warning(np.mojo(x = x, G =  G, threshold = "manual", threshold.val = 0.1, reps = 200),
                 "reps is only used with threshold=bootstrap")
})

test_that("Euclidean kernel has appropriate kern.par value",{
  expect_error(np.mojo(x = x, G =  G, kernel.f = "euclidean", kern.par = 10, data.driven.kern.par = FALSE),
               "For the 'euclidean' kernel, the kernel parameter must be in the interval (0,2).",
               fixed = TRUE)
})

test_that("warning when Euclidean kernel used with data.driven.par = TRUE",{
  expect_warning(np.mojo(x = x, G =  G, kernel.f = "euclidean", data.driven.kern.par = TRUE),
               "Data driven parameter choice not suited for Euclidean kernel. Parameter kern.par reset to 1.")
})

test_that("G is not too large",{
  expect_error(np.mojo(x = x, G =  400),
                 "Bandwidth is too large for the length of time series.")
})

