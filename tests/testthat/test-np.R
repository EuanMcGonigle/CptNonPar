#Tests for np.mojo function

set.seed(1)
x <- rnorm(500)+c(rep(0,200),rep(2,300))
G <- 80

x.c <- np.mojo(x,G)


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



