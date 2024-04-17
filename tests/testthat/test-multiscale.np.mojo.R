test_that("multiscale NP-MOJO executes", {
  n <- 500
  noise <- c(rep(1, 300), rep(0.4, 200)) * stats::rnorm(n)
  signal <- c(rep(0, 100), rep(2, 400))
  x <- signal + noise
  x.c <- multiscale.np.mojo(x, G = c(50,80), lags = c(0, 1), reps = 10)
  expect_equal(class(x.c), "list")
})

test_that("multiscale NP-MOJO rejects incorrect bandwidths", {
  expect_error(multiscale.np.mojo(stats::rnorm(100), G = c(-10,30)),
               "Bandwidth parameter G must be numeric positive integer.")
})
