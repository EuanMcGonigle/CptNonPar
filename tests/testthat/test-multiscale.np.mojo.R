test_that("multiscale NP-MOJO executes", {
  n <- 500
  noise <- c(rep(1, 300), rep(0.4, 200)) * stats::rnorm(n)
  signal <- c(rep(0, 100), rep(2, 400))
  x <- signal + noise
  x.c <- multiscale.np.mojo(x, G = c(50, 80), lags = c(0, 1), reps = 10)
  expect_equal(class(x.c), "list")
})

test_that("multiscale NP-MOJO rejects incorrect bandwidths", {
  expect_error(
    multiscale.np.mojo(stats::rnorm(100), G = c(-10, 30)),
    "Bandwidth parameter G must be numeric positive integer."
  )
})

test_that("multiscale NP-MOJO rejects non-numeric lags", {
  expect_error(
    multiscale.np.mojo(stats::rnorm(100), G = c(10, 30), lags = c("1", "2")),
    "The set of lags must be a numeric vector of positive integer values."
  )
})

test_that("multiscale NP-MOJO rejects negative lags", {
  expect_error(
    multiscale.np.mojo(stats::rnorm(100), G = c(10, 30), lags = c(-1, 0, 3)),
    "The set of lags must be a numeric vector of positive integer values."
  )
})

test_that("multiscale NP-MOJO executes with manual threshold and single numeric threshold.val specification", {
  x.c <- multiscale.np.mojo(stats::rnorm(100), G = c(20, 30), lags = c(0, 1),
                            threshold = "manual", threshold.val = 0.1)
  expect_equal(class(x.c), "list")
})

test_that("multiscale NP-MOJO executes with manual threshold and list threshold.val specification", {
  x.c <- multiscale.np.mojo(stats::rnorm(100), G = c(20, 30), lags = c(0, 1),
                            threshold = "manual", threshold.val = list(c(0.1, 0.2), c(0.15,0.25)))
  expect_equal(class(x.c), "list")
})

test_that("multiscale NP-MOJO rejects incorrect threshold.val specification", {
  expect_error(
    multiscale.np.mojo(stats::rnorm(100), G = c(10, 30), threshold = "manual",
                       threshold.val = c(1, 2)),
    "With threshold = manual, threshold.val must be either a single numeric value, or list of vectors, with list length length(G), and each vector of length length(lags).",
    fixed = TRUE
  )
})




