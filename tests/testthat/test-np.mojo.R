#Tests for np.mojo.multilag function

set.seed(1)
x <- rnorm(500)+c(rep(0,200),rep(2,300))
G <- 80

x.c <- np.mojo.multilag(x, G = G)

test_that("np.mojo.multilag executes", {
  skip_on_cran()
  expect_equal(class(x.c), "list")
})

test_that("Change point merging type is recognised",{
  expect_error(np.mojo.multilag(x, G = G, merge.type = "pruning"),
               "Error: change point merging type must be either 'sequential' or 'bottom-up'.")
})

test_that("lags parameter is numeric",{
  expect_error(np.mojo.multilag(x, G = G, lags = c("0", "1")),
               "The set of lags must be a numeric vector of positive integer values.")
})

test_that("lags parameter contains positive values",{
  expect_error(np.mojo.multilag(x, G = G, lags = c(0,-1)),
               "The set of lags must be a numeric vector of positive integer values.")
})
