#Tests for multilag.cpts.merge function

n <- 1000
set.seed(123)

noise1 <- stats::arima.sim(model = list(ar = -0.5), n = n, sd = sqrt(1 - 0.5^2))
noise2 <- stats::arima.sim(model = list(ar = 0.5), n = n, sd = sqrt(1 - 0.5^2))

noise <- c(noise1[1:650], noise2[651:n])
signal <- c(rep(0, 300), rep(0.7, 700))
x <- signal + noise

set.seed(1)
x.nocp <- stats:::rnorm(1000)

G <- 166

x.c0 <- np.mojo(x,G,0)
x.c1 <- np.mojo(x,G,1)

test_that("multilag.cpts.merge executes with sequential merging", {
  skip_on_cran()
  expect_equal(class(multilag.cpts.merge(x.c = list(x.c0,x.c1), merge.type = "sequential")), "list")
})

test_that("multilag.cpts.merge executes with bottom-up merging", {
  skip_on_cran()
  expect_equal(class(multilag.cpts.merge(x.c = list(x.c0, x.c1), merge.type = "bottom-up")), "list")
})

test_that("Change point merging type is recognised",{
  expect_error(multilag.cpts.merge(x.c = list(x.c0, x.c1), merge.type = "pruning"),
               "Error: change point merging type must be either 'sequential' or 'bottom-up'.")
})

test_that("multilag.cpts.merge executes when no change points", {
  skip_on_cran()
  x.nocp0 <- np.mojo(x.nocp,G,0, threshold = "manual", threshold.val = 5)
  x.nocp1 <- np.mojo(x.nocp,G,1, threshold = "manual", threshold.val = 5)
  expect_equal(class(multilag.cpts.merge(x.c = list(x.nocp0,x.nocp1), merge.type = "sequential")), "list")
})
