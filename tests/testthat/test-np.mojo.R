#Tests for np.mojo.multilag function

set.seed(1)
x <- rnorm(500)+c(rep(0,200),rep(2,300))
G <- 80

test_that("Change point merging type is recognised",{
  expect_error(np.mojo.multilag(x, merge.type = "pruning"),
               "Error: change point merging type must be either 'sequential' or 'bottom-up'.")
})
