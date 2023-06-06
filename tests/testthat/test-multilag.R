#Tests for multilag.cpts.merge function

set.seed(1)
x <- rnorm(500)+c(rep(0,200),rep(2,300))
G <- 80

x.c0 <- np.mojo(x,G,0)
x.c1 <- np.mojo(x,G,1)

test_that("Change point merging type is recognised",{
  expect_error(multilag.cpts.merge(x.c = list(x.c0, x.c1), merge.type = "pruning"),
               "Error: change point merging type must be either 'sequential' or 'bottom-up'.")
})
