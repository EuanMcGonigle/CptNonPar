#Tests for multilag.cpts.merge function

set.seed(1)
x <- rnorm(500)+c(rep(0,200),rep(2,300))
G <- 80

x.c0 <- np.mojo(x,G,0)
x.c1 <- np.mojo(x,G,1)

x.c <- multilag.cpts.merge(x.c = list(x.c0, x.c1))

test_that("multilag.cpts.merge executes", {
  skip_on_cran()
  expect_equal(class(x.c), "list")
})

test_that("Change point merging type is recognised",{
  expect_error(multilag.cpts.merge(x.c = list(x.c0, x.c1), merge.type = "pruning"),
               "Error: change point merging type must be either 'sequential' or 'bottom-up'.")
})
