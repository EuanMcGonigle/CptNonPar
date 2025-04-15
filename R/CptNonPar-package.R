#' @keywords internal
#' @examples
#' set.seed(1)
#' n <- 500
#' noise <- c(rep(1, 300), rep(0.4, 200)) * stats::arima.sim(model = list(ar = 0.3), n = n)
#' signal <- c(rep(0, 100), rep(2, 400))
#' x <- signal + noise
#' x.c <- np.mojo.multilag(x, G = 83, lags = c(0, 1))
#' x.c$cpts
#' x.c$cpt.clusters
#' @seealso \link{np.mojo}, \link{np.mojo.multilag}, \link{multilag.cpts.merge}
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL


