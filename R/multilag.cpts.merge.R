#' @title Merge Change Point Estimators from Multiple Lags
#' @description Merges change point estimators from different lagged values into a final set of overall change point estimators.
#' @details See McGonigle and Cho (2025) for further details.
#' @param x.c A \code{list} object, where each element of the list is the output of the \link{np.mojo} function computed at a
#' different lag.
#' @param eta.merge A positive numeric value for the minimal mutual distance of
#' changes, relative to bandwidth, used to merge change point estimators across different lags.
#' @param merge.type String indicating the method used to merge change point estimators from different lags. Possible choices are
#'  \itemize{
#'    \item \code{"sequential"}: starting from the left-most change point estimator and proceeding forward in time, estimators
#'    are grouped into clusters based on mutual distance. The estimator yielding the largest corresponding importance score is
#'    chosen as the change point estimator for that cluster. See McGonigle and Cho (2025) for details.
#'        \item \code{"bottom-up"}: starting with the largest importance score, the change points are merged using bottom-up merging (Messer
#'        et al. (2014)).
#' }
#'
#' @return A \code{list} object which contains the following fields
#' \item{cpts}{A matrix with rows corresponding to final change point estimators, with estimated change point location and associated lag and importance score given in columns.}
#'    \item{cpt.clusters}{A \code{list} object of length given by the number of detected change points. Each field contains a matrix of all
#'    change point estimators that are declared to be associated to the corresponding change point in the \code{cpts} field.}
#' @references McGonigle, E.T., Cho, H. (2025). Nonparametric data segmentation in multivariate time series via joint characteristic functions.
#' \emph{Biometrika} (to appear).
#' @references Messer M., Kirchner M., Schiemann J., Roeper J., Neininger R., Schneider G. (2014). A Multiple Filter Test for
#' the Detection of Rate Changes in Renewal Processes with Varying Variance. \emph{The Annals of Applied Statistics}, 8(4), 2027-2067.
#' @export
#'
#' @examples
#' set.seed(1)
#' n <- 500
#' noise <- c(rep(1, 300), rep(0.4, 200)) * stats::arima.sim(model = list(ar = 0.3), n = n)
#' signal <- c(rep(0, 100), rep(2, 400))
#' x <- signal + noise
#' x.c0 <- np.mojo(x, G = 83, lag = 0)
#' x.c1 <- np.mojo(x, G = 83, lag = 1)
#' x.c <- multilag.cpts.merge(list(x.c0, x.c1))
#' x.c
#' @seealso \link{np.mojo}, \link{np.mojo.multilag}
multilag.cpts.merge <- function(x.c, eta.merge = 1, merge.type = c("sequential", "bottom-up")[1]) {
  stopifnot(
    "Error: change point merging type must be either 'sequential' or 'bottom-up'." =
      merge.type == "sequential" || merge.type == "bottom-up"
  )

  cpts <- init.cpts <- matrix(NA, nrow = 0, ncol = 3)
  dimnames(init.cpts)[[2]] <- c("cpt", "lag", "score")
  dimnames(cpts)[[2]] <- c("cpt", "lag", "score")

  cpt.clusters <- list()

  L <- length(x.c)

  G <- x.c[[1]]$G

  for (l in 1:L) {
    lag.cpts <- x.c[[l]]
    if(lag.cpts$threshold == "manual"){
      lag.cpts$scores <- lag.cpts$test.stat[lag.cpts$cpts]
    }

    if (length(lag.cpts$cpts) > 0) {
      new.cpts <- cbind(lag.cpts$cpts, rep(lag.cpts$lag, length(lag.cpts$cpts)), lag.cpts$scores)
      init.cpts <- rbind(init.cpts, new.cpts)
    }
  }

  if (nrow(init.cpts) <= 1) {
    return(list(cpts = init.cpts, cpt.clusters = init.cpts))
  }

  if (merge.type == "sequential") {
    cpt.order <- order(init.cpts[, 1])

    init.cpts <- init.cpts[cpt.order, ]

    j <- 1

    while (nrow(init.cpts) != 0) {
      if (nrow(init.cpts) == 1) {
        cpts <- rbind(cpts, init.cpts)
        cpt.clusters[[j]] <- init.cpts
      } else {
        cpt.left <- init.cpts[1, 1]
        linked.cpts <- init.cpts[2:nrow(init.cpts), 1]
        linked.cpts <- c(cpt.left, linked.cpts[linked.cpts - cpt.left < eta.merge * G])

        num.linked.cpts <- length(linked.cpts)

        linked.lags <- init.cpts[1:num.linked.cpts, 2]
        linked.scores <- init.cpts[1:num.linked.cpts, 3]

        final.cpt <- linked.cpts[which.max(linked.scores)]
        final.lag <- linked.lags[which.max(linked.scores)]
        final.pval <- linked.scores[which.max(linked.scores)]

        cpts <- rbind(cpts, c(final.cpt, final.lag, final.pval))

        cpt.clusters[[j]] <- init.cpts[1:num.linked.cpts, , drop = FALSE][order(init.cpts[1:num.linked.cpts, 2]), , drop = FALSE]
      }
      init.cpts <- init.cpts[-(1:num.linked.cpts), , drop = FALSE]
      j <- j + 1
    }
  } else {
    cpt.order <- order(-init.cpts[, 3]) #larger importance score is better

    init.cpts <- init.cpts[cpt.order, ]

    j <- 1

    while (nrow(init.cpts) != 0) {
      if (nrow(init.cpts) == 1) {
        cpts <- rbind(cpts, init.cpts)
        cpt.clusters[[j]] <- init.cpts
      } else {
        cpt.min <- init.cpts[1, 1]

        nearby.cpts <- which(abs(init.cpts[, 1] - cpt.min) < eta.merge * G)
        linked.cpts <- init.cpts[nearby.cpts, , drop = FALSE]

        cpts <- rbind(cpts, init.cpts[1, ])
        cpt.clusters[[j]] <- linked.cpts[order(linked.cpts[, 2]), ]
      }
      init.cpts <- init.cpts[-nearby.cpts, , drop = FALSE]
      j <- j + 1
    }

    final.cpt.order <- order(cpts[, 1])
    cpts <- cpts[final.cpt.order, ]
    cpt.clusters <- cpt.clusters[final.cpt.order]
  }

  return(list(cpts = cpts, cpt.clusters = cpt.clusters))
}
