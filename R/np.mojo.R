#' @title Nonparametric Single Lag Change Point Detection
#' @description For a given lagged value of the time series, performs nonparametric change point detection of a possibly multivariate
#' time series. If \code{lag} \eqn{\ell = 0}, then only marginal changes are detected.
#' If \code{lag} \eqn{\ell \neq 0}, then changes in the pairwise distribution of \eqn{(X_t , X_{t+\ell})} are detected.
#' @details The single-lag NP-MOJO algorithm for nonparametric change point detection is described in McGonigle, E. T. and Cho, H. (2023)
#' Nonparametric data segmentation in multivariate time series via joint characteristic functions.  \emph{arXiv preprint arXiv:2305.07581}.
#' @param x Input data (a \code{numeric} vector or an object of classes \code{ts} and \code{timeSeries},
#' or a \code{numeric} matrix with rows representing observations and columns representing variables).
#' @param G An integer value for the moving sum bandwidth;
#' \code{G} should be less than half the length of the time series.
#' @param lag The lagged values of the time series used to detect changes. If \code{lag} \eqn{\ell = 0}, then only marginal changes are detected.
#' If \code{lag} \eqn{\ell \neq 0}, then changes in the pairwise distribution of \eqn{(X_t , X_{t+\ell})} are detected.
#' @param kernel.f String indicating which kernel function to use when calculating the NP-MOJO detectors statistics; with \code{kern.par} \eqn{= a}, possible values are
#'  \itemize{
#'    \item \code{"quad.exp"}: kernel \eqn{h_2} in McGonigle and cho (2025), kernel 5 in Fan et al. (2017):
#'    \deqn{h (x,y) = \prod_{i=1}^{2p} \frac{ (2a - (x_i - y_i)^2) \exp (-\frac{1}{4a} (x_i - y_i)^2 )}{2a} .}
#'    \item \code{"gauss"}: kernel \eqn{h_1} in McGonigle and cho (2025), the standard Gaussian kernel:
#'    \deqn{h (x,y) = \exp ( - \frac{a^2}{2} \Vert x - y  \Vert^2) .}
#'    \item \code{"euclidean"}: kernel \eqn{h_3} in McGonigle and cho (2025), the Euclidean distance-based kernel:
#'    \deqn{h (x, y ) = \Vert x - y \Vert^a  .}
#'    \item \code{"laplace"}: kernel 2 in Fan et al. (2017), based on a Laplace weight function:
#'      \deqn{h (x, y ) = \prod_{i=1}^{2p} \left( 1+ a^2 (x_i - y_i)^2  \right)^{-1}. }
#'    \item \code{"sine"}: kernel 4 in Fan et al. (2017), based on a sinusoidal weight function:
#'      \deqn{h (x, y ) = \prod_{i=1}^{2p} \frac{-2 | x_i - y_i |  + | x_i - y_i - 2a|  + | x_i - y_i +2a| }{4a} .}
#' }
#' @param kern.par The tuning parameter that appears in the expression for the kernel function, which acts as a scaling parameter,
#' only to be used if \code{data.driven.kern.par = FALSE}. If \code{kernel.f = "euclidean"}, then \code{kern.par} \eqn{\in (0,2)},
#' otherwise \code{kern.par} \eqn{> 0}.
#' @param data.driven.kern.par A \code{logical} variable, if set to \code{TRUE}, then the kernel tuning parameter is calculated
#'  using the median heuristic, if \code{FALSE} it is given by \code{kern.par}.
#' @param alpha A numeric value for the significance level with
#' \code{0 <= alpha <= 1}; use iff \code{threshold = "bootstrap"}.
#' @param reps An integer value for the number of bootstrap replications performed, if \code{threshold = "bootstrap"}.
#' @param boot.dep A positive value for the strength of dependence in the multiplier bootstrap sequence, if \code{threshold = "bootstrap"}.
#' @param parallel A \code{logical} variable, if set to \code{TRUE}, then parallel computing is used in the bootstrapping procedure
#'  if bootstrapping is performed.
#' @param boot.method A string indicating the method for creating bootstrap replications. It is not recommended to change this. Possible choices are
#' \itemize{
#'    \item \code{"mean.subtract"}: the default choice, as described in McGonigle and cho (2025).
#'    Empirical mean subtraction is performed to the bootstrapped replicates, improving power.
#'        \item \code{"no.mean.subtract"}: empirical mean subtraction is not performed, improving size control.
#' }
#' @param criterion String indicating how to determine whether each point \code{k} at which NP-MOJO statistic
#' exceeds the threshold is a change point; possible values are
#'  \itemize{
#'    \item \code{"epsilon"}: \code{k} is the maximum of its local exceeding environment,
#'    which has at least size \code{epsilon*G}.
#'        \item \code{"eta"}: there is no larger exceeding in an \code{eta*G} environment of \code{k}.
#'        \item \code{"eta.and.epsilon"}: the recommended default option; \code{k} satisfies both
#'        the eta and epsilon criterion.
#'        Recommended to use with the standard value of eta that would be used if \code{criterion = "eta"} (e.g. 0.4),
#'        but much smaller value of epsilon than would be used if \code{criterion = "epsilon"}, e.g. 0.02.
#' }
#' @param eta A positive numeric value for the minimal mutual distance of
#' changes, relative to bandwidth (if \code{criterion = "eta"} or \code{criterion = "eta.and.epsilon"}).
#' @param epsilon a numeric value in (0,1] for the minimal size of exceeding
#' environments, relative to moving sum bandwidth (if \code{criterion = "epsilon"} or \code{criterion = "eta.and.epsilon"}).
#' @param use.mean \code{Logical variable}, only to be used if \code{data.drive.kern.par=TRUE}. If set to \code{TRUE}, the mean
#' of pairwise distances is used to set the kernel function tuning parameter, instead of the median. May be useful for binary data,
#' not recommended to be used otherwise.
#' @param threshold String indicating how the threshold is computed. Possible values are
#'  \itemize{
#'    \item \code{"bootstrap"}: the threshold is calculated using the bootstrap method
#'    with significance level \code{alpha}.
#'        \item \code{"manual"}: the threshold is set by the user and must be
#'        specified using the \code{threshold.val} parameter.
#' }
#' @param threshold.val The value of the threshold used to declare change points, only to be used if \code{threshold = "manual"}.
#'
#' @return A \code{list} object that contains the following fields:
#'    \item{x}{Input data}
#'    \item{G}{Moving window bandwidth}
#'    \item{lag}{Lag used to detect changes}
#'    \item{kernel.f, data.driven.kern.par, use.mean}{Input parameters}
#'    \item{kern.par}{The value of the kernel tuning parameter}
#'    \item{threshold, alpha, reps, boot.dep, boot.method, parallel}{Input parameters}
#'    \item{threshold.val}{Threshold value for declaring change points}
#'    \item{criterion, eta, epsilon}{Input parameters}
#'    \item{test.stat}{A vector containing the NP-MOJO detector statistics computed from the input data}
#'    \item{cpts}{A vector containing the estimated change point locations}
#'    \item{p.vals}{The corresponding p values of the change points, if the bootstrap method was used}
#' @references McGonigle, E.T., Cho, H. (2025). Nonparametric data segmentation in multivariate time series via joint characteristic functions. \emph{Biometrika}.
#' @references Fan, Y., de Micheaux, P.L., Penev, S. and Salopek, D. (2017). Multivariate nonparametric test of independence. \emph{Journal of Multivariate Analysis},
#' 153, pp.189-210.
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' n <- 500
#' noise <- c(rep(1, 300), rep(0.4, 200)) * stats::arima.sim(model = list(ar = 0.3), n = n)
#' signal <- c(rep(0, 100), rep(2, 400))
#' x <- signal + noise
#' x.c <- np.mojo(x, G = 83, lag = 0)
#' x.c$cpts
#' x.c$p.vals
#' @importFrom Rcpp evalCpp
#' @useDynLib CptNonPar, .registration = TRUE
#' @importFrom foreach %dopar%
#' @seealso \link{np.mojo.multilag}
np.mojo <- function(x, G, lag = 0, kernel.f = c("quad.exp", "gauss", "euclidean", "laplace", "sine")[1],
                    kern.par = 1, data.driven.kern.par = TRUE, alpha = 0.1, threshold = c("bootstrap", "manual")[1],
                    threshold.val = NULL, reps = 200, boot.dep = 1.5 * (nrow(as.matrix(x))^(1 / 3)), parallel = FALSE,
                    boot.method = c("mean.subtract", "no.mean.subtract")[1],
                    criterion = c("eta", "epsilon", "eta.and.epsilon")[3], eta = 0.4, epsilon = 0.02, use.mean = FALSE) {
  mojo.error.checks(
    x = x, G = G, lag = lag, kernel.f = kernel.f, kern.par = kern.par, data.driven.kern.par = data.driven.kern.par,
    alpha = alpha, threshold = threshold, threshold.val = threshold.val, reps = reps, boot.dep = boot.dep,
    parallel = parallel, boot.method = boot.method, criterion = criterion, eta = eta, epsilon = epsilon,
    use.mean = use.mean
  )

  if (is.null(dim(x))) {
    data.len <- length(x)
    x <- matrix(x)
  } else {
    data.len <- nrow(x)
  }

  if (G > data.len / 2) {
    stop("Bandwidth is too large for the length of time series.")
  }

  if (lag != 0) {
    y <- cbind(x[1:(data.len - lag), ], x[(lag + 1):data.len, ])
  } else {
    y <- x
  }

  # calculate the Euclidean distances required for the test statistic calculation:

  if (data.driven.kern.par == TRUE) {
    if (kernel.f == "euclidean") {
      warning("Data driven parameter choice not suited for Euclidean kernel. Parameter kern.par reset to 1.")
      kern.par <- 1
    } else {
      D.par <- mosum_dist_calc(y = y, G = G, n = data.len - lag, kern = "euc.dist", kern_par = kern.par)
      if (use.mean == TRUE) {
        delta <- mean(D.par[D.par > 0])
      } else {
        delta <- stats::median(D.par[D.par > 0])
      }
      if (kernel.f == "gauss") {
        kern.par <- 1 / sqrt(delta)
      } else if (kernel.f == "quad.exp") {
        kern.par <- 0.5 * delta
      }
    }
  }


  if (kernel.f == "gauss" && data.driven.kern.par == TRUE) {
    D <- D.par
  } else {
    D <- mosum_dist_calc(y = y, G = G, n = data.len - lag, kern = kernel.f, kern_par = kern.par)

    if (kernel.f != "euclidean" && kernel.f != "gauss") {
      diag(D) <- 1
    }
  }

  # calculate the entries of the kernel statistic:

  h.mat <- h.matcalc(D = D, lag = lag, G = G, kernel.f = kernel.f, kern.par = kern.par)

  h.mat.ncol <- ncol(h.mat)

  init.val <- sum(h.mat[1:(G - lag), 1:(G - lag)])

  # calculate the MOSUM test statistic to find change points:

  test.stat <- rolling_matrix_sum(stat_mat = h.mat, G = G, lag = lag, init_val = init.val, n = h.mat.ncol)

  test.stat <- c(rep(0, G - 1), test.stat[1:(data.len - 2 * G + 1)] / ((G - lag)^2), rep(0, G))

  # perform wild multiplier bootstrap to assess significance of changes:

  if (threshold == "bootstrap") {
    if (parallel == TRUE) {
      cl <- parallel::makeCluster(parallelly::availableCores())
      doParallel::registerDoParallel(cl)

      Tstar <- foreach::foreach(iterators::icount(reps), .combine = cbind, .packages = c("foreach", "Rcpp", "CptNonPar")) %dopar% bootstrap.char(h.mat, test.stat, G, lag, reps = 1, boot.dep, parallel, boot.method, data.len, h.mat.ncol)

      parallel::stopCluster(cl)
    } else {
      Tstar <- bootstrap.char(
        h.mat = h.mat, test.stat = test.stat, G = G, lag = lag, reps = reps, boot.dep = boot.dep,
        parallel = parallel, boot.method = boot.method, data.len = data.len, h.mat.ncol = h.mat.ncol
      )
    }

    threshold.val <- stats::quantile(Tstar, probs = 1 - alpha)
  }

  cpt.locs <- numeric(0)
  p.vals <- numeric(0)

  if (max(test.stat) > threshold.val) {
    exceedings <- (test.stat > threshold.val)

    if (criterion == "epsilon") {
      exceedingsCount <- (exceedings) * unlist(lapply(
        rle(exceedings)$lengths,
        seq_len
      ))
      minIntervalSize <- max(1, G * epsilon)
      intervalEndPoints <- which(diff(exceedingsCount) <= -minIntervalSize)
      intervalBeginPoints <- intervalEndPoints - exceedingsCount[intervalEndPoints] + 1

      numChangePoints <- length(intervalBeginPoints)
      if (numChangePoints > 0) {
        for (i in 1:numChangePoints) {
          changePoint <- intervalBeginPoints[i] + which.max(test.stat[(intervalBeginPoints[i]):(intervalEndPoints[i])]) -
            1
          cpt.locs <- c(cpt.locs, changePoint)
        }
      }
    } else {
      localMaxima <- (c((diff.default(test.stat) < 0), NA) & c(NA, diff.default(test.stat) > 0))
      localMaxima[data.len - G] <- TRUE

      if (criterion == "eta.and.epsilon") {
        stat.exceed <- which(test.stat > threshold.val)
        r <- rle(diff(stat.exceed))
        end.r <- cumsum(r$lengths)
        start.r <- end.r - r$lengths + 1

        epsilon.satisfied.start <- stat.exceed[start.r[r$lengths > epsilon * G]]
        epsilon.satisfied.end <- stat.exceed[end.r[r$lengths > epsilon * G]] + 1

        epsilon.exceedings <- rep(FALSE, data.len)
        for (i in seq_len(length(epsilon.satisfied.start))) {
          epsilon.exceedings[epsilon.satisfied.start[i]:epsilon.satisfied.end[i]] <- TRUE
        }

        p.candidates <- which(epsilon.exceedings & localMaxima)
      } else {
        p.candidates <- which(localMaxima)
      }

      cpt.locs <- mojo_eta_criterion_help(p.candidates, test.stat, eta, G, G)
    }

    if (threshold == "bootstrap") {
      p.vals <- numeric(0)

      for (i in seq_len(length(cpt.locs))) {
        p.vals <- c(p.vals, sum(Tstar >= test.stat[cpt.locs[i]]) / reps)
      }
    }
  } else {
    if (threshold == "bootstrap") {
      p.vals <- sum(Tstar >= max(test.stat)) / reps
    }
  }


  ret <- list(
    x = x,
    G = G,
    lag = lag,
    kernel.f = kernel.f,
    kern.par = kern.par,
    data.driven.kern.par = data.driven.kern.par,
    threshold = threshold,
    threshold.val = threshold.val,
    boot.dep = boot.dep,
    boot.method = boot.method,
    reps = reps,
    parallel = parallel,
    alpha = alpha,
    criterion = criterion,
    eta = eta,
    epsilon = epsilon,
    use.mean = use.mean,
    test.stat = test.stat,
    cpts = cpt.locs,
    p.vals = p.vals
  )

  return(ret)
}
