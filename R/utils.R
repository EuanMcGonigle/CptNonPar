#' @title Compute Running Mean
#' @description Internal function for performing running mean used in the bootstrapping procedure
#' @keywords internal
#' @noRd

runmean <- function(x, n = 5) {
  stats::filter(x, rep(1 / n, n), sides = 1)
}

#' @title Performing multiplier bootstrapping procedure
#' @description Internal function for performing bootstrapping procedure with implementation using Rcpp
#' @importFrom Rcpp evalCpp
#' @useDynLib CptNonPar, .registration = TRUE
#' @keywords internal
#' @noRd

bootstrap.tstat <- function(h.mat, test.stat, G, lag, boot.dep, boot.method, data.len, h.mat.ncol) {
  # set.seed(1)

  Wtstar <- stats::arima.sim(
    model = list(ar = exp(-1 / boot.dep)), sd = sqrt(1 - exp(-2 / boot.dep)),
    n = h.mat.ncol, n.start = 1, start.innov = stats::rnorm(1)
  )

  t.stat <- rep(0, h.mat.ncol)

  if (boot.method == "mean.subtract") {
    t.stat <- rep(0, h.mat.ncol + G + lag)

    W.mean <- runmean(Wtstar, G - lag)

    Wtstar.out.a <- Rfast::Outer(Wtstar, Wtstar)

    h.mat.a <- C_matmatprod_elwise_inplace(Wtstar.out.a, h.mat) # FIRST ARGUMENT GETS OVER-WRITTEN

    h.mat.b <- C_matvecprod_elwise(h.mat, Wtstar)

    init.val.a <- sum(h.mat.a[1:(G - lag), 1:(G - lag)])

    init.val.b <- sum(h.mat.b[1:(G - lag), 1:(G - lag)])

    t.stat.a <- rolling_matrix_sum(stat_mat = h.mat.a, G = G, lag = lag, init_val = init.val.a, n = h.mat.ncol)
    t.stat.b <- rolling_matrix_sum(stat_mat = h.mat.b, G = G, lag = lag, init_val = init.val.b, n = h.mat.ncol)

    t.stat.a <- c(rep(0, G - 1), t.stat.a[1:(data.len - 2 * G + 1)] / ((G - lag)^2), rep(0, G))

    t.stat.b <- c(rep(0, G - 1), t.stat.b[1:(data.len - 2 * G + 1)] / ((G - lag)^2), rep(0, G))

    t.stat[1:(h.mat.ncol)] <- t.stat.a[1:(h.mat.ncol)] - 2 * W.mean * t.stat.b[1:(h.mat.ncol)] +
      W.mean^2 * test.stat[1:(h.mat.ncol)]
  } else {
    Wtstar.out <- outer(Wtstar, Wtstar)

    h.mat2 <- Wtstar.out * h.mat

    init.val <- sum(h.mat2[1:(G - lag), 1:(G - lag)])

    t.stat <- rolling_matrix_sum(stat_mat = h.mat2, G = G, lag = lag, init_val = init.val, n = h.mat.ncol)

    t.stat <- t.stat / ((G - lag)^2)
  }

  max(t.stat, na.rm = TRUE)
}

#' @title Compute Entries of the h Kernel Matrix
#' @description Internal function for computing the entries of the h kernel matrix required to compute the NP-MOJO detector statistics
#' @keywords internal
#' @noRd

h.matcalc <- function(D, lag, G, kernel.f, kern.par) {

  n <- dim(D)[2]

  h.dim <- n - G

  if (kernel.f == "gauss") {
    h.mat <- exp(-kern.par^2 * (D[1:h.dim, 1:h.dim]) / 2) + exp(-kern.par^2 * (D[(G + 1):(n), (G + 1):(n)]) / 2) -
      exp(-kern.par^2 * (D[1:h.dim, (G + 1):(n)]) / 2) - exp(-kern.par^2 * (D[(G + 1):(n), 1:h.dim]) / 2)
  } else if (kernel.f == "euclidean") {
    h.mat <- (D[1:h.dim, (G + 1):(n)]) + (D[(G + 1):n, 1:h.dim]) - (D[1:h.dim, 1:h.dim]) -
      (D[(G + 1):(n), (G + 1):(n)])
  } else {
    h.mat <- -(D[1:h.dim, (G + 1):(n)]) - D[(G + 1):n, 1:h.dim] +
      (D[1:h.dim, 1:h.dim]) + (D[(G + 1):(n), (G + 1):(n)])
  }

  h.mat
}

#' @title Generate Bootstrap Replicates
#' @description Internal helper function for performing bootstrapping procedure
#' @keywords internal
#' @noRd

bootstrap.char <- function(h.mat, test.stat, G, lag, reps, boot.dep, parallel, boot.method, data.len, h.mat.ncol) {
  if (parallel == FALSE) {
    d <- replicate(reps, bootstrap.tstat(
      h.mat = h.mat, test.stat = test.stat, G = G, lag = lag, boot.dep = boot.dep,
      boot.method = boot.method, data.len = data.len, h.mat.ncol = h.mat.ncol
    ))
  } else {
    d <- bootstrap.tstat(h.mat, test.stat, G, lag, boot.dep, boot.method, data.len = data.len, h.mat.ncol)
  }

  d
}

#' @title Perform Error Checks
#' @description Internal function for performing error checks for the inputs to the np,mojo function
#' @keywords internal
#' @noRd

mojo.error.checks <- function(x, G, lag, kernel.f, kern.par, data.driven.kern.par, alpha, threshold,
                              threshold.val, reps, boot.dep, parallel, boot.method, criterion, eta, epsilon,
                              use.mean) {
  stopifnot("Error: alpha must be a number between 0 and 1." = alpha >= 0 && alpha <= 1)

  stopifnot("Error: change point detection criterion must be one of 'eta', 'epsilon', or 'eta.and.epsilon'." = criterion == "epsilon" || criterion == "eta" || criterion == "eta.and.epsilon")

  stopifnot("Error: epsilon must be a positive number." = criterion != "epsilon" || epsilon >= 0)
  stopifnot("Error: eta must be a positive number." = criterion != "eta" || eta >= 0)
  stopifnot("Error: epsilon must be a positive number." = criterion != "eta.and.epsilon" || epsilon >= 0)
  stopifnot("Error: eta must be a positive number." = criterion != "eta.and.epsilon" || eta >= 0)

  stopifnot("Error: 'parallel' argument must be logical variable." = is.logical(parallel))
  stopifnot("Error: 'data.driven.kern.par' argument must be logical variable." = is.logical(data.driven.kern.par))
  stopifnot("Error: 'use.mean' argument must be logical variable." = is.logical(use.mean))

  if (!is.numeric(lag)) {
    stop("The lag parameter should be a single positive integer.")
  }
  if ((length(lag) != 1) || (lag %% 1 != 0) || (lag < 0)) {
    stop("The lag parameter should be a single positive integer.")
  }

  if(!is.numeric(G) || (G < 0)) {
    stop("Bandwidth parameter G must be numeric positive integer.")
  }

  if(G <= lag) {
    stop("Bandwidth parameter G must be larger than the chosen lag.")
  }

  if (!is.numeric(x)) {
    stop("Data must be numeric.")
  }
  if (any(is.na(x))) {
    stop("Missing values in data: NA is not allowed in the data.")
  }

  if (threshold != "bootstrap" && threshold != "manual") {
    stop("Threshold type parameter 'threshold' must be set to be either 'bootstrap' or 'manual'.")
  }
  if ((threshold != "bootstrap") && (reps != 199)) {
    warning("reps is only used with threshold=bootstrap")
  }
  if (is.null(threshold.val) && threshold == "manual") {
    stop("Threshold type has been set to 'manual', but threshold.val has not been set.")
  }

  if (!is.numeric(reps)) {
    stop("Number of bootstrap replications should be a single positive integer.")
  }
  if ((length(reps) != 1) || (reps %% 1 != 0) || (reps < 1)) {
    stop("Number of bootstrap replications should be a single positive integer.")
  }

  if (kernel.f != "gauss" && kernel.f != "euclidean" && kernel.f != "laplace" && kernel.f != "sine" && kernel.f != "quad.exp") {
    stop("The kernel.f function must be either 'quad.exp', 'gauss', 'laplace', 'sine', or 'euclidean'")
  }
  if ((kernel.f == "euclidean") && (kern.par <= 0 || kern.par >= 2)) {
    stop("For the 'euclidean' kernel function, the kernel parameter must be in the interval (0,2)")
  }

  if (kern.par < 0) {
    stop("The kernel parameter must be a positive value.")
  }

  if (boot.dep < 0) {
    stop("The bootstrap dependence parameter 'boot.dep' must be a positive value.")
  }
  if (boot.method != "mean.subtract" && boot.method != "no.mean.subtract") {
    stop("Parameter 'boot.method' must be set to be either 'mean.subtract' or 'no.mean.subtract'. Highly
         recommended to set boot.method = 'mean.subtract'.")
  }

}
