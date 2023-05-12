#' @title Nonparametric Single Lag Change Point Detection for Multivariate Time Series
#' @description For a given lagged value of the time series, performs nonparametric change point detection of a possibly multivariate
#' time series. If \code{lag = 0}, then only marginal changes are detected.
#' If \code{lag} \eqn{= \ell \neq 0}, then changes in the pairwise distribution of \eqn{(X_t , X_{t+\ell})} are detected.
#' @details See McGonigle and Cho (2023) for further details.
#' @param x Input data (a \code{numeric} vector or an object of classes \code{ts} and \code{timeSeries},
#' or a \code{numeric} matrix with rows representing variables)
#' @param G An integer value for the moving sum bandwidth;
#' \code{G} should be less than \code{length(n)/2}.
#' @param lag The lagged values of the time series used to detect changes. If \code{lag = 0}, then only marginal changes are detected.
#' If \code{lag} \eqn{= \ell \neq 0}, then changes in the pairwise distribution of \eqn{(X_t , X_{t+\ell})} are detected.
#' @param kernel.f String indicating which kernel function to use when calculating the NP-MOJO statistic; possible values are
#'  \itemize{
#'    \item{\code{"quad.exp"}}{: kernel \eqn{h_2} in McGonigle and Cho (2023), kernel 5 in Fan et al. (2017):
#'    \deqn{h (x,y) = \prod_{i=1}^{2p} \frac{ (2a - (x_i - y_i)^2) \exp (-\frac{1}{4a} (x_i - y_i)^2 )}{2a} .}  }
#'    \item{\code{"gauss"}}{: kernel \eqn{h_1} in McGonigle and Cho (2023), the standard Gaussian kernel:
#'    \deqn{h (x,y) = \exp ( - \frac{a^2}{2} \Vert x - y  \Vert^2) .} }
#'    \item{\code{"euclidean"}}{: kernel \eqn{h_3} in McGonigle and Cho (2023), the Euclidean distance-based kernel:
#'    \deqn{h (x, y ) = \Vert x - y \Vert^a  .}}
#'    \item{\code{"laplace"}}{: kernel 2 in Fan et al. (2017).}
#'    \item{\code{"sine"}}{: kernel 4 in Fan et al. (2017).}
#' }
#' @param kern.par The tuning parameter that appears in the expression for the kernel function, which acts as a scaling parameter,
#' only to be used if \code{data.driven.kern.par = FALSE}. If \code{kernel.f = "euclidean"}, then \code{kern.par} \eqn{\in (0,2)},
#' otherwise \code{kern.par} \eqn{> 0}.
#' @param data.driven.kern.par A \code{logical} variable, if set to \code{TRUE}, then the kernel tuning parameter is calculated
#'  using the median heuristic, if \code{FALSE} it is given by \code{kern.par}.
#' @param alpha A numeric value for the significance level with
#' \code{0 <= alpha <= 1}; use iff \code{threshold = "bootstrap"}
#' @param reps An integer value for the number of bootstrap replications performed, if \code{threshold = "bootstrap"}.
#' @param boot.dep A positive value for the strength of dependence in the multiplier bootstrap sequence, if \code{threshold = "bootstrap"}.
#' @param parallel A \code{logical} variable, if set to \code{TRUE}, then parallelisation is used if bootstrapping is performed,
#'  if \code{FALSE} no parallelisation is performed.
#' @param boot.method A string indicating the method for creating bootstrap replications. It is not recommended to change this. Possible choices are
#' #'  \itemize{
#'    \item{\code{"1"}}{: empirical mean subtraction is performed to the bootstrapped replicates, improving power.}
#'        \item{\code{"2"}}{: empirical mean subtraction is not performed, improving size control.}
#' }
#' @param criterion String indicating how to determine whether each point \code{k} at which NP-MOJO statistic
#' exceeds the threshold is a change point; possible values are
#'  \itemize{
#'    \item{\code{"epsilon"}}{\code{k} is the maximum of its local exceeding environment, which has at least size \code{epsilon*G}}
#'        \item{\code{"eta.and.ep"}}{there is no larger exceeding in an \code{eta*G} environment of \code{k},
#'    and \code{epsilon}-criteria is satisfied}
#' }
#' @param eta A positive numeric value for the minimal mutual distance of
#' changes, relative to bandwidth (if \code{criterion = "eta"}).
#' @param epsilon a numeric value in (0,1] for the minimal size of exceeding
#' environments, relative to moving sum bandwidth (iff \code{criterion = "epsilon"})
#' @param use.mean \code{Logical variable}, only to be used if \code{data.drive.kern.par=TRUE}. If set to \code{TRUE}, the mean
#' of pairwise distances is used to set the kernel function tuning parameter, instead of the median. May be useful for binary data,
#' not recommended to be used otherwise.
#' @param threshold String indicating how the threshold is computed. Possible values are
#'  \itemize{
#'    \item{\code{"bootstrap"}}{: the threshold is calculated using the bootstrap method with significance level \code{alpha}.}
#'        \item{\code{"manual"}}{: the threshold is set by the user and must be specified using the \code{threshold.val} parameter.}
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
#'    \item{test.stat}{A vector containing the NP-MOSTAT detector statistics computed from the input data}
#'    \item{cpts}{A vector containing the estimated change point locations}
#'    \item{p.vals}{The corresponding p values of the change points, if the bootstrap method was used}
#' @references McGonigle, E.T., Cho, H. (2023). Nonparametric change point analysis for multivariate time series. \emph{arXiv preprint} \href{https://arxiv.org/abs/2108.07550}{arXiv:2108.07550}.
#' @references Fan, Y., de Micheaux, P.L., Penev, S. and Salopek, D. (2017). Multivariate nonparametric test of independence. \emph{Journal of Multivariate Analysis},
#' 153, pp.189-210.
#'
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' n <- 1000
#' noise <- c(rep(1,300),rep(0.4,700))*stats::arima.sim(model = list(ar =0.3), n = 1000)
#' signal <- c(rep(0,700),rep(0.5,300))
#' x <- signal + noise
#' x.c <- np.mojo(x, G = 166, lag = 0)
#' x.c$cpts
#' x.c$p.vals
#' }
#' @importFrom Rcpp evalCpp
#' @useDynLib CptNonPar, .registration = TRUE
#' @importFrom foreach %dopar%
#' @seealso \link{np.mojo.multilag}
np.mojo = function(x, G,lag = 0, kernel.f = c("quad.exp", "gauss", "euclidean","laplace","sine")[1],
                      kern.par = 1, data.driven.kern.par = TRUE, alpha = 0.1, threshold = c("bootstrap", "manual")[1],
                      threshold.val = NULL, reps = 199,  boot.dep = 1.5*(dim(as.matrix(x))[1]^(1/3)), parallel = FALSE,
                      boot.method = 1, criterion = "eta", eta = 0.4, epsilon = 0.02, use.mean = FALSE){



  #A faster implementation of cpt.nonpar using calls to Rcpp for matrix and vector operations in the bootstrap procedure

  #Finds nonparametric changes using the integrated squared distance between
  #joint characteristic function to left and right of change within a MOSUM procedure

  #Error handling:

  stopifnot(criterion == "epsilon" || criterion == "eta")
  stopifnot(criterion != "epsilon" || epsilon >= 0)
  stopifnot(criterion != "eta" || eta >= 0)

  if (!is.numeric(x)) {
    stop("Data must be numeric")
  }
  if (any(is.na(x))) {
    stop("Missing values in data: NA is not allowed in the data.")
  }


  if (!is.numeric(reps)) {
    stop("Number of bootstrap replications should be a single positive integer")
  }
  if ((length(reps) != 1) | (reps%%1 != 0) | (reps < 1)) {
    stop("Number of bootstrap replications should be a single positive integer")
  }
  if(kernel.f != "gauss" & kernel.f != "euclidean" & kernel.f != "laplace" & kernel.f != "sine" & kernel.f != "quad.exp"){
    stop("The kernel.f function must be either 'quad.exp', 'gauss', 'laplace', 'sine', or 'euclidean'")
  }
  if ((threshold != "bootstrap") & (reps != 199)) {
    warning("reps is only used with threshold=bootstrap")
  }
  if (kern.par<0) {
    warning("The kernel parameter must be a positive value.")
  }
  if((kernel.f == "euclidean") & (kern.par <=0 || kern.par>=2)){
    stop("For the 'euclidean' kernel function, the kernel parameter must be in the interval (0,2)")
  }

  if(is.null(dim(x))){
    data.len <- length(x)
    x <- matrix(x)
  }else{
    data.len <- dim(x)[1]
  }

  if (G>data.len/2) {
    stop("Bandwidth is too large for the length of time series.")
  }

  if(lag!=0){
    y <- cbind(x[1:(data.len-lag),], x[(lag+1):data.len,])
  } else{
    y <- x
  }

  p.vals <- numeric(0)

  #calculate the Euclidean distances required for the test statistic calculation:

  if(data.driven.kern.par == TRUE){
    if(kernel.f == "euclidean"){
      warning("Data driven parameter choice not suited for this kernel. Parameter kern.par reset to 1.")
      kern.par <- 1
    }
    else{
      D.par <- mosum_dist_calc(y = y, G = G, n = data.len-lag, kern = "euc.dist", kern_par = kern.par)
      if(use.mean == TRUE){
        delta <- mean(D.par[D.par>0])
      } else{
        delta <- stats::median(D.par[D.par>0])
      }
      if(kernel.f == "gauss"){
        kern.par = 1/sqrt(delta)
      } else if (kernel.f == "quad.exp"){
        kern.par = 0.5*delta
      }
    }
  }


  if(kernel.f == "gauss" & data.driven.kern.par == TRUE){
    D <- D.par
  } else{
    D <- mosum_dist_calc(y = y, G = G, n = data.len-lag, kern = kernel.f, kern_par = kern.par)

    if(kernel.f != "euclidean" && kernel.f != "gauss"){
      diag(D) <- 1
    }
  }



  h.matcalc = function(D, lag, G, kernel.f, kern.par){

    #calculates the entries of the h kernel matrix required to compute the MOSUM test statistic

    n <- dim(D)[2]

    h.dim <- n-G

    if(kernel.f == "gauss"){
      h.mat <- exp(-kern.par^2*(D[1:h.dim,1:h.dim])/2)+exp(-kern.par^2*(D[(G+1):(n),(G+1):(n)])/2)-
        exp(-kern.par^2*(D[1:h.dim, (G+1):(n)])/2) - exp(-kern.par^2*(D[(G+1):(n), 1:h.dim])/2)
    }else if(kernel.f == "euclidean"){
      h.mat <- (D[1:h.dim, (G+1):(n)])+ (D[(G+1):n, 1:h.dim]) -(D[1:h.dim,1:h.dim])-
        (D[(G+1):(n),(G+1):(n)])
    }else{
      h.mat <- -(D[1:h.dim, (G+1):(n)]) - D[(G+1):n, 1:h.dim] +
        (D[1:h.dim,1:h.dim]) +(D[(G+1):(n),(G+1):(n)])
    }

    h.mat
  }

  #calculate the entries of the kernel statistic:

  h.mat <- h.matcalc(D = D, lag = lag, G = G, kernel.f = kernel.f, kern.par = kern.par)

  h.mat.ncol <- ncol(h.mat)

  init.val = sum(h.mat[1:(G-lag),1:(G-lag)])

  #calculate the MOSUM test statistic to find change points:

  test.stat = rolling_matrix_sum(stat_mat = h.mat, G = G, lag = lag, init_val = init.val, n = h.mat.ncol)

  test.stat = c(rep(0,G-1),test.stat[1:(data.len-2*G+1)]/((G-lag)^2),rep(0,G))

  bootstrap.char = function(h.mat, test.stat, G,lag, reps, boot.dep, parallel, boot.method, data.len, h.mat.ncol){

    if (parallel==FALSE){
      d = replicate(reps, bootstrap.tstat.faster(h.mat = h.mat, test.stat = test.stat, G = G, lag = lag, boot.dep = boot.dep,
                                          boot.method = boot.method, data.len = data.len, h.mat.ncol = h.mat.ncol))
    }
    else{
      d = bootstrap.tstat.faster(h.mat, test.stat, G, lag, boot.dep, boot.method, data.len = data.len, h.mat.ncol)
    }

    d

  }

  #perform wild multiplier bootstrap to assess significance of changes:

  if(threshold=="bootstrap"){

      if (parallel == TRUE) {
        closeAllConnections()
        cl <- parallel::makeCluster(parallel::detectCores())
        doParallel::registerDoParallel(cl)
        parallel::clusterSetRNGStream(cl = cl, iseed = 9182)

        Tstar <- foreach::foreach(iterators::icount(reps), .combine = cbind, .packages = c("foreach", "Rcpp", "CptNonPar")) %dopar% bootstrap.char(h.mat,test.stat,G,lag,reps=1,boot.dep,parallel,boot.method,data.len,h.mat.ncol)

        parallel::stopCluster(cl)
      }
      else{
        Tstar <- bootstrap.char(h.mat = h.mat, test.stat = test.stat, G = G, lag = lag, reps = reps, boot.dep = boot.dep,
                                parallel = parallel, boot.method = boot.method, data.len = data.len, h.mat.ncol = h.mat.ncol)
      }

      threshold.val <- stats::quantile(Tstar, probs = 1-alpha)

  }

  cpt.locs <- numeric(0)

  if(max(test.stat)>threshold.val){

    exceedings <- (test.stat>threshold.val)

    if (criterion == "epsilon") {
      exceedingsCount <- (exceedings) * unlist(lapply(rle(exceedings)$lengths,
                                                      seq_len))
      minIntervalSize <- max(1, G * epsilon)
      intervalEndPoints <- which(diff(exceedingsCount) <= -minIntervalSize)
      intervalBeginPoints <- intervalEndPoints - exceedingsCount[intervalEndPoints] + 1

      if (exceedings[data.len - G] && !((data.len - G) %in% intervalEndPoints)) {

        lastBeginPoint <- data.len - G - exceedingsCount[data.len - G] + 1
        stopifnot(exceedings[seq(lastBeginPoint, data.len - G)])
        stopifnot(!(lastBeginPoint %in% intervalBeginPoints))
        highestStatPoint <- which.max(test.stat[seq(lastBeginPoint, data.len - G)]) + lastBeginPoint - 1

        if (highestStatPoint - lastBeginPoint >= minIntervalSize/2) {
          intervalEndPoints <- c(intervalEndPoints, data.len - G)
          intervalBeginPoints <- c(intervalBeginPoints,lastBeginPoint)
        }
      }
      if (exceedings[G] && !(G %in% intervalBeginPoints)) {

        firstEndPoint <- which(diff(exceedingsCount) < 0)[1]
        stopifnot(exceedings[seq(G, firstEndPoint)])
        stopifnot(!(firstEndPoint %in% intervalEndPoints))
        highestStatPoint <- which.max(test.stat[seq(G, firstEndPoint)]) + G - 1

        if (firstEndPoint - highestStatPoint >= minIntervalSize/2) {
          intervalEndPoints <- c(firstEndPoint, intervalEndPoints)
          intervalBeginPoints <- c(G, intervalBeginPoints)
        }
      }

      numChangePoints <- length(intervalBeginPoints)
      if (numChangePoints > 0) {
        for (i in 1:numChangePoints) {
          changePoint <- intervalBeginPoints[i] + which.max(test.stat[(intervalBeginPoints[i]):(intervalEndPoints[i])]) -
            1
          cpt.locs <- c(cpt.locs, changePoint)
        }
      }
    }
    else {
      stat.exceed <- which(test.stat > threshold.val)
      r = rle(diff(stat.exceed))
      end.r = cumsum(r$lengths)
      start.r = end.r - r$lengths + 1

      epsilon.satisfied.start <- stat.exceed[start.r[r$lengths>epsilon*G]]
      epsilon.satisfied.end <- stat.exceed[end.r[r$lengths>epsilon*G]]

      epsilon.exceedings <- rep(FALSE, data.len)
      for(i in seq_len(length(epsilon.satisfied.start))){
        epsilon.exceedings[epsilon.satisfied.start[i]:epsilon.satisfied.end[i]] <- TRUE
      }

      localMaxima <- (c((diff.default(test.stat) < 0), NA) & c(NA, diff.default(test.stat) > 0))
      localMaxima[data.len - G] <- TRUE
      p.candidates <- which(epsilon.exceedings & localMaxima)
      cpt.locs <- mojo_eta_criterion_help(p.candidates, test.stat, eta, G, G)
    }

    if(threshold == "bootstrap"){
      p.vals <- numeric(0)
      if(length(cpt.locs)==0){
        p.vals <- sum(Tstar >= max(test.stat))/(reps + 1)
      }
      else{
        for(i in 1:length(cpt.locs)){
          p.vals <- c(p.vals, sum(Tstar >= test.stat[cpt.locs[i]])/(reps + 1))
        }
      }
    }
  }
  else{
    if(threshold == "bootstrap"){
      p.vals <- sum(Tstar >= max(test.stat))/(reps + 1)
    }
  }


  ret <- list(x = x,
              G = G,
              lag = lag,
              kernel.f = kernel.f,
              kern.par = kern.par,
              data.driven.kern.par = data.driven.kern.par,
              test.stat = test.stat,
              cpts = cpt.locs,
              p.vals = p.vals,
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
              use.mean = use.mean)

  return(ret)

}


