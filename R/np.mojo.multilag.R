#' @title Nonparametric Multiple Lag Change Point Detection for Multivariate Time Series
#' @description For a given set of lagged values of the time series, performs nonparametric change point detection of a possibly multivariate
#' time series.
#' @details See McGonigle and Cho (2023) for further details.
#' @param x Input data (a \code{numeric} vector or an object of classes \code{ts} and \code{timeSeries},
#' or a \code{numeric} matrix with rows representing variables)
#' @param G An integer value for the moving sum bandwidth;
#' \code{G} should be less than \code{length(n)/2}.
#' @param lags A \code{numeric} vector giving the range of lagged values of the time series that will be used to detect changes. See
#' \link{np.mojo} for further details.
#' @param kernel.f String indicating which kernel function to use when calculating the NP-MOJO statistic; possible values are
#'  \itemize{
#'    \item{\code{"quad.exp"}}{: kernel \eqn{h_2} in McGonigle and Cho (2023), kernel 5 in Fan et al. (2017):
#'    \deqn{h (x,y) = \prod_{i=1}^{2p} \frac{ (2a - (x_i - y_i)^2) \exp (-\frac{1}{4a} (x_i - y_i)^2 )}{2a} .}  }
#'    \item{\code{"gauss"}}{: kernel \eqn{h_1} in McGonigle and Cho (2023), the standard Gaussian kernel:
#'    \deqn{h (x,y) = \exp ( - \frac{a^2}{2} \Vert x - y  \Vert^2) .} }
#'    \item{\code{"euclidean"}}{: kernel \eqn{h_3} in McGonigle and Cho (2023), the Euclidean distance-based kernel:
#'    \deqn{h (x, y ) = \Vert x - y \Vert^a  .}}
#'    \item{\code{"laplace"}}{: kernel 2 in Fan et al. (2017), based on a Laplace weight function:
#'      \deqn{h (x, y ) = \prod_{i=1}^{2p} \left( 1+ a^2 (x_i - y_i)^2  \right)^{-1}. }}
#'    \item{\code{"sine"}}{: kernel 4 in Fan et al. (2017), based on a sinusoidal weight function:
#'      \deqn{h (x, y ) = \prod_{i=1}^{2p} \frac{-2 | x_i - y_i |  + | x_i - y_i - 2a|  + | x_i - y_i +2a| }{4a} .}}
#' }
#' @param kern.par The tuning parameter that appears in the expression for the kernel function, which acts as a scaling parameter.
#' @param data.driven.kern.par A \code{logical} variable, if set to \code{TRUE}, then the kernel tuning parameter is calculated
#'  using the median heuristic, if \code{FALSE} it is given by \code{kern.par}.
#' @param alpha a numeric value for the significance level with
#' \code{0 <= alpha <= 1}; use iff \code{threshold = "bootstrap"}
#' @param reps An integer value for the number of bootstrap replications performed, if \code{threshold = "bootstrap"}.
#' @param boot.dep A positive value for the strength of dependence in the multiplier bootstrap sequence, if \code{threshold = "bootstrap"}.
#' @param parallel A \code{logical} variable, if set to \code{TRUE}, then computation is performed in parallel is used if bootstrapping is performed,
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
#'    and \code{epsilon}-criterion is satisfied}
#' }
#' @param eta A positive numeric value for the minimal mutual distance of
#' changes, relative to bandwidth (if \code{criterion = "eta"}).
#' @param epsilon a numeric value in (0,1] for the minimal size of exceeding
#' environments, relative to moving sum bandwidth (iff \code{criterion = "epsilon"})
#' @param use.mean \code{Logical variable}, only to be used if \code{data.drive.kern.par=TRUE}. If set to \code{TRUE}, the mean
#' of pairwise distances is used to set the kernel function tuning parameter, instead of the median. May be useful for binary data,
#' not recommended to be used otherwise.
#' @param threshold String indicating how th
#' @param eta.merge A positive numeric value for the minimal mutual distance of
#' changes, relative to bandwidth, used to merge change point estimators across different lags.
#' @param merge.type String indicating the method used to merge change point estimators from different lags. Possible choices are
#'  \itemize{
#'    \item{\code{"sequential"}}{:  Starting from the left-most change point estimator and proceeding forward in time, estimators
#'    are grouped into clusters based on mutual distance. The estimator yielding the smallest corresponding p-value is
#'    chosen as the change point estimator for that cluster. See McGonigle and Cho (2023) for details.}
#'        \item{\code{"bottom-up"}}{: starting with the smallest p-value, the change points are merged using bottom-up merging (Messer
#'        et al. (2014)).}
#' }
#' @param threshold String indicating how the threshold is computed. Possible values are
#'  \itemize{
#'    \item{\code{"bootstrap"}}{: the threshold is calculated using the bootstrap method with significance level \code{alpha}.}
#'        \item{\code{"manual"}}{: the threshold is set by the user and must be specified using the \code{threshold.val} parameter.}
#' }
#' @param threshold.val The value of the threshold used to delcare change points, only to be used if \code{threshold = "manual"}.
#'
#' @return A \code{list} object that contains the following fields:
#'    \item{G}{Moving window bandwidth}
#'    \item{lags}{Lags used to detect changes}
#'    \item{kernel.f, data.driven.kern.par, use.mean}{Input parameters}
#'    \item{threshold, alpha, reps, boot.dep, boot.method, parallel}{Input parameters}
#'    \item{criterion, eta, epsilon}{Input parameters}
#'    \item{merged.cpts}{A \code{list} object which contains the following fields:
#'    \itemize{
#'    \item{cpts}{: a matrix with rows corresponding to final change point estimators, with associate lag and p-value given in columns}
#'    \item{cpt.clusters}{: a \code{list} object of length given by the number of detected change points. Each field contains a matrix of all
#'    change point estimators that are declared to be associated to the corresponding change point in the \code{cpts} field.}}
#'    }
#' @references McGonigle, E.T., Cho, H. (2023). Nonparametric data segmentation in multivariate time series via joint characteristic functions. \emph{arXiv preprint \href{https://doi.org/10.48550/arXiv.2305.07581}{arXiv:2305.07581}.}
#' @references Fan, Y., de Micheaux, P.L., Penev, S. and Salopek, D. (2017). Multivariate nonparametric test of independence. \emph{Journal of Multivariate Analysis},
#' 153, pp.189-210.
#' @references Messer M., Kirchner M., Schiemann J., Roeper J., Neininger R., Schneider G. (2014). A Multiple Filter Test for
#' the Detection of Rate Changes in Renewal Processes with Varying Variance. \emph{The Annals of Applied Statistics}, 8(4), 2027-2067.
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' n <- 1000
#' noise <- c(rep(1,300),rep(0.4,700))*stats::arima.sim(model = list(ar =0.3), n = 1000)
#' signal <- c(rep(0,700),rep(0.5,300))
#' x <- signal + noise
#' x.c <- np.mojo.multilag(x, G = 166, lags = c(0,1,2))
#' x.c$merged.cpts
#' }
#'
#' @seealso \link{np.mojo}, \link{multilag.cpts.merge}
np.mojo.multilag <- function(x, G, lags = c(0), kernel.f = c("quad.exp","gauss", "euclidean", "laplace","sine")[1],
                             kern.par = 1, data.driven.kern.par = TRUE, threshold = c("bootstrap", "manual")[1], threshold.val = NULL,
                             alpha = 0.1,  reps = 199, boot.dep = 1.5*(dim(as.matrix(x))[1]^(1/3)), parallel = FALSE, boot.method = 1,
                              criterion = "eta", eta = 0.4, epsilon = 0.02, use.mean = FALSE, eta.merge = 1,
                             merge.type = c("sequential", "bottom-up")[1]){

  if (!is.numeric(lags)) {
    stop("The set of lags must be a numeric vector.")
  }

  lag.cpts <- vector(mode = "list", length = length(lags))

  cpts <- init.cpts <- matrix(NA, nrow = 0, ncol = 3)
  dimnames(init.cpts)[[2]] <- c('cp', 'lag', 'p.val')
  dimnames(cpts)[[2]] <- c('cp', 'lag', 'p.val')

  for(l in 1:length(lag.cpts)){

    lag.cpts[[l]] <- np.mojo(x = x, G = G, lag = lags[l], kernel.f = kernel.f,  data.driven.kern.par = data.driven.kern.par,
                                       alpha = alpha, kern.par = kern.par, reps = reps, boot.dep = boot.dep, parallel = parallel,
                                       boot.method = boot.method, criterion = criterion, eta=eta, epsilon = epsilon, use.mean = use.mean, threshold = threshold,
                               threshold.val = threshold.val)

    if(length(lag.cpts[[l]]$cpts)>0){
      new.cpts <- cbind(lag.cpts[[l]]$cpts, rep(lags[l], length(lag.cpts[[l]]$cpts)), lag.cpts[[l]]$p.vals)
      init.cpts <- rbind(init.cpts, new.cpts)
    }

  }

  merged.cpts <- multilag.cpts.merge(lag.cpts,  eta.merge = eta.merge, merge.type = merge.type)

  ret <- list(G = G,
              lags = lags,
              kernel.f = kernel.f,
              data.driven.kern.par = data.driven.kern.par,
              merged.cpts = merged.cpts,
              threshold = threshold,
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


