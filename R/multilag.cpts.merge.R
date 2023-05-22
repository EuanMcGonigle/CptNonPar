#' @title Merging Change Point Estimators from Multiple Lags Into Final Set of Change Point Estimators
#' @description Merges change point estimators from different lagged values into a final set of overall change point estimators.
#' @details See McGonigle and Cho (2023) for further details.
#' @param x.c A \code{list} object, where each element of the list is the output of the \link{np.mojo} function computed at a
#' different lag.
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
#'
#' @return A \code{list} object which contains the following fields
#' \item{cpts}{A matrix with rows corresponding to final change point estimators, with associate lag and p-value given in columns}
#'    \item{cpt.clusters}{A \code{list} object of length given by the number of detected change points. Each field contains a matrix of all
#'    change point estimators that are declared to be associated to the corresponding change point in the \code{cpts} field.}
#' @references McGonigle, E.T., Cho, H. (2023). Nonparametric data segmentation in multivariate time series via joint characteristic functions. \emph{arXiv preprint \href{https://doi.org/10.48550/arXiv.2305.07581}{arXiv:2305.07581}.}
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' n <- 1000
#' noise <- c(rep(1,300),rep(0.4,700))*stats::arima.sim(model = list(ar =0.3), n = 1000)
#' signal <- c(rep(0,700),rep(0.5,300))
#' x <- signal + noise
#' x.c0 <- np.mojo(x, G = 166, lag = 0)
#' x.c1 <- np.mojo(x, G = 166, lag = 1)
#' x.c2 <- np.mojo(x, G = 166, lag = 2)
#' x.c <- multilag.cpts.merge(list(x.c0,x.c1,x.c2))
#' x.c
#' }
#' @seealso \link{np.mojo}, \link{np.mojo.multilag}
multilag.cpts.merge <- function(x.c, eta.merge = 1, merge.type = c("sequential", "bottom-up")[1]){


  cpts <- init.cpts <- matrix(NA, nrow = 0, ncol = 3)
  dimnames(init.cpts)[[2]] <- c('cp', 'lag', 'p.val')
  dimnames(cpts)[[2]] <- c('cp', 'lag', 'p.val')

  cpt.clusters <- list()

  L <- length(x.c)

  G <- x.c[[1]]$G

  for(l in 1:L){

    lag.cpts <- x.c[[l]]

    if(length(lag.cpts$cpts)>0){
      new.cpts <- cbind(lag.cpts$cpts, rep(lag.cpts$lag, length(lag.cpts$cpts)), lag.cpts$p.vals)
      init.cpts <- rbind(init.cpts, new.cpts)
    }

  }

  if(nrow(init.cpts)<=1){
    return(list(cpts = init.cpts, cpt.clusters = init.cpts))
  }

  if(merge.type == "sequential"){

    cpt.order <- order(init.cpts[,1])

    init.cpts <- init.cpts[cpt.order,]

    j <- 1

    while(nrow(init.cpts)!=0){

      if(nrow(init.cpts)==1){
        cpts <- rbind(cpts, init.cpts)
        cpt.clusters[[j]] <- init.cpts
      } else{
        cpt.left <- init.cpts[1,1]
        linked.cpts <- init.cpts[2:nrow(init.cpts),1]
        linked.cpts <- c(cpt.left,linked.cpts[linked.cpts-cpt.left < eta.merge*G])

        num.linked.cpts <- length(linked.cpts)

        linked.lags <- init.cpts[1:num.linked.cpts,2]
        linked.p.vals <- init.cpts[1:num.linked.cpts,3]

        final.cpt <- linked.cpts[which.min(linked.p.vals)]
        final.lag <- linked.lags[which.min(linked.p.vals)]
        final.pval <- linked.p.vals[which.min(linked.p.vals)]

        cpts <- rbind(cpts, c(final.cpt, final.lag, final.pval))

        cpt.clusters[[j]] <- init.cpts[1:num.linked.cpts, , drop = FALSE][order(init.cpts[1:num.linked.cpts,2]), , drop = FALSE]
      }
      init.cpts <- init.cpts[-(1:num.linked.cpts),, drop = FALSE]
      j <- j + 1
    }

  }else{

    cpt.order <- order(init.cpts[,3])

    init.cpts <- init.cpts[cpt.order,]

    j <- 1


    while(nrow(init.cpts)!=0){

      if(nrow(init.cpts)==1){
        cpts <- rbind(cpts, init.cpts)
        cpt.clusters[[j]] <- init.cpts
      } else{

        cpt.min <- init.cpts[1,1]

        nearby.cpts <- which(abs(init.cpts[,1]-cpt.min) < eta.merge*G)
        linked.cpts <- init.cpts[nearby.cpts,, drop = FALSE]

        cpts <- rbind(cpts,  init.cpts[1,])
        cpt.clusters[[j]] <- linked.cpts[order(linked.cpts[,2]),]
      }
      init.cpts <- init.cpts[-nearby.cpts,, drop = FALSE]
      j <- j + 1
    }

    cpts <- cpts[order(cpts[,1]),]

  }

  return(list(cpts = cpts, cpt.clusters = cpt.clusters))

}

