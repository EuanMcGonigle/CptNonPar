#' @keywords internal

runmean <- function(x, n = 5){stats::filter(x, rep(1 / n, n), sides = 1)}

#' @keywords internal

bootstrap.tstat.faster = function(h.mat, test.stat, G, lag, boot.dep, boot.method, data.len, h.mat.ncol) {

  #set.seed(1)

  Wtstar <- stats::arima.sim(model = list(ar = exp(-1/boot.dep)), sd = sqrt(1-exp(-2/boot.dep)),
                      n = h.mat.ncol, n.start = 1, start.innov = stats::rnorm(1))

  t.stat = rep(0, h.mat.ncol)

  if(boot.method==1){

    t.stat <- rep(0, h.mat.ncol+G+lag)

    W.mean = runmean(Wtstar, G-lag)

    Wtstar.out.a <- Rfast::Outer(Wtstar,Wtstar)

    h.mat.a <- C_matmatprod_elwise_inplace(Wtstar.out.a, h.mat) #FIRST ARGUMENT GETS OVER-WRITTEN

    h.mat.b <- C_matvecprod_elwise(h.mat, Wtstar)

    init.val.a <- sum(h.mat.a[1:(G-lag),1:(G-lag)])

    init.val.b <- sum(h.mat.b[1:(G-lag),1:(G-lag)])

    t.stat.a = rolling_matrix_sum(stat_mat = h.mat.a, G = G, lag = lag, init_val = init.val.a, n = h.mat.ncol)
    t.stat.b = rolling_matrix_sum(stat_mat = h.mat.b, G = G, lag = lag, init_val = init.val.b, n = h.mat.ncol)

    t.stat.a = c(rep(0,G-1),t.stat.a[1:(data.len-2*G+1)]/((G-lag)^2),rep(0,G))

    t.stat.b = c(rep(0,G-1),t.stat.b[1:(data.len-2*G+1)]/((G-lag)^2),rep(0,G))

    t.stat[1:(h.mat.ncol)] <- t.stat.a[1:(h.mat.ncol)]- 2*W.mean*t.stat.b[1:(h.mat.ncol)] +
      W.mean^2*test.stat[1:(h.mat.ncol)]


  }
  else{

    Wtstar.out = outer(Wtstar,Wtstar)

    h.mat2 = Wtstar.out*h.mat

    init.val = sum(h.mat2[1:(G-lag),1:(G-lag)])

    t.stat = rolling_matrix_sum(stat_mat = h.mat2,G = G, lag = lag, init_val = init.val, n = h.mat.ncol)

    t.stat = t.stat/((G-lag)^2)

  }

  #print(max(t.stat, na.rm=TRUE))

  max(t.stat, na.rm = TRUE)



}
