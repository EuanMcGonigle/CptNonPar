# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' @keywords internal
mosum_dist_calc <- function(y, G, n, kern, kern_par) {
    .Call(`_CptNonPar_mosum_dist_calc`, y, G, n, kern, kern_par)
}

rolling_matrix_sum <- function(stat_mat, G, lag, init_val, n) {
    .Call(`_CptNonPar_rolling_matrix_sum`, stat_mat, G, lag, init_val, n)
}

#' @keywords internal
mojo_eta_criterion_help <- function(candidates, m_values, eta, G_left, G_right) {
    .Call(`_CptNonPar_mojo_eta_criterion_help`, candidates, m_values, eta, G_left, G_right)
}

C_matvecprod_elwise <- function(X, y) {
    .Call(`_CptNonPar_C_matvecprod_elwise`, X, y)
}

C_matmatprod_elwise_inplace <- function(X, Y) {
    .Call(`_CptNonPar_C_matmatprod_elwise_inplace`, X, Y)
}

