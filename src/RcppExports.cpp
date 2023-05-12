// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// mosum_dist_calc
NumericMatrix mosum_dist_calc(Rcpp::NumericMatrix y, unsigned int G, unsigned int n, String kern, double kern_par);
RcppExport SEXP _CptNonPar_mosum_dist_calc(SEXP ySEXP, SEXP GSEXP, SEXP nSEXP, SEXP kernSEXP, SEXP kern_parSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type y(ySEXP);
    Rcpp::traits::input_parameter< unsigned int >::type G(GSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type n(nSEXP);
    Rcpp::traits::input_parameter< String >::type kern(kernSEXP);
    Rcpp::traits::input_parameter< double >::type kern_par(kern_parSEXP);
    rcpp_result_gen = Rcpp::wrap(mosum_dist_calc(y, G, n, kern, kern_par));
    return rcpp_result_gen;
END_RCPP
}
// rolling_matrix_sum
/* function to compute rolling matrix sums needed to compute the test statistic. */  NumericVector rolling_matrix_sum(NumericMatrix stat_mat, int G, int lag, double init_val, int n);
RcppExport SEXP _CptNonPar_rolling_matrix_sum(SEXP stat_matSEXP, SEXP GSEXP, SEXP lagSEXP, SEXP init_valSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type stat_mat(stat_matSEXP);
    Rcpp::traits::input_parameter< int >::type G(GSEXP);
    Rcpp::traits::input_parameter< int >::type lag(lagSEXP);
    Rcpp::traits::input_parameter< double >::type init_val(init_valSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(rolling_matrix_sum(stat_mat, G, lag, init_val, n));
    return rcpp_result_gen;
END_RCPP
}
// mojo_eta_criterion_help
IntegerVector mojo_eta_criterion_help(const IntegerVector& candidates, const NumericVector& m_values, double eta, double G_left, double G_right);
RcppExport SEXP _CptNonPar_mojo_eta_criterion_help(SEXP candidatesSEXP, SEXP m_valuesSEXP, SEXP etaSEXP, SEXP G_leftSEXP, SEXP G_rightSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector& >::type candidates(candidatesSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type m_values(m_valuesSEXP);
    Rcpp::traits::input_parameter< double >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< double >::type G_left(G_leftSEXP);
    Rcpp::traits::input_parameter< double >::type G_right(G_rightSEXP);
    rcpp_result_gen = Rcpp::wrap(mojo_eta_criterion_help(candidates, m_values, eta, G_left, G_right));
    return rcpp_result_gen;
END_RCPP
}
// C_matvecprod_elwise_inplace
NumericMatrix C_matvecprod_elwise_inplace(NumericMatrix& X, const NumericVector& y);
RcppExport SEXP _CptNonPar_C_matvecprod_elwise_inplace(SEXP XSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(C_matvecprod_elwise_inplace(X, y));
    return rcpp_result_gen;
END_RCPP
}
// C_matvecprod_elwise
SEXP C_matvecprod_elwise(const NumericMatrix& X, const NumericVector& y);
RcppExport SEXP _CptNonPar_C_matvecprod_elwise(SEXP XSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(C_matvecprod_elwise(X, y));
    return rcpp_result_gen;
END_RCPP
}
// C_matmatprod_elwise
SEXP C_matmatprod_elwise(const NumericMatrix& X, const NumericMatrix& Y);
RcppExport SEXP _CptNonPar_C_matmatprod_elwise(SEXP XSEXP, SEXP YSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Y(YSEXP);
    rcpp_result_gen = Rcpp::wrap(C_matmatprod_elwise(X, Y));
    return rcpp_result_gen;
END_RCPP
}
// C_matmatprod_elwise_inplace
NumericMatrix C_matmatprod_elwise_inplace(NumericMatrix& X, const NumericMatrix& Y);
RcppExport SEXP _CptNonPar_C_matmatprod_elwise_inplace(SEXP XSEXP, SEXP YSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Y(YSEXP);
    rcpp_result_gen = Rcpp::wrap(C_matmatprod_elwise_inplace(X, Y));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_CptNonPar_mosum_dist_calc", (DL_FUNC) &_CptNonPar_mosum_dist_calc, 5},
    {"_CptNonPar_rolling_matrix_sum", (DL_FUNC) &_CptNonPar_rolling_matrix_sum, 5},
    {"_CptNonPar_mojo_eta_criterion_help", (DL_FUNC) &_CptNonPar_mojo_eta_criterion_help, 5},
    {"_CptNonPar_C_matvecprod_elwise_inplace", (DL_FUNC) &_CptNonPar_C_matvecprod_elwise_inplace, 2},
    {"_CptNonPar_C_matvecprod_elwise", (DL_FUNC) &_CptNonPar_C_matvecprod_elwise, 2},
    {"_CptNonPar_C_matmatprod_elwise", (DL_FUNC) &_CptNonPar_C_matmatprod_elwise, 2},
    {"_CptNonPar_C_matmatprod_elwise_inplace", (DL_FUNC) &_CptNonPar_C_matmatprod_elwise_inplace, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_CptNonPar(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
