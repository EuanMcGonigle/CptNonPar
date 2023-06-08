#define restrict __restrict__ // gcc/clang
//#define restrict __restrict   // MS Visual Studio

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <cstdlib>
using std::size_t;

#include <R.h>
#include <Rdefines.h>


/* function to compute distance matrix needed to compute the test statistic. */
//' @keywords internal
// [[Rcpp::export]]

NumericMatrix mosum_dist_calc(Rcpp::NumericMatrix y, unsigned int G, unsigned int n, String kern, double kern_par) {
  Rcpp::NumericMatrix res(n);
  double d;
  unsigned int i = 0, j = 0;

  if (kern == "quad.exp"){
    double d2 = 2*kern_par;
    double d3 = 2*d2;
    for (i=0; i<n-1; ++i) {
      Rcpp::NumericVector v1 = y.row(i);
      for(j = i+1; j< std::min(i + 2*G, n); ++j){
        Rcpp::NumericVector v2 = -pow(v1-y.row(j), 2.0);
        Rcpp::NumericVector v3 = exp(v2/d3)*(d2 + v2)/d2;
        d = std::accumulate(v3.begin(), v3.end(), 1.0, std::multiplies<double>());
        res(i,j) = d;
        res(j,i) = d;
      }
    }
  } else if (kern == "sine"){
    double d2 = 2*kern_par;
    double d3 = 2*d2;
    for (i=0; i<n-1; ++i) {
    Rcpp::NumericVector v1 = y.row(i);
      for(j = i+1; j< std::min(i + 2*G, n); ++j){
        Rcpp::NumericVector v2 = (-2*abs(v1-y.row(j))+abs(v1-y.row(j) -d2)+abs(v1-y.row(j)+d2))/d3  ;
        d = std::accumulate(v2.begin(), v2.end(), 1.0, std::multiplies<double>());
        res(i,j) = d;
        res(j,i) = d;
      }
    }
  } else if (kern == "laplace"){
    double d2 = pow(kern_par,2.0);
    for (i=0; i<n-1; ++i) {
    Rcpp::NumericVector v1 = y.row(i);
      for(j = i+1; j< std::min(i + 2*G, n); ++j){
        Rcpp::NumericVector v2 = 1/(1+d2*pow(v1-y.row(j), 2.0));
        d = std::accumulate(v2.begin(), v2.end(), 1.0, std::multiplies<double>());
        res(i,j) = d;
        res(j,i) = d;
      }
    }
  }else if (kern == "gauss"){
    for (i=0; i<n-1; ++i) {
      Rcpp::NumericVector v1 = y.row(i);
      for(j = i+1; j< std::min(i + 2*G, n); ++j){
        d = sum(pow(v1-y.row(j), 2.0));
        res(i,j) = d;
        res(j,i) = d;
      }
    }
  } else if (kern == "euclidean"){
    for (i=0; i<n-1; ++i) {
      Rcpp::NumericVector v1 = y.row(i);
      for(j = i+1; j< std::min(i + 2*G, n); ++j){
        d = pow(sum(pow(v1-y.row(j), 2.0)), kern_par/2);
        res(i,j) = d;
        res(j,i) = d;
      }
    }
  } else if (kern == "euc.dist"){
    for (i=0; i<n-1; ++i) {
      Rcpp::NumericVector v1 = y.row(i);
      for(j = i+1; j< std::min(i + 2*G, n); ++j){
        d = sum(pow(v1-y.row(j), 2.0));
        res(i,j) = d;
        res(j,i) = d;
      }
    }
  }
  return res;
}


// [[Rcpp::export]]

/* function to compute rolling matrix sums needed to compute the test statistic. */

NumericVector rolling_matrix_sum(NumericMatrix stat_mat, int G, int lag, double init_val, int n) {
  NumericVector res(n, NA_REAL);
  double currentSum = init_val;
  res[0] = currentSum;
  for (int t=1; t<n-G+lag+1; ++t) {
      for(int j = t; j<t-1+G-lag; ++j){
        currentSum += stat_mat(j,t-1+G-lag) + stat_mat(t+G-lag-1,j) - stat_mat(t-1,j) - stat_mat(j-1,t-1);
        /* currentSum -= stat_mat(t-1,j) + stat_mat(j-1,t-1); */
      }
      currentSum += stat_mat(t-1+G-lag,t-1+G-lag) - stat_mat(t-2+G-lag,t-1);
      /* currentSum -= stat_mat(t-2+G-lag,t-1); */
    res[t] = currentSum;
  }
  return res;
}

/* extract changepoints from candidates with eta criterion  */
//' @keywords internal
// [[Rcpp::export]]
IntegerVector mojo_eta_criterion_help(const IntegerVector &candidates,
                                 const NumericVector &m_values,
                                 double eta, double G_left, double G_right) {
  const int n = m_values.length();
  IntegerVector res(0);
  const int left_length = std::floor(eta*G_left);
  const int right_length = std::floor(eta*G_right);
  for (int j=0; j<candidates.length(); ++j) {
    const int k_star = candidates[j];
    // Careful: k_star 1-indexed, vectors here are 0-indexed
    const double m_star = m_values[k_star-1];
    const int left_thresh = std::max(1, k_star-left_length) - 1; // see above
    const int right_thresh = std::min(n, k_star+right_length) - 1; // ""
    bool largest = true;
    for (int l=left_thresh; (l<=right_thresh) && largest; ++l) {
      if (m_values[l] > m_star) {
        largest = false;
      }
    }
    if (largest) {
      // k_star is maximum in its eta*G environment --> accept as changepoint
      res.push_back(k_star);
    }
  }
  return res;
}

/* functions to compute element-wise matrix times vector, faster than R base version, with no over-writing. */

void hadamardMultiplyMatrixByVector(const double* restrict x,
                                    size_t numRows, size_t numCols,
                                    const double* restrict y,
                                    double* restrict z)
{
  if (numRows == 0 || numCols == 0) return;

  for (size_t row = 0; row < numRows; ++row) {
    const double* restrict x_row = x + row * numCols;
    double* restrict z_row = z + row * numCols;

    for (size_t col = 0; col < numCols; ++col) {
      z_row[col] = x_row[col] * y[row];
    }
  }
}

// [[Rcpp::export]]
SEXP C_matvecprod_elwise(const NumericMatrix& X, const NumericVector& y)
{
  size_t numRows = X.nrow();
  size_t numCols = X.ncol();



  SEXP Z = PROTECT(Rf_allocVector(REALSXP, (int) (numRows * numCols)));
  SEXP dimsExpr = PROTECT(Rf_allocVector(INTSXP, 2));
  int* dims = INTEGER(dimsExpr);
  dims[0] = (int) numRows;
  dims[1] = (int) numCols;
  Rf_setAttrib(Z, R_DimSymbol, dimsExpr);

  hadamardMultiplyMatrixByVector(X.begin(), X.nrow(), X.ncol(), y.begin(), REAL(Z));

  UNPROTECT(2);

  return Z;
}


/* function to compute element-wise matrix times matrix, much faster than R base version, first argument is over-written. */

void hadamardMultiplyMatrixByMatrixInPlace(double* restrict X, size_t numRows, size_t numCols, const double* restrict Y) {
  if (numRows == 0 || numCols == 0) return;

  for (size_t row = 0; row < numRows; ++row) {
    double* restrict X_row = X + row * numCols;
    const double* restrict Y_row = Y + row * numCols;

    for (size_t col = 0; col < numCols; ++col) {
      X_row[col] *= Y_row[col];
    }
  }
}

// [[Rcpp::export]]

NumericMatrix C_matmatprod_elwise_inplace(NumericMatrix& X, const NumericMatrix& Y)
{


  hadamardMultiplyMatrixByMatrixInPlace(X.begin(), X.nrow(), X.ncol(), Y.begin());

  return X;
}
