// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
void do_it() {
  mat A = arma::randu<mat>(5,5);
  cout << A << endl;
  A.shed_row(2);
  cout << A << endl;
  A.shed_cols(2,4);
  cout << A << endl;
  return;
}
