// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

double kreg(mat Y1,
			mat X1,
			rowvec x1,
			double b1) {

	int n_obs = Y1.n_rows;
	int n_features = X1.n_cols;

	mat xx1 = repmat(x1, n_obs, 1);
	mat psi = (X1 - xx1)/b1;
	mat K = prod(1.0/std::sqrt(2*PI)*arma::exp(-arma::pow(psi,2)/2),1);

	return accu(Y1 % K)/accu(K);
}


// [[Rcpp::export]]
NumericVector nw_loo(NumericMatrix Y,
			         NumericMatrix X,
			         double b1) {

    mat Y1 = as<mat>(Y);
	mat X1 = as<mat>(X);

	int n_obs = X1.n_rows;
	int n_features = X1.n_cols;

	mat x = mat(1, n_features);
    mat y = mat(1, 1);

	colvec regression_results = colvec(n_obs);

	for (int i = 0; i < n_obs; i++) {
		x = X1.row(i);
        y = Y1.row(i);
        X1.shed_row(i);
        Y1.shed_row(i);

		regression_results(i) = kreg(Y1, X1, x, b1);
	    
        X1.insert_rows(i, x);
        Y1.insert_rows(i, y); 
    }

	return wrap(regression_results);
}



