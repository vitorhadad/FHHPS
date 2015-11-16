// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// G1-inverse matrix
mat g1inv(rowvec x) {
	mat g1 = mat(2,2, fill::none);
	g1(0,0) = 1; g1(0,1) = x(0);
	g1(1,0) = 1; g1(1,1) = x(1);
	mat g1inv = inv(g1);
	return(g1inv);
}


// G3 matrix
mat g3(rowvec x) {
    mat g2 = mat(2,2);
    g2(0,0) = 0; g2(0,1) = 0;
    g2(1,0) = 1; g2(1,1) = x(1);
    return g1inv(x) * g2;
}

// Kernel regression
cx_double kreg(cx_mat Y1,
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



// numerator_var[i] = exp(i * s * G1^{-1}(x[i]) * Y[i])
cx_mat numerator_var_creator(mat Y1,
			   		         rowvec x1,
			   		         mat s1)  {

	mat tmp = s1.t() * g1inv(x1) * Y1.t();

    // Using Euler's formula: exp(iz) = cos(z) + isin(z) 
    cx_mat numerator_var = cx_mat(arma::cos(tmp), arma::sin(tmp));
	return numerator_var.t();
}


// Characteristic function of numerator 
cx_colvec cf_numerator(mat y1,
			           mat x1,
			           mat s1,
			           double b1) {

	int n_obs = x1.n_rows;
	int n_features = x1.n_cols;

	mat x = mat(1, n_features);
	cx_mat numerator_var = cx_mat(n_obs, 1);
	cx_colvec results = cx_colvec(n_obs);


	for (int i = 0; i < n_obs; i++) {
		x = x1.row(i);
		numerator_var = numerator_var_creator(y1, x, s1);
		results(i) = kreg(numerator_var, x1, x, b1);
	}
	return results;
}


// increases the modulus to lower_bound
cx_colvec trim(cx_colvec v, double lower_bound) {

    for (int i = 0; i < v.n_elem; i++) {
        if (abs(v[i]) < lower_bound) {
            v[i] = v[i]/abs(v[i]) * lower_bound;
        }
    }
    return v;
} 


// Characteristic fct of denominator
cx_colvec cf_denominator(mat muW,
                        mat SigmaW,
			           mat X1,
			           mat s1,
			           double b1,
                       double t) {

	int n_obs = X1.n_rows;
	int n_features = X1.n_cols;

	mat x = mat(1, n_features);
	cx_colvec results = cx_colvec(n_obs);
    
	for (int i = 0; i < n_obs; i++) {
		x = X1.row(i);
        results(i) = as_scalar(exp(
                    1i * s1.t() * g3(x) * muW 
                    - .5 * s1.t() * g3(x) * SigmaW * g3(x).t() * s1));
    }
	return trim(results, t);
}



// Computes a matrix of characteristic function values over s_grid
//[[Rcpp::export]]
ComplexVector target_cf(NumericMatrix Y,
						NumericMatrix X,
					    NumericMatrix muW,
                        NumericMatrix SigmaW,
                        NumericVector s_grid,
						double b1,
                        double t)  {

	mat Y1 = as<mat>(Y);
	mat muW1 = as<mat>(muW);
    mat SigmaW1 = as<mat>(SigmaW);
	mat X1 = as<mat>(X);
	vec s_grid1 = as<vec>(s_grid);

    int s_pts = s_grid1.n_elem;
    mat s = zeros(2,1);

    cx_mat target_cf = cx_mat(s_pts, s_pts);
    
    for (int i = 0; i < s_pts; i++) {
        for (int j = 0; j < s_pts; j++) {
            s(0) = s_grid1(i); s(1) = s_grid1(j);
            target_cf(i,j) = arma::mean(arma::mean(cf_numerator(Y1, X1, s, b1)/
                                        cf_denominator(muW1, SigmaW1, X1, s, b1, t)));
        }
    }  
	return wrap(target_cf);
}



// [[Rcpp::export]]
ComplexVector kreg_R(ComplexMatrix Y,
			         NumericMatrix X,
			         NumericMatrix x,
			         double b)   {

	cx_mat Y1 = as<cx_mat>(Y);
	mat X1 = as<mat>(X);    
	rowvec x1 = as<rowvec>(x);
    
    cx_double output = kreg(Y1, X1, x1, b);
    
    return wrap(output);
}






