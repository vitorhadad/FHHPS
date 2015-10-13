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



// Creates the dependent variable vector:
// numerator_var[i] = exp(i * s * G1^{-1}(x[i]) * Y[i])
cx_mat numerator_var_creator(mat Y1,
			   		   mat X1,
			   		   rowvec x1,
			   		   mat s1)  {


    mat g1inv_mat = g1inv(x1);
	mat tmp = s1.t() * g1inv_mat * Y1.t();

    // Using Euler's formula: exp(iz) = cos(z) + isin(z) 
    cx_mat numerator_var = cx_mat(arma::cos(tmp), arma::sin(tmp));
	return numerator_var.t();
}


// Characteristic function of G1^{-1}(x)'Y 
cx_colvec cf_numerator(mat Y1,
			           mat X1,
			           mat s1,
			           double b1) {

	int n_obs = X1.n_rows;
	int n_features = X1.n_cols;

	mat x = mat(1, n_features);
	cx_mat numerator_var = cx_mat(n_obs,1);
	cx_colvec results = cx_colvec(n_obs);


	for (int i = 0; i < n_obs; i++) {
		x = X1.row(i);
		numerator_var = numerator_var_creator(Y1, X1, x, s1);
		results(i) = kreg(numerator_var, X1, x, b1);
	}

	return results;
}


// Characteristic function of denominator variable
// denominator variable = i * s * G_3(x) * W
cx_colvec cf_denominator(mat W1,
			   		    mat X1,
			   		    mat s1) {

	int n_obs = X1.n_rows;
	cx_colvec results = cx_colvec(n_obs);
	cx_mat tmp = cx_mat(1, n_obs);
	rowvec x = rowvec(2);

	for (int i = 0; i < n_obs; i++) {
		x = X1.row(i);
	    // Euler's formula
        tmp.set_real(arma::cos(s1.t() * g3(x) * W1.t()));
		tmp.set_imag(arma::sin(s1.t() * g3(x) * W1.t()));
		results(i) = mean(mean(tmp, 1));
	}
	return results;
}


// [[Rcpp::export]]
ComplexVector target_cf(NumericMatrix Y,
						NumericMatrix W,
					    NumericMatrix X,
                        NumericVector s_grid,
						double b1)  {

	mat Y1 = as<mat>(Y);
	mat W1 = as<mat>(W);
	mat X1 = as<mat>(X);
	vec s_grid1 = as<vec>(s_grid);

    int s_pts = s_grid1.n_elem;
    mat s = zeros(2,1);

    cx_vec target_cf = cx_vec(s_pts*s_pts);
    
    int k = 0;
    for (int i = 0; i < s_pts; i++) {
        for (int j = 0; j < s_pts; j++) {
            s(0) = s_grid1(i); s(1) = s_grid1(j);
            target_cf[k++] = arma::mean(arma::mean(cf_numerator(Y1, X1, s, b1)/cf_denominator(W1, X1, s)));
        }
    }  

	return wrap(target_cf);
}



// Divides the char. fcts pointwise and takes mean
// [[Rcpp::export]]
ComplexVector rhs_variable(NumericMatrix Y,
						   NumericMatrix W,
					       NumericMatrix X,
						   NumericMatrix s,
						   double b1)  {

	mat Y1 = as<mat>(Y);
	mat W1 = as<mat>(W);
	mat X1 = as<mat>(X);
	mat s1 = as<mat>(s);

	cx_colvec numerator = cf_numerator(Y1, X1, s1, b1);
	cx_colvec denominator = cf_denominator(W1, X1, s1);

	return wrap(arma::mean(numerator/denominator));
}



// This is just a wrapper to call function from R
// [[Rcpp::export]]
ComplexVector cf_numerator_R(NumericMatrix Y,
        				     NumericMatrix X,
        					 NumericMatrix s,
                             double b1)  {

	mat Y1 = as<mat>(Y);
	mat X1 = as<mat>(X);    
	mat s1 = as<mat>(s);

	cx_colvec output = cf_numerator(Y1, X1, s1, b1);

  return wrap(output);
}


// This is just a wrapper to call function from R
// [[Rcpp::export]]
ComplexVector cf_denominator_R(NumericMatrix W,
                               NumericMatrix X,
        					   NumericMatrix s)  {

	mat W1 = as<mat>(W);
	mat X1 = as<mat>(X);    
	mat s1 = as<mat>(s);

	cx_colvec output = cf_denominator(W1, X1, s1);

  return wrap(output);
}





// [[Rcpp::export]]
ComplexVector kreg_R(NumericMatrix Y,
			         NumericMatrix X,
			         NumericMatrix x,
			         double b)   {

	cx_mat Y1 = as<cx_mat>(Y);
	mat X1 = as<mat>(X);    
	rowvec x1 = as<rowvec>(x);


    return wrap(kreg(Y1, X1, x1, b));
}






