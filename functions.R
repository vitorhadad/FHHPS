library(MASS)
library(cubature)
library(caret)
library(stats)
library(Rcpp)
library(RcppArmadillo)
library(fastGHQuad)

#Get C++ functions (takes about 2 minutes to compile)
sourceCpp("cpp_functions.cpp")


# Definitions
i1 <- complex(imaginary = 1)


# Functions
create_data <- function(n_obs, error_dist = "normal") {

    # Draw coefficients for first period
    mAB <- c(1,2)
    covAB <- matrix(c(2,1,1,2),nrow = 2)
    AB <- mvrnorm(n = n_obs, mu = mAB, Sigma = covAB)
    A <- AB[,1,drop = F]; B <- AB[,2,drop = F]

    # Draw shocks in second period
    if (error_dist == "normal") {
        mUV <- c(0,0)
        covUV <- 2*matrix(c(1,0,0,1), nrow = 2)
        W <- mvrnorm(n = n_obs, mu = mUV, Sigma = covUV)
        U <- W[,1]; V <- W[,2]
    } else if (error_dist == "laplace") {
        U <- rexp(n_obs)*(1 - 2*rbinom(n_obs, 1, 1/2))
        V <- rexp(n_obs)*(1 - 2*rbinom(n_obs, 1, 1/2))
        W <- cbind(U,V);
    } else {
        stop("unknown distribution")
    }
    

    # Draw regressors
    # X1 correlated with A,B
    X_sd <- sqrt(5)
    X1 <- matrix(0,nrow = n_obs)
    for (i in seq(n_obs))
        X1[i] <- rnorm(n = 1,mean = AB[i,1] + AB[i,2],sd = X_sd)

    # X2 is not correlated with A,B, but has same support
    X2 <- rnorm(n = n_obs, mean = mean(X1), sd = X_sd)
    X <- cbind(X1,X2)

    # Compute dependent variable
    Y1 <- A + B*X1
    Y2 <- (A + U) + (B + V)*X2
    Y <- cbind(Y1,Y2)

    # Bind everything in a list to return
    dataset <- list()
    dataset$X <- X
    dataset$Y <- Y
    dataset$A <- A
    dataset$B <- B
    dataset$W <- W

    # Assert correlations are appropriate
   # stopifnot(cor(X, A)[1] > .2)
   # stopifnot(cor(X, B)[1] > .2)
   # stopifnot(all(cor(X, W) < .05))
   # stopifnot(all(cor(W,A) < .1))
   # stopifnot(all(cor(W,B) < .1))

    return(dataset)
}



# Gamma(x) matrices
g1inv <- function(x) solve(matrix(c(1,1,x[1],x[2]), 2, 2))
g3 <- function(x)  g1inv(x) %*% matrix(c(0, 1, 0, x[2]), 2, 2)



# Computes the dependent variable
rhs <- function(s,Y,X,W, b = .1) {

	n_obs <- dim(Y)[1]
	n_features <- dim(Y)[2]

	# Make sure everything is conformable
	Y <- matrix(Y, n_obs, n_features)
	X <- matrix(X, n_obs, n_features)
	W <- matrix(W, n_obs, n_features)
	s <- matrix(s, n_features, 1)

	# Now call C++ functions
	rhs_variable(Y, W, X, s, b)

}


# Normal characteristic function
normal_cf <- function(s, mu, S) {
	s <- matrix(s, nrow = 2)
  mu <- matrix(mu, nrow = 2)
  output <- exp(-.5*t(s) %*% S %*% s)*(cos(t(s) %*% mu) + 1i*sin(t(s) %*% mu))
  return(output)
}


squared_loss <- function(v,s,Y,X,W, bw = 1) {

 	# Build the parameters
 	mu <- matrix(v[1:2], nrow = 2)
 	L <- matrix(c(v[3],0,v[4],v[5]), nrow = 2) # Cholesky of covariance
    s <- matrix(s, 2, 1)

    a <- normal_cf(s, mu, t(L) %*% L)
    b <- rhs(s,Y,X,W, bw)
 	loss <- (a - b)
 	sq_loss <- Re(loss*Conj(loss))
 	return (sq_loss)
}


ISE <- function(v, Y, X, W, lim, rule = gaussHermiteData(1000)) {

    # Prints current set of parameters
	print(v)
    
    # Vectorized square loss (necessary to use ghQuad)
    v_squared_loss <- Vectorize(function(s1,s2) squared_loss(v, c(s1,s2), Y, X, W))

    # Integrated squared error (ISE: squared loss integrated over s)
    integ <- ghQuad(function(s2)
                    ghQuad(function(s1)
                      v_squared_loss(s1,s2),
                    rule),
                  rule)
  return(integ)
}


minimize_ISE <- function(v_init, Y, X, W, lim) {

  # Create rule for Gaussian quadrature
  rule <- gaussHermiteData(1000)

  # Define objective function
  objective_fct <- function(v) ISE(v, Y, X, W, lim, rule)

  # Minimize it!
  return(optim(v_init, objective_fct))
}
