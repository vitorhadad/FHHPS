library(MASS)
library(Rcpp)
library(RcppArmadillo)


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
    mUV <- c(.5,.5)
    covUV <- matrix(c(1,0,0,1), nrow = 2)
    W <- mvrnorm(n = n_obs, mu = mUV, Sigma = covUV)
    U <- W[,1]; V <- W[,2]
   
    # Draw regressors
    # X1 correlated with A,B
    X_sd <- sqrt(5)
    X1 <-  .2*A + .5*B + 
        .5*A^2 - .3*B^2 + 
        rnorm(n = 1,mean = 0, sd = X_sd)
    # [[De-mean!]]
    X1 <- X1 - mean(X1)

    # X2 is not correlated with A,B, but has same support
    X2 <- rnorm(n = n_obs, mean = 0, sd = X_sd)
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

    return(dataset)
}

# Gamma matrices
g1 <- function(x) matrix(c(1,1,x[1],x[2]), 2, 2)
g1inv <- function(x)  solve(g1(x))
g2 <- function(x) matrix(c(0, 1, 0, x[2]), 2, 2)
g3 <- function(x)  g1inv(x) %*% g2(x)


# Normal characteristic function
normal_cf <- function(s, mu, Sig) {
    s <- matrix(c(s[1], s[2]), nrow = 2)
    return(exp(1i * t(s) %*% mu - .5 * t(s) %*% Sig %*% s))
}

# Shock pameters
get_shock_params <- function(Y, X, b1) {

    n_obs1 <- dim(Y)[1]
    Y <- matrix(Y, n_obs1, 2) # Ensure it's in matrix form
    X <- matrix(X, n_obs1, 2)

    # First moments
    DY <- Y[,2] - Y[,1]
    shock_means <- lp(DY, 
                      cbind(1, X[,2]), 
                      X[,2] - X[,1], 
                      b = b1)
    EU <- shock_means[1]
    EV <- shock_means[2]
    muW <- matrix(c(EU, EV), 2, 1)

    # Second moments
    shock_cov <- lp((DY - EU - EV*X[,2])^2, 
                    cbind(1, 2*X[,2], X[,2]^2), 
                    X[,2] - X[,1],
                    b = b1)
    VarU <- shock_cov[1]
    CovUV <- shock_cov[2]
    VarV <- shock_cov[3]

    SigmaW <- matrix(c(VarU, CovUV, CovUV, VarV), 2, 2)

    return(list(muW = muW, SigmaW = SigmaW))
}

# Characteristic function of (A,B)
chf_AB <- function(s, X, Y, muW, SigmaW, d, b2, t) {

    # Keep if |X(2)-X(1)| > d
    valid <- abs(X[,2] - X[,1]) > d
    X <- X[valid,]
    Y <- Y[valid,]
    n_obs <- sum(valid)

    # Conformable s
    s <- matrix(c(s[1], s[2]), 2, 1) 

    # Initialize ratio
    ratio <- rep(0 + 0i, n_obs)

    for (i in seq(n_obs)) {
        # Get evaluation points
        x <- matrix(X[i,], 1, 2)
        # Numerator
        num_var <- matrix(exp(1i * t(s) %*% g1inv(x) %*% t(Y[-i,])), n_obs-1, 1)
        num <- kreg_R(num_var, X[-i,], x, b2)
       
        # Denominator
        denom <- normal_cf(s, g3(x) %*% muW, g3(x) %*% SigmaW %*% t(g3(x)))  
        denom <- ifelse(abs(denom) > t, denom, denom/abs(denom)*t) 
        
        # Compute ratio
        ratio[i] <- num/denom
    }
    return(mean(ratio, na.rm = TRUE))
}


# Second step objective function
obj <- function(z, chf2, s_grid, X, Y, muW, SigmaW, d, b2, t) {
  
    muAB <- matrix(z[1:2], 2, 1)
    LAB <- matrix(c(z[3], 0, z[4], z[5]), 2, 2) 
    SigmaAB <- t(LAB) %*% LAB

    loss <- rep(0, dim(s_grid)[1])
    for (i in seq(dim(s_grid)[1])) {
        s <- matrix(c(s_grid[i,1], s_grid[i,2]), 2, 1)
        chf1 <- normal_cf(s, muAB, SigmaAB)
        loss[i] <- chf1 - chf2[i]
    }
    return(sum(loss*Conj(loss)))
}


estimate <- function(Y, X, b1, b2, t, d, lim, n_pts) {

    # Shock moments
    shock_params <- get_shock_params(Y, X, b1) 
    muW <- shock_params$muW 
    SigmaW <- shock_params$SigmaW 

    # Evaluation grid
    s_grid <- expand.grid(runif(n_pts, min = -lim, max = +lim),
                          runif(n_pts, min = -lim, max = +lim))

    chf2 <- rep(0, n_pts^2)
    for (i in seq(dim(s_grid)[1])) {
        s <- matrix(c(s_grid[i,1], s_grid[i,2]), 2, 1)
        chf2[i] <- chf_AB(s, X, Y, muW, SigmaW, d, b2, t) 
    }

    z_guess <- c(1,1,1,0,1)
    result <- optim(z_guess,
                    function(z) obj(z, chf2, s_grid, X, Y, muW, SigmaW, d, b2, t))

    z <- result$par

    muAB <- matrix(z[1:2], 2, 1)
    LAB <- matrix(c(z[3], 0, z[4], z[5]), 2, 2) 
    SigmaAB <- t(LAB) %*% LAB

    estimates <- c(n_obs, b1, b2, t, d, lim, n_pts,
                   c(c(muW), 
                   c(SigmaW[1,1], SigmaW[1,2], SigmaW[2,2]),
                   c(muAB), 
                   c(SigmaAB[1,1], SigmaAB[1,2], SigmaAB[2,2])))

    return(estimates)
}

# Local polynomial regression 
lp <- function(YY, XX, xx, b) {
    K <- diag(dnorm((xx)/b))
    B <- solve(t(XX) %*% K %*% XX) %*% (t(XX) %*% K %*% YY)
    return(B)
}

