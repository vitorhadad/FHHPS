source("functions.R")
sourceCpp("cpp_functions.cpp")


# Data
n_obs =  1000
ds <- create_data(n_obs)
Y <- matrix(ds$Y, n_obs, 2)
X <- matrix(ds$X, n_obs, 2)
W <- matrix(ds$W, n_obs, 2)


# Gamma matrices
g1 <- function(x) matrix(c(1,1,x[1],x[2]), 2, 2)
g1inv <- function(x) solve(g1(x))
g2 <- function(x) matrix(c(0, 1, 0, x[2]), 2, 2)
g3 <- function(x)  g1inv(x) %*% g2(x)

# 2D Conditional Kernel CF Estimator of G1inv(x)'Y|X=x
ft_num <- function(Y, X, x, s, H) {
    n <- dim(Y)[1]
    a <- exp(-1/2 * t(s) %*% H %*% s)
    k1 <- dnorm(X[, 1, drop = F], mean = x[1], sd = H[1,1])  
    k2 <- dnorm(X[, 2, drop = F], mean = x[2], sd = H[2,2])
    b <- t(exp(1i * t(s) %*% t(Y)))
    return (a * mean(k1 * k2 * b) / mean(k1 * k2)) 
}

# 2D Unconditional Kernel CF Estimator of G3(x)'W
ft_denom <- function(W, x, s, H) {
    a <- exp(-1/2 * t(s) %*% H %*% s)
    b <- exp(1i * t(s) %*% t(W))
    return (a %*% mean(b)) 
}

impose_lower_bound <- function(z, bound) {
    if (abs(z) < bound) {
        if (abs(z) > 0) { 
            z <- bound*z/abs(z)
        } else {
            z <- 1e-6 + 1e-6i
        }
    }
    return(z)
}



# Run through grid of s
s_pts <- 9
s_grid <- seq(-.5, .5, length.out = s_pts)
estimated_cf <- matrix(0, s_pts, s_pts)
numerator_cf <- matrix(0, s_pts, s_pts)
denominator_cf <- matrix(0, s_pts, s_pts)

mu_AB <- matrix(apply(cbind(ds$A, ds$B), 2, mean), 2, 1)
Sigma_AB <- cov(cbind(ds$A, ds$B))
target_cf <- matrix(0, s_pts, s_pts)


X <- X[abs(X[,1] - X[,2]) > .05,]
W <- X[abs(X[,1] - X[,2]) > .05,]
Y <- X[abs(X[,1] - X[,2]) > .05,]


# Numerators (Works!)
par(mfrow = c(3,2))
for (c in c(0.01, .05, .1, .5, 1, 1.5)) {
    Hy = c * diag(c(1, 1))
    bound <- .1

    for (i in seq(s_pts)) {
        print(i)
        for (j in seq(s_pts)) {
            s <- matrix(c(s_grid[i], s_grid[j]), 2, 1)
            numerator_cf[i,j] <- mean(apply(X, 1, function(x) ft_num(t(g1inv(x) %*% t(Y)), X, x, s, Hy)))
       }
    }
    persp(abs(numerator_cf), phi = 45, theta = 45)
}



