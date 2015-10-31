rm(list = ls())

source("functions.R")
sourceCpp("cpp_functions.cpp")

# Data
n_obs =  10000
ds <- create_data(n_obs)
Y <- matrix(ds$Y, n_obs, 2)
X <- matrix(ds$X, n_obs, 2)
W <- matrix(ds$W, n_obs, 2)

# Rule out small determinants
lower_bound = 1 
valid <- sapply(seq(n_obs), function(i) (det(g1(X[i,])) > lower_bound))
Y <- Y[valid,]
X <- X[valid,]
W <- W[valid,]
n_obs <- sum(valid)

# Numerator
numerator_variable <- function(s, x) t(exp(1i * t(s) %*% g1inv(x) %*% t(Y)))
numerator_regression <- function(s) apply(X, 1,
                             function(x) kreg_R(Y = numerator_variable(s, x),
                                                X = X,
                                                matrix(x, 1, 2),
                                                b = diff(range(X))*.01))
# Denominator
denominator_variable <- function(s, x) t(exp(1i * t(s) %*% g3(x) %*% t(W)))
trim <- function(z) return(if (abs(z) < .1) { .1*z/abs(z) } else {z}) 
denominator_regression <- function(s) apply(X, 1,
                             function(x) trim(kreg_R(Y = denominator_variable(s, x),
                                                X = X,
                                                matrix(x, 1, 2),
                                                b = diff(range(X))*.01)))
# Mean of ratio
ratio <- function(s) mean(numerator_regression(s) / denominator_regression(s))


# Normal cf (estimation target)
mean_AB <- matrix(c(mean(ds$A), mean(ds$B)), 2, 1)
cov_AB <- cov(cbind(ds$A, ds$B))
normal_cf <- function(s) mean(exp(1i * t(s) %*% mean_AB 
                                  - .5 * t(s) %*% cov_AB %*% s))



# Evaluate over grid
n_pts <- 13
s_grid <- seq(-.5, .5, length = n_pts)
estimated <- matrix(0, n_pts, n_pts)
target <- matrix(0, n_pts, n_pts)
denominator <- matrix(0, n_pts, n_pts)
numerator <- matrix(0, n_pts, n_pts)

for (i in seq(n_pts)) {
    print(i)
    for (j in seq(n_pts)) {
        s <- matrix(c(s_grid[i], s_grid[j]), 2, 1)
        estimated[i,j] <- ratio(s) 
        target[i,j] <- normal_cf(s)
        denominator[i,j] <- mean(denominator_regression(s))
        numerator[i,j] <- mean(numerator_regression(s))
    }
}


dev.new()
par(mfrow = c(2,2))
persp(s_grid, s_grid, abs(estimated), theta = 45, phi = 45, main = "estimated")
persp(s_grid, s_grid, abs(target), theta = 45, phi = 45, main = "target")
persp(s_grid, s_grid, abs(numerator), theta = 45, phi = 45, main = "mean(numerator)")
persp(s_grid, s_grid, abs(denominator), theta = 45, phi = 45, main = "mean(denominator)")

save.image("minimal_implementation.RData")

dev.new()
par(mfrow = c(2,2))
contour(s_grid, s_grid, abs(estimated), main = "estimated")
contour(s_grid, s_grid, abs(target), main = "target")
contour(s_grid, s_grid, abs(numerator), main = "mean(numerator)")
contour(s_grid, s_grid, abs(denominator), main = "mean(denominator)")

save.image("minimal_implementation.RData")




