library(plot3D)

# Load functions
source("functions.r")


# Data
n_obs =  2000
ds <- create_data(n_obs, "normal")
Y <- matrix(ds$Y, n_obs, 2)
X <- matrix(ds$X, n_obs, 2)
W <- matrix(ds$W, n_obs, 2)


# Drop obs that have x[1] approx x[2]
valid_obs <- (abs(X[,1] - X[,2]) > 1)
Y <- Y[valid_obs,]
X <- X[valid_obs,]
W <- W[valid_obs,]
n_obs <- sum(valid_obs)

# Parameters
lim <- .5
b <- .01
mu_AB <- matrix(c(mean(ds$A), mean(ds$B)), 2, 1) # Oracle mean
Sigma_AB <- cov(cbind(ds$A, ds$B)) # Oracle covariance

# Get grid
n_pts <- 7
s_grid <- seq(-lim, +lim, length = n_pts)

# Look at chf contour difference
oracle_normal <- normal_cf_matrix(s_grid, mu_AB, Sigma_AB)

oracle_ratio <- target_cf(Y, W, X, s_grid, b)
l2loss <- Re((oracle_normal - oracle_ratio) * Conj(oracle_normal - oracle_ratio))
l1loss <- abs(oracle_normal) - abs(oracle_ratio) 

# Plot
dev.new()
par(mfrow = c(2,2))
contour(abs(oracle_normal), main = "oracle (normal)")
contour(abs(oracle_ratio), main = "oracle (ratio)")
contour(l1loss, main = "difference btw absolute values")
contour(l2loss, main = "difference btw squared values")

# Plot
dev.new()
par(mfrow = c(1,4))
persp(abs(oracle_normal), main = "oracle (normal)")
persp(abs(oracle_ratio), main = "oracle (ratio)")

persp(l2loss, main = "squared difference")



v_true <- c(mu_AB[1], mu_AB[2], L[1,1], L[1,2], L[2,2])
v_perturb <- v_true + 0.05*rnorm(5)

loss <- function(mu, Sig) abs(normal_cf_matrix(s_grid, mu, Sig)) - abs(target_cf(Y, W, X, s_grid, b))
sqloss <- function(mu, Sig) { 
    L <- loss(mu, Sig)
    return(Re(sum(L * Conj(L))))
}
