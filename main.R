library(plot3D)

# Load functions
source("functions.r")

# Data
n_obs =  500
ds <- create_data(n_obs, "normal")
Y <- matrix(ds$Y, n_obs, 2)
X <- matrix(ds$X, n_obs, 2)
W <- matrix(ds$W, n_obs, 2)
lim = 2


# Drop obs that have x[1] approx x[2]
valid_obs <- (abs(X[,1] - X[,2]) > 1e-3)
Y <- Y[valid_obs,]
X <- X[valid_obs,]
W <- W[valid_obs,]


# Plot the variables
s_grid <- seq(-1, 1, by = .2)
n_pts <- length(s_grid)
matching_cf <- matrix(0, n_pts, n_pts)
true_cf <- matrix(0, n_pts, n_pts)
numerator <- matrix(0, n_pts, n_pts)
denominator <- matrix(0, n_pts, n_pts)
ratio <- matrix(0, n_pts, n_pts)
b <- 20
mu_AB = matrix(c(mean(ds$A), mean(ds$B)), 2, 1)
cov(cbind(ds$A, ds$B))

for (i in seq(n_pts)) {
	for (j in seq(n_pts)) {
		s <- matrix(c(s_grid[i], s_grid[j]), 2, 1)
		true_cf[i,j] <- normal_cf(s, mu_AB, Sigma_AB)
		matching_cf[i,j] <- rhs(s, Y, X, W, b = b)
     }
}


par(mfrow = c(1,2))
persp(abs(true_cf), main = "True CF", theta = 30, phi = 30, col = "lightblue")
persp(abs(matching_cf), main = "Estimated CF", theta = 30, phi = 30,col = "lightblue")







 # num <- cf_numerator_R(Y, X, s, b)
 # denom <- cf_denominator_R(Y, X, s)
 #      // numerator[i,j] <- mean(num)
 #      // denominator[i,j] <- mean(denom)
 #      // ratio[i,j] <- mean(num/denom)

