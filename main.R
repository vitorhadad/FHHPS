library(plot3D)

# Load functions
source("functions.r")

# Data
n_obs =  1500
ds <- create_data(n_obs, "normal")
Y <- matrix(ds$Y, n_obs, 2)
X <- matrix(ds$X, n_obs, 2)
W <- matrix(ds$W, n_obs, 2)


# Drop obs that have x[1] approx x[2]
valid_obs <- (abs(X[,1] - X[,2]) > 1e-3)
Y <- Y[valid_obs,]
X <- X[valid_obs,]
W <- W[valid_obs,]




# Plot the variables
matching_cf <- matrix(0, n_pts, n_pts)
true_cf <- matrix(0, n_pts, n_pts)
numerator <- matrix(0, n_pts, n_pts)
denominator <- matrix(0, n_pts, n_pts)
ratio <- matrix(0, n_pts, n_pts)
b <- .3*diff(range(X[,2]))
mu_AB <- matrix(c(mean(ds$A), mean(ds$B)), 2, 1)
Sigma_AB <- cov(cbind(ds$A, ds$B))

s_grid <- seq(-1/b, 1/b, by = .05)
n_pts <- length(s_grid)


for (i in seq(n_pts)) {
	for (j in seq(n_pts)) {
		s <- matrix(c(s_grid[i], s_grid[j]), 2, 1)
		#true_cf[i,j] <- normal_cf(s, mu_AB, Sigma_AB)
		#matching_cf[i,j] <- rhs(s, Y, X, W, b = b)
        
        num <- cf_numerator_R(Y, X, s, b)
        denom <- cf_denominator_R(Y, X, s)
        numerator[i,j] <- mean(num)
        denominator[i,j] <- mean(denom)
     }
}


par(mfrow = c(2,2))
persp(abs(true_cf), main = "True CF", 
      theta = 30, phi = 30, col = "lightblue")
persp(abs(matching_cf), main = "Estimated CF", 
      theta = 30, phi = 30, col = "lightblue")

persp(abs(numerator), main = "E[Numerator]", 
      theta = 30, phi = 30, col = "lightblue")
persp(abs(denominator), main = "E[Denominator]", 
      theta = 30, phi = 30, col = "lightblue")



s = matrix(c(1,1), 2, 1)
denom <- cf_denominator_R(Y, X, s)
num <- cf_numerator_R(Y, X, s, b)




 #      // ratio[i,j] <- mean(num/denom)

