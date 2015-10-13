source("functions.r")

# Data
n_obs =  1000
ds <- create_data(n_obs, "normal")
Y <- matrix(ds$Y, n_obs, 2)
X <- matrix(ds$X, n_obs, 2)
W <- matrix(ds$W, n_obs, 2)

# Local means
E_Y1_x = nw_loo(Y[,1, drop = FALSE], X, b = 1)
E_Y2_x = nw_loo(Y[,2, drop = FALSE], X, b = 1)
E_Y_x = cbind(E_Y1_x, E_Y2_x)

Y1sq_demeaned = matrix((Y[,1]-mean(Y[,1]))^2, n_obs, 1)
Y2sq_demeaned = matrix((Y[,2]-mean(Y[,2]))^2, n_obs, 1)
Y1Y2_demeaned = matrix((Y[,1]-mean(Y[,1]))*(Y[,2]-mean(Y[,2])), n_obs, 1)

E_Y1sq_x = nw_loo(Y1sq_demeaned, X, b = 10)
E_Y2sq_x = nw_loo(Y2sq_demeaned, X, b = 10)
E_Y1Y2_x = nw_loo(Y1Y2_demeaned, X, b = 10)


Cov_W = cov(W)

E_AB_x = matrix(0, n_obs, 2)
Var_A_x = matrix(0, n_obs)
Var_B_x = matrix(0, n_obs)
Cov_AB_x = matrix(0, n_obs)

for (i in seq(n_obs)) {
    # Means
    E_AB_x[i,] = g1inv(X[i,]) %*% matrix(E_Y_x[i,],2,1)
    
    # Covariances
    Cov_Y_x = matrix(c(E_Y1sq_x[i],
                       E_Y1Y2_x[i],
                       E_Y1Y2_x[i],
                       E_Y2sq_x[i]),
                       2,2)

    tmp = g1inv(X[i,]) %*% Cov_Y_x %*% t(g1inv(X[i,])) +
          g3(X[i,]) %*% Cov_W %*% t(g3(X[i,]))

    Var_A_x[i] = tmp[1,1]
    Var_B_x[i] = tmp[2,2]
    Cov_AB_x[i] = tmp[1,2]
}


E_AB = apply(E_AB_x, 2, mean)


E_Cov_AB_x = matrix(c(mean(Var_A_x), mean(Cov_AB_x), 
                      mean(Cov_AB_x), mean(Var_B_x)), 2, 2)

Cov_E_AB_x = cov(E_AB_x)
Var_AB = E_Cov_AB_x  + Cov_E_AB_x 








