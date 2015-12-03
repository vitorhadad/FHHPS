n_obs <- 500
b1 <- 0.1
b2 <- 0.1
t <- 0.2
d <- 0.5
lim <- 0.1
n_pts <- 10

source("functions.r")

# Data
n_iter <- 8
results <- matrix(c(0), n_iter, 17)

for (k in seq(n_iter)) {
    ds <- create_data(n_obs)
    Y <- matrix(ds$Y, n_obs, 2)
    X <- matrix(ds$X, n_obs, 2)
    results[k,] <- estimate(Y, X, b1, b2, t, d, lim, n_pts) 
}

write.table(x = matrix(apply(results, 2, mean), 1, 17), 
            file = "results1.csv", append = TRUE, 
            quote = FALSE, sep = ",", eol = "\n", 
            col.names = FALSE, row.names = FALSE)
