source("functions.r")

set.seed(12345)

n_obs <- 10000
b1 <- 0.1
b2 <- 0.1
t <- 0.2
d <- 1
lim <- 0.5
n_pts <- 10


# Data
n_iter <- 200
results <- matrix(c(0), n_iter, 17)

write.table(x = matrix(c("n_obs", "b1", "b2", "t", "d", "lim", "n_pts", 
                     "EU", "EV", "VarU", "CovUV", "VarV", 
                     "EA", "EB", "VarA", "CovAB", "VarB"), 1, 17), 
            file = "results_5000.csv", append = FALSE, 
            quote = FALSE, sep = ",", eol = "\n", 
            col.names = FALSE, row.names = FALSE) 
for (k in seq(n_iter)) {
    ds <- create_data(n_obs)
    Y <- matrix(ds$Y, n_obs, 2)
    X <- matrix(ds$X, n_obs, 2)
    results[k,] <- estimate(Y, X, b1, b2, t, d, lim, n_pts)
    write.table(x = matrix(results[k,], 1, 17), 
            file = "results_5000.csv", append = TRUE, 
            quote = FALSE, sep = ",", eol = "\n", 
            col.names = FALSE, row.names = FALSE) 
}


