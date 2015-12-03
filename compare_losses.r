
table <- read.csv("results_10000.csv")
colnames(table) <- c("n_obs", "b1", "b2", "t", "d", "lim", "n_pts", 
                     "EU", "EV", "VarU", "CovUV", "VarV", 
                     "EA", "EB", "VarA", "CovAB", "VarB")

truth <- c(.5, .5, 1, 0, 1, 1, 2, 2, 1, 2)
names(truth) <- c("EU", "EV", "VarU", "CovUV", "VarV", 
                     "EA", "EB", "VarA", "CovAB", "VarB")


vars <- c("EU", "EV", "VarU", "CovUV", "VarV", 
                     "EA", "EB", "VarA", "CovAB", "VarB")
cond <- which(table[,"n_obs"] == 10000)

tab <- table[cond,]

absloss <- apply(tab, 1, function(x) mean(abs(x[vars] - truth[vars])))
rmsloss <- apply(tab, 1, function(x) sqrt(mean((x[vars] - truth[vars])^2)))

id_rms <- which.min(rmsloss)
id_abs <- which.min(absloss)

tab[id_rms,]
tab[id_abs,]

