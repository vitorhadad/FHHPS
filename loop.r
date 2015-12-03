rm(list = ls())

# Test client

# Variations on parameters
params <- NULL
params$n_obs <- c(500, 2000, 5000, 10000)
params$b1 <- c(.1, .2, .3)
params$b2 <- c(.1, .2, .3)
params$t <- c(.2, .4, .8)
params$d <- c(.5, 1, 1.5, 2)
params$lim <- c(.1, .5)
params$n_pts <- c(10)
variations <- expand.grid(params)

for (k in seq(dim(variations)[1])) {
    
    # Open main file
    lines <- readLines("main.r")  
    # Open file to be written
    main_name <- paste("new_main_",k,".r",sep="")
    new_main <- file(main_name, "w")
    # Write parameters
    for (param in names(variations)) {
        writeLines(paste(param, variations[k, param], sep = " <- "),
               con = new_main)
    }
    # Write rest of file
    writeLines(lines, con = new_main)
    close(new_main)

    # Create new job file
    job_name <- paste("j_main_",k,".pbs",sep="")
    new_job <- file(job_name, "w")
    lines <- readLines("j_main.pbs")
    writeLines(lines, con = new_job)
    writeLines(paste("R < ", main_name, " --no-save", sep = ""),
               con = new_job)
    close(new_job)



    # Submit to queue    
    system(paste("qsub ", job_name, sep = ""))    

    Sys.sleep(5)

}
