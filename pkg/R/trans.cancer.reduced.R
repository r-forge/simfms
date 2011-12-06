trans.cancer.reduced <-
function (names){
    tmat <- matrix(NA, 4, 4)
    tmat[1, 2:4] <- 1:3
    tmat[2:3, 4] <- 4:5
    if (missing(names)) 
        names <- c("NED", "LR", "DM", "De")    else {
        if (length(names) != 4)
            stop("incorrect length of \"names\" argument")
    }
    dimnames(tmat) <- list(from = names, to = names)
    return(tmat)
}
