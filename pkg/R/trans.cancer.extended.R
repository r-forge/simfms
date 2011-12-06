trans.cancer.extended <-
function (names){
    tmat <- matrix(NA, 5, 5)
    tmat[1, c(2:5)] <- 1:4
    tmat[2, 4:5] <- 5:6
    tmat[3, 4:5] <- 7:8
    tmat[4, 5] <- 9
    if (missing(names)) 
        names <- c("NED", "LR", "DM", "LRDM", "De")   else {
        if (length(names) != 5) 
            stop("incorrect length of \"names\" argument")
    }
    dimnames(tmat) <- list(from = names, to = names)
    return(tmat)
}
