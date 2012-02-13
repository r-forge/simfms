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


trans.cancer <-
  function (names){
    tmat <- matrix(NA, 5, 5)
    tmat[1, c(2:3, 5)] <- 1:3
    tmat[2, 4:5] <- 4:5
    tmat[3, 4:5] <- 6:7
    tmat[4, 5] <- 8
    if (missing(names)) 
      names <- c("NED", "LR", "DM", "LRDM", "De")   else {
        if (length(names) != 5) 
          stop("incorrect length of \"names\" argument")
      }
    dimnames(tmat) <- list(from = names, to = names)
    return(tmat)
  }


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
