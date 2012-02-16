################################################################################
#  Transition matrices for cancer mulit-state structure                        #
################################################################################
#                                                                              #
#   Date: February, 14, 2012                                                   #
#   Last modification on: February, 16, 2012                                   #
################################################################################


## - REDUCED MODEL, without LRDM - #############################################
trans.cancer.reduced <- function (names=c("NED", "LR", "DM", "De")) {
  tmat <- matrix(NA, 4, 4)
  tmat[1, 2:4] <- 1:3
  tmat[2:3, 4] <- 4:5
  
  if (length(names) != 4)
    stop("incorrect length of \"names\" argument")

  dimnames(tmat) <- list(from = names, to = names)
  return(tmat)
}
#################################################### - END of REDUCED MODEL - ##


## - STANDARD MODEL - ##########################################################
trans.cancer <- function (names=c("NED", "LR", "DM", "LRDM", "De")) {
  tmat <- matrix(NA, 5, 5)
  tmat[1, c(2:3, 5)] <- 1:3
  tmat[2, 4:5] <- 4:5
  tmat[3, 4:5] <- 6:7
  tmat[4, 5] <- 8
  
  if (length(names) != 5)
    stop("incorrect length of \"names\" argument")

  dimnames(tmat) <- list(from = names, to = names)
  return(tmat)
}
################################################### - END of STANDARD MODEL - ##


## - EXTENDED MODEL, with direct transition into LRDM - ########################
trans.cancer.extended <- function (names=c("NED", "LR", "DM", "LRDM", "De")) {
  tmat <- matrix(NA, 5, 5)
  tmat[1, c(2:5)] <- 1:4
  tmat[2, 4:5] <- 5:6
  tmat[3, 4:5] <- 7:8
  tmat[4, 5] <- 9
  
  if (length(names) != 5)
    stop("incorrect length of \"names\" argument")

  dimnames(tmat) <- list(from = names, to = names)
  return(tmat)
}
################################################### - END of EXTENDED MODEL - ##
