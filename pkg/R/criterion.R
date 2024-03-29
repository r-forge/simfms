################################################################################
#  Criterion function for simulation parameters of a comepeting risks block    #
################################################################################
#                                                                              #
#  Computes the criterion function sum{ log(p/hat p)^2 + log(m/hat m)^2 }      #
#  given target values of                                                      #
#   - competing events probabilities and                                       #
#   - medians of uncensored times                                              #
#                                                                              #
#  Its parameters are                                                          #
#   - data      : the dataframe with data simulated up to now                  #
#   - atState   : the state from which new transitions are considered,         #
#                 irrespective, for all its possible incoming ones             #
#   - eta       : the linear predictors matrix, with                           #
#                 as many rows as data                                         #
#                 as many columns as the number of transitions in 'tmat'       #
#   - tmat      : the trnasitions matrix                                       #
#   - clock     : either 'forward' or 'reset'                                  #
#   - marg      : the marginal baseline hazards. A list with components        #
#                 dist    : the name of the baseline hazard distribution       #
#                            (one value)                                       #
#                 eachpar : initial values of each baseline parameter          #
#                           (either one value or as many as the number         #
#                            of transitions in 'tmat')                         #
#   - cens      : the censoring time distributions. A list with components     #
#                 dist : the name of the censoring distributions (one value)   #
#                 eachpar : each censoring distribution parameter              #
#                           (either one value or as many as the number of      #
#                            possible starting states in 'tmat')               #
#                 admin: the time of administrative censoring                  #
#   - copula    : the copula model. A list with components                     #
#                 name : the name of the copula                                #
#                 par  : the copula parameter                                  #
#   - target    : target values for probabilities of competing events and for  #
#                 medians of uncensored times of each transition.              #
#                 A list with elements 'prob' and 'meds', both vector of the   #
#                 same length as the number of transitions in 'tmat'           #
#                                                                              #
#                                                                              #
#   Date: February, 20, 2012                                                   #
#   Last modification on: March, 30, 2012                                      #
################################################################################


criterion <- function(data,
                      atState,
                      eta,
                      tmat,
                      clock,
                      marg,
                      cens,
                      copula,
                      target,
                      verbose=FALSE
                      ) {
  # All possible conditioning transitions
  inTrans <- tmat[which(!is.na(tmat[, atState])), atState]
  
  ### - Simulation of data with current parameter values
  if (!length(inTrans)) {
#     subjs <- 1:nrow(data)
    data <- scan.tmat(data=data, inTrans=NULL,
                      eta=eta,   tmat=tmat,  
                      clock=clock, marg=marg, 
                      cens=cens,  copula=copula,
                      iterative=FALSE,
                      verbose=verbose) 
  } else for (it in inTrans) {
    it.subjs <- which(data[, paste("tr", it, ".status", sep="")] == 1)
    data[it.subjs, ] <- scan.tmat(data=data[it.subjs, ], inTrans=it,
                                  eta=eta[it.subjs, ],   tmat=tmat,  
                                  clock=clock, marg=marg, 
                                  cens=cens,  copula=copula,
                                  iterative=FALSE,
                                  verbose=verbose)
  }
  
  # All possible outgoing transitions
  outTrans <- tmat[atState, which(!is.na(tmat[atState, ]))]
  
  ### - Computation of criterion 
  # sum{ log( prob / \hat prob )^2}
  crit <- sum(log(
    colSums(data[, paste("tr", outTrans, ".status", sep=""), drop=FALSE]) / (
      nrow(data) *
        target$prob[which(tmat %in% outTrans)] ))^2)
  # sum{ log( meds / \hat meds )^2}
  crit <- crit + sum(log(
    sapply(outTrans, function(x)
      median(data[, paste("tr", x, c(".time"), sep="")][
        data[, paste("tr", x, c(".status"), sep="")] == 1])) / (
          target$meds[which(tmat %in% outTrans)] ))^2)
  
  return(crit)
}
