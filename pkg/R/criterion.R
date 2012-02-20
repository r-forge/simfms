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
#   - inTrans   : the id number of the incoming transition                     #
#   - subjs     : the id numbers of the subjects concerned                     #
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
#   Last modification on: February, 20, 2012                                   #
################################################################################


criterion <- function(data,
                      atState,
                      subjs,
                      eta,
                      tmat,
                      clock,
                      marg,
                      cens,
                      copula,
                      target
  ) {
  # All possible conditioning transitions
  inTrans <- tmat[which(!is.na(tmat[, atState])), atState]
  if (!length(inTrans))
    inTrans <- NULL

  ### - Simulation of data with current parameter values
  for (it in inTrans) {
    if (is.null(it))
      it.subjs <- subjs else
        it.subjs <- which(data[subjs, paste("tr.", it, ".trans", sep="")] == 1)
    
    data <- scan.tmat(data=data, inTrans=it, subjs=subjs,
                      eta=eta,   tmat=tmat,  clock=clock,
                      marg=marg, cens=cens,  copula=copula)
  }
  
  
  outTrans <- tmat[atState, which(!is.na(tmat[atState, ]))]
  ### - Computation of criterion 
  # sum{ log( prob / \hat prob )^2}
  crit <- sum(log(
    colSums(data[subjs, paste("tr", outTrans, ".status", sep="")]) / (
      length(subjs) *
        target$prob[which(tmat %in% outTrans)] ))^2)
  # sum{ log( meds / \hat meds )^2}
  crit <- crit + sum(log(
    sapply(outTrans, function(x)
      median(data[subjs, paste("tr", x, c(".time"), sep="")][
        data[subjs, paste("tr", x, c(".status"), sep="")] == 1])) / (
          target$meds[which(tmat %in% outTrans)] ))^2)
  
  return(crit)
}
