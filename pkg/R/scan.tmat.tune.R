################################################################################
#  Tuning of simulation parameters for a comepeting risks block                #
################################################################################
#                                                                              #
#  Tunes iteratively the simulation parameters of a competing risks blcock,    #
#     according to given target values for probabilities of competing events   #
#     and of medians of uncensored times                                       #
#                                                                              #
#  Its parameters are                                                          #
#   - pars      : the list with parameters values chosen till now              #
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

scan.tmat.tune <- function(pars,
                           data,
                           inTrans,
                           subjs,
                           eta,
                           # Other parameters from tune.simfms()
                           tmat,
                           clock,
                           marg,
                           cens,
                           copula
                           ){
  ### - PREPARATION - ##########################################################
  # Present state and Conditioning transition infos
  if (is.null(inTrans)){ # from the starting state
    atState <- colnames(tmat)[which(colSums(tmat, na.rm=TRUE) == 0)]
    condTime <- condMarg <- NULL
  } else { # from all the other states
    atState <- colnames(tmat)[which(tmat == inTrans, arr.ind=TRUE)[2]]
    condTime <- data[, paste("tr", inTrans, ".time", sep="")]
    condMarg <- extractMargs(as.data.frame(marg)[inTrans,])
  }
  outTrans <- tmat[atState, which(!is.na(tmat[atState, ]))]
  # if ending state, then return results
  if (length(outTrans) == 0)
    return(pars=pars)
  ################################################### - END of PREPARATION - ###

  
  ### - COMPETING RISKS PARAMETERS TUNING - #####################################
  
  ######################################################### - END of TUNING - ###
  
  
  ### - UPDATE PARAMETERS - #####################################################
#   data[subjs, sapply(c(".time", ".status"), function(x)
#     paste("tr", outTrans, x, sep=""))] <-
#       t(apply(cbind(data[subjs, paste("tr", outTrans, ".time", sep="")],
#                     C.time=C.time), 1, function(x)
#                       c(rep(min(x), length(x)-1), 
#                         1:(length(x)-1) == which.min(x))))
  ############################################## - END of UPDATE PARAMETERS - ###
  
  
  ### - NEXT CRs BLOCKS - #######################################################
#   for (ot in outTrans) { # ot, the number of the transition in tmat
#     # find out concerned subjects
#     subjs <- data[data[[paste("tr", ot, ".status", sep="")]] > 0, "ID"]
#     # call scan.tmat on them
#     if (length(subjs))
#       data <- scan.tmat(pars=pars, 
#                         data=data, inTrans=ot, subjs=subjs,
#                         eta=eta,   tmat=tmat,  clock=clock,
#                         marg=marg, cens=cens,  copula=copula)
#   }
  ################################################ - END of NEXT CRs BLOCKS - ###
  
  return(pars=pars)
}
