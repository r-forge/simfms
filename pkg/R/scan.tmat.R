################################################################################
#  Simualation of times for a comepeting risks block                           #
################################################################################
#                                                                              #
#  Computes iteratively the times of a competing risks blcock, according       #
#     to a given copula model, with a possible conditioning event into the     #
#     present state and conditional, too, to the previously simulated          #
#     competing events in the block                                            #
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
#   - iterative : shall the simulation continue on children transitions?       #
#                                                                              #
#                                                                              #
#   Date: February, 13, 2012                                                   #
#   Last modification on: March, 29, 2012                                      #
################################################################################

scan.tmat <- function(data,
                      inTrans,
#                       subjs,
                      eta,
                      # Other parameters from simfms()
                      tmat,
                      clock,
                      marg,
                      cens,
                      copula,
                      iterative = TRUE
                      ) {
  ### - PREPARATION - ##########################################################
  # Present state and Conditioning transition infos
  if (is.null(inTrans)){ # from the starting state
    atState <- colnames(tmat)[which(colSums(tmat, na.rm=TRUE) == 0)]
    condTime <- condMarg <- NULL
  } else { # from all the other states
    atState <- colnames(tmat)[which(tmat == inTrans, arr.ind=TRUE)[2]]
    condTime <- data[, paste("tr", inTrans, ".time", sep="")]
    condMarg <- marg[[paste(inTrans)]]
  }
  outTrans <- tmat[atState, which(!is.na(tmat[atState, ]))]
  
  # if ending state, then return results
  if (length(outTrans) == 0)
    return(data)
  ################################################### - END of PREPARATION - ###

  
  ### - COMPETING EVENTS TIMES - ###############################################
  for (ot in outTrans) { # ot, the number of the transition in tmat!!!!!!!!!!!!!
    ot.N <- which(outTrans == ot) # ot.N its rank in the CRs block!!!!!!!!!!!!!!
    # Previous transition(s) infos
    if (ot.N == 1) {
      prevOTs <- prevTimes <- prevMargs <- NULL
    } else {
      prevOTs <- outTrans[1:(ot.N - 1)]
      prevTimes <- data[, paste("tr", prevOTs, ".time", sep=""),
                        drop=FALSE]
      prevMargs <- marg[paste(prevOTs)]
    }
    
    for (subj in 1:nrow(data)) {
      data[subj, paste("tr", outTrans[ot.N], ".time", sep="")] <-
        eval(parse(text=copula$name))(
          par=copula$par,
          condTime=condTime[subj],
          condMarg=condMarg,
          trans=ot, marg=marg[[paste(ot)]],
          prevTimes=prevTimes[subj, ], prevMargs=prevMargs,
          eta=eta[subj, c(inTrans, prevOTs, ot)], tmat=tmat,
          clock=clock)
    }
  }
  ######################################## - END of COMPETING EVENTS TIMES - ###
  
  
  ### - CENSORING - ############################################################
  C.time <- sapply(cens$f[[paste(atState)]](runif(nrow(data)), inv=TRUE), 
                   function(x) min(x, cens$admin))
    
  ##################################################### - END of CENSORING - ###
    
  
  ### - UPDATE DATASET - #######################################################
  data[, sapply(c(".time", ".status"), function(x)
    paste("tr", outTrans, x, sep=""))] <-
      t(apply(cbind(data[, paste("tr", outTrans, ".time", sep="")],
                    C.time=C.time), 1, function(x)
                      c(rep(min(x), length(x)-1), 
                        1:(length(x)-1) == which.min(x))))
  ################################################ - END of UPDATE DATASET - ###
  
  
  ### - NEXT EVENTS TIMES - ####################################################
  if (iterative) {
    for (ot in outTrans) { # ot, the number of the transition in tmat
      # find out concerned subjects
      subjs <- which(data[[paste("tr", ot, ".status", sep="")]] > 0)
      # call scan.tmat on them
      if (length(subjs)) {
        data[subjs, ] <- scan.tmat(data=data[subjs, ], inTrans=ot, #subjs=subjs,
                                   eta=eta[subjs, ],   tmat=tmat,  
                                   clock=clock,        marg=marg,
                                   cens=cens,          copula=copula)
      }
    }
  }
  ############################################# - END of NEXT EVENTS TIMES - ###
  
  return(data)
}
