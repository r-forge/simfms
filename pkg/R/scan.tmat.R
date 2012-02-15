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
#   - nsim      : the number of subjects to simulate                           #
#   - tmat      : the trnasitions matrix                                       #
#   - clock     : either 'forward' or 'reset'                                  #
#   - frailty   : the frailty term specifications. A list with components      #
#                 dist: the name of the frailty distribution                   #
#                 par : the frailty parameter(s) value                         #
#   - nclus     : the number of clusters to simulate                           #
#   - csize     : the size(s) of cluster                                       #
#   - covs      : the covariates to simulate. A list with components           #
#                 nameOfCovariate = simulationfunction                         #
#   - beta      : the regression coefficients for covariates. A list of the    #
#                 same length as 'covs', with elements                         #
#                 nameOfCovariate = a vector of the same length as the number  #
#                 of transitions in 'tmat'                                     #
#   - marg      : the marginal baseline hazards. A list with components        #
#                 dist    : the name of the baseline hazard distribution       #
#                           (either one value or as many as the number         #
#                            of transitions in 'tmat')                         #
#                 eachpar : each baseline parameter                            #
#                           (either one value or as many as the number         #
#                            of transitions in 'tmat')                         #
#   - cens      : the censoring time distributions. A list with components     #
#                 dist : the name of the censoring distributions               #
#                           (either one value or as many as the number of      #
#                            possible starting states in 'tmat')               #
#                 eachpar : each censoring distribution parameter              #
#                           (either one value or as many as the number of      #
#                            possible starting states in 'tmat')               #
#   - copula    : the copula model. A list with components                     #
#                 name : the name of the copula                                #
#                 par  : the copula parameter                                  #
#                                                                              #
#                                                                              #
#   Date: February, 13, 2012                                                   #
#   Last modification on: February, 15, 2012                                   #
################################################################################

# scan.tmat <- function(data,
#                       inTrans,
#                       subjs,
#                       # Other parameters
#                       nsim,
#                       tmat,
#                       clock,
#                       frailty,
#                       nclus,
#                       csize,
#                       covs,
#                       beta,
#                       marg,
#                       cens,
#                       copula
#                       ){
  ### - PREPARATION - ##########################################################
  # Present state and Conditioning transition infos
  if (is.null(inTrans)){ # from the starting state
    atState <- colnames(tmat)[which(colSums(tmat, na.rm=TRUE) == 0)]
    condTime <- condMarg <- NULL
  } else { # from all the other states
    atState <- colnames(tmat)[which(tmat == inTrans, arr.ind=TRUE)[2]]
    condTime <- data[subjs, paste(atState, "time", sep=".")]
    condMarg <- extractMargs(as.data.frame(marg)[inTrans,])
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
    if (ot.N == 1)
      prevTimes <- prevMargs <- NULL else {
        prevOTs <- outTrans[1:(ot.N - 1)]
        prevTimes <- data[subjs, paste(names(prevOTs), "time", sep="."),
                          drop=FALSE]
        prevMargs <- apply(as.data.frame(marg)[prevOTs, ], 1, extractMargs)
      }
    
    data[subjs, paste(names(outTrans[outTrans == ot]), "time", sep=".")] <-
      sapply(subjs, function(x) eval(parse(text=copula$name))(
        par=copula$par,
        condTime=condTime[x], condMarg=condMarg,
        prevTimes=prevTimes[x,, drop=FALSE], prevMargs=prevMargs,
        marg=extractMargs(as.data.frame(marg)[ot,]),
        eta=eta[x, c(inTrans, outTrans[1:ot.N]), drop=FALSE], 
        clock=clock))
  }
  ######################################## - END of COMPETING EVENTS TIMES - ###
  
  
  ### - CENSORING - ############################################################
  C.time <- extractMargs(as.data.frame(cens[names(cens)!="admin"])[atState,])(
    runif(length(subjs)), inv=TRUE)
  ##################################################### - END of CENSORING - ###
    
  
  ### - UPDATE DATASET - #######################################################
  data[subjs, sapply(c("time", "status"), function(x)
    paste(names(outTrans), x, sep="."))] <-
      t(apply(cbind(data[subjs, paste(names(outTrans), "time", sep=".")],
                    C.time=C.time), 1, function(x)
                      c(rep(min(x), length(x)-1), 
                        1:(length(x)-1) == which.min(x))))
  ################################################ - END of UPDATE DATASET - ###
  
  
  ### - NEXT EVENTS TIMES - ####################################################
  for (ot in outTrans) { # ot, the number of the transition in tmat
    # find out concerned subjects
    # call scan.tmat on them
  }
  ############################################# - END of NEXT EVENTS TIMES - ###
  
# }
  
  