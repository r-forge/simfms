scan.tmat <- function(data,
                      inTrans,
                      subjs,
                      # Other parameters
                      nsim,
                      tmat,
                      clock,
                      frailty,
                      nclus,
                      csize,
                      covs,
                      beta,
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
    condTime <- data[subjs, paste(atState, "time", sep=".")]
    condMarg <- extractMargs(as.data.frame(marg)[inTrans,])
#       eval(parse(text=marg[inTrans]$dist))(
#       pars=eval(parse(text=paste("list(",
#                                  paste(names(marg)[names(marg) != "dist"],
#                                        as.data.frame(marg)[inTrans, 
#                                                            names(marg) !="dist"],
#                                        sep="=", collapse=", "),
#                                  ")"))))        
  }
  outTrans <- tmat[atState, which(!is.na(tmat[atState, ]))]
  ################################################### - END of PREPARATION - ###
      
  
  ### - RECURSION on the COMPETING RISKS BLOCK - ###############################
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
  ########################################## - END of RECURSION on the CRs - ###
  
  # if no  children, then return results
  if (length(outTrans == 0))
    return(data)
  # else
  else
    {}
  ## pass results to its children
  ## merge their results
  ## return results
}