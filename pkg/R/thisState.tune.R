################################################################################
#  Tuning of simulation parameters for a present state                         #
################################################################################
#                                                                              #
#  Tunes the simulation parameters for a ''present state'', that is for its    #
#     CRs block and for all possible incoming transitions,                     #
#     according to given target values for probabilities of competing events   #
#     and of medians of uncensored times                                       #
#                                                                              #
#  Its parameters are                                                          #
#   - target    : target values for probabilities of competing events and for  #
#                 medians of uncensored times of each transition.              #
#                 A list with elements 'prob' and 'meds', both vector of the   #
#                 same length as the number of transitions in 'tmat'           #
#   - data      : the dataframe with data simulated up to now                  #
#   - atState   : the state from which new transitions are considered,         #
#                 irrespective, for all its possible incoming ones             #
#   - subjs     : the id numbers of the subjects concerned                     #
#   - eta       : the linear predictors matrix, with                           #
#                 as many rows as data                                         #
#                 as many columns as the number of transitions in 'tmat'       #
#   - tmat      : the trnasitions matrix                                       #
#   - clock     : either 'forward' or 'reset'                                  #
#   - marg      : the marginal baseline hazards. A list with components        #
#                 dist    : the name of the baseline hazard distribution       #
#                            (one value)                                       #
#                 eachpar : initial or present values of each                  #
#                           baseline parameter (either one value or as many    #
#                           as the number of transitions in 'tmat')            #
#   - cens      : the censoring time distributions. A list with components     #
#                 dist : the name of the censoring distributions (one value)   #
#                 eachpar : each censoring distribution parameter              #
#                           (either one value or as many as the number of      #
#                            possible starting states in 'tmat')               #
#                 admin: the time of administrative censoring                  #
#   - copula    : the copula model. A list with components                     #
#                 name : the name of the copula                                #
#                 par  : the copula parameter                                  #
#                                                                              #
#                                                                              #
#   Date: February, 20, 2012                                                   #
#   Last modification on: March, 30, 2012                                      #
################################################################################

thisState.tune <- function(target,
                           data,
                           atState,
                           eta,
                           tmat,
                           clock,
                           marg,
                           cens,
                           copula
                           ){
  ### - PREPARATION - ##########################################################
  if (is.null(atState))
    atState <- colnames(tmat)[which(colSums(tmat, na.rm=TRUE) == 0)]
  
  # All possible conditioning transitions
  inTrans <- tmat[which(!is.na(tmat[, atState])), atState]
  if (!length(inTrans))
    inTrans <- NULL
  
  # All possible outgoing transitions
  outTrans <- tmat[atState, which(!is.na(tmat[atState, ]))]
  # if ending state, then return results
  if (length(outTrans) == 0)
    return(list(marg=marg, cens=cens))

  # Reparameterization on all R
  optimPars <- list(
    marg = t(apply(marg[outTrans, names(marg) != "dist", drop=FALSE], 1, 
                   function(x) attr(eval(parse(text=as.character(marg[1, "dist"]))),
                                    "optimPars")(x))),
    cens = t(apply(cens$f[atState, names(cens$f) != "dist", drop=FALSE], 1, 
                   function(x) attr(eval(parse(text=as.character(cens$f[1, "dist"]))),
                                    "optimPars")(x))))    
  ################################################### - END of PREPARATION - ###

  parnames <- unique(c(colnames(optimPars$marg), colnames(optimPars$cens)))
  
  for (par in parnames) {
    inipar <- c(optimPars$marg[, par], optimPars$cens[, par])
    
    tomin <- function(x, par.name, outTrans, data, atState, eta,
                      tmat, clock, marg, cens, copula, target) {
      margdist <- attr(eval(parse(text=as.character(marg[1, "dist"]))), 
                       "optimPars")
      marg[outTrans, par.name] <-
        margdist(x[1:length(outTrans)], inv=TRUE)
      cmarg <- apply(marg, 1, extractMargs)        
        names(cmarg) <- 1:length(cmarg)
      ccens <- list(f = apply(cens$f, 1, extractMargs),
                    admin = cens$admin)
      
      criterion(data=data, atState=atState, eta=eta,
                tmat=tmat, clock=clock, 
                marg=cmarg,
                cens=ccens,
                copula=copula, target=target,
                verbose=FALSE)
    }
    
    cat(paste("Starting optimisation of parameter ",
              par,
              " at ", format(Sys.time(), "%X"),
              " (", Sys.Date(), ")", sep=""))
    optime <- system.time({par.res <- optim(inipar, 
                                  fn=tomin, 
                                  par.name=par, outTrans=outTrans,
                                  data=data, atState=atState, eta=eta,
                                  tmat=tmat, clock=clock, marg=marg, cens=cens,
                                  copula=copula, target=target,
                                  method="SANN")})
    cat(paste("\n Execution time: ", sec2ext(optime[1]), "\n", sep=""))
    
    marg[outTrans, par] <- attr(eval(parse(text=as.character(marg[1, "dist"]))), 
                                "optimPars")(par.res$par[outTrans], inv=TRUE)
    cens$f[atState, par] <- attr(eval(parse(text=as.character(cens$f[1, "dist"]))), 
                                  "optimPars")(par.res$par[-outTrans], inv=TRUE)
  }
  
  ### - COMPETING RISKS PARAMETERS TUNING - ####################################
  
  ######################################################## - END of TUNING - ###
  
  
  ### - UPDATE PARAMETERS - ####################################################
#   data[, sapply(c(".time", ".status"), function(x)
#     paste("tr", outTrans, x, sep=""))] <-
#       t(apply(cbind(data[, paste("tr", outTrans, ".time", sep="")],
#                     C.time=C.time), 1, function(x)
#                       c(rep(min(x), length(x)-1), 
#                         1:(length(x)-1) == which.min(x))))
  ############################################# - END of UPDATE PARAMETERS - ###
  
  
   return(list(marg=marg, cens=cens))
 }
