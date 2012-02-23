################################################################################
#  Clayton copula conditional quantile function                                #
################################################################################
#                                                                              #
#  Computes the simulated time for the k-th time variable, given               #
#     a possible conditioning transition time and distribution                 #
#     a possible set of previously simulated competing transitions             #
#       times and distributions                                                #
#                                                                              #
#  Its parameters are                                                          #
#   - par       : the copula parameter par                                     #
#   - condTime  : the time of the possible conditioning transition,            #
#                 with the transition number as name                           #
#   - condMarg  : the marginal baseline function of the possible               #
#                 conditioning transition,                                     #
#                 with the transition number as name                           #
#   - prevTimes : the times of the possible previous transitions,              #
#                 with the transition numbers as names                         #
#   - prevMargs : the marginal baseline functions of the possible              #
#                 previous transitions,                                        #
#                 with the transition numbers as names                         #
#   - marg      : the marginal baseline function of the transition to simulate #
#   - eta       : the vector of the linear predictors, of length k,            #
#                 in the order:                                                #
#                   1 (possibly) for the conditioning transition               #
#                   those of the (possible) previous transitions               #
#                   1 for the present transition                               #
#   - clock     : either 'forward' or 'reset'                                  #
#                                                                              #
#                                                                              #
#   Date: February, 14, 2012                                                   #
#   Last modification on: February, 16, 2012                                   #
################################################################################

clayton <- function(par,
                    subj,
                    atState,
                    inTrans,
                    outTrans,
                    trans,
                    data,
                    eta,
                    tmat,
                    clock,
                    marg,
                    cens) {
  # Present state and Conditioning transition infos
  if (is.null(inTrans)){ # from the starting state
    condTime <- 
      condMarg <- NULL
  } else { # from all the other states
    condTime <- data[subj, paste("tr", inTrans, ".time", sep="")]
    condMarg <- extractMargs(as.data.frame(marg)[inTrans,])
  }
  
  ot.N <- which(outTrans == trans) # ot.N its rank in the CRs block!!!!!!!!!!!!!!
  # Previous transition(s) infos
  if (ot.N == 1) {
    prevTimes <- 
      prevMargs <- NULL
  } else {
    prevOTs <- outTrans[1:(ot.N - 1)]
    prevTimes <- data[subj, paste("tr", prevOTs, ".time", sep=""),
                      drop=FALSE]
    prevMargs <- apply(as.data.frame(marg)[prevOTs, ], 1, extractMargs)
  }
  
  # Marginal of present transition
  this.marg <- extractMargs(as.data.frame(marg)[trans, ])
  
  
# clayton <- function(par,
#                     condTime  = NULL,
#                     condMarg  = NULL,
#                     prevTimes = NULL,
#                     prevMargs = NULL,
#                     marg = NULL,
#                     eta = NULL,
#                     clock = NULL) {
  k <- 1 + length(condTime) + length(prevTimes)
  prevTimes <- c(condTime, prevTimes)
  prevMargs <- c(condMarg, prevMargs)
  
  ### DENOMINATOR: 1 + sum_{j=1}^{k-1}[ S_j(t_j)^(-th exp(eta_j)) - 1] #########
  denom <- 2 - k 
  if (length(prevTimes)) {
    denom <- denom +
      sum(sapply(1:length(prevMargs),
                 function(x) prevMargs[[x]](prevTimes[[x]])^(
                   - par * exp(eta[x])),
                 USE.NAMES=FALSE))
  }
  ####################################################### END of DENOMINATOR ###
  
  ### CLOCK FORWARD CORRECTION #################################################
  if (clock == "forward" && !is.null(condTime)) {
#     clock <- (1 + marg(condTime)^(-par * exp(eta[ncol(eta)])) / denom)^(
    clock <- (1 + this.marg(condTime)^(-par * exp(eta[ncol(eta)])) / denom)^(
        1 - k - 1 / par)
  }
  else
    clock <- 1
  ########################################## END of CLOCK FORWARD CORRECTION ###
      
  u <- runif(n=1, min=0, max=1)
  arg <- (1 + denom * ((u * clock)^(1 / (1 - k - 1 / par)) - 1))^(
    - 1 / (par * exp(eta[ncol(eta)])))
#   T <- marg(arg, inv=TRUE)
  T <- this.marg(arg, inv=TRUE)

  ### CLOCK RESET CORRECTION ###################################################
  if (clock == "reset" && !is.null(condTime)) {
    T <- condTime +  T
  }
  ############################################ END of CLOCK RESET CORRECTION ###
  
#   return(marg(arg, inv=TRUE))
  return(T)
}


