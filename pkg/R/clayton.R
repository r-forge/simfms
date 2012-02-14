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
#   - theta     : the copula  parameter                                        #
#   - condTime  : the time of the possible conditioning transition,            #
#                 with the transition number as name                           #
#   - condMarg  : the marginal baseline function of the possible               #
#                 conditioning transition,                                     #
#                 with the transition number as name                           #
#   - prevTime  : the times of the posssible previous transitions,             #
#                 with the transition numbers as names                         #
#   - prevMarg  : the marginal baseline functions of the possible              #
#                 conditioning transitions,                                    #
#                 with the transition numbers as names                         #
#   - marg      : the marginal baseline function of the transition to simulate #
#   - eta       : the vector of the linear predictors, of length k,            #
#                 k = 1 + the number of previous transitions                   #
#                       + 1 if there is a conditioning transition              #
#   - clock     : either 'forward' or 'reset'                                  #
#                                                                              #
#                                                                              #
#   Date: February, 14, 2012                                                   #
#   Last modification on: February, 14, 2012                                   #
################################################################################

clayton <- function(theta = 1,
                    condTime  = NULL,
                    condMarg  = NULL,
                    prevTimes = NULL,
                    prevMargs = NULL,
                    marg = NULL,
                    eta = NULL,
                    clock = NULL) {
  k <- 1 + length(condTime) + length(prevTimes)
  prevTimes <- c(condTime, prevTimes)
  prevMargs <- c(condMarg, prevMargs)
  
  ### DENOMINATOR: 1 + sum_{j=1}^{k-1}[ S_j(t_j)^(-th exp(eta_j)) - 1] #########
  denom <- 2 - k 
  if (length(prevTimes))
    denom <- denom +
      sum(sapply(names(prevMargs), 
                 function(x) prevMargs[[x]](prevTimes[x])^(
                   - theta * exp(eta[x])),
                 USE.NAMES=FALSE))
  ####################################################### END of DENOMINATOR ###
  
  ### CLOCK FORWARD CORRECTION #################################################
  if (clock == "reset")
    clock <- 1
  else if (clock == "forward")
    clock <- (1 + marg(condTime)^(-theta * exp(eta[names(marg)])) / denom)^(
      1 - k - 1 / theta)
  ########################################## END of CLOCK FORWARD CORRECTION ###
      
  T <- function(u) {
    arg <- (1 + denom * ((u * clock)^(k - 1 + 1 / theta) - 1))^(
      - 1 / (theta * exp(eta[names(marg)])))
    return(marg(arg, inv=TRUE))
  }
  return(T)
}


