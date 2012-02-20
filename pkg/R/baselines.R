################################################################################
#  Marginal baseline parametric distributions                                  #
################################################################################
#                                                                              #
#  Gives the survival function and its converse (quantile function)            #
#  for a given parametric marginal haard                                       #
#  Possible distributions:                                                     #
#   - gompertz                                                                 #
#   - loglogistic                                                              #
#   - lognormal                                                                #
#   - weibull                                                                  #
#                                                                              #
#  The argument is                                                             #
#   - pars : list of names parameters                                          #
#                                                                              #
#  The returned value is a function with arguments                             #
#   - x   : either the percentile (if inv=TRUE) or the time value (inv=FALSE)  #
#   - inv : inverted survival function? if TRUE, the quantile function         #
#           is obtained, if FALSE, the Survival function is.                   #
#                                                                              #
#                                                                              #
#   Date: February, 14, 2012                                                   #
#   Last modification on: February, 14, 2012                                   #
################################################################################


## - GOMPERTZ - ################################################################
gompertz <- function(pars=list(lambda=1, gamma=1)) {
  resf <- function(x, inv=FALSE) {
    if (!is.logical(inv)) 
      stop("Parameter inv must be logical.\n")
    if (inv) {
      if (x < 0 || x > 1)
        stop("With inv=TRUE, argument 'x' must be between 0 and 1!\n")
      c(T = log(1 - log(x) * pars[[2]] / pars[[1]]) / pars[[2]])
    } else {
      if (x < 0) 
        stop("With inv=FALSE, argument 'x' must be greater than 0!\n")
      c(S = exp(-pars[[1]] / pars[[2]] * (exp(x * pars[[2]]))))
    }
  }
  return(resf)
}

attr(gompertz, "optimPars") <- function(pars, inv=FALSE) {
  if(inv) {
    pars <- exp(pars)
  } else
    pars <- log(pars)
  return(pars)
}
######################################################### - END of GOMPERTZ - ##


## - LOGLOGISTIC - #############################################################
loglogistic <- function(pars=list(alpha=1, kappa=1)) {
  resf <- function(x, inv=FALSE) {
    if (!is.logical(inv)) 
      stop("Parameter inv must be logical.\n")
    if (inv){
      if (x < 0 || x > 1) 
        stop("With inv=TRUE, argument 'x' must be between 0 and 1!\n")
      c(T = (-(1 - 1 / x) * exp(-pars[[1]]))^(1 / pars[[2]]))
    } else {
      if (x<0) 
        stop("With inv=FALSE, argument 'x' must be greater than 0!\n")
      c(S = 1 / (1 + exp(pars[[1]]) * x^pars[[2]]))
    }
  }
  return(resf)
}

attr(loglogistic, "optimPars") <- function(pars, inv=FALSE) {
  if(inv) {
    pars$kappa <- exp(pars$kappa)
  } else
    pars$kappa <- log(pars$kappa)
  return(pars)
}
###################################################### - END of LOGLOGISTIC - ##



## - LOGNORMAL - ###############################################################
lognormal <- function(pars=list(mu=1, sigma=1)) {
  resf <- function(x, inv=FALSE) {
    if (!is.logical(inv)) 
      stop("Parameter inv must be logical.\n")
    if (inv) {
      if (x < 0 || x > 1) 
        stop("With inv=TRUE, argument 'x' must be between 0 and 1!\n")
      c(T = exp( pars[[1]] + pars[[2]] * qnorm(1 - x)))
    } else {
      if (x < 0)
        stop("With inv=FALSE, argument 'x' must be greater than 0!\n")
      c(S = 1 - pnorm((log(x) - pars[[1]]) / pars[[2]]))
    }
  }
  return(resf)
}

attr(lognormal, "optimPars") <- function(pars, inv=FALSE) {
  if(inv) {
    pars$sigma <- exp(pars$sigma)
  } else
    pars$sigma <- log(pars$sigma)
  return(pars)
}
######################################################## - END of LOGNORMAL - ##


## - WEIBULL - #################################################################
weibull <- function(pars=list(lambda=1, rho=1)) {
  if (length(pars) != 2) 
    stop("Weibull distribution needs 2 parameters!")
  resf <- function(x, inv=FALSE) {
    if (!is.logical(inv)) 
      stop("Parameter inv must be logical.\n")
    if (inv) {
      if (x < 0 || x > 1)
        stop("With inv=TRUE, argument 'x' must be between 0 and 1!\n")
      c(T = (-log(x) / pars[[1]])^(1 / pars[[2]]))
    } else {
      if (x < 0)
        stop("With inv=FALSE, argument 'x' must be greater than 0!\n")
      c(S = exp(-pars[[1]] * x^pars[[2]]))
    }
  }
  return(resf)
}

attr(weibull, "optimPars") <- function(pars, inv=FALSE) {
  if(inv) {
    pars <- exp(pars)
  } else
    pars <- log(pars)
  return(pars)
}
########################################################## - END of WEIBULL - ##

