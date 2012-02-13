clayton <- function(theta=1,
                    cond =NULL,
                    prev =NULL,
                    nsim =NULL,
                    marg =NULL,
                    pars =NULL,
                    eta  =NULL,
                    clock=NULL) {
  if (is.null(clock)) {
    clock <- "r"
    warning("Parameter 'clock' is not set! Fixed to 'reset'.")
  }
  else if (!(substring(clock, 1, 1) %in% c("r","f")))
    stop("Clock parameter 'clock' must be either 'f' (forward) or 'r' (reset)!")
  
  # NEITHER CONDITIONING NOR PREVIOUS TIMES
  if (is.null(cond) && is.null(prev)){
    if (is.null(nsim)) # nsim NEEDED!
      stop(paste("\nProvide either the number of subjects to simulate 'nsim'",
                 "or \n the conditioning times 'cond' and/or",
                 "the previous times 'prev'!\n"))
    else
      k <- 1
  }
  # EITHER CONDITIONING OR PREVIOUS TIMES
  else {
    if (!is.null(nsim)) # nsim CANNOT BE ALTERED!
      stop(paste("\nProvide only one between",
                 "the number of subjects to simulate 'nsim'",
                 "and \n the conditioning times 'cond' and/or",
                 "the previous times 'prev'!\n"))
    # WITH CONDITIONING PREVIOUS TIMES
    if (!is.null(cond)) {
      cond <- as.matrix(cond)
        if (ncol(cond)!=1)
          stop(paste("\nThe conditioning times 'cond' must be either a vector",
                     "or a one-column matrix!"))
      nsim <- nrow(cond)
    }
    # WITH PREVIOUS TIMES
    if (!is.null(prev)) {
      prev <- as.matrix(prev)
      k    <- ncol(prev) + 1
      nsim <- nrow(prev)      
    }
    # BOTH CONDITIONING AND PREVIOUS TIMES
    if (!is.null(cond) && !is.null(prev)) {
      if (nrow(cond) != nrow(prev))
        stop(paste("The number of previous times 'prev' (",
                   nrow(prev),   ") and that of conditioning times 'cond' (",
                   length(cond), ") do not match!", sep=""))
    }
  }
  
  marg <- eval(parse(text=paste(marg, "(pars=c(", 
                                paste(pars, collapse=","), "))",sep=""  )))
  # problem: we would need parameters of all previous and conditioning
  # transitions to be used in the denominator of the conditional survival

  
  # in case of covariates and/or frailties
  if (!(is.null(nsim) || is.null(eta)))
    if (nsim!=length(eta))
      stop (paste("The number of simulations 'nsim' does not match",
                  "the number of simulations of the covariates in 'eta'!"))

  if (is.null(eta))
    eta <- rep(0, nsim)
  margEta <- margCovFrail(marg, lp=eta)
  
  # simulation of uniforms
  u <- runif(nsim)

# - here problems!!! - #
  if (substring(clock, 1, 1) == "f" && !is.null(cond))
    u <- u * ( 1 + (margEta(cond, inv=FALSE)^(-theta)-1) / (
      1 + apply(cbind(cond, eta), 1,
                function(x) sum(margEta(x, inv=FALSE)^(-theta)-1))
      ))^(1-k-1/theta)

margEta(cond, inv=FALSE) #########
  
 
  # simulation of new times
  if (k==1) { tk <- margEta(u, inv=TRUE) } else {
    tk <- margEta(
      (1+(u^(theta/(theta*(1-k)-1))-1)*(
        1+apply(cbind(cond, eta), 1,
                function(x) sum(margCovFrail(marg, lp=x[2])(x[1], inv=FALSE)^(-theta)-1))
 # not correct! marginal marg, should use parameters of each previous marginal
      ))^(-1/theta),
    inv=TRUE)
  }
  
  if (substring(clock, 1, 1) == "r" && !is.null(cond))
    tk <- cond + tk
  
    tk <- matrix(as.numeric(tk), 
    dimnames=list(NULL,paste("T", k, sep="")))
  return(tk)
}
