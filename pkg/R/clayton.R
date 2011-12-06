clayton <-
function(theta=1,
                    prev =NULL,
                    nsim =NULL,
                    marg =NULL,
                    pars =NULL,
                    eta  =NULL,
                    clock=NULL) {
  if (is.null(clock))
    clock <- "r"
  if (!(substring(clock, 1, 1) %in% c("r","f")))
    stop("Clock parameter 'clock' must be either 'f' (forward) or 'r' (reset)!")
  
  # starting state or conditioned?
  if (is.null(prev)){
    if (is.null(nsim)) 
      stop(paste("\nProvide either the conditioning times 'prev'",
                 "or the number of subject to simulate 'nsim'!\n"))
    k <- 1
  } else {
    if (!is.null(nsim))
      stop(paste("\nProvide only one betwen the conditioning times 'prev'",
                 "and the number of subject to simulate 'nsim'!\n"))
    prev <- as.matrix(prev)
    k<-ncol(prev)+1
    nsim<-nrow(prev)
  }
  
  marg <- eval(parse(text=paste(marg, "(pars=c(", 
                                paste(pars, collapse=","), "))",sep=""  )))
  
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
  if (substring(clock, 1, 1) == "f" && !is.null(prev))
    u <- u * margEta(prev, inv=FALSE)
  
 
  # simulation of new times
  if (k==1) { tk <- margEta(u, inv=TRUE)} else {
    tk <- margEta(
      (1+(u^(theta/(theta*(1-k)-1))-1)*(
        1+apply(cbind(prev, eta), 1,
                function(x) sum(margCovFrail(marg, lp=x[2])(x[1], inv=FALSE)^(-theta)-1))
      ))^(-1/theta),
    inv=TRUE)
  }
  
  if (substring(clock, 1, 1) == "r" && !is.null(prev))
    tk <- prev + tk
  
    tk <- matrix(as.numeric(tk), 
    dimnames=list(NULL,paste("T", k, sep="")))
  return(tk)
}
