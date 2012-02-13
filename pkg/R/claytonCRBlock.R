claytonCRBlock <- function(theta  =1,
                           cond  =NULL,
                           nsim  =NULL, 
                           marg  =NULL,
                           pars  =NULL,
                           cens  =NULL, 
                           cpars =NULL, 
                           adcens=NULL,
                           eta   =NULL,
                           names =NULL,
                           clock =NULL) {
  if (is.null(clock)) {
    clock <- "r"
    warning("Parameter 'clock' is not set! Fixed to 'reset'.")
  }
  if (!(substring(clock, 1, 1)%in%c("r","f")))
    stop("Clock parameter 'clock' must be either 'f' (forward) or 'r' (reset)!")
 
  if (substring(clock, 1, 1) == "f" && !is.null(cond) && !is.null(adcens))
    adcens <- adcens - cond
  
  # starting CRs block or conditioned?
  if (is.null(cond)){
    if (is.null(nsim))
      stop(paste("\nProvide either the conditioning times 'cond'",
                 "or the number of subject to simulate 'nsim'!\n"))
  } else {
    if (!is.null(nsim))
      stop(paste("\nProvide only one betwen the conditioning times 'cond'",
                 "and the number of subject to simulate 'nsim'!\n"))
    else {
      cond <- as.matrix(cond)
      if (ncol(cond)!=1)
        stop(paste("\nThe conditioning times 'cond' must be either a vector",
                   "or a one-column matrix!"))
    }
  }
  # Administrative Censoring
  if (!is.null(adcens)) if(min(adcens) <= 0)
    stop("Administrative censoring time(s) 'adcens' must be positive!")
  # number of time variables in the block
  if (!is.matrix(pars)) {
    if (is.vector(pars)) pars <- matrix(pars) else
      stop(paste("The parameters object 'pars' must be a matrix",
        "with one column for each time variable",
        "and a row for each distribution parameter!"))
  } 
  nt <- ncol(pars)
  
  # inizialisation of 'rest'
  rest <- matrix(NA, max(nrow(cond),nsim), nt)
                   
  # variables' names
#   if (is.null(cond))
#     k <- 0
#   else
#     k <- max(as.numeric(substr(colnames(cond), 2, 5)), na.rm=TRUE)
  if (is.null(names)) 
#     names <- k+(1:nt) 
    names <- 1:nt
  colnames(rest) = paste("T", names, sep="")
                   
  # Censoring
  if (!is.null(adcens) || !is.null(cens)) {
    rest <- cbind(rest, C=NA)
    if (!is.null(cens)) {
      cens <- eval(parse(text=paste(cens, "(pars=c(",
                   paste(cpars, collapse=","), "))",sep=""  )))
      uC <- runif(max(nrow(cond),nsim))
      if (substring(clock, 1, 1) == "f" && !is.null(cond))
        uC <- uC * cens(cond, inv=FALSE)
      rest[,"C"] <- cens(uC, inv=TRUE )
    }

    if (substring(clock, 1, 1) == "r" && !is.null(cond))
      rest[,"C"] <- cond + rest[,"C"]
   
    if (!is.null(adcens))
      rest[,"C"] <- apply(cbind(rest[,"C"], adcens), 1, min, na.rm=TRUE)
  }
  
  # first variable of the block
    rest[, 1] <- clayton(theta=theta, nsim=nsim, cond=cond,
                         marg=marg, pars=pars[, 1], eta=eta[, 1], clock=clock)
    
  # following variables of the block
  if (nt>1)
    for (i in 2:nt) {
      rest[, i] <- clayton(theta=theta, cond=cond, prev=rest[, 1:(i-1)],
                           marg=marg, pars=pars[, i], eta=eta[, i], clock=clock)
    }

  return(rest)
}
