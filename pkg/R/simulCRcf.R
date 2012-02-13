simulCRcf <- function(nsim  = NULL,
                      trans = NULL,
                      clock = "reset",
                    # Frailty
                      fdist = "gamma", 
                      ftheta= .5,
                      nclus = NULL, 
                      csize = NULL,
                    # Covariates
                      covs,
                      beta,
                    # Marginals
                      ctheta= 1, 
#                       prev = NULL, 
                      marg  = "weibull", 
                      pars  = c(lambda=1, rho=1), 
                      cens  = "weibull", 
                      cpars = c(lambda=1, rho=1), 
                      adcens= 72) {
  res <- NULL

  # number of simulations
  if (is.null(nsim))
    stop("Number of simulations 'nsim' not defined!")
  if (!(is.numeric(nsim) && length(nsim) == 1))
    stop(paste("The number of simulations 'nsim' must be",
               "a one-dimensional numeric object!"))
  
  # number of time variables
  if (!is.matrix(pars)) {
    if (is.vector(pars)) pars <- matrix(pars) else
      stop(paste("The parameters object 'pars' must be a matrix",
        "with one column for each time variable",
        "and a row for each distribution parameter!"))
  }
  nt <- ncol(pars)
  
  # number of transitions and CRs blocks
  if (max(trans, na.rm=TRUE) != nt)
    stop(paste("The number of transitions in the transition matrix 'trans'",
               "does not match with the number of columns of 'pars'!"))
  if (length(which(rowSums(trans, na.rm=TRUE) > 0)) != ncol(cpars))
    stop(paste("The number of competing events blocks",
               "(number of rows of 'trans' with non NA elements)",
               "does not match with the number of columns of 'cpars'!"))
 
  
  # ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## #
  # Simulation of frailties + covariates and computation of linear predictors  #
  # ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## #
 
  # Frailty simulation
  if (is.null(fdist) && is.null(ftheta)) {
    X.frail <- NULL
    eta.frail <- rep(0, nsim)
  }
  else {
    if (is.null(fdist))
      stop(paste("No frailty distribution 'fdist' is specified,",
                 "while its dispersion 'ftheta' is!"))
    if (is.null(ftheta))
      stop(paste("No frailty dispersion 'ftheta' is specified,",
                 "while frailty distribution 'fdist' is!"))
    X.frail <- simFrail(dist=fdist, theta=ftheta,
                      nsim=nsim, nclus=nclus, csize=csize)
    eta.frail <- log(X.frail[, "z"])
    attributes(res)$frailty <- list(dist = fdist, frail = X.frail)
  }

  # Covariates simulation
  if (is.null(beta) && is.null(covs)) {
    X.cov <- NULL
    eta.cov <- rep(0, nsim)
  }
  else {
    if (is.null(beta))
      stop(paste("No regression coefficients 'beta' are specified,",
                 "while covariates 'covs' are!"))
    if (is.null(covs))
      stop(paste("No covariates 'covs' are specified,",
                 "while regression coefficients 'beta' are!"))
    if (length(beta) != length(covs))
      stop(paste("\nThe number of covariates (", length(covs),
                 ") is different from the number of regression coefficients (",
                 length(beta), ")!\n", sep=""))
    if (diff(range((sapply(beta, length)))) != 0)
      stop("All beta coefficients in 'beta' must have the same length!")
    if (length(beta[[1]]) != nt)
      stop(paste("The length of each beta coefficient in 'beta' must be",
                 "the same as the length of each parameter in 'pars'!"))
    X.cov <- simCov(covs=covs, nsim=nsim)
    eta.cov <- as.matrix(X.cov) %*%
      t(matrix(unlist(beta), nt, dimnames=list(NULL, names(beta))))
    attributes(res)$covariates <- list(dist = covs, beta = beta)
  }
 
  # Merge frailties and covariates
  eta <- eta.cov + eta.frail
  if (!(is.null(X.frail) || is.null(X.cov)))
    Xdata <- cbind(X.frail, X.cov)
  else if (is.null(X.frail))
    Xdata <- X.cov
  else Xdata <- X.frail


  # ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## #
  # Simulation of first blocks                                                 #
  # ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## #

  ss <- which(colSums(trans, na.rm=TRUE) == 0)
  if (length(ss) > 1)
    stop(paste("This method is implmented for multi-state structures",
               "with only one starting state!"))
  
  dest.states <- which(!is.na(trans[ss, ]))
  transitions <- trans[ss, dest.states]
  res <- compTimes(claytonCRBlock(nsim=nsim, theta=ctheta, 
                                  marg=marg, pars=pars[, transitions],
                                  cens=cens,  cpars=cpars[, ss],
                                  adcens=adcens, eta=eta,
                                  names=colnames(trans)[dest.states],
                                  clock=clock))
  # Present states after the first transition
  pres.state <- apply(res[, (ncol(res)/2+1):ncol(res)], 1, function(x) {
    ifelse(sum(x) != 0, names(x)[which(x==1)], NA)
  })
  
  # ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## #
  # Simulation of following blocks                                             #
  # ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## #
  
  repeat {
    # the states in which subjects are 
    # and from which a new transition is possible
    startstates <- which(sapply(rownames(trans), function(x) {
        x %in% names(table(pres.state)) && sum(trans[x,], na.rm=TRUE)>0}))
    if (length(startstates) == 0) 
      break
    for (ss in startstates) {
      # ID numbers of patients who are at state ss
      whichss <- which(pres.state == rownames(trans)[ss])
      # ID numbers of states to which they can move
      dest.states <- which(!is.na(trans[ss, ]))
      # ID numbers of transitions from ss to dest.states
      transitions <- trans[ss, dest.states]
      
      res2 <- compTimes(claytonCRBlock(prev=res[whichss, paste("T",
                                         rownames(trans)[ss], sep="")],
                                       theta=ctheta, 
                                       marg=marg, pars=pars[, transitions],
                                       cens=cens,  cpars=cpars[, ss],
                                       adcens=adcens, 
                                       eta=matrix(eta[whichss,], length(whichss)),
                                       names=colnames(trans)[dest.states],
                                       clock=clock))
      for (ds in dest.states) {
        # possible addition of new state variables
        if (! colnames(trans)[ds] %in% colnames(res))
          res <- cbind(res, matrix(cbind(res[, paste("T",
                                         rownames(trans)[ss], sep="")], 0),
                                   nrow(res), dimnames=list(NULL, 
            c(paste("T", colnames(trans)[ds], sep=""), colnames(trans)[ds]))))
        
        res[whichss, paste(c("T",""), colnames(trans)[ds], sep="")] <-
          res2[    , paste(c("T",""), colnames(trans)[ds], sep="")]
      }

      # New present states
      pres.state[whichss] <- apply(
        matrix(res2[, (ncol(res2)/2+1):ncol(res2)], length(whichss)), 
        1, function(x) {
          ifelse(sum(x) != 0, names(dest.states)[which(x==1)], NA) })
    }
  }
  
  
  return(as.data.frame(res))
  
}
