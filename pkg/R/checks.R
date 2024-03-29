################################################################################
#  Preliminary controls on simulation parameters                               #
################################################################################
#                                                                              #
#  Performs checks on simuation parameters                                     #
#                                                                              #
#  Its parameters are                                                          #
#   - nsim      : the number of subjects to simulate                           #
#   - tmat      : the trnasitions matrix                                       #
#   - clock     : either 'forward' or 'reset'                                  #
#   - frailty   : the frailty term specifications. A list with components      #
#                 dist: the name of the frailty distribution                   #
#                 par : the frailty parameter(s) value                         #
#                 type: the frailty type: 'shared', 'iid' or 'nested'          #
#   - nclus     : the number of clusters to simulate                           #
#   - csize     : the size(s) of cluster                                       #
#   - covs      : the covariates to simulate. A list with components           #
#                 nameOfCovariate = simulationfunction                         #
#   - beta      : the regression coefficients for covariates. A list of the    #
#                 same length as 'covs', with elements                         #
#                 nameOfCovariate = a vector of the same length as the number  #
#                 of transitions in 'tmat'                                     #
#   - marg      : the marginal baseline hazards. A list with components        #
#                 dist : the name of the censoring distributions (one value)   #
#                 eachpar : each baseline parameter                            #
#                           (either one value or as many as the number         #
#                            of transitions in 'tmat')                         #
#   - cens      : the censoring time distributions. A list with components     #
#                 dist : the name of the censoring distributions (one value)   #
#                 eachpar : each censoring distribution parameter              #
#                           (either one value or as many as the number of      #
#                            possible starting states in 'tmat')               #
#   - copula    : the copula model. A list with components                     #
#                 name : the name of the copula                                #
#                 par  : the copula parameter                                  #
#                                                                              #
#                                                                              #
#   Date: February, 13, 2012                                                   #
#   Last modification on: March, 29, 2012                                      #
################################################################################

checks <- function(nsim,
                   tmat,
                   clock,
                   # Frailty
                   frailty,
                   nclus,
                   csize,
                   # Covariates
                   covs,
                   beta,
                   # Marginals
                   marg,
                   cens,
                   # Copula
                   copula,
                   # Data format
                   format
                   ) {
  ### Transition matrix ###
  if (is.null(tmat))
    stop("No transition matrix 'tmat' defined!")
  if (sum(colSums(tmat, na.rm=TRUE) == 0) != 1)
    stop("Simulation method implemented only for one starting state!")
  if (sum(1 - is.na(tmat[lower.tri(tmat)])) > 1)
    stop("Only for acyclic multi-state structures can be used!")
  
  ### Clock ###
  if (is.null(clock)) {
    clock <- "f"
    cat("\nParameter 'clock' is not set! Fixed to 'forward'.")
  }
  if (tolower(substring(clock, 1, 1)) == "r")
    clock <- "reset"
  else if (tolower(substring(clock, 1, 1)) == "f")
    clock <- "forward"
  else
    stop("Clock parameter 'clock' must be either 'f' (forward) or 'r' (reset)!")
  
  ### Frailties ###
  if (!is.list(frailty))
    stop("The frailty object 'frailty' must be a list!")
  if (is.null(frailty$dist))
    stop("No frailty distribution is specified!")
  if (is.null(frailty$par))
    stop("No frailty parameter is specified!")
  validfrailties <- c("gamma", "none")
  if (!(frailty$dist %in% validfrailties))
    stop(paste("Frailty distribution not valid! \nIt must be one of '",
               paste(validfrailties, collapse="', '"),
               "'.", sep=""))
  if (is.null(frailty$type)) {
    cat("\nNo frailty type specified! Set to 'shared'")
    Ftype <- "shared"
  } else {
    Ftype <- frailty$type
  }
  if (! substr(Ftype, 1, 1) %in% c("s", "i", "n")) {
    stop(paste("The frailties type must be",
               "'shared', 'iid' or 'nested'!"))
  }
  nFpars <- length(frailty$par)
  nFpars.possible <- unlist(list(s = 1,
                                 i = c(1, max(tmat, na.rm=TRUE)),
                                 n = c(2, max(tmat, na.rm=TRUE) + 1))[
                                   substr(Ftype, 1, 1)])
  if (! nFpars %in% nFpars.possible) {
    stop(paste("Wrong length of frailty parameters vector.",
               "It is", nFpars, "instead of", 
               paste(nFpars.possible, collapse=" or ")))
  }
  rm(nFpars.possible)
  
  ### Covariates and regression parameters ###
  if (!(is.list(covs) || is.null(covs)))
    stop("The covariates object 'covs' must be a list!")
  if (!(is.list(beta) || is.null(covs)))
    stop("The regression coefficient object 'beta' must be a list!")
  if (length(beta) != length(covs))
    stop(paste("\nThe number of covariates (", length(covs),
               ") is different from the number of regression coefficients (",
               length(beta), ")!\n", sep=""))
  if (!is.null(beta)) {
    if (any(sort(names(beta)) != sort(names(covs))))
      stop(paste("\nThe covariates names in 'covs' ",
                 "are different from those in 'beta'!\n",
                 " names(covs):  ", paste(sort(names(covs)), collapse=", "),
                 "\n names(beta):  ", paste(sort(names(beta)), collapse=", "),
                 sep=""))
    if (diff(range((sapply(beta, length)))) != 0)
      stop("All beta coefficients in 'beta' must have the same length!")
    if (length(beta[[1]]) != max(tmat, na.rm=TRUE))
      stop(paste("The length of each beta coefficient in 'beta' must be",
                 "the number of transitions in 'tmat'!"))
  }
  
  ### Marginals ###
  if (!is.list(marg) || is.null(marg$dist))
    stop(paste("The marginal baselines object 'marg' must be a list with",
               "at least an element 'dist'!"))
  if (length(marg$dist) != 1)
    stop (paste("The method is implemented only for marginal baselines with the",
                "same distribution!",
                "Please, set a single value for 'marg$dist'."))
  validbaselines <- c("gompertz", "loglogistic", "lognormal", "weibull")
  if (!(marg$dist %in% validbaselines))
    stop(paste("Baseline distributions not valid! \nThey must be one of '",
               paste(validbaselines, collapse="', '"),
               "'.", sep=""))
  for (marEl in names(marg)[names(marg) != "dist"]) {
    if (length(marg[[marEl]]) == 1)
      marg[[marEl]] <- rep(marg[[marEl]], max(tmat, na.rm=TRUE))
    else if (length(marg[[marEl]]) != max(tmat, na.rm=TRUE))
      stop(paste("Invalid length of marg$", marEl,
                 ": ", length(marg[[marEl]]),
                 " instead of ", max(tmat, na.rm=TRUE),
                 sep=""))
  }
  
  ### Censoring ###
  if (!is.list(cens) || is.null(cens$dist))
    stop(paste("The censoring distributions object 'cens' must be a list with",
               "at least an element 'dist'!"))
  if (length(cens$dist) != 1)
    stop (paste("The method is implemented only for censoring times with the",
                "same parametric distribution!",
                "Please, set a single value for 'cens$dist'."))
  if (!all(cens$dist %in% validbaselines))
    stop(paste("Censoring distributions not valid! \nThey must be one of '",
               paste(validbaselines, collapse="', '"),
               "'.", sep=""))
  if (is.null(cens$admin))
    cens$admin <- Inf
  
  censNames <- rownames(tmat)[rowSums(tmat, na.rm=TRUE) > 0]
  for (censEl in names(cens)[! names(cens) %in% c("dist", "admin")]) {
    if (length(cens[[censEl]]) == 1) {
      cens[[censEl]] <- rep(cens[[censEl]], length(censNames))
      names(cens[[censEl]]) <- censNames
    } else
      if (length(cens[[censEl]]) == length(censNames)) {
        if (is.null(names(cens[[censEl]]))) {
          names(cens[[censEl]]) <- censNames
        } else 
          if (!all(names(cens[[censEl]]) %in% censNames))
            stop(paste("The names of elements of 'cens$", censEl, 
                       "' do not match the names",
                       "of the non-final states in 'tmat'!",
                       "\nNon-final states in 'tmat':\n  ", 
                       paste(sort(censNames), collapse=", "),
                       "\nNames of cens$", censEl, "':\n  ", 
                       paste(sort(names(cens[[censEl]])), collapse=", "),
                       sep=""))                       
    } else
      stop(paste("Invalid length of cens$", censEl,
                 ": ", length(cens[[censEl]]),
                 " instead of ", sum(rowSums(tmat, na.rm=TRUE) > 0),
                 sep=""))
  }
  
  ### Copula ###
  validcopulas <- c("clayton")
  if (!(copula$name %in% validcopulas))
    stop(paste("Copula name not valid! \nIt must be one of '",
               paste(validcopulas, collapse="', '"),
               "'.", sep=""))
  
  
  
  ##############################################################################
  # --------------------------------- Sizes ---------------------------------- #
  ##############################################################################
  #[1]# ALL SIZES (nsim, nclus, csize) ARE GIVEN #
  if (is.null(nsim) + is.null(nclus) + is.null(csize) == 0) {
    if (!(length(csize) %in% c(1, nclus))) 
      stop("Number of clusters 'nclus' and cluster size(s) 'csize' incoherent!")
    if (length(csize) == 1)
      csize <- rep(csize, nclus)
    if (nsim != sum(csize)) 
      stop(paste("Number of subjects 'nsim' incoherent to the number",
                 "of clusters 'nclus' and cluster size(s) 'csize'!"))
  }
  ###############################################
  #[2]# TWO SIZES (nsim, nclus, csize) ARE GIVEN #
  else if (is.null(nsim) + is.null(nclus) + is.null(csize) == 1) {
    #[2.a]  MISSING nsim - GIVEN nclus AND csize
    if (is.null(nsim)) {
      if (length(csize) == 1) 
        csize <- rep(csize, nclus) 
      else if (length(csize) != nclus)
        stop("Number of clusters 'nclus' and cluster size(s) 'csize' incoherent!")
      nsim <- sum(csize)
      #[2.b] MISSING nclus - GIVEN nsim AND csize
    } else if (is.null(nclus)) {
      if (length(csize) == 1)
        stop(paste("No number of cluster 'nclus' provided, with a single value ",
                   "for the clusters length (", csize, ")!", sep=""))
      if (nsim / sum(csize) - nsim %/% sum(csize) != 0) 
        stop("Number of subjects 'nsim' and cluster size(s) 'csize' incoherent!")
      nclus <- length(csize)
    }
    #[2.c] MISSING csize - GIVEN nsim AND nclus
    # see later on
  }
  ##############################
  #[3]# ONLY ONE SIZE IS GIVEN #
  else if (is.null(nsim) + is.null(nclus) + is.null(csize) == 2) {
    #[3.a] ONLY csize VECTOR GIVEN
    if (!is.null(csize)) {
      if (length(csize) == 1)
        stop("Too few information (nsim, nclus, csize) to simulate frailties!")
      nsim <- sum(csize)
      nclus <- length(csize)
    }
  #[3.b] ONLY nsim GIVEN
  else if (!is.null(nsim)) {
    nclus <- round((runif(1)*nsim)^.5)
    cat(paste("\n Randomly generated number of clusters 'nclus':", 
                  nclus, "\n"))
  }
  #[3.c] ONLY nclus GIVEN
  else if (!is.null(nclus)) 
    stop("At least the number of simulations 'nsim' must be provided!")
  }
  ####################
  # NO INFO IS GIVEN #
  else
    stop("No information (nsim, nclus, csize) to simulate frailties is present!")
  ##############################################################################
  
  
  # In the cases [2.c] and [3.b]
  if (is.null(csize)) {
    csize <- runif(nclus)
    csize <- round(csize / sum(csize) * nsim)
    csize[nclus] <- nsim - sum(csize[-nclus])
    cat("\n\n Randomly generated clusters' sizes 'cszie'")
  }
  # last very last control
  if(nsim != sum(csize) || length(csize) != nclus)
    stop("\n --- something went wrong :( ---\n")
  ##############################################################################
  # ------------------------------ end of Sizes ------------------------------ #
  ##############################################################################

  if (! substr(format, 1, 1) %in% c("l", "w")) {
    cat("\nInvalid 'format' value. Set to 'long'.")
  }
  
  return(list(nsim  = nsim, 
              nclus = nclus,
              csize = csize,
              Ftype = Ftype,
              clock = clock,
              marg  = marg,
              cens  = cens))
}

