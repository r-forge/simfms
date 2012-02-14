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
#   - nclus     : the number of clusters to simulate                           #
#   - csize     : the size(s) of cluster                                       #
#   - covs      : the covariates to simulate. A list with components           #
#                 nameOfCovariate = simulationfunction                         #
#   - beta      : the regression coefficients for covariates. A list of the    #
#                 same length as 'covs', with elements                         #
#                 nameOfCovariate = a vector of the same length as the number  #
#                 of transitions in 'tmat'                                     #
#   - marg      : the marginal baseline hazards. A list with components        #
#                 dist    : the name of the baseline hazard distribution       #
#                           (either one value or as many as the number         #
#                            of transitions in 'tmat')                         #
#                 eachpar : each baseline parameter                            #
#                           (either one value or as many as the number         #
#                            of transitions in 'tmat')                         #
#   - cens      : the censoring time distribution. A list with components      #
#                 dist : the name of the censoring distribution                #
#                 par  : the vector of the censoring distribution parameters   #
#                 admin: the time of administrative censoring                  #
#   - copula    : the copula model. A list with components                     #
#                 name : the name of the copula                                #
#                 par  : the copula parameter                                  #
#                                                                              #
#                                                                              #
#   Date: February, 13, 2012                                                   #
#   Last modification on: February, 13, 2012                                   #
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
                   copula= list(name="clayton",
                                par= 1)
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
    warning("Parameter 'clock' is not set! Fixed to 'forward'.")
  }
  if (!(substring(clock, 1, 1) %in% c("r","f")))
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
  
  ### Covariates and regression parameters ###
  if (!is.list(covs))
    stop("The covariates object 'covs' must be a list!")
  if (!is.list(beta))
    stop("The regression coefficient object 'beta' must be a list!")
  if (length(beta) != length(covs))
    stop(paste("\nThe number of covariates (", length(covs),
               ") is different from the number of regression coefficients (",
               length(beta), ")!\n", sep=""))
  if (sum(sort(names(beta)) != sort(names(covs))) > 0)
    stop(paste("\nThe covariates names in 'covs' ",
               "are different from those in 'beta'!\n",
               " names(covs):  ", paste(sort(names(covs)), collapse=", "),
               "\n names(beta):  ", paste(sort(names(beta)), collapse=", "),
               sep=""))
  if (!is.null(beta)){
    if (diff(range((sapply(beta, length)))) != 0)
      stop("All beta coefficients in 'beta' must have the same length!")
    if (length(beta[[1]]) != max(tmat, na.rm=TRUE))
      stop(paste("The length of each beta coefficient in 'beta' must be",
                 "the number of transitions in 'tmat'!"))
  }
  
  ### Marginals ###
  validbaselines <- c("gompertz", "loglogistic", "lognormal", "weibull")
  if (!all(marg$dist %in% validbaselines))
    stop(paste("Baseline distributions not valid! \nThey must be one of '",
               paste(validbaselines, collapse="', '"),
               "'.", sep=""))
  
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
    warning(paste("\n Randomly generated number of clusters 'nclus':", 
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
    warning("\n Randomly generated clusters' sizes 'cszie'")
  }
  # last very last control
  if(nsim != sum(csize) || length(csize) != nclus)
    stop("\n --- something went wrong :( ---\n")
  ##############################################################################
  # ------------------------------ end of Sizes ------------------------------ #
  ##############################################################################

  return(list(nsim  = nsim, 
              nclus = nclus,
              csize = csize,
              clock = clock))
}

