################################################################################
#  Simulation of frailty term                                                  #
################################################################################
#                                                                              #
#  Simulates the frailty term for a frailty multi-state simulation model       #
#                                                                              #
#  Its parameters are                                                          #
#   - dist      : the name of the fraty distribution                           #
#   - theta     : the frailty parameter                                        #
#   - nsim      : the number of subjects to simulate                           #
#   - nclus     : the number of clusters to simulate                           #
#   - csize     : the size(s) of cluster                                       #
#                                                                              #
#                                                                              #
#   Date: February, 13, 2012                                                   #
#   Last modification on: February, 13, 2012                                   #
################################################################################

simFrail <-function(dist="gamma", 
                    theta=.5, 
                    nsim=NULL, 
                    nclus=NULL, 
                    csize=NULL) {
  ##############################################################################
  ### * CONTROLS * #############################################################
  ##############################################################################
  ################################################
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
    if (!is.null(csize))) {
      if (length(csize) == 1)
        stop("Too few information (nsim, nclus, csize) to simulate frailties!")
      nsim <- sum(csize)
      nclus <- length(csize)
    }
    #[3.b] ONLY nsim GIVEN
    else if (!is.null(nsim))) {
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

  # In the cases [2.c] and [3.b]
  if (is.null(csize)) {
    csize <- runif(nclus)
    csize <- round(csize / sum(csize) * nsim)
    csize[nclus] <- nsim - sum(csize[-nclus])
    warning("\n Randomly generated clusters' sizes 'cszie'")
  }
  # last extreeeeme control
  if(nsim ! =sum(csize) || length(csize) != nclus)
    stop("\n --- something went wrong :( ---\n")
  ##############################################################################
  ###################################################### * END of CONTROLS * ###
  ##############################################################################

  res <- NULL
  attributes(res)$nsim  <- nsim
  attributes(res)$nclus <- nclus
  attributes(res)$csize <- csize
  
  if (substr(dist, 1, 3)=="gam") {
    z <- rgamma(nclus,  shape=1/theta, scale=theta)
    z <- as.vector(unlist(apply(cbind(z,csize), 1, function(x) rep(x[1],x[2]))))
    res$z <- z
    res$Cluster <- as.factor(unlist(apply(cbind(1:nclus,csize), 1, 
                                          function(x) rep(x[1],x[2]))))
  } 
  else
    stop(paste("Unknown frailty distribution '", dist, "'!\n", sep=""))
  
  return(as.data.frame(res))
}
