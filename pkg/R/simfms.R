################################################################################
#  Simulation of frailty multi-state data                                      #
################################################################################
#                                                                              #
#  Simulates multi-state data                                                  #
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

simfms <- function(nsim  = NULL,
                   tmat  = NULL,
                   clock = "forward",
                   # Frailty
                   frailty = list(dist="gamma",
                                  par= .5),
                   nclus = NULL, 
                   csize = NULL,
                   # Covariates
                   covs = NULL,
                   beta = NULL,
                   # Marginals
                   marg  = list(dist="weibull",
                                lambda=1, rho=1), 
                   cens  = list(dist="weibull", 
                                par = c(lambda=1, rho=1), 
                                admin= 72),
                   # Copula
                   copula= list(name="clayton",
                                par= 1)
                   ) {
  ### - CONTROLS - #############################################################
  checks <- checks(nsim=nsim, tmat=tmat, clock=clock,
                   frailty=frailty, nclus=nclus,  csize=csize,
                   covs=covs, beta=beta, marg=marg, 
                   cens=cens, copula=copula)
  nsim  <- checks$nsim
  nclus <- checks$nclus
  csize <- checks$csize
  clock <- checks$clock
  marg  <- checks$marg
  rm("checks")
  data <- data.frame(ID=1:nsim)
  ###################################################### - END of CONTROLS - ###

  
  ### - FRAILTIES - ############################################################
  data <- cbind(data, simFrail(Fdist=frailty$dist, 
                               Fpar =frailty$par,
                               nsim=nsim, nclus=nclus, csize=csize))
  ##################################################### - END of FRAILTIES - ###
  
  
  ### - COVARIATES - ###########################################################
  if (!is.null(covs))
    data <- cbind(data, simCov(covs=covs, nsim=nsim))
  #################################################### - END of COVARIATES - ###
  
  if (is.null(beta))
    eta <- 0
  else
    eta <- as.matrix(data[,-(1:3)]) %*% t(as.data.frame(beta))
  
  eta <- eta +
         matrix(rep(log(data$z), max(tmat, na.rm=TRUE)), nrow(data), 
                dimnames=list(ID=data$ID, trans=1:max(tmat, na.rm=TRUE)))

  
#   # Starting state
#   startState <- which(colSums(tmat, na.rm=TRUE) == 0)
#   if (length(startState) > 1)
#     stop(paste("This method is implmented for multi-state structures",
#                "with only one starting state!"))
#   which(!is.na(tmat[startState, ]))
#   
#   cat(c(nsim, nclus, csize))
  
  return(eta)
  #return(data)
}
