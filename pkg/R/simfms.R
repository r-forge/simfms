################################################################################
#  Simulation of frailty multi-state data                                      #
################################################################################
#                                                                              #
#  Simulates multi-state data                                                  #
#                                                                              #
#  Its parameters are                                                          #
#   - nsim      : the number of subjects to simulate                           #
#   - tmat      : the transitions matrix                                       #
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
#                 dist    : the name of the baseline hazard distributions      #
#                           (either one value or as many as the number         #
#                            of transitions in 'tmat')                         #
#                 eachpar : each baseline parameter                            #
#                           (either one value or as many as the number         #
#                            of transitions in 'tmat')                         #
#   - cens      : the censoring time distributions. A list with components     #
#                 dist : the name of the censoring distributions               #
#                           (either one value or as many as the number of      #
#                            possible starting states in 'tmat',               #
#                            possibly with states' names)                      #
#                 eachpar : each censoring distribution parameter              #
#                           (either one value or as many as the number of      #
#                            possible starting states in 'tmat',               #
#                            possibly with states' names)                      #
#                 admin: the time of administrative censoring                  #
#   - copula    : the copula model. A list with components                     #
#                 name : the name of the copula                                #
#                 par  : the copula parameter                                  #
#                                                                              #
#                                                                              #
#   Date: February, 13, 2012                                                   #
#   Last modification on: February, 16, 2012                                   #
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
                                lambda=1, rho=1, 
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
  cens  <- checks$cens
  rm("checks")
  ###################################################### - END of CONTROLS - ###

  
  ### - DATA OBJECT INITIALIZATION #############################################
  data <- as.data.frame(cbind(1:nsim, matrix(c(NA, 0), 
                                             nsim, 2 * max(tmat, na.rm=TRUE),
                                             byrow=TRUE)))
  names(data) <- c("ID", sapply(1:max(tmat, na.rm=TRUE), function(x)
    paste("tr", x, c(".time", ".status"), sep="")))
  
  #################################### - END of DATA OBJECT INITIALIZATION - ###

  
  ### - FRAILTIES - ############################################################
  data <- cbind(data, simFrail(Fdist=frailty$dist, 
                               Fpar =frailty$par,
                               nsim=nsim, nclus=nclus, csize=csize))
  ##################################################### - END of FRAILTIES - ###
  
  
  ### - COVARIATES - ###########################################################
  if (!is.null(covs))
    data <- cbind(data, simCov(covs=covs, nsim=nsim))
  #################################################### - END of COVARIATES - ###
  
  ### - LINEAR PREDICTORS - ####################################################
  # Covariates constributions
  notCov <- 1:(3 +                      # ID, z, Cluster
    2 * max(tmat, na.rm=TRUE))  # tr1.time, tr1.status, ...
  if (is.null(beta))
    eta <- 0  else
      eta <- as.matrix(data[,-notCov]) %*% t(as.data.frame(beta))
  
  # Frailty contribution
  eta <- eta +
         matrix(rep(log(data$z), max(tmat, na.rm=TRUE)), nrow(data), 
                dimnames=list(ID=data$ID, trans=1:max(tmat, na.rm=TRUE)))
  ############################################# - END of LINEAR PREDICTORS - ###

  ### - COMPUTATION of TRANSITION TIMES - ######################################
  # Detailed data for each transition
  data <- scan.tmat(data=data, inTrans=NULL, subjs=1:nrow(data),
                    eta=eta,   tmat=tmat,    clock=clock,
                    marg=marg, cens=cens,    copula=copula)
  
  # Possible arrival states
  events <- colnames(tmat)[colSums(tmat, na.rm=TRUE) > 0]
  
  # Summarization into possible arrival states
  resdata <- as.data.frame(lapply(events, function(x) {
    res <- cbind(apply(data[, paste("tr", sort(tmat[, x]), ".time", sep=""),
                            drop=FALSE], 
                       1, function(x) {
                         if (all(is.na(x))) NA else
                         max(x, na.rm=TRUE)}),
                 apply(data[, paste("tr", sort(tmat[, x]), ".status", sep=""),
                            drop=FALSE], 
                       1, max, na.rm=TRUE))
    colnames(res) <- paste(x, c("time", "status"), sep=".")
    return(res)
  }))
  ########################################## - END of COMPUTATION of TIMES - ###
  
  return(resdata)
}
