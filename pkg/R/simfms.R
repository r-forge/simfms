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
#                 dist    : the name of the baseline hazard distribution       #
#                            (one value)                                       #
#                 eachpar : initial values of each baseline parameter          #
#                           (either one value or as many as the number         #
#                            of transitions in 'tmat')                         #
#   - cens      : the censoring time distributions. A list with components     #
#                 dist : the name of the censoring distributions (one value)   #
#                 eachpar : each censoring distribution parameter              #
#                           (either one value or as many as the number of      #
#                            possible starting states in 'tmat')               #
#                 admin: the time of administrative censoring                  #
#   - copula    : the copula model. A list with components                     #
#                 name : the name of the copula                                #
#                 par  : the copula parameter                                  #
#                                                                              #
#                                                                              #
#   Date: February, 13, 2012                                                   #
#   Last modification on: February, 20, 2012                                   #
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
  checked <- checks(nsim=nsim, tmat=tmat, clock=clock,
                   frailty=frailty, nclus=nclus,  csize=csize,
                   covs=covs, beta=beta, marg=marg, 
                   cens=cens, copula=copula)
  nsim  <- checked$nsim
  nclus <- checked$nclus
  csize <- checked$csize
  clock <- checked$clock
  marg  <- checked$marg
  cens  <- checked$cens
  rm("checked")
  ###################################################### - END of CONTROLS - ###
  
  ### - INITIALIZATION - #######################################################
  initial <- initialize.fms(nsim    = nsim,
                            tmat    = tmat,
                            clock   = clock,
                            frailty = frailty,
                            nclus   = nclus,
                            csize   = csize,
                            covs    = covs,
                            beta    = beta)
  data  <- initial$data
  eta   <- initial$eta
  rm("initial")
  ################################################ - END of INITIALIZATION - ###
  
  
  ### - COMPUTATION of TRANSITION TIMES - ######################################
  # Detailed data for each transition
  data <- scan.tmat(data=data, inTrans=NULL, subjs=1:nrow(data),
                    eta=eta,   tmat=tmat,    clock=clock,
                    marg=marg, cens=cens,    copula=copula)
  ########################################### - END of COMPUTATION of TIMES - ###
  
  
  if (fulldata)
    resdata <- data else {
      
  ### - DATA in WIDE FORMAT - ##################################################
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
    }
  ################################################### - END of WIDE FORMAT - ###
  
  return(resdata)
}
