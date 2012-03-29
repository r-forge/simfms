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
#                 type: the type of frailty: 'shared', 'iid' or 'nested'       #
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
#   Last modification on: March, 29, 2012                                      #
################################################################################

simfms <- function(nsim  = NULL,
                   tmat  = NULL,
                   clock = "forward",
                   # Frailty
                   frailty = list(dist="gamma",
                                  par= .5,
                                  type="shared"),
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
                                par= 1),
                   format="long"
                   ) {
  exectime <- system.time({
    cat("#####################################################################\n")
    cat("######### Simulation of clustered multi-state survival data #########\n")
    cat("#####################################################################\n")
    
    ### - CONTROLS - #############################################################
    cat("\n# 1 # Checks on input parameters") 
    checked <- checks(nsim=nsim, tmat=tmat, clock=clock,
                      frailty=frailty, nclus=nclus,  csize=csize,
                      covs=covs, beta=beta, marg=marg, 
                      cens=cens, copula=copula, format=format)
    nsim  <- checked$nsim
    nclus <- checked$nclus
    csize <- checked$csize
    frailty$type <- checked$Ftype
    clock <- checked$clock
    marg  <- apply(as.data.frame(checked$marg), 1, extractMargs)
    names(marg) <- 1:length(marg)
    cens  <- list(f = apply(subset(as.data.frame(checked$cens), 
                                   select=-admin), 1, extractMargs),
                  admin = checked$cens$admin)
    rm("checked")
    ###################################################### - END of CONTROLS - ###
    
    ### - INITIALIZATION - #######################################################
    cat("\n\n# 2 # Initialisation of the dataframe object\n")
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
    cat("\n# 3 # Simulation of the time data\n")
    # Detailed data for each transition
    data <- scan.tmat(data=data, inTrans=NULL, #subjs=1:nrow(data),
                      eta=eta,   tmat=tmat,    clock=clock,
                      marg=marg, cens=cens,    copula=copula)
    ########################################### - END of COMPUTATION of TIMES - ###
    
    
    ### - DATA in WIDE FORMAT - ##################################################
    if (substr(format, 1, 1) == "l") {
      cat("\n# 4 # Tranformation of data from wide to long format\n")
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
      ################################################### - END of WIDE FORMAT - ###
    } else {
      resdata <- data
    }
  })
  cat("\n#####################################################################")
  cat(paste("\n# Total elapsed time:", round(exectime[1],2), "sec(s)"))
  cat("\n#####################################################################")
  cat("\n\n\n")
  return(resdata)
}
