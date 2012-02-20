################################################################################
#  Initialization of the data matrix and of linear predictors                  #
################################################################################
#                                                                              #
#  Builds the data matrix, simulates the frailties and the covariates,         #
#  computes the linear predictors 'eta'
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
#                                                                              #
#                                                                              #
#   Date: February, 20, 2012                                                   #
#   Last modification on: February, 20, 2012                                   #
################################################################################

initialize.fms <- function(nsim,
                           tmat,
                           clock,
                           frailty,
                           nclus,
                           csize,
                           covs,
                           beta
                           ) {
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
  
  return(list(data=data, 
              eta=eta))
}