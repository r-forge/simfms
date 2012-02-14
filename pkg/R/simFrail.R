################################################################################
#  Simulation of frailty term                                                  #
################################################################################
#                                                                              #
#  Simulates the frailty term for a frailty multi-state simulation model       #
#                                                                              #
#  Its parameters are                                                          #
#   - Fdist     : the name of the frailty distribution                         #
#   - Fpar      : the frailty parameter                                        #
#   - nsim      : the number of subjects to simulate                           #
#   - nclus     : the number of clusters to simulate                           #
#   - csize     : the size(s) of cluster                                       #
#                                                                              #
#                                                                              #
#   Date: February, 13, 2012                                                   #
#   Last modification on: February, 14, 2012                                   #
################################################################################

simFrail <-function(Fdist=NULL, 
                    Fpar=NULL, 
                    nsim=NULL, 
                    nclus=NULL, 
                    csize=NULL) {
  
  if (substr(Fdist, 1, 3)=="gam") 
    z <- rgamma(nclus,  shape=1/Fpar, scale=Fpar)
  else if (substr(Fdist, 1, 2)=="no") {
    nclus <- 1
    csize <- nsim
    z <- rep(1, nclus)
  }
  else
    stop(paste("Unknown frailty distribution '", Fdist, "'!\n", sep=""))
  
  res <- NULL
  res$z <- as.vector(unlist(mapply(rep, x=z, times=csize)))
  res$Cluster <- as.factor(unlist(mapply(rep, x=1:nclus, times=csize)))
  attributes(res)$nsim  <- nsim
  attributes(res)$nclus <- nclus
  attributes(res)$csize <- csize

  return(as.data.frame(res))
}
