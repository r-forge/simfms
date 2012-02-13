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
#   Last modification on: February, 13, 2012                                   #
################################################################################

simFrail <-function(Fdist="gamma", 
                    Fpar=.5, 
                    nsim=NULL, 
                    nclus=NULL, 
                    csize=NULL) {
  res <- NULL
  attributes(res)$nsim  <- nsim
  attributes(res)$nclus <- nclus
  attributes(res)$csize <- csize
  
  if (substr(Fdist, 1, 3)=="gam") {
    z <- rgamma(nclus,  shape=1/Fpar, scale=Fpar)
    z <- as.vector(unlist(apply(cbind(z,csize), 1, function(x) rep(x[1],x[2]))))
    res$z <- z
    res$Cluster <- as.factor(unlist(apply(cbind(1:nclus,csize), 1, 
                                          function(x) rep(x[1],x[2]))))
  } 
  else
    stop(paste("Unknown frailty distribution '", Fdist, "'!\n", sep=""))
  
  return(as.data.frame(res))
}
