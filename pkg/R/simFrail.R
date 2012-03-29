################################################################################
#  Simulation of frailty term                                                  #
################################################################################
#                                                                              #
#  Simulates the frailty term for a frailty multi-state simulation model       #
#                                                                              #
#  Its parameters are                                                          #
#   - Fdist     : the name of the frailty distribution                         #
#   - Ftype     : the type of frailty: 'shared', 'iid' or 'nested'             #
#   - Fpar      : the frailty parameter, with dimension                        #
#                 1               if Ftype is 'shared'                         #
#                 1 or ntrans     if Ftype is 'iid'                            #
#                 2 or ntrans + 1 if Ftype is 'nested                          #
#   - nclus     : the number of clusters                                       #
#   - ntrans    : the number of transtion types                                #
#                                                                              #
#                                                                              #
#   Date: February, 13, 2012                                                   #
#   Last modification on: March, 29, 2012                                      #
################################################################################

simFrail <-function(Fdist=NULL,
                    Ftype="shared",
                    Fpar=NULL, 
                    nclus=NULL, 
                    ntrans=NULL) {
  
  Ftype <- substr(Ftype, 1, 1)  
  
  z <- matrix(1, nrow = nclus, ncol = ntrans)
  
  if (Ftype != "i") {
    cat("Simultaion of shared frailty...\n")
    if (substr(Fdist, 1, 3) == "gam") {
      z <- z * matrix(rep(rgamma(nclus,  shape=1/Fpar[1], scale=Fpar[1]), 
                          ntrans), ncol = ntrans)
    } else if (substr(Fdist, 1, 2) != "no") {
      stop(paste("Unknown frailty distribution '", Fdist, "'!\n", sep=""))
    }
  }
  
  if (Ftype != "s") {
    cat("Simultaion of transition-specific frailty...\n")
    if (substr(Fdist, 1, 3) == "gam") {
      Rgamma <- function(nclus, Fpar, ntrans) {
        rgamma(nclus, shape=1/Fpar, scale=Fpar)
      }
      VRgamma <- function(nclus, Fpar, ntrans) {
        if (length(Fpar) == 1) {
          Fpar <- rep(Fpar, ntrans)
        }
        Vectorize(Rgamma, "Fpar")(nclus, Fpar, ntrans)
      }
      
      if (Ftype == "n") {pars <- Fpar[-1]} else if (Ftype == "i") {pars <- Fpar}
      z <- z * VRgamma(nclus,  pars, ntrans)
      rm(pars)
    } else if (substr(Fdist, 1, 2) != "no") {
      stop(paste("Unknown frailty distribution '", Fdist, "'!\n", sep=""))
    }
  }

  z <- as.data.frame(cbind(1:nclus, z))
  colnames(z) <- c("Cluster", paste("frail", 1:ntrans, sep="."))
  attributes(z)$Fdist  <- Fdist
  attributes(z)$Ftype  <- Ftype
  attributes(z)$Fpar   <- Fpar
  attributes(z)$nclus  <- nclus
  attributes(z)$ntrans <- ntrans

  return(z)
}
