simFrail <-
function(dist="gamma", 
                     theta=.5, 
                     nsim=NULL, 
                     nclus=NULL, 
                     csize=NULL) {
  if (is.null(nsim)+is.null(nclus)+is.null(csize)==0) {
    if (!(length(csize)%in%c(1,nclus))) stop(
      "Number of clusters 'nclus' and cluster size(s) 'csize' incoherent!")
    if (length(csize)==1) csize <- rep(csize, nclus)
    nsim2 <- sum(csize)
    if (nsim!=nsim2) stop("Number of subjects 'nsim' incoherent to the number of clusters 'nclus' and cluster size(s) 'csize'!")
  }
  if (is.null(nsim)+is.null(nclus)+is.null(csize)==1) {
    if(is.null(nsim)){
      if (length(csize)==1) csize <- rep(csize, nclus) else
        if (length(csize)!=nclus) stop("Number of clusters 'nclus' and cluster size(s) 'csize' incoherent!")
      nsim <- sum(csize)
    } else if (is.null(nclus)){
      if (nsim/sum(csize)-nsim%/%sum(csize)!=0) stop("Number of subjects 'nsim' and cluster size(s) 'csize' incoherent!")
      nclus <- length(csize)
    }
  }
  if (is.null(nsim)+is.null(nclus)+is.null(csize)>=2) {
    if (is.null(nsim)) stop("At least the number of simulations 'nsim' must be provided!")
    nclus <- round((runif(1)*nsim)^.5)
    warning(paste("\n Randomly generated number of clusters 'nclus':", nclus, "\n"))
  }
  if (is.null(csize)){
    csize <- runif(nclus)
    csize <- round(csize/sum(csize) * nsim)
    csize[nclus] <- nsim - sum(csize[-nclus])
    warning("\n Randomly generated clusters' sizes 'cszie'")
  }
  if(nsim!=sum(csize)||length(csize)!=nclus) stop("\n --- something went wrong :( ---\n")
  
  res <- NULL
  attributes(res)$nsim=nsim
  attributes(res)$nclus=nclus
  attributes(res)$csize=csize
  
  if (substr(dist,1,3)=="gam"){
    z <- rgamma(nclus,  shape=1/theta, scale=theta)
    z <- as.vector(unlist(apply(cbind(z,csize), 1, function(x) rep(x[1],x[2]))))
    res$z <- z
    res$Cluster <- as.factor(unlist(apply(cbind(1:nclus,csize), 1, function(x) rep(x[1],x[2]))))
  } else(stop("unknown frailty distribution!\n"))
  
  return(as.data.frame(res))
}
