simCov <-
function(covs, nsim=3000){
  if (is.null(nsim)) stop("Please, set the number of subjects 'nsim'!")
  if (is.null(cov)) stop("Please, set the covariates object 'covs'!")
  if (!is.list(covs)) stop("The covariates object 'covs' must be a list!")
  if (is.null(names(covs))||""%in%names(covs))
    stop("All the names of covariates must be specified as
         covs=list(name1=..., ..., nameN=...)")
  res <- NULL
  for (i in 1:length(covs)){
    if (length(covs[[i]])!=1) stop(paste("The covariate '", names(covs)[i],
      "' has length ", length(covs[[i]]), " instead of 1!", sep=""))
    if (!is.function(covs[[i]])) stop(paste("The element '", names(covs)[i],
      "' of the covariates list 'covs' is not a function!", sep=""))
    res[[names(covs)[[i]]]] <- covs[[i]](nsim)
  }
  return(as.data.frame(res))
}
