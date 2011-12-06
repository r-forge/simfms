lognormal <-
function(pars=list(mu=1, sigma=1)){
  resf <- Vectorize(function(x, inv=FALSE){
    if (!is.logical(inv)) stop("Parameter inv must be logical.\n")
      if (inv){
        if (x<0 || x>1) stop("With inv=TRUE, argument 'x' must be between 0 and 1!\n")
        c(T = exp( pars[[1]]+pars[[2]]*qnorm(1-x) ) )
      } else {
        if (x<0) stop("With inv=FALSE, argument 'x' must be greater than 0!\n")
        c( S = 1-pnorm((log(x)-pars[[1]])/pars[[2]]) )
      }
    })
  attributes(resf)$pars <- pars
  return(resf)
}
