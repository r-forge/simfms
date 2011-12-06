weibull <-
function(pars=list(lambda=1, rho=1)){
  if (length(pars)!=2) stop("Weibull distribution needs 2 parameters!")
  resf <- Vectorize(function(x, inv=FALSE){
    if (!is.logical(inv)) 
      stop("Parameter inv must be logical.\n")
    if (inv){
      if (x<0 || x>1)
        stop("With inv=TRUE, argument 'x' must be between 0 and 1!\n")
      c(T = (-log(x)/pars[[1]])^(1/pars[[2]]))
    } else {
      if (x<0)
        stop("With inv=FALSE, argument 'x' must be greater than 0!\n")
      c(S = exp(-pars[[1]] * x^pars[[2]]))
    }
  })
  attributes(resf)$pars <- pars
  return(resf)
}
