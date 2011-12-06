loglogistic <-
function(pars=list(alpha=1, kappa=1)){
  resf <- Vectorize(function(x, inv=FALSE){
    if (!is.logical(inv)) stop("Parameter inv must be logical.\n")
      if (inv){
        if (x<0 || x>1) stop("With inv=TRUE, argument 'x' must be between 0 and 1!\n")
        c(T = (-(1-1/x)*exp(-pars[[1]]))^(1/pars[[2]]) )
      } else {
        if (x<0) stop("With inv=FALSE, argument 'x' must be greater than 0!\n")
        c( S = 1/(1+exp(pars[[1]])*x^pars[[2]]) )
      }
    })
  attributes(resf)$pars <- pars
  return(resf)
}
