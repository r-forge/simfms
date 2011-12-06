margCovFrail <-
function(f, lp=0){
  resf <- function(x, inv=FALSE) {
    if (length(x)!=length(lp)) 
      warning("'x' and 'eta' have different lengths!")
    if (inv)
      c(T = as.numeric(f(x^exp(-lp), inv=TRUE)))
    else
      c(S = as.numeric(f(x, inv=FALSE))^exp(lp) )
  }
  attributes(resf)$nsim <- length(lp)
  return(resf)
}
