scan.tmat <- function(tmat,
                     incoming,
                     data=NULL){
  presentstate <- which(tmat == incoming, arr.ind=TRUE)[2]
  # simulate for its CRsBlock
  
  children <- which(!is.na(tmat[presentstate, ]))
  
  # if no  children return results
  if (length(children == 0))
    return(data)
  # else
  else
    {}
  ## pass results to its children
  ## merge their results
  ## return their results
}