compTimes <-
function(widedata) {
  nCRs <- length(colnames(widedata)) - 1
  tnames <- colnames(widedata)[1:nCRs]
  
  res <- t(apply(widedata, 1, function(x) {
    c(rep(min(x), nCRs),
      as.numeric(x == min(x))[1:nCRs])      
  }))
  
  colnames(res) <- c(tnames,
                     substr(colnames(widedata), 2, 100)[1:nCRs] )
  
  return(res)
}
