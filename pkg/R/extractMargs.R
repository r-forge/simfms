extractMargs <- function(x) {
  resf <- eval(parse(text=as.character(x[["dist"]])))(
    pars=eval(parse(text=paste("list(",
                               paste(names(x)[names(x) != "dist"],
                                     x[names(x) !="dist"],
                                     sep="=", collapse=", "),
                               ")"))))
  return(resf)
}