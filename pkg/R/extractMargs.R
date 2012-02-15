################################################################################
#  Building of a marginal survival function from its name and parameters       #
################################################################################
#                                                                              #
#  Creates the marginal survival function for given values of the name         #
#  and the names and values of parameters                                      #
#                                                                              #
#  Its parameters are                                                          #
#   - x  : a named vector with element 'dist' the names of the distribution    #
#          and the others which are values of parameters and whose namese are  #
#           the parameters' names                                              #
#                                                                              #
#                                                                              #
#   Date: February, 15, 2012                                                   #
#   Last modification on: February, 15, 2012                                   #
################################################################################

extractMargs <- function(x) {
  resf <- eval(parse(text=as.character(x[["dist"]])))(
    pars=eval(parse(text=paste("list(",
                               paste(names(x)[names(x) != "dist"],
                                     x[names(x) !="dist"],
                                     sep="=", collapse=", "),
                               ")"))))
  return(resf)
}