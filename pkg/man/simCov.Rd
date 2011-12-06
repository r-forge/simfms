\name{simCov}
\alias{simCov}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
simCov(covs, nsim = 3000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{covs}{
%%     ~~Describe \code{covs} here~~
}
  \item{nsim}{
%%     ~~Describe \code{nsim} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (covs, nsim = 3000) 
{
    if (is.null(nsim)) 
        stop("Please, set the number of subjects 'nsim'!")
    if (is.null(cov)) 
        stop("Please, set the covariates object 'covs'!")
    if (!is.list(covs)) 
        stop("The covariates object 'covs' must be a list!")
    if (is.null(names(covs)) || "" \%in\% names(covs)) 
        stop("All the names of covariates must be specified as\n         covs=list(name1=..., ..., nameN=...)")
    res <- NULL
    for (i in 1:length(covs)) {
        if (length(covs[[i]]) != 1) 
            stop(paste("The covariate '", names(covs)[i], "' has length ", 
                length(covs[[i]]), " instead of 1!", sep = ""))
        if (!is.function(covs[[i]])) 
            stop(paste("The element '", names(covs)[i], "' of the covariates list 'covs' is not a function!", 
                sep = ""))
        res[[names(covs)[[i]]]] <- covs[[i]](nsim)
    }
    return(as.data.frame(res))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
