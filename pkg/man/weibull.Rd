\name{weibull}
\alias{weibull}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
weibull(pars = list(lambda = 1, rho = 1))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pars}{
%%     ~~Describe \code{pars} here~~
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
function (pars = list(lambda = 1, rho = 1)) 
{
    if (length(pars) != 2) 
        stop("Weibull distribution needs 2 parameters!")
    resf <- Vectorize(function(x, inv = FALSE) {
        if (!is.logical(inv)) 
            stop("Parameter inv must be logical.\n")
        if (inv) {
            if (x < 0 || x > 1) 
                stop("With inv=TRUE, argument 'x' must be between 0 and 1!\n")
            c(T = (-log(x)/pars[[1]])^(1/pars[[2]]))
        }
        else {
            if (x < 0) 
                stop("With inv=FALSE, argument 'x' must be greater than 0!\n")
            c(S = exp(-pars[[1]] * x^pars[[2]]))
        }
    })
    attributes(resf)$pars <- pars
    return(resf)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
