\name{compTimes}
\alias{compTimes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
compTimes(widedata)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{widedata}{
%%     ~~Describe \code{widedata} here~~
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
function (widedata) 
{
    nCRs <- length(colnames(widedata)) - 1
    tnames <- colnames(widedata)[1:nCRs]
    res <- t(apply(widedata, 1, function(x) {
        c(rep(min(x), nCRs), as.numeric(x == min(x))[1:nCRs])
    }))
    colnames(res) <- c(tnames, substr(colnames(widedata), 2, 
        100)[1:nCRs])
    return(res)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
