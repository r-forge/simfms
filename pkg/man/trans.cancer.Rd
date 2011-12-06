\name{trans.cancer}
\alias{trans.cancer}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
trans.cancer(names)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{names}{
%%     ~~Describe \code{names} here~~
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
function (names) 
{
    tmat <- matrix(NA, 5, 5)
    tmat[1, c(2:3, 5)] <- 1:3
    tmat[2, 4:5] <- 4:5
    tmat[3, 4:5] <- 6:7
    tmat[4, 5] <- 8
    if (missing(names)) 
        names <- c("NED", "LR", "DM", "LRDM", "De")
    else {
        if (length(names) != 5) 
            stop("incorrect length of \"names\" argument")
    }
    dimnames(tmat) <- list(from = names, to = names)
    return(tmat)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
