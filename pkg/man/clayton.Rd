\name{clayton}
\alias{clayton}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
clayton(theta = 1, prev = NULL, nsim = NULL, marg = NULL, pars = NULL, eta = NULL, clock = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{
%%     ~~Describe \code{theta} here~~
}
  \item{prev}{
%%     ~~Describe \code{prev} here~~
}
  \item{nsim}{
%%     ~~Describe \code{nsim} here~~
}
  \item{marg}{
%%     ~~Describe \code{marg} here~~
}
  \item{pars}{
%%     ~~Describe \code{pars} here~~
}
  \item{eta}{
%%     ~~Describe \code{eta} here~~
}
  \item{clock}{
%%     ~~Describe \code{clock} here~~
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
function (theta = 1, prev = NULL, nsim = NULL, marg = NULL, pars = NULL, 
    eta = NULL, clock = NULL) 
{
    if (is.null(clock)) 
        clock <- "r"
    if (!(substring(clock, 1, 1) \%in\% c("r", "f"))) 
        stop("Clock parameter 'clock' must be either 'f' (forward) or 'r' (reset)!")
    if (is.null(prev)) {
        if (is.null(nsim)) 
            stop(paste("\nProvide either the conditioning times 'prev'", 
                "or the number of subject to simulate 'nsim'!\n"))
        k <- 1
    }
    else {
        if (!is.null(nsim)) 
            stop(paste("\nProvide only one betwen the conditioning times 'prev'", 
                "and the number of subject to simulate 'nsim'!\n"))
        prev <- as.matrix(prev)
        k <- ncol(prev) + 1
        nsim <- nrow(prev)
    }
    marg <- eval(parse(text = paste(marg, "(pars=c(", paste(pars, 
        collapse = ","), "))", sep = "")))
    if (!(is.null(nsim) || is.null(eta))) 
        if (nsim != length(eta)) 
            stop(paste("The number of simulations 'nsim' does not match", 
                "the number of simulations of the covariates in 'eta'!"))
    if (is.null(eta)) 
        eta <- rep(0, nsim)
    margEta <- margCovFrail(marg, lp = eta)
    u <- runif(nsim)
    if (substring(clock, 1, 1) == "f" && !is.null(prev)) 
        u <- u * margEta(prev, inv = FALSE)
    if (k == 1) {
        tk <- margEta(u, inv = TRUE)
    }
    else {
        tk <- margEta((1 + (u^(theta/(theta * (1 - k) - 1)) - 
            1) * (1 + apply(cbind(prev, eta), 1, function(x) sum(margCovFrail(marg, 
            lp = x[2])(x[1], inv = FALSE)^(-theta) - 1))))^(-1/theta), 
            inv = TRUE)
    }
    if (substring(clock, 1, 1) == "r" && !is.null(prev)) 
        tk <- prev + tk
    tk <- matrix(as.numeric(tk), dimnames = list(NULL, paste("T", 
        k, sep = "")))
    return(tk)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
