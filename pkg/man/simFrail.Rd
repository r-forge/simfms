\name{simFrail}
\alias{simFrail}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
simFrail(dist = "gamma", theta = 0.5, nsim = NULL, nclus = NULL, csize = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dist}{
%%     ~~Describe \code{dist} here~~
}
  \item{theta}{
%%     ~~Describe \code{theta} here~~
}
  \item{nsim}{
%%     ~~Describe \code{nsim} here~~
}
  \item{nclus}{
%%     ~~Describe \code{nclus} here~~
}
  \item{csize}{
%%     ~~Describe \code{csize} here~~
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
function (dist = "gamma", theta = 0.5, nsim = NULL, nclus = NULL, 
    csize = NULL) 
{
    if (is.null(nsim) + is.null(nclus) + is.null(csize) == 0) {
        if (!(length(csize) \%in\% c(1, nclus))) 
            stop("Number of clusters 'nclus' and cluster size(s) 'csize' incoherent!")
        if (length(csize) == 1) 
            csize <- rep(csize, nclus)
        nsim2 <- sum(csize)
        if (nsim != nsim2) 
            stop("Number of subjects 'nsim' incoherent to the number of clusters 'nclus' and cluster size(s) 'csize'!")
    }
    if (is.null(nsim) + is.null(nclus) + is.null(csize) == 1) {
        if (is.null(nsim)) {
            if (length(csize) == 1) 
                csize <- rep(csize, nclus)
            else if (length(csize) != nclus) 
                stop("Number of clusters 'nclus' and cluster size(s) 'csize' incoherent!")
            nsim <- sum(csize)
        }
        else if (is.null(nclus)) {
            if (nsim/sum(csize) - nsim\%/\%sum(csize) != 0) 
                stop("Number of subjects 'nsim' and cluster size(s) 'csize' incoherent!")
            nclus <- length(csize)
        }
    }
    if (is.null(nsim) + is.null(nclus) + is.null(csize) >= 2) {
        if (is.null(nsim)) 
            stop("At least the number of simulations 'nsim' must be provided!")
        nclus <- round((runif(1) * nsim)^0.5)
        warning(paste("\n Randomly generated number of clusters 'nclus':", 
            nclus, "\n"))
    }
    if (is.null(csize)) {
        csize <- runif(nclus)
        csize <- round(csize/sum(csize) * nsim)
        csize[nclus] <- nsim - sum(csize[-nclus])
        warning("\n Randomly generated clusters' sizes 'cszie'")
    }
    if (nsim != sum(csize) || length(csize) != nclus) 
        stop("\n --- something went wrong :( ---\n")
    res <- NULL
    attributes(res)$nsim = nsim
    attributes(res)$nclus = nclus
    attributes(res)$csize = csize
    if (substr(dist, 1, 3) == "gam") {
        z <- rgamma(nclus, shape = 1/theta, scale = theta)
        z <- as.vector(unlist(apply(cbind(z, csize), 1, function(x) rep(x[1], 
            x[2]))))
        res$z <- z
        res$Cluster <- as.factor(unlist(apply(cbind(1:nclus, 
            csize), 1, function(x) rep(x[1], x[2]))))
    }
    else (stop("unknown frailty distribution!\n"))
    return(as.data.frame(res))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
