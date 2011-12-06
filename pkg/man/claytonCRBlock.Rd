\name{claytonCRBlock}
\alias{claytonCRBlock}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
claytonCRBlock(theta = 1, prev = NULL, nsim = NULL, marg = NULL, pars = NULL, cens = NULL, cpars = NULL, adcens = NULL, eta = NULL, names = NULL, clock = NULL)
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
  \item{cens}{
%%     ~~Describe \code{cens} here~~
}
  \item{cpars}{
%%     ~~Describe \code{cpars} here~~
}
  \item{adcens}{
%%     ~~Describe \code{adcens} here~~
}
  \item{eta}{
%%     ~~Describe \code{eta} here~~
}
  \item{names}{
%%     ~~Describe \code{names} here~~
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
    cens = NULL, cpars = NULL, adcens = NULL, eta = NULL, names = NULL, 
    clock = NULL) 
{
    if (is.null(clock)) 
        clock <- "r"
    if (!(substring(clock, 1, 1) \%in\% c("r", "f"))) 
        stop("Clock parameter 'clock' must be either 'f' (forward) or 'r' (reset)!")
    if (substring(clock, 1, 1) == "f" && !is.null(prev) && !is.null(adcens)) 
        adcens <- adcens - prev
    if (is.null(prev)) {
        if (is.null(nsim)) 
            stop(paste("\nProvide either the conditioning times 'prev'", 
                "or the number of subject to simulate 'nsim'!\n"))
    }
    else {
        if (!is.null(nsim)) 
            stop(paste("\nProvide only one betwen the conditioning times 'prev'", 
                "and the number of subject to simulate 'nsim'!\n"))
        prev <- as.matrix(prev)
        if (ncol(prev) != 1) 
            stop(paste("\nThe conditioning times 'prev' must be either a vector", 
                "or a one-column matrix!"))
    }
    if (!is.null(adcens)) 
        if (min(adcens) <= 0) 
            stop("Administrative censoring time 'adcens' must be positive!")
    if (!is.matrix(pars)) {
        if (is.vector(pars)) 
            pars <- matrix(pars)
        else stop(paste("The parameters object 'pars' must be a matrix", 
            "with one column for each time variable", "and a row for each distribution parameter!"))
    }
    nt <- ncol(pars)
    rest <- matrix(NA, max(nrow(prev), nsim), nt)
    k <- 0
    if (is.null(names)) 
        names <- k + (1:nt)
    colnames(rest) = paste("T", names, sep = "")
    if (!is.null(adcens) || !is.null(cens)) {
        rest <- cbind(rest, C = NA)
        if (!is.null(cens)) {
            cens <- eval(parse(text = paste(cens, "(pars=c(", 
                paste(cpars, collapse = ","), "))", sep = "")))
            uC <- runif(max(nrow(prev), nsim))
            if (substring(clock, 1, 1) == "f" && !is.null(prev)) 
                uC <- uC * cens(prev, inv = FALSE)
            rest[, "C"] <- cens(uC, inv = TRUE)
        }
        if (substring(clock, 1, 1) == "r" && !is.null(prev)) 
            rest[, "C"] <- prev + rest[, "C"]
        if (!is.null(adcens)) 
            rest[, "C"] <- apply(cbind(rest[, "C"], adcens), 
                1, min, na.rm = TRUE)
    }
    rest[, 1] <- clayton(theta = theta, nsim = nsim, prev = prev, 
        marg = marg, pars = pars[, 1], eta = eta[, 1], clock = clock)
    if (nt > 1) {
        for (i in 2:nt) {
            rest[, i] <- clayton(theta = theta, prev = rest[, 
                i - 1], marg = marg, pars = pars[, i], eta = eta[, 
                i], clock = clock)
        }
    }
    return(rest)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
