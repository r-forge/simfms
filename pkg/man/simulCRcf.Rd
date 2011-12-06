\name{simulCRcf}
\alias{simulCRcf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
simulCRcf(nsim = NULL, trans = NULL, clock = "reset", fdist = "gamma", ftheta = 0.5, nclus = NULL, csize = NULL, covs, beta, ctheta = 1, marg = "weibull", pars = c(lambda = 1, rho = 1), cens = "weibull", cpars = c(lambda = 1, rho = 1), adcens = 72)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nsim}{
%%     ~~Describe \code{nsim} here~~
}
  \item{trans}{
%%     ~~Describe \code{trans} here~~
}
  \item{clock}{
%%     ~~Describe \code{clock} here~~
}
  \item{fdist}{
%%     ~~Describe \code{fdist} here~~
}
  \item{ftheta}{
%%     ~~Describe \code{ftheta} here~~
}
  \item{nclus}{
%%     ~~Describe \code{nclus} here~~
}
  \item{csize}{
%%     ~~Describe \code{csize} here~~
}
  \item{covs}{
%%     ~~Describe \code{covs} here~~
}
  \item{beta}{
%%     ~~Describe \code{beta} here~~
}
  \item{ctheta}{
%%     ~~Describe \code{ctheta} here~~
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
function (nsim = NULL, trans = NULL, clock = "reset", fdist = "gamma", 
    ftheta = 0.5, nclus = NULL, csize = NULL, covs, beta, ctheta = 1, 
    marg = "weibull", pars = c(lambda = 1, rho = 1), cens = "weibull", 
    cpars = c(lambda = 1, rho = 1), adcens = 72) 
{
    res <- NULL
    if (is.null(nsim)) 
        stop("Number of simulations 'nsim' not defined!")
    if (!(is.numeric(nsim) && length(nsim) == 1)) 
        stop(paste("The number of simulations 'nsim' must be", 
            "a one-dimensional numeric object!"))
    if (!is.matrix(pars)) {
        if (is.vector(pars)) 
            pars <- matrix(pars)
        else stop(paste("The parameters object 'pars' must be a matrix", 
            "with one column for each time variable", "and a row for each distribution parameter!"))
    }
    nt <- ncol(pars)
    if (max(trans, na.rm = TRUE) != nt) 
        stop(paste("The number of transitions in the transition matrix 'trans'", 
            "does not match with the number of columns of 'pars'!"))
    if (length(which(rowSums(trans, na.rm = TRUE) > 0)) != ncol(cpars)) 
        stop(paste("The number of competing events blocks", "(number of rows of 'trans' with non NA elements)", 
            "does not match with the number of columns of 'cpars'!"))
    if (is.null(fdist) && is.null(ftheta)) {
        X.frail <- NULL
        eta.frail <- rep(0, nsim)
    }
    else {
        if (is.null(fdist)) 
            stop(paste("No frailty distribution 'fdist' is specified,", 
                "while its dispersion 'ftheta' is!"))
        if (is.null(ftheta)) 
            stop(paste("No frailty dispersion 'ftheta' is specified,", 
                "while frailty distribution 'fdist' is!"))
        X.frail <- simFrail(dist = fdist, theta = ftheta, nsim = nsim, 
            nclus = nclus, csize = csize)
        eta.frail <- log(X.frail[, "z"])
        attributes(res)$frailty <- list(dist = fdist, frail = X.frail)
    }
    if (is.null(beta) && is.null(covs)) {
        X.cov <- NULL
        eta.cov <- rep(0, nsim)
    }
    else {
        if (is.null(beta)) 
            stop(paste("No regression coefficients 'beta' are specified,", 
                "while covariates 'covs' are!"))
        if (is.null(covs)) 
            stop(paste("No covariates 'covs' are specified,", 
                "while regression coefficients 'beta' are!"))
        if (length(beta) != length(covs)) 
            stop(paste("\nThe number of covariates (", length(covs), 
                ") is different from the number of regression coefficients (", 
                length(beta), ")!\n", sep = ""))
        if (diff(range((sapply(beta, length)))) != 0) 
            stop("All beta coefficients in 'beta' must have the same length!")
        if (length(beta[[1]]) != nt) 
            stop(paste("The length of each beta coefficient in 'beta' must be", 
                "the same as the length of each parameter in 'pars'!"))
        X.cov <- simCov(covs = covs, nsim = nsim)
        eta.cov <- as.matrix(X.cov) \%*\% t(matrix(unlist(beta), 
            nt, dimnames = list(NULL, names(beta))))
        attributes(res)$covariates <- list(dist = covs, beta = beta)
    }
    eta <- eta.cov + eta.frail
    if (!(is.null(X.frail) || is.null(X.cov))) 
        Xdata <- cbind(X.frail, X.cov)
    else if (is.null(X.frail)) 
        Xdata <- X.cov
    else Xdata <- X.frail
    ss <- which(colSums(trans, na.rm = TRUE) == 0)
    if (length(ss) > 1) 
        stop(paste("This method is implmented for multi-state structures", 
            "with only one starting state!"))
    dest.states <- which(!is.na(trans[ss, ]))
    transitions <- trans[ss, dest.states]
    res <- compTimes(claytonCRBlock(nsim = nsim, theta = ctheta, 
        marg = marg, pars = pars[, transitions], cens = cens, 
        cpars = cpars[, ss], adcens = adcens, eta = eta, names = colnames(trans)[dest.states], 
        clock = clock))
    pres.state <- apply(res[, (ncol(res)/2 + 1):ncol(res)], 1, 
        function(x) {
            ifelse(sum(x) != 0, names(x)[which(x == 1)], NA)
        })
    repeat {
        startstates <- which(sapply(rownames(trans), function(x) {
            x \%in\% names(table(pres.state)) && sum(trans[x, ], 
                na.rm = TRUE) > 0
        }))
        if (length(startstates) == 0) 
            break
        for (ss in startstates) {
            whichss <- which(pres.state == rownames(trans)[ss])
            dest.states <- which(!is.na(trans[ss, ]))
            transitions <- trans[ss, dest.states]
            res2 <- compTimes(claytonCRBlock(prev = res[whichss, 
                paste("T", rownames(trans)[ss], sep = "")], theta = ctheta, 
                marg = marg, pars = pars[, transitions], cens = cens, 
                cpars = cpars[, ss], adcens = adcens, eta = eta[whichss, 
                  ], names = colnames(trans)[dest.states], clock = clock))
            for (ds in dest.states) {
                if (!colnames(trans)[ds] \%in\% colnames(res)) 
                  res <- cbind(res, matrix(NA, nrow(res), 2, 
                    dimnames = list(NULL, c(paste("T", colnames(trans)[ds], 
                      sep = ""), colnames(trans)[ds]))))
                res[whichss, paste(c("T", ""), colnames(trans)[dest.states], 
                  sep = "")] <- res2[, paste(c("T", ""), colnames(trans)[dest.states], 
                  sep = "")]
            }
            pres.state[whichss] <- apply(res2[, (ncol(res2)/2 + 
                1):ncol(res2)], 1, function(x) {
                ifelse(sum(x) != 0, names(x)[which(x == 1)], 
                  NA)
            })
        }
    }
    return(res)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
