\name{simfms}
\alias{simfms}

\title{Simulation of clustered multi-state data}
\description{
Random generation of survival data with a clustered and multi-state structure.

The simulation procedure is based on a copula model for each 
competing events block, allowing to specify the marginal distributions
of time variables.

The effect of simulated frailties and covariates can be added 
in a proportional hazards way.
%
%The tuning of parameters is done by numerical minimization of a criterion function based on the ratios of target and observed values of median times and of probabilities of competing events.
}
\usage{
simfms(nsim = NULL, tmat = NULL, clock = "forward", 
       frailty = list(dist="gamma", par= .5, type="shared"), 
       nclus = NULL, csize = NULL, 
       covs = NULL, beta = NULL, 
       marg = list(dist = "weibull", lambda = 1, rho = 1), 
       cens = list(dist = "weibull", lambda = 1, rho = 1, admin = 72), 
       copula = list(name = "clayton", par = 1),
       format = "long")
}

\arguments{
  \item{nsim}{The number of subjects to simulate}
  \item{tmat}{The transitions matrix}
  \item{clock}{The approach to be used for tractation of baseline hazards and generation of times: 
    \tabular{ll}{
    \code{reset} \tab
      generates times at each step ranging from 0 to infinity 
      and adds them to those of previous transitions\cr
    \code{forward} \tab
      generates times at each step, conditional on the time of entering the present state
      }
    }
  \item{frailty}{The frailty term specifications. A list with components:
    \itemize{
      \item \code{dist}:
        the name of the frailty distribution 
        (only \code{gamma} implemented at now)
      \item \code{par}:
        the frailty parameter(s) value
      \item \code{type}: the type of frailty
        \itemize{
          \item \code{shared}:
            iid and
            shared by the times of transitions of any type,
            for all subjects in a group
          \item \code{iid}:
            iid and
            shared by the times of transitions of the same type,
            for all subjects in a group
            (possible different frailty parameters 
            for different transition types)
          \item \code{nested}:
            id and correlated by means of nested frailties            
            (possible different frailty parameters 
            for different transition types)        
        }
    }
  }
  \item{nclus}{The number of clusters to simulate}
  \item{csize}{The size(s) of cluster}
  \item{covs}{The covariates to simulate. A list with components 
    \code{nameOfCovariate = simulationfunction}}
  \item{beta}{The regression coefficients for covariates. 
    A list of the same length as \code{covs}, with elements
    \code{nameOfCovariate}, 
    vectors of the same length as the number of transitions 
    \code{max(tmat, na.rm=TRUE)}}
  \item{marg}{The marginal baseline hazards. A list with components
    \itemize{
      \item \code{dist}:
        the name of the baseline hazard distributions
        (one value for all transitions)
      \item \code{eachpar}:
      each baseline parameter (either one value
      or as many as the number of transitions in code{tmat})
      }
    }
  \item{cens}{The censoring time distributions. A list with components
    \itemize{
      \item \code{dist}:
        the name of the censoring distributions
        (one value for all transitions)
      \item \code{eachpar}:
        each censoring distribution parameter (either one value
        or as many as the number of possible starting states
        in \code{tmat}, possibly with states' names)
       }
     }
  \item{copula}{The copula model. A list with components
    \itemize{
      \item \code{name}: the name of the copula 
        (only \code{clayton} available at now)
      \item \code{par}: the copula parameter
    }
  }
  \item{format}{the data format, either \code{long} or \code{wide}
  }
}

%%\details{}

\value{A \code{\link{data.frame}} containing the simulated data, with columns}

\references{Rotolo, F, Legrand, Catherine, Van Keilegom, I (2012)
  Simulation of clustered multi-state survival data based on a copula model.
  Submitted}
  
\author{Federico Rotolo <federico.rotolo@stat.unipd.it>}

%%\note{}

\seealso{\code{\link{trans.cancer}}}

\examples{
################################################################################
### - Cancer reduced multi-state structure - ###################################
################################################################################
trans.cancer.reduced()
simfms(nsim  = NULL,
       tmat  = trans.cancer.reduced(),
       clock = "forward",
       frailty = list(dist="gamma", par= .5),
       nclus = 5,
       csize = 2,
       covs = list(age=function(x) rnorm(x, mean=60, sd=7),
                   treat=function(x) rbinom(x, 1, .5)),
       beta = list(age=rep(.02,5), treat=rep(2,5)),
       marg  = list(dist="weibull", lambda=1, rho=1), 
       cens  = list(dist="weibull", lambda=1, rho=1, admin= 72),
       copula= list(name="clayton", par= 1))


################################################################################
### - Cancer standard multi-state structure - ##################################
################################################################################
trans.cancer()
simfms(nsim  = NULL,
       tmat  = trans.cancer(),
       clock = "forward",
       frailty = list(dist="gamma", par= .5),
       nclus = 5,
       csize = 2,
       covs = list(age=function(x) rnorm(x, mean=60, sd=7),
                   treat=function(x) rbinom(x, 1, .5)),
       beta = list(age=rep(.02, 8), treat=rep(2, 8)),
       marg  = list(dist="weibull", lambda=1, rho=1), 
       cens  = list(dist="weibull", lambda=1, rho=1,  admin= 72),
       copula= list(name="clayton", par= 1))


################################################################################
### - Cancer extended multi-state structure - ##################################
################################################################################
trans.cancer.extended()
simfms(nsim  = NULL,
       tmat  = trans.cancer.extended(),
       clock = "forward",
       frailty = list(dist="gamma", par= .5),
       nclus = 5,
       csize = 2,
       covs = list(age=function(x) rnorm(x, mean=60, sd=7),
                   treat=function(x) rbinom(x, 1, .5)),
       beta = list(age=rep(.02, 9), treat=rep(2, 9)),
       marg  = list(dist="weibull", lambda=1, rho=1), 
       cens  = list(dist="weibull", lambda=1, rho=1, admin= 72),
       copula= list(name="clayton", par= 1))
}

\keyword{frailty}
\keyword{multi-state}
\keyword{multistate}
