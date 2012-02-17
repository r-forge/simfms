\name{trans.cancer}
\alias{trans.cancer}
\alias{trans.cancer.reduced}
\alias{trans.cancer.extended}

\title{Transition matrices for cancer clinical trials}

\description{
  Transition matrices for standard multi-state structures
  in cancer clinical trials.
}

\usage{
trans.cancer.reduced( names = c("NED", "LR", "DM",         "De"))
trans.cancer(         names = c("NED", "LR", "DM", "LRDM", "De"))
trans.cancer.extended(names = c("NED", "LR", "DM", "LRDM", "De"))
}

\arguments{
  \item{names}{The names of the states of the multi-state structure}
}

%%\details{}
\value{A transition matrix (see examples)}

\references{
Putter H, Fiocco M, Geskus RB (2007). Tutorial in biostatistics: 
Competing risks and multi-state models. Statistics in Medicine 26, 2389--2430.

L. C. de Wreede, M. Fiocco, and H. Putter (2010). 
The mstate package for estimation and prediction in non- and semi-parametric multi-state and competing risks models. Computer Methods and Programs in Biomedicine 99: 261--274.
}

\author{Federico Rotolo <federico.rotolo@stat.unipd.it>}

%%\note{}

\seealso{\code{\link{simfms}},
          \code{\link[mstate]{mstate}}}
          
\examples{
  trans.cancer.reduced()
  trans.cancer()
  trans.cancer.extended()
}

\keyword{transition}
\keyword{multi-state}
\keyword{multistate}
\keyword{acyclic}
\keyword{cancer}
\keyword{survival}
\keyword{clinical trial}
\keyword{stochastic process}
