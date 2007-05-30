\name{mulnos}
\alias{mulnos}
\title{Multiple Unoise}
\description{
  Compute relative power contributions in differential and integrated form, assuming the orthogonality between noise sources.
}
\usage{
  mulnos(y, max.order=NULL, ncon=NULL, nman=0, h, inw=NULL)
}
\arguments{
  \item{y}{a multivariate time series.}
  \item{max.order}{upper limit of model order. Default is 2*sqrt(n), where n is the length of time series y.}
  \item{ncon}{number of controlled variables. Default is d, where d is the dimension of the time series y.}
  \item{nman}{number of maninpulated variables.}
  \item{h}{specify frequencies i/2h (i=0,...,h).}
  \item{inw}{indicator; inw[i] (i=1,...,ncon) indicate the controlled variables and

    inw[i+ncon] (i=1,...,nman) indicate the manipulate variables.}
}
\value{
  \item{nperr}{a normalized prediction error covaiance matrix.}
  \item{diffr}{differential relative power contribution.}
  \item{integr}{integrated relative power contribution.}
}

\references{
  H.Akaike and T.Nakagawa (1988) \emph{Statistical Analysis and Control
    of Dynamic Systems.} Kluwer Academic publishers.
}

\examples{
  ar <- array(0,dim=c(3,3,2))
  ar[,,1] <- matrix(c(0.4,  0,   0.3,
                      0.2, -0.1, -0.5,
                      0.3,  0.1, 0),3,3,byrow=TRUE)
  ar[,,2] <- matrix(c(0,  -0.3,  0.5,
                      0.7, -0.4,  1,
                      0,   -0.5,  0.3),3,3,byrow=TRUE)
  x <- matrix(rnorm(200*3),200,3)
  y <- mfilter(x,ar,"recursive")
  mulnos(y, max.order=10, ncon=3, nman=0, h=20)
}

\keyword{ts}
