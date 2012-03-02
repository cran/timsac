\name{mulnos}
\alias{mulnos}
\title{Relative Power Contribution}
\description{
  Compute relative power contributions in differential and integrated form, assuming the orthogonality between noise sources.
}
\usage{
  mulnos(y, max.order=NULL, control=NULL, manip=NULL, h)
}
\arguments{
  \item{y}{a multivariate time series.}
  \item{max.order}{upper limit of model order. Default is \eqn{2 \sqrt{n}}{2*sqrt(n)}, where \eqn{n} is the length of time series \code{y}.}
  \item{control}{controlled variables. Default is \eqn{c(1:d)}, where \eqn{d} is the dimension of the time series \code{y}.}
  \item{manip}{manipulated variables. Default number of manipulated variable is \eqn{0}.}
  \item{h}{specify frequencies \eqn{i/2}\code{h} (\eqn{i=0,...,}\code{h}).}
}
\value{
  \item{nperr}{a normalized prediction error covariance matrix.}
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
  mulnos(y, max.order=10, h=20)
}

\keyword{ts}
