\name{mulfrf}
\alias{mulfrf}
\title{Frequency Response Function (Multiple Channel)}
\description{
  Compute multiple frequency response function, gain, phase, multiple coherency, partial coherency and relative error statistics.
}
\usage{
 mulfrf(y, lag=NULL, niv, iovar=c(1:(niv+1)))
}
\arguments{
  \item{y}{a multivariate time series.}
  \item{lag}{maximum lag. Default is 2*sqrt(n), where n is the length of the time series y.}
  \item{niv}{number of input variables.}
  \item{iovar}{input and output variables. iovar(i)-th variables (1<=i<=niv) are input variables and iovar(niv+1)-th variable is output variable.}
}
\value{
  \item{cospec}{spectrum (complex).}
  \item{freqr}{frequency response function : real part.}
  \item{freqi}{frequency response function : imaginary part.}
  \item{gain}{gain.}
  \item{phase}{phase.}
  \item{pcoh}{partial coherency.}
  \item{errstat}{relative error statistics.}
  \item{mcoh}{multiple coherency.}
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
  mulfrf(y, lag=20, niv=2)
}

\keyword{ts}
