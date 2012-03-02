\name{mulfrf}
\alias{mulfrf}
\title{Frequency Response Function (Multiple Channel)}
\description{
  Compute multiple frequency response function, gain, phase, multiple coherency, partial coherency and relative error statistics.
}
\usage{
 mulfrf(y, lag=NULL, iovar=NULL)
}
\arguments{
  \item{y}{a multivariate time series.}
  \item{lag}{maximum lag. Default is \eqn{2 \sqrt{n}}{2*sqrt(n)}, where \eqn{n} is the number of rows in \code{y}.}
  \item{iovar}{input variables (\code{iovar[i], i=1,k}) and output variable (\code{iovar[k+1]}) (\eqn{1 \le k \le d}{1<=i<=d}),
 where \eqn{d} is the number of columns in \code{y}. Default is \eqn{c(1:d)}.}
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
  mulfrf(y, lag=20)
}

\keyword{ts}
