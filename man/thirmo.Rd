\name{thirmo}
\alias{thirmo}
\title{Third Order Monents}
\description{
  Compute the third order moments.
}
\usage{
thirmo(y, lag=NULL, plot=TRUE)
}
\arguments{
  \item{y}{a univariate time series.}
  \item{lag}{maximum lag. Default is 2*sqrt(n), where n is the length of the time series y.}
  \item{plot}{logical. If TRUE (default) autocovariance acor is plotted.}
}
\value{
  \item{mean}{mean.}
  \item{acov}{autocovariance.}
  \item{acor}{normalized covariance.}
  \item{tmomnt}{third order moments.}
}

\references{
  H.Akaike, E.Arahata and T.Ozaki (1975) \emph{Computer Science Monograph, No.6,
   Timsac74, A Time Series Analysis and Control Program Package (2)}.
  The Institute of Statistical Mathematics.
}

\examples{
  data(bispecData)
  z <- thirmo(bispecData, lag=30)
  z$tmomnt
}

\keyword{ts}
