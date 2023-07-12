\name{thirmo}
\alias{thirmo}
\title{Third Order Moments}
\description{
  Compute the third order moments.
}
\usage{
thirmo(y, lag = NULL, plot = TRUE)
}
\arguments{
  \item{y}{a univariate time series.}
  \item{lag}{maximum lag. Default is \eqn{2 \sqrt{n}}{2*sqrt(n)}, where \eqn{n}
    is the length of the time series \code{y}.}
  \item{plot}{logical. If \code{TRUE} (default), autocovariance \code{acor} is
    plotted.}
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
z <- thirmo(bispecData, lag = 30)
z$tmomnt
}

\keyword{ts}
