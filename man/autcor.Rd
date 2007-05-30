\name{autcor}
\alias{autcor}
\title{Autocorrelation}
\description{
  Estimate autocovariances and autocorrelations.
}
\usage{
autcor(y, lag=NULL, plot=TRUE, lag_axis=TRUE)
}
\arguments{
  \item{y}{a univariate time series.}
  \item{lag}{maximum lag. Default is 2*sqrt(n), where n is the length of the time series y.}
  \item{plot}{logical. If TRUE (default) autocorrelations are plotted.}
  \item{lag_axis}{logical. If TRUE (default) with plot=TRUE, x\_axis is drawn.}
}
\value{
  \item{acov}{autocovariances.}
  \item{acor}{autocorrelations (normalized covariances).}
  \item{mean}{mean of y.}
}

\references{
  H.Akaike and T.Nakagawa (1988) \emph{Statistical Analysis and Control
    of Dynamic Systems.} Kluwer Academic publishers.
}

\examples{
  # Example 1 for the normal distribution 
  y <- rnorm(200)
  autcor(y, lag_axis=FALSE)

  # Example 2 for the ARIMA model
  y <- arima.sim(list(order=c(2,0,0), ar=c(0.64,-0.8)), n=200)
  autcor(y, lag=20)
}

\keyword{ts}
