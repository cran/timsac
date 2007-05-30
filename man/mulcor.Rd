\name{mulcor}
\alias{mulcor}
\alias{print.mulcor}
\title{Multiple Correlation}
\description{
  Estimate multiple correlation.
}
\usage{
  mulcor(y, lag=NULL, plot=TRUE, lag_axis=TRUE)
}
\arguments{
  \item{y}{a multivariate time series.}
  \item{lag}{maximum lag. Default is 2*sqrt(n), where n is the length of the time series y.}
  \item{plot}{logical. If TRUE (default) correlations cor are plotted.}
  \item{lag_axis}{logical. If TRUE (default) with plot=TRUE, x\_axis is drawn.}
}
\value{
  \item{cov}{covariances.}
  \item{cor}{correlations (normalized covariances).}
  \item{mean}{mean.}
}

\references{
  H.Akaike and T.Nakagawa (1988) \emph{Statistical Analysis and Control
    of Dynamic Systems.} Kluwer Academic publishers.
}

\examples{
  # Example 1 
  y <- rnorm(1000)
  dim(y) <- c(500,2)
  mulcor(y, lag_axis=FALSE)

  # Example 2
  xorg <- rnorm(1003)
  x <- matrix(0,1000,2)
  x[,1] <- xorg[1:1000]
  x[,2] <- xorg[4:1003]+0.5*rnorm(1000)
  mulcor(x, lag=20)
}

\keyword{ts}
