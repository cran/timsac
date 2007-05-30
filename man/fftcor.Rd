\name{fftcor}
\alias{fftcor}
\title{Auto And/Or Cross Correlations via FFT}
\description{
  Compute auto and/or cross covariances and correlations via FFT.
}
\usage{
  fftcor(y, lag=NULL, isw=4, plot=TRUE, lag_axis=TRUE)
}
\arguments{
  \item{y}{data of channel X and Y (data of channel Y is given for isw=2 or 4 only).}
  \item{lag}{maximum lag. Default is 2*sqrt(n), where n is the length of the time series y.}
  \item{isw}{numerical flag giving the type of computation.

    isw = 1 : autocorrelation of X (one-channel)

    isw = 2 : autocorrelations of X and Y (two-channel)

    isw = 4 : auto- and cross- correlations of X and Y (two-channel)}
  \item{plot}{logical. If TRUE (default) crosscorrelations are plotted.}
  \item{lag_axis}{logical. If TRUE (default) with plot=TRUE, x\_axis is drawn.}
}
\value{
  \item{acov}{autocovariance.}
  \item{ccov12}{crosscovariance.}
  \item{ccov21}{crosscovariance.}
  \item{acor}{autocorrelation.}
  \item{ccor12}{crosscorrelation.}
  \item{ccor21}{crosscorrelation.}
  \item{mean}{mean.}
}

\references{
  H.Akaike and T.Nakagawa (1988) \emph{Statistical Analysis and Control
    of Dynamic Systems.} Kluwer Academic publishers.
}

\examples{
  # Example 1
  x <- rnorm(200)
  y <- rnorm(200)
  xy <- array(c(x,y), dim=c(200,2))
  fftcor(xy, lag_axis=FALSE)

  # Example 2
  xorg <- rnorm(1003)
  x <- matrix(0,1000,2)
  x[,1] <- xorg[1:1000]
  x[,2] <- xorg[4:1003]+0.5*rnorm(1000)
  fftcor(x, lag=20)
}

\keyword{ts}
