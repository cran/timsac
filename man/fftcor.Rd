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
  \item{y}{data of channel X and Y (data of channel Y is given for \code{isw} = 2 or 4 only).}
  \item{lag}{maximum lag. Default is \eqn{2 \sqrt{n}}{2*sqrt(n)}, where \eqn{n} is the length of the time series \code{y}.}
  \item{isw}{numerical flag giving the type of computation.
    \tabular{rl}{
      1 : \tab auto-correlation of X (one-channel)\cr
      2 : \tab auto-correlations of X and Y (two-channel)\cr
      4 : \tab auto- and cross- correlations of X and Y (two-channel)
    }}
  \item{plot}{logical. If TRUE (default) cross-correlations are plotted.}
  \item{lag_axis}{logical. If TRUE (default) with plot=TRUE, \eqn{x}-axis is drawn.}
}
\value{
  \item{acov}{auto-covariance.}
  \item{ccov12}{cross-covariance.}
  \item{ccov21}{cross-covariance.}
  \item{acor}{auto-correlation.}
  \item{ccor12}{cross-correlation.}
  \item{ccor21}{cross-correlation.}
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
