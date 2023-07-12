\name{sglfre}
\alias{sglfre}
\title{Frequency Response Function (Single Channel)}
\description{
  Compute 1-input,1-output frequency response function, gain, phase, coherency
  and relative error statistics.
}
\usage{
  sglfre(y, lag = NULL, invar, outvar)
}
\arguments{
  \item{y}{a multivariate time series.}
  \item{lag}{maximum lag. Default \eqn{2 \sqrt{n}}{2*sqrt(n)}, where \eqn{n} is
    the length of the time series \code{y}.}
  \item{invar}{within \eqn{d} variables of the spectrum, \code{invar}-th
    variable is taken as an input variable.}
  \item{outvar}{within \eqn{d} variables of the spectrum, \code{outvar}-th
    variable is taken as an output variable .}
}
\value{
  \item{inspec}{power spectrum (input).}
  \item{outspec}{power spectrum (output).}
  \item{cspec}{co-spectrum.}
  \item{qspec}{quad-spectrum.}
  \item{gain}{gain.}
  \item{coh}{coherency.}
  \item{freqr}{frequency response function : real part.}
  \item{freqi}{frequency response function : imaginary part.}
  \item{errstat}{relative error statistics.}
  \item{phase}{phase.}
}

\references{
  H.Akaike and T.Nakagawa (1988) \emph{Statistical Analysis and Control
    of Dynamic Systems.} Kluwer Academic publishers.
}
\examples{
ar <- array(0, dim = c(3,3,2))
ar[, , 1] <- matrix(c(0.4,  0,   0.3,
                      0.2, -0.1, -0.5,
                      0.3,  0.1,  0), nrow = 3, ncol = 3, byrow = TRUE)
ar[, , 2] <- matrix(c(0,  -0.3,  0.5,
                      0.7, -0.4,  1,
                      0,   -0.5,  0.3), nrow = 3, ncol = 3, byrow = TRUE)
x <- matrix(rnorm(200*3), nrow = 200, ncol = 3)
y <- mfilter(x, ar, "recursive")
sglfre(y, lag = 20, invar = 1, outvar = 2)
}

\keyword{ts}
