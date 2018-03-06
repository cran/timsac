\name{covgen}
\alias{covgen}
\title{Covariance Generation}
\description{
  Produce the Fourier transform of a power gain function in the form of an
  autocovariance
  sequence.
}
\usage{
  covgen(lag, f, gain, plot = TRUE)
}
\arguments{
  \item{lag}{desired maximum lag of covariance.}
  \item{f}{frequency \code{f[i]} \eqn{(i=1,...,k)}, where \eqn{k} is the number
    of data points. By definition \code{f[1]} = 0.0 and \code{f[k]} = 0.5,
    \code{f[i]}'s are arranged in increasing order.}
  \item{gain}{power gain of the filter at the frequency \code{f[i]}.}
  \item{plot}{logical. If \code{TRUE} (default), autocorrelations are plotted.}
}
\value{
  \item{acov}{autocovariance.}
  \item{acor}{autocovariance normalized.}
}

\references{
  H.Akaike, E.Arahata and T.Ozaki (1975) \emph{Computer Science Monograph, No.5,
   Timsac74, A Time Series Analysis and Control Program Package (1)}.
  The Institute of Statistical Mathematics.
}

\examples{
spec <- raspec(h = 100, var = 1, arcoef = c(0.64,-0.8), plot = FALSE)
covgen(lag = 100, f = 0:100/200, gain = spec)
}

\keyword{ts}
