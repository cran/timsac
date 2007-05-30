\name{bispec}
\alias{bispec}
\title{Bispectrum}
\description{
  Compute bi-spectrum using the direct Fourier transform of sample third order monents.
}
\usage{
bispec(y, lag=NULL, window="Akaike", log=FALSE, plot=TRUE)
}
\arguments{
  \item{y}{a univariate time series.}
  \item{lag}{maximum lag. Default is 2*sqrt(n), where n is the length of the time series y.}
  \item{window}{character string giving the definition of smoothing window. Allowed values are "Akaike" (default) or "Hanning".}
  \item{log}{logical. If TRUE the spectrum pspec is plotted as log(pspec).}
  \item{plot}{logical. If TRUE (default) the spectrum pspec is plotted.}
}
\details{
  Hanning Window : a1(0)=0.5, a1(1)=a1(-1)=0.25, a1(2)=a1(-2)=0

  Akaike Window  : a2(0)=0.625, a2(1)=a2(-1)=0.25, a2(2)=a2(-2)=-0.0625
}
\value{
  \item{pspec}{power spectrum smoothed by "window".}
  \item{sig}{significance.}
  \item{cohe}{coherence.}
  \item{breal}{real part of bispectrum.}
  \item{bimag}{imaginary part of bispectrum.}
  \item{exval}{aproximate expected value of coherence under Gaussian assumption.}
}

\references{
  H.Akaike, E.Arahata and T.Ozaki (1975) \emph{Computer Science Monograph, No.6,
   Timsac74, A Time Series Analysis and Control Program Package (2)}.
  The Institute of Statistical Mathematics.
}

\examples{
  data(bispecData)
  bispec(bispecData, lag=30)
}

\keyword{ts}
