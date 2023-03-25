\name{bispec}
\alias{bispec}
\title{Bispectrum}
\description{
  Compute bi-spectrum using the direct Fourier transform of sample third order
  moments.
}
\usage{
bispec(y, lag = NULL, window = "Akaike", log = FALSE, plot = TRUE)
}
\arguments{
  \item{y}{a univariate time series.}
  \item{lag}{maximum lag. Default is \eqn{2 \sqrt{n}}{2*sqrt(n)}, where \eqn{n}
    is the length of the time series \code{y}.}
  \item{window}{character string giving the definition of smoothing window.
    Allowed strings are "Akaike" (default) or "Hanning".}
  \item{log}{logical. If \code{TRUE}, the spectrum \code{pspec} is plotted as
    \eqn{log}(\code{pspec}).}
  \item{plot}{logical. If \code{TRUE} (default), the spectrum \code{pspec} is
    plotted.}
}
\details{
  \tabular{rlll}{
    Hanning Window : \tab a1(0)=0.5,\tab a1(1)=a1(-1)=0.25,\tab
    a1(2)=a1(-2)=0\cr
    Akaike Window : \tab a2(0)=0.625,\tab a2(1)=a2(-1)=0.25,\tab
    a2(2)=a2(-2)=-0.0625
  }
}
\value{
  \item{pspec}{power spectrum smoothed by '\code{window}'.}
  \item{sig}{significance.}
  \item{cohe}{coherence.}
  \item{breal}{real part of bispectrum.}
  \item{bimag}{imaginary part of bispectrum.}
  \item{exval}{approximate expected value of coherence under Gaussian assumption.}
}

\references{
  H.Akaike, E.Arahata and T.Ozaki (1975) \emph{Computer Science Monograph, No.6,
   Timsac74, A Time Series Analysis and Control Program Package (2)}.
  The Institute of Statistical Mathematics.
}

\examples{
data(bispecData)
bispec(bispecData, lag = 30)
}

\keyword{ts}
