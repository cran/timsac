\name{mulspe}
\alias{mulspe}
\alias{ptint.mulspe}
\title{Multiple Spectrum}
\description{
  Compute multiple spectrum estimates using Akaike window or Hanning window.
}
\usage{
  mulspe(y, lag=NULL, window="Akaike", plot=TRUE, plot.scale=FALSE)
}
\arguments{
  \item{y}{a multivariate time series with d variables and n observations.
           (y[n,d])}
  \item{lag}{maximum lag. Default is 2*sqrt(n), where n is the number of observations.}
  \item{window}{character string giving the definition of smoothing window. Allowed values are "Akaike" (default) or "Hanning".}
  \item{plot}{logical. If TRUE (default) spectrums are plotted as (d,d) matrix.
              Diagonal parts: Auto spectrums for each series.
              Lower triangular parts: Amplitude spectrums.
              Upper triangular part: Pahse spectrums.}
  \item{plot.scale}{logical. IF TRUE the common range of the y-axisis is used.}
}
\details{
  Hanning Window : a1(0)=0.5, a1(1)=a1(-1)=0.25, a1(2)=a1(-2)=0

  Akaike Window  : a2(0)=0.625, a2(1)=a2(-1)=0.25, a2(2)=a2(-2)=-0.0625
}
\value{
  \item{spec}{spectrum smoothing by "window".
              Lower triangular parts: Real parts.
              Upper triangular parts: Imaginary parts.}
  \item{stat}{test statistics.}
  \item{coh}{simple coherence by "window".}
}

\references{
  H.Akaike and T.Nakagawa (1988) \emph{Statistical Analysis and Control
    of Dynamic Systems.} Kluwer Academic publishers.
}

\examples{
  sgnl <- rnorm(1003)
  x <- matrix(0,1000,2)
  x[,1] <- sgnl[4:1003]
  x[,2] <- 0.9*sgnl[1:1000]+0.2*rnorm(1000)  #x[i,2]=0.9*x[i-3,1]+0.2*N(0,1)
  mulspe(x, 100, window="Hanning", plot=TRUE, plot.scale=TRUE)
}

\keyword{ts}
