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
  \item{y}{a multivariate time series with \eqn{d} variables and \eqn{n} observations.
           (\eqn{y[n,d])}}
  \item{lag}{maximum lag. Default is \eqn{2 \sqrt{n}}{2*sqrt(n)}, where \eqn{n} is the number of observations.}
  \item{window}{character string giving the definition of smoothing window. Allowed values are "Akaike" (default) or "Hanning".}
  \item{plot}{logical. If TRUE (default) spectrums are plotted as \eqn{(d,d)} matrix.
    \tabular{ll}{
      Diagonal parts : \tab Auto spectrums for each series.\cr
      Lower triangular parts : \tab Amplitude spectrums.\cr
      Upper triangular part : \tab Phase spectrums.
    }}
  \item{plot.scale}{logical. IF TRUE the common range of the \eqn{y}-axis is used.}
}
\details{
  \tabular{rlll}{
    Hanning Window : \tab a1(0)=0.5, \tab a1(1)=a1(-1)=0.25, \tab a1(2)=a1(-2)=0\cr
    Akaike Window : \tab a2(0)=0.625, \tab a2(1)=a2(-1)=0.25, \tab a2(2)=a2(-2)=-0.0625
  }
}
\value{
  \item{spec}{spectrum smoothing by "\code{window}".
    \tabular{rl}{
      Lower triangular parts : \tab Real parts\cr
      Upper triangular parts : \tab Imaginary parts
    }}
  \item{stat}{test statistics.}
  \item{coh}{simple coherence by "\code{window}".}
}

\references{
  H.Akaike and T.Nakagawa (1988) \emph{Statistical Analysis and Control
    of Dynamic Systems.} Kluwer Academic publishers.
}

\examples{
  sgnl <- rnorm(1003)
  x <- matrix(0,1000,2)
  x[,1] <- sgnl[4:1003]
  #x[i,2]=0.9*x[i-3,1]+0.2*N(0,1)
  x[,2] <- 0.9*sgnl[1:1000]+0.2*rnorm(1000)
  mulspe(x, 100, "Hanning", plot.scale=TRUE)
}

\keyword{ts}
