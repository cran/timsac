\name{auspec}
\alias{auspec}
\title{Power Spectrum}
\description{
  Compute power spectrum estimates for two trigonometric windows of Blackman-Tukey type by goertzel method.
}
\usage{
auspec(y, lag=NULL, window="Akaike", log=FALSE, plot=TRUE)
}
\arguments{
  \item{y}{a univariate time series.}
  \item{lag}{maximum lag. Default is 2*sqrt(n), where n is the length of time series y.}
  \item{window}{character string giving the definition of smoothing window. Allowed values are "Akaike" (default) or "Hanning".}
  \item{log}{logical. If TRUE, the spectrum spec is plotted as log(spec).}
  \item{plot}{logical. If TRUE (default) the spectrum is plotted.}
}
\details{
  Hanning Window : a1(0)=0.5, a1(1)=a1(-1)=0.25, a1(2)=a1(-2)=0

  Akaike Window  : a2(0)=0.625, a2(1)=a2(-1)=0.25, a2(2)=a2(-2)=-0.0625
}
\value{
  \item{spec}{spectrum smoothing by "window"}
  \item{stat}{test statistics.}
}

\references{
  H.Akaike and T.Nakagawa (1988) \emph{Statistical Analysis and Control
    of Dynamic Systems.} Kluwer Academic publishers.
}

\examples{
y <- arima.sim(list(order=c(2,0,0), ar=c(0.64,-0.8)), n=200)
auspec(y, log=TRUE)
}

\keyword{ts}
