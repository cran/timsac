\name{auspec}
\alias{auspec}
\title{Power Spectrum}
\description{
  Compute power spectrum estimates for two trigonometric windows of
  Blackman-Tukey type by Goertzel method.
}
\usage{
auspec(y, lag = NULL, window = "Akaike", log = FALSE, plot = TRUE)
}
\arguments{
  \item{y}{a univariate time series.}
  \item{lag}{maximum lag. Default is \eqn{2 \sqrt{n}}{2*sqrt(n)}, where \eqn{n}
    is the length of time series \code{y}.}
  \item{window}{character string giving the definition of smoothing window.
    Allowed strings are "Akaike" (default) or "Hanning".}
  \item{log}{logical. If \code{TRUE}, the spectrum \code{spec} is plotted as
    log(\code{spec}).}
  \item{plot}{logical. If \code{TRUE} (default), the spectrum \code{spec} is
    plotted.}
}
\details{
  \tabular{rlll}{
    Hanning Window : \tab a1(0)=0.5, \tab a1(1)=a1(-1)=0.25, \tab
    a1(2)=a1(-2)=0\cr
    Akaike Window  : \tab a2(0)=0.625, \tab a2(1)=a2(-1)=0.25, \tab
    a2(2)=a2(-2)=-0.0625
  }
}
\value{
  \item{spec}{spectrum smoothing by '\code{window}'}
  \item{stat}{test statistics.}
}

\references{
  H.Akaike and T.Nakagawa (1988) \emph{Statistical Analysis and Control
    of Dynamic Systems.} Kluwer Academic publishers.
}

\examples{
y <- arima.sim(list(order=c(2,0,0), ar=c(0.64,-0.8)), n = 200)
auspec(y, log = TRUE)
}

\keyword{ts}
