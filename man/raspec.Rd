\name{raspec}
\alias{raspec}
\title{Rational Spectrum}
\description{
  Compute power spectrum of ARMA process.
}
\usage{
  raspec(h, var, arcoef=NULL, macoef=NULL, log=FALSE, plot=TRUE)
}
\arguments{
  \item{h}{specify frequencies \eqn{i/2h (i=0,1,...,h)}.}
  \item{var}{variance.}
  \item{arcoef}{AR coefficients.}
  \item{macoef}{MA coefficients.}
  \item{log}{logical. If TRUE the spectrum is plotted as log(sepc).}
  \item{plot}{logical. If TRUE (default) the spectrum is plotted.}
}
\value{
  \code{raspec} gives the rational spectrum.
}
\details{
  ARMA process :

  \eqn{y(t) - a(1)y(t-1) -...- a(p)y(t-p) = u(t) - b(1)u(t-1) -...- b(q)u(t-q)}

  where p is AR order, q is MA order and u(t) is a white noise with zero mean and variance equal to \code{var}.
}
\references{
  H.Akaike and T.Nakagawa (1988) \emph{Statistical Analysis and Control
    of Dynamic Systems.} Kluwer Academic publishers.
}

\examples{
  # Example 1 for the AR model
  raspec(h=100, var=1, arcoef=c(0.64,-0.8))

  # Example 2 for the MA model
  raspec(h=20, var=1, macoef=c(0.64,-0.8))
}

\keyword{ts}
