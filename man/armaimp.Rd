\name{armaimp}
\alias{armaimp}
\title{Calculate Characteristics of Scalar ARMA Model}
\description{
  Calculate impulse, autocovariance, partial autocorrelation function and characteristic roots of scalar ARMA model for given AR and MA coefficients.
}
\usage{
  armaimp(arcoef, macoef, v, n=1000, lag=NULL, nf=200, plot=TRUE)
}
\arguments{
  \item{arcoef}{AR coefficients.}
  \item{macoef}{MA coefficients.}
  \item{v}{innovation variance.}
  \item{n}{data length.}
  \item{lag}{maximum lag of autocovariance function. Default is \eqn{2 \sqrt{n}}{2*sqrt(n)}.}
  \item{nf}{number of frequencies in evaluating spectrum.}
  \item{plot}{logical. If TRUE (default) impulse response function, autocovariance, power spectrum and characteristic roots are plotted.}
}
\details{
  The ARMA model is given by

    \deqn{y(t) - a(1)y(t-1) - ... - a(p)y(t-p) = u(t) - b(1)u(t-1) - ... - b(q)u(t-q),}

  where \eqn{p} is AR order, \eqn{q} is MA order and \eqn{u(t)} is a zero mean white noise.
}
\value{
  \item{impuls}{impulse response function.}
  \item{acov}{autocovariance function.}
  \item{parcor}{partial autocorrelation function.}
  \item{spec}{power spectrum.}
  \item{croot.ar}{characteristic roots of AR operator. Characteristic root is a list with components named real (real part \eqn{R}),
       image (imaginary part \eqn{I}), amp (\eqn{=\sqrt{R^2+I^2}}{=sqrt(R**2+I**2)}), atan(\eqn{=atan(I/R)}) and degree.}
  \item{croot.ma}{characteristic roots of MA operator.}
}

\references{
  G.Kitagawa (1993) \emph{Time series analysis programing (in Japanese)}. The Iwanami Computer Science Series.
}

\examples{
  # ARMA model :  y(n) = 0.9sqrt(3)y(n-1) - 0.81y(n-2)
  #                      + v(n) -0.9sqrt(2)v(n-1) + 0.81v(n-2)
  a <- c(0.9*sqrt(3), -0.81)
  b <- c(0.9*sqrt(2), -0.81)
  z <- armaimp(arcoef=a, macoef=b, v=1.0, n=1000, lag=20)
  z$croot.ar
  z$croot.ma

  # AR model : y(n) = 0.9sqrt(3)y(n-1) - 0.81y(n-2) + v(n)
  z <- armaimp(arcoef=a, v=1.0, n=1000, lag=20)
  z$croot.ar

  # MA model : y(n) = v(n) -0.9sqrt(2)v(n-1) + 0.81v(n-2)
  z <- armaimp(macoef=b, v=1.0, n=1000, lag=20)
  z$croot.ma
}

\keyword{ts}
