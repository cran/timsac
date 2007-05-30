\name{xsarma}
\alias{xsarma}
\title{Exact Maximum Likelihood Method of Scalar ARMA Model Fitting}
\description{
  Produce exact maximum likelihood estimates of the parameters of a scalar ARMA model.
}
\usage{
  xsarma(y, arcoefi, macoefi)
}
\arguments{
  \item{y}{a univariate time series.}
  \item{arcoefi}{initial estimates of AR coefficients.}
  \item{macoefi}{initial estimates of MA coefficients.}
}
\value{
  \item{gradi}{initial gradient.}
  \item{lkhoodi}{initial (-2)log likelihood.}
  \item{arcoef}{final estimates of AR coefficients.}
  \item{macoef}{final estimates of MA coefficients.}
  \item{grad}{final gradient.}
  \item{alph.ar}{final ALPH (AR part) at subroutine ARCHCK.}
  \item{alph.ma}{final ALPH (MA part) at subroutine ARCHCK.}
  \item{lkhood}{final (-2)log likelihood.}
  \item{wnoise.var}{white noise variance.}
}
\details{
  The ARMA model is given by

    \eqn{y(t) - a(1)y(t-1) - ... - a(p)y(t-p) = u(t) - b(1)u(t-1) - ... - b(q)u(t-q)},

  where \code{p} is AR order, \code{q} is MA order and \code{u(t)} is a zero mean white noise.
}

\references{
  H.Akaike (1978)
  Covariance matrix computation of the state variable of a stationary Gaussian process.
   Research Memo. No.139. The Institute of Statistical Mathematics.

  H.Akaike, G.Kitagawa, E.Arahata and F.Tada (1979)
  \emph{Computer Science Monograph, No.11, Timsac78.}
  The Institute of Statistical Mathematics.
}

\examples{
  # "arima.sim" is a function in "stats".
  # Note that the sign of MA coefficient is opposite from that in "timsac".
  arcoef <- c(1.45, -0.9)
  macoef <- c(-0.5)
  y <- arima.sim(list(order=c(2,0,1),ar=arcoef,ma=macoef),n=100)
  arcoefi <- c(1.5, -0.8)
  macoefi <- c(0.0)
  z <- xsarma(y, arcoefi, macoefi)
  z$arcoef
  z$macoef
}

\keyword{ts}
