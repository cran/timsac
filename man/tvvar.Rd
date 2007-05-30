\name{tvvar}
\alias{tvvar}
\title{Time Varying Variance}
\description{
  Estimate time-varying variance.
}
\usage{
  tvvar(y, trend.order, tau20, delta, plot=TRUE)
}

\arguments{
  \item{y}{univariate time series.}
  \item{trend.order}{trend order.}
  \item{tau20}{initial estimate of tau2.}
  \item{delta}{search width.}
  \item{plot}{logical. If TRUE (default) normdat, ts, trend and noise are plotted.}
}

\details{
    A chi-square distribution with digree 2 is given by

    \eqn{s(m) = y(2m-1)**2 + y(2m)**2}

    where \eqn{y(n)} is original scalar time series and \eqn{\sigma(2m-1)**2 = \sigma(2m)**2}.

    \eqn{z(m) = log(s(m)/2)}.

    \eqn{z(m) = log(\sigma**2) + w(m)},

    where \eqn{w(m)} is a double exponential distribution with density \eqn{h(w) = exp{w-e**w}}.

    The space state model is given by

    \eqn{z(m) = t(m) + w(m)}.
}

\value{
  \item{tvvar}{time varying variance.}
  \item{normdat}{normalized data.}
  \item{ts}{tranceformed time series s(m).}
  \item{trend}{trend.}
  \item{noise}{residuals.}
  \item{tau2}{variance of the system noise tau2.}
  \item{sigma2}{variance of the observational noise.}
  \item{lkhood}{log-likelihood of the mode.}
  \item{aic}{AIC.}
}

\references{
  Kitagawa, G. (1993) \emph{Time series analysis programing (in Japanese)}. The Iwanami Computer Science Senes.

  Kitagawa, G. and Gersch, W. (1996) \emph{Smoothness Priors Analysis of Time Series}. Lecture Notes in Statistics, No.116, Springer-Verlag.

  Kitagawa, G. and Gersch, W. (1985) \emph{A smoothness priors time varying AR coefficient modeling of nonstationary time series}. IEEE trans. on Automatic Controle, AC-30, 48-56.
}

\examples{
  data(MYE1F) # an earthquake wave data
  z <- tvvar(MYE1F, trend.order=2, tau20= 6.6e-06, delta=1.0e-06)
  z$lkhood
  z$aic
}

\keyword{ts}
