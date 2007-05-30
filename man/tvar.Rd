\name{tvar}
\alias{tvar}
\title{Time Varying Coefficients AR model}
\description{
  Estimate time varying coefficients AR model.
}
\usage{
  tvar(y, ar.order, trend.order=2, span, outlier, tau20=NULL, delta=NULL, plot=TRUE)
}

\arguments{
  \item{y}{a univariate time series.}
  \item{ar.order}{AR order.}
  \item{trend.order}{trend order (=1 or 2).}
  \item{span}{local stationary span.}
  \item{outlier}{positions of outliers.}
  \item{tau20}{initial value for computing variance of the system noise tau2.}
  \item{delta}{delta for computing variance of the system noise tau2. If tau20 is NULL or delta is NULL, tau2 is computed automatically.}
  \item{plot}{logical. If TRUE (default) parcor is plotted.}
}
\details{
  The time-varying coefficients  AR model is given by

    \eqn{y(t) = a(1,t)y(t-1) + ... + a(p,t)y(t-p) + u(t)}

    where \eqn{a(i,t)} is i-lag AR coefficient at time t and \eqn{u(t)} is a zwro mean white noise.
}

\value{
  \item{tau2}{variance of the system noise.}
  \item{sigma2}{variance of the observational noise.}
  \item{lkhood}{log-likelihood.}
  \item{aic}{AIC.}
  \item{arcoef}{time varying AR coefficients.}
  \item{parcor}{partial autocorrelation coefficient.}
  \item{spec}{time varying spectrum.}
}

\references{
  Kitagawa, G. (1993) \emph{Time series analysis programing (in Japanese)}. The Iwanami Computer Science Senes.

  Kitagawa, G. and Gersch, W. (1996) \emph{Smoothness Priors Analysis of Time Series}. Lecture Notes in Statistics, No.116, Springer-Verlag.

  Kitagawa, G. and Gersch, W. (1985) \emph{A smoothness priors time varying AR coefficient modeling of nonstationary time series}. IEEE trans. on Automatic Controle, AC-30, 48-56.
}

\examples{
  data(MYE1F) # an earthquake wave data
  z <- tvvar(MYE1F, trend.order=2, tau20= 6.6e-06, delta=1.0e-06, plot=FALSE)
  zz <- tvar(z$normdat, ar.order=4, trend.order=2, span=20, tau20=6.6e-06, delta=1.0e-06, outlier=c(630,1026))
  zz$tau2
  zz$sigma2
  zz$lkhood
  zz$aic
}

\keyword{ts}
