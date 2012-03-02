\name{tvar}
\alias{tvar}
\title{Time Varying Coefficients AR model}
\description{
  Estimate time varying coefficients AR model.
}
\usage{
  tvar(y, ar.order, trend.order=2, span, outlier=NULL, tau20=NULL,
       delta=NULL, plot=TRUE)
}

\arguments{
  \item{y}{a univariate time series.}
  \item{ar.order}{AR order.}
  \item{trend.order}{trend order (1 or 2).}
  \item{span}{local stationary span.}
  \item{outlier}{positions of outliers.}
  \item{tau20}{initial estimate of variance of the system noise \code{tau2}.}
  \item{delta}{search width. If \code{tau2} is NULL or \code{delta} is NULL, \code{tau2} is computed automatically.}
  \item{plot}{logical. If TRUE (default) \code{parcor} is plotted.}
}
\details{
  The time-varying coefficients  AR model is given by

    \deqn{y_t = a_{1,t}y_{t-1} + \ldots + a_{p,t}y_{t-p} + u_t}{y(t) = a(1,t)y(t-1) + ... + a(p,t)y(t-p) + u(t)}

    where \eqn{a_{i,t}}{a(i,t)} is \eqn{i}-lag AR coefficient at time \eqn{t} and \eqn{u_t}{u(t)} is a zero mean white noise.
}

\value{
  \item{tau2max}{variance of the system noise for maximum log-likelihood.}
  \item{sigma2}{variance of the observational noise.}
  \item{lkhood}{log-likelihood.}
  \item{aic}{AIC.}
  \item{arcoef}{time varying AR coefficients.}
  \item{parcor}{partial autocorrelation coefficient.}
}

\references{
  Kitagawa, G. (1993) \emph{Time series analysis programing (in Japanese)}. The Iwanami Computer Science Series.

  Kitagawa, G. and Gersch, W. (1996) \emph{Smoothness Priors Analysis of Time Series}. Lecture Notes in Statistics, No.116, Springer-Verlag.

  Kitagawa, G. and Gersch, W. (1985) \emph{A smoothness priors time varying AR coefficient modeling of nonstationary time series}. IEEE trans. on Automatic Control, AC-30, 48-56.
}

\examples{
  data(MYE1F) # an earthquake wave data
  z <- tvar(MYE1F, 4, 2, 20, c(630,1026), 6.6e-06, 1.0e-06)
  z$tau2max
  z$sigma2
  z$lkhood
  z$aic
}

\keyword{ts}
