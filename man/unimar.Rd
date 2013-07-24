\name{unimar}
\alias{unimar}
\title{Univariate Case of Minimum AIC Method of AR Model Fitting}
\description{
  This is the basic program for the fitting of autoregressive models of successively higher
  by the method of least squares realized through householder transformation. 
}
\usage{
unimar(y, max.order=NULL, plot=FALSE)
}
\arguments{
  \item{y}{a univariate time series.}
  \item{max.order}{upper limit of AR order. Default is \eqn{2 \sqrt{n}}{2*sqrt(n)}, where \eqn{n} is the length of the time series \eqn{y}.}
  \item{plot}{logical. If TRUE \code{daic} is plotted.}
}
\value{
  \item{mean}{mean.}
  \item{var}{variance.}
  \item{v}{innovation variance.}
  \item{aic}{AIC.}
  \item{aicmin}{minimum AIC.}
  \item{daic}{AIC-\code{aicmin}.}
  \item{order.maice}{order of minimum AIC.}
  \item{v.maice}{innovation variance attained at \code{order.maice}.}
  \item{arcoef}{AR coefficients.}
}
\details{
  The AR model is given by
    \deqn{y(t) = a(1)y(t-1) + \ldots + a(p)y(t-p) + u(t),}{y(t) = a(1)y(t-1) + ... + a(p)y(t-p) + u(t),}
  where \eqn{p} is AR order and \eqn{u(t)} is Gaussian white noise with mean \eqn{0} and variance \eqn{v}.
  AIC is defined by
    \deqn{AIC = n\log(det(v)) + 2k,}{AIC = nlog(det(v)) + 2k,}
  where \eqn{n} is the length of data, \eqn{v} is the estimates of the innovation variance
  and \eqn{k} is the number of parameter.
}

\references{
  G.Kitagawa and H.Akaike (1978) A Procedure For The Modeling of Non-Stationary Time Series.
  Ann. Inst. Statist. Math.,30, B, 351--363.

  H.Akaike, G.Kitagawa, E.Arahata and F.Tada (1979)
  \emph{Computer Science Monograph, No.11, Timsac78.}
  The Institute of Statistical Mathematics.
}

\examples{
  data(Canadianlynx)
  z <- unimar(Canadianlynx, max.order=20)
  z$arcoef
}

\keyword{ts}
