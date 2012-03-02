\name{fpeaut}
\alias{fpeaut}
\title{FPE Auto}
\description{
  Perform FPE(Final Prediction Error) computation for one-dimensional AR model.}
\usage{
  fpeaut(y, max.order=NULL)
}
\arguments{
  \item{y}{a univariate time series.}
  \item{max.order}{upper limit of model order. Default is \eqn{2 \sqrt{n}}{2*sqrt(n)}, where \eqn{n} is the length of the time series \code{y}.}
}
\value{
  \item{ordermin}{order of minimum FPE.}
  \item{best.ar}{AR coefficients with minimum FPE.}
  \item{sigma2m}{= \code{sigma2}(\code{ordermin}).}
  \item{fpemin}{minimum FPE.}
  \item{rfpemin}{minimum RFPE.}
  \item{ofpe}{OFPE.}
  \item{arcoef}{AR coefficients.}
  \item{sigma2}{\eqn{\sigma^2.}{sigma^2.}}
  \item{fpe}{FPE (Final Prediction Error).}
  \item{rfpe}{RFPE.}
  \item{parcor}{partial correlation.}
  \item{chi2}{chi-squared.}
}
\details{
  The AR model is given by
    \deqn{y(t) = a(1)y(t-1) + .... + a(p)y(t-p) + u(t)}
  where \eqn{p} is AR order and \eqn{u(t)} is a zero mean white noise.
}

\references{
  H.Akaike and T.Nakagawa (1988) \emph{Statistical Analysis and Control
    of Dynamic Systems.} Kluwer Academic publishers.
}

\examples{
  y <- arima.sim(list(order=c(2,0,0), ar=c(0.64,-0.8)), n=200)
  fpeaut(y, max.order=20)
}

\keyword{ts}
