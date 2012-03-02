\name{exsar}
\alias{exsar}
\title{Exact Maximum Likelihood Method of Scalar AR Model Fitting}
\description{
  Produce exact maximum likelihood estimates of the parameters of a scalar AR model.
}
\usage{
  exsar(y, max.order=NULL, plot=FALSE)
}
\arguments{
  \item{y}{a univariate time series.}
  \item{max.order}{upper limit of AR order. Default is \eqn{2 \sqrt{n}}{2*sqrt(n)}, where \eqn{n} is the length of the time series \code{y}.}
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
  \item{v.maice}{MAICE innovation variance.}
  \item{arcoef.maice}{MAICE AR coefficients.}
  \item{v.mle}{maximum likelihood estimates of innovation variance.}
  \item{arcoef.mle}{maximum likelihood estimates of AR coefficients.}
}
\details{
  The AR model is given by

    \deqn{y(t) = a(1)y(t-1) + .... + a(p)y(t-p) + u(t)}

  where \eqn{p} is AR order and \eqn{u(t)} is a zero mean white noise.
}

\references{
  H.Akaike, G.Kitagawa, E.Arahata and F.Tada (1979)
  \emph{Computer Science Monograph, No.11, Timsac78.}
  The Institute of Statistical Mathematics.
}

\examples{
  data(Canadianlynx)
  z <- exsar(Canadianlynx, max.order=14)
  z$arcoef.maice
  z$arcoef.mle
}

\keyword{ts}
