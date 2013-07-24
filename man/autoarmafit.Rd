\name{autoarmafit}
\alias{autoarmafit}
\alias{print.autoarmafit}
\title{Automatic ARMA Model Fitting}
\description{
  Provide an automatic ARMA model fitting procedure.
  Models with various orders are fitted and the best choice is determined with the aid of
 the statistics AIC.
}
\usage{
  autoarmafit(y, max.order=NULL)
}
\arguments{
  \item{y}{a univariate time series.}
  \item{max.order}{upper limit of AR order and MA order. Default is \eqn{2 \sqrt{n}}{2*sqrt(n)}, where \eqn{n} is the length of the time series \code{y}.}
}
\details{
  The maximum likelihood estimates of the coefficients of a scalar ARMA model
    \deqn{y(t) - a(1)y(t-1) -...- a(p)y(t-p) = u(t) - b(1)u(t-1) -...- b(q)u(t-q)}
  of a time series \eqn{y(t)} are obtained by using DAVIDON's variance algorithm.
  Where \eqn{p} is AR order, \eqn{q} is MA order and \eqn{u(t)} is a zero mean white noise.
  Pure autoregression is not allowed.
}
\value{
  \item{best.model}{the best choice of ARMA coefficients.}
  \item{model}{a list with components \code{arcoef} (Maximum likelihood estimates of AR coefficients), \code{macoef} (Maximum likelihood estimates of MA coefficients),
  \code{arstd} (AR standard deviation), \code{mastd} (MA standard deviation), \code{v} (Innovation variance), \code{aic} (AIC \eqn{= n \log(det(v))+2(p+q)}{= n log(det(v))+2(p+q)})
  and \code{grad} (Final gradient) in AIC increasing order.}
}

\references{
  H.Akaike, E.Arahata and T.Ozaki (1975) \emph{Computer Science Monograph, No.5,
   Timsac74, A Time Series Analysis and Control Program Package (1)}.
  The Institute of Statistical Mathematics.
}

\examples{
  # "arima.sim" is a function in "stats".
  # Note that the sign of MA coefficient is opposite from that in "timsac".
  y <- arima.sim(list(order=c(2,0,1),ar=c(0.64,-0.8),ma=-0.5),n=1000)
  autoarmafit(y)
}

\keyword{ts}
