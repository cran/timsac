\name{armafit}
\alias{armafit}
\title{ARMA Model Fitting}
\description{
  Fit an ARMA model with specified order by using DAVIDON's algorithm.
}
\usage{
  armafit(y, model.order)
}
\arguments{
  \item{y}{a univariate time series.}
  \item{model.order}{a numerical vector of the form c(ar, ma) which gives the order to be fitted successively. }
}
\details{
  The maximum likelihood estimates of the coefficients of a scalar ARMA model
    \deqn{y(t) - a(1)y(t-1) -...- a(p)y(t-p) = u(t) - b(1)u(t-1) -...- b(q)u(t-q)}
  of a time series \eqn{y(t)} are obtained by using DAVIDON's algorithm.
  Pure autoregression is not allowed.
}
\value{
  \item{arcoef}{maximum likelihood estimates of AR coefficients.}
  \item{macoef}{maximum likelihood estimates of MA coefficients.}
  \item{arstd}{standard deviation (AR).}
  \item{mastd}{standard deviation (MA).}
  \item{v}{innovation variance.}
  \item{aic}{AIC.}
  \item{grad}{final gradient.}
}

\references{
  H.Akaike, E.Arahata and T.Ozaki (1975) \emph{Computer Science Monograph, No.5,
   Timsac74, A Time Series Analysis and Control Program Package (1)}.
  The Institute of Statistical Mathematics.
}

\examples{
  # "arima.sim" is a function in "stats".
  # Note that the sign of MA coefficient is opposite from that in "timsac".
  y <- arima.sim(list(order=c(2,0,1), ar = c(0.64,-0.8), ma=-0.5), 1000)
  z <- armafit(y, model.order=c(2,1))
  z$arcoef
  z$macoef
}

\keyword{ts}
