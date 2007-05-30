\name{canarm}
\alias{canarm}
\title{Canonical Correlation Analysis of Scalar Time Series}
\description{
  Fit an ARMA model to stationary scalar time series through the analysis of
  canonical correlations between the future and past sets of observations.
}
\usage{
canarm(y, max.order=NULL, plot=TRUE)
}
\arguments{
  \item{y}{a univariate time series.}
  \item{max.order}{upper limit of AR order and MA order. Default is 2*sqrt(n), where n is the length of the time series y.}
  \item{plot}{logical. If TRUE (default) parcor is plotted.}
}
\details{
  The ARMA model of stationary scalar time series y(t) (t=1,...,n) is given by

  \eqn{y(t) - a(1)y(t-1) -...- a(p)y(t-p) = u(t) - b(1)u(t-1) -...- b(q)u(t-q)},

  where p is AR order and q is MA order.
}
\value{
  \item{arinit}{AR coefficients of initial AR model fitting by the minimum AIC procedure.}
  \item{v}{innovation vector.}
  \item{aic}{AIC.}
  \item{aicmin}{minimum AIC.}
  \item{order.maice}{order of minimum AIC.}
  \item{parcor}{partial autocorrelation.}
  \item{nc}{total number of case.}
  \item{future}{number of present and future variables.}
  \item{past}{number of present and past variables.}
  \item{cweight}{future set canonical weight.}
  \item{canocoef}{canonical R.}
  \item{canocoef2}{R-squared.}
  \item{chisquar}{chi-square.}
  \item{ndf}{N.D.F.}
  \item{dic}{DIC.}
  \item{dicmin}{minimum DIC.}
  \item{order.dicmin}{order of minimum DIC.}
  \item{arcoef}{AR coefficients a(i) (i = 1,...,p).}
  \item{macoef}{MA coefficients b(i) (i = 1,...,q).}
}

\references{
  H.Akaike, E.Arahata and T.Ozaki (1975) \emph{Computer Science Monograph, No.5,
   Timsac74, A Time Series Analysis and Control Program Package (1)}.
  The Institute of Statistical Mathematics.
}

\examples{
  # "arima.sim" is a function in "stats".
  # Note that the sign of MA coefficient is opposite from that in "timsac".
  y <- arima.sim(list(order=c(2,0,1), ar=c(0.64,-0.8), ma=c(-0.5)), n=1000)
  z <- canarm(y, max.order=30)
  z$arcoef
  z$macoef
}

\keyword{ts}
