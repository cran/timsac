\name{decomp}
\alias{decomp}
\title{Time Series Decomposition (Seasonal Adjustment) by Square-Root Filter}
\description{
  Decompose a nonstationary time series into several possible components by square-root filter.
}
\usage{
  decomp(y, trend.order=2, ar.order=2, frequency=12,
         seasonal.order=1, log=FALSE, trade=FALSE, diff=1,
         year=1980, month=1, miss=0, omax=99999.9, plot=TRUE)
}
\arguments{
  \item{y}{a univariate time series.}
  \item{trend.order}{trend order (0, 1, 2 or 3).}
  \item{ar.order}{AR order (less than 11, try 2 first).}
  \item{frequency}{number of seasons in one period.}
  \item{seasonal.order}{seasonal order (0, 1 or 2).}
  \item{log}{log transformation of data (if \code{log} = TRUE).}
  \item{trade}{trading day adjustment (if \code{trade} = TRUE).}
  \item{diff}{numerical differencing (1 sided or 2 sided).}
  \item{year}{the first year of the data.}
  \item{month}{the first month of the data.}
  \item{miss}{missing data flag.
    \tabular{rl}{
      = 0 : \tab no consideration\cr
      > 0 : \tab values which are greater than \code{omax} are treated as missing data\cr
      < 0 : \tab values which are less than \code{omax} are treated as missing data
    }}
  \item{omax}{maximum or minimum data value (if \code{miss} > 0 or \code{miss} < 0).}
  \item{plot}{logical. If TRUE (default) \code{trend}, \code{seasonal}, \code{ar} and \code{trad} are plotted.}
}

\details{
  The Basic Model\cr
    \eqn{y(t) = T(t) + AR(t) + S(t) + TD(t) + W(t)}\cr
  where \eqn{T(t)} is trend component, \eqn{AR(t)} is AR process, \eqn{S(t)} is seasonal component,
  \eqn{TD(t)} is trading day factor and \eqn{W(t)} is observational noise. \cr

  Component Models

  Trend component  (trend.order m1)\cr
    \eqn{m1= 1 : T(t) = T(t-1) + V1(t)}\cr
    \eqn{m1= 2 : T(t) = 2T(t-1) - T(t-2) + V1(t)}\cr
    \eqn{m1= 3 : T(t) = 3T(t-1) -3T(t-2) + T(t-2) + V1(t)}\cr\cr

  AR component  (ar.order m2)\cr
    \eqn{AR(t) = a(1)AR(t-1) + \ldots + a(m2)AR(t-m2) + V2(t)}{AR(t) = a(1)AR(t-1) + ... + a(m2)AR(t-m2) + V2(t)}\cr\cr

  Seasonal component  (seasonal.order k, frequency f)\cr
    \eqn{k=1 : S(t) = -S(t-1) - \ldots - S(t-f+1) + V3(t)}{k=1 : S(t) = -S(t-1) - ... - S(t-f+1) + V3(t)}\cr
    \eqn{k=2 : S(t) = -2S(t-1) - \ldots -f\ S(t-f+1) - \ldots - S(t-2f+2) + V3(t)}{k=2 : S(t) = -2S(t-1) - ... -f S(t-f+1) - ... - S(t-2f+2) + V3(t)}\cr\cr

  Trading day effect\cr
    \eqn{TD(t) = b(1) TRADE(t,1) + \ldots + b(7) TRADE(t,7)}{TD(t) = b(1) TRADE(t,1) + ... + b(7) TRADE(t,7)}\cr
  where \eqn{TRADE(t,i)} is the number of \eqn{i}-th days of the week in \eqn{t}-th data and \eqn{b(1) + \ldots + b(7) = 0}{b(1) + ... + b(7) = 0}.
}

\value{
  \item{trend}{trend component.}
  \item{seasonal}{seasonal component.}
  \item{ar}{AR process.}
  \item{trad}{trading day factor.}
  \item{noise}{observational noise.}
  \item{aic}{AIC.}
  \item{lkhd}{likelihood.}
  \item{sigma2}{sigma^2.}
  \item{tau1}{system noise variances tau2(1).}
  \item{tau2}{system noise variances tau2(2).}
  \item{tau3}{system noise variances tau2(3).}
  \item{arcoef}{vector of AR coefficients.}
  \item{tdf}{trading day factor \code{tdf(i)} (i=1,7).}
}

\references{
  G.Kitagawa (1981) \emph{A Nonstationary Time Series Model and Its Fitting by a Recursive Filter}
  Journal of Time Series Analysis, Vol.2, 103-116.

  W.Gersch and G.Kitagawa (1983) \emph{The prediction of time series with Trends and Seasonalities}
  Journal of Business and Economic Statistics, Vol.1, 253-264.

  G.Kitagawa (1984)
  \emph{A smoothness priors-state space modeling of Time Series with Trend and Seasonality}
  Journal of American Statistical Association, VOL.79, NO.386, 378-389.
}

\examples{
  data(Blsallfood)
  z <- decomp(Blsallfood, trade=TRUE, year=1973)
  z$aic
  z$lkhd
  z$sigma2
  z$tau1
  z$tau2
  z$tau3
}

\keyword{ts}
