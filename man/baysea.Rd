\name{baysea}
\alias{baysea}
\title{Bayesian Seasonal Adjustment Procedure}
\description{
    Decompose a nonstationary time series into several possible components.
}

\usage{
  baysea(y, period=12, span=4, shift=1, forecast=0, trend.order=2, seasonal.order=1, year=0, month=1, out=0,
         rigid=1, zersum=1, delta=7, alpha=0.01, beta=0.01, gamma=0.1, spec=TRUE, plot=TRUE, separate.graphics=FALSE)
}

\arguments{
  \item{y}{a univariate time series.}
  \item{period}{number of seasonals within a period.}
  \item{span}{number of periods to be processed at one time.}
  \item{shift}{number of periods to be shifted to define the new span of data.}
  \item{forecast}{length of forecast at at the end of data.}
  \item{trend.order}{order of differencing of trend.}
  \item{seasonal.order}{order of differencing of seasonal. seasonal.order is smaller than or equal to span.}
  \item{year}{trading-day adjustment option. 

    \eqn{=0} : without trading-day adjustment

    \eqn{>0} : with trading-day adjustment (the series is supposed to start at this "year") }
  \item{month}{number of the month in which the series starts. If \eqn{year=0} this parameter is ignored.}
  \item{out}{outlier correction option.

    \eqn{=0} : without outlier detection

    \eqn{=1} : with outlier detection by marginal probability

    \eqn{=2} : with outlier detection by model selection }
  \item{rigid}{controls the rigidity of the seasonal component. more rigid seasonal with larger than rigid.}
  \item{zersum}{controls the sum of the seasonals within a period.}
  \item{delta}{controls the leap year effect.}
  \item{alpha}{controls prior variance of initial trend.}
  \item{beta}{controls prior variance of initial seasonal.}
  \item{gamma}{controls prior variance of initial sum of seasonal.}
  \item{spec}{logical. If TRUE (default) estimate spectra of irregular and differenced adjusted.}
  \item{plot}{logical. If TRUE (default) plot trend, adjust, smoothed, season and irregular.}
  \item{separate.graphics}{logical. If TRUE a graphic device is opened for each graphics display.}
}

\details{
  This function realized a decomposition of time series Y into the form

    \eqn{y(t) = T(t) + S(t) + I(t) + TDC(t) + OCF(t)}

  where \eqn{T(t)} is trend component, \eqn{S(t)} is seasonal component, \eqn{I(t)} is irregular, \eqn{TDC(t)} is trading day factor and \eqn{OCF(t)} is outlier correction factor.

  For the purpose of comparison of models the criterion ABIC is defined

    \eqn{ABIC = -2(log maximum likelihood of the model)}

  Smaller value of ABIC represents better fit.
}

\value{
  \item{outlier}{outlier correction factor.}
  \item{trend}{trend.}
  \item{season}{seasonal.}
  \item{tday}{trading-day component if \eqn{year>0}.}
  \item{irregular}{irregular = data - trend - season - tday - ootlier.}
  \item{adjust}{adjusted = trend - irreguar.}
  \item{smoothed}{smoothed = trend + season + tday.}
  \item{aveABIC}{averaged ABIC.}
  \item{irregular.spec}{a list of acov(autocovariances), acor(normalized covariances), mean, v(innovation variance), aic(AIC), parcor(partial autocorrelation)
                        and rspec(rational spectrum) of irregular if \eqn{spec=TRUE}.}
  \item{adjusted.spec}{a list of acov(autocovariances), acor(normalized covariances), mean, v(innovation variance), aic(AIC), parcor(partial autocorrelation)
                       and rspec(rational spectrum) of differenced adjusted series if \eqn{spec=TRUE}.}
  \item{differenced.trend}{a list of acov(autocovariances), acor(normalized covariances), mean, v(innovation variance), aic(AIC) and parcor(partial autocorrelation)
                           of differenced trend series if \eqn{spec=TRUE}.}
  \item{differenced.season}{a list of acov(autocovariances), acor(normalized covariances), mean, v(innovation variance), aic(AIC) and parcor(partial autocorrelation)
                             of differenced seasonal series if \eqn{spec=TRUE}.}
}

\references{
  H.Akaike, T.Ozaki, M.Ishiguro, Y.Ogata, G.Kitagawa, Y-H.Tamura, E.Arahata, K.Katsura and Y.Tamura (1985) \emph{Computer Science Monograph, No.22,
   Timsac84 Part 1}.
  The Institute of Statistical Mathematics.
}

\examples{
  data(LaborData)
  baysea(y=LaborData, forecast=12)
}

\keyword{ts}
