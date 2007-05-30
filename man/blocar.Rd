\name{blocar}
\alias{blocar}
\title{Bayesian Method of Locally Stationary AR Model Fitting; Scalar Case}
\description{
  Locally fit autoregressive models to non-stationary time series by a Bayesian procedure.
}
\usage{
blocar(y, max.order=NULL, span, plot=TRUE)
}
\arguments{
  \item{y}{a univariate time series.}
  \item{max.order}{upper limit of the order of AR model. Default is 2*sqrt(n), where n is the length of the time series y.}
  \item{span}{length of basic local span.}
  \item{plot}{logical. If TRUE (default) spectrums pspec are plotted.}
}
\value{
  \item{var}{variance.}
  \item{aic}{AIC.}
  \item{bweight}{Bayesian weight.}
  \item{pacoef}{partial autocorrelation.}
  \item{arcoef}{coefficients ( average by the Bayesian weights ).}
  \item{v}{innovation variance.}
  \item{init}{initial point of the data fitted to the current model.}
  \item{end}{end point of the data fitted to the current model.}
  \item{pspec}{power spectrum.}
}

\details{
  The basic AR model of scalar time series y(t) (t=1,...,n) is given by

    \eqn{y(t) = a(1)y(t-1) + a(2)y(t-2) +...+ a(p)y(t-p) + u(t)},

  where p is order of the model and u(t) is Gaussian white noise with mean 0 and variance v.

  At each stage of modeling of locally AR model, a two-step Bayesian procedure is applied

    1.  Averaging of the models with different orders fitted to the newly obtained data.

    2.  Averaging of the models fitted to the present and preceding spans.

  AIC of the model fitted to the new span is defined by

  \eqn{AIC = ns log( det(v) ) + 2k},

  where ns is the length of new data, v is the innovation variance
  and k is the equivalent number of parameters, defined as the sum of squares of the Bayesian weights.

  AIC of the model fitted to the preceding spans are defined by

    \eqn{AIC( j+1 ) = ns log( det(v(j) ) + 2}

  where v(j) is the prediction error variance by the model fitted to j periods former span.
}

\references{
  G.Kitagawa and H.Akaike (1978)
  A Procedure for The Modeling of Non-Stationary Time Series.
  Ann. Inst. Statist. Math., 30, B, 351--363.

  H.Akaike (1978)
  A Bayesian Extension of the Minimin MIC Procedure of Autoregressive Model Fitting.
  Reseach Meno. NO.126. The Institute of The Statistical Mathematics.

  H.Akaike, G.Kitagawa, E.Arahata and F.Tada (1979)
  \emph{Computer Science Monograph, No.11, Timsac78.}
  The Institute of Statistical Mathematics.
}

\examples{
  data(locarData)
  z <- blocar(locarData, max.order=10, span=300)
  z$arcoef
}

\keyword{ts}
