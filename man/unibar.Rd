\name{unibar}
\alias{unibar}
\title{Univariate Bayesian Method of AR Model Fitting}
\description{
  This program fits an autoregressive model by a Bayesian procedure.
  The least squares estimates of the parameters are obtained by
  the householder transformation.
}
\usage{
  unibar(y, ar.order=NULL, plot=TRUE)
}
\arguments{
  \item{y}{a univariate time series.}
  \item{ar.order}{order of the AR model. Default is 2*sqrt(n), where n is the length of the time series y.}
  \item{plot}{logical. If TRUE (default) daic, pacoef and pspec are plotted.}
}
\value{
  \item{mean}{mean.}
  \item{var}{variance.}
  \item{v}{innovation variance.}
  \item{aic}{AIC(m) (m = 0,...,ar.order).}
  \item{aicmin}{minimum AIC.}
  \item{daic}{AIC(m)-aicmin (m = 0,...,ar.order).}
  \item{order.maice}{order of minimum AIC.}
  \item{v.maice}{innovation variance attained at m=order.maice.}
  \item{pacoef}{partial autocorrelation coefficients (least squares estimate).}
  \item{bweight}{Bayesian Weight.}
  \item{integra.bweight}{integrated Bayesian weights.}
  \item{v.bay}{innovation variance of Bayesian model.}
  \item{aic.bay}{AIC of Bayesian model.}
  \item{np}{equivalent number of parameters.}
  \item{pacoef.bay}{partial autocorrelation coefficients of Bayesian model.}
  \item{arcoef}{AR coefficients of Bayesian model.}
  \item{pspec}{power spectrum.}
}
\details{
  The AR model is given by

    \eqn{y(t) = a(1)y(t-1) + .... + a(p)y(t-p) + u(t)}

  where p is AR order and u(t) is Gaussian white noise with mean 0 and variance v(p).

  The basic statistic AIC is defined by

    \eqn{AIC = nlog(det(v)) + 2m},

  where n is the length of data, v is the estimate of innovation variance, and m is the order of the model.

  Bayesian weight of the m-th order model is defined by

    \eqn{W(m) = CONST*C(m) / (m+1)},

  where CONST is the normalizing constant and \eqn{C(m)=exp(-0.5AIC(m))}.

  The equivalent number of free parameter for the Bayesian model is defined by

    \eqn{ek = D(1)^2 +...+ D(k)^2 +1}

  where D(j) is defined by \eqn{D(j)=W(j) +...+ W(k)}.

  m in the definition of AIC is replaced by ek to be define an equivalent AIC for a Bayesian model.
}

\references{
  H.Akaike (1978) A Bayesian Extension of The Minimum AIC Procedure of Autoregressivemodel Fitting.
  Research memo. No.126. The Institute of Statistical Mathematics.

  G.Kitagawa and H.Akaike (1978) A Procedure for The Modeling of Non-Stationary Time Series.
  Ann. Inst. Statist. Math., 30, B, 351--363.

  H.Akaike, G.Kitagawa, E.Arahata and F.Tada (1979)
  \emph{Computer Science Monograph, No.11, Timsac78.}
  The Institute of Statistical Mathematics.
}

\examples{
  data(Canadianlynx)
  z <- unibar(Canadianlynx, ar.order=20)
  z$arcoef
}

\keyword{ts}
