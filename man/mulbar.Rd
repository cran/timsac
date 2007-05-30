\name{mulbar}
\alias{mulbar}
\title{Multivariate Bayesian Method of AR Model Fitting}
\description{
  Determine multivariate autoregressive models by a Bayesian procedure.
  The basic least squares estimates of the parameters are obtained by the householder transformation.
}
\usage{
  mulbar(y, max.order=NULL, plot=FALSE)
}
\arguments{
  \item{y}{a multivariate time series.}
  \item{max.order}{upper limit of the order of AR model. Default is 2*sqrt(n), where n is the length of the time series y.}
  \item{plot}{logical. If TRUE daic is plotted.}

}
\value{
  \item{mean}{mean.}
  \item{var}{variance.}
  \item{v}{innovation variance.}
  \item{aic}{AIC(m), (m = 0,...,max.order).}
  \item{aicmin}{minimum AIC.}
  \item{daic}{AIC(m)-aicmin (m = 0,...,max.order).}
  \item{order.maice}{order of minimum AIC.}
  \item{v.maice}{innovation variance attained at m=order.maice.}
  \item{bweight}{Bayesian weights.}
  \item{integra.bweight}{integrated Bayesian Weights.}
  \item{arcoef.for}{AR coefficients (forward model). arcoef.for[i,j,k] shows the value of i-th row, j-th column, k-th order.}
  \item{arcoef.back}{AR coefficients (backward model). arcoef.back[i,j,k] shows the value of i-th row, j-th column, k-th order.}
  \item{pacoef.for}{partial autoregression coefficients (forward model).}
  \item{pacoef.back}{partial autoregression coefficients (backward model).}
  \item{v.bay}{innovation variance of the Bayesian model.}
  \item{aic.bay}{equivalent AIC of the Bayesian (forward) model.}
}
\details{
  The statistic AIC is defined by

    \eqn{AIC = nlog(det(v)) + 2k},

  where n is the number of data, v is the estimate of innovation variance matrix,
  det is the determinant and k is the number of free parameters.

  Bayesian weight of the m-th order model is defined by

    \eqn{W(n) = const * C(m) / (m+1)},

  where const is the normalizing constant and \eqn{C(m)=exp(-0.5AIC(m))}.

  The Bayesian estimates of partial autoregression coefficient matrices of forward and backward models are
  obtained by (m = 1,...,lag)

    \eqn{G(m) = G(m) D(m)},

    \eqn{H(m) = H(m) D(m)},

  where the original G(m) and H(m) are the (conditional) maximum likelihood estimates
  of the highest order coefficient matrices of forward and backward AR models of order m
  and D(m) is defined by

    \eqn{D(m) = W(m) + ... + W(lag)}.

  The equivalent number of parameters for the Bayesian model is defined by

    \eqn{ek = (D(1)^2 + ... + D(lag)^2) id + id (id+1)/2}

  where id denotes dimension of the process.
}

\references{
  H.Akaike (1978)
  A Bayesian Extension of The Minimum AIC Procedure of Autoregressive Model Fitting.
   Research Memo. NO.126, The Institute of Statistical Mathematics.

  G.Kiagawa and H.Akaike (1978)
  A Procedure for The Modeling of Non-stationary Time Series.
   Ann. Inst. Statist. Math., 30, B, 351--363.

  H.Akaike, G.Kitagawa, E.Arahata and F.Tada (1979)
  \emph{Computer Science Monograph, No.11, Timsac78.}
   The Institute of Statistical Mathematics.
}

\examples{
  data(Powerplant)
  z <- mulbar(Powerplant, max.order=10)
  z$pacoef.for
  z$pacoef.back
}

\keyword{ts}
