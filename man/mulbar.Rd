\name{mulbar}
\alias{mulbar}
\title{Multivariate Bayesian Method of AR Model Fitting}
\description{
  Determine multivariate autoregressive models by a Bayesian procedure.
  The basic least squares estimates of the parameters are obtained by the
  householder transformation.
}
\usage{
  mulbar(y, max.order = NULL, plot = FALSE)
}
\arguments{
  \item{y}{a multivariate time series.}
  \item{max.order}{upper limit of the order of AR model, less than or equal to
    \eqn{n/2d} where \eqn{n} is the length and \eqn{d} is the dimension of the
    time series \code{y}.
 Default is \eqn{min(2 \sqrt{n}, n/2d)}{min(2*sqrt(n), n/2d)}.}
  \item{plot}{logical. If \code{TRUE}, \code{daic} is plotted.}

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
  \item{bweight}{Bayesian weights.}
  \item{integra.bweight}{integrated Bayesian Weights.}
  \item{arcoef.for}{AR coefficients (forward model). \code{arcoef.for[i,j,k]}
    shows the value of \eqn{i}-th row, \eqn{j}-th column, \eqn{k}-th order.}
  \item{arcoef.back}{AR coefficients (backward model). \code{arcoef.back[i,j,k]}
    shows the value of \eqn{i}-th row, \eqn{j}-th column, \eqn{k}-th order.}
  \item{pacoef.for}{partial autoregression coefficients (forward model).}
  \item{pacoef.back}{partial autoregression coefficients (backward model).}
  \item{v.bay}{innovation variance of the Bayesian model.}
  \item{aic.bay}{equivalent AIC of the Bayesian (forward) model.}
}
\details{
  The statistic AIC is defined by
    \deqn{AIC = n \log(det(v)) + 2k,}{AIC = n log(det(v)) + 2k,}
  where \eqn{n} is the number of data, \eqn{v} is the estimate of innovation
  variance matrix, \eqn{det} is the determinant and \eqn{k} is the number of
  free parameters.

  Bayesian weight of the \eqn{m}-th order model is defined by
    \deqn{W(n) = const \times \frac{C(m)}{m+1},}{%
          W(n) = const * C(m) / (m+1),}
  where \eqn{const} is the normalizing constant and
  \eqn{C(m)=\exp(-0.5 AIC(m))}{C(m)=exp(-0.5 AIC(m))}. The Bayesian estimates of
  partial autoregression coefficient matrices of forward and backward models are
  obtained by \eqn{(m = 1,\ldots,lag)}{(m = 1,...,lag)}
    \deqn{G(m) = G(m) D(m),}
    \deqn{H(m) = H(m) D(m),}
  where the original \eqn{G(m)} and \eqn{H(m)} are the (conditional) maximum
  likelihood estimates of the highest order coefficient matrices of forward and
   backward AR models of order \eqn{m} and \eqn{D(m)} is defined by
    \deqn{D(m) = W(m) + \ldots + W(lag).}{D(m) = W(m) + ... + W(lag).}
  The equivalent number of parameters for the Bayesian model is defined by
    \deqn{ek = \{ D(1)^2 + \ldots + D(lag)^2 \} id + \frac{id(id+1)}{2}}{%
          ek = {D(1)^2 + ... + D(lag)^2} * id + id * (id+1)/2}
  where \eqn{id} denotes dimension of the process.
}

\references{
  H.Akaike (1978)
  A Bayesian Extension of The Minimum AIC Procedure of Autoregressive Model
  Fitting. Research Memo. NO.126, The Institute of Statistical Mathematics.

  G.Kitagawa and H.Akaike (1978)
  A Procedure for The Modeling of Non-stationary Time Series.
   Ann. Inst. Statist. Math., 30, B, 351--363.

  H.Akaike, G.Kitagawa, E.Arahata and F.Tada (1979)
  \emph{Computer Science Monograph, No.11, Timsac78.}
   The Institute of Statistical Mathematics.
}

\examples{
data(Powerplant)
z <- mulbar(Powerplant, max.order = 10)
z$pacoef.for
z$pacoef.back
}

\keyword{ts}
