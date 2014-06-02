\name{perars}
\alias{perars}
\alias{print.perars}
\title{Periodic Autoregression for a Scalar Time Series}
\description{
  This is the program for the fitting of periodic autoregressive models by the method of least squares realized
  through householder transformation.
}
\usage{
  perars(y, ni, lag=NULL, ksw=0)
}
\arguments{
  \item{y}{a univariate time series.}
  \item{ni}{number of instants in one period.}
  \item{lag}{maximum lag of periods. Default is \eqn{2 \sqrt{ni}}{2*sqrt(ni)}.}
  \item{ksw}{integer. \eqn{0} denotes constant vector is not included as a regressor
    and \eqn{1} denotes constant vector is included as the first regressor.}
}
\value{
  \item{mean}{mean.}
  \item{var}{variance.}
  \item{subset}{specification of i-th regressor (i=1,...,\code{ni}).}
  \item{regcoef}{regression coefficients.}
  \item{rvar}{residual variances.}
  \item{np}{number of parameters.}
  \item{aic}{AIC.}
  \item{v}{innovation variance matrix.}
  \item{arcoef}{AR coefficient matrices. \code{arcoef[i,,k]} shows \eqn{i}-th regressand of \eqn{k}-th period former.}
  \item{const}{constant vector.}
  \item{morder}{order of the MAICE model.}
}
\details{
  Periodic autoregressive model (\eqn{i=1,\ldots,nd, j=1,\ldots,}{i=1,...,nd, j=1,...,}\code{ni}) is defined by\cr
    \eqn{z(i,j) = y(ni(i-1)+j)},\cr
    \eqn{z(i,j) = c(j) + A(1,j,0)z(i,1) + \ldots + A(j-1,j,0)z(i,j-1) + A(1,j,1)z(i-1,1) + \ldots + A(ni,j,1)z(i-1,ni) + \ldots + u(i,j)}{z(i,j) = c(j) + A(1,j,0)z(i,1) + ... + A(j-1,j,0)z(i,j-1) + A(1,j,1)z(i-1,1) + ... + A(ni,j,1)z(i-1,ni) + ... + u(i,j)},\cr
  where \eqn{nd} is the number of periods, \eqn{ni} is the number of instants in one period and \eqn{u(i,j)} is the Gaussian white noise.
  When \code{ksw} is set to \eqn{0}, the constant term \eqn{c(j)} is excluded.
     
  The statistics AIC is defined by \eqn{AIC = n \log(det(v)) + 2k}{AIC = nlog(det(v)) + 2k},
  where \eqn{n} is the length of data, \eqn{v} is the estimate of the innovation variance matrix and \eqn{k} is the number of parameters.
  The outputs are the estimates of the regression coefficients and innovation variance of the periodic AR model for each instant.
}

\references{
  M.Pagano (1978)
  On Periodic and Multiple Autoregressions.
   Ann. Statist., 6, 1310--1317.

  H.Akaike, G.Kitagawa, E.Arahata and F.Tada (1979)
  \emph{Computer Science Monograph, No.11, Timsac78.}
  The Institute of Statistical Mathematics.
}

\examples{
  data(Airpollution)
  perars(Airpollution, ni=6, lag=2, ksw=1)
}

\keyword{ts}
