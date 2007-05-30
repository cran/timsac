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
  \item{lag}{maximum lag of periods. Default is 2*sqrt(ni).}
  \item{ksw}{integer. 0 denotes constant vector is not included as a regressor
    and 1 denotes constant vector is included as the first regressor.}
}
\value{
  \item{mean}{mean.}
  \item{var}{variance.}
  \item{ord}{specification of i-th regressor (i=1,...,ni).}
  \item{regcoef}{regression coefficients.}
  \item{rvar}{residual variances.}
  \item{np}{number of parameters.}
  \item{aic}{AIC.}
  \item{v}{innovation variance matrix.}
  \item{arcoef}{AR coefficient matrices. arcoef[i,,k] shows i-th regressand of k-th period formar.}
  \item{const}{constant vector.}
  \item{morder}{order of the MAICE model.}
}
\details{
  Periodic autoregressive model (i=1,...,nd, j=1,...,ni) is defined by

    \eqn{z(i,j) = y(ni(i-1)+j)},

    \eqn{z(i,j) = c(j) + A(1,j,0)z(i,1) + ... + A(j-1,j,0)z(i,j-1) + A(1,j,1)z(i-1,1) + ... + A(ni,j,1)z(i-1,ni) + ... + u(i,j)},

  where nd is the number of periods, ni is the number of instants in one period and u(i,j) is the Gaussian white noise.
  When ksw is set to 0, the constant term c(j) is excluded.
     
  The statistics AIC is defined by

    \eqn{AIC = nlog(det(v)) + 2k},

  where n is the length of data, v is the estimate of the innovation variance matrix and k is the number of parameters.

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
  data(Airpolution)
  z <- perars(Airpolution, ni=6, lag=2, ksw=1)
  z$regcoef
  z$v
}

\keyword{ts}
