\name{mulmar}
\alias{mulmar}
\title{Multivariate Case of Minimum AIC Method of AR Model Fitting}
\description{
  Fit a multivariate autoregressive model by the minimum AIC procedure.
  Only the possibilities of zero coefficients at the beginning and end of the model are considered.
  The least squares estimates of the parameters are obtained by the householder transformation.
}
\usage{
mulmar(y, max.order=NULL, plot=FALSE)
}
\arguments{
  \item{y}{a multivariate time series.}
  \item{max.order}{upper limit of the order of AR model. Default is \eqn{2 \sqrt{n}}{2*sqrt(n)}, where \eqn{n} is the length of the time series \code{y}.}
  \item{plot}{logical. If TRUE \code{daic[[1]],...,daic[[d]]} are plotted, where \eqn{d} is the dimension of the multivariate time series.}
}
\value{
  \item{mean}{mean.}
  \item{var}{variance.}
  \item{v}{innovation variance.}
  \item{aic}{AIC.}
  \item{aicmin}{minimum AIC.}
  \item{daic}{AIC-aicmin.}
  \item{order.maice}{order of minimum AIC.}
  \item{v.maice}{MAICE innovation variance.}
  \item{np}{number of parameters.}
  \item{jnd}{specification of \eqn{i}-th regressor.}
  \item{subregcoef}{subset regression coefficients.}
  \item{rvar}{residual variance.}
  \item{aicf}{final estimate of AIC (\eqn{=n\log}{=nlog}(\code{rvar})\eqn{+2}\code{np}).}
  \item{respns}{instantaneous response.}
  \item{matv}{innovation variance matrix.}
  \item{morder}{order of the MAICE model.}
  \item{arcoef}{AR coefficients. \code{arcoef[i,j,k]} shows the value of \eqn{i}-th row, \eqn{j}-th column, \eqn{k-}th order.}
  \item{aicsum}{the sum of aicf.}
}
\details{
  Multivariate autoregressive model is defined by
    \deqn{y(t) = A(1)y(t-1) + A(2)y(t-2) +\ldots+ A(p)y(t-p) + u(t),}{y(t) = A(1)y(t-1) + A(2)y(t-2) +...+ A(p)y(t-p) + u(t),}
  where \eqn{p} is order of the model and \eqn{u(t)} is Gaussian white noise with mean \eqn{0} and variance matrix \code{matv}.
  AIC is defined by
    \deqn{AIC = n \log(det(v)) + 2k,}{AIC = n log(det(v)) + 2k,}
  where \eqn{n} is the number of data, \eqn{v} is the estimate of innovation variance matrix,
  \eqn{det} is the determinant and \eqn{k} is the number of free parameters.
}

\references{
  G.Kitagawa and H.Akaike (1978) A Procedure for The Modeling of Non-stationary Time Series.
  Ann. Inst. Statist. Math., 30, B, 351--363.

  H.Akaike, G.Kitagawa, E.Arahata and F.Tada (1979)
  \emph{Computer Science Monograph, No.11, Timsac78.}
  The Institute of Statistical Mathematics.
}

\examples{
  # Example 1
  data(Powerplant)
  z <- mulmar(Powerplant, max.order=10)
  z$arcoef

  # Example 2
  ar <- array(0,dim=c(3,3,2))
  ar[,,1] <- matrix(c(0.4,  0,   0.3,
                      0.2, -0.1, -0.5,
                      0.3,  0.1, 0),3,3,byrow=TRUE)
  ar[,,2] <- matrix(c(0,  -0.3,  0.5,
                      0.7, -0.4,  1,
                      0,   -0.5,  0.3),3,3,byrow=TRUE)
  x <- matrix(rnorm(200*3),200,3)
  y <- mfilter(x,ar,"recursive")
  z <- mulmar(y, max.order=10)
  z$arcoef
}

\keyword{ts}
