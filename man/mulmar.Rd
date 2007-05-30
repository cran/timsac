\name{mulmar}
\alias{mulmar}
\title{Multivariate Case of Minimum AIC Method of AR Model Fitting}
\description{
  Fit a multivariate autoregressive model by the minimum AIC procedure.
  Only the possiblities of zero coefficients at the begining and end of the model are considered.
  The least squares estimates of the parameters are obtained by the householder transformation.
}
\usage{
mulmar(y, max.order=NULL, plot=FALSE, tmp.file=NULL)
}
\arguments{
  \item{y}{a multivariate time series.}
  \item{max.order}{upper limit of the order of AR model. Default is 2*sqrt(n), where n is the length of the time series y.}
  \item{plot}{logical. If TRUE daic[[1]],...,daic[[d]] are plotted, where d is the dimention of the multivariate time series.}
  \item{tmp.file}{a character string naming a file written intermediate results of AIC minimazation.
        If NULL (default) output no file.}
}
\value{
  \item{mean}{mean.}
  \item{var}{variance.}
  \item{v}{innovation variance.}
  \item{aic}{AIC(m) (m = 0,...,max.order).}
  \item{aicmin}{minimum AIC.}
  \item{daic}{AIC(m)-aicmin (m = 0,...,max.order).}
  \item{order.maice}{order of minimum AIC.}
  \item{v.maice}{innovation variance attained at m = order.maice.}
  \item{np}{number of parameters.}
  \item{jnd}{specification of i-th regressor.}
  \item{subregcoef}{subset regression coefficients.}
  \item{rvar}{residual variance.}
  \item{aicf}{final estimate of AIC (=nlog(rvar)+2np).}
  \item{respns}{instantaneous response.}
  \item{matv}{innovation variance matrix.}
  \item{morder}{order of the MAICE model.}
  \item{arcoef}{AR coefficients. arcoef[i,j,k] shows the value of i-th row, j-th column, k-th order.}
  \item{aicsum}{the sum of aicf.}
}
\details{
  Multivariate autoregressive model is defined by

    \eqn{y(t) = A(1)y(t-1) + A(2)y(t-2) +...+ A(p)y(t-p) + u(t)},

  where p is order of the model and u(t) is Gaussian white noise with mean 0 and variance matrix matv.

  AIC is defined by

    \eqn{AIC = nlog(det(v)) + 2k}

  where n is the number of data, v is the estimate of innovation variance matrix,
  det is the determinant and k is the number of free parameters.
}

\references{
  G.Kitagawa and H.Akaike (1978) A Procedure for The Modeling of Non-stationary Time Series.
  Ann. Inst. Statist. Math., 30, B, 351--363.

  H.Akaike, G.Kiragawa, E.Arahata and F.Tada (1979)
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
