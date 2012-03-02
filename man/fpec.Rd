\name{fpec}
\alias{fpec}
\alias{print.fpec}
\title{AR model Fitting for Control}
\description{
  Perform AR model fitting for control.
}
\usage{
  fpec(y, max.order=NULL, control=NULL, manip=NULL)
}
\arguments{
  \item{y}{a multivariate time series.}
  \item{max.order}{upper limit of model order. Default is \eqn{2 \sqrt{n}}{2*sqrt(n)}, where \eqn{n} is the length of time series \code{y}.}
  \item{control}{controlled variables. Default is \eqn{c(1:d)}, where \eqn{d} is the dimension of the time series \code{y}.}
  \item{manip}{manipulated variables. Default number of manipulated variable is \eqn{0}.}
}
\value{
  \item{cov}{covariance matrix rearrangement by \code{inw}.}
  \item{fpec}{FPEC (AR model fitting for control).}
  \item{rfpec}{RFPEC.}
  \item{aic}{AIC.}
  \item{ordermin}{order of minimum FPEC.}
  \item{fpecmin}{minimum FPEC.}
  \item{rfpecmin}{minimum RFPEC.}
  \item{aicmin}{minimum AIC.}
  \item{perr}{prediction error covariance matrix.}
  \item{arcoef}{a set of coefficient matrices. \code{arcoef[i,j,k]} shows the value of \eqn{i}-th row, \eqn{j}-th column, \eqn{k}-th order.}
}

\references{
  H.Akaike and T.Nakagawa (1988) \emph{Statistical Analysis and Control
    of Dynamic Systems.} Kluwer Academic publishers.
}

\examples{
  ar <- array(0,dim=c(3,3,2))
  ar[,,1] <- matrix(c(0.4,  0,   0.3,
                      0.2, -0.1, -0.5,
                      0.3,  0.1, 0),3,3,byrow=TRUE)
  ar[,,2] <- matrix(c(0,  -0.3,  0.5,
                      0.7, -0.4,  1,
                      0,   -0.5,  0.3),3,3,byrow=TRUE)
  x <- matrix(rnorm(200*3),200,3)
  y <- mfilter(x,ar,"recursive")
  fpec(y, max.order=10)
}

\keyword{ts}
