\name{optdes}
\alias{optdes}
\title{Optimal Controller Design}
\description{
  Compute optimal controller gain matrix for a quadratic criterion defined by two positive definite matrices Q and R.
}
\usage{
  optdes(y, max.order=NULL, ns, q, r)
}
\arguments{
  \item{y}{a multivariate time series.}
  \item{max.order}{upper limit of model order. Default is \eqn{2 \sqrt{n}}{2*sqrt(n)}, where \eqn{n} is the length of the time series \code{y}.}
  \item{ns}{number of D.P. stages.}
  \item{q}{positive definite \eqn{(m, m)} matrix \eqn{Q}, where \eqn{m} is the number of controlled variables.
   
    A quadratic criterion is defined by \eqn{Q} and \eqn{R}.}
  \item{r}{positive definite \eqn{(l, l)} matrix \eqn{R}, where \eqn{l} is the number of manipulated variables.}

}
\value{
  \item{perr}{prediction error covariance matrix.}
  \item{trans}{first \eqn{m} columns of transition matrix, where \eqn{m} is the number of controlled variables.}
  \item{gamma}{gamma matrix.}
  \item{gain}{gain matrix.}
}

\references{
  H.Akaike and T.Nakagawa (1988) \emph{Statistical Analysis and Control
    of Dynamic Systems.} Kluwer Academic publishers.
}

\examples{
# Multivariate Example Data
  ar <- array(0,dim=c(3,3,2))
  ar[,,1] <- matrix(c(0.4,  0,   0.3,
                      0.2, -0.1, -0.5,
                      0.3,  0.1, 0),3,3,byrow=TRUE)
  ar[,,2] <- matrix(c(0,  -0.3,  0.5,
                      0.7, -0.4,  1,
                      0,   -0.5,  0.3),3,3,byrow=TRUE)
  x <- matrix(rnorm(200*3),200,3)
  y <- mfilter(x,ar,"recursive")
  q <- matrix(c(0.16,0,0,0.09), 2, 2)
  r <- matrix(0.001, 1, 1)
  optdes(y,, ns=20, q, r)
}

\keyword{ts}
