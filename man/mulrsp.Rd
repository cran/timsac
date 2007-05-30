\name{mulrsp}
\alias{mulrsp}
\title{Multiple Rational Spectrum}
\description{
  Compute rational spectrum for d-dimensional ARMA process.
}
\usage{
mulrsp(h, d, cov, ar=NULL, ma=NULL, log=FALSE, plot=TRUE, plot.scale=FALSE)
}
\arguments{
  \item{h}{specify frequencies i/2h (i=0,1,...,h).}
  \item{d}{dimension of the observation vector.}
  \item{cov}{covariance matrix.}
  \item{ar}{coefficient matrix of autoregressive model. ar[i,j,k] shows the value of i-th row, j-th column, k-th order.}
  \item{ma}{coefficient matrix of moving average model. ma[i,j,k] shows the value of i-th row, j-th column, k-th order.}
  \item{log}{logical. If TRUE rational spectrums rspec are plotted as log(rspec).}
  \item{plot}{logical. If TRUE rational spectrums rspec are plotted.}
  \item{plot.scale}{logical. IF TRUE the common range of the y-axisis is used.}
}
\value{
  \item{rspec}{rational spectrum.}
  \item{scoh}{simple coherence.}
}
\details{
  ARMA process :

  \eqn{y(t) - A(1)y(t-1) -...- A(p)y(t-p) = u(t) - B(1)u(t-1) -...- B(q)u(t-q)}

  where u(t) is a white noise with zero mean vector and covariance matrix \code{cov}.
}

\references{
  H.Akaike and T.Nakagawa (1988) \emph{Statistical Analysis and Control
    of Dynamic Systems.} Kluwer Academic publishers.
}

\examples{
  # Example 1 for the normal distribution
  xorg <- rnorm(1003)
  x <- matrix(0,1000,2)
  x[,1] <- xorg[1:1000]
  x[,2] <- xorg[4:1003]+0.5*rnorm(1000)
  aaa <- ar(x)
  mulrsp(h=20, d=2, cov=aaa$var.pred, ar=aaa$ar, plot=TRUE, plot.scale=TRUE)

  # Example 2 for the AR model
  ar <- array(0,dim=c(3,3,2))
  ar[,,1] <- matrix(c(0.4,  0,   0.3,
                      0.2, -0.1, -0.5,
                      0.3,  0.1, 0),3,3,byrow=TRUE)
  ar[,,2] <- matrix(c(0,  -0.3,  0.5,
                      0.7, -0.4,  1,
                      0,   -0.5,  0.3),3,3,byrow=TRUE)
  x <- matrix(rnorm(200*3),200,3)
  y <- mfilter(x,ar,"recursive")
  z <- fpec(y, max.order=10, ncon=3, nman=0)
  mulrsp(h=20, d=3, cov=z$perr, ar=z$arcoef)   # d=ncon+nman
}

\keyword{ts}
