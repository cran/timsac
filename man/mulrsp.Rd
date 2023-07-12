\name{mulrsp}
\alias{mulrsp}
\title{Multiple Rational Spectrum}
\description{
  Compute rational spectrum for d-dimensional ARMA process.
}
\usage{
mulrsp(h, d, cov, ar = NULL, ma = NULL, log = FALSE, plot = TRUE, ...)
}
\arguments{
  \item{h}{specify frequencies \eqn{i/2}\code{h} (\eqn{i=0,1,...,}\code{h}).}
  \item{d}{dimension of the observation vector.}
  \item{cov}{covariance matrix.}
  \item{ar}{coefficient matrix of autoregressive model. \code{ar[i,j,k]} shows
    the value of \eqn{i}-th row, \eqn{j}-th column, \eqn{k}-th order.}
  \item{ma}{coefficient matrix of moving average model. \code{ma[i,j,k]} shows
    the value of \eqn{i}-th row, \eqn{j}-th column, \eqn{k}-th order.}
  \item{log}{logical. If \code{TRUE}, rational spectrums \code{rspec} are
    plotted as \eqn{log(}\code{rspec}\eqn{)}.}
  \item{plot}{logical. If \code{TRUE}, rational spectrums \code{rspec} are
    plotted.}
  \item{\dots}{graphical arguments passed to \code{\link{plot.specmx}}.}
}
\value{
  \item{rspec}{rational spectrum. An object of class \code{"specmx"}.}
  \item{scoh}{simple coherence.}
}
\details{
  ARMA process :
  \deqn{y(t) - A(1)y(t-1) -...- A(p)y(t-p) = u(t) - B(1)u(t-1) -...- B(q)u(t-q)}
  where \eqn{u(t)} is a white noise with zero mean vector and covariance matrix
  \code{cov}.
}

\references{
  H.Akaike and T.Nakagawa (1988) \emph{Statistical Analysis and Control
    of Dynamic Systems.} Kluwer Academic publishers.
}

\examples{
# Example 1 for the normal distribution
xorg <- rnorm(1003)
x <- matrix(0, nrow = 1000, ncol = 2)
x[, 1] <- xorg[1:1000]
x[, 2] <- xorg[4:1003] + 0.5*rnorm(1000)
aaa <- ar(x)
mulrsp(h = 20, d = 2, cov = aaa$var.pred, ar = aaa$ar)

# Example 2 for the AR model
ar <- array(0, dim = c(3,3,2))
ar[, , 1] <- matrix(c(0.4,  0,   0.3,
                      0.2, -0.1, -0.5,
                      0.3,  0.1, 0), nrow = 3, ncol = 3, byrow = TRUE)
ar[, , 2] <- matrix(c(0,  -0.3,  0.5,
                      0.7, -0.4,  1,
                      0,   -0.5,  0.3), nrow = 3, ncol = 3, byrow = TRUE)
x <- matrix(rnorm(200*3), nrow = 200, ncol = 3)
y <- mfilter(x, ar, "recursive")
z <- fpec(y, max.order = 10)
mulrsp(h = 20, d = 3, cov = z$perr, ar = z$arcoef)
}

\keyword{ts}
