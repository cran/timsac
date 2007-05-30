\name{ngsmth}
\alias{ngsmth}
\title{Non-Gaussian Smoothing}
\description{
  Trend estimation by Non-Gaussian smoothing.
}
\usage{
  ngsmth(y,noisev=2,tau2,bv=1.0,noisew=1,sig2,bw=1.0,initd=1,k=200,plot=TRUE)
}

\arguments{
  \item{y}{a univariate time series.}
  \item{noisev}{type of system noise density. (1:Gaussian (normal), 2:Pearson family, 3:two-sides exponential)}
  \item{tau2}{variance of dispersion of system noise.}
  \item{bv}{shape parameter of system noise (for noisev=2).}
  \item{noisew}{type of observation noise density. (1:Gaussian (normal), 2:Pearson family, 3:two-sided exponential, 4:double exponential)}
  \item{sig2}{variance of dispersion of observation noise.}
  \item{bw}{shape parameter of observation noise (for noisew=2).}
  \item{initd}{type of density function. (0:two-sided exponential, 1:Gaussian (normal), 2:uniform)}
  \item{k}{number of intervals}
  \item{plot}{logical. If TRUE (default) trend and smoothed density are plotted.}
}

\details{
  Consider a one dimensional state space model

    \eqn{x(n) = x(n-1) + v(n)},

    \eqn{y(n) = x(n) + w(n)},

  where the observation noise \eqn{w(n)} is assumed to be Gaussian distributed and
  the system noise \eqn{v(n)} is assumed to be distributed as the Pearson system

    \eqn{q(v(n)) = c/(\tau^2 + v(n)^2)^b}  

  with \eqn{1/2 < b < \infty} and \eqn{c = \tau^{2b-1}\Gamma(b) / \Gamma(1/2)\Gamma(b - 1/2)}.

  This broad family of distributions includes the Cauchy distribution (b=1).
}

\value{
  \item{trend}{trend.}
  \item{smt}{smoothed density.}
  \item{lkhood}{log-likelihood.}
}

\references{
  Kitagawa, G., (1993) \emph{Time series analysis programing (in Japanese)}. The Iwanami Computer Science Senes.

  Kitagawa, G. and Gersch, W., (1996) \emph{Smoothness Priors Analysis of Time Series}. Lecture Notes in Statistics, No.116, Springer-Verlag.
}

\examples{

  # trend model
  x <- rep(0,400)
  x[101:200] <- 1
  x[201:300] <- -1
  y <- rnorm(600, mean=0, sd=0.5)
  y <- y[101:500]
  z1 <- ngsmth(x+y, noisev=2, tau2=0.211e-09, bv=1, noisew=2, sig2=1.042, bw=1 )

  z2 <- ngsmth(x+y, noisev=1, tau2=0.14e-01, bv=1, noisew=2, sig2=1.048, bw=1 )


  #  an earthquake wave data
  data(MYE1F)
  n <- length(MYE1F)
  m <- n/2
  y <- rep(0, n)
  for( i in 2:n )  y[i] <- MYE1F[i] - 0.5*MYE1F[i-1]
  yy <- rep(0, m)
  for( i in 1:m ) yy[i] <- y[i*2]
  z <- tvvar(yy, trend.order=2, tau20= 6.6e-06, delta=1.0e-06, plot=FALSE)
  z1 <- ngsmth( z$ts, noisev=2, tau2=0.00026, bv=1, noisew=2, sig2=1.644934, bw=1, k=190 )

  z2 <- ngsmth( z$ts, noisev=1, tau2=0.04909, bv=1, noisew=2, sig2=1.644934, bw=1, k=190 )

}

\keyword{ts}
