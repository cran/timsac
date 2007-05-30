\name{tsmooth}
\alias{tsmooth}
\title{Kalman Filter}
\description{
  State estimation of user-defined state space model by Kalman filter.
}
\usage{
  tsmooth(y,f,g,h,q,r,x0=NULL,v0=NULL,filter.end=NULL,predict.end=NULL,outmin=-10.0e+30,outmax=10.0e+30,missed=NULL,np=NULL,plot=FALSE)
}

\arguments{
  \item{y}{\eqn{l}-dimensional time series y(n).}
  \item{f}{\eqn{m * m} state transition matrix F(n), where \eqn{m} is the dimension of the state vector x(n).}
  \item{g}{\eqn{m * k} matrix G(n), where \eqn{k} is the demension of the system noise.}
  \item{h}{\eqn{l * m} matrix H(n).}
  \item{q}{system noise variance Q(n).}
  \item{r}{observational noise variance R(n).}
  \item{x0}{initial state vector X(0|0).}
  \item{v0}{initial state covariance matrix V(0|0).}
  \item{filter.end}{end point of filtering.}
  \item{predict.end}{end point of prediction.}
  \item{outmin}{lower limits of observations.}
  \item{outmax}{upper limits of observations.}
  \item{missed}{start position of missed intervals.}
  \item{np}{number of missed observations.}
  \item{plot}{logical. If TRUE estimated smoothed state is plotted.}
}

\details{
  The linear Gaussian state space model is

    \eqn{x(n) = F(n)x(n-1) + G(n)v(n)},

    \eqn{y(n) = H(n)x(n) + w(n)},

  where \eqn{y(n)} is an \eqn{l}-dimensional time series, \eqn{x(n)} is \eqn{m}-dimensional state vector.
  \eqn{v(n)} and \eqn{w(n)} are \eqn{k}- and \eqn{l}-dimensional white noise sequences. 
  \eqn{F(n)}, \eqn{G(n)} and \eqn{H(n)} are \eqn{m * m}, \eqn{m * k} and \eqn{l * m} matrices, respectively.
  \eqn{R(n)} and \eqn{Q(n)} are \eqn{k * k} and \eqn{l * l} matrices, respectively.
  We asusume that \eqn{E(v(n),w(n)) = 0}, \eqn{v(n) ~ N(0,Q(n))} and \eqn{w(n) ~ N(0,R(n))}.
  User should give all the matrices of a state space model and its parameters.
  In current version, F(n), G(n), H(n), Q(n), R(n) should be time invariant.
}

\value{
  \item{mean.smooth}{mean vectors of the smoother.}
  \item{cov.smooth}{variance of the smoother.}
  \item{esterr}{estimation error.}
  \item{lkhood}{log-likelihood.}
  \item{aic}{AIC.}
}

\references{
  Kitagawa, G., (1993) \emph{Time series analysis programing (in Japanese)}. The Iwanami Computer Science Senes.

  Kitagawa, G. and Gersch, W., (1996) \emph{Smoothness Priors Analysis of Time Series}. Lecture Notes in Statistics, No.116, Springer-Verlag.
}

\examples{
## AR model (l=1, m=10, k=1)
#  m <- 5
  m <- 10
  k <- 1
  data(Blsallfood)
  z1 <- exsar(Blsallfood, max.order=m)
  var <- z1$var
  tau2 <- z1$v.mle
  arcoef <- z1$arcoef.mle
  f <- matrix( 0.0e0, m, m )
  f[1,] <- arcoef
  for( i in 2:m ) f[i,i-1] <- 1
  g <- c(1, rep(0.0e0, m-1))
  h <- c(1, rep(0.0e0, m-1))
  q <- tau2
  r <- 0.0e0
  x0 <- rep(0.0e0, m)
  v0 <- matrix( 0.0e0, m, m )
  for( i in 1:m ) v0[i,i] <- var
  z <- tsmooth(Blsallfood, f, g, h, q, r, x0, v0, filter.end=156, predict.end=170, missed=c(41,101), np=c(30,20))

# plot mean vector and estimation error
  xss <- z$mean.smooth[1,] + mean(Blsallfood)
  cov <- z$cov.smooth
  c1 <- xss + sqrt(cov[1,])
  c2 <- xss - sqrt(cov[1,])
  err <- z$esterr
  par(mfcol=c(2,1))
  ymax <- as.integer(max(xss,c1,c2)+1)
  ymin <- as.integer(min(xss,c1,c2)-1)
  plot(c1, type='l', ylim=c(ymin,ymax), col=2, xlab="Mean vectors of the smoother XSS(1,) +/- standard deviation", ylab="")
  par(new=TRUE)
  plot(c2, type='l', ylim=c(ymin,ymax), col=3, xlab="", ylab="")
  par(new=TRUE)
  plot(xss, type='l', ylim=c(ymin,ymax), xlab="", ylab="")
  plot(err[,1,1], type='h', xlim=c(1,length(xss)),xlab="estimation error", ylab="")

## Trend model (l=3, m=2, k=2)
  l <- 3
#  n <- 400
  n <- 500
  m <- 2
  k <- 2
  f <- matrix( c(1, 0,   0, 1), nrow=m, ncol=m, byrow=TRUE )
  g <- matrix( c(1, 0,   0, 1), nrow=m, ncol=k, byrow=TRUE )
  h <- matrix( c(0.1, -0.1,   -0.05, 0.05,   0.2, 0.15), nrow=l, ncol=m, byrow=TRUE)
  q <- matrix( c(0.2*0.2, 0,   0, 0.3*0.3), nrow=k, ncol=k, byrow=TRUE)
  r <- matrix( c(0.2*0.2, 0, 0,   0, 0.1*0.1, 0,   0, 0, 0.15*0.15), nrow=l, ncol=l, byrow=TRUE)
  Xn <- matrix( 0, nrow=m, ncol=n )
  x1 <- rnorm(n+100, 0, 0.2)
  x2 <- rnorm(n+100, 0, 0.3)
  x1 <- cumsum(x1)[101:(n+100)]
  x2 <- cumsum(x2)[101:(n+100)]
  Xn[1,] <- x1-mean(x1)
  Xn[2,] <- x2-mean(x2)
  Yn <- matrix( 0, nrow=l, ncol=n )
  Wn <- matrix( 0, nrow=l, ncol=n )
  Wn[1,] <- rnorm(n, 0, 0.2)
  Wn[2,] <- rnorm(n, 0, 0.1)
  Wn[3,] <- rnorm(n, 0, 0.15)
  Yn <- h \%*\% Xn + Wn
  Yn <- aperm(Yn, c(2,1))
  x0 <- c(Xn[1,1], Xn[2,1])
  v0 <- matrix( c(var(Yn[,1]), 0,   0, var(Yn[,2])), 2, 2, byrow=TRUE)
  npe <- n+20
  z <- tsmooth(Yn, f, g, h, q, r, x0, v0, filter.end=n, predict.end=npe, missed=n/2, np=n/20)
 
# plot mean vector and state vector
  xss <- z$mean.smooth
  par(mfcol=c(m,1))
  for( i in 1:m ) {
    ymax <- as.integer(max(xss[i,],Xn[i,])+1)
    ymin <- as.integer(min(xss[i,],Xn[i,])-1)
    plot(Xn[i,], type='l', xlim=c(1,npe), ylim=c(ymin,ymax), xlab=paste(" red : mean.smooth[",i,",]  /  black : Xn[",i,",]"), ylab="")
    par(new=TRUE)
    plot(xss[i,], type='l', ylim=c(ymin,ymax), xlab="", ylab="", col=2)
  }

}

\keyword{ts}
