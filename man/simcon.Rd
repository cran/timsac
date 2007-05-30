\name{simcon}
\alias{simcon}
\title{Optimal Controller Design and Simulation}
\description{
  Produce optimal controller gain and simulate the controlled process.
}
\usage{
simcon(span, len, r, arcoef, impuls, v, weight)
}
\arguments{
  \item{span}{span of control performance evaluation.}
  \item{len}{length of experimental observation.}
  \item{r}{dimension of control input, less than or equal to d (dimension of a vector).}
  \item{arcoef}{matrices of autoregressive coefficients. arcoef[i,j,k] shows the value of i-th row, j-th column, k-th order.}
  \item{impuls}{impulse responce matrices.}
  \item{v}{covariance matrix of innovation.}
  \item{weight}{weighting matrix of performance.}
}
\details{
  The basic state space model is obtained from the autoregressive moving average model of
  a vector process y(t);

  \eqn{y(t) - A(1)y(t-1) -...- A(p)y(t-p) = u(t) - B(1)u(t-1) -...- B(p-1)u(t-p+1)},

  where A(i) (i=1,p)  are the autoregressive coefficients of the ARMA representation of y(t).
}

\value{
  \item{gain}{controller gain.}
  \item{ave}{average value of i-th component of y.}
  \item{var}{variance.}
  \item{std}{standard deviation.}
  \item{bc}{(p*d,r) sub matrices of impulse response matrices, where p is the order of the process, d is the dimension of the vector and
            r is the dimension of the control input.}
  \item{bd}{(p*d,d-r) sub matrices of impulse response matrices.}
}

\references{
  H.Akaike, E.Arahata and T.Ozaki (1975) \emph{Computer Science Monograph, No.6,
   Timsac74, A Time Series Analysis and Control Program Package (2)}.
  The Institute of Statistical Mathematics.
}

\examples{
  x <- matrix(rnorm(1000*2),1000,2)
  ma <- array(0,dim=c(2,2,2))
  ma[,,1] <- matrix(c( -1.0,  0.0,
                        0.0, -1.0), 2,2,byrow=TRUE)
  ma[,,2] <- matrix(c( -0.2,  0.0,
                       -0.1, -0.3), 2,2,byrow=TRUE)
  y <- mfilter(x,ma,"convolution")
  ar <- array(0,dim=c(2,2,3))
  ar[,,1] <- matrix(c( -1.0,  0.0,
                        0.0, -1.0), 2,2,byrow=TRUE)
  ar[,,2] <- matrix(c( -0.5, -0.2,
                       -0.2, -0.5), 2,2,byrow=TRUE)
  ar[,,3] <- matrix(c( -0.3, -0.05,
                       -0.1, -0.30), 2,2,byrow=TRUE)
  y <- mfilter(y,ar,"recursive")
  z <- markov(y)
  weight <-  matrix(c(0.0002,  0.0,
                      0.0,    2.9 ), 2,2,byrow=TRUE)
  simcon(span=50, len=700, r=1, z$arcoef, z$impuls, z$v, weight)
}

\keyword{ts}
