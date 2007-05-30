\name{markov}
\alias{markov}
\title{Maximum Likelihood Computation of Markovian Model}
\description{
  Compute maximum likelihood estimates of Markovian model.
}
\usage{
markov(y, tmp.file=NULL)
}
\arguments{
  \item{y}{a multivariate time series.}
  \item{tmp.file}{temporary file name. If NULL (default) output no file.}
}
\details{
  This function is usually used with "simcon".
}
\value{
  \item{id}{id[i]=1 means that the i-th row of F contains free parameters.}
  \item{ir}{ir[i] denotes the position of the last non-zero element within the i-th row of F.}
  \item{ij}{ij[i] denotes the position of the i-th non-trivial row within F.}
  \item{ik}{ik[i] denotes the number of free parameters within the i-th non-trivial row of F.}
  \item{grad}{gradient vector.}
  \item{matFi}{initial estimate of the transition matrix (F).}
  \item{matF}{transition matrix (F).}
  \item{matG}{input matrix (G).}
  \item{davvar}{DAVIDON variance.}
  \item{arcoef}{AR coefficient matrices. arcoef[i,j,k] shows the value of i-th row, j-th column, k-th order.}
  \item{impuls}{impulse response matrices.}
  \item{macoef}{MA coefficient matrices. macoef[i,j,k] shows the value of i-th row, j-th column, k-th order.}
  \item{v}{inovation variance.}
  \item{aic}{AIC.}
}

\references{
  H.Akaike, E.Arahata and T.Ozaki (1975) \emph{Computer Science Monograph, No.5,
   Timsac74, A Time Series Analysis and Control Program Package (1)}.
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
  z <- mfilter(y,ar,"recursive")
  markov(z)
}

\keyword{ts}
