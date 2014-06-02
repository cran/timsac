\name{markov}
\alias{markov}
\title{Maximum Likelihood Computation of Markovian Model}
\description{
  Compute maximum likelihood estimates of Markovian model.
}
\usage{
markov(y)
}
\arguments{
  \item{y}{a multivariate time series.}
}
\details{
  This function is usually used with \code{\link{simcon}}.
}
\value{
  \item{id}{\code{id[i]}\eqn{=1} means that the \eqn{i}-th row of \eqn{F} contains free parameters.}
  \item{ir}{\code{ir[i]} denotes the position of the last non-zero element within the \eqn{i}-th row of \eqn{F}.}
  \item{ij}{\code{ij[i]} denotes the position of the \eqn{i}-th non-trivial row within \eqn{F}.}
  \item{ik}{\code{ik[i]} denotes the number of free parameters within the \eqn{i}-th non-trivial row of \eqn{F}.}
  \item{grad}{gradient vector.}
  \item{matFi}{initial estimate of the transition matrix \eqn{F}.}
  \item{matF}{transition matrix \eqn{F}.}
  \item{matG}{input matrix \eqn{G}.}
  \item{davvar}{DAVIDON variance.}
  \item{arcoef}{AR coefficient matrices. \code{arcoef[i,j,k]} shows the value of \eqn{i}-th row, \eqn{j}-th column, \eqn{k}-th order.}
  \item{impulse}{impulse response matrices.}
  \item{macoef}{MA coefficient matrices. \code{macoef[i,j,k]} shows the value of \eqn{i}-th row, \eqn{j}-th column, \eqn{k}-th order.}
  \item{v}{innovation variance.}
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
