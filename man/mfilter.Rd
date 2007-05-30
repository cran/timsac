\name{mfilter}
\alias{mfilter}
\title{Linear Filtering on a Multivariate Time Series}
\description{
  Applies linear filtering to a multivariate time series.
}
\usage{
mfilter(x, filter, method=c("convolution","recursive"), init)
}
\arguments{
  \item{x}{a multivariate (m-dimensinal,n length) time series x[n,m].}
  \item{filter}{a array of filter coefficients. filter[i,j,k] shows the value of i-th row, j-th column, k-th order}
  \item{method}{either "convolution" or "recursive" (and can be abbreviated).
   If "convolution" a moving average is used: if "recursive" an autoregression is used.
   For convolution filters, the filter coefficients are foe past value only.}
  \item{init}{specifies the initial values of the time series just prior to the start value,
   in reverse time order. The default is a set of zeros.}
}
\details{
  This is a malutivariate version of "filter" function.

  Missing values are allowed in 'x' but not in 'filter'
  (where they would lead to missing values everywhere in the output).

  Note that there is an implied coefficient 1 at lag 0 in the recursive filter, which gives

     \eqn{y[i,]' =x[,i]' + f[,,1]*y[i-1,]' + ... +f[,,p]*y[i-p,]'},

  No check is made to see if recursive filter is invertible:
  the output maydiverge if it is not.

  The convolution filter is

     \eqn{y[i,]' = f[,,1]*x[i,]' + ... + f[,,p]*x[i-p+1,]'}
}
\value{
  \item{a time series object.}
}

\note{
  'convolve(, type="folter")' uses the FFT for computations and so \emph{may} be faster foe long filters on 
  univariate time series (and so the time alignment is unclear), nor does it handle missing values.
  'filter' is faster for a filter of length 100 on a series 1000, for examples.
}

\seealso{
  'convolve','arima.sim'
}

\examples{
#AR model simulation
  ar <- array(0,dim=c(3,3,2))
  ar[,,1] <- matrix(c(0.4,  0,   0.3,
                      0.2, -0.1, -0.5,
                      0.3,  0.1, 0),3,3,byrow=TRUE)
  ar[,,2] <- matrix(c(0,  -0.3,  0.5,
                      0.7, -0.4,  1,
                      0,   -0.5,  0.3),3,3,byrow=TRUE)
  x <- matrix(rnorm(100*3),100,3)
  y <- mfilter(x,ar,"recursive")

#Back to white noise
  ma <- array(0,dim=c(3,3,3))
  ma[,,1] <- diag(3)
  ma[,,2] <- -ar[,,1]
  ma[,,3] <- -ar[,,2]
  z <- mfilter(y,ma,"convolution")
  mulcor(z)

#AR-MA model simulation
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
}

\keyword{ts}
