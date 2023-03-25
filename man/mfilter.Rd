\name{mfilter}
\alias{mfilter}
\title{Linear Filtering on a Multivariate Time Series}
\description{
  Applies linear filtering to a multivariate time series.
}
\usage{
mfilter(x, filter, method = c("convolution","recursive"), init)
}
\arguments{
  \item{x}{a multivariate (\eqn{m}-dimensional, \eqn{n} length) time series
    \eqn{x[n,m]}.}
  \item{filter}{an array of filter coefficients. \code{filter[i,j,k]} shows the
    value of \eqn{i}-th row, \eqn{j}-th column, \eqn{k}-th order}
  \item{method}{either "convolution" or "recursive" (and can be abbreviated).
    If "convolution" a moving average is used: if "recursive" an autoregression
    is used. For convolution filters, the filter coefficients are for past value
    only.}
  \item{init}{specifies the initial values of the time series just prior to the
    start value, in reverse time order. The default is a set of zeros.}
}
\details{
  This is a multivariate version of "filter" function.
  Missing values are allowed in '\code{x}' but not in '\code{filter}'
  (where they would lead to missing values everywhere in the output).
  Note that there is an implied coefficient \eqn{1} at lag \eqn{0} in the
  recursive filter, which gives
    \deqn{y[i,]' =x[,i]' + f[,,1] \times y[i-1,]' + ... +f[,,p] \times
 y[i-p,]',}{y[i,]' =x[,i]' + f[,,1]*y[i-1,]' + ... +f[,,p]*y[i-p,]',}

  No check is made to see if recursive filter is invertible: the output may
  diverge if it is not. The convolution filter is
     \deqn{y[i,]' = f[,,1] \times x[i,]' + ... + f[,,p] \times
 x[i-p+1,]'.}{y[i,]' = f[,,1]*x[i,]' + ... + f[,,p]*x[i-p+1,]'.}

}
\value{
  \code{mfilter} returns a time series object.
}

\note{
  '\code{convolve(, type="filter")}' uses the FFT for computations and so may be
  faster for long filters on univariate time series (and so the time alignment
  is unclear), nor does it handle missing values. 'filter' is faster for a
  filter of length 100 on a series 1000, for examples.
}

\seealso{
  \code{\link{convolve}}, \code{\link{arima.sim}}
}

\examples{
#AR model simulation
ar <- array(0, dim = c(3,3,2))
ar[, , 1] <- matrix(c(0.4,  0,   0.3,
                      0.2, -0.1, -0.5,
                      0.3,  0.1, 0), nrow = 3, ncol = 3, byrow = TRUE)
ar[, , 2] <- matrix(c(0,  -0.3,  0.5,
                      0.7, -0.4,  1,
                      0,   -0.5,  0.3), nrow = 3, ncol = 3, byrow = TRUE)
x <- matrix(rnorm(100*3), nrow = 100, ncol = 3)
y <- mfilter(x, ar, "recursive")

#Back to white noise
ma <- array(0, dim = c(3,3,3))
ma[, , 1] <- diag(3)
ma[, , 2] <- -ar[, , 1]
ma[, , 3] <- -ar[, , 2]
z <- mfilter(y, ma, "convolution")
mulcor(z)

#AR-MA model simulation
x <- matrix(rnorm(1000*2), nrow = 1000, ncol = 2)
ma <- array(0, dim = c(2,2,2))
ma[, , 1] <- matrix(c( -1.0,  0.0,
                        0.0, -1.0), nrow = 2, ncol = 2, byrow = TRUE)
ma[, , 2] <- matrix(c( -0.2,  0.0,
                       -0.1, -0.3), nrow = 2, ncol = 2, byrow = TRUE)
y <- mfilter(x, ma, "convolution")

ar <- array(0, dim = c(2,2,3))
ar[, , 1] <- matrix(c( -1.0,  0.0,
                        0.0, -1.0), nrow = 2, ncol = 2, byrow = TRUE)
ar[, , 2] <- matrix(c( -0.5, -0.2,
                       -0.2, -0.5), nrow = 2, ncol = 2, byrow = TRUE)
ar[, , 3] <- matrix(c( -0.3, -0.05,
                       -0.1, -0.30), nrow = 2, ncol = 2, byrow = TRUE)
z <- mfilter(y, ar, "recursive")
}

\keyword{ts}
