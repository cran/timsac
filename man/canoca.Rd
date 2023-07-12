\name{canoca}
\alias{canoca}
\title{Canonical Correlation Analysis of Vector Time Series}
\description{
  Analyze canonical correlation of a d-dimensional multivariate time series.
}
\usage{
canoca(y)
}
\arguments{
  \item{y}{a multivariate time series.}
}
\details{
  First AR model is fitted by the minimum AIC procedure. The results are used to
  ortho-normalize the present and past variables. The present and future
  variables are tested successively to decide on the dependence of their
  predictors. When the last DIC (=chi-square - 2.0*N.D.F.) is negative the
  predictor of the variable is decided to be linearly dependent on the
  antecedents.
}
\value{
  \item{aic}{AIC.}
  \item{aicmin}{minimum AIC.}
  \item{order.maice}{MAICE AR model order.}
  \item{v}{innovation variance.}
  \item{arcoef}{autoregressive coefficients. \code{arcoef[i,j,k]} shows the
    value of \eqn{i}-th row, \eqn{j}-th column, \eqn{k}-th order.}
  \item{nc}{number of cases.}
  \item{future}{number of variable in the future set.}
  \item{past}{number of variables in the past set.}
  \item{cweight}{future set canonical weight.}
  \item{canocoef}{canonical R.}
  \item{canocoef2}{R-squared.}
  \item{chisquar}{chi-square.}
  \item{ndf}{N.D.F.}
  \item{dic}{DIC.}
  \item{dicmin}{minimum DIC.}
  \item{order.dicmin}{order of minimum DIC.}
  \item{matF}{the transition matrix \eqn{F}.}
  \item{vectH}{structural characteristic vector \eqn{H} of the canonical
    Markovian representation.}
  \item{matG}{the estimate of the input matrix \eqn{G}.}
  \item{vectF}{matrix \eqn{F} in vector form.}
}

\references{
  H.Akaike, E.Arahata and T.Ozaki (1975) \emph{Computer Science Monograph, No.5,
   Timsac74, A Time Series Analysis and Control Program Package (1)}.
  The Institute of Statistical Mathematics.
}

\examples{
ar <- array(0, dim = c(3,3,2))
ar[, , 1] <- matrix(c(0.4,  0,   0.3,
                      0.2, -0.1, -0.5,
                      0.3,  0.1, 0), nrow = 3, ncol = 3, byrow= TRUE)
ar[, , 2] <- matrix(c(0,  -0.3,  0.5,
                      0.7, -0.4,  1,
                      0,   -0.5,  0.3), nrow = 3, ncol = 3, byrow = TRUE)
x <- matrix(rnorm(1000*3), nrow = 1000, ncol = 3)
y <- mfilter(x, ar, "recursive")
z <- canoca(y)
z$arcoef
}

\keyword{ts}
