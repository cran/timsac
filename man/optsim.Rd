\name{optsim}
\alias{optsim}
\title{Optimal Control Simulation}
\description{
  Perform optimal control simulation and evaluate the means and variances of the
  controlled and manipulated variables X and Y.
}
\usage{
  optsim(y, max.order = NULL, ns, q, r, noise = NULL, len, plot = TRUE)
}
\arguments{
  \item{y}{a multivariate time series.}
  \item{max.order}{upper limit of model order. Default is \eqn{2 \sqrt{n}}{2*sqrt(n)}, where \eqn{n} is the length of the time series \code{y}.}
  \item{ns}{number of steps of simulation.}
  \item{q}{positive definite matrix \eqn{Q}.}
  \item{r}{positive definite matrix \eqn{R}.}
  \item{noise}{noise. If not provided, Gaussian vector white noise with the
    length \code{len} is generated.}
  \item{len}{length of white noise record.}
  \item{plot}{logical. If \code{TRUE} (default), controlled variables \eqn{X}
    and  manipulated variables \eqn{Y} are plotted.}
}
\value{
  \item{trans}{first \eqn{m} columns of transition matrix, where \eqn{m} is the
    number of controlled variables.}
  \item{gamma}{gamma matrix.}
  \item{gain}{gain matrix.}
  \item{convar}{controlled variables \eqn{X}.}
  \item{manvar}{manipulated variables \eqn{Y}.}
  \item{xmean}{mean of \eqn{X}.}
  \item{ymean}{mean of \eqn{Y}.}
  \item{xvar}{variance of \eqn{X}.}
  \item{yvar}{variance of \eqn{Y}.}
  \item{x2sum}{sum of \eqn{X^2}.}
  \item{y2sum}{sum of \eqn{Y^2}.}
  \item{x2mean}{mean of \eqn{X^2}.}
  \item{y2mean}{mean of \eqn{Y^2}.}
}

\references{
  H.Akaike and T.Nakagawa (1988) \emph{Statistical Analysis and Control of
  Dynamic Systems.} Kluwer Academic publishers.
}

\examples{
# Multivariate Example Data
ar <- array(0, dim = c(3,3,2))
ar[, , 1] <- matrix(c(0.4,  0,    0.3,
                      0.2, -0.1, -0.5,
                      0.3,  0.1, 0), nrow = 3, ncol = 3, byrow = TRUE)
ar[, , 2] <- matrix(c(0,  -0.3,  0.5,
                      0.7, -0.4,  1,
                      0,   -0.5,  0.3), nrow = 3, ncol = 3, byrow = TRUE)
x <- matrix(rnorm(200*3), nrow = 200, ncol = 3)
y <- mfilter(x, ar, "recursive")
q.mat <- matrix(c(0.16,0,0,0.09), nrow = 2, ncol = 2)
r.mat <- as.matrix(0.001)
optsim(y, max.order = 10, ns = 20, q = q.mat, r = r.mat, len = 20)
}

\keyword{ts}
