\name{optsim}
\alias{optsim}
\title{Optimal Control Simulation}
\description{
  Perform optimal control simulation and evaluate the means and variances of the controlled and manipulated variables X and Y.
}
\usage{
  optsim(y, max.order=NULL, ns, q, r, noise=NULL, len, plot=TRUE)
}
\arguments{
  \item{y}{a multivariate time series.}
  \item{max.order}{upper limit of model order. Default is 2*sqrt(n), where n is the length of the time series y.}
  \item{ns}{number of steps of simulation.}
  \item{q}{positive definite matrix \eqn{Q}.}
  \item{r}{positive definite matrix \eqn{R}.}
  \item{noise}{noise. If not provided, Gaussian vector white noise with the length len is generated.}
  \item{len}{length of white noise record.}
  \item{plot}{logical. If TRUE (default) controlled variables \eqn{X} and  manipulated variables \eqn{Y} are plotted.}
}
\value{
  \item{trans}{first ir columns of transition matrix, where ir is the number of controlled variables.}
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
  H.Akaike and T.Nakagawa (1988) \emph{Statistical Analysis and Control
    of Dynamic Systems.} Kluwer Academic publishers.
}

\examples{
# Multivariate Example Data
  ar <- array(0,dim=c(3,3,2))
  ar[,,1] <- matrix(c(0.4,  0,   0.3,
                      0.2, -0.1, -0.5,
                      0.3,  0.1, 0),3,3,byrow=TRUE)
  ar[,,2] <- matrix(c(0,  -0.3,  0.5,
                      0.7, -0.4,  1,
                      0,   -0.5,  0.3),3,3,byrow=TRUE)
  x <- matrix(rnorm(200*3),200,3)
  y <- mfilter(x,ar,"recursive")
  q <- matrix(c(0.16,0,0,0.09), 2, 2)
  r <- matrix(0.001, 1, 1)
  optsim(y, max.order=10, ns=20, q, r, len=20)
}

\keyword{ts}
