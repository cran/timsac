\name{exsar}
\alias{exsar}
\title{Exact Maximum Likelihood Method of Scalar AR Model Fitting}
\description{
  Produce exact maximum likelihood estimates of the parameters of a scalar AR model.
}
\usage{
  exsar(y, max.order=NULL, plot=FALSE, tmp.file=NULL)
}
\arguments{
  \item{y}{a univariate time sries.}
  \item{max.order}{upper limit of AR order. Default is 2*sqrt(n), where n is the length of the time series y.}
  \item{plot}{logical. If TRUE daic is plotted.}
  \item{tmp.file}{a character string naming a file written intermediate results of minimization by DAVIDON-FLETCHER-POWELL procedure.
        If NULL (default) output no file.}
}
\value{
  \item{mean}{mean.}
  \item{var}{variance.}
  \item{v}{innovation variance.}
  \item{aic}{AIC.}
  \item{aicmin}{minimum AIC.}
  \item{daic}{AIC-aicmin.}
  \item{order.maice}{order of minimum AIC.}
  \item{v.maice}{MAICE innovation variance.}
  \item{arcoef.maice}{MAICE AR coefficients.}
  \item{v.mle}{maximum likelihood estimates of innovation variance.}
  \item{arcoef.mle}{maximum likelihood estimates of AR coefficients.}
}
\details{
  The AR model is given by

    \eqn{y(t) = a(1)y(t-1) + .... + a(p)y(t-p) + u(t)}

  where p is AR order and u(t) is a zero mean white noise.
}

\references{
  H.Akaike, G.Kitagawa, E.Arahata and F.Tada (1979)
  \emph{Computer Science Monograph, No.11, Timsac78.}
  The Institute of Statistical Mathematics.
}

\examples{
  data(Canadianlynx)
  z <- exsar(Canadianlynx, max.order=14)
  z$arcoef.maice
  z$arcoef.mle
}

\keyword{ts}
