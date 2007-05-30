\name{blomar}
\alias{blomar}
\alias{print.blomar}
\title{Bayesian Method of Locally Stationary Multivariate AR Model Fitting}
\description{
  Locally fit multivariate autoregressive models to non-stationary time series by a Bayesian procedure.
}
\usage{
  blomar(y, max.order=NULL, span)
}
\arguments{
  \item{y}{A multivariate time series.}
  \item{max.order}{upper limit of the order of AR model. Default is 2*sqrt(n), where n is the length of the time series y.}
  \item{span}{length of basic local span.}
}
\value{
  \item{mean}{mean.}
  \item{var}{variance.}
  \item{bweight}{Bayesian weight.}
  \item{aic}{AIC with respect to the present data.}
  \item{arcoef}{AR coefficients. arcoef[[m]][i,j,k] shows the value of i-th row, j-th column, k-th order of m-th model.}
  \item{v}{innovation variance.}
  \item{eaic}{equivalent AIC of Bayesian model.}
  \item{init}{start point of the data fitted to the current model.}
  \item{end}{end point of the data fitted to the current model.}
}

\details{
  The basic AR model is given by

    \eqn{y(t) = A(1)y(t-1) + A(2)y(t-2) +...+ A(p)y(t-p) + u(t)},

  where p is order of the AR model and u(t) is innovation variance v.
}

\references{
  G.Kitagawa and H.Akaike (1978)
  A Procedure for the Modeling of Non-stationary Time Series.
  Ann. Inst. Statist. Math., 30, B, 351--363.

  H.Akaike (1978)
  A Bayesian Extension of The Minimum AIC Procedure of Autoregressive Model Fitting.
  Research Memo. NO.126. The institute of Statistical Mathematics.

  H.Akaike, G.Kitagawa, E.Arahata and F.Tada (1979)
  \emph{Computer Science Monograph, No.11, Timsac78.}
  The Institute of Statistical Mathematics.
}

\examples{
  data(Amerikamaru)
  blomar(Amerikamaru, max.order=10, span=300)
}

\keyword{ts}
