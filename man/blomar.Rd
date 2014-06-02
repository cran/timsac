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
  \item{max.order}{upper limit of the order of AR model, less than or equal to \eqn{n/(2*d)} where \eqn{n} is the length and \eqn{d} is the dimension of the time series \code{y}.
 Default is \eqn{min(2 \sqrt{n}, n/(2*d))}{min(2*sqrt(n), n/(2*d))}.}
  \item{span}{length of basic local span.
 If \eqn{n}-\code{max.order}-1 is less than or equal to \code{span} or \eqn{n}-\code{max.order}-1-\code{span} is less than 2*\code{max.order}*\eqn{d}, \code{span} is \eqn{n}-\code{max.order}.}
}
\value{
  \item{mean}{mean.}
  \item{var}{variance.}
  \item{bweight}{Bayesian weight.}
  \item{aic}{AIC with respect to the present data.}
  \item{arcoef}{AR coefficients. \code{arcoef[[m]][i,j,k]} shows the value of \eqn{i}-th row, \eqn{j}-th column, \eqn{k}-th order of \eqn{m}-th model.}
  \item{v}{innovation variance.}
  \item{eaic}{equivalent AIC of Bayesian model.}
  \item{init}{start point of the data fitted to the current model.}
  \item{end}{end point of the data fitted to the current model.}
}

\details{
  The basic AR model is given by
    \deqn{y(t) = A(1)y(t-1) + A(2)y(t-2) +...+ A(p)y(t-p) + u(t),}
  where \eqn{p} is order of the AR model and \eqn{u(t)} is innovation variance \code{v}.
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
