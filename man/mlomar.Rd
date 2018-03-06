\name{mlomar}
\alias{mlomar}
\alias{print.mlomar}
\title{Minimum AIC Method of Locally Stationary Multivariate AR Model Fitting}
\description{
  Locally fit multivariate autoregressive models to non-stationary time series
  by the minimum AIC procedure using the householder transformation.
}
\usage{
  mlomar(y, max.order = NULL, span, const = 0)
}
\arguments{
  \item{y}{a multivariate time series.}
  \item{max.order}{upper limit of the order of AR model, less than or equal to
    \eqn{n/2d} where \eqn{n} is the length and \eqn{d} is the dimension of the
    time series \code{y}. Default is
    \eqn{min(2 \sqrt{n}, n/2d)}{min(2*sqrt(n), n/2d)}.}
  \item{span}{length of basic local span. Let \eqn{m} denote \code{max.order},
    if \eqn{n-m-1} is less than or equal to \code{span} or \eqn{n-m-1-}\code{span}
    is less than \eqn{2md+}\code{const}, \code{span} is \eqn{n-m}.}
  \item{const}{integer. '\eqn{0}' denotes constant vector is not included as a
    regressor and '\eqn{1}' denotes constant vector is included as the first
    regressor.}
}
\value{
  \item{mean}{mean.}
  \item{var}{variance.}
  \item{ns}{the number of local spans.}
  \item{order}{order of the current model.}
  \item{aic}{AIC of the current model.}
  \item{arcoef}{AR coefficient matrices of the current model.
    \code{arcoef[[m]][i,j,k]} shows the value of \eqn{i}-th row, \eqn{j}-th
    column, \eqn{k}-th order of \eqn{m}-th model.}
  \item{v}{innovation variance of the current model.}
  \item{init}{initial point of the data fitted to the current model.}
  \item{end}{end point of the data fitted to the current model.}
  \item{npre}{data length of the preceding stationary block.}
  \item{nnew}{data length of the new block.}
  \item{order.mov}{order of the moving model.}
  \item{aic.mov}{AIC of the moving model.}
  \item{order.const}{order of the constant model.}
  \item{aic.const}{AIC of the constant model.}
}
\details{
  The data of length \eqn{n} are divided into \eqn{k} locally stationary spans,

    \deqn{|<-- n_1 -->|<-- n_2 -->|<-- n_3 -->| ..... |<-- n_k -->|}{%
          |<-- n1 -->|<-- n2 -->|<-- n3 -->| ..... |<-- nk -->|}

  where \eqn{n_i}{ni} \eqn{(i=1,\ldots,k)}{(i=1,...,k)} denoted the number of
  basic spans, each of length span, which constitute the \eqn{i}-th locally
  stationary span. At each local span, the process is represented by a
  stationary autoregressive model.
}

\references{
  G.Kitagawa and H.Akaike (1978)
  A Procedure for The Modeling of Non-Stationary Time Series.
  Ann. Inst. Statist. Math., 30, B, 351--363.

  H.Akaike, G.Kitagawa, E.Arahata and F.Tada (1979)
  \emph{Computer Science Monograph, No.11, Timsac78.}
  The Institute of Statistical Mathematics.
}

\examples{
data(Amerikamaru)
mlomar(Amerikamaru, max.order = 10, span = 300, const = 0)
}

\keyword{ts}
