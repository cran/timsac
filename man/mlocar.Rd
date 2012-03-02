\name{mlocar}
\alias{mlocar}
\title{Minimum AIC Method of Locally Stationary AR Model Fitting; Scalar Case}
\description{
  Locally fit autoregressive models to non-stationary time series by minimum AIC procedure.
}
\usage{
  mlocar(y, max.order=NULL, span, const=0, plot=TRUE)
}
\arguments{
  \item{y}{a univariate time series.}
  \item{max.order}{upper limit of the order of AR model. Default is \eqn{2 \sqrt{n}}{2*sqrt(n)}, where \eqn{n} is the length of the time series \code{y}.}
  \item{span}{length of the basic local span.}
  \item{const}{integer. \code{0} denotes constant vector is not included as a regressor
        and \code{1} denotes constant vector is included as the first regressor.}
  \item{plot}{logical. If TRUE (default) spectrums \code{pspec} are plotted.}
}
\value{
  \item{mean}{mean.}
  \item{var}{variance.}
  \item{ns}{the number of local spans.}
  \item{order}{order of the current model.}
  \item{arcoef}{AR coefficients of current model.}
  \item{v}{innovation variance of the current model.}
  \item{init}{initial point of the data fitted to the current model.}
  \item{end}{end point of the data fitted to the current model.}
  \item{pspec}{power spectrum.}
  \item{npre}{data length of the preceding stationary block.}
  \item{nnew}{data length of the new block.}
  \item{order.mov}{order of the moving model.}
  \item{v.mov}{innovation variance of the moving model.}
  \item{aic.mov}{AIC of the moving model.}
  \item{order.const}{order of the constant model.}
  \item{v.const}{innovation variance of the constant model.}
  \item{aic.const}{AIC of the constant model.}
}
\details{
  The data of length \eqn{n} are devided into \eqn{k} locally stationary spans,

    \deqn{|<-- n_1 -->|<-- n_2 -->|<-- n_3 -->| ..... |<-- n_k -->|}{  |<-- n1 -->|<-- n2 -->|<-- n3 -->| ..... |<-- nk -->|}

  where \eqn{n_i}{ni} (\eqn{i=1,\ldots,k}{i=1,...,k}) denotes the number of basic spans, each of length span,
  which constitute the \eqn{i}-th locally stationary span. At each local span,
  the process is represented by a stationary autoregressive model.
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
  data(locarData)
  z <- mlocar(locarData, max.order=10, span=300, const=0)
  z$arcoef
}

\keyword{ts}
