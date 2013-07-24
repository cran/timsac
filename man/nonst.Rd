\name{nonst}
\alias{nonst}
\title{Non-stationary Power Spectrum Analysis}
\description{
  Locally fit autoregressive models to non-stationary time series by AIC criterion.
}
\usage{
nonst(y, span, max.order=NULL, plot=TRUE)
}
\arguments{
  \item{y}{a univariate time series.}
  \item{span}{length of the basic local span.}
  \item{max.order}{highest order of AR model. Default is \eqn{2 \sqrt{n}}{2*sqrt(n)}, where \eqn{n} is the length of the time series \code{y}.}
  \item{plot}{logical. If TRUE (the default) spectrums are plotted.}
}
\details{
  The basic AR model is given by
    \deqn{y(t) = A(1)y(t-1) + A(2)y(t-2) +...+ A(p)y(t-p) + u(t),}
  where \eqn{p} is order of the AR model and \eqn{u(t)} is innovation variance.
  AIC is defined by
    \deqn{AIC = n \log(det(sd)) + 2k,}{AIC = nlog(det(sd)) + 2k,}
  where \eqn{n} is the length of data, \eqn{sd} is the estimates of the innovation variance
  and \eqn{k} is the number of parameter.
}
\value{
  \item{ns}{the number of local spans.}
  \item{arcoef}{AR coefficients.}
  \item{v}{innovation variance.}
  \item{aic}{AIC.}
  \item{daic21}{= AIC2-AIC1.}
  \item{daic}{= \code{daic21}\eqn{/n} (\eqn{n} is the length of the current model).}
  \item{init}{start point of the data fitted to the current model.}
  \item{end}{end point of the data fitted to the current model.}
  \item{pspec}{power spectrum.}
}

\references{
  H.Akaike, E.Arahata and T.Ozaki (1976) \emph{Computer Science Monograph, No.6,
   Timsac74 A Time Series Analysis and Control Program Package (2)}.
  The Institute of Statistical Mathematics.
}

\examples{
# Non-stationary Test Data
  data(nonstData)
  nonst(nonstData, span=700, max.order=49)
}

\keyword{ts}
