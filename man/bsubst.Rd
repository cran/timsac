\name{bsubst}
\alias{bsubst}
\title{Bayesian Type All Subset Analysis}
\description{
  Produce Bayesian estimates of time series models such as pure AR models, AR models
  with non-linear terms, AR models with polynomial type mean value functions, etc.
  The goodness of fit of a model is checked by the analysis of several steps ahead
  prediction errors.
}
\usage{
bsubst(y, mtype, lag=NULL, nreg, reg=NULL, term.lag=NULL, cstep=5,
       plot=TRUE)
}

\arguments{
  \item{y}{a univariate time series.}

  \item{mtype}{model type. Allowed values are
    \tabular{rl}{
      1 : \tab autoregressive model,\cr
      2 : \tab polynomial type non-linear model (lag's read in),\cr
      3 : \tab polynomial type non-linear model (lag's automatically set),\cr
      4 : \tab AR-model with polynomial mean value function,\cr
      5,6,7 : \tab originally defined but omitted here.
    }}

  \item{lag}{maximum time lag. Default is \eqn{2 \sqrt(n)}{2*sqrt(n)}, where \eqn{n} is the length of the time series \code{y}.}

  \item{nreg}{number of regressors.}

  \item{reg}{specification of regressor (\code{mtype} = 2).\cr
    \eqn{i}-th regressor is defined by \eqn{z(n-L1(i)) \times z(n-L2(i)) \times z(n-L3(i))}{z(n-L1(i)) * z(n-L2(i)) * z(n-L3(i))},
    where \eqn{L1(i)} is \code{reg(1,i)}, \eqn{L2(i)} is \code{reg(2,i)} and \eqn{L3(i)} is \code{reg(3,i)}.
    \eqn{0}-lag term \eqn{z(n-0)} is replaced by the constant \eqn{1}.}

  \item{term.lag}{maximum time lag specify the regressors (\eqn{L1(i),L2(i),L3(i)}) (i=1,\ldots,\code{nreg}) (\code{mtype} = 3).
    \tabular{rl}{
      term.lag(1) : \tab maximum time lag of linear term\cr
      term.lag(2) : \tab maximum time lag of squared term\cr
      term.lag(3) : \tab maximum time lag of quadratic crosses term\cr
      term.lag(4) : \tab maximum time lag of cubic term\cr
      term.lag(5) : \tab maximum time lag of cubic cross term.
    }}

  \item{cstep}{prediction errors checking (up to \code{cstep}-steps ahead) is requested. (\code{mtype} = 1,2,3).}

  \item{plot}{logical. If TRUE (default) \code{daic}, \code{perr} and \code{peautcor} are plotted.}
}
\details{
  The AR model is given by ( \code{mtype} = 2 )\cr
    \deqn{y(t) = a(1)y(t-1) + ... + a(p)y(t-p) + u(t).}
  The non-linear model is given by ( \code{mtype} = 2,3 )\cr
    \deqn{y(t) = a(1)z(t,1) + a(2)z(t,2) + ... + a(p)z(t,p) + u(t).}
  Where \eqn{p} is AR order and \eqn{u(t)} is Gaussian white noise with mean \eqn{0} and variance \eqn{v(p)}.
}
\value{
  \item{ymean}{mean of \code{y}.}
  \item{yvar}{variance of \code{y}.}
  \item{v}{innovation variance.}
  \item{aic}{AIC(m), (m=0,...,\code{nreg}).}
  \item{aicmin}{minimum AIC.}
  \item{daic}{AIC(m)-\code{aicmin} (m=0,...,\code{nreg}).}
  \item{order.maice}{order of minimum AIC.}
  \item{v.maice}{innovation variance attained at \code{order.maice}.}
  \item{arcoef.maice}{AR coefficients attained at \code{order.maice}.}
  \item{v.bay}{residual variance of Bayesian model.}
  \item{aic.bay}{AIC of Bayesian model.}
  \item{np.bay}{equivalent number of parameters.}
  \item{arcoef.bay}{AR coefficients of Bayesian model.}
  \item{ind.c}{index of \code{parcor2} in order of increasing magnitude.}
  \item{parcor2}{square of partial correlations (normalized by multiplying N).}
  \item{damp}{binomial type damper.}
  \item{bweight}{final Bayesian weights of partial correlations.}
  \item{parcor.bay}{partial correlations of the Bayesian model.}
  \item{eicmin}{minimum EIC.}
  \item{esum}{whole subset regression models.}
  \item{npmean}{mean of number of parameter.}
  \item{npmean.nreg}{=\code{npmean}/\code{nreg}.}
  \item{perr}{prediction error.}
  \item{mean}{mean.}
  \item{var}{variance.}
  \item{skew}{skewness.}
  \item{peak}{peakedness.}
  \item{peautcor}{autocorrelation function of 1-step ahead prediction error.}
  \item{pspec}{power spectrum (\code{mtype} = 1).}
}

\references{
  H.Akaike, G.Kitagawa, E.Arahata and F.Tada (1979)
  \emph{Computer Science Monograph, No.11, Timsac78.}
  The Institute of Statistical Mathematics.
}

\examples{
  data(Canadianlynx)
  Regressor <- matrix(
       c( 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 1, 2, 1, 3, 1, 2, 3,
          0, 0, 0, 0, 0, 0, 0, 0, 0,  0,  0,  0, 1, 2, 2, 3, 1, 2, 3,
          0, 0, 0, 0, 0, 0, 0, 0, 0,  0,  0,  0, 0, 0, 0, 0, 1, 2, 3 ),
        3,19, byrow=TRUE)

  z <- bsubst(Canadianlynx, 2, 12, 19, Regressor)

  z$arcoef.bay
}

\keyword{ts}
