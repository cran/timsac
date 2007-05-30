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
bsubst(y,mtype,lag=NULL,nreg,reg=NULL,term.lag=NULL,cstep=5,plot=TRUE)
}

\arguments{
  \item{y}{a univariate time series.}

  \item{mtype}{model type. Allowed values are

    \code{1} : (autoregressive model),

    \code{2} : (polinomial type non-linear model, lag's read in),

    \code{3} : (polinomial type non-linear model, lag's automatically set) and

    \code{4} : (AR-model with polinomial mean value function).

    \code{5},\code{6} and \code{7} are originaly defined but omitted here.}

  \item{lag}{maximum time lag. Default is 2*sqrt(n), where n is the length of the time series y.}

  \item{nreg}{number of regressors.}

  \item{reg}{specification of regressor (mtype = 2).

    i-th regressor is defined by \eqn{z(n-L1(i)) * z(n-L2(i)) * z(n-L3(i))},
    where L1(i) is reg(1,i), L2(i) is reg(2,i) and L3(i) is reg(3,i).

    0-lag term z(n-0) is replaced by the constant 1.}

  \item{term.lag}{maximum time lag specify the regressors (L1(i),L2(i),L3(i)) (i=1,...,nreg) (mtype = 3).

    i-th regressor is defined by \eqn{z(n-L1(i)) * z(n-L2(i)) * z(n-L3(i))},
    where 0-lag term z(n-0) is replaced by the constant 1.

    term.lag(1) : maximum time lag of linear term

    term.lag(2) : maximum time lag of squared term

    term.lag(3) : maximum time lag of quadratic cross term

    term.lag(4) : maximum time lag of cubic term

    term.lag(5) : maximum time lag of cubic cross term.}

  \item{cstep}{prediction errors checking (up to cstep-steps ahead) is requested.}

  \item{plot}{logical. If TRUE (default) daic, pre.err and peautcor are plotted.}
}
\details{
  The AR model is given by (mtype = 2)

    \eqn{y(t) = a(1)y(t-1) + .... + a(p)y(t-p) + u(t)}.

  The non-linear model is given by ( mtype = 2,3 )

  \eqn{y(t) = a(1)z(t,1) + a(2)z(t,2) +...+ a(p)z(t,p) + u(t)}.

  Where p is AR order and u(t) is Gaussian white noise with mean 0 and variance v(p).
}
\value{
  \item{ymean}{mean of y.}
  \item{yvar}{variance of y.}
  \item{v}{innovation variance.}
  \item{aic}{AIC(m), (m=0,...,nreg).}
  \item{aicmin}{minimum AIC.}
  \item{daic}{AIC(m)-aicmin (m=0,...,nreg).}
  \item{order.maice}{order of minimum AIC.}
  \item{v.maice}{innovation variance attained at order.maice.}
  \item{arcoef.maice}{AR coefficients attained at order.maice.}
  \item{v.bay}{residual variance of Bayesian model.}
  \item{aic.bay}{AIC of Bayesian model.}
  \item{np.bay}{equivalent number of parameters.}
  \item{arcoef.bay}{AR coefficients of Bayesian model.}
  \item{ind.c}{index of parcor2 in order of increasing magnitude.}
  \item{parcor2}{square of partial correlations (normalisedby multiplying N).}
  \item{damp}{binomial type damper.}
  \item{bweight}{final Bayesian weights of partial correlations.}
  \item{parcor.bay}{partial correlations of the Bayesian model.}
  \item{eicmin}{minimum EIC.}
  \item{esum}{whole subset regression models.}
  \item{npmean}{mean of number of parameter.}
  \item{npmean.nreg}{(=npmean/nreg).}
  \item{perr}{prediction error.}
  \item{mean}{mean.}
  \item{var}{variance.}
  \item{skew}{skewness.}
  \item{peak}{peakedness.}
  \item{peautcor}{autocorrelation function of 1-step ahead prediction error.}
  \item{pspec}{power spectrum (mtype = 1).}
}

\references{
  H.Akaike, G.Kitagawa, E.Arahata and F.Tada (1979)
  \emph{Computer Science Monograph, No.11, Timsac78.}
  The Institute of Statistical Mathematics.
}

\examples{
  data(Canadianlynx)
  Regressor <- matrix(c( 1, 0, 0, 2, 0, 0, 3, 0,  0, 4, 0,  0, 5, 0,  0, 6, 0, 0, 7,
                         0, 0, 8, 0, 0, 9, 0, 0, 10, 0, 0, 11, 0, 0, 12, 0, 0, 1, 1,
                         0, 2, 2, 0, 1, 2, 0, 3,  3, 0, 1,  1, 1, 2,  2, 2, 3, 3, 3 ), 3,19)
  z <- bsubst(Canadianlynx, mtype=2, lag=12, nreg=19, reg=Regressor, cstep=5 )
  z$arcoef.bay
}

\keyword{ts}
