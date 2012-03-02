\name{tvvar}
\alias{tvvar}
\title{Time Varying Variance}
\description{
  Estimate time-varying variance.
}
\usage{
  tvvar(y, trend.order, tau20=NULL, delta=NULL, plot=TRUE)
}

\arguments{
  \item{y}{univariate time series.}
  \item{trend.order}{trend order.}
  \item{tau20}{initial estimate of \code{tau2}.}
  \item{delta}{search width.}
  \item{plot}{logical. If TRUE (default) \code{normdat}, \code{ts}, \code{trend} and \code{noise} are plotted.}
}

\details{
    Assuming that \eqn{\sigma_{2m-1}^2 = \sigma_{2m}^2}{sigma_(2m-1)^2 = sigma_(2m)^2}, we define a transformed time series \eqn{s_1,...,s_{N/2}}{s_1,...,s_(N/2)} by 

      \deqn{s_m = y_{2m-1}^2 + y_{2m}^2,}{s_m = y_(2m-1)^2 + y_(2m)^2,}

    where \eqn{y_n} is a Gaussian white noise with mean \eqn{0} and variance \eqn{\sigma_n^2}{sigma_n^2}.

    \eqn{s_m} is distributed as a \eqn{\chi^2}{chi-square} distribution with \eqn{2} degrees of freedom, so the probability density function of \eqn{s_m} is given by

      \deqn{f(s) = \frac{1}{2\sigma^2} e^{-s/2\sigma^2}.}{f(s) = 1/(2sigma^2)  e^(-s/2sigma^2).}
 
    By further transformation

        \deqn{z_m = log(\frac{s_m}{2}),}{z_m = log(s_m/2),}

    the probability density function of \eqn{z_m} is given by

      \deqn{g(z) = \frac{1}{\sigma^2} exp{\{z-\frac{e^z}{\sigma^2}\}} = exp{\{(z-log\sigma^2)} - e^{(z-log\sigma^2)}\}.}{%
            g(z) = 1/sigma^2 exp{z - e^z/sigma^2} = exp{(z-log(sigma^2)) - e^(z-log(sigma^2))}.}

    Therefore the transformed time series is given by

      \deqn{z_m = log \sigma^2 + w_m,}{z_m = log(sigma^2) + w_m,}

    where \eqn{w_m} is a double exponential distribution with probability density function 

      \deqn{h(w) = exp{\{w-e^w\}}.}{h(w) = exp{w-e^w}.}

    In the space state model 

    \deqn{z_m = t_m + w_m}

    by identifying trend components of \eqn{z_m}, the log variance of original time series \eqn{y_n} is obtained.
}

\value{
  \item{tvvar}{time varying variance.}
  \item{normdat}{normalized data.}
  \item{ts}{transformed time series \eqn{s_m}.}
  \item{trend}{trend.}
  \item{noise}{residuals.}
  \item{tau2}{variance of the system noise tau2.}
  \item{sigma2}{variance of the observational noise.}
  \item{lkhood}{log-likelihood of the mode.}
  \item{aic}{AIC.}
}

\references{
  Kitagawa, G. (1993) \emph{Time series analysis programing (in Japanese)}. The Iwanami Computer Science Series.

  Kitagawa, G. and Gersch, W. (1996) \emph{Smoothness Priors Analysis of Time Series}. Lecture Notes in Statistics, No.116, Springer-Verlag.

  Kitagawa, G. and Gersch, W. (1985) \emph{A smoothness priors time varying AR coefficient modeling of nonstationary time series}. IEEE trans. on Automatic Control, AC-30, 48-56.
}

\examples{
  data(MYE1F) # an earthquake wave data
  z <- tvvar(MYE1F, 2, 6.6e-06, 1.0e-06)
  z$lkhood
  z$aic
}

\keyword{ts}
