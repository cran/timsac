\name{tvspc}
\alias{tvspc}
\title{Time Evolution of Power Spectra of Time Varying AR model}
\description{
  Estimate the time evolution of the power spectra of time varying AR model.
}
\usage{
  tvspc(ar.order, sigma2, arcoef, var=NULL, span, nf=200)
}

\arguments{
  \item{ar.order}{AR order.}
  \item{sigma2}{variance of the observational noise.}
  \item{arcoef}{time varying AR coefficients.}
  \item{var}{time varying variance.}
  \item{span}{local stationary span.}
  \item{nf}{number of frequencies in evaluating spectrum.}
}

\value{
  \item{spec}{time varying spectrum.}
}

\references{
  Kitagawa, G. (1993) \emph{Time series analysis programing (in Japanese)}. The Iwanami Computer Science Series.

  Kitagawa, G. and Gersch, W. (1996) \emph{Smoothness Priors Analysis of Time Series}. Lecture Notes in Statistics, No.116, Springer-Verlag.

  Kitagawa, G. and Gersch, W. (1985) \emph{A smoothness priors time varying AR coefficient modeling of nonstationary time series}. IEEE trans. on Automatic Control, AC-30, 48-56.
}

\examples{
  data(MYE1F) # an earthquake wave data
  z <- tvar(MYE1F, 4, 2, 20, c(630,1026), 6.6e-06, 1.0e-06)
  spec <- tvspc(4, z$sigma2, z$arcoef,, 20)
  persp(spec$x, spec$y, spec$z, expand=0.5, theta=20, col = "lightblue", ticktype="detailed",
        xlab="f", ylab="n", zlab="log p(f)")
}

\keyword{ts}
