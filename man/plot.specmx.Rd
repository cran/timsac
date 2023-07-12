\name{plot.specmx}
\alias{plot.specmx}
\title{Plot Spectrum}
\description{
  Plot spectrum returned by \code{\link{mulspe}} and \code{\link{mulrsp}}.
  On and lower diagonal are real parts, and upper diagonal are imaginary parts. 
}
\usage{
\method{plot}{specmx}(x, plot.scale = TRUE, \dots)
}
\arguments{
  \item{x}{An object of class \code{"specmx"}.}
  \item{plot.scale}{logical. IF \code{TRUE}, the common range of the
    \eqn{y}-axis is used.}
  \item{\dots}{further graphical parameters may also be supplied as arguments.}
}
\keyword{ts}
