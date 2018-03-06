\name{wnoise}
\alias{wnoise}
\title{White Noise Generator}
\description{
  Generate approximately Gaussian vector white noise.
}
\usage{
  wnoise(len, perr, plot = TRUE)
}
\arguments{
  \item{len}{length of white noise record.}
  \item{perr}{prediction error.}
  \item{plot}{logical. If \code{TRUE} (default), white noises are plotted.}
}
\value{
  \code{wnoise} gives white noises.
}

\references{
  H.Akaike and T.Nakagawa (1988) \emph{Statistical Analysis and Control
    of Dynamic Systems.} Kluwer Academic publishers.
}

\examples{
# Example 1
wnoise(len = 100, perr = 1)

# Example 2
v <- matrix(c(1,  0,  0,
              0,  2,  0,
              0,  0,  3), nrow = 3, ncol = 3, byrow = TRUE)
wnoise(len = 20, perr = v)
}

\keyword{ts}
