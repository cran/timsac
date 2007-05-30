\name{prdctr}
\alias{prdctr}
\alias{print.prdctr}
\title{Prediction Program}
\description{
  Operate on a real record of a vector process and compute predicted values.
}
\usage{
prdctr(y, r, s, h, arcoef, macoef=NULL, impuls=NULL, v, plot=TRUE)
}
\arguments{
  \item{y}{a univariate time series or A multivariate time series.}
  \item{r}{one step ahead prediction starting position \eqn{R}.}
  \item{s}{long range forecast starting position \eqn{S}.}
  \item{h}{maximum span of long range forecast \eqn{H}.}
  \item{arcoef}{AR coefficient matrices.}
  \item{macoef}{MA coefficient matrices.}
  \item{impuls}{impulse response matrices.}
  \item{v}{inovation variance.}
  \item{plot}{logical. If TRUE (default) the real data and predicted values are plotted.}
}
\details{
  One step ahead Prediction starts at time \eqn{R} and ends at time \eqn{S}.
  Prediction is continued without new observations until time \eqn{S+H}.
  Basic model is the autoregressive moving average model of \eqn{y(t)} which is given by

    \eqn{y(t) - A(t)y(t-1) -...- A(p)y(t-p) = u(t) - B(1)u(t-1) -...- B(q)u(t-q)},

  where p is AR order and q is MA order.
}
\value{
  \item{predct}{predicted values : predct(i) \eqn{(r<=i<=s+h)}.}
  \item{ys}{predct(i) - y(i) \eqn{(r<=i<=n)}.}
  \item{pstd}{predct(i) + (standard deviation) \eqn{(s<=i<=s+h)}.}
  \item{p2std}{predct(i) + 2*(standard deviation) \eqn{(s<=i<=s+h)}.}
  \item{p3std}{predct(i) + 3*(standard deviation) \eqn{(s<=i<=s+h)}.}
  \item{mstd}{predct(i) - (standard deviation) \eqn{(s<=i<=s+h)}.}
  \item{m2std}{predct(i) - 2*(standard deviation) \eqn{(s<=i<=s+h)}.}
  \item{m3std}{predct(i) - 3*(standard deviation) \eqn{(s<=i<=s+h)}.}
}

\references{
  H.Akaike, E.Arahata and T.Ozaki (1975) \emph{Computer Science Monograph, No.6,
   Timsac74, A Time Series Analysis and Control Program Package (2)}.
  The Institute of Statistical Mathematics.
}

\examples{
  # "arima.sim" is a function in "stats".
  # Note that the sign of MA coefficient is opposite from that in "timsac".
  yy <- arima.sim(list(order=c(2,0,1), ar = c(0.64,-0.8), ma=c(-0.5)), n=350)
  y1 <- yy[51:300]
  z <- autoarmafit(y1)
  ar <- z$model[[1]]$arcoef
  ma <- z$model[[1]]$macoef
  v <- z$model[[1]]$v
  y2 <- yy[301:350]
  prdctr(y2, r=30, s=50, h=10, arcoef=ar, macoef=ma, v=v)
}

\keyword{ts}
