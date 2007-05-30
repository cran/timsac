.packageName <- "timsac"

mfilter <-
function (x, filter, method = c("convolution", "recursive"), init = NULL) 
{
    method <- match.arg(method)
    x <- as.ts(x)
    xtsp <- tsp(x)
    x <- as.matrix(x)
    n <- nrow(x)
    nser <- ncol(x)
    p <- dim(filter)[3]
    nd <- n+p
    if (any(is.na(filter))) 
        stop("missing values in 'filter'")
    y <- matrix(0, nd, nser)
    if (method == "convolution") {
        if (p > n) 
            stop("'filter' is longer than time series")
        if (missing(init)) {
            init <- matrix(0, p, nser)
        }
        else {
            ni <- NROW(init)
            if (ni != p) 
                stop("length of 'init' must equal length of 'filter'")
            if (NCOL(init) != 1 && NCOL(init) != nser) 
                stop(gettextf("'init'; must have 1 or %d cols", 
                  nser), domain = NA)
            if (!is.matrix(init)) 
                init <- matrix(init, p, nser)
        }
        xx <- matrix(NA,nd,nser)
        for (i in 1:p) xx[i, ] <- init[i, ]
        tar <- array(0,dim=c(p,nser,nser))
        for (i in 1:p) tar[i, ,] <- t(filter[, ,i])
        i <- p + 1
        while (i <= nd) {
          xx[i,] = x[i-p,]
          y[i, ] <- x[i-p, ]
          for (j in 1:p) y[i, ] <- y[i, ] - xx[i - j, ] %*% tar[j, ,]
             i <- i + 1
          }
    }
    else {
        if (missing(init)) {
            init <- matrix(0, p, nser)
        }
        else {
            ni <- NROW(init)
            if (ni != p) 
                stop("length of 'init' must equal length of 'filter'")
            if (NCOL(init) != 1 && NCOL(init) != nser) 
                stop(gettextf("'init'; must have 1 or %d cols", 
                  nser), domain = NA)
            if (!is.matrix(init)) 
                init <- matrix(init, p, nser)
        }
        for (i in 1:p) y[i, ] <- init[i, ]
        tar <- array(0,dim=c(p,nser,nser))
        for (i in 1:p) tar[i, ,] <- t(filter[, ,i])
        i <- p + 1
        while (i <= nd) {
          y[i, ] <- x[i-p, ]
          for (j in 1:p) y[i, ] <- y[i, ] + 
                     y[i - j, ] %*% tar[j, ,]
             i <- i + 1
          }
     }
    y <- y[(p+1):nd, ]
    y <- drop(y)
    tsp(y) <- xtsp
    class(y) <- "mts"
    return(y)
}

#####   TIMSAC72   #####

autcor <- function (y, lag=NULL, plot=TRUE, lag_axis=TRUE)
{
    n <- length(y)
    if( is.null(lag) )  lag <- as.integer(2*sqrt(n))    # maximum lag
    lag1 <- lag+1
    acov <- rep(0,lag1)  # auto covariance
    acor <- rep(0,lag1)  # auto correlation
    mean <- 0            # mean of y

    z <- .C("autcor",
	     as.double(y),
	     as.integer(n),
	     acov = as.double(acov),
	     acor = as.double(acor),
	     as.integer(lag1),
	     mean = as.double(mean))

    if( plot == TRUE ) {
      plot((0:lag), z$acor, type="h", ylab="Autocorrelation", xlab="Lag")
      if( lag_axis == TRUE ) abline(h=0, lty=1)
    }

    autcor.out <- list( acov=z$acov, acor=z$acor, mean=z$mean )
    return( autcor.out )
}


fftcor <- function (y, lag=NULL, isw=4, plot=TRUE, lag_axis=TRUE)
{
    ld <- nrow(y)                # length of data
    d <- ncol(y)                 # dimension of the observation vector
    x1 <- y[,1]                  # data of channel X
    y1 <- rep(0,ld)
    if( d==2 ) y1 <- y[,2]       # data of channel Y
    if( is.null(lag) )  lag <- as.integer(2*sqrt(ld))    # maximum lag
    lag1 <- lag+1

#  n2p, n : definition    n=2**n2p >= ld+lag1 >= 2**(n2p-1)
    nd <- ld+lag1
    n2p <- 1
    n <- 2**n2p
    while( (n-nd) < 0 ) {
	n2p <- n2p+1
	n <- 2**n2p
    }

    acov <- array(0, dim=c(n,2))      # auto covariance
    ccov21 <- rep(0,n)                # cross covariance
    ccov12 <- rep(0,n)                # cross covariance
    acor <- array(0, dim=c(lag1,2))   # auto correlation
    ccor21 <- rep(0,lag1)             # cross correlation
    ccor12 <- rep(0,lag1)             # cross correlation
    mean <- rep(0,2)                  # mean

    z <- .C("fftcor",
	as.integer(ld),
	as.integer(lag1),
	as.integer(n),
	as.integer(n2p),
	as.integer(isw),
	as.double(x1),
	as.double(y1),
	acov = as.double(acov),
	ccov21 = as.double(ccov21),
	ccov12 = as.double(ccov12),
	acor = as.double(acor),
	ccor21 = as.double(ccor21),
	ccor12 = as.double(ccor12),
	mean = as.double(mean))
    
    acv <- array(z$acov, dim=c(n,2))
    acov <- array(, dim=c(lag1,2))
    for( i in 1:lag1 ) acov[i,1] <- acv[i,1]
    for( i in 1:lag1 ) acov[i,2] <- acv[i,2]
    
    if( isw == 4 ) {
      fftcor.out  <- list( acov=acov, ccov12=z$ccov12[1:lag1], ccov21=z$ccov21[1:lag1],
	                   acor=array(z$acor, dim=c(lag1,2)), ccor12=z$ccor12, ccor21=z$ccor21, mean=z$mean )
      if( plot == TRUE ) {
        par(mfrow=c(2,1))
        plot((0:lag), z$ccor12, type="h", xlab="Lag", ylab="Crosscorrelation ccor12")
        if( lag_axis == TRUE ) abline(h=0, lty=1)
        plot((0:lag), z$ccor21, type="h", xlab="Lag", ylab="Crosscorrelation ccor21")
        if( lag_axis == TRUE ) abline(h=0, lty=1) }
        par(mfrow=c(1,1))
    }

    if( isw != 4 ) fftcor.out  <- list( acov=acov, acor=array(z$acor, dim=c(lag1,2)), mean=z$mean )

    return( fftcor.out )
}


mulcor <-
function (y, lag=NULL, plot=TRUE, lag_axis=TRUE)
{
    n <- nrow(y)        # length of data
    d <- ncol(y)        # dimension of the observation vector
    if( is.null(lag) )  lag <- as.integer(2*sqrt(n))  # maximum lag
    lag1 <- lag+1
    mean <- rep(0,d)    # mean
    cov <- array(0,dim=c(lag1,d,d))   # covariance
    cor <- array(0,dim=c(lag1,d,d))   # normalized covariance

    z <- .C("mulcor",
	as.double(y),
	as.integer(n),
	as.integer(d),
	as.integer(lag1),
        mean = as.double(mean),
	cov = as.double(cov),
	cor = as.double(cor))

    cor=array(z$cor,dim=c(lag1,d,d))
    if( plot == TRUE ) {
      x <- rep(0,lag1)
      for( i in 1:lag1 ) x[i] <- i-1
      par(mfrow=c(d,d))
      for( j in 1:d )
        for( k in 1:d )  {
           cy <- cor[,j,k]
           plot(x, cy, type="h", xlab="lag", ylab=paste("cor (",j,",",k,")"),ylim=c(-1,1))
           if( lag_axis == TRUE ) abline(h=0, lty=1)
        }
      par(mfrow=c(1,1))
      }

    mulcor.out <- list( cov=array(z$cov,dim=c(lag1,d,d)), cor=array(z$cor,dim=c(lag1,d,d)), mean=z$mean )
    class( mulcor.out ) <- "mulcor"
    return( mulcor.out )
}

print.mulcor <- function(x, ...)
{
  lag <- dim(x$cov)[1]-1
  ip <- dim(x$cov)[2]

  for( i in 1:ip ) {
    cat(sprintf("\nAUTOCOVARIANCE \t\tMEAN=%f\n", x$mean[i]))
    cat(sprintf("L\t C(%i , %i)\tNORMALIZED\n",i,i))
    for( j in 0:lag ) cat(sprintf("%i\t%f\t%f\n", j,x$cov[j+1,i,i],x$cor[j+1,i,i]))
    if( i != 1 ) for( j in 1:(i-1) ) {
      cat("\nCROSS COVARIANCE\n")
      cat(sprintf("L\t C(%i , %i)\tNORMALIZED\t C(%i, %i)\tNORMALIZED\n",j,i,i,j))
      for( k in 0:lag ) cat(sprintf("%i\t%f\t%f\t%f\t%f\n", k,x$cov[k+1,j,i],x$cor[k+1,j,i],x$cov[k+1,i,j],x$cor[k+1,i,j]))
    }
  }
}

auspec <-
function (y, lag=NULL, window="Akaike", log=FALSE, plot=TRUE)
{
    n <- length(y)
    if( is.null(lag) ) lag <- as.integer(2*sqrt(n))  # maximum lag
    lag1 <- lag+1

    spec1 <- rep(0,lag1)   # spectrum smoothing by window W1
    spec2 <- rep(0,lag1)   # spectrum smoothing by window W2
    stat <- rep(0,lag1)    # test statistics

    z1 <- autcor(y, lag, plot=FALSE)

    z <- .C("auspec",
	     as.integer(n),
	     as.integer(lag1),
	     as.double(z1$acov),
	     spec1 = as.double(spec1),
	     spec2 = as.double(spec2),
	     stat = as.double(stat))

    if( window == "Akaike" ) spec <- z$spec2
    if( window == "Hanning" ) spec <- z$spec1

    if( plot == TRUE ) {
      x <- rep(0,lag1)
      for( i in 1:lag1 ) x[i] <- (i-1)/(2*lag)
      if ( log == TRUE ) plot(x, spec, type="l", log="y", xlab="Frequency", ylab="Spectrum", )
      if ( log == FALSE ) plot(x, spec, type="l", xlab="Frequency", ylab="Spectrum", ) }

    auspec.out <- list( spec=spec, stat=z$stat )
    return( auspec.out )
}


mulspe <- 
function (y, lag = NULL, window = "Akaike", plot = TRUE, plot.scale = FALSE) 
{
    n <- nrow(y)
    d <- ncol(y)
    if (is.null(lag)) 
        lag <- as.integer(2 * sqrt(n))
    lag1 <- lag + 1
    spec1 <- array(0, dim = c(lag1, d, d))
    spec2 <- array(0, dim = c(lag1, d, d))
    stat <- array(0, dim = c(lag1, d))
    coh1 <- array(0, dim = c(lag1, d, d))
    coh2 <- array(0, dim = c(lag1, d, d))
    z1 <- mulcor(y, lag, plot = FALSE)
    z <- .C("mulspe", as.integer(n), as.integer(d), as.integer(lag1), 
        as.integer(lag1), as.double(z1$cov), spec1 = as.double(spec1), 
        spec2 = as.double(spec2), stat = as.double(stat), coh1 = as.double(coh1), 
        coh2 = as.double(coh2))
    if (window == "Akaike") {
        spec <- array(z$spec2, dim = c(lag1, d, d))
        coh <- array(z$coh2, dim = c(lag1, d, d))
    }
    if (window == "Hanning") {
        spec <- array(z$spec1, dim = c(lag1, d, d))
        coh <- array(z$coh1, dim = c(lag1, d, d))
    }
    cspec <- array(0, dim = c(d, d, lag1))
    if (plot == TRUE) {
        for (j in 1:d) {
            for (k in 1:j) {
                for (l in 1:lag1) {
                  cspec[j, k, l] <- spec[l, j, k] + (0+1i) * spec[l, k, j]
                  cspec[k, j, l] <- spec[l, j, k] - (0+1i) * spec[l, k, j]
                }
            }
        }
        plot.mulspec(cspec, d, lag, plot.scale)
    }
    mulspe.out <- list(spec = spec, stat = array(z$stat, dim = c(lag1, d)), coh = coh)
    class(mulspe.out) <- "mulspe"
    return(mulspe.out)
}


print.mulspe <- function(x, ...)
{
  lag <- dim(x$spec)[1]-1
  ip <- dim(x$spec)[2]

  for( i in 1:ip ) {
    cat(sprintf("\nPOWER SPECTRUM P(%i,%i)\tSIGNIFICANCE\n",i,i))
    for( j in 0:lag ) cat(sprintf("%i\t%f\t%f\n", j,x$spec[j+1,i,i],x$stat[j+1,i]))
    if( i != 1 ) for( j in 1:(i-1) ) {
      cat(sprintf("\nCROSS SPECTRUM P(%i,%i)\n",i,j))
      cat("I\t CO-SPECTRUM\tQUAD-SPECTRUM\tSIMPLE COHERENCE\n")
      for( k in 0:lag ) cat(sprintf("%i\t%f\t%f\t%f\n", k,x$spec[k+1,i,j],x$spec[k+1,j,i],x$coh[k+1,i,j]))
    }
  }
}


plot.mulspec <- function(spec, d, lag, plot.scale)
{
  par(mfrow = c(d,d))
  lag1 <- lag+1
  x <- rep(0,lag1)
  for (i in 1:lag1) x[i] <- (i-1)/(2*lag)
  dspec <- array(0, dim=c(d,d,lag1))
  for (j in 1:d) for (k in 1:d) {
    if ( j >= k ) {
      dspec[j,k,] = Mod(spec[j,k,])
    } else {
      dspec[j,k,] = Arg(spec[j,k,])
    }
  }

  if(plot.scale == TRUE) {
    mx  <- max(dspec[1,1,])
    for(j in 1:d) for (k in 1:j) if (mx < max(dspec[j,k,])) mx <- max(dspec[j,k,])
    for(j in 1:d) for (k in 1:d)
      if ( j >= k ) {
      ylabs = paste("AmpSp (", j, ",", k, ")")
      plot(x, dspec[j, k, ], type = "l", xlab = "Frequency", ylab = ylabs, ylim=c(0,mx))
    } else {
      ylabs = paste("PhaseSp (", j, ",", k, ")")
      plot(x, dspec[j, k, ], type = "l", xlab = "Frequency", ylab = ylabs, ylim=c(-pi,pi))
    }
  } else {
    for(j in 1:d) for (k in 1:d)
      if ( j >= k ) {
      ylabs = paste("AmpSp (", j, ",", k, ")")
      plot(x, dspec[j, k, ], type = "l", xlab = "Frequency", ylab = ylabs)
    } else {
      ylabs = paste("PhaseSp (", j, ",", k, ")")
      plot(x, dspec[j, k, ], type = "l", xlab = "Frequency", ylab = ylabs)
    }
  }
  par(mfrow=c(1,1))
}


sglfre <-
function (y, lag=NULL, invar, outvar)
{
    n <- nrow(y)        # length of data
    d <- ncol(y)        # dimension of the observation vector
    if( is.null(lag) ) lag <- as.integer(2*sqrt(n))
    lag1 <- lag+1

    spec1 <- rep(0,lag1)   # power spectrum
    spec2 <- rep(0,lag1)   # power spectrum
    cspec <- rep(0,lag1)   # co-spectrum
    qspec <- rep(0,lag1)   # quad-spectrum
    gain <- rep(0,lag1)    # gain
    coh <- rep(0,lag1)     # coherncy
    freqr <- rep(0,lag1)   # frequency response function : real part
    freqi <- rep(0,lag1)   # frequency response function : imaginary part
    errstat <- rep(0,lag1) # relative error statistics
    phase <- rep(0,lag1)   # phase

    z1 <- mulspe(y, lag, window="Akaike", plot=FALSE)

    z <- .C("sglfre",
	     as.integer(invar),
	     as.integer(outvar),
	     as.integer(n),
	     as.integer(lag1),
	     as.integer(d),
	     as.double(z1$spec),
	     spec1 = as.double(spec1),
	     spec2 = as.double(spec2),
	     cspec = as.double(cspec),
	     qspec = as.double(qspec),
	     gain = as.double(gain),
	     coh = as.double(coh),
	     freqr = as.double(freqr),
	     freqi = as.double(freqi),
	     errstat = as.double(errstat),
	     phase = as.double(phase))

    sglfre.out <- list( inspec=z$spec1, outspec=z$spec2, cspec=z$cspec, qspec=z$qspec, gain=z$gain, coh=z$coh,
			freqr=z$freqr, freqi=z$freqi, errstat=z$errstat, phase=z$phase )
    return( sglfre.out )
}


mulfrf <-
function (y, lag=NULL, niv, iovar=c(1:(niv+1)))
{
    n <- nrow(y)        # length of data
    d <- ncol(y)        # dimension of the observation vector
    if( is.null(lag) ) lag <- as.integer(2*sqrt(n))  # maximum lag
    lag1 <- lag+1
    nv <- niv           # number of input variables

    cospec <- array(complex(real=0,imag=0),dim=c(d,d,lag1))  # spectrum
    freqr <- array(0,dim=c(nv,lag1))     # frequency response function : real part
    freqi <- array(0,dim=c(nv,lag1))     #                             : imag. part
    gain <- array(0,dim=c(nv,lag1))      # gain
    phase <- array(0,dim=c(nv,lag1))     # phase
    pcoh <- array(0,dim=c(nv,lag1))      # partial coherency
    errstat <- array(0,dim=c(nv,lag1))   # relative error statistics
    mcoh <- rep(0,lag1)                  # multiple coherency

    z1 <- mulspe(y, lag, plot=FALSE)

    z <- .C("mulfrf",
	as.integer(nv),
	as.integer(iovar),
	as.integer(n),
	as.integer(lag1),
	as.integer(d),
	as.double(z1$spec),
	cospec = as.complex(cospec),
	freqr = as.double(freqr),
	freqi = as.double(freqi),
	gain = as.double(gain),
	phase = as.double(phase),
	pcoh = as.double(pcoh),
	errstat = as.double(errstat),
	mcoh = as.double(mcoh))

    mulfrf.out <- list( cospec=array(z$cospec, dim=c(d,d,lag1)), freqr=array(z$freqr, dim=c(nv,lag1)),
			freqi=array(z$freqi, dim=c(nv,lag1)), gain=array(z$gain, dim=c(nv,lag1)),
			phase=array(z$phase, dim=c(nv,lag1)), pcoh=array(z$pcoh, dim=c(nv,lag1)),
			errstat=array(z$errstat, dim=c(nv,lag1)), mcoh=z$mcoh )
    return( mulfrf.out )
}


fpeaut <-
function (y, max.order=NULL)
{
    n <- length(y)
    if( is.null(max.order) ) max.order <- as.integer(2*sqrt(n))   # upper limit of model order
    morder <- max.order
    sigma2 <- rep(0,morder)            # sigma**2
    fpe <- rep(0,morder)               # FPE (Final Prediction Error)
    rfpe <- rep(0,morder)              # RFPE
    parcor <- rep(0,morder)            # parcor
    chi2 <- rep(0,morder)              # chi-squared
    ofpe <- 0                          # ofpe
    fpemin <- 0                        # minimum fpe
    rfpemin <- 0                       # minimum rfpe
    ordermin <- 0                      # order of minimum fpe
    sigma2m <- 0                       # sigma**2 attained at m=order.fpemin
    a <- array(0,dim=c(morder,morder)) # coefficients of ar-process
    ao <- rep(0,morder)                # coefficients of ar-process attained at m=order.fpemin

    lag <- morder
    lag1 <- lag+1
    z1 <- autcor(y, lag, plot=FALSE)

    z <- .C("fpeaut",
	as.integer(morder),
	as.integer(n),
	as.double(z1$acov[1]),
	as.double(z1$acov[2:(lag1)]),
	sigma2 = as.double(sigma2),
	fpe = as.double(fpe),
	rfpe = as.double(rfpe),
	parcor = as.double(parcor),
	chi2 = as.double(chi2),
	ofpe = as.double(ofpe),
	fpemin = as.double(fpemin),
	rfpemin = as.double(rfpemin),
	ordermin = as.integer(ordermin),
	sigma2m = as.double(sigma2m),
	a = as.double(a),
	ao = as.double(ao))

    a <- array(z$a,dim=c(morder,morder))
    arcoef <- list()
    for( i in 1:morder ) arcoef[[i]] <- a[1:i,i]

    fpeaut.out <- list( ordermin=z$ordermin, best.ar=arcoef[[z$ordermin]], sigma2m=z$sigma2m, fpemin=z$fpemin,
                        rfpemin=z$rfpemin, ofpe=z$ofpe, arcoef=arcoef, sigma2=z$sigma2, fpe=z$fpe, rfpe=z$rfpe,
			parcor=z$parcor, chi2=z$chi2 )
    return( fpeaut.out )
}


fpec <-
function (y, max.order=NULL, ncon=NULL, nman=0, inw=NULL)
{
    n <- nrow(y)        # length of data
    d <- ncol(y)        # dimension of the observation vector
    if( is.null(max.order) ) max.order <- as.integer(2*sqrt(n))  # upper limit of model order
    if( is.null(ncon) ) ncon <- d                    # number of controlled variables
    if( is.null(nman) ) nman <- 0                    # number of manipulated variables
    ip <- ncon+nman
    if( is.null(inw) ) inw <- c(1:ip)                # indicator

    morder <- max.order
    morder1 <- morder+1
    ccv <- array(0,dim=c(morder1,ip,ip))      # matrix rearrangement by inw
    fpec <- rep(0,morder1)                    # FPEC (AR-MODEL FITTING FOR CONTROL)
    rfpec <- rep(0,morder1)                   # RFPEC
    aic <- rep(0,morder1)                     # aic
    ordermin <- 0                             # order of minimum fpec
    fpecmin <- 0                              # minimum fpec
    rfpecmin <- 0                             # minimum rfpec
    aicmin <- 0                               # minimum aic
    perr <- array(0,dim=c(ncon,ncon))         # prediction error covariance matrix
    arcoef <- array(0,dim=c(morder,ncon,ip))  # the set of coefficient matrices

    morder1 <- morder+1
    z1 <- mulcor(y, morder, plot=FALSE)
    cov <- z1$cov[1:morder1,,]

    z <- .C("fpec7",
	as.integer(n),
	as.integer(morder),
	as.integer(ncon),
	as.integer(ip),
	as.integer(d),
	as.integer(inw),
	as.double(cov),
	ccv = as.double(ccv),
	fpec = as.double(fpec),
	rfpec = as.double(rfpec),
	aic = as.double(aic),
	ordermin = as.integer(ordermin),
	fpecmin = as.double(fpecmin),
	rfpecmin = as.double(rfpecmin),
	aicmin = as.double(aicmin),
	perr = as.double(perr),
	arcoef = as.double(arcoef))

    cov <- array(z$ccv, dim=c(morder1,ip,ip))
    cov <- aperm(cov, c(2,3,1))
    arcoef <- array(z$arcoef,dim=c(morder,ncon,ip))
    arcoef <- arcoef[1:z$ordermin,,]
    arcoef <- aperm(arcoef, c(2,3,1))

    fpec.out  <- list( cov=cov, fpec=z$fpec, rfpec=z$rfpec, aic=z$aic, ordermin=z$ordermin,
                       fpecmin=z$fpecmin, rfpecmin=z$rfpecmin, aicmin=z$aicmin,
                        perr=array(z$perr,dim=c(ncon,ncon)), arcoef=arcoef )
    class( fpec.out ) <- "fpec"
    return( fpec.out )
}

print.fpec <- function(x, ...)
{
  m1 <- dim(x$arcoef)[1]
  m2 <- dim(x$arcoef)[2]
  m3 <- dim(x$arcoef)[3]
  lag <- dim(x$cov)[3]-1

  cat("$cov\nCOVARIANCE MATRIX\n\n")
  print(x$cov)

  cat("M\tFPEC\t\tRFPEC\t\tAIC\n")
  for( i in 0:lag ) cat(sprintf("%i\t%f\t%f\t%f\n",i, x$fpec[i+1],x$rfpec[i+1],x$aic[i+1]))
  cat(sprintf("\n  MINIMUM FPEC = %f\tMINIMUM RFPEC = %f  ATTANED AT M = %i\n", x$fpecmin, x$rfpecmin,  x$ordermin))
  cat(sprintf("  MINIMUM AIC = %f\n", x$aicmin))

  cat("\n\n$perr\nPREDICTION ERROR COVARIANCE MATRIX\n")
  print(x$perr)

  cat("\n\n$arcoef\nAR COEFFICIENT MATRIX\n\n")
  print(x$arcoef)

}

mulnos <-
function ( y, max.order=NULL, ncon=NULL, nman=0, h, inw=NULL )
{
    n <- nrow(y)        # length of data
    d <- ncol(y)        # dimension of the observation vector
    if( is.null(max.order) ) max.order <- as.integer(2*sqrt(n))  # upper limit of model order
    if( is.null(ncon) ) ncon <- d           # number of controlled variables
    if( is.null(nman) ) nman <- 0           # number of manipulated variables
    ip <- ncon+nman
    if( is.null(inw) ) inw <- c(1:ip)       # indicator

    nperr <- array(0,dim=c(ip,ip))          # normalized prediction error covaiance matrix
    diffr <- array(0,dim=c(ip,ip,(h+1)))    # differential relative power contribution
    integr <- array(0,dim=c(ip,ip,(h+1)))   # integrated relative power contribution

    z1 <- fpec(y, max.order, ncon, nman, inw)
    arcoef <- aperm(z1$arcoef, c(3,1,2))

    z <- .C("mulnos",
	as.integer(h),
	as.integer(z1$ordermin),
	as.integer(ip),
	as.double(z1$perr),
	as.double(arcoef),
	nperr = as.double(nperr),
	diffr = as.double(diffr),
	integr = as.double(integr))

    mulnos.out <- list( nperr=array(z$nperr,dim=c(ip,ip)), diffr=array(z$diffr,dim=c(ip,ip,(h+1))),
			integr=array(z$integr,dim=c(ip,ip,(h+1))) )
    return( mulnos.out )
}


raspec <-
function (h, var, arcoef=NULL, macoef=NULL, log = FALSE, plot=TRUE)
{
    l <- length(arcoef)
    k <- length(macoef)
    if( is.null(arcoef) )  arcoef <- 0   # coefficient matrix of autoregressive model
    if( is.null(macoef) )  macoef <- 0   # coefficient matrix of moving average model
    macoef <- -macoef
    rspec <- rep(0,(h+1))                # rational spectrum

    z <- .C("raspec",
	as.integer(h),
	as.integer(l),
	as.integer(k),
	as.double(var),
	as.double(arcoef),
	as.double(macoef),
	rspec = as.double(rspec))

    if( plot == TRUE ) {
      x <- rep(0,(h+1))
      for( i in 1:(h+1) ) x[i] <- (i-1)/(2*h)
      if( log == TRUE ) plot(x, z$rspec, type="l", log="y", xlab="Frequency", ylab="Rational Spectrum")
      if( log == FALSE ) plot(x, z$rspec, type="l", xlab="Frequency", ylab="Rational Spectrum") }

    return( raspec=z$rspec )
}

mulrsp <-
function (h, d, cov, ar = NULL, ma = NULL, log = FALSE, plot = TRUE, 
    plot.scale = FALSE) 
{
    if (is.null(ar)) {
        l <- 0
        ar <- 0
    }
    else {
        l <- dim(ar)[3]
    }
    if (is.null(ma)) {
        k <- 0
        ma <- 0
    }
    else {
        k <- dim(ma)[3]
    }
    arcoef <- array(0, dim = c(d, d, l))
    macoef <- array(0, dim = c(d, d, k))
    if (l != 0) 
        arcoef <- aperm(ar, c(3, 1, 2))
    if (k != 0) 
        macoef <- aperm(ma, c(3, 1, 2))
    macoef <- -macoef
    rspec <- array(0, dim = c(d, d, (h + 1)))
    spec <- array(0, dim = c(d, d, (h + 1)))
    scoh <- array(0, dim = c(d, d, (h + 1)))
    z <- .C("mulrsp", as.integer(h), as.integer(l), as.integer(d), 
        as.integer(k), as.double(cov), as.double(arcoef), as.double(macoef), 
        rspec = as.complex(rspec), scoh = as.double(scoh))
    mulrsp.out <- list(rspec = array(z$rspec, dim = c(d, d, (h + 
        1))), scoh = array(z$scoh, dim = c(d, d, (h + 1))))
    if (plot == TRUE) 
        plot.mulspec(mulrsp.out$rspec, d, h, plot.scale)
    return(mulrsp.out)
}


optdes <-
function (y, max.order=NULL, ns, q, r)
{
    n <- nrow(y)       # length of data
    d <- ncol(y)       # dimension of the observation vector
    if( is.null(max.order) ) max.order <- as.integer(2*sqrt(n))  # upper limit of model order
    ncon <- dim(q)[1]  # number of controlled variables
    nman <- dim(r)[1]  # number of manipulated variables

    z1 <- fpec(y, max.order, ncon, nman)
    order <- z1$ordermin
    ao <- z1$arcoef
    ao <- aperm(ao, c(3,1,2))
    ao <- aperm(ao, c(2,1,3))
    dim(ao) <- c(ncon*order,(ncon+nman))
    a <- ao[,(1:ncon)]
    b <- ao[,(ncon+1):(ncon+nman)]

    gain <- array(0,c(nman,order*ncon))   # gain

    z <- .C("optdes",
	as.integer(ncon),
	as.integer(nman),
	as.integer(ns),
	as.integer(order),
	as.double(q),
	as.double(r),
	as.double(z1$perr),
	as.double(a),
	as.double(b),
	gain = as.double(gain))

    optdes.out <- list( perr=z1$perr, trans=a, gamma=b, gain=array(z$gain, dim=c(nman,order*ncon)) )
    return( optdes.out )
}


optsim <-
function (y, max.order=NULL, ns, q, r, noise=NULL, len, plot=TRUE)
{
    n <- nrow(y)        # length of data
    d <- ncol(y)        # dimension of the observation vector
    if( is.null(max.order) ) max.order <- as.integer(2*sqrt(n))  # upper limit of model order
    ir <- dim(q)[1]     # number of controlled variables
    il <- dim(r)[1]     # number of manipulated variables

    xx <- array(0,dim=c(ir,ns))  # controlled variables X
    yy <- array(0,dim=c(il,ns))  # manipulated variables Y
    xmean <- rep(0,ir)       # mean of X
    ymean <- rep(0,il)       # mean of Y
    x2sum <- rep(0,ir)       # sum of X**2
    y2sum <- rep(0,il)       # sum of Y**2
    x2mean <- rep(0,ir)      # mean of X**2
    y2mean <- rep(0,il)      # mean of Y**2
    xvar <- rep(0,ir)        # variance of X
    yvar <- rep(0,il)        # variance of Y

    z1 <- optdes(y, max.order, ns, q, r)
    trans <- z1$trans
    gamma <- z1$gamma
    gain <- z1$gain
    order <- dim(gain)[2]/ir
    if( is.null(noise) ) noise <- wnoise(len, z1$perr, plot=FALSE)     # white noise

    z <- .C("optsim",
	as.integer(ns),
	as.integer(order),
	as.integer(ir),
	as.integer(il),
	as.double(trans),
	as.double(gamma),
	as.double(gain),
	as.double(noise),
	xx = as.double(xx),
	yy = as.double(yy),
	xmean = as.double(xmean),
	ymean = as.double(ymean),
	x2sum = as.double(x2sum),
	y2sum = as.double(y2sum),
	x2mean = as.double(x2mean),
	y2mean = as.double(y2mean),
	xvar = as.double(xvar),
	yvar = as.double(yvar))

    convar <- array(z$xx,dim=c(ir,ns))
    manvar <- array(z$yy,dim=c(il,ns))

    if( plot == TRUE ) {
      nc <- max(ir,il)
      par(mfcol=c(2,nc))
      for( i in 1:ir ) {
        plot(convar[i,], type="h", xlab="Step",ylab=paste("Controlled variables X(",i,")"))
        abline(h=0, lty=1) }
      for( i in 1:il ) {
        plot(manvar[i,], type="h", xlab="Step",ylab=paste("Manipulated variables Y(",i,")"))
        abline(h=0, lty=1) }
      par(mfrow=c(1,1))
    }

    optsim.out <- list( trans=trans, gamma=gamma, gain=gain, convar=convar, manvar=manvar,
			xmean=z$xmean, ymean=z$ymean, xvar=z$xvar, yvar=z$yvar,
			x2sum=z$x2sum, y2sum=z$y2sum, x2mean=z$x2mean, y2mean=z$y2mean )
    return( optsim.out )
}


wnoise <-
function (len, perr, plot=TRUE)
{
    ir <- 1                           # number of controlled variables
    if( is.array(perr) ) ir <- ncol(perr)
    w <- array(0,dim=c(ir,len))       # white noise

    z <- .C("wnoise",
	as.integer(len),
	as.integer(ir),
	as.double(perr),
	w = as.double(w))

    if( ir != 1 ) wnoise <- array(z$w,dim=c(ir,len))
    if( ir == 1 ) wnoise <- z$w[1:len]

    if( plot == TRUE ) {
#      nc <- as.integer((ir+1)/2)
#      if( ir != 1 ) par(mfrow=c(2,nc))
      par(mfcol=c(ir,1))
      if( ir == 1 ) {
        plot(wnoise, type="h", ylab="White Noise")
        abline(h=0, lty=1)
      } else {
        for( i in 1:ir ) {
          plot(wnoise[i,], type="h", ylab=paste("White Noise (",i,")"))
          abline(h=0, lty=1)
        }
      }
      par(mfrow=c(1,1))
    }

    return( wnoise=wnoise )
}

#####   TIMSAC74   #####

autoarmafit <-
function (y, max.order=NULL, tmp.file=NULL)
{
    n <- length(y)
    if( is.null(max.order) ) max.order <- as.integer(2*sqrt(n))
    if( is.null(tmp.file) )  tmp.file <- " "

    morder <- max.order
    autcv <- autcor(y, morder, plot=FALSE)$acov   # covariance sequence
    z1 <- canarm(y, morder, plot=FALSE)
    inc <- 1               # total number of case
    arcoef1 <- z1$arcoef   # initial estimates of AR-coefficients
    macoef1 <- z1$macoef   # initial estimates of MA-coefficients
    p1 <- length(arcoef1)  # initial AR order
    q1 <- length(macoef1)  # initial MA order

    icst <- 190
    mmax <- morder
    lmax <- morder+icst+mmax
    nmax <- 25

    newn <- 0                          # number of cases
    p <- rep(0,nmax)                   # AR order
    a <- array(0, dim=c(mmax,nmax))    # maximum lokelihood estimates of AR coefficients
    q <- rep(0,nmax)                   # MA order
    b <- array(0, dim=c(mmax,nmax))    # maximum lokelihood estimates of MA coefficients
    std <- array(0, dim=c(mmax,nmax))  # standard deviation
    v <- rep(0,nmax)                   # innovation variance
    gr <- array(0, dim=c(mmax,nmax))   # final gradient
    aic <- rep(0,nmax)                 # n*log(v)+2*(p+q)
    aicm <- rep(0,nmax)                # minimum AIC
    pbest <- 0                         # AR order of best choice
    qbest <- 0                         # MA order of best choice

    z <- .C("autarm",
	as.integer(n),
	as.integer(morder+1),
	as.double(autcv),
	as.integer(inc),
	as.integer(p1),
	as.double(arcoef1),
	as.integer(q1),
	as.double(macoef1),
	newn = as.integer(newn),
	p = as.integer(p),
	a = as.double(a),
	q = as.integer(q),
	b = as.double(b),
	std = as.double(std),
	v = as.double(v),
        gr = as.double(gr),
	aic = as.double(aic),
        aicm = as.double(aicm),
	qbest = as.integer(qbest),
	pbest = as.integer(pbest),
        as.character(tmp.file),
	as.integer(lmax),
	as.integer(mmax),
	as.integer(nmax) )

    nc <- z$newn
    p <- z$p[1:nc]
    q <- z$q[1:nc]
    a <- array(z$a, dim=c(mmax,nmax))
    b <- array(z$b, dim=c(mmax,nmax))
    std <- array(z$std, dim=c(mmax,nmax))
    v <- rep(0,nc)
    gr <- array(z$gr, dim=c(mmax,nmax))
    aaic <- z$aic[1:nc]
    arcoef <- list()
    macoef <- list()
    arstd <- list()
    mastd <- list()
    grad <- list()

    aicorder <- sort(aaic, index.return=TRUE)

    for ( i in 1:nc ) {
        j <- aicorder$ix[i]
        arcoef[[i]] <- -a[(1:p[j]),j]
        macoef[[i]] <- -b[(1:q[j]),j]
        arstd[[i]] <- std[(q[j]+1):(p[j]+q[j]),j]
        mastd[[i]] <- std[1:q[j],j]
        v[i] <- z$v[j]
        grad[[i]] <- gr[(1:(p[j]+q[j])),j]
    }

    best.model <- list( ar=arcoef[[1]], ma=macoef[[1]] )
    best.order <- list( arorder=z$qbest, maorder=z$pbest )

    model <- list()
    for ( i in 1:nc )
      model[[i]] <- list( arcoef=arcoef[[i]], macoef=macoef[[i]], arstd=arstd[[i]], mastd=mastd[[i]], v=v[i],
                 aic=aicorder$x[i], grad=grad[[i]] )

    autoarmafit.out <- list( best.order=best.order, best.model=best.model, model=model )
    class( autoarmafit.out ) <- "autoarmafit"
    return( autoarmafit.out )
}

print.autoarmafit <- function(x, ...)
{
  cat("\nBest ARMA model")
  cat(sprintf("\n AR coefficient (order =%i)",x$best.order$arorder))
  for( j in 1:x$best.order$arorder ) cat(sprintf("\t%f",x$best.model$ar[j]))
  cat(sprintf("\n MA coefficient (order =%i)",x$best.order$maorder))
  for( j in 1:x$best.order$maorder ) cat(sprintf("\t%f",x$best.model$ma[j]))

  nc <- length(x$model)
  for( i in 1:nc) {
    cat(sprintf("\n\nCase No. %i\n", i))
    cat("\nI\tAR(I)\tSTANDARD DEVIATION\n")
    arorder <- length(x$model[[i]]$arcoef)
    for( j in 1:arorder )
      cat(sprintf("%i\t%f\t%f\n", j,x$model[[i]]$arcoef[j],x$model[[i]]$arstd[j]))
    cat("\nI\tMA(I)\tSTANDARD DEVIATION\n")
    maorder <- length(x$model[[i]]$macoef)
    for( j in 1:maorder )
      cat(sprintf("%i\t%f\t%f\n", j,x$model[[i]]$macoef[j],x$model[[i]]$mastd[j]))
    cat(sprintf("\nAIC\t%f\n", x$model[[i]]$aic))
    cat(sprintf("Innovation variance\t%f\n",x$model[[i]]$v))
    cat("Final gradient")
    for( j in 1:(arorder+maorder) ) cat(sprintf("\t%e", x$model[[i]]$grad[j]))
    cat("\n")
  }
}

armafit <-
function (y, model.order, tmp.file=NULL)
{
    n <- length(y)
    lag <- as.integer(2*sqrt(n))     # maximum lag
    lag1 <- lag+1

    inc <- 1             # total number of case
    p1 <- rep(0,2)
    q1 <- rep(0,2)

    autcv <- autcor(y, lag, plot=FALSE)$acov		# covariance sequence
    z1 <- canarm(y, lag, plot=FALSE)
    arcoef1 <- z1$arcoef	# initial estimates of AR-coefficients
    macoef1 <- z1$macoef	# initial estimates of MA-coefficients
    p1[1] <- length(arcoef1)	# initial AR order
    q1[1] <- length(macoef1)	# initial MA order

    if( length(model.order) == 2 ) {
      inc <- inc+1
      p1[2] <- model.order[1]     # AR order to be fitted successively
      q1[2] <- model.order[2]     # MA order to be fitted successively
    }

    icst <- 190
    mmax <- 50
    lmax <- lag+icst+mmax
    nmax <- 25	         # the limit of the total cases

    newn <- 0                          # number of cases
    p <- rep(0,nmax)                   # AR order
    a <- array(0, dim=c(mmax,nmax))    # maximum lokelihood estimates of AR coefficients
    q <- rep(0,nmax)                   # MA order
    b <- array(0, dim=c(mmax,nmax))    # maximum lokelihood estimates of MA coefficients
    std <- array(0, dim=c(mmax,nmax))  # standard deviation
    v <- rep(0,nmax)                   # innovation variance
    gr <- array(0, dim=c(mmax,nmax))   # final gradient
    aic <- rep(0,nmax)                 # n*log(cxx)+2*(l+k)
    aicm <- rep(0,nmax)                # minimum AIC
    pbest <- 0                         # AR order of best choice
    qbest <- 0                         # MA order of best choice
    if(is.null(tmp.file))  tmp.file <- " "
    
    z <- .C("autarm",
	as.integer(n),
	as.integer(lag1),
	as.double(autcv),
	as.integer(inc),
	as.integer(p1),
	as.double(arcoef1),
	as.integer(q1),
	as.double(macoef1),
	newn = as.integer(newn),
	p = as.integer(p),
	a = as.double(a),
	q = as.integer(q),
	b = as.double(b),
	std = as.double(std),
	v = as.double(v),
        gr = as.double(gr),
	aic = as.double(aic),
        aicm = as.double(aicm),
	pbest = as.integer(pbest),
	qbest = as.integer(qbest),
        as.character(tmp.file),
	as.integer(lmax),
	as.integer(mmax),
	as.integer(nmax) )

    nc <- z$newn
    p <- z$p[1:nc]
    q <- z$q[1:nc]
    a <- array(z$a, dim=c(mmax,nmax))
    b <- array(z$b, dim=c(mmax,nmax))
    std <- array(z$std, dim=c(mmax,nmax))
    gr <- array(z$gr, dim=c(mmax,nmax))
    aaic <- z$aic[1:nc]
    aic <- NULL
    arcoef <- NULL
    macoef <- NULL
    arstd <- NULL
    mastd <- NULL
    grad <- NULL

    aicorder <- sort(aaic, index.return=TRUE)
    order.maice <- aicorder$ix[1]

    arorder <- model.order[1]
    maorder <- model.order[2]

    for( i in 1:nc ) {
        if ( is.null(aic) == FALSE ) break
        j <- aicorder$ix[i]
        if( p[j]==arorder  &&  q[j]==maorder) {
          arcoef <- -a[(1:p[j]),j]
          macoef <- -b[(1:q[j]),j]
          arstd <- std[(q[j]+1):(p[j]+q[j]),j]
          mastd <- std[1:q[j],j]
          v <- z$v[j]
          aic <- z$aic[j]
          grad <- gr[(1:(p[j]+q[j])),j]
        }
    }

    armafit.out <- list( arcoef=arcoef, macoef=macoef, arstd=arstd, mastd=mastd, v=v, aic=aic, grad=grad )
    return( armafit.out )
}


bispec <-
function (y, lag=NULL, window="Akaike", log=FALSE, plot=TRUE)
{
    n <- length(y)
    if( is.null(lag) ) lag <- as.integer(2*sqrt(n))  # maximum lag
    lag1 <- lag+1

    z1 <- thirmo(y, lag, plot=FALSE)
    cv <- z1$acov               # autocovariance

    tmnt <- array(0,dim=c(lag1,lag1))
    for( i in 1:lag1 ) tmnt[i,1:i] <- z1$tmomnt[[i]]   # third order moments

    pspec1 <- rep(0,lag1)	# power spectrum smoothed by window w1
    pspec2 <- rep(0,lag1)	# power spectrum smoothed by window w2
    sig <- rep(0,lag1)		# significance
    ch <- array(0, dim=c(lag1,lag1))	# coherence
    br <- array(0, dim=c(lag1,lag1))	# real part of bispectrum
    bi <- array(0, dim=c(lag1,lag1))	# imaginary part of bispectrum
    rat <- 0

    z <- .C("bispec",
	as.integer(n),
	as.integer(lag),
	as.double(cv),
	as.double(tmnt),
	pspec1 = as.double(pspec1),
	pspec2 = as.double(pspec2),
	sig = as.double(sig),
	ch = as.double(ch),
	br = as.double(br),
	bi = as.double(bi),
	rat = as.double(rat))

    if( window == "Akaike" ) pspec <- z$pspec2
    if( window == "Hanning" ) pspec <- z$pspec1

    if( plot == TRUE ) {
      x <- rep(0,lag1)
      for( i in 1:lag1 ) x[i] <- (i-1)/(2*lag)
      ylab = paste("Spectrum smoothing by ",window," window")
      if( log == TRUE ) plot(x, pspec, type="l", log="y", ylab=ylab, xlab="Frequency")
      if( log == FALSE ) plot(x, pspec, type="l", ylab=ylab, xlab="Frequency") }

    bispec.out <- list( pspec=pspec, sig=z$sig, cohe=array(z$ch, dim=c(lag1,lag1)),
			breal=array(z$br, dim=c(lag1,lag1)), bimag=array(z$bi, dim=c(lag1,lag1)), exval=z$rat )
    return( bispec.out )
}


canarm <-
function (y, max.order=NULL, plot=TRUE)
{
    n <- length(y)
    if( is.null(max.order) ) max.order <- as.integer(2*sqrt(n))
    morder <- max.order

#    mmax <- 50
#    nmax <- 101
    mmax <- morder
    nmax <- 2*morder+1

    arcoef <- rep(0,nmax)	# AR-coefficients
    l1 <- 0			# upper limit of the model order +1
    v <- rep(0,mmax+1)		# innovation vector
    aic <- rep(0,mmax+1)	# AIC
    oaic <- 0              	# minimum AIC
    mo <- 0                	# order of AR
    parcor <- rep(0,mmax)	# partial auto-correlation

    nc <- 0                # total number of case
    m1 <- rep(0,mmax)      # number of present and future variables
    m2 <- rep(0,mmax)      # number of present and past variables
    w <- array(0, dim=c(mmax,mmax,mmax))  # future set canonical weight
    z <- array(0, dim=c(mmax,mmax))       # canonical R
    Rs <- array(0, dim=c(mmax,mmax))      # R-squared
    chi <- array(0, dim=c(mmax,mmax))     # chi-square
    ndt <- array(0, dim=c(mmax,mmax))     # N.D.F
    dic <- array(0,dim=c(mmax,mmax))      # DIC
    dicm <- rep(0,mmax)                   # minimum DIC
    po <- rep(0,mmax)                     # order of minimum DIC

    k <- 0                 # order of AR
    b <- rep(0,mmax)       # AR-coefficients
    l <- 0                 # order of MA
    a <- rep(0,mmax)       # MA-coefficients

    z2 <- autcor(y, morder, plot=FALSE)
    autcv <- z2$acov
    tmp.file <- " "


    z1 <- .C("canarm",
	as.integer(n),
	as.integer(morder+1),
	as.double(autcv),
	arcoef = as.double(arcoef),
	l1 = as.integer(l1),
	v = as.double(v),
	aic = as.double(aic),
	oaic = as.double(oaic),
	mo = as.integer(mo),
	parcor = as.double(parcor),
	nc = as.integer(nc),
	m1 = as.integer(m1),
	m2 = as.integer(m2),
	w = as.double(w),
	z = as.double(z),
	Rs = as.double(Rs),
	chi = as.double(chi),
	ndt = as.integer(ndt),
	dicp = as.double(dic),
	dicm = as.double(dicm),
	po = as.integer(po),
	k = as.integer(k),
	b = as.double(b),
	l = as.integer(l),
	a = as.double(a),
        as.character(tmp.file),
	as.integer(mmax),
	as.integer(nmax) )

    l1 <- z1$l1
    mo <- z1$mo
    parcor <- z1$parcor[1:(l1-1)]
    nc <- z1$nc
    m1 <- z1$m1[1:nc]
    m2 <- z1$m2[1:nc]
    w <- array(z1$w, dim=c(mmax,mmax,nc))
    z <- array(z1$z, dim=c(mmax,nc))
    Rs <- array(z1$Rs, dim=c(mmax,nc))
    chi <- array(z1$chi, dim=c(mmax,nc))
    ndt <- array(z1$ndt, dim=c(mmax,nc))
    dicp <- array(z1$dicp, dim=c(mmax,nc))
    cweight <- list()
    cR <- list()
    Rsquar <- list()
    chisquar <- list()
    ndf <- list()
    dic <- list()
    for( i in 1:nc )  cweight[[i]] <- w[(1:m1[i]),(1:m1[i]),i]
    for( i in 1:nc )  cR[[i]] <- z[(1:m1[i]),i]
    for( i in 1:nc )  Rsquar[[i]] <- Rs[(1:m1[i]),i]
    for( i in 1:nc )  chisquar[[i]] <- chi[(1:m1[i]),i]
    for( i in 1:nc )  ndf[[i]] <- ndt[(1:m1[i]),i]
    for( i in 1:nc )  dic[[i]] <- dicp[(1:m1[i]),i]
    k <- z1$k
    l <- z1$l

    if( plot == TRUE ) {
      plot(parcor, type="h", ylab="Partial Autocorrelation", xlab="Lag")
      abline(h=0, lty=1)
#      abline(h=1/sqrt(n), lty=3)
#      abline(h=-1/sqrt(n), lty=3)
      abline(h=2/sqrt(n), lty=3)
      abline(h=-2/sqrt(n), lty=3) }

    canarm.out <- list( arinit=-z1$arcoef[1:mo], v=z1$v[1:l1], aic=z1$aic[1:l1], aicmin=z1$oaic,
                        order.maice=mo, parcor=parcor, nc=nc, future=m1, past=m2, cweight=cweight,
                        canocoef=cR, canocoef2=Rsquar, chisquar=chisquar, ndf=ndf, dic=dic, dicmin=z1$dicm[1:nc],
                        order.dicmin=z1$po[1:nc], arcoef=-z1$b[1:k], macoef=-z1$a[1:l] )
    return( canarm.out )
}


canoca <-
function (y)
{
    n <- nrow(y)        # length of data
    d <- ncol(y)        # dimension of Y(I)
    lag <- as.integer(2*sqrt(n))     # maximum lag
    lag1 <- lag+1
    inw <- c(1:d)       # inw(k)=j means that the k-th component of y(i) is the j-th component of the original record z(i)

    z1 <- mulcor(y, lag, plot=FALSE)
    cov <- z1$cov       # covariance matrix

    lmax <- 12
    mj0 <- lmax+1
    mj1 <- mj0*d

    l <- 0		# upper bound of AR-order
    aic <- rep(0,mj0)	# AIC
    oaic <- 0		# minimum AIC
    mo <- 0		# MAICE AR-model order                 
    v <- array(0, dim=c(d,d))		# innovation variance
    ac <- array(0, dim=c(mj0,d,d))	# autoregressive coefficients
    nc <- 0				# number of cases
    m1 <- rep(0,mj1)			# number of variable in the future set
    m2 <- rep(0,mj1)			# number of variables in the past set
    w <- array(0,dim=c(mj1,mj1,mj1))	# future set canonical weight
    z <- array(0,dim=c(mj1,mj1))	# canonical R
    Rs <- array(0,dim=c(mj1,mj1))	# R-squared
    chi <- array(0,dim=c(mj1,mj1))	# chi-square
    ndt <- array(0,dim=c(mj1,mj1))	# N.D.F
    dic <- array(0,dim=c(mj1,mj1))	# DIC(=CHI**2-2*D.F)
    dicm <- rep(0,mj1)			# minimum DIC
    po <- rep(0,mj1)			# order of minimum DIC
    f <- array(0,dim=c(mj1,mj1))	# transition matrix F
    k <- 0				# number of structual characteristic vector
    nh <- rep(0,mj1)			# structual characteristic vector
    g <- array(0, dim=c(mj1,d))		# input matrix
    ivf <- 0				# number of vector vf
    vf <- rep(0,mj1*mj1)		# F matrix in vector form

    z1 <- .C("canoca",
	as.integer(d),
	as.integer(inw),
	as.integer(n),
	as.integer(lag1),
	as.integer(d),
	as.double(cov),
	l = as.integer(l),
	aic = as.double(aic),
        oaic = as.double(oaic),
	mo = as.integer(mo),
	v = as.double(v),
	ac = as.double(ac),
	nc = as.integer(nc),
	m1 = as.integer(m1),
	m2 = as.integer(m2),
	w = as.double(w),
	z = as.double(z),
	Rs = as.double(Rs),
	chi = as.double(chi),
	ndt = as.integer(ndt),
	dic = as.double(dic),
	dicm = as.double(dicm),
	po = as.integer(po),
	f = as.double(f),
	k = as.integer(k),
	nh = as.integer(nh),
	g = as.double(g),
	ivf = as.integer(ivf),
	vf = as.double(vf),
	as.integer(lmax),
	as.integer(mj0),
	as.integer(mj1) )

    l <- z1$l
    mo <- z1$mo
    ac <- array(z1$ac, dim=c(mj0,d,d))
    arcoef <- array(,dim=c(d,d,mo))
    for( i in 1:mo ) arcoef[,,i] <- -ac[i,(1:d),(1:d)]
    nc <- z1$nc
    m1 <- z1$m1[1:nc]
    m2 <- z1$m2[1:nc]
    w <- array(z1$w, dim=c(mj1,mj1,mj1))
    z <- array(z1$z, dim=c(mj1,mj1))
    Rs <- array(z1$Rs, dim=c(mj1,mj1))
    chi <- array(z1$chi, dim=c(mj1,mj1))
    ndt <- array(z1$ndt, dim=c(mj1,mj1))
    dicp <- array(z1$dic, dim=c(mj1,mj1))
    f <- array(z1$f, dim=c(mj1,mj1))
    cweight <- list()
    cR <- list()
    Rsquar <- list()
    chisquar <- list()
    ndf <- list()
    dic <- list()
    matF <- list()
    for( i in 1:nc )  cweight[[i]] <- w[(1:m1[i]),(1:m1[i]),i]
    for( i in 1:nc )  cR[[i]] <- z[(1:m1[i]),i]
    for( i in 1:nc )  Rsquar[[i]] <- Rs[(1:m1[i]),i]
    for( i in 1:nc )  chisquar[[i]] <- chi[(1:m1[i]),i]
    for( i in 1:nc )  ndf[[i]] <- ndt[(1:m1[i]),i]
    for( i in 1:nc )  dic[[i]] <- dicp[(1:m1[i]),i]
    for( i in 1:nc )  matF[[i]] <- f[1:(m1[i]-1),i]
    k <- z1$k
    g <- array(z1$g, dim=c(mj1,d))
    ivf <- z1$ivf

    canoca.out <- list( aic=z1$aic[1:(l+1)], aicmin=z1$oaic, order.maice=mo, v=array(z1$v, dim=c(d,d)),
			arcoef=arcoef, nc=nc, future=m1, past=m2, cweight=cweight, canocoef=cR, canocoef2=Rsquar,
			chisquar=chisquar, ndf=ndf, dic=dic, dicmin=z1$dicm[1:nc], order.dicmin=z1$po[1:nc],
                        matF=matF, vectH=z1$nh[1:k], matG=g[(1:k),(1:d)], vectF=z1$vf[1:ivf] )
    return( canoca.out )
}


covgen <-
function (lag, f, gain, plot=TRUE)
{
    k <- length(f)	# number of data points

    acov <- rep(0,(lag+1))	# auto covariance
    acor <- rep(0,(lag+1))	# auto covariance normalized
 
    z <- .C("covgen",
	as.integer(lag),
	as.integer(k),
	as.double(f),
	as.double(gain),
	acov = as.double(acov),
	acor = as.double(acor))

    if( plot == TRUE ) {
      plot((0:lag), z$acor, type="h", ylab="Auto Covariance Normalized", xlab="Lag")
      abline(h=0, lty=1) }

    covgen.out <- list( acov=z$acov, acor=z$acor )
    return( covgen.out )
}


markov <-
function (y, tmp.file=NULL)
{
    n <- nrow(y)        # length of data
    d <- ncol(y)        # dimension of the observation vector
    lag <- as.integer(2*sqrt(n))  # maximum lag
    lag1 <- lag+1
    icont <- 2		# output control (0:ARMA coefficients, 1:for SIMCON input, 2: for both)

    z1 <- mulcor(y, lag, plot=FALSE)
    cov <- z1$cov		# covariance matrix

    z2 <- canoca(y)
    nh <- z2$vectH		# structural characteristic vector
    k <- length(nh)		# dimension of the state vector
    vectF <- z2$vectF		# initial estimate of the vector of free parameters in F
    nvf <- length(vectF)	# length of vector vf
    matGi <- rep(z2$matG,1)     # initial estimates of the free parameters in G

    mj3 <- max(lag1,100)
    mj4 <- nvf+k*d
    mj6 <- 2*(k+d)-1
    mj7 <- as.integer((k-1)/d+1)

    id <- rep(0,k)      # id(i)=1 means that the i-th row of a contains free parameters
    ir <- rep(0,k)      # denotes the position of the last non-zero element within the i-th row in F
    ij <- rep(0,d)      # denotes the position of the i-th non-trivial row within F
    ik <- rep(0,d)      # denotes the number of free parameters within the i-th non-trivial row in F
    ngr <- 0				# number of gradient vector
    gr <- rep(0, mj4)			# gradient vector
    a1 <- array(0, dim=c(k,k))		# initial estimate of the transition matrix (F)
    a <- array(0, dim=c(k,k))		# transition matrix (F)
    b <- array(0, dim=c(k,d))		# input matrix (G)
    vd <- array(0, dim=c(mj4,mj4))	# DAVIDON variance
    iqm <- 0				# AR-order
    bm <- array(0, dim=c(d,d,mj7))	# AR-coefficient matrices
    au <- array(0, dim=c(d,d,mj7))	# impulse response matrices
    zz <- array(0, dim=c(d,d,mj7))	# MA-coefficient matrices
    v <- array(0, dim=c(d,d))		# inovation variance
    aic <- 0				# AIC
    if( is.null(tmp.file) )  tmp.file <- " "

    z <- .C("markov",
	as.integer(n),
	as.integer(lag1),
	as.integer(d),
	as.double(cov),
	as.integer(k),
	as.integer(nh),
	as.integer(nvf),
	as.double(vectF),
	as.double(matGi),
	as.integer(icont),
	idd = as.integer(id),
	ir = as.integer(ir),
	ij = as.integer(ij),
	ik = as.integer(ik),
	ngr = as.integer(ngr),
	gr = as.double(gr),
	a1 = as.double(a1),
	a = as.double(a),
	b = as.double(b),
	vd = as.double(vd),
	iqm = as.integer(iqm),
	bm = as.double(bm),
	au = as.double(au),
	zz = as.double(zz),
	v = as.double(v),
	aic = as.double(aic),
        as.character(tmp.file),
	as.integer(mj3),
	as.integer(mj4),
	as.integer(mj6),
	as.integer(mj7) )

    ngr <- z$ngr
    vd <- array(z$vd, dim=c(mj4,mj4))
    iqm <- z$iqm
    bm <- array(z$bm, dim=c(d,d,mj7))
    au <- array(z$au, dim=c(d,d,mj7))
    zz <- array(z$zz, dim=c(d,d,mj7))

    arcoef <- array(-bm, dim=c(d,d,iqm))
    macoef <- array(-zz, dim=c(d,d,iqm-1))

    markov.out <- list( id=z$id, ir=z$ir, ij=z$ij, ik=z$ik, grad=z$gr[1:ngr], matFi=array(z$a1, dim=c(k,k)),
			matF=array(z$a, dim=c(k,k)), matG=array(z$b, dim=c(k,d)), davvar=vd[(1:ngr),(1:ngr)],
			arcoef=arcoef, impuls=array(au, dim=c(d,d,iqm-1)),
			macoef=macoef, v=array(z$v, dim=c(d,d)), aic=z$aic )
    return( markov.out )
}


nonst <-
function (y, span, max.order=NULL, plot=TRUE)
{
    n <- length(y)	# length of data
    if( is.null(max.order) ) max.order <- as.integer(2*sqrt(n))  # highest order of AR model
    tmp.file <- " "

    morder <- max.order
    ns <- as.integer(n/span)
    p <- rep(0,ns)			# AR order
    coef <- array(0, dim=c(morder,ns))	# AR-coefficients
    v <- rep(0,ns) 			# innovation variance
    aic <- rep(0,ns)			# AIC
    daic21 <- rep(0,ns)			# AIC2-AIC1
    daic <- rep(0,ns)			# DAIC/N
    ks <- rep(0,ns)			# start point of the current model
    ke <- rep(0,ns)			# end point of the current model
    pspec <- array(0, dim=c(121,ns))

    z <- .C("nonst",
	as.integer(n),
	as.integer(span),
	as.double(y),
	as.integer(ns),
	as.integer(morder),
	p = as.integer(p),
	coef = as.double(coef),
	v = as.double(v),
	aic = as.double(aic),
	daic21 = as.double(daic21),
	daic = as.double(daic),
	ks = as.integer(ks),
	ke = as.integer(ke),
	pspec = as.double(pspec),
        as.character(tmp.file) )

    coef <- array(z$coef, dim=c(morder,ns))
    arcoef <- list() 
    for (i in 1:ns )  arcoef[[i]] <- coef[1:z$p[i],i]
    pspec <- array(z$pspec, dim=c(121,ns))

    if( plot == TRUE ) {
      x <- rep(0,121)
      for( i in 1:121 ) x[i] <- (i-1)/240
      par(mfrow=c(ns,1))
      for( i in 1:ns ) {
        plot(x, pspec[,i], type="l", main=paste("y(", z$ks[i], "),...,y(", z$ke[i], ")"),
             xlab="Frequency", ylab="Power Spectrum") }
      par(mfrow=c(1,1))
    }

    nonst.out <- list( ns=ns, arcoef=arcoef, v=z$v, aic=z$aic, daic21=z$daic21, daic=z$daic,
		       init=z$ks, end=z$ke, pspec=array(z$pspec, dim=c(121,ns)) )
    return( nonst.out )
}


prdctr <-
function (y, r, s, h, arcoef, macoef=NULL, impuls=NULL, v, plot=TRUE)
{
    if (is.array(y)) {
      n <- nrow(y)      # length of data
      d <- ncol(y)	# dimension of vector y(i)
    } else {
      n <- length(y)
      d <- 1
      y <- array(y, dim=c(n,1))
    }

    if(is.array(arcoef)) {	# AR-coefficient matrices
        p <- length(arcoef)[3]
        arcoef <- -arcoef
    } else {
        p <- length(arcoef)
        arcoef <- array(-arcoef, dim=c(1,1,p))
    }

    jsw <- 0
    if( is.null(macoef) ) {
      jsw <- 1
      if( is.array(impuls) ) {	# impulse response matrices
        q <- dim(impuls)[3]
      } else {
        q <- length(impuls)
        impuls <- array(impuls,dim=c(1,1,q))
      }
      macoef <- array(0, dim=c(d,d,q))
    } else {
      if( is.array(macoef) ) {	# MA-coefficient matrices
        q <- dim(macoef)[3]
        macoef <- -macoef
      } else {
        q <- length(macoef)
        macoef <- array(-macoef,dim=c(1,1,q))
      }
      if( is.null(impuls) ) impuls <- array(0, dim=c(d,d,q))
    }

    yreal <- array(0, dim=c(s+h,d))	# real data
    yori <- array(0, dim=c(h+1,d))	#
    ypre <- array(0, dim=c(s+h,d))	# predicted values
    ys <- array(0, dim=c(n,d))		# predicted - (real data)
    z1 <- array(0, dim=c(s+h,d))	# predicted + (standard deviation)
    z1 <- array(0, dim=c(s+h,d))	# predicted + (standard deviation)
    z2 <- array(0, dim=c(s+h,d))	# predicted + 2*(standard deviation)
    z3 <- array(0, dim=c(s+h,d))	# predicted + 3*(standard deviation)
    zz1 <- array(0, dim=c(s+h,d))	# predicted - (standard deviation)
    zz2 <- array(0, dim=c(s+h,d))	# predicted - 2*(standard deviation)
    zz3 <- array(0, dim=c(s+h,d))	# predicted - 3*(standard deviation)
    tmp.file <- " "

    z <- .C("prdctr",
	as.integer(n),
	as.integer(r),
	as.integer(s),
	as.integer(h),
	as.integer(d),
	as.integer(p),
	as.integer(q),
	as.integer(jsw),
	as.double(y),
	as.double(arcoef),
	as.double(macoef),
	as.double(impuls),
	as.double(v),
	yreal = as.double(yreal),
	yori = as.double(yori),
	ypre = as.double(ypre),
	ys = as.double(ys),
	z1 = as.double(z1),
	z2 = as.double(z2),
	z3 = as.double(z3),
	zz1 = as.double(zz1),
	zz2 = as.double(zz2),
	zz3 = as.double(zz3),
        as.character(tmp.file) )

#    yreal <- array(z$yreal, dim=c(s+h,d))
#    yori <- array(z$yori, dim=c(h+1,d))
#    for(i in 1:(n-s+1)) yreal[s+i-1,] <- yori[i,]
#    for(j in 1:d) yreal[,j] <- yreal[(1:n),j]

    predct <- array(z$ypre, dim=c(s+h,d))
    for(i in 1:(r-1)) predct[i,] <- NA

#    ys <- array(z$ys, dim=c(n,d))
#    for(i in 1:(r-1)) ys[i,] <- NA
#    for(i in s:n) ys[i,] <- NA
    ys <- array(NA, dim=c(n,d))
    for(j in 1:d)
      for(i in r:n) ys[i,j] <- y[i,j]-predct[i,j]

    pstd <- array(z$z1, dim=c(s+h,d))
    pstd[1:s-1,] <- NA
    p2std <- array(z$z2, dim=c(s+h,d))
    p2std[1:s-1,] <- NA
    p3std <- array(z$z3, dim=c(s+h,d))
    p3std[1:s-1,] <- NA
    mstd <- array(z$zz1, dim=c(s+h,d))
    mstd[1:s-1,] <- NA
    m2std <- array(z$zz2, dim=c(s+h,d))
    m2std[1:s-1,] <- NA
    m3std <- array(z$zz3, dim=c(s+h,d))
    m3std[1:s-1,] <- NA

    if( plot == TRUE ) {
      par(mfrow=c(1,d))
      for( i in 1:d ) {
        ymin <- min(y[(1:n),i], predct[r:(s+h),i])
        ymax <- max(y[(1:n),i], predct[r:(s+h),i])
        plot(y[,i], type="l", xlim=c(0,(s+h)), ylim=c(ymin,ymax), xlab="Time", ylab="real data / predicted values")
        par(new=TRUE)
        plot(predct[,i], type="l", xlim=c(0,(s+h)), ylim=c(ymin,ymax), col="red", xlab="", ylab="") }
      par(mfrow=c(1,1))
    }

    prdctr.out <- list( predct=predct, ys=ys, pstd=pstd, p2std=p2std, p3std=p3std,
                        mstd=mstd, m2std=m2std, m3std=m3std )

    class( prdctr.out ) <- "prdctr"
    return( prdctr.out )
}

print.prdctr <- function(x, ...)
{
  n <- dim(x$ys)[1]
  d <- dim(x$ys)[2]
  sh <- dim(x$predct)[1]
  for (i in n:1 ) {
    if ( is.na(x$ys[i,1]) == FALSE ) r <- i
    if ( is.na(x$pstd[i,1]) == FALSE ) s <- i
  }

  for ( id in 1:d ) {
    cat(sprintf("\n\n d = %i\n\n", id))
    cat("\tPREDICTED-REAL\tPREDICTED\tPREDICTED\tPREDICTED\tPREDICTED\n")
    cat("\t\t\tVALUES\t\t+/-STANDARD\t+/-2STANDARD\t+/-3STANDARD\n")
    cat("\t\t\t\t\tDEVIATION\tDEVIATION\tDEVIATION\n")
    for( i in r:(s-1) ) cat(sprintf("N= %i\t%f\t%f\n", i,x$ys[i,id],x$predct[i,id]))
    for( i in s:n ) {
      cat(sprintf("N= %i\t%f\t%f\t%f\t%f\t%f\n", i,x$ys[i,id],x$predct[i,id],x$pstd[i,id],x$p2std[i,id],x$p3std[i,id]))
      cat(sprintf("\t\t\t\t\t%f\t%f\t%f\n", x$mstd[i,id],x$m2std[i,id],x$m3std[i,id])) }
    for( i in (n+1):sh ) {
      cat(sprintf("N= %i\t\t\t%f\t%f\t%f\t%f\n", i,x$predct[i,id],x$pstd[i,id],x$p2std[i,id],x$p3std[i,id]))
      cat(sprintf("\t\t\t\t\t%f\t%f\t%f\n", x$mstd[i,id],x$m2std[i,id],x$m3std[i,id])) }
  }
}


simcon <-
function (span, len, r, arcoef, impuls, v, weight)
{
    arcoef <- -arcoef	# matrices of autoregressive coefficients
    d <- dim(arcoef)[1]		# dimension of Y(I)
    k <- dim(arcoef)[3]		# order of the process

    bc <- array(0, dim=c(k*d,r))	# kd*(d-r) sub matrices of BX
    bd <- array(0, dim=c(k*d,d-r))	# kd*r sub matrices of BX
    g <- array(0, dim=c(r,k*d))		# controller gain G
    av <- rep(0, d)	# average value of i-th componentof Y
    si <- rep(0, d)	# variance
    s2 <- rep(0, d)	# standard deviation
    
    z <- .C("simcon",
	as.integer(d),
	as.integer(k),
	as.integer(span),
	as.integer(len),
	as.integer(r),
	as.double(arcoef),
	as.double(impuls),
	as.double(v),
	as.double(weight),
	bc = as.double(bc),
	bd = as.double(bd),
	g = as.double(g),
	av = as.double(av),
	si = as.double(si),
	s2 = as.double(s2))

    simcon.out <- list( gain=array(z$g, dim=c(r,k*d)), ave=z$av[1:d], var=z$si[1:d], std=z$s2[1:d],
			bc=array(z$bc, dim=c(k*d,r)), bd=array(z$bd, dim=c(k*d,d-r)) )
    return( simcon.out )
}


thirmo <-
function (y, lag=NULL, plot=TRUE)
{
    n <- length(y)
    if( is.null(lag) ) lag <- as.integer(2*sqrt(n))  # maximum lag
    lag1 <- lag+1

    mean <- 0		# mean
    acov <- rep(0,lag1)	# autocovariance
    acor <- rep(0,lag1)		# normalized covariance
    mnt <- array(0, dim=c(lag1,lag1))	# third order moments

    z <- .C("thirmo",
	as.integer(n),
	as.integer(lag),
	as.double(y),
	mean = as.double(mean),
	acov = as.double(acov),
	acor = as.double(acor),
	mnt = as.double(mnt))

    mnt <- array(z$mnt, dim=c(lag1,lag1))
    tmomnt <- list()
    for( i in 1:lag1 ) tmomnt[[i]] <- mnt[i,1:i]

    if( plot == TRUE ) {
      plot((0:lag), z$acor, type="h", ylab="Auto Covariance Normalized", xlab="Lag")
      abline(h=0, lty=1)
    }

    thirmo.out <- list( mean=z$mean, acov=z$acov, acor=z$acor, tmomnt=tmomnt )
    return( thirmo.out )
}


#####   TIMSAC78   #####

blocar <-
function (y, max.order=NULL, span, plot=TRUE)
{
    n <- length(y)
    if( is.null(max.order) ) max.order <- as.integer(2*sqrt(n))  # upper limit of the order of AR model

    morder <- max.order
    ns <- as.integer((n-morder+span-1)/span)
    mean <- 0
    var <- 0
    aic <- array(0, dim=c(ns,ns))
    bw <- array(0, dim=c(ns,ns))  # Bayesian Weight
    b <- array(0, dim=c(morder,ns))  # Partial Autocorrelation
    a <- array(0, dim=c(morder,ns))  # Coefficients ( average by the Bayesian weights )
    v  <- rep(0,ns)                 # Innovation Variance
    ks <- rep(0,ns)                 # initial point of data
    ke <- rep(0,ns)                 # end point of data
    pxx <- array(0, dim=c(121,ns))  # power spectrum

    z <- .C("blocar",
	as.double(y),
	as.integer(n),
	as.integer(morder),
	as.integer(span),
	as.integer(ns),
	mean = as.double(mean),
	var = as.double(var),
	aic = as.double(aic),
	bw = as.double(bw),
	b = as.double(b),
	a = as.double(a),
	v = as.double(v), 
	ks = as.integer(ks),
	ke = as.integer(ke),
	pxx = as.double(pxx) )

    aic <- list()
    bweight <- list()
    pacoef <- list()
    arcoef <- list()
    aic[[1]] <- NA
    bweight[[1]] <- NA
    for( i in 2:ns ) {
        j <- ((i-1)*ns+1):((i-1)*ns+i)
        aic[[i]] <- z$aic[j]
        bweight[[i]] <- z$bw[j]
    }
    for( i in 1:ns ) pacoef[[i]] <- z$b[((i-1)*morder+1):(i*morder)]
    for( i in 1:ns ) arcoef[[i]] <- z$a[((i-1)*morder+1):(i*morder)]
    pspec <- array(z$pxx, dim=c(121,ns))

    if( plot == TRUE ) {
      x <- rep(0,121)
      for( i in 1:121 ) x[i] <- (i-1)/240
      par(mfrow=c(ns,1))
      for( i in 1:ns ) plot(x, pspec[,i], type="l", main=paste("y(", z$ks[i], "),...,y(", z$ke[i], ")"),
                            xlab="Frequency", ylab="Power Spectrum")
      par(mfrow=c(1,1))
    }

    blocar.out <- list( mean=z$mean, var=z$var, aic=aic, bweight=bweight, pacoef=pacoef, arcoef=arcoef,
			v=z$v, init=z$ks, end=z$ke, pspec=pspec )
    return( blocar.out )
}


blomar <-
function (y, max.order=NULL, span)
{
    n <- nrow(y)
    d <- ncol(y)
    if( is.null(max.order) ) max.order <- as.integer(2*sqrt(n))  # upper limit of the order of AR-model
    morder <- max.order

    calb<-rep(1,d)   # calibration constant for channel j (j=1,d)
    ns <- as.integer((n-morder+span-1)/span)
    mean <- rep(0,d)
    var <- rep(0,d)
    bw <- array(0, dim=c(ns,ns))        # Bayesian Weight
    raic <- array(0, dim=c(ns,ns))      # AIC
    a <- array(0, dim=c(d,d,morder,ns))    # AR-coefficient matrices
    e <- array(0, dim=c(d,d,ns))        # innovation variance
    aic <- rep(0,ns)                    # equivalent AIC of Bayesian model
    ks <- rep(0,ns)                     # initial point of data
    ke <- rep(0,ns)                     # end point of data

    z <- .C("blomar",
	as.double(y),
	as.integer(n),
	as.integer(d),
	as.double(calb),
	as.integer(morder),
	as.integer(span),
	as.integer(ns),
	mean = as.double(mean),
	var = as.double(var),
	bw = as.double(bw),
	raic = as.double(raic),
	a = as.double(a),
	e = as.double(e),
	aic = as.double(aic),
	ks = as.integer(ks),
	ke = as.integer(ke) )

    bw <- array(z$bw, dim=c(ns,ns))
    bweight <- list()
    bweight[[1]] <- NA
    for( i in 2:ns ) bweight[[i]] <- bw[1:i,i]

    raic <- array(z$raic, dim=c(ns,ns))
    aic <- list()
    aic[[1]] <- NA
    for( i in 2:ns ) aic[[i]] <- raic[1:i,i]

    a <- array(z$a, dim=c(d,d,morder,ns))
    arcoef <- list()
    for( i in 1:ns ) arcoef[[i]] <- array(a[,,,i], dim=c(d,d,morder))

    e <- array(z$e, dim=c(d,d,ns))
    v <- list()
    for( i in 1:ns ) v[[i]] <- array(e[,,i], dim=c(d,d))

    blomar.out <- list( mean=z$mean, var=z$var, bweight=bweight, aic=aic, arcoef=arcoef, v=v, eaic=z$aic,
			init=z$ks, end=z$ke )
    class( blomar.out ) <- "blomar"
    return( blomar.out )
}

print.blomar <- function(x, ...)
{

  cat("\n\n I\tMEAN\t\tVARIANCE\n")
  id <- length(x$mean)
  for( i in 1:id ) cat(sprintf(" %i\t%f\t%f\n", i,x$mean[i],x$var[i]))

  ns <- length(x$bweight)
  for( i in 1:ns ) {
    if( i != 1 ) {
      cat("\n\nAR-MODEL FITTED TO\t! BAYESIAN WEIGHTS\t! AIC WITH RESPECT TO THE PRESENT DATA\n")
      cat("--------------------------------------------------------------------------------------\n")
      cat(sprintf("CURRENT BLOCK\t\t! %f\t\t! %f\n", x$bweight[[i]][1], x$aic[[i]][1]))
      for( k in 1:(i-1) ) cat(sprintf("%i PERIOD FORMER BLOCK\t! %f\t\t! %f\n", k,x$bweight[[i]][k+1],x$aic[[i]][k+1]))
    }
    cat("\n..........  CURRENT MODEL  ..........\n\n")
    cat(sprintf(" M\tAM(I,J)\t\t\tDATA  Z(K,.); K= %i,%i\n", x$init[i],x$end[i]))
    id <- dim(x$arcoef[[i]])[1]
    mf <- dim(x$arcoef[[i]])[3]
    for( j in 1:mf )
      for( k1 in 1:id ) {
        if( k1 == 1 ) cat(sprintf(" %i",j))
        for( k2 in 1:id ) cat(sprintf("\t%f", x$arcoef[[i]][k1,k2,j]))
        cat("\n")
      }
    cat(sprintf("\n\nORDER = %i\nAIC = %f\n", mf,x$eaic[i]))
    cat("\nINNOVATION VARIANCE MATRIX\n")
    print(x$v[[i]])
  }
}


bsubst <-
function (y, mtype, lag=NULL, nreg, reg=NULL, term.lag=NULL, cstep=5, plot=TRUE)
{
    n <- length(y)
    if( is.null(lag) ) lag <- as.integer(2*sqrt(n))  # maximum time lag use in the model
    k <- nreg                                        # number of regressors
    if(is.null(reg)) reg <- array(0, dim=c(3,nreg))  # specification of regressor (mtype=2,6)
    if(is.null(term.lag)) term.lag <- rep(0,5)       # specification of regressor (mtype=3)
    f <- ""                                          # specification of regressor (mtype=5)
    cnst <- 0                                        # constant value (mtype=6)

    ymean <- 0
    yvar <- 0
    m <- 0
    aicm <- 0
    vm <- 0
    a1 <- rep(0,k)       # AR-coefficients
    v  <- rep(0,k+1)
    aic <- rep(0,k+1)
    daic <- rep(0,k+1)
    aicb <- 0            # AIC of Bayesian model
    vb <- 0              # residual variance of Bayesian model
    ek <- 0              # equivalent number of parameters
    a2 <- rep(0,k)       # AR-coefficients of Bayesian model
    ind <- rep(0,k)      # index of c(i) in order of increasing magnitude
    c <- rep(0,k)        # square of partial correlations (c(i)= n*b(i)**2)
    c1 <- rep(0,k+1)     # binomial type damper
    c2 <- rep(0,k)       # final Bayesian weights of partial correlations
    b <- rep(0,k)        # partial correlations of the Bayesian model b(i)
    eicmin <- 0          # minimum EIC
    esum <- rep(0,k+1)
    npmean <- 0          # mean of number of parameter
    npmean.nreg <- 0     # (=npmean/nreg)
    e <- array(0, dim=c(n,cstep))  # prediction error
    mean <- rep(0,cstep) # mean
    var <- rep(0,cstep)  # variance
    skew <- rep(0,cstep) # skewness
    peak <- rep(0,cstep) # peakedness
    cov <- rep(0,101)    # autocorrelation function
    pxx <- rep(0,121)    # power spectrum

    z <- .C("bsubst",
	as.double(y),
	as.integer(n),
	as.integer(mtype),
	as.integer(lag),
	as.integer(nreg),
	as.integer(cstep),
 	as.integer(reg),
	as.integer(term.lag),
	as.character(f),
	as.double(cnst),
	ymean = as.double(ymean),
	yvar = as.double(yvar),
	m = as.integer(m),
	aicm = as.double(aicm),
	vm = as.double(vm),
	a1 = as.double(a1),
	v = as.double(v),
	aic = as.double(aic),
	daic = as.double(daic),
	aicb = as.double(aicb),
	vb = as.double(vb),
	ek = as.double(ek),
	a2 = as.double(a2),
	ind = as.integer(ind),
	c = as.double(c),
	c1 = as.double(c1),
	c2 = as.double(c2),
	b = as.double(b),
	eicmin = as.double(eicmin),
	esum = as.double(esum),
	npmean = as.double(npmean),
	npmean.nreg = as.double(npmean.nreg),
	e = as.double(e),
	mean = as.double(mean),
	var = as.double(var),
	skew = as.double(skew),
	peak = as.double(peak),
	cov = as.double(cov),
	pxx = as.double(pxx) )

    perr <- array(z$e, dim=c(n,cstep))

    if( plot == TRUE ) {
      mm <- 2
      if( mtype==1 ) mm <- mm+1
      nc <- (mm+cstep+1)/2
      par(mfcol=c(nc,2))

      for( i in 1:cstep ) {
        sig <- sqrt(z$var[i])*0.5
        hist( perr[(lag+1):n,i], breaks=c(-sig*10:1,0,sig*1:10), main="", xlab=paste(i,"-step ahead prediction error") )
      }

      plot((0:nreg), z$daic, ylim=c(0,40), type="l", xlab="Lag", ylab="AIC(M)-AICMIN" )
      abline(h=0, lty=1)

      if( mtype==1 ) {
        x <- rep(0,121)
        for( i in 1:121 ) x[i] <- (i-1)/240
        plot(x, pxx, type="l", xlab="Frequency", ylab=paste("Power Spectrum"))
      }

      plot( order<-c(0:100),z$cov,ylim=c(-1,1),bty="l",type="h", main="Autocorrelation of\n1-step ahead prediction error",
            xlab="Lines show +/-2sd\n ( sd = sqrt(1/n) )", ylab="peautcor" )
      abline( h=0, lty=1 )
      abline( h=2*sqrt(1/(n-lag)), lty=3 )
      abline( h=-2*sqrt(1/(n-lag)), lty=3 )
      par(mfrow=c(1,1)) }

    if(mtype == 1) 
	bsubst.out <- list( ymean=z$ymean, yvar=z$yvar, v=z$v, aic=z$aic, aicmin=z$aicm, daic=z$daic, order.maice=z$m,
			    v.maice=z$vm, arcoef.maice=z$a1, v.bay=z$vb, aic.bay=z$aicb, np.bay=z$ek, arcoef.bay=z$a2,
			    ind.c=z$ind, parcor2=z$c, damp=z$c1, bweight=z$c2, parcor.bay=z$b, eicmin=z$eicmin,
			    esum=z$esum, npmean=z$npmean, npmean.nreg=z$npmean.nreg, perr=perr, mean=z$mean,
			    var=z$var, skew=z$skew, peak=z$peak, peautcor=z$cov, pspec=z$pxx )
    if(mtype != 1)
	bsubst.out <- list( ymean=z$ymean, yvar=z$yvar, v=z$v, aic=z$aic, aicmin=z$aicm, daic=z$daic, order.maice=z$m,
			    v.maice=z$vm, arcoef.maice=z$a1, v.bay=z$vb, aic.bay=z$aicb, np.bay=z$ek, arcoef.bay=z$a2,
			    ind.c=z$ind, parcor2=z$c, damp=z$c1, bweight=z$c2, parcor.bay=z$b, eicmin=z$eicmin,
			    esum=z$esum, npmean=z$npmean, npmean.nreg=z$npmean.nreg, perr=perr, mean=z$mean,
			    var=z$var, skew=z$skew, peak=z$peak, peautcor=z$cov )
    return( bsubst.out )
}


exsar <-
function (y, max.order=NULL, plot=FALSE, tmp.file=NULL)
{
    n <- length(y)
    if( is.null(max.order) ) max.order <- as.integer(2*sqrt(n))
    if( is.null(tmp.file) )  tmp.file <- " "

    morder <- max.order
    mean <-0
    var <-0
    v  <- rep(0,morder+1)
    aic <- rep(0,morder+1)
    daic <- rep(0,morder+1)
    m <- 0               # MAICE order
    aicm <- 0            # minimum AIC
    sdm1 <- 0            # MAICE innovation variance
    a1 <- rep(0,morder)  # MAICE AR-coefficients
    sdm2 <- 0            # maximum likelihood estimates of innovation variance
    a2 <- rep(0,morder)  # maximum likelihood estimates of AR-coefficients

    z <- .C("exsar",
	as.double(y),
	as.integer(n),
	as.integer(morder),
	mean = as.double(mean),
	var = as.double(var),
	v = as.double(v),
	aic = as.double(aic),
	daic = as.double(daic),
	m = as.integer(m),
	aicm = as.double(aicm),
	sdm1 = as.double(sdm1),
	a1 = as.double(a1), 
	sdm2 = as.double(sdm2),
	a2 = as.double(a2),
	as.character(tmp.file) )

    if( plot == TRUE ) {
      plot((0:morder), z$daic, ylim=c(0,40), type="l", xlab="Lag", ylab="AIC(M)-AICMIN (Truncated at 40.0)")
      abline(h=0, lty=1) }

    exsar.out  <- list( mean=z$mean, var=z$var, v=z$v, aic=z$aic, aicmin=z$aicm, daic=z$daic, order.maice=z$m, 
			v.maice=z$sdm1, arcoef.maice=z$a1[1:z$m], v.mle=z$sdm2, arcoef.mle=z$a2[1:z$m] )
    return( exsar.out )
}


mlocar <-
function (y, max.order=NULL, span, const=0, plot=TRUE)
{
    n <- length(y)
    if( is.null(max.order) ) max.order <- as.integer(2*sqrt(n))
    morder <- max.order

    ns <- as.integer((n-morder+span-1)/span)
    mean <- 0
    var <- 0
    a <- array(0, dim=c(morder+const,ns)) # current model : AR-coefficients
    mf <- rep(0,ns)                  #               : order
    sdf  <- rep(0,ns)                #               : innovation variance
    ks <- rep(0,ns)                  #               : initial point of data
    ke <- rep(0,ns)                  #               : end point of data
    pxx <- array(0, dim=c(121,ns))   #               : power spectrum
    ld1 <- rep(0,ns)                 # data length of the preceding stationary block
    ld2 <- rep(0,ns)                 # data length of new block
    ms <- rep(0,ns)                  # moving model   : order
    sdms <- rep(0,ns)                #                : innovation variance
    aics <- rep(0,ns)                #                : AIC 
    mp <- rep(0,ns)                  # constant model : order
    sdmp <- rep(0,ns)                #                : innovation variance
    aicp <- rep(0,ns)                #                : AIC 

    z <- .C("mlocar",
	as.double(y),
	as.integer(n),
	as.integer(morder),
	as.integer(span),
	as.integer(const),
	as.integer(ns),
	mean = as.double(mean),
	var = as.double(var),
	a = as.double(a),
	mf = as.integer(mf),
	sdf = as.double(sdf),
	ks = as.integer(ks),
	ke = as.integer(ke),
	pxx = as.double(pxx),
	ld1 = as.integer(ld1),
	ld2 = as.integer(ld2),
	ms = as.integer(ms),
	sdms = as.double(sdms),
	aics = as.double(aics),
	mp = as.integer(mp),
	sdmp = as.double(sdmp),
	aicp = as.double(aicp) )

    a <- array(z$a, dim=c(morder+const,ns))
    arcoef <- list()
    for(i in 1:ns) arcoef[[i]] <- a[1:z$mf[i],i]
    pspec <- array(z$pxx, dim=c(121,ns))
    npre <- z$ld1
    order.const=z$mp
    v.const=z$sdmp
    aic.const=z$aicp
    npre[1] <- NA
    order.const[1] <- NA
    v.const[1] <- NA
    aic.const[1] <- NA

    if( plot == TRUE ) {
      x <- rep(0,121)
      for( i in 1:121 ) x[i] <- (i-1)/240
      par(mfrow=c(ns,1))
      for(i in 1:ns) plot(x, pspec[,i], type="l", main=paste("y(", z$ks[i], "),...,y(", z$ke[i], ")"),
                          xlab="Frequency", ylab="Power Spectrum")
      par(mfrow=c(1,1))
    }

    mlocar.out <- list( mean=z$mean, var=z$var, ns=ns, order=z$mf, arcoef=arcoef, v=z$sdf, init=z$ks, end=z$ke, 
			pspec=pspec, npre=npre, nnew=z$ld2, order.mov=z$ms, v.mov=z$sdms, aic.mov=z$aics,
			order.const=order.const, v.const=v.const, aic.const=aic.const )
    return( mlocar.out )
}


mlomar <-
function (y, max.order=NULL, span, const=0)
{
    n <- nrow(y)
    d <- ncol(y)
    if( is.null(max.order) ) max.order <- as.integer(2*sqrt(n))
    morder <- max.order

    calb <- rep(1,d)   # calibration for channel j (j=1,d)
    ns <- as.integer((n-morder+span-1)/span)
    mean <- rep(0,d)
    var <- rep(0,d)
    ld1<- rep(0,ns)    # data length of the preceding stationary block
    ld2 <- rep(0,ns)   # data length of new block
    ms <- rep(0,ns)    # moving model   : AR-order
    aicm <- rep(0,ns)  #                : aic
    mp <- rep(0,ns)    # constant model : AR-order
    aicc <- rep(0,ns)  #                : aic
    mf <- rep(0,ns)                   # current model : order
    aic <- rep(0,ns)                  #               : aic
    a <- array(0, dim=c(d,d,morder,ns))  #               : AR-coefficient matrices
    e <- array(0, dim=c(d,d,ns))      #               : innovation variance
    ks <- rep(0,ns)                   #               : initial point of data
    ke <- rep(0,ns)                   #               : end point of data

    z <- .C("mlomar",
	as.double(y),
	as.integer(n),
	as.integer(d),
	as.double(calb),
	as.integer(morder),
	as.integer(span),
	as.integer(const),
	as.integer(ns),
	mean = as.double(mean),
	var = as.double(var),
	ld1 = as.integer(ld1),
	ld2 = as.integer(ld2),
	ms = as.integer(ms),
	aicm = as.double(aicm),
	mp = as.integer(mp),
	aicc = as.double(aicc),
	mf = as.integer(mf),
	aic = as.double(aic),
	a = as.double(a),
	e = as.double(e),
	ks = as.integer(ks),
	ke = as.integer(ke) )

    npre=z$ld1
    order.mov=z$ms
    aic.mov=z$aicm
    order.const=z$mp
    aic.const=z$aicc
    npre[1] <- NA
    order.mov[1] <- NA
    aic.mov[1] <- NA
    order.const[1] <- NA
    aic.const[1] <- NA

    a <- array(z$a, dim=c(d,d,morder,ns))
    arcoef <- list()
    for( i in 1:ns ) arcoef[[i]] <- array(a[,,,i],dim=c(d,d,z$mf[i]))

    e <- array(z$e, dim=c(d,d,ns))
    v <- list()
    for( i in 1:ns ) v[[i]] <- array(e[,,i],dim=c(d,d))

    mlomar.out <- list( mean=z$mean, var=z$var, ns=ns, order=z$mf, aic=z$aic, arcoef=arcoef, v=v,
			init=z$ks, end=z$ke, npre=npre, nnew=z$ld2, order.mov=order.mov, aic.mov=aic.mov,
			order.const=order.const, aic.const=aic.const )
    class( mlomar.out ) <- "mlomar"
    return( mlomar.out )
}

print.mlomar <- function(x, ...)
{

  cat("\n\n I\tMEAN\t\tVARIANCE\n")
  id <- length(x$mean)
  for( i in 1:id ) cat(sprintf(" %i\t%f\t%f\n", i,x$mean[i],x$var[i]))

  ns <- x$ns
  for( i in 1:ns ) {
    if( i == 1 )
      cat(sprintf("\n\n INITIAL LOCAL MODEL:  NS = %i\t\tAIC = %f\n", x$nnew[i],x$aic[i]))
    if( i != 1 ) {
      cat("\n\n ---  THE FOLLOWING TWO MODELS ARE COMPARED  ---\n\n")
      np <- x$npre[i]+x$nnew[i]
      cat(sprintf("  MOVING MODEL:   (NF = %i\tNS = %i)\tMS = %i\t\tAIC = %f\n", x$npre[i], x$nnew[i], x$order.mov[i], x$aic.mov[i]))
      cat(sprintf("  CONSTANT MODEL: (NP = %i)\t\t\tMP = %i\t\tAIC = %f\n", np,x$order.const[i],x$aic.const[i]))
      if( x$aic.mov[i] < x$aic.const[i] ) {
        cat("\n *****     NEW MODEL ADOPTED     *****\n")
      } else {
        cat("\n *****  CONSTANT MODEL ADOPTED  *****\n")
      }
    }
    cat("\n\n..........  CURRENT MODEL  ..........\n\n")
    cat(sprintf(" M\tAM(I,J)\t\t\tDATA  Z(K,.); K= %i,%i\n", x$init[i],x$end[i]))
    id <- dim(x$arcoef[[i]])[1]
    mf <- x$order[i]
    for( j in 1:mf )
      for( k1 in 1:id ) {
        if( k1 == 1 ) cat(sprintf(" %i",j))
        for( k2 in 1:id ) cat(sprintf("\t%f", x$arcoef[[i]][k1,k2,j]))
        cat("\n")
      }
    cat(sprintf("\n\nORDER = %i\nAIC = %f\n", mf,x$aic[i]))
    cat("\nINNOVATION VARIANCE MATRIX\n")
    print(x$v[[i]])
  }

}


mulbar <-
function (y, max.order=NULL, plot=FALSE)
{
    n <- nrow(y)
    d <- ncol(y)
    if( is.null(max.order) ) max.order <- as.integer(2*sqrt(n))  # upper limit of the order of AR-model
    morder <- max.order

    calb<-rep(1,d)   # calibration of channel i (i=1,d)
    mean <- rep(0,d)
    var <- rep(0,d)
    v <- rep(0,morder+1)
    aic <- rep(0,morder+1)
    daic <- rep(0,morder+1)
    m <- 0
    aicm <- 0
    vm <- 0
    w1 <- rep(0,morder+1)             # Bayesian weights
    w2 <- rep(0,morder)               # integrated Bayesian Weights
    a <- array(0, dim=c(d,d,morder))  # AR-coefficients (forward model)
    b <- array(0, dim=c(d,d,morder))  # AR-coefficients (backward model)
    g <- array(0, dim=c(d,d,morder))  # partial autoregression coefficients (forward model)
    h <- array(0, dim=c(d,d,morder))  # partial autoregression coefficients (backward model)
    e <- array(0, dim=c(d,d))         # innovation variance matrix
    aicb <- 0                         # equivalent AIC of the Bayesian (forward) model

    z <- .C("mulbar",
	as.double(y),
	as.integer(n),
	as.integer(d),
	as.double(calb),
	as.integer(morder),
	mean = as.double(mean),
	var = as.double(var),
	v = as.double(v),
	aic = as.double(aic),
	daic = as.double(daic),
	m = as.integer(m),
	aicm = as.double(aicm),
	vm = as.double(vm),
	w1 = as.double(w1),
	w2 = as.double(w2),
	a = as.double(a),
	b = as.double(b),
	g = as.double(g),
	h = as.double(h),
	e = as.double(e),
	aicb = as.double(aicb) )

    if( plot == TRUE ) {
      plot( (0:morder), z$daic, ylim=c(0,40), type="l", xlab="Lag", ylab="AIC(M)-AICMIN (Truncated at 40.0)" )
      abline(h=0, lty=1) }

    mulbar.out <- list( mean=z$mean, var=z$var, v=z$v, aic=z$aic, aicmin=z$aicm, daic=z$daic,
			order.maice=z$m, v.maice=z$vm, bweight=z$w1, integra.bweight=z$w2,
			arcoef.for=array(z$a, dim=c(d,d,morder)), arcoef.back=array(z$b, dim=c(d,d,morder)),
			pacoef.for=array(z$g, dim=c(d,d,morder)), pacoef.back=array(z$h, dim=c(d,d,morder)),
			v.bay=array(z$e, dim=c(d,d)), aic.bay=z$aicb )
    return( mulbar.out )
}


mulmar <-
function (y, max.order=NULL, plot=FALSE, tmp.file=NULL) 
{
    n <- nrow(y)
    d <- ncol(y)
    if( is.null(max.order) ) max.order <- as.integer(2*sqrt(n))
    lag1 <- max.order+1

    calb <- rep(1, d)
    mean <- rep(0, d)
    var <- rep(0, d)
    v <- array(0, dim = c(lag1, d))
    aic <- array(0, dim = c(lag1, d))
    daic <- array(0, dim = c(lag1, d))
    m <- rep(0, d)
    aicm <- rep(0, d)
    vm <- rep(0, d)
    npr <- rep(0, d)
    jnd <- array(0, dim = c(lag1 * d, d))
    a <- array(0, dim = c(lag1 * d, d))
    rv <- rep(0, d)
    aicf <- rep(0, d)
    ei <- array(0, dim = c(d, d))
    bi <- array(0, dim = c(d, d, lag1))
    matv <- array(0, dim = c(d, d))
    arcoef <- array(0, dim = c(d, d, lag1))
    morder <- 0
    aics <- 0
    if (is.null(tmp.file))  tmp.file <- " "


    z <- .C("mulmar", as.double(y), as.integer(n), as.integer(d), 
        as.double(calb), as.integer(max.order), mean = as.double(mean), 
        var = as.double(var), v = as.double(v), aic = as.double(aic), 
        daic = as.double(daic), m = as.integer(m), aicm = as.double(aicm), 
        vm = as.double(vm), npr = as.integer(npr), jnd = as.integer(jnd), 
        a = as.double(a), rv = as.double(rv), aicf = as.double(aicf), 
        ei = as.double(ei), bi = as.double(bi), matv = as.double(matv), 
        arcoef = as.double(arcoef), morder = as.integer(morder), 
        aics = as.double(aics), as.character(tmp.file))
    v <- list()
    aic <- list()
    daic <- list()
    for (i in 1:d) {
        j <- ((i - 1) * lag1 + 1):(i * lag1)
        v[[i]] <- z$v[j]
        aic[[i]] <- z$aic[j]
        daic[[i]] <- z$daic[j]
    }
    jnd <- list()
    subregcoef <- list()
    ind <- array(z$jnd, dim = c(lag1 * d, d))
    a <- array(z$a, dim = c(lag1 * d, d))
    for (i in 1:d) jnd[[i]] <- ind[(1:z$npr[[i]]), i]
    for (i in 1:d) subregcoef[[i]] <- a[(1:z$npr[[i]]), i]

    if (plot == TRUE) {
        par(mfrow=c(d,1))  
        for (i in 1:d) {
            plot((0:max.order), daic[[i]], ylim = c(0, 40), type = "l", 
                main = paste(" d=", i), xlab = "Lag", ylab = "AIC(M)-AICMIN (Truncated at 40.0)")
            abline(h = 0, lty = 1)
        }
        par(mfrow=c(1,1))
    }
    mulmar.out <- list(mean = z$mean, var = z$var, v = v, aic = aic, aicmin = z$aicm,
        daic = daic, order.maice = z$m, v.maice = z$vm, np = z$npr, jnd = jnd,
        subregcoef = subregcoef, rvar = z$rv, aicf = z$aicf, respns = array(z$ei, dim = c(d, d)),
#        regcoef = array(z$bi, dim = c(d, d, z$morder)), matv = array(z$matv, dim = c(d, d)),
        matv = array(z$matv, dim = c(d, d)),
        morder = z$morder, arcoef = array(z$arcoef, dim = c(d, d, z$morder)), aicsum = z$aics)
    return(mulmar.out)
}


perars <-
function (y, ni, lag=NULL, ksw=0)
{
    n <- length(y)
    if( is.null(lag) ) lag <- as.integer(2*sqrt(ni))
    lag1 <- lag+1

    mean <-0
    var <-0
    np  <- rep(0,ni)                        # number of parameter
    jnd <- array(0, dim=c(lag1*ni+ksw,ni))  # specification of i-th regressor (i=1,...,ip)
    a <- array(0, dim=c(lag1*ni+ksw,ni))    # regression coefficients
    aic <- rep(0,ni)
    b <- array(0, dim=c(ni,ni,lag))  # AR-coefficient matrices
    v <- array(0, dim=c(ni,ni))      # innovation variance matrix
    c <- rep(0,ni)                   # constant vector
    osd <- rep(0,ni)                 # residual variances
    morder <- 0                      # order of the maice model

    z <- .C("perars",
	as.double(y),
	as.integer(n),
	as.integer(ni),
	as.integer(lag),
	as.integer(ksw),
	mean = as.double(mean),
	var = as.double(var),
	np = as.integer(np),
	jnd = as.integer(jnd),
	a = as.double(a),
	aic = as.double(aic),
	b = as.double(b),
	v = as.double(v),
	c = as.double(c),
	osd = as.double(osd),
	morder = as.integer(morder) )

    ind <- array(z$jnd, dim=c(lag1*ni+ksw,ni))
    jnd <- list()
    for( i in 1:ni )  jnd[[i]] <- ind[1:z$np[i],i]

    a <- array(z$a, dim=c(lag1*ni+ksw,ni))
    regcoef <- list()
    for( i in 1:ni )  regcoef[[i]] <- a[1:z$np[i],i]

    perars.out <- list( mean=z$mean, var=z$var, ord=jnd, regcoef=regcoef, rvar=z$osd, np=z$np, aic=z$aic,
                        v=array(z$v, dim=c(ni,ni)), arcoef=array(z$b, dim=c(ni,ni,z$morder)), const=z$c, morder=z$morder )

    class( perars.out ) <- "perars"
    return( perars.out )

}

print.perars <- function(x, ...)
{
  cat(sprintf("\n\n MEAN = %f\n", x$mean))
  cat(sprintf(" VARIANCE = %f\n", x$var))

  ni <- nrow(x$v)
  for( i in 1:ni ) {
    cat(sprintf("\n REGRESSION MODEL FOR THE REGRESSAND I = %i  ............\n", i))
    cat("\n SUBSET\t REGRESSION COEFFICIENTS\n")
    cat("   J\t\t  A(J)\n")
    sorder <- length(x$ord[[i]])
    for( j in 1:sorder )  cat(sprintf("   %i\t\t%f\n", x$ord[[i]][j],x$regcoef[[i]][j]))
    cat(sprintf("\n RVAR = RESIDUAL VARIANCE = %f\n", x$rvar[i]))
    cat(sprintf(" NP = NUMBER OF PARAMETER = %i\n", x$np[i]))
    cat(sprintf(" AIC = N*LOG(RVAR) + 2*NP = %f\n", x$aic[i]))
  }

  cat("\n\nmatrix of regression coefficients\n")
  print(x$v)
  cat("\n\nregression coefficients within the present period\n")
  print(x$arcoef)
  cat("\n\nconstants within the regression models\n")
  print(x$const)
}


unibar <-
function (y, ar.order=NULL, plot=TRUE)
{
    n <- length(y)
    if( is.null(ar.order) ) ar.order <- as.integer(2*sqrt(n))
    ar.order1 <- ar.order+1

    mean <- 0
    var <- 0
    v  <- rep(0,ar.order1)
    aic <- rep(0,ar.order1)
    daic <- rep(0,ar.order1)
    m <- 0
    aicm <- 0
    v.maice <- 0
    pa <- rep(0,ar.order)     # partial autocorrelation coefficients (AR-model)
    bw <- rep(0,ar.order1)    # Bayesian Weight
    sbw <- rep(0,ar.order)    # integrated Bayesian Weights
    pab <- rep(0,ar.order)    # partial autocorrelation coefficients (Bayesian model)
    aicb <- 0                 # AIC of Bayesian model
    vb <- 0                   # innovation variance of Bayesian model
    np <- 0                   # equivalent number of parameters
    a <- rep(0,ar.order)      # ar-coefficients (Bayesian model)
    pxx <- rep(0,121)         # power spectrum

    z <- .C("unibar",
	as.double(y),
	as.integer(n),
	as.integer(ar.order),
	mean = as.double(mean),
	var = as.double(var),
	v = as.double(v),
	aic = as.double(aic),
	daic = as.double(daic),
	m = as.integer(m),
	aicm = as.double(aicm),
	v.maice = as.double(v.maice),
	pa = as.double(pa),
	bw = as.double(bw),
	sbw = as.double(sbw),
	pab = as.double(pab),
	aicb = as.double(aicb),
	vb = as.double(vb),
	np = as.double(np),
	a = as.double(a), 
	pxx = as.double(pxx) )

    if( plot == TRUE ) {
      par(mfrow=c(3,1))
      plot((0:ar.order), z$daic, ylim=c(0,40), type="l", xlab="Lag", ylab="AIC(M)-AICMIN (Truncated at 40.0)")
      abline(h=0, lty=1)
      plot(z$pa, type="h", xlab="Lag", ylab="Partial autocorrelation")
      abline(h=0, lty=1)
#      abline(h=1/sqrt(n), lty=3)
#      abline(h=-1/sqrt(n),lty=3)
      abline(h=2/sqrt(n), lty=3)
      abline(h=-2/sqrt(n), lty=3)
      x <- rep(0,121)
      for( i in 1:121 ) x[i] <- (i-1)/240
      plot(x, z$pxx, type="l", xlab="Frequency", ylab="Power Spectral Density")
      par(mfrow=c(1,1))
    }

    unibar.out <- list( mean=z$mean, var=z$var, v=z$v, aic=z$aic, aicmin=z$aicm, daic=z$daic, 
			order.maice=z$m, v.maice=z$v.maice, pacoef=z$pa, bweight=z$bw[2:(ar.order1)], integra.bweight=z$sbw,
			v.bay=z$vb, aic.bay=z$aicb, np=z$np, pacoef.bay=z$pab, arcoef=z$a, pspec=z$pxx )
    return( unibar.out )
}


unimar <-
function (y, max.order=NULL, plot=FALSE, tmp.file=NULL)
{
    n <- length(y)
    if( is.null(max.order) ) max.order <- as.integer(2*sqrt(n))  # upper limit of AR-order
    morder <- max.order
    if( is.null(tmp.file) )  tmp.file <- " "

    mean <- 0
    var <- 0
    v  <- rep(0,morder+1)    # estimate of the innovation variance
    aic <- rep(0,morder+1)   # AIC
    daic <- rep(0,morder+1)  # AIC(M)-AICM
    m <- 0
    aicm <- 0                # minimum AIC
    v.maice <- 0             # innovation variance attained at m
    a <- rep(0,morder)       # AR-coefficients

    z <- .C("unimar",
	as.double(y),
	as.integer(n),
	as.integer(morder),
	mean = as.double(mean),
	var = as.double(var),
	v = as.double(v),
	aic = as.double(aic),
	daic = as.double(daic),
	m = as.integer(m),
	aicm = as.double(aicm),
	v.maice = as.double(v.maice),
	a = as.double(a), 
	as.character(tmp.file) )

    if( plot == TRUE ) {
      plot((0:morder), z$daic, ylim=c(0,40), type="l", xlab="Lag", ylab="AIC(M)-AICMIN (Truncated at 40.0)")
      abline(h=0, lty=1)
    }

    unimar.out <- list( mean=z$mean, var=z$var, v=z$v, aic=z$aic, aicmin=z$aicm, daic=z$daic, order.maice=z$m,
			v.maice=z$v.maice, arcoef=z$a[1:z$m] )
    return( unimar.out )
}


xsarma <-
function (y, arcoefi, macoefi)
{
    n <- length(y)
    arcoefi <- -arcoefi   # Initial estimates of AR coefficients
    macoefi <- -macoefi   # Initial estimates of MA coefficients
    p <- length(arcoefi)  # AR-ORDER
    q <- length(macoefi)  # MA-ORDER
    p01 <- c(arcoefi,macoefi)    # INITIAL ESTIMATES OF AR- AND  MA-COEFFICIENTS
    g1  <- rep(0,p+q)     # INITIAL GRADIENT
    tl1 <- 0              # INITIAL (-2)LOG LIKELIHOOD
    p02 <- rep(0,p+q)     # AR- AND MA-COEFFICIENTS
    g2  <- rep(0,p+q)     # FINAL GRADIENT
    alph.ar <- rep(0,p)   # FINAL ALPH (AR-PART)
    alph.ma <- rep(0,q)   #            (MA-PART)
    tl2 <- 0              # FINAL (-2)LOG LIKELIHOOD
    sigma2 <- 0           # WHITE NOISE VARIANCE

    z <- .C( "xsarma",
	as.double(y),
	as.integer(n),
	as.integer(p),
	as.integer(q),
	p01 = as.double(p01),
	g1 = as.double(g1),
	tl1 = as.double(tl1),
	p02 = as.double(p02),
	g2 = as.double(g2),
	alph.ar = as.double(alph.ar),
	alph.ma = as.double(alph.ma),
	tl2 = as.double(tl2), 
	sigma2 = as.double(sigma2) )

    p02 <- -z$p02

    xsarma.out <- list( gradi=z$g1, lkhoodi=z$tl1, arcoef=p02[1:p], macoef=p02[(p+1):(p+q)],
			grad=z$g2, alph.ar=z$alph.ar, alph.ma=z$alph.ma, lkhood=z$tl2, wnoise.var=z$sigma2)
    return( xsarma.out )
}
#####   DECOMP     #####

decomp <- function(y, trend.order=2, ar.order=2, frequency=12, seasonal.order=1, log=FALSE, trade=FALSE,
	    idif=1, year=1980, month=1, imiss=1, omax=99999.9, plot=TRUE)
{
    m1 <- trend.order
    m2 <- ar.order
    ilog <- 0
    if( log == TRUE ) ilog <- 1
    itrade <- 0
    if( trade == TRUE ) itrade <- 1

    n <- length(y)
    ipar <- rep(0, 9)
    ipar[1] <- trend.order
    ipar[2] <- ar.order
    ipar[3] <- frequency
    ipar[4] <- seasonal.order
    ipar[5] <- ilog
    ipar[6] <- itrade
    ipar[7] <- idif
    ipar[8] <- year
    ipar[9] <- month

    trend <- rep(0, n)
    seasonal <- rep(0, n)
    ar <- rep(0, n)
    trad <- rep(0, n)
    noise <- rep(0, n)
#    para <- rep(0, 13+m2)
    para <- rep(0, 26)

    z <- .C("decomp",
             as.double(y),
             as.integer(n),
             as.integer(ipar),
             trend = as.double(trend),
             seasonal = as.double(seasonal),
             ar = as.double(ar),
             trad = as.double(trad),
             noise = as.double(noise),
             para = as.double(para),
             as.integer(imiss),
             as.double(omax))	

    aic=z$para[1]
    lkhd=z$para[2]
    sigma2=z$para[3]
    tau1=z$para[4]
    tau2=z$para[5]
    tau3=z$para[6]
    arcoef=z$para[7:(6+m2)]
    tdf=z$para[(7+m2):(13+m2)]

    if( plot == TRUE ) {
      if(ilog == 1) y <- log(y)
      nc <- 2
      if( ar.order != 0 ) nc <- nc+1
      if( seasonal.order != 0 ) nc <- nc+1
      if( itrade == 1 ) nc <- nc+1
      if( nc > 3 ) par(mfrow=c((nc+1)/2,2))
      if( nc <= 3 ) par(mfrow=c(nc,1))
      matD <- array(c(y,z$trend),dim=c(n,2))
      matplot(matD, pch=0, type="l", col=1:2, main="Original and Trend", xlab="", ylab="")
      ymax <- max(z$season, z$ar, z$trad)
      ymin <- min(z$season, z$ar, z$trad)
      my <- max(ymax, abs(ymin))*1.5
      if( seasonal.order != 0 ) plot(z$seasonal, type="l", main= "Seasonal", xlab="", ylab="", ylim=c(-my,my))
      plot(z$noise, type="l", main= "Noise", xlab="", ylab="", ylim=c(-my,my))
      if( ar.order != 0 ) plot(z$ar, type="l", main="AR component", xlab="", ylab="", ylim=c(-my,my))
      if( itrade == 1) plot(z$trad, type="l", main="Trading Day Effect", xlab="", ylab="", ylim=c(-my,my))
      par(mfrow=c(1,1))
    }

    decomp.out <- list(trend=z$trend, seasonal=z$seasonal, ar=z$ar, trad=z$trad, noise=z$noise,
                       aic=aic, lkhd=lkhd, sigma2=sigma2, tau1=tau1, tau2=tau2, tau3=tau3, arcoef=arcoef, tdf=tdf)
    return(decomp.out)
}

#####   IWANAMI     #####

armaimp <- function( arcoef=NULL, macoef=NULL, v, n=1000, lag=NULL, nf=200, plot=TRUE )     # PROGRAM 6.1
{
    if( is.null(arcoef) ) {       # AR coefficients
      arorder <- 0
      arcoef <- 0.0
    } else {
      arorder <- length(arcoef)   # AR order
    }
    if( is.null(macoef) ) {       # MA coefficients
      maorder <- 0
      macoef <- 0.0
    } else {
      maorder <- length(macoef)   # MA order
    }
#   v                             # innovation variance
#   n                             # original data length
    if( is.null(lag) )  lag <- as.integer(2*sqrt(n))    # maximum lag of autocovariance function
#   nf                            # number of frequencies in evaluating spectrum

    g <- rep(0,(lag+1))           # impulse response function
    acov <- rep(0,(lag+1))        # autocovariance function
    parcor <- rep(0,lag)          # partial autocorrelation coefficient
    spec <- rep(0,(nf+1))         # power spectrum
    roota <- array(0, dim=c(arorder,2))   # characteristic roots of AR operator
    rootb <- array(0, dim=c(maorder,2))   # characteristic roots of MA operator
    ier <- 0
    jer <- 0

    z <- .C("arma",
	     as.integer(arorder),
	     as.integer(maorder),
	     as.double(arcoef),
	     as.double(macoef),
	     as.double(v),
	     as.integer(n),
	     as.integer(lag),
	     as.integer(nf),
	     g = as.double(g),
	     acov = as.double(acov),
	     parcor = as.double(parcor),
	     spec = as.double(spec),
	     roota = as.double(roota),
	     rootb = as.double(rootb),
             ier = as.integer(ier),
             jer = as.integer(jer))

    impuls <- z$g
    acov <- z$acov
    parcor <- z$parcor
    spec <- z$spec
    roota <- array(z$roota, dim=c(arorder,2))
    rootb <- array(z$rootb, dim=c(maorder,2))

    croot.ar <- list()
    croot.ma <- list()
    if( arorder != 0 ) {
    for ( i in 1:arorder ) {
      re <- roota[i,1]
      im <- roota[i,2]
      amp <- sqrt(re**2 + im**2)
      atan <- atan2(im,re)
      croot.ar[[i]] <- list(real=re, image=im, amp=amp, atan=atan, degree=atan*57.29577951) } }
    if( maorder != 0 ) {
    for ( i in 1:maorder ) {
      re <- rootb[i,1]
      im <- rootb[i,2]
      amp <- sqrt(re**2 + im**2)
      atan <- atan2(im,re)
      croot.ma[[i]] <- list(real=re, image=im, amp=amp, atan=atan, degree=atan*57.29577951) } }

    if( z$ier == 1 ) cat(" ***** ERROR : MATRIX WITH ZERO ROW IN DECOMPOSE\n" )
    if( z$ier == 2 ) cat(" ***** ERROR : SINGULAR MATRIX IN DECOMPOSE.ZERO DIDIVIDE IN SOLVE\n" )
    if( z$ier == 3 ) cat(" ***** ERROR : CONVERGENCE IN IMPRUV.MATRIX IS NEARLY SINGULAR\n" )
    if( z$jer == 1 ) cat(" ***** ERROR : NON-CONVERGENCE AT POLYRT\n" )

    if( plot == TRUE ) {
      par(mfrow=c(3,2))
      x <- c(0:lag)
      ymin <- as.integer(min(impuls)-1)
      ymax <- as.integer(max(impuls)+1)
      plot(x,impuls, type='l', xlim=c(0,lag), ylim=c(ymin,ymax), xlab='', ylab='')
      par(new=TRUE)
      plot(x,impuls, type='h', xlim=c(0,lag), ylim=c(ymin,ymax), xlab='lag', ylab='impulse')

      ymin <- as.integer(min(acov)-1)
      ymax <- as.integer(max(acov)+1)
      plot(x, acov, type='l', xlim=c(0,lag), ylim=c(ymin,ymax), xlab='', ylab='')
      par(new=TRUE)
      plot(x, acov, type='h', xlim=c(0,lag), ylim=c(ymin,ymax), xlab='lag', ylab='autocovariance')

      ymin <- as.integer(min(parcor)-1)
      ymax <- as.integer(max(parcor)+1)
      plot(parcor, type='l', xlim=c(0,lag), ylim=c(ymin,ymax), xlab='', ylab='')
      par(new=TRUE)
      plot(parcor, type='h', xlim=c(0,lag), ylim=c(ymin,ymax), xlab='lag', ylab='parcor')

      k1 <- length(spec)
      k <- k1-1
      x <- rep(0,k1)
      for( i in 1:k1 ) x[i] <- (i-1)/(2*k)
      plot(x, spec, type="l", xlab="frequency", ylab="log spectrum", )

      if( z$jer == 0 ) {
        par(pty="s")
        plot(x=c(-1.0,1.0), y=c(0,0), type='l', xlim=c(-1.1,1.1), ylim=c(-1.1,1.1), xlab="ARMA characteristic roots\n(square:AR, triangle:MA)", ylab="", axes=FALSE)
        par(new=TRUE)
        plot(x=c(0,0), y=c(-1.0,1.0), type='l', xlim=c(-1.1,1.1), ylim=c(-1.1,1.1), xlab="", ylab="", axes=FALSE)
        symbols(x=0, y=0, circles=1, xlim=c(-1.0,1.0), ylim=c(-1.0,1.0), inches=FALSE, add=TRUE)
        if( arorder != 0 ) {
        for( i in 1:arorder ) {     # characteristic roots of AR operator
          x1 <- croot.ar[[i]]$re
          y1 <- croot.ar[[i]]$im
          points(x1,y1,pch=19) }}
        if( maorder != 0 ) {
        for( i in 1:maorder ) {     # characteristic roots of MA operator
          x2 <- croot.ma[[i]]$re
          y2 <- croot.ma[[i]]$im
          points(x2,y2,pch=24) } } }
      par(mfrow=c(1,1), pty="m")
    }

    if( arorder!=0 && maorder!=0 ) armaimp.out <- list( impuls=impuls, acov=acov, parcor=parcor, spec=spec, croot.ar=croot.ar, croot.ma=croot.ma )
    if( arorder==0 && maorder!=0 ) armaimp.out <- list( impuls=impuls, acov=acov, parcor=parcor, spec=spec, croot.ma=croot.ma )
    if( arorder!=0 && maorder==0 ) armaimp.out <- list( impuls=impuls, acov=acov, parcor=parcor, spec=spec, croot.ar=croot.ar )

    return(armaimp.out)
}


tvvar <- function(y, trend.order, tau20, delta, plot=TRUE)     # PROGRAM 13.1
{
    n <- length(y)             # length of data
    n1 <- n/2
    m <- trend.order           # trend order
    iopt <- 1                  # search method
    if( is.null(tau20) || is.null(delta) ) iopt <- 0
    if( is.null(tau20) ) tau20 <- 0    # initial estimate of TAU2
    if( is.null(delta) ) delta <- 0    # search width

    tvvar <- rep(0,n1)
    normdat <- rep(0,n)
    y1 <- rep(0,n1)
    trend <- array(0, dim=c(n1,3))
    noise <- rep(0,n1)
    taumax <- 0
    sig2m <- 0
    ffmax <- 0
    aic <- 0

    z <- .C("tvvar",
	     as.double(y),
	     as.integer(n),
	     as.integer(m),
	     as.double(tau20),
	     as.integer(iopt),
	     as.double(delta),
	     tvvar = as.double(tvvar),
	     normdat = as.double(normdat),
	     y1 = as.double(y1),
	     n1 = as.integer(n1),
	     trend = as.double(trend),
	     noise = as.double(noise),
	     taumax = as.double(taumax),
	     sig2m = as.double(sig2m),
	     ffmax = as.double(ffmax),
	     aic = as.double(aic))

    normdat <- z$normdat
    ts <- z$y1
    trend <- array(z$trend,dim=c(z$n1,3))
    noise <- z$noise

    if( plot == TRUE ) {
      par(mfcol=c(4,1))
      ymin <- as.integer(min(normdat)-1)
      ymax <- as.integer(max(normdat)+1)
      plot(normdat, type='h', ylim=c(ymin,ymax), xlab="m", ylab="normalized data")
      par(new=TRUE)
      abline(h=0)

      ymin <- as.integer(min(ts)-1)
      ymax <- as.integer(max(ts)+1)
      plot(ts, type='l', xlab="m", ylab="s(m)=y(2m-1)**2 + y(2m)**2", ylim=c(ymin,ymax))

      plot(trend[,1], type="l", ylim=c(ymin,ymax), xlab="m", ylab="trend  t(m)") 
      par(new=TRUE) 
      plot(trend[,2], type="l", ylim=c(ymin,ymax), xlab="", ylab="", col=2)
      par(new=TRUE)
      plot(trend[,3], type="l", ylim=c(ymin,ymax), xlab="", ylab="") 

      ymin <- as.integer(min(noise)-1)
      ymax <- as.integer(max(noise)+1)
      plot(noise,type='h', ylim=c(ymin,ymax), xlab="m", ylab="noise")
      par(new=TRUE)
      abline(h=0)
      par(mfcol=c(1,1), new=FALSE)
    }

    tvvar.out <- list(tvvar=z$tvvar, normdat=normdat, ts=ts, trend=trend, noise=noise, tau2=z$taumax, sigma2=z$sig2m, lkhood=z$ffmax, aic=z$aic)
    return(tvvar.out)
}


tvar <-
function (y,ar.order,trend.order=2,span,outlier,tau20=NULL,delta=NULL,plot=TRUE)     # PROGRAM 13.2
{
#    y                    # original data
    n <- length(y)        # data length
#    ar.order             # AR order
#    trend.order          # Trend order
    if( trend.order!=1 && trend.order!=2 ) stop( " ***** ERROR : 'trend.order' is 1 or 2." )
#    span                 # local stationary span
#    outlier              # position of i-th outlier
    nout <- length(outlier)   # number of outliers
#    method               # search method
#    tau20                # initial variance of systen noise
#    delta                # delta for computing variance of system noise
    method <- 1
    if( is.null(tau20) || is.null(delta) ) method <- 0
    if( is.null(tau20) ) tau20 <- 0
    if( is.null(delta) ) delta <- 0
    
    nn <- n/span
    tau2 <- 0.0
    sigma2 <- 0.0
    lkhood <- 0.0
    aic <- 0.0
    arcoef <- array(0, dim=c(ar.order,nn))
    parcor <- array(0, dim=c(ar.order,nn))

    z <- .C("tvar",
	     as.double(y),
	     as.integer(n),
	     as.integer(ar.order),
	     as.integer(trend.order),
	     as.integer(span),
	     as.integer(method),
	     as.integer(nout),
	     as.integer(outlier),
	     as.double(tau20),
	     as.double(delta),
	     tau2 = as.double(tau2),
	     sigma2 = as.double(sigma2),
	     lkhood = as.double(lkhood),
	     aic = as.double(aic),
	     arcoef = as.double(arcoef),
	     parcor = as.double(parcor))

    sigma2 <- z$sigma2
    arcoef <- array(z$arcoef, dim=c(ar.order,nn))
    parcor <- array(z$parcor, dim=c(ar.order,nn))

    if( plot == TRUE ) {
      x <- span
      for( i in 2:nn ) x <- c(x, i*span)
      if( ar.order < 6 ) par(mfrow=c(ar.order,1))
      if( ar.order > 5 ) par(mfrow=c(5,1))
      for( i in 1:ar.order ) {
          if( (i%%5 == 1) & i > 1 ) par(ask=TRUE)
          plot(x, parcor[i,], type="l", xlab="", ylim=c(-1.0,1.0),ylab=paste("parcor( i=",i,")")) } 
      par(mfrow=c(1,1)) }

#  PROGRAM 13.3  TVSPC   ...  CHANGING SPECTRUM ...  

    nf <- 200          # parameter ( number of frequencies )
    ivar <- 0          # =1: for variance correction
    var <- rep(0, n)   # time varying variance
    spec <- array(0, dim=c(nf+1,nn))

    z1 <- .C("tvspc",
	     as.integer(nn),
	     as.integer(ar.order),
	     as.integer(span),
	     as.integer(nf),
	     as.integer(ivar),
	     as.double(sigma2),
	     as.double(arcoef),
	     as.double(var),
	     spec = as.double(spec))

    spec <- array(z1$spec, dim=c(nf+1,nn))

    if( plot == TRUE ) {
      par(ask=TRUE)
      x <- seq(0, 0.5, length=nf+1)
      y <- seq(0, n, length=nn)
      zmin <- as.integer(min(spec)-1)
      zmax <- as.integer(max(spec)+1)
      persp(x,y,z=spec,zlim=c(zmin,zmax),theta=10,phi=20,expand=0.5,col="lightblue",xlab="f",zlab="log p(f)",ticktype="detail")
      par(ask=FALSE) }

    tvar.out <- list( tau2=z$tau2, sigma2=sigma2, lkhood=z$lkhood, aic=z$aic, arcoef=arcoef, parcor=parcor, spec=spec )
    return( tvar.out )
}


ngsmth <- function(y,noisev=2,tau2,bv=1.0,noisew=1,sig2,bw=1.0,initd=1,k=200,plot=TRUE)     # PROGRAM 14.1
{
#    y                    # original data
    n <- length(y)        # data length
#   noisev                # type of system noise density (0,1,2,3)
                          # 1: Gaussian (normal) / 2: Pearson family / 3: two-sides exponential
    if( noisev!=0 && noisev!=1 && noisev!=2 && noisev!=3 ) stop( " ***** ERROR : 'noisev' is numeric in {1,2,3}" )
#   tau2                  # variance of dispersion of system noise
#   bv                    # shape parameter of system noise (for noisev=2)
#   noisew                # type of observation noise density (0,1,2,3,4)
                          # 1: Gaussian (normal) / 2: Pearson family / 3: two-sided exponential / 4:double exponential
    if( noisew!=0 && noisew!=1 && noisew!=2 && noisew!=3 && noisew!=4 ) stop( " ***** ERROR : 'noisew' is numeric in {1,2,3,4}" )
#   sig2                  # variance of dispersion of observation noise
#   bw                    # shape parameter of observation noise (for noisew=2)
#   initd                 # type of density function
                          # 0: two-sided exponential / 1: Gaussian (normal) / 2: uniform
    if( initd!=0 && initd!=1 && initd!=2 ) stop( " ***** ERROR : 'initd' is numeric in {0,1,2}" )
#   k                     # number of intervals

    ns <- 1
    nfe <- n
    npe <- n
    k1 <- k+1

    trend <- array(0, dim=c(npe,7))      # trend
    lkhood <- 0.0                        # log-likelihood
    smt <- array(0, dim=c(k1,npe))       # smoothed density
    loc <- rep(0,npe)                    # location of the center of the interval

    z <- .C("ngsmth",
	     as.double(y),
	     as.integer(n),
	     as.integer(noisev),
	     as.double(tau2),
	     as.double(bv),
	     as.integer(noisew),
	     as.double(sig2),
	     as.double(bw),
	     as.integer(initd),
	     trend = as.double(trend),
	     smt = as.single(smt),
	     lkhood = as.double(lkhood),
	     as.integer(ns),
	     as.integer(nfe),
	     as.integer(npe),
	     as.integer(k1))

    trend <- array(z$trend, dim=c(npe,7))
    smt <- array(z$smt, dim=c(k1,npe))

    if( plot == TRUE ) {
      ymin <- 0
      ymax <- 0
      for( i in 1:7 )
        for( j in 1:npe ) {
          if( is.na(trend[j,i]) == FALSE ) {
          if( trend[j,i] < ymin ) ymin <- trend[j,i]
          if( trend[j,i] > ymax ) ymax <- trend[j,i] }}
      ymin <- as.integer(ymin-1)
      ymax <- as.integer(ymax+1)
      for( i in 1:7 ) {
        if( i != 1 ) par(new=TRUE) 
        if( i == 4 ) {
          plot(trend[,i], ylim=c(ymin,ymax), type="l", xlab="n",ylab="trend  tn", col=2)
        } else {
          plot(trend[,i], ylim=c(ymin,ymax), type="l", xlab="",ylab="") }
      }
    }

#  plot smoothed density  : subroutine post3d

    ndif <- 1
    if( n >= 100 ) ndif <- 2
    if( n >= 200 ) ndif <- 4
    if( n >= 300 ) ndif <- 6
    if( n >= 500 ) ndif <- as.integer(n/50)
    n0 <- as.integer(ndif/2)+1
    nn <- as.integer((n-n0)/ndif)+1

    ss <- array(0, dim=c(k+1,nn))
    jj <- 0
    for( j in n0:npe ) {
      if( (j-n0)%%ndif == 0 ) {
	ss[,jj] <- smt[,j]
        jj <- jj+1
      }
    }

    if( plot == TRUE ) {
      par(ask=TRUE)
      xs <- min(y)
      xe <- max(y)
      x <- seq(xs, xe, length=k+1)
      y <- seq(1, n, length=nn)
      zmin <- 0
      zmax <- as.integer(max(ss)+1)
      persp(x,y,ss,zlim=c(zmin,zmax),theta=30,phi=20,expand=0.5,col="lightblue",xlab="tn",ylab="n",zlab="p(tn)",ticktype="detaile") 

      par(ask=FALSE)
    }

    ngsmth.out <- list( trend=trend, smt=smt, lkhood=z$lkhood )
    return( ngsmth.out )

}

tsmooth <- function( y, f, g, h, q, r, x0=NULL, v0=NULL, filter.end=NULL, predict.end=NULL, outmin=-10.0e+30, outmax=10.0e+30, missed=NULL, np=NULL, plot=FALSE)     # PROGRAM 9.1
{
#   y                    # time series
    yy <- y
    if( is.matrix(yy) == FALSE ) yy <- matrix(yy, length(y), 1)
    n <- dim(yy)[1]      # data length
    l <- dim(yy)[2]      # dimension

    ff <- f
    if( is.matrix(ff) == FALSE ) ff <- as.matrix(ff)
    m <- dim(ff)[1]      # dimension of the state vector

    qq <- q
    if( is.matrix(qq) == FALSE ) qq <- as.matrix(qq)
    k <- dim(qq)[1]      # dimension of the system noise
    
    rr <- r
    if( is.matrix(rr) == FALSE ) rr <- matrix(rr, l, l)

    gg <- g
    if( is.matrix(gg) == FALSE ) gg <- matrix(gg, m, k)

    hh <- h
    if( is.matrix(hh) == FALSE ) hh <- matrix(hh, l, m)

    if( is.null(x0) ) x0 <- rep(0.0e0, m)
    if( is.null(v0) ) {
      v0 <- matrix( 0.0e0, m, m )
      for( i in 1:m ) v0[i,i] <- 10.0e-5
    }

   if( is.null(filter.end) )   filter.end <- n      # end point of filtering
   if( is.null(predict.end) )  predict.end <- n     # end point of prediction
#   outmin               # lower limits of observations
#   outmax               # upper limits of observations

    startp <- missed
    if( is.null(startp) || is.null(np) ) {
      nmiss <- 0
      startp <- 0
      np <- 0
    } else {
      nmiss <- min(length(startp),length(np))   # number of missed intervals
      startp <- startp[1:nmiss]                 # start position of missed intervals
      np <- np[1:nmiss]                         # number of missed observations
    }
 
    npe <- predict.end
    xss <- array(0, dim=c(m,npe))     # mean vectors of the smoother
    vss <- array(0, dim=c(m,m,npe))   # covariance matrices of the smoother
    lkhood <- 0.0                     # log-likelihood
    aic <- 0.0                        # AIC

    z <- .C("smooth",
	     as.double(yy),
	     as.integer(n),
	     as.integer(l),
	     as.integer(m),
	     as.integer(k),
	     as.double(ff),
	     as.double(gg),
	     as.double(hh),
	     as.double(qq),
	     as.double(rr),
	     as.double(x0),
	     as.double(v0),
	     as.integer(filter.end),
	     as.integer(predict.end),
	     as.double(outmin),
	     as.double(outmax),
	     as.integer(nmiss),
	     as.integer(startp),
	     as.integer(np),
	     xss = as.double(xss),
	     vss = as.double(vss),
	     lkhood = as.double(lkhood),
	     aic = as.double(aic))

    xss <- array(z$xss, dim=c(m,npe))
    vss <- array(z$vss, dim=c(m,m,npe))
    cov.smooth <- array(,dim=c(m,npe))
    for( i in 1:m ) cov.smooth[i,] <- vss[i,i,]

    err <- array(0, dim=c(n,m,l))
    for( j in 1:l )
      for( ij in 1:m )
        for( i in 1:n ) err[i,ij,j] <- yy[i,j] - (mean(yy[,j]) + xss[ij,i])

    if( plot == TRUE ) {
      par(mfcol=c(m,1))
      ymin <- as.integer(min(xss)-1)
      ymax <- as.integer(max(xss)+1)
      plot(xss[1,], type="l", ylim=c(ymin,ymax), main=paste("Mean vectors of the smoother XSS(i,j) ( i=1:",m,", j=1:",npe,")"), xlab="i = 1", ylab="")
      if( m > 1 ) for( i in 2:m )
        plot(xss[i,], type="l", ylim=c(ymin,ymax), xlab=paste("i =",i), ylab="")
    }

    tsmooth.out <- list(mean.smooth=xss, cov.smooth=cov.smooth, esterr=err, lkhood=z$lkhood, aic=z$aic)
    return(tsmooth.out)
}

.noGenerics <- TRUE

options(warn.FPU=FALSE)
