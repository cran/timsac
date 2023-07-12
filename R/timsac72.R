#####   TIMSAC72   #####

autcor <- function (y, lag = NULL, plot = TRUE, lag_axis = TRUE)
{
    n <- length(y)
    if (is.null(lag))
        lag <- as.integer(2*sqrt(n))
    lag1 <- lag + 1

    z <- .Fortran(C_autcorf,
                  as.double(y),
                  as.integer(n),
                  acov = double(lag1),
                  acor = double(lag1),
                  as.integer(lag1),
                  mean = double(1))


    if (plot == TRUE) {
        plot((0:lag), z$acor, type = "h", ylab = "Autocorrelation",
             xlab = "Lag")
        if (lag_axis == TRUE)
            abline(h = 0, lty = 1)
    }

    autcor.out <- list(acov = z$acov, acor = z$acor, mean = z$mean)
    return(autcor.out)
}


fftcor <- function (y, lag = NULL, isw = 4, plot = TRUE, lag_axis = TRUE)
{ 
    ld <- nrow(y)                # length of data
    d <- ncol(y)                 # dimension of the observation vector
    x1 <- y[,1]                  # data of channel X
    y1 <- rep(0, ld)
    if (d == 2)
        y1 <- y[, 2]       # data of channel Y
    if (is.null(lag))           # maximum lag
        lag <- as.integer(2*sqrt(ld))
    lag1 <- lag + 1

#  n2p, n : definition    n=2**n2p >= ld+lag1 >= 2**(n2p-1)
    nd <- ld + lag1
    n2p <- 1
    n <- 2**n2p
    while((n-nd) < 0) {
	n2p <- n2p + 1
	n <- 2**n2p
    }

    acov <- array(0, dim=c(n, 2))      # auto covariance
    ccov21 <- rep(0, n)                # cross covariance
    ccov12 <- rep(0, n)                # cross covariance
    acor <- array(0, dim = c(lag1, 2))   # auto correlation
    ccor21 <- rep(0, lag1)             # cross correlation
    ccor12 <- rep(0, lag1)             # cross correlation
    mean <- rep(0, 2)                  # mean

    z <- .Fortran(C_fftcorf,
                  as.integer(ld),
                  as.integer(lag1),
                  as.integer(n),
                  as.integer(n2p),
                  as.integer(isw),
                  as.double(x1),
                  as.double(y1),
                  acov = double(2*n),
                  cov21 = double(n),
                  cov12 = double(n),
                  acor = double(2*lag1),
                  cor21 = double(lag1),
                  cor12 = double(lag1),
                  mean = double(2))

    acv <- array(z$acov, dim = c(n, 2))
    cna1 <- array(z$acor, dim = c(lag1, 2))
    if (isw == 1 ) {
        acov <-  acv[1:lag1, 1]
        cna1 <- cna1[1:lag1, 1]
        amean <- z$mean[1]
    } else {
        acov <- array(0, dim = c(lag1, 2))
        for(i in 1:lag1) {
            acov[i, 1] <- acv[i, 1]
            acov[i, 2] <- acv[i, 2]
        }
        amean <- z$mean
    }

    if(isw == 4) {
        fftcor.out  <- list(acov = acov, ccov12 = z$cov12[1:lag1],
                            ccov21 = z$cov21[1:lag1], acor = cna1,
                            ccor12 = z$cor12, ccor21 = z$cor21, mean = amean)
        if(plot == TRUE) {
            par(mfrow=c(2,1))
            plot((0:lag), z$cor12, type = "h", xlab = "Lag",
                 ylab = "Crosscorrelation ccor12")
            if(lag_axis == TRUE)
                abline(h = 0, lty = 1)
            plot((0:lag), z$cor21, type = "h", xlab = "Lag",
                 ylab = "Crosscorrelation ccor21")
            if(lag_axis == TRUE)
                abline(h = 0, lty = 1)
            par(mfrow = c(1, 1))
        }
    } else if(isw != 4) {
        fftcor.out  <- list(acov = acov, acor = cna1, mean = amean)
    }
    return(fftcor.out)
}


mulcor <- function (y, lag = NULL, plot = TRUE, lag_axis = TRUE)
{
    n <- nrow(y)        # length of data
    d <- ncol(y)        # dimension of the observation vector
    if (is.null(lag))
         lag <- as.integer(2*sqrt(n))  # maximum lag
    lag1 <- lag + 1

    z <- .Fortran(C_mulcorf,
                  as.double(y),
                  as.integer(n),
                  as.integer(d),
                  as.integer(lag1),
                  mean = double(d),
                  cov = double(lag1*d*d),
                  cor = double(lag1*d*d))

    cov = array(z$cov, dim = c(lag1, d, d))
    cor = array(z$cor, dim = c(lag1, d, d))

    if (plot == TRUE) {
        oldpar <- par(no.readonly = TRUE)
        x <- rep(0, lag1)
        for(i in 1:lag1)
            x[i] <- i - 1
        par(mfrow = c(d, d))
        for(j in 1:d)
            for(k in 1:d) {
                cy <- cor[, j, k]
                plot(x, cy, type = "h", xlab = "lag",
                     ylab = paste("cor [,",j,",",k,"]"), ylim = c(-1, 1))
                if(lag_axis == TRUE)
                    abline(h = 0, lty = 1)
            }
        par(oldpar)
    }

    mulcor.out <- list(cov = cov, cor = cor, mean = z$mean)
    class(mulcor.out) <- "mulcor"
    return(mulcor.out)

}


print.mulcor <- function(x, ...)
{
    lag <- dim(x$cov)[1] - 1
    ip <- dim(x$cov)[2]

    cat("\n Mean =")
    for(i in 1:ip)
        cat(sprintf(" %f", x$mean[i]))

    cat("\n\n Autocovariance\n")
    for(i in 1:ip) {
        cat(sprintf("\n lag\tcov[,%i,%i]\tNormalized\n", i, i))
        for(j in 0:lag)
            cat(sprintf(" %i\t%f\t%f\n", j, x$cov[j+1, i, i], x$cor[j+1, i, i]))
    }

    cat("\n\n Cross-covariance\n")
    for(i in 2:ip)
        for(j in 1:(i-1)) {
            cat(sprintf("\n lag\tcor[,%i,%i]\tNormalized", j, i))
            cat(sprintf("\tcor[,%i,%i]\tNormalized\n", i, j))
            for(k in 0:lag)
                cat(sprintf(" %i\t%f\t%f\t%f\t%f\n", k, x$cov[k+1, j, i],
                    x$cor[k+1, j, i], x$cov[k+1, i, j], x$cor[k+1, i, j]))
        }
}


auspec <-
function (y, lag=NULL, window = "Akaike", log = FALSE, plot = TRUE)
{
    if (window != "Akaike" && window != "Hanning")
        stop("allowed windows are 'Akaike' or 'Hanning'")

    n <- length(y)
    if (is.null(lag))
        lag <- as.integer(2*sqrt(n))  # maximum lag
    lag1 <- lag + 1

    spec1 <- rep(0, lag1)   # spectrum smoothing by window W1
    spec2 <- rep(0, lag1)   # spectrum smoothing by window W2
    stat <- rep(0, lag1)    # test statistics

    z1 <- autcor(y, lag, plot = FALSE)

    z <- .Fortran(C_auspecf,
                  as.integer(n),
                  as.integer(lag1),
                  as.double(z1$acov),
                  spec1 = double(lag1),
                  spec2 = double(lag1),
                  stat = double(lag1))

    if (window == "Akaike")
        spec <- z$spec2
    if (window == "Hanning")
        spec <- z$spec1

    if (plot == TRUE) {
        x <- rep(0, lag1)
        for(i in 1:lag1) 
            x[i] <- (i - 1) / (2 * lag)
        if (log == TRUE)
            plot(x, spec, type = "l", log = "y", xlab = "Frequency",
                 ylab = "Spectrum")
        if (log == FALSE)
            plot(x, spec, type = "l", xlab = "Frequency", ylab = "Spectrum")
        auspec.out <- list(stat=z$stat)
    } else {
        auspec.out <- list(spec = spec, stat = z$stat)
    }
    return(auspec.out)
}


mulspe <- 
function (y, lag = NULL, window = "Akaike", plot = TRUE, ...) 
{
    if (window != "Akaike" && window != "Hanning")
        stop("allowed windows are 'Akaike' or 'Hanning'")
 
    n <- nrow(y)
    d <- ncol(y)
    if (is.null(lag)) 
        lag <- as.integer(2 * sqrt(n))
    lag1 <- lag + 1

    z1 <- mulcor(y, lag, plot = FALSE)

    z <- .Fortran(C_mulspef,
                  as.integer(n),
                  as.integer(d),
                  as.integer(lag1), 
                  as.integer(lag1),
                  as.double(z1$cov),
                  spec1 = double(lag1*d*d),
                  spec2 = double(lag1*d*d),
                  stat = double(lag1*d),
                  coh1 = double(lag1*d*d),
                  coh2 = double(lag1*d*d))

    if (window == "Akaike") {
        spec <- array(z$spec2, dim = c(lag1, d, d))
        coh <- array(z$coh2, dim = c(lag1, d, d))
    }
    if (window == "Hanning") {
        spec <- array(z$spec1, dim = c(lag1, d, d))
        coh <- array(z$coh1, dim = c(lag1, d, d))
    }

    specmx <- array(0, dim = c(d, d, lag1))
    for(j in 1:d)
        for(k in 1:j)
            for(l in 1:lag1) {
              specmx[j, k, l] <- spec[l, j, k] + (0 + 1i) * spec[l, k, j]
              specmx[k, j, l] <- spec[l, j, k] - (0 + 1i) * spec[l, k, j]
            }

    mulspe.out <- list(spec = spec, specmx = specmx,
                       stat = array(z$stat, dim = c(lag1, d)), coh = coh)
    class(mulspe.out$specmx) <- "specmx"

    if (plot == TRUE) {
        plot.specmx(mulspe.out$specmx, ...)
        return(invisible(mulspe.out))
    } else {
        return(mulspe.out)
    }
}


print.mulspe <- function(x, ...)
{
    lag <- dim(x$spec)[1] - 1
    ip <- dim(x$spec)[2]

    for(i in 1:ip) {
        cat(sprintf("\nPower spectrum P(%i,%i)\tSignificance\n", i, i))
        for(j in 0:lag)
            cat(sprintf("%i\t%f\t%f\n", j, x$spec[j+1, i, i], x$stat[j+1, i]))
        if (i != 1) 
            for(j in 1:(i-1)) {
                cat(sprintf("\nCross spectrum P(%i,%i)\n", i, j))
                cat("I\tCo-spectrum\tQuad-spectrum\tSimple coherence\n")
                for(k in 0:lag)
                    cat(sprintf("%i\t%f\t%f\t%f\n", k, x$spec[k+1, i, j],
                        x$spec[k+1, j, i], x$coh[k+1, i, j]))
            }
    }
}


plot.specmx <- function(x, plot.scale = TRUE, ...)
{
    spec <- x
    lag1 <- dim(spec)[3]
    lag <- lag1 - 1
    d <- dim(spec)[1]

    par(mfrow = c(d, d))
    xx <- rep(0, lag1)
    for(i in 1:lag1)
        xx[i] <- (i - 1) / (2 * lag)

    dspec <- array(0, dim = c(d, d, lag1))
    for(j in 1:d)
        for(k in 1:d) {
            if (j >= k) {
                dspec[j, k, ] = Mod(spec[j, k, ])
            } else {
                dspec[j, k, ] = Arg(spec[j, k, ])
            }
        }
    for(j in 1:d)
        for(k in 1:d) {
            if (j >= k) {
                dspec[j, k, ] = Mod(spec[j, k, ])
            } else {
                dspec[j, k, ] = Arg(spec[j, k, ])
            }
        }

    if(plot.scale == TRUE) {
        mx  <- max(dspec[1, 1, ])
        for(j in 1:d)
            for(k in 1:j)
                if(mx < max(dspec[j, k, ]))
                    mx <- max(dspec[j, k, ])
        for(j in 1:d)
            for(k in 1:d)
                if (j >= k) {
                    ylabs = paste("AmpSp (", j, ",", k, ")")
                    plot(xx, dspec[j, k, ], type = "l", xlab = "Frequency",
                         ylab = ylabs, ylim=c(0,mx), ...)
                } else {
                    ylabs = paste("PhaseSp (", j, ",", k, ")")
                    plot(xx, dspec[j, k, ], type = "l", xlab = "Frequency",
                         ylab = ylabs, ylim=c(-pi,pi), ...)
                }
    } else {
        for(j in 1:d)
            for(k in 1:d)
                if (j >= k) {
                    ylabs = paste("AmpSp (", j, ",", k, ")")
                    plot(xx, dspec[j, k, ], type = "l", xlab = "Frequency",
                         ylab = ylabs, ...)
                } else {
                    ylabs = paste("PhaseSp (", j, ",", k, ")")
                    plot(xx, dspec[j, k, ], type = "l", xlab = "Frequency",
                         ylab = ylabs, ...)
                }
    }
    par(mfrow = c(1, 1))
}


sglfre <-
function (y, lag = NULL, invar, outvar)
{
    n <- nrow(y)        # length of data
    d <- ncol(y)        # dimension of the observation vector
    if (is.null(lag))
        lag <- as.integer(2*sqrt(n))
    lag1 <- lag + 1
    z1 <- mulspe(y, lag, window = "Akaike", plot = FALSE)

    z <- .Fortran(C_sglfref,
                  as.integer(invar),
                  as.integer(outvar),
                  as.integer(n),
                  as.integer(lag1),
                  as.integer(d),
                  as.double(z1$spec),
                  spec1 = double(lag1),
                  spec2 = double(lag1),
                  cspec = double(lag1),
                  qspec = double(lag1),
                  gain = double(lag1),
                  coh = double(lag1),
                  freqr = double(lag1),
                  freqi = double(lag1),
                  err = double(lag1),
                  phase = double(lag1))

    sglfre.out <- list(inspec = z$spec1, outspec = z$spec2, cspec = z$cspec,
                       qspec = z$qspec, gain = z$gain, coh = z$coh,
                       freqr = z$freqr, freqi = z$freqi, errstat = z$err,
                       phase = z$phase)
    return(sglfre.out)
}


mulfrf <-
function (y, lag = NULL, iovar = NULL)
{
    n <- nrow(y)        # length of data
    d <- ncol(y)        # dimension of the observation vector
    if (is.null(lag))
        lag <- as.integer(2*sqrt(n))  # maximum lag
    lag1 <- lag + 1
    if (is.null(iovar))
        iovar <- c(1:d)
    nv <- length(iovar) - 1          # number of input variables
    if (nv == d || nv > d)
        stop("number of input variable is smaller than d") 
    for(i in 1: (nv+1)) {
        if (iovar[i] < 1)
            stop("control variable is greater than or equal to 1")
        if (iovar[i] > d)
            stop("control variable is smaller than or equal to d")
    }

    z1 <- mulspe(y, lag, plot = FALSE)

    z <- .Fortran(C_mulfrff,
                  as.integer(nv),
                  as.integer(iovar),
                  as.integer(n),
                  as.integer(lag1),
                  as.integer(d),
                  as.double(z1$spec),
                  cospec = complex(d*d*lag1),
                  freqr = double(nv*lag1),
                  freqi = double(nv*lag1),
                  gain = double(nv*lag1),
                  phase = double(nv*lag1),
                  pcoh = double(nv*lag1),
                  errstat = double(nv*lag1),
                  mcoh = double(lag1))

    csp <- array(z$cospec, dim = c(d, d, lag1))
    fr <- array(z$freqr, dim = c(nv, lag1))
    fi <- array(z$freqi, dim = c(nv, lag1))
    g <- array(z$gain, dim = c(nv, lag1))
    ph <- array(z$phase, dim = c(nv, lag1))
    pch <- array(z$pcoh, dim = c(nv, lag1))
    err <- array(z$errstat, dim = c(nv, lag1))
    mch <- z$mcoh

    mulfrf.out <- list(cospec = csp, freqr = fr, freqi = fi, gain = g,
                       phase = ph, pcoh = pch, errstat = err, mcoh = mch)
    return(mulfrf.out)
}


fpeaut <- function (y, max.order = NULL)
{
    n <- length(y)
    if (is.null(max.order))
        max.order <- as.integer(2*sqrt(n))   # upper limit of model order
    morder <- max.order

    lag <- morder
    lag1 <- lag + 1
    z1 <- autcor(y, lag, plot = FALSE)
	sd <- z1$acov[1]
	cxx <- z1$acov[2:(lag1)]

    z <- .Fortran(C_fpeautf,
                  as.integer(morder),
                  as.integer(n),
                  as.double(sd),
                  as.double(cxx),
                  sig2 = double(morder),
                  fpe = double(morder),
                  rfpe = double(morder),
                  parcor = double(morder),
                  chi2 = double(morder),
                  ofpe = double(1),
                  fpem = double(1),
                  rfpem = double(1),
                  mo = integer(1),
                  sig2m = double(1),
                  a = double(morder*morder),
                  ao = double(morder))

    a <- array(z$a, dim = c(morder, morder))
    arcoef <- list()
    for(i in 1:morder)
        arcoef[[i]] <- a[1:i, i]

    mo <- z$mo
    fpeaut.out <- list(ordermin = mo, best.ar = arcoef[[mo]], sigma2m = z$sig2m,
                       fpemin = z$fpem, rfpemin = z$rfpem, ofpe = z$ofpe,
                       arcoef = arcoef, sigma2 = z$sig2, fpe = z$fpe,
                       rfpe = z$rfpe, parcor = z$parcor, chi2 = z$chi2)
    return(fpeaut.out)
}


fpec <- function (y, max.order = NULL, control = NULL, manip = NULL)
{
    n <- nrow(y)        # length of data
    d <- ncol(y)        # dimension of the observation vector
    if (is.null(max.order))
        max.order <- as.integer(2*sqrt(n))  # upper limit of model order
    ncon <- length(control)       # number of controlled variables
    nman <- length(manip)         # number of manipulated variables
    if (is.null(control)) {
        ncon <- d
        control <- c(1:d)
        inw <- control
    } else {
        inw <- control
        for(i in 1:ncon) {
            if (inw[i] < 1)
                stop("control variable is greater than or equal to 1")
            if (inw[i] > d)
                stop("control variable is smaller than or equal to d")
        }
    }
    if (is.null(manip)) {
        nman <- 0
    } else {
        for(i in 1:nman)  {
            if (manip[i] < 1)
                stop("manipulate variable is greater than or equal to 1")
            if (manip[i] > d)
                stop("manipulate variable is smaller than or equal to d")
        }
        inw <- c(inw, manip)
    }
    ip <- ncon + nman
    if (ip > d)
        stop("length of control and manipulate variables is smaller than or
 equal to d")

    morder <- max.order
    morder1 <- morder + 1
    z1 <- mulcor(y, morder, plot = FALSE)
    cov <- z1$cov[1:morder1, , ]

	
    z <- .Fortran(C_fpec7,
                  as.integer(n),
                  as.integer(morder),
                  as.integer(ncon),
                  as.integer(ip),
                  as.integer(d),
                  as.integer(inw),
                  as.double(cov),
                  ccv = double(morder1*ip*ip),
                  fpec = double(morder1),
                  rfpec = double(morder1),
                  aic = double(morder1),
                  mo = integer(1),
                  fpecm = double(1),
                  rfpecm = double(1),
                  aicm = double(1),
                  perr = double(ncon*ncon),
                  arcoef = double(morder*ncon*ip))

    cov <- array(z$ccv, dim = c(morder1, ip, ip))
    cov <- aperm(cov, c(2, 3, 1))
    arcoef <- array(z$arcoef, dim = c(morder, ncon, ip))
    arcoef <- arcoef[1:z$mo, , , drop=F]
    arcoef <- aperm(arcoef, c(2, 3, 1))

    fpec.out  <- list(cov = cov, fpec = z$fpec, rfpec = z$rfpec, aic = z$aic,
                      ordermin = z$mo, fpecmin = z$fpecm, rfpecmin = z$rfpecm,
                      aicmin = z$aicm, perr=array(z$perr, dim = c(ncon,ncon)),
                      arcoef = arcoef)
    class(fpec.out) <- "fpec"
    return(fpec.out)
}

print.fpec <- function(x, ...)
{
    m1 <- dim(x$arcoef)[1]
    m2 <- dim(x$arcoef)[2]
    m3 <- dim(x$arcoef)[3]
    lag <- dim(x$cov)[3] - 1

    cat("\n\nCovariance matrix\n")
    print(x$cov)

    cat("m\tFPEC\t\tRFPEC\t\tAIC\n")
    for(i in 0:lag)
        cat(sprintf("%i\t%f\t%f\t%f\n",i, x$fpec[i+1], x$rfpec[i+1], x$aic[i+1]))
    cat(sprintf("\n Minimum FPEC = %f\n", x$fpecmin))
    cat(sprintf(" Minimum RFPEC = %f\n", x$rfpecmin))
    cat(sprintf(" Minimum AIC = %f attained at m = %i\n", x$aicmin, x$ordermin))

    cat("\n\nPrediction error covariance matrix\n")
    print(x$perr)

    cat("\n\nAR coefficient matrix\n\n")
    print(x$arcoef)
}

mulnos <- function (y, max.order = NULL, control = NULL, manip = NULL, h)
{
    n <- nrow(y)                    # length of data
    d <- ncol(y)                    # dimension of the observation vector
    if (is.null(max.order))
        max.order <- as.integer(2*sqrt(n))  # upper limit of model order
    ncon <- length(control)       # number of controlled variables
    nman <- length(manip)         # number of manipulated variables
    n <- nrow(y)        # length of data
    d <- ncol(y)        # dimension of the observation vector
    if (is.null(max.order))
        max.order <- as.integer(2*sqrt(n))  # upper limit of model order
    ncon <- length(control)       # number of controlled variables
    nman <- length(manip)         # number of manipulated variables
    if (is.null(control)) {
        ncon <- d
        control <- c(1:d)
        inw <- control
    } else {
        inw <- control
        for(i in 1:ncon) {
            if (inw[i] < 1)
                stop("control variable is greater than or equal to 1")
            if (inw[i] > d)
                stop("control variable is smaller than or equal to d")
        }
    }
    if (is.null(manip)) {
        nman <- 0
    } else {
        for(i in 1:nman)  {
            if (manip[i] < 1)
                stop("manipulate variable is greater than or equal to 1")
            if (manip[i] > d)
                stop("manipulate variable is smaller than or equal to d")
        }
        inw <- c(inw, manip)
    }
    ip <- ncon + nman
    if (ip > d)
        stop("length of control and manipulate variables is smaller than or
 equal to d")

    z1 <- fpec(y, max.order, control, manip)
    arcoef <- aperm(z1$arcoef, c(3, 1, 2))
    h1 <- h + 1

    z <- .Fortran(C_mulnosf,
                  as.integer(h),
                  as.integer(z1$ordermin),
                  as.integer(ip),
                  as.double(z1$perr),
                  as.double(arcoef),
                  nsd = double(ip*ip),
                  drpc = double(ip*ip*h1),
                  irpc = double(ip*ip*h1))

    mulnos.out <- list(nperr = array(z$nsd, dim = c(ip, ip)),
                       diffr = array(z$drpc, dim = c(ip, ip, h1)),
                       integr = array(z$irpc, dim = c(ip, ip, h1)))
    return(mulnos.out)
}


raspec <-
function (h, var, arcoef = NULL, macoef = NULL, log = FALSE, plot = TRUE)
{
    l <- length(arcoef)
    k <- length(macoef)
    if (is.null(arcoef))
        arcoef <- 0   # coefficient matrix of autoregressive model
    if (is.null(macoef))
        macoef <- 0   # coefficient matrix of moving average model
    macoef <- -macoef
    h1 <- h + 1

    z <- .Fortran(C_raspecf,
                  as.integer(h),
                  as.integer(l),
                  as.integer(k),
                  as.double(var),
                  as.double(arcoef),
                  as.double(macoef),
                  rspec = double(h1))

    rspec <- z$rspec

    if (plot == TRUE) {
        x <- rep(0, h1)
        for(i in 1:h1)
            x[i] <- (i - 1) / (2 * h)
        if (log == TRUE)
            plot(x, rspec, type = "l", log = "y", xlab = "Frequency",
                 ylab = "Rational Spectrum")
        if (log == FALSE)
            plot(x, rspec, type = "l", xlab = "Frequency",
                 ylab = "Rational Spectrum")
        return(invisible(rspec))
    } else {
        return(raspec = rspec)
    }
}

mulrsp <-
function (h, d, cov, ar = NULL, ma = NULL, log = FALSE, plot = TRUE, ...)
{
    if (is.null(ar)) {
        l <- 0
        ar <- 0
    } else {
        l <- dim(ar)[3]
    }
    if (is.null(ma)) {
        k <- 0
        ma <- 0
    } else {
        k <- dim(ma)[3]
    }
    arcoef <- array(0, dim = c(d, d, l))
    macoef <- array(0, dim = c(d, d, k))
    if (l != 0) 
        arcoef <- aperm(ar, c(3, 1, 2))
    if (k != 0) 
        macoef <- aperm(ma, c(3, 1, 2))
    macoef <- -macoef
    h1 <- h + 1

    z <- .Fortran(C_mulrspf,
                  as.integer(h),
                  as.integer(l),
                  as.integer(d),
                  as.integer(k),
                  as.double(cov),
                  as.double(arcoef),
                  as.double(macoef),
                  rspec = complex(d*d*h1),
                  scoh = double(d*d*h1))

    rspec <- array(z$rspec, dim = c(d, d, h1))
    mulrsp.out <- list(rspec = rspec, scoh = array(z$scoh, dim = c(d, d, h1)))
    class(mulrsp.out$rspec) <- "specmx"

    if (plot == TRUE) {
        plot.specmx(mulrsp.out$rspec, ...)
        return(invisible(mulrsp.out))
    } else {
        return(mulrsp.out)
    }
}


optdes <- function (y, max.order = NULL, ns, q, r)
{
    n <- nrow(y)       # length of data
    d <- ncol(y)       # dimension of the observation vector
    if (is.null(max.order))
        max.order <- as.integer(2*sqrt(n))  # upper limit of model order
    ncon <- dim(q)[1]  # number of controlled variables
    nman <- dim(r)[1]  # number of manipulated variables
    control <- c(1:ncon)
    manip <- c((ncon+1):(ncon+nman))
    z1 <- fpec(y, max.order, control, manip)
    order <- z1$ordermin
    ao <- z1$arcoef
    ao <- aperm(ao, c(3, 1, 2))
    ao <- aperm(ao, c(2, 1, 3))
    dim(ao) <- c(ncon*order, (ncon+nman))
    a <- ao[, (1:ncon)]
    b <- ao[, (ncon+1):(ncon+nman)]

    z <- .Fortran(C_optdesf,
                  as.integer(ncon),
                  as.integer(nman),
                  as.integer(ns),
                  as.integer(order),
                  as.double(q),
                  as.double(r),
                  as.double(z1$perr),
                  as.double(a),
                  as.double(b),
                  gain = double(ncon*nman*order))

    optdes.out <- list(perr = z1$perr, trans = a, gamma = b,
                       gain = array(z$gain, dim = c(nman, order*ncon)))
    return(optdes.out)
}


optsim <-
function (y, max.order =NULL, ns, q, r, noise = NULL, len, plot = TRUE)
{
    n <- nrow(y)        # length of data
    d <- ncol(y)        # dimension of the observation vector
    if (is.null(max.order))
        max.order <- as.integer(2*sqrt(n))  # upper limit of model order
    ir <- dim(q)[1]     # number of controlled variables
    il <- dim(r)[1]     # number of manipulated variables

    z1 <- optdes(y, max.order, ns, q, r)
    trans <- z1$trans
    gamma <- z1$gamma
    gain <- z1$gain
    order <- dim(gain)[2] / ir
    if (is.null(noise))
        noise <- wnoise(len, z1$perr, plot = FALSE)     # white noise

    z <- .Fortran(C_optsimf,
                  as.integer(ns),
                  as.integer(order),
                  as.integer(ir),
                  as.integer(il),
                  as.double(trans),
                  as.double(gamma),
                  as.double(gain),
                  as.double(noise),
                  x = double(ir*ns),
                  y = double(il*ns),
                  xmean = double(ir),
                  ymean = double(il),
                  x2sum = double(ir),
                  y2sum = double(il),
                  x2mean = double(ir),
                  y2mean = double(il),
                  xvar = double(ir),
                  yvar = double(il))

    convar <- array(z$x, dim = c(ir, ns))
    manvar <- array(z$y, dim = c(il, ns))

    if (plot == TRUE) {
        nc <- max(ir, il)
        par(mfcol = c(2, nc))
        for(i in 1:ir) {
            plot(convar[i, ], type = "h", xlab = "Step",
                 ylab = paste("Controlled variables X(",i,")"))
            abline(h = 0, lty = 1)
        }
        for(i in 1:il) {
            plot(manvar[i, ], type = "h", xlab = "Step",
                 ylab = paste("Manipulated variables Y(",i,")"))
            abline(h = 0, lty = 1)
        }
        par(mfrow = c(1, 1))
        optsim.out <- list(trans = trans, gamma = gamma, gain = gain,
                           xmean = z$xmean, ymean = z$ymean, xvar = z$xvar,
                           yvar = z$yvar, x2sum = z$x2sum, y2sum = z$y2sum,
                           x2mean = z$x2mean, y2mean = z$y2mean)
    } else {
        optsim.out <- list(trans = trans, gamma = gamma, gain = gain,
                           convar = convar, manvar = manvar, xmean = z$xmean,
                           ymean = z$ymean, xvar = z$xvar, yvar = z$yvar,
                           x2sum = z$x2sum, y2sum = z$y2sum, x2mean = z$x2mean,
                           y2mean = z$y2mean)
    }
    return(optsim.out)
}


wnoise <- function (len, perr, plot = TRUE)
{
    ir <- 1                           # number of controlled variables
    if (is.array(perr))
        ir <- ncol(perr)

    z <- .Fortran(C_wnoisef,
                  as.integer(len),
                  as.integer(ir),
                  as.double(perr),
                  wn = double(len*ir))

    if (ir != 1)
        wnoise <- array(z$wn, dim = c(ir, len))
    if (ir == 1)
        wnoise <- z$wn[1:len]

    if (plot == TRUE) {
        par(mfcol = c(ir, 1))
        if (ir == 1) {
            plot(wnoise, type = "h", main = "White Noise", xlab = "n",
                 ylab = expression(y[n]))
            abline(h = 0, lty = 1)
        } else {
            plot(wnoise[1, ], type = "h",
                 main = expression(paste("White noise ", y[n])), xlab = "n",
                 ylab = "d=1")
            abline(h = 0, lty = 1)
            for(i in 2:ir) {
                plot(wnoise[i, ], type = "h", xlab = "n", ylab = paste("d=", i))
                abline(h=0, lty=1)
            }
        }
        par(mfrow = c(1, 1))
    } else {
        return(wnoise=wnoise)
    }
}

