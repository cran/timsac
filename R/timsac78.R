#####   TIMSAC78   #####

blocar <- function (y, max.order = NULL, span, plot = TRUE)
{
    n <- length(y)
    if (is.null(max.order))
        max.order <- as.integer(2*sqrt(n))  
    morder <- max.order

    if (span < 0)
        span <- n
    ns <- as.integer((n-morder+span-1)/span)

    z <- .Fortran(C_blocarf,
                  as.double(y),
                  as.integer(n),
                  as.integer(morder),
                  as.integer(span),
                  as.integer(ns),
                  mean = double(1),
                  var = double(1),
                  aic = double(ns*ns), 
                  bw = double(ns*ns),
                  b = double(ns*morder),
                  a = double(ns*morder), 
                  v = double(ns),
                  ks = integer(ns),
                  ke = integer(ns),
                  pxx = double(ns*121))

    aic <- list()
    bweight <- list()
    pacoef <- list()
    arcoef <- list()
    aic[[1]] <- NA
    bweight[[1]] <- NA
    for(i in 2:ns) {
        j <- ((i-1)*ns+1):((i-1)*ns+i)
        aic[[i]] <- z$aic[j]
        bweight[[i]] <- z$bw[j]
    }
    for(i in 1:ns)
        pacoef[[i]] <- z$b[((i-1)*morder+1):(i*morder)]
    for(i in 1:ns)
        arcoef[[i]] <- z$a[((i-1)*morder+1):(i*morder)]
    pspec <- array(z$pxx, dim = c(121, ns))

    if (plot == TRUE) {
        x <- rep(0,121)
        for(i in 1:121)
            x[i] <- (i - 1) / 240
        ymin <- min(pspec)
        ymax <- max(pspec)
        par(mfrow = c(ns, 1))
        nk <- ns
        if (ns == 4) {
            par(mfrow = c(2, 2))
            nk <- 4
        }
        if (ns > 4) {
            par(mfrow = c(3, 2))
            nk <- 6
        }
        if (ns > 6) {
            par(mfrow = c(3, 3))
            nk <- 9
        }
        nn <- 0
        for(i in 1:ns) {
            if (nn == nk) {
                par(ask = TRUE)
                nn <- 0
            }
            plot(x, pspec[,i], type = "l", ylim = c(ymin,ymax),
                 main = paste("y(", z$ks[i], "),...,y(", z$ke[i], ")"),
                 xlab = "Frequency", ylab = "Power Spectrum")
            nn <- nn + 1
        }
        par(mfrow = c(1, 1))
    }

    blocar.out <- list(mean = z$mean, var = z$var, aic = aic,
                       bweight = bweight, pacoef = pacoef, arcoef = arcoef,
                       v = z$v, init = z$ks, end = z$ke, pspec = pspec)
    class(blocar.out) <- "blocar"
    return(blocar.out)
}


print.blocar <- function(x, ...)
{
    cat(sprintf("\n Mean\t%f\n", x$mean))
    cat(sprintf(" Variance\t%f\n", x$var))

    ns <- length(x$bweight)
    for(i in 1:ns) {
        if (i != 1) {
            cat("\n\nAR-Model fitted to\t! Bayesian weights\t")
            cat("! AIC with respect to the present data\n")
            cat("----------------------------------------")
            cat("----------------------------------------\n")
            cat(sprintf("Current block\t\t! %f\t\t! %f\n", x$bweight[[i]][1],
                x$aic[[i]][1]))
            for(k in 1:(i-1))
                cat(sprintf("%i period former block\t! %f\t\t! %f\n", k,
                    x$bweight[[i]][k+1], x$aic[[i]][k+1]))
        }

        morder <- length(x$arcoef[[i]])
        cat("\n\n..........  Current model  (Average by the Bayesian weight)")
        cat("   ..........\n\n")
        cat(sprintf(" This model was fitted to the data  y( %i ),...,y( %i )\n",
            x$init[i], x$end[i]))
        cat(sprintf(" Innovation variance = %e\n\n", x$v[i]))
        cat(sprintf(" AR coefficients ( order %i ) \n", morder ))
        print(x$arcoef[[i]])
    }
    cat("\n")
}


blomar <- function (y, max.order = NULL, span)
{
    n <- nrow(y)
    d <- ncol(y)
# upper limit of the order of AR-model
    if (is.null(max.order))
        max.order <- as.integer(2*sqrt(n))

    icflag <- 0
    if (max.order > n/(2*d))
        icflag <- -1
    max.order <- as.integer(min(max.order, n/(2*d)))
    morder <- max.order

    calb <- rep(1, d)   # calibration constant for channel j (j=1,d)
    if (span < 1)
        span <- n
    ns <- as.integer((n-morder+span-1)/span)

    z <- .Fortran(C_blomarf,
                  as.double(y),
                  as.integer(n),
                  as.integer(d),
                  as.double(calb),
                  as.integer(morder),
                  as.integer(span),
                  as.integer(ns),
                  mean = double(d),
                  var = double(d),
                  bw = double(ns*ns),
                  raic = double(ns*ns),
                  a = double(d*d*morder*ns),
                  e = double(d*d*ns),
                  aic = double(ns),
                  ks = integer(ns),
                  ke = integer(ns),
                  nns = integer(1))

    nns <- z$nns
    bw <- array(z$bw, dim = c(ns, ns))
    bweight <- list()
    bweight[[1]] <- NA
    if (nns > 1)
        for(i in 2:nns)
            bweight[[i]] <- bw[1:i, i]

    raic <- array(z$raic, dim = c(ns, ns))
    aic <- list()
    aic[[1]] <- NA
    if (nns > 1)
        for(i in 2:nns)
            aic[[i]] <- raic[1:i, i]

    a <- array(z$a, dim = c(d, d, morder, ns))
    arcoef <- list()
    for(i in 1:nns)
        arcoef[[i]] <- array(a[, , , i], dim = c(d, d, morder))

    e <- array(z$e, dim = c(d, d, ns))
    v <- list()
    for(i in 1:nns)
        v[[i]] <- array(e[, , i], dim = c(d, d))

    aicbay <- z$aic[1:nns]
    start <- z$ks[1:nns]
    end <- z$ke[1:nns]

    blomar.out <- list(mean = z$mean, var = z$var, bweight = bweight,
                       aic = aic, arcoef = arcoef, v = v, eaic = aicbay,
                       init = start, end = end)

    class(blomar.out) <- "blomar"

    if (icflag == -1)
		warning(gettextf("max.order is corrected n/(2*d) = %d\n\n", max.order),
                domain=NA)

    return(blomar.out)
}

print.blomar <- function(x, ...)
{
    id <- length(x$mean)
    cat("\n\n Mean ")
    for(i in 1:id)
        cat(sprintf("  %f", x$mean[i]))
    cat("\n Variance ")
    for(i in 1:id)
        cat(sprintf("  %f", x$var[i]))
    cat("\n")

    ns <- length(x$bweight)
    for(i in 1:ns) {
        if (i != 1) {
            cat("\n\nAR-model fitted to\t! Bayesian weights\t")
            cat("! AIC with respect to the present data\n")
            cat("----------------------------------------")
            cat("----------------------------------------\n")
            cat(sprintf("Current block\t\t! %f\t\t! %f\n", x$bweight[[i]][1],
                x$aic[[i]][1]))
            for(k in 1:(i-1))
                cat(sprintf("%i Period former block\t! %f\t\t! %f\n", k,
                    x$bweight[[i]][k+1], x$aic[[i]][k+1]))
        }

        mf <- dim(x$arcoef[[i]])[3]
        cat("\n\n..........  Current model  ..........\n\n")
        cat(sprintf(" This model was fitted to the data  y( %i ),...,y( %i )\n",
            x$init[i], x$end[i]))
        cat(sprintf(" aic = %f\n", x$eaic[i]))
        cat("\n Innovation variance matrix\n")
        print(x$v[[i]])
        cat(sprintf("\n AR coefficient matrix ( order %i )\n", mf))
        print(x$arcoef[[i]])
    }
}


bsubst <- function (y, mtype, lag = NULL, nreg = NULL, reg = NULL,
                    term.lag = NULL, cstep = 5, plot = TRUE)
{
    if (mtype<0 || mtype>4) 
        stop("allowed model type are
             1 : autoregressive model,
             2 : polynomial type non-linear model (lag's read in),
             3 : polynomial type non-linear model (lag's automatically set),
             4 : AR-model with polynomial mean value function\n")

# specification of number of regressors (mtype=1,2,4)
    if (is.null(nreg))
        if (mtype != 3)
            stop("number of regressors is not specified")

# specification of maximum time lag  (mtype=2)
    if (is.null(reg))
        if (mtype == 2)
            stop("regressors are not specified")

# specification of regressor (mtype=3)
    if (is.null(term.lag))
    if (mtype == 3)
        stop("maximum time lags are not specified")

    n <- length(y)
    if (is.null(lag))
        lag <- as.integer(2*sqrt(n))   # maximum time lag use in the model
# number of regressors
    if (mtype == 3) {
        k <- 0
        lag5 <- term.lag[5]
        for(i in 1:lag5)
            k <- k + i * (i + 1) / 2 - 1
        nreg <- term.lag[1] + term.lag[2] + term.lag[3] + term.lag[4] + k
    }

    if (is.null(reg))
        reg <- array(0, dim = c(3, nreg)) 
    if (is.null(term.lag))
        term.lag <- rep(0, 5)
#    f <- ""                           # specification of regressor (mtype=5)
#    cnst <- 0                         # constant value (mtype=6)

    z <- .Fortran(C_bsubstf,
                  as.double(y),
                  as.integer(n),
                  as.integer(mtype),
                  as.integer(lag),
                  as.integer(nreg),
                  as.integer(cstep),
                  as.integer(reg),
                  as.integer(term.lag),
                  ymean = double(1),
                  yvar = double(1),
                  m = integer(1),
                  aicm = double(1),
                  vm = double(1),
                  a1 = double(nreg),
                  v = double(nreg+1),
                  aic = double(nreg+1),
                  daic = double(nreg+1),
                  aicb = double(1),
                  vb = double(1),
                  ek = double(1),
                  a2 = double(nreg),
                  ind = integer(nreg),
                  c = double(nreg),
                  c1 = double(nreg+1),
                  c2 = double(nreg),
                  b = double(nreg),
                  eicmin = double(1),
                  esum = double(nreg+1),
                  npm = double(1),
                  npmnreg = double(1),
                  e = double(n*cstep),
                  mean = double(cstep),
                  var = double(cstep),
                  skew = double(cstep),
                  peak = double(cstep),
                  cov = double(101),
                  pxx = double(121))

    mo <- z$m
    arcoefm=z$a1[1:mo]
    if (mtype != 4) {
        nps <- lag + 1
        perr <- array(z$e, dim = c(n, cstep))
        perr <- perr[nps:n, ]
        lagh <- 100
        nn <- n - nps - 1
        if ((lagh > nn) || (lagh == nn))
            lagh <- nn - 1
        lag1 <- lagh + 1
        autcor1=z$cov[1:lag1]
    }

    if (plot == TRUE) {

        if (mtype != 4)  {
            mm <- 2
            if (mtype == 1)
                mm <- mm + 1
            nc <- (mm + cstep + 1) / 2
            par(mfcol = c(nc, 2))
            for(i in 1:cstep) {
                sig <- sqrt(z$var[i]) * 0.5
                hist(perr[,i], breaks = c(-sig*10:1, 0, sig*1:10),
                    main = "", xlab = paste(i, "-step ahead prediction error"))
            }
        }

        if (mtype == 1) {
            k <- lag
            k1 <- k + 1
            plot((0:k), z$daic[1:k1], ylim = c(0,40), type = "l", xlab = "Lag",
                 ylab = "AIC(m)-aicmin")
            abline(h=0, lty=1)
        } else {
            plot((0:nreg), z$daic, ylim = c(0,40), type = "l", xlab = "Lag",
                 ylab = "AIC(m)-aicmin")
            abline(h=0, lty=1)
        }

        if (mtype == 1) {
            x <- rep(0, 121)
            for(i in 1:121)
                x[i] <- (i - 1) / 240
            plot(x, z$pxx, type = "l", xlab = "Frequency",
                 ylab = paste("Power Spectrum"))
        }


        if (mtype != 4) {
            plot(order <- c(0:lagh), autcor1, ylim = c(-1, 1), bty = "l",
                 type = "h",
                 main = "Autocorrelation of\n1-step ahead prediction error", 
                 xlab = "Lines show +/-2sd\n ( sd = sqrt(1/n) )",
                 ylab = "peautcor")
            abline(h = 0, lty = 1)
            abline(h = 2*sqrt(1/(n-lag)), lty = 3)
            abline(h = -2*sqrt(1/(n-lag)), lty = 3)
        }
        par(mfrow = c(1, 1))
    }

    if (mtype == 1) {
        k <- lag
        k1 <- k + 1
	bsubst.out <- list(ymean = z$ymean, yvar = z$yvar, v = z$v[1:k1],
                    aic = z$aic[1:k1], aicmin = z$aicm, daic = z$daic[1:k1],
                    order.maice = mo, v.maice = z$vm, arcoef.maice = arcoefm,
                    v.bay = z$vb, aic.bay = z$aicb, np.bay = z$ek,
                    arcoef.bay = z$a2[1:k], ind.c = z$ind[1:k],
                    parcor2 = z$c[1:k], damp = z$c1[1:k1],
                    bweight = z$c2[1:k], parcor.bay = z$b[1:k],
                    eicmin = z$eicmin, esum = z$esum[1:k1], npmean = z$npm,
                    npmean.nreg = z$npmnreg, perr = perr, mean = z$mean,
                    var = z$var, skew = z$skew, peak = z$peak,
                    peautcor = autcor1, pspec = z$pxx)
    }

    if (mtype==2 || mtype==3)
	bsubst.out <- list(ymean = z$ymean, yvar = z$yvar, v = z$v,
                    aic = z$aic, aicmin = z$aicm, daic = z$daic,
                    order.maice = mo, v.maice = z$vm, arcoef.maice = arcoefm,
                    v.bay = z$vb, aic.bay = z$aicb, np.bay = z$ek,
                    arcoef.bay = z$a2, ind.c = z$ind, parcor2 = z$c,
                    damp = z$c1, bweight = z$c2, parcor.bay = z$b,
                    eicmin = z$eicmin, esum = z$esum, npmean = z$npm,
                    npmean.nreg = z$npmnreg, perr = perr, mean = z$mean,
                    var = z$var, skew = z$skew, peak = z$peak,
                    peautcor = autcor1)

    if (mtype == 4)
	bsubst.out <- list(ymean = z$ymean, yvar = z$yvar, v = z$v,
                    aic = z$aic, aicmin = z$aicm, daic = z$daic,
                    order.maice = mo, v.maice = z$vm, arcoef.maice = arcoefm,
                    v.bay = z$vb, aic.bay = z$aicb, np.bay = z$ek,
                    arcoef.bay = z$a2, ind.c = z$ind, parcor2 = z$c,
                    damp = z$c1, bweight = z$c2, parcor.bay = z$b,
                    eicmin = z$eicmin, esum = z$esum, npmean = z$npm,
                    npmean.nreg = z$npmnreg)

    return(bsubst.out)
}

exsar <- function (y, max.order = NULL, plot = FALSE)
{
    n <- length(y)
    if (is.null(max.order))
        max.order <- as.integer(2*sqrt(n))
    morder <- max.order
	morder1 <- morder + 1

    z <- .Fortran(C_exsarf,
                  as.double(y),
                  as.integer(n),
                  as.integer(morder),
                  mean = double(1),
                  var = double(1),
                  v = double(morder1),
                  aic = double(morder1),
                  daic = double(morder1),
                  m = integer(1),
                  aicm = double(1),
                  sdm1 = double(1),
                  a1 = double(morder),
                  sdm2 = double(1),
                  a2 = double(morder),
                  ier = integer(1))

    if (z$ier != 0)
        stop("in FUNCT : SD is less than or equal to 0")

    if (plot == TRUE) {
        plot((0:morder), z$daic, ylim = c(0,40), type = "l", xlab = "Lag",
              ylab = "AIC(m)-aicmin (Truncated at 40.0)")
        abline(h = 0, lty = 1)
    }

    m <- z$m
    exsar.out  <- list(mean = z$mean, var = z$var, v = z$v, aic = z$aic,
                       aicmin = z$aicm, daic = z$daic, order.maice = m,
                       v.maice = z$sdm1, arcoef.maice = z$a1[1:m],
                       v.mle = z$sdm2, arcoef.mle = z$a2[1:m])

    return(exsar.out)
}


mlocar <- function (y, max.order = NULL, span, const = 0, plot = TRUE)
{
    n <- length(y)
    if (is.null(max.order))
        max.order <- as.integer(2*sqrt(n))
    morder <- max.order

    if (span < 1)
        span <- n
    ns <- as.integer((n-morder+span-1)/span)
    k <- morder + const

    z <- .Fortran(C_mlocarf,
                  as.double(y),
                  as.integer(n),
                  as.integer(morder),
                  as.integer(span),
                  as.integer(const),
                  as.integer(ns),
                  mean = double(1),
                  var = double(1),
                  a = double(k*ns),
                  mf = integer(ns),
                  sdf = double(ns),
                  ks = integer(ns),
                  ke = integer(ns),
                  pxx = double(121*ns),
                  ld1 = integer(ns),
                  ld2 = integer(ns),
                  ms = integer(ns),
                  sdms = double(ns),
                  aics = double(ns),
                  mp = integer(ns),
                  sdmp = double(ns),
                  aicp = double(ns))

    a <- array(z$a, dim = c(k, ns))
    mf <- z$mf
    mj.max <- 0
    for(i in 1:ns)
        if (mf[i] != 0) {
            mj.max <- i
        }
    ns <- mj.max
    arcoef <- list()
    for(i in 1:ns) {
        arcoef[[i]] <- a[1:mf[i], i]
    }

    pspec <- array(z$pxx, dim=c(121,ns))
    pspec <- pspec[, 1:ns]
    npre <- z$ld1
    order.const <- z$mp
    v.const <- z$sdmp
    aic.const <- z$aicp
    npre[1] <- NA
    order.const[1] <- NA
    v.const[1] <- NA
    aic.const[1] <- NA

    if (plot == TRUE) {
        oldpar <- par(no.readonly = TRUE)
        x <- rep(0, 121)
        for(i in 1:121)
            x[i] <- (i - 1) / 240
        ymin <- min(pspec)
        ymax <- max(pspec)
        par(mfrow=c(ns, 1))
        nk <- ns
        if (ns == 4) {
            par(mfrow = c(2, 2))
            nk <- 4
        }
        if (ns > 4) {
            par(mfrow = c(3, 2))
            nk <- 6
        }
        if (ns > 6) {
            par(mfrow = c(3, 3))
            nk <- 9
        }
        nn <- 0
        for(i in 1:ns) {
            if (nn == nk) {
                par(ask = TRUE)
                nn <- 0
            }
            plot(x, pspec[, i], type = "l", ylim = c(ymin, ymax),
                 main = paste("y(", z$ks[i], "),...,y(", z$ke[i], ")"),
                 xlab = "Frequency", ylab = "Power Spectrum")
            nn <- nn + 1
        }
        par(oldpar)
    }

    mlocar.out <- list(mean = z$mean, var = z$var, ns = ns, order = mf,
                       arcoef = arcoef, v = z$sdf, init = z$ks,
                       end = z$ke, pspec = pspec, npre = npre,
                       nnew = z$ld2, order.mov = z$ms, v.mov = z$sdms,
                       aic.mov = z$aics, order.const = order.const,
                       v.const = v.const, aic.const = aic.const)

    class(mlocar.out) <- "mlocar"
    return(mlocar.out)
}

print.mlocar <- function(x, ...)
{
    cat(sprintf("\n\n Mean\t%f\n", x$mean))
    cat(sprintf(" Variance\t%f\n", x$var))

    ns <- x$ns
    for(i in 1:ns) {
        if (i == 1) {
            cat(sprintf("\n Initial local model: (nnew = %i)", x$nnew[i]))
            cat(sprintf("\tvariance = %e\taic = %f\n", x$v.mov[i],
                x$aic.mov[i]))
        }
        if (i != 1) {
            cat("\n\n >>>  The following two models are compared  <<<\n")
            np <- x$npre[i] + x$nnew[i]
            cat(sprintf(" Moving model: (npre = %i, nnew = %i)\t", x$npre[i],
                x$nnew[i]))
            cat(sprintf("variance = %e\taic = %f\n", x$v.mov[i], x$aic.mov[i]))
            cat(sprintf(" Constant model: (npre+nnew = %i)\tvariance = %e",
                np, x$v.const[i]))
            cat(sprintf("\taic = %f\n", x$aic.const[i]))
            if (x$aic.mov[i] < x$aic.const[i]) {
                cat(" --->  New model adopted\n")
            } else {
                cat(" --->  Constant model adopted\n")
            }
        }
        cat("\n\n..........  Current model  ..........\n\n")
        cat(sprintf(" This model was fitted to the data  y( %i ),...,y( %i )\n",
            x$init[i], x$end[i]))
        cat(sprintf(" Innovation variance = %e\n\n", x$v[i]))
        cat(sprintf(" AR coefficients ( order %i ) \n", x$order[i]))
        print(x$arcoef[[i]])
    }
}


mlomar <- function (y, max.order = NULL, span, const = 0)
{
    n <- nrow(y)
    d <- ncol(y)
    if (is.null(max.order))
        max.order <- as.integer(2*sqrt(n))
    icflag <- 0
    if (max.order > n/(2*d))
        icflag <- -1
    max.order <- as.integer(min(max.order, n/(2*d)))
    morder <- max.order

    calb <- rep(1, d)   # calibration for channel j (j=1,d)
    if (span < 1)
        span <- n
    ns <- as.integer((n-morder+span-1)/span)

    z <- .Fortran(C_mlomarf,
                  as.double(y),
                  as.integer(n),
                  as.integer(d),
                  as.double(calb),
                  as.integer(morder),
                  as.integer(span),
                  as.integer(const),
                  as.integer(ns),
                  mean = double(d),
                  var = double(d),
                  ld1 = integer(ns),
                  ld2 = integer(ns),
                  ms = integer(ns),
                  aicm = double(ns),
                  mp = integer(ns),
                  aicc = double(ns),
                  mf = integer(ns),
                  aic = double(ns),
                  a = double(d*d*morder*ns),
                  e = double(d*d*ns),
                  ks = integer(ns),
                  ke = integer(ns),
                  nns = integer(1))

    nns <- z$nns
    npre <- z$ld1[1:nns]
    nnew <- z$ld2[1:nns]
    order.mov <- z$ms[1:nns]
    aic.mov <- z$aicm[1:nns]
    order.const <- z$mp[1:nns]
    aic.const <- z$aicc[1:nns]
    order <- z$mf[1:nns]
    aic <- z$aic[1:nns]
    start <- z$ks[1:nns]
    end <- z$ke[1:nns]
    npre[1] <- NA
    order.mov[1] <- NA
    aic.mov[1] <- NA
    order.const[1] <- NA
    aic.const[1] <- NA

    a <- array(z$a, dim = c(d, d, morder, ns))
    arcoef <- list()
    for(i in 1:nns)
        arcoef[[i]] <- array(a[, , , i], dim = c(d, d, z$mf[i]))

    e <- array(z$e, dim = c(d, d, ns))
    v <- list()
    for(i in 1:nns)
        v[[i]] <- array(e[, , i], dim = c(d, d))

    mlomar.out <- list(mean = z$mean, var = z$var, ns = nns, order = order,
                       aic = aic, arcoef = arcoef, v = v, init = start,
                       end = end, npre = npre, nnew = nnew,
                       order.mov = order.mov, aic.mov = aic.mov,
                       order.const = order.const, aic.const = aic.const)

    class(mlomar.out) <- "mlomar"

    if (icflag == -1)
		warning(gettextf("max.order is corrected n/(2*d) = %d\n\n", max.order),
                domain=NA)

    return(mlomar.out)
}

print.mlomar <- function(x, ...)
{
    id <- length(x$mean)
    cat("\n\n Mean ")
    for(i in 1:id)
        cat(sprintf("  %f", x$mean[i]))
    cat("\n Variance ")
    for(i in 1:id)
        cat(sprintf("  %f", x$var[i]))

    ns <- x$ns
    for(i in 1:ns) {
        if (i == 1)
            cat(sprintf("\n\n Initial local model: (nnew = %i), aic = %f\n",
                x$nnew[i], x$aic[i]))

        if (i != 1) {
            cat("\n\n >>>  The following two models are compared  <<<\n")
            np <- x$npre[i] + x$nnew[i]
            cat(sprintf(" Moving model: (npre = %i, nnew = %i)\taic = %f\n",
                x$npre[i], x$nnew[i], x$aic.mov[i]))
            cat(sprintf(" Constant model: (npre+nnew = %i)\taic = %f\n",
                np,x$aic.const[i]))
            if (x$aic.mov[i] < x$aic.const[i]) {
                cat(" ---> New model adopted \n")
            } else {
                cat(" ---> Constant model adopted \n")
            }
        }
        cat("\n\n..........  Current model  ..........\n\n")
        cat(sprintf(" This model was fitted to the data  y( %i ),...,y( %i )\n",
            x$init[i], x$end[i]))
        cat(sprintf(" aic = %f\n", x$aic[i]))
        cat("\n Innovation variance matrix\n")
        print(x$v[[i]])
        cat(sprintf("\n AR coefficient matrix ( order %i )\n", x$order[i]))
        print(x$arcoef[[i]])
    }
}


mulbar <- function (y, max.order = NULL, plot = FALSE)
{
    n <- nrow(y)
    d <- ncol(y)
    calb <- rep(1, d)   # calibration of channel i (i=1,d)
    if (is.null(max.order))   # upper limit of the order of AR-model
        max.order <- as.integer(2*sqrt(n))
    icflag <- 0
    if (max.order > n/(2*d))
        icflag <- -1
    max.order <- as.integer(min(max.order, n/(2*d)))
    morder <- max.order

	morder1 <- morder + 1
	d2 <- d * d

    z <- .Fortran(C_mulbarf,
                  as.double(y),
                  as.integer(n),
                  as.integer(d),
                  as.double(calb),
                  as.integer(morder),
                  mean = double(d),
                  var = double(d),
                  v = double(morder1),
                  aic = double(morder1),
                  daic = double(morder1),
                  m = integer(1),
                  aicm = double(1),
                  vm = double(1),
                  w1 = double(morder1),
                  w2 = double(morder),
                  a = double(d2*morder),
                  b = double(d2*morder),
                  g = double(d2*morder),
                  h = double(d2*morder),
                  e = double(d2),
                  aicb = double(1))

    if (plot == TRUE) {
        plot((0:morder), z$daic, ylim = c(0,40), type = "l", xlab = "Lag",
             ylab = "AIC - aicmin (Truncated at 40.0)")
        abline(h = 0, lty = 1)
    }

    mulbar.out <- list(mean = z$mean, var = z$var, v = z$v,
                       aic = z$aic, aicmin = z$aicm, daic = z$daic,
                       order.maice = z$m, v.maice = z$vm,
                       bweight = z$w1, integra.bweight = z$w2,
                       arcoef.for = array(z$a, dim = c(d, d, morder)),
                       arcoef.back = array(z$b, dim = c(d, d, morder)),
                       pacoef.for = array(z$g, dim = c(d, d, morder)),
                       pacoef.back = array(z$h, dim = c(d, d, morder)),
                       v.bay = array(z$e, dim = c(d, d)),
                       aic.bay = z$aicb)

    if (icflag == -1)
		warning(gettextf("max.order is corrected n/(2*d) = %d\n\n", max.order),
                domain=NA)

    return(mulbar.out)
}


mulmar <- function (y, max.order = NULL, plot = FALSE) 
{
    n <- nrow(y)
    d <- ncol(y)
    calb <- rep(1, d)
    if (is.null(max.order))
        max.order <- as.integer(2*sqrt(n))
    icflag <- 0
    if (max.order > n/(2*d))
        icflag <- -1
    max.order <- as.integer(min(max.order, n/(2*d)))
	lag <- max.order
    lag1 <- lag + 1

    z <- .Fortran(C_mulmarf,
                  as.double(y),
                  as.integer(n),
                  as.integer(d),
                  as.double(calb),
                  as.integer(max.order),
                  mean = double(d),
                  var = double(d),
                  v = double(lag1*d),
                  aic = double(lag1*d),
                  daic = double(lag1*d),
                  m = integer(d),
                  aicm = double(d),
                  vm = double(d),
                  npr = integer(d),
                  jnd = integer(lag1*d*d),
                  a = double(lag1*d*d),
                  rv = double(d),
                  aicf = double(d),
                  ei = double(d*d),
                  bi = double(d*d*lag),
                  matv = double(d*d),
                  arcoef = double(d*d*lag),
                  morder = integer(1),
                  aics = double(1))

    v <- list()
    aic <- list()
    daic <- list()
    for(i in 1:d) {
        j <- ((i - 1) * lag1 + 1):(i * lag1)
        v[[i]] <- z$v[j]
        aic[[i]] <- z$aic[j]
        daic[[i]] <- z$daic[j]
    }
    jnd <- list()
    subregcoef <- list()
    ind <- array(z$jnd, dim = c(lag1*d, d))
    a <- array(z$a, dim = c(lag1*d, d))
    for(i in 1:d)
        jnd[[i]] <- ind[(1:z$npr[[i]]), i]
    for(i in 1:d)
        subregcoef[[i]] <- a[(1:z$npr[[i]]), i]

    if (plot == TRUE) {
        par(mfrow = c(d, 1))  
        for(i in 1:d) {
            plot((0:max.order), daic[[i]], ylim = c(0, 40), type = "l", 
                 main = paste(" d=", i), xlab = "Lag",
                 ylab = "AIC(m)-aicmin (Truncated at 40.0)")
            abline(h = 0, lty = 1)
        }
        par(mfrow = c(1, 1))
    }

    morder <- z$morder
    mulmar.out <- list(mean = z$mean, var = z$var, v = v, aic = aic,
                       aicmin = z$aicm, daic = daic, order.maice = z$m,
                       v.maice = z$vm, np = z$npr, jnd = jnd,
                       subregcoef = subregcoef, rvar = z$rv, aicf = z$aicf,
                       respns = array(z$ei, dim = c(d, d)),
                       regcoef = array(z$bi, dim = c(d, d, morder)),
                       matv = array(z$matv, dim = c(d, d)), morder = morder,
                       arcoef = array(z$arcoef, dim = c(d, d, morder)),
                       aicsum = z$aics)

    if (icflag == -1)
		warning(gettextf("max.order is corrected n/(2*d) = %d\n\n", max.order),
                domain=NA)

    return(mulmar.out)
}


perars <- function (y, ni, lag = NULL, ksw = 0)
{
    n <- length(y)
    if (is.null(lag))
        lag <- as.integer(2*sqrt(ni))
    lag1 <- lag + 1
    k <- (lag1*ni + ksw)*ni
	
    z <- .Fortran(C_perarsf,
                  as.double(y),
                  as.integer(n),
                  as.integer(ni),
                  as.integer(lag),
                  as.integer(ksw),
                  mean = double(1),
                  var = double(1),
                  np = integer(ni),
                  jnd = integer(k),
                  a = double(k),
                  aic = double(ni),
                  b = double(ni*ni*lag),
                  v = double(ni*ni),
                  c = double(ni),
                  osd = double(ni),
                  mo = integer(1))

    npara <- z$np
    ind <- array(z$jnd, dim = c(lag1*ni+ksw, ni))
    jnd <- list()
    for(i in 1:ni)
        jnd[[i]] <- ind[1:npara[i], i]

    a <- array(z$a, dim = c(lag1*ni+ksw, ni))
    regcoef <- list()
    for(i in 1:ni)
        regcoef[[i]] <- a[1:npara[i], i]

    perars.out <- list(mean = z$mean, var = z$var, subset = jnd,
                       regcoef = regcoef, rvar = z$osd, np = npara,
                       aic = z$aic, v = array(z$v, dim = c(ni, ni)),
                       arcoef = array(z$b, dim = c(ni, ni, z$morder)),
                       const = z$c, morder = z$mo)

    class(perars.out) <- "perars"
    return(perars.out)

}

print.perars <- function(x, ...)
{
    cat(sprintf("\n\n Mean = %f\n", x$mean))
    cat(sprintf(" Variance = %f\n", x$var))

    ni <- nrow(x$v)
    for(i in 1:ni) {
        cat("\n ........  Regression Model for the regressand")
        cat(sprintf(" i = %i  ........\n", i))
        cat("\n subset")
        sorder <- length(x$subset[[i]])
        for(j in 1:sorder)
            cat(sprintf("\t\t%i", x$subset[[i]][j]))
        cat("\n regcoef")
        for(j in 1:sorder)
            cat(sprintf("\t%f", x$regcoef[[i]][j]))
        cat(sprintf("\n\n residual variance (rvar) = %f\n", x$rvar[i]))
        cat(sprintf(" number of parameter (np) = %i\n", x$np[i]))
        cat(sprintf(" AIC = n*log(rvar) + 2*np = %f\n\n", x$aic[i]))
    }

    cat("\n\nmatrix of regression coefficients\n")
    print(x$v)
    cat("\n\nregression coefficients within the present period\n")
    print(x$arcoef)
    cat("\n\nconstants within the regression models\n")
    print(x$const)
}


unibar <- function (y, ar.order = NULL, plot = TRUE)
{
    n <- length(y)
    if (is.null(ar.order))
        ar.order <- as.integer(2*sqrt(n))
    ar.order1 <- ar.order + 1

    z <- .Fortran(C_unibarf,
                  as.double(y),
                  as.integer(n),
                  as.integer(ar.order),
                  mean = double(1),
                  var = double(1),
                  v = double(ar.order1),
                  aic = double(ar.order1),
                  daic = double(ar.order1),
                  m = integer(1),
                  aicm = double(1),
                  vm = double(1),
                  pa = double(ar.order),
                  bw = double(ar.order1),
                  sbw = double(ar.order),
                  pab = double(ar.order),
                  aicb = double(1),
                  vb = double(1),
                  pn = double(1),
                  a = double(ar.order),
                  pxx = double(121))


    if (plot == TRUE) {
        oldpar <- par(no.readonly = TRUE)
        par(mfrow = c(3, 1))
        plot((0:ar.order), z$daic, ylim = c(0,40), type = "l", xlab = "Lag",
             ylab = "AIC(m)-aicmin (Truncated at 40.0)")
        abline(h = 0, lty = 1)
        plot(z$pa, type = "h", xlab = "Lag", ylab = "Partial autocorrelation")
        abline(h = 0, lty = 1)
#        abline(h = 1/sqrt(n), lty = 3)
#        abline(h = -1/sqrt(n),lty = 3)
        abline(h = 2/sqrt(n), lty = 3)
        abline(h = -2/sqrt(n), lty = 3)
        x <- rep(0, 121)
        for(i in 1:121)
            x[i] <- (i - 1) / 240
        plot(x, z$pxx, type = "l", xlab = "Frequency",
             ylab = "Power Spectral Density")
        par(oldpar)
        unibar.out <- list(mean = z$mean, var = z$var, v = z$v,
                           aic = z$aic, aicmin = z$aicm, order.maice = z$m,
                           v.maice = z$vm, bweight = z$bw[2:(ar.order1)],
                           integra.bweight = z$sbw, v.bay = z$vb,
                           aic.bay = z$aicb, np = z$pn,
                           pacoef.bay = z$pab, arcoef = z$a)
    } else {
        unibar.out <- list(mean = z$mean, var = z$var, v = z$v,
                           aic = z$aic, aicmin = z$aicm, daic = z$daic,
                           order.maice = z$m, v.maice = z$vm,
                           pacoef = z$pa, bweight = z$bw[2:(ar.order1)],
                           integra.bweight = z$sbw, v.bay = z$vb,
                           aic.bay = z$aicb, np = z$pn,
                           pacoef.bay = z$pab, arcoef = z$a,
                           pspec = z$pxx)
    }
    return(unibar.out)
}


unimar <- function (y, max.order = NULL, plot = FALSE)
{
    n <- length(y)
    if (is.null(max.order))
        max.order <- as.integer(2*sqrt(n))  # upper limit of AR-order
    morder <- max.order
    morder1 <- morder + 1

    z <- .Fortran(C_unimarf,
                  as.double(y),
                  as.integer(n),
                  as.integer(morder),
                  mean = double(1),
                  var = double(1),
                  v = double(morder1),
                  aic = double(morder1),
                  daic = double(morder1),
                  m = integer(1),
                  aicm = double(1),
                  vm = double(1),
                  a = double(morder))

    m <- z$m
    if (plot == TRUE) {
        plot((0:morder), z$aicm, ylim = c(0,40), type = "l", xlab = "Lag",
             ylab = "AIC(m)-aicmin (Truncated at 40.0)")
      abline(h = 0, lty = 1)
      unimar.out <- list(mean = z$mean, var = z$var, v = z$v,
                         aic = z$aic, aicmin = z$aicm, order.maice = m,
                         v.maice = z$vm, arcoef = z$a[1:m])
    } else { 
        unimar.out <- list(mean = z$mean, var = z$var, v = z$v,
                           aic = z$aic, aicmin = z$aicm, daic = z$daic,
                           order.maice = m, v.maice = z$vm,
                           arcoef = z$a[1:m])
    }
    return(unimar.out)
}


xsarma <- function (y, arcoefi, macoefi)
{
    n <- length(y)
    arcoefi <- -arcoefi   # Initial estimates of AR coefficients
    macoefi <- -macoefi   # Initial estimates of MA coefficients
    p <- length(arcoefi)  # AR-ORDER
    q <- length(macoefi)  # MA-ORDER
    p01 <- c(arcoefi, macoefi)   # INITIAL ESTIMATES OF AR- AND  MA-COEFFICIENTS

    z <- .Fortran(C_xsarmaf,
                  as.double(y),
                  as.integer(n),
                  as.integer(p),
                  as.integer(q),
                  as.double(p01),
                  g1 = double(p+q),
                  tl1 = double(1),
                  p02 = double(p+q),
                  g2 = double(p+q),
                  alpar = double(p),
                  alpma = double(q),
                  tl2 = double(1),
                  sigma2 = double(1))

    coef <- -z$p02

    xsarma.out <- list(gradi = z$g1, lkhoodi = z$tl1, arcoef = coef[1:p],
                       macoef = coef[(p+1):(p+q)], grad = z$g2,
                       alph.ar = z$alpar, alph.ma = z$alpma, lkhood = z$tl2,
                       wnoise.var = z$sigma2)
    return(xsarma.out)
}

