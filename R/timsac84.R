#####   TIMSAC84   #####

baysea <- function(y, period = 12, span = 4, shift = 1, forecast = 0,
                   trend.order = 2, seasonal.order = 1, year = 0, month = 1,
                   out = 0, rigid = 1, zersum = 1, delta = 7, alpha = 0.01,
                   beta = 0.01, gamma = 0.1, spec = TRUE, plot = TRUE,
                   separate.graphics = FALSE)
{
    if (seasonal.order > span)
        stop("seasonal.order is smaller than or equal to span")
    if (span < 1)
        stop("span is greater than or equal to 1")
    if (trend.order < 1)
        stop("trend.order is greater than or equal to 1")

    ndata <- length(y)
    npf <- ndata + forecast

    ipara <- rep(0, 12)
    ipara[1] <- period
    ipara[2] <- span
    ipara[3] <- shift
    ipara[4] <- trend.order
    ipara[5] <- seasonal.order
    ipara[6] <- 0            # logt
    ipara[7] <- year
    ipara[8] <- month
    ipara[9] <- 1            # nday
    if (spec == TRUE)
        ipara[10] <- 1       # spectrum estimation option
    if (spec == FALSE)
        ipara[10] <- 0
    ipara[11] <- out         # ioutd : outlier correction option 

    para <- rep(0, 8)
    para[1] <- rigid     # controls the rigidity of the seasonal component
    para[2] <- 1         # wtrd
    para[3] <- 1         # dd
    para[4] <- zersum    # controls the sum of the seasonals within a period
    para[5] <- delta     # controls the leap year effect
    para[6] <- alpha     # controls prior variance of initial trend
    para[7] <- beta      # controls prior variance of initial seasonal
    para[8] <- gamma     # controls prior variance of initial sum of seasonal

    arft <- rep(0, 3)
    arfs <- rep(0, 3)
    arfn <- rep(0, 3)

#   subroutine arcoef ---> subroutine partar
    iart <- 0
#    for(i in 1:3) if(arft[i] != 0.0) iord <- i
#    if(iord != 0) {
#      arm <- matrix(0, dim = c(3,3))
#      for(i in 1:3) arm[i,i] <- arft[i]
#      arm[2,1] <- arm[1,1]-arft[2]*arm[1,1]
#      arm[3,1] <- arm[2,1]-arft[3]*arm[2,2]
#      arm[3,2] <- arm[2,2]-arft[3]*arm[2,1]
#      for(i in 1:iord) arft[i] <- arm[iord,i]
#    }
    iars <- 0
    iarn <- 0

    is <- period * seasonal.order
    iprd <- 2
    if (period == 1)
        iprd <- 1
    lftrn <- trend.order + iart
    lfsea <- (seasonal.order + iars) * period + iarn
    idc <- lftrn * iprd + 1
    idcx <- lfsea * 2 + 1
    if ((period > 1) && (idc < idcx))
        idc <- idcx
    if ((period > 1) && (idc < period*2-1))
        idc <- period * 2 - 1
    ipara[12] <- idc

    z <- .Fortran(C_bayseaf,
                  as.double(y),
                  as.integer(ndata),
                  as.integer(forecast),
                  outlier = double(ndata),
                  dmoi = double(ndata),
                  trend = double(npf),
                  season = double(npf),
                  tdcmp = double(npf),
                  irreg = double(ndata),
                  adjust = double(ndata),
                  est = double(npf),
                  psds = double(npf),
                  psdt = double(npf),
                  avabic = double(1),
                  as.integer(ipara),
                  as.double(para),
                  as.double(arft),
                  as.double(arfs),
                  as.double(arfn),
                  as.integer(iart),
                  as.integer(iars),
                  as.integer(iarn))	

    outlier <- NULL
    if (out != 0)
        outlier <- z$outlier
    tre <- z$trend
    sea <- z$season
    tday <- NULL
    if (year != 0)
        tday <- z$tdcmp
    irr <- z$irreg
    adj <- z$adjust
    est <- z$est
    psds <- z$psds
    psdt <- z$psdt
    abic <- z$avabic

    baysea.out <- list(outlier = outlier, trend = tre, season = sea, tday = tday,
                       irregular = irr, adjust = adj, smoothed = est,
                       aveABIC = abic)

    if (spec == TRUE) {
        z1 <- spec.baysea(y, period, tre, trend.order, sea, seasonal.order, irr,
                          adj)
        spec1 <- z1$irregular.spec
        spec2 <- z1$adjusted.spec
        spec3 <- z1$differenced.trend
        spec4 <- z1$differenced.season

        if (plot == TRUE)
            plot.baysea(y, period, outlier, tre, trend.order, sea,
                        seasonal.order, tday, irr, adj, est, psdt, psds, spec1,
                        spec2, spec3, spec4, spec, separate.graphics)

        ir.spec <- list(acov = spec1$acov, acor = spec1$acor, mean = spec1$mean,
                        v = spec1$v, aic = spec1$aic, parcor = spec1$parcor,
                        rspec = spec1$rspec) 
        ad.spec <- list(acov = spec2$acov, acor = spec2$acor, mean = spec2$mean,
                        v = spec2$v, aic = spec2$aic, parcor = spec2$parcor,
                        rspec = spec2$rspec) 
        diff.trend <- list(acov = spec3$acov, acor = spec3$acor,
                           mean = spec3$mean, v = spec3$v, aic = spec3$aic,
                           parcor = spec3$parcor)
        diff.season <- list(acov = spec4$acov, acor = spec4$acor,
                            mean = spec4$mean, v = spec4$v, aic = spec4$aic,
                            parcor = spec4$parcor)

        baysea.out <- c(baysea.out, irregular.spec = ir.spec,
                        adjusted.spec = ad.spec, differenced.trend = diff.trend,
                        differenced.season = diff.season)
    } else {
        if (plot == TRUE) {
            ir.spec <- NULL
            ad.spec <- NULL
            diff.trend <- NULL
            diff.season <- NULL
            plot.baysea(y, period, outlier, tre, trend.order, sea,
                        seasonal.order, tday, irr, adj, est, psdt, psds, ir.spec,
                        ad.spec, diff.trend, diff.season, spec,
                        separate.graphics)
        }
    }

    return(baysea.out)
}

spec.baysea <- function(y, period, trend, trend.order, season, seasonal.order,
                        irregular, adjust)
{
    ndata <- length(y)
    npf <- length(trend)

    lag <- min((ndata-1), 60)
    lag1 <- lag + 1
#    ifpl <- min(3*sqrt(ndata), 50, lag)
    ifpl <- min(30, ndata-1)
    ifpl1 <- ifpl + 1

 # SPECTRUM OF IRREGULAR                               
    mode <- 1
    z <- .Fortran(C_spgrh,
                  as.double(irregular),
                  as.integer(ndata),
                  as.integer(lag1),
                  as.integer(ifpl1),
                  as.integer(mode),
                  as.integer(period),
                  cxx = double(lag1),
                  cn = double(lag1),
                  xmean = double(1),
                  sd = double(ifpl1),
                  aic = double(ifpl1),
                  parcor = double(ifpl),
                  pxx = double(lag1),
                  ier = integer(1))
	
    pxx <- z$pxx
    sxx <- rep(0, lag1)
    for(i in 1:lag1) {
        t <- abs(pxx[i])
        sxx[i] <- log10(t)
    }

    if (z$ier == 2600)
        warning("accuracy of computation lost")

    irregular.spec <- list(n = ndata, acov = z$cxx, acor = z$cn,
                           mean = z$xmean, v = z$sd, aic = z$aic,
                           parcor = z$parcor, rspec = pxx, rpspec = sxx)

# SPECTRUM OF DIFFERENCED ADJUSTED SERIES
    n1 <- ndata - 1
    dadj <- adjust
    for(i in 1:n1)
        dadj[i] <- dadj[i+1] - dadj[i]
    dadj <- dadj[1:n1]
    lag <- min((n1-1), 60)
    lag1 <- lag + 1

    z <- .Fortran(C_spgrh,
                  as.double(dadj),
                  as.integer(n1),
                  as.integer(lag1),
                  as.integer(ifpl1),
                  as.integer(mode),
                  as.integer(period),
                  cxx = double(lag1),
                  cn = double(lag1),
                  xmean = double(1),
                  sd = double(ifpl1),
                  aic = double(ifpl1),
                  parcor = double(ifpl),
                  pxx = double(lag1),
                  ier = integer(1))

    pxx <- z$pxx
    sxx <- rep(0, lag1)
    for(i in 1:lag1) {
        t <- abs(pxx[i])
        sxx[i] <- log10(t)
    }

    if (z$ier == 2600)
        cat(" ***** Warnning : Accuracy of computation lost\n")

    adjusted.spec <- list(n = n1, acov = z$cxx, acor = z$cn, mean = z$xmean,
                          v = z$sd, aic = z$aic, parcor = z$parcor,
                          rspec = pxx, rpspec = sxx)

#  PARCOR OF TREND.ORDER TIME(S) DIFFERENCED TREND SERIES
    n1 <- ndata		
    trendd <- trend
    for(j in 1:trend.order) {
        n1 <- n1 - 1
        for(i in 1:n1)
            trendd[i] <- trendd[i+1] - trendd[i]
    }
    trendd <- trendd[1:n1]
    lag <- min((n1-1), 60)
    lag1 <- lag + 1
    mode <- 0

    z <- .Fortran(C_spgrh,
                  as.double(trendd),
                  as.integer(n1),
                  as.integer(lag1),
                  as.integer(ifpl1),
                  as.integer(mode),
                  as.integer(period),
                  cxx = double(lag1),
                  cn = double(lag1),
                  xmean = double(1),
                  sd = double(ifpl1),
                  aic = double(ifpl1),
                  parcor = double(ifpl),
                  pxx = double(lag1),
                  ier = integer(1))

    if (z$ier == 2600)
        warning("accuracy of computation lost")

    differenced.trend <- list(n = n1, acov = z$cxx, acor = z$cn,
                              mean = z$xmean, v = z$sd, aic = z$aic,
                              parcor = z$parcor)

#  PARCOR OF SEASONAL.ORDER TIME(S) DIFFERENCED SEASONAL SERIES
    n1 <- ndata
    seasond <- season
    for(j in 1:seasonal.order) {
        n1 <- n1 - period
        for(i in 1:n1)
            seasond[i] <- seasond[i+period] - seasond[i]
    }
    seasond <- seasond[1:n1]

    lag <- min((n1-1), 60)
    lag1 <- lag + 1

    z <- .Fortran(C_spgrh,
                  as.double(seasond),
                  as.integer(n1),
                  as.integer(lag1),
                  as.integer(ifpl1),
                  as.integer(mode),
                  as.integer(period),
                  cxx = double(lag1),
                  cn = double(lag1),
                  xmean = double(1),
                  sd = double(ifpl1),
                  aic = double(ifpl1),
                  parcor = double(ifpl),
                  pxx = double(lag1),
                  ier = integer(1))

    if (z$ier == 2600)
        warning("accuracy of computation lost")

    differenced.season <- list(n = n1, acov = z$cxx, acor = z$cn,
                               mean = z$xmean, v = z$sd, aic = z$aic,
                               parcor = z$parcor)

    spec.baysea.out <- list(irregular.spec = irregular.spec,
                            adjusted.spec = adjusted.spec,
                            differenced.trend = differenced.trend,
                            differenced.season = differenced.season)
}


plot.baysea <- function(y, period, outlier, trend, trend.order, season,
                        seasonal.order, tday, irregular, adjust, smoothed, psdt,
                        psds, irregular.spec, adjusted.spec, differenced.trend,
                        differenced.season, spec, separate.graphics)

{
    oldpar <- par(no.readonly = TRUE)
    ndata <- length(y)
    npf <- length(trend)

    ymax1 <- max(y, trend, adjust, smoothed)
    ymin1 <- min(y, trend, adjust, smoothed)

    plot(y, type = "l", main = "Original Data", xlab = "time", xlim = c(0,npf),
         ylim = c(ymin1,ymax1))
    nw <- 1

    if (separate.graphics == TRUE) {
        dev.new()
    } else {
        par(ask = TRUE)
    }

    plot(trend, type = "l", main = "Trend and 2*(post SD)", cex.main = 0.9,
         xlab = "time", ylab = "", xlim = c(0, npf), ylim = c(ymin1, ymax1))
    par(new = TRUE)
    xtem <- trend + psdt
    plot(xtem, type = "l", lty = 3, main = "", xlab = "", ylab = "",
         xlim = c(0, npf), ylim = c(ymin1, ymax1))
    par(new = TRUE)
    xtem <- trend - psdt
    plot(xtem, type = "l", lty = 3, main = "", xlab = "", ylab = "",
         xlim = c(0, npf), ylim = c(ymin1, ymax1))

    if (separate.graphics == TRUE)
        dev.new()
    plot(adjust, pch = 0, type = "l",
       main = "Adjusted = Original Data - Seasonal - Trading.Day.Comp - Outlier",
        cex.main = 0.9, xlab = "time", ylab = "", xlim = c(0, npf),
        ylim = c(ymin1, ymax1))

    if (separate.graphics == TRUE)
        dev.new()
    plot(smoothed, pch = 0, type = "l",
        main = "Smoothed = Trend + Seasonal + Trading.Day.Comp", cex.main = 0.9,
        xlab = "time", ylab = "", xlim = c(0, npf), ylim = c(ymin1, ymax1))

    ymax2 <- max(irregular)
    ymin2 <- min(irregular)
    if (seasonal.order != 0) {
        ymax2 <- max(season, ymax2)
        ymin2 <- min(season, ymin2)
    }
    if (is.null(tday) == FALSE) {
        ymax2 <- max(tday, ymax2)
        ymin2 <- min(tday, ymin2)
    }
    my <- max(ymax2, abs(ymin2)) * 1.5
    if (seasonal.order != 0) {
        if (separate.graphics == TRUE)
            dev.new()
        plot(season, type = "l", main = "Seasonal and 2*(post SD)",
             cex.main = 0.9,  xlab = "time", ylab = "", ylim = c(-my, my))
        par(new = TRUE)
        xtem <- season + psds
        plot(xtem, type = "l", lty = 3, main = "", xlab = "", ylab = "",
             ylim = c(-my, my))
        par(new = TRUE)
        xtem <- season - psds
        plot(xtem, type = "l", lty = 3, main = "", xlab = "", ylab = "",
             ylim = c(-my, my))
    }

    if (separate.graphics == TRUE)
        dev.new()
    par(mfrow = c(2, 1))
    plot(irregular, type = "l",
       main = "Irregular = Original Data - Trend - Seasonal - Trading.Day.Comp",
         cex.main = 0.9, xlab = "time", ylab = "", ylim = c(-my, my))
    vy <- sd(irregular) * 5
    plot(irregular, type = "l",
         main = "Irregular ( Scaled by the Standard Deviation)", cex.main = 0.9,
         xlab = "time", ylab = "", ylim = c(-vy, vy))
    par(mfrow = c(1, 1))

 # SPECTRUM OF IRREGULAR                               
    if (spec == TRUE) {
        lag <- length(irregular.spec$acov) - 1
        n <- irregular.spec$n
        if (separate.graphics == TRUE)
            dev.new()
        par(mfrow = c(3, 1))
        plot((0:lag), irregular.spec$acor, type = "h",
             main = "Autocorrelation & Parcor of Irregular",
             ylab = "Autocorrelation", xlab = "lag", ylim = c(-1, 1))
        plot(irregular.spec$parcor, type = "h", ylab = "Parcor", xlab = "lag",
             ylim = c(-1, 1))
        abline(h = 0, lty = 1)
        abline(h = 2/sqrt(n), lty = 3)
        abline(h = -2/sqrt(n), lty = 3)

        it <- irregular.spec$rpspec * 10
        ymin <- (min(irregular.spec$rpspec) - 1) * 10
        ymax <- (max(irregular.spec$rpspec) + 1) * 10
        plot(it, type = 'l',
             main = "High Order AR-Spectrum as an approximation
 to periodgram ( Order is fixed at 30 )",
             ylab = "Rational Spectrum", xlab = "", ylim = c(ymin, ymax))
        par(mfrow = c(1, 1))

# SPECTRUM OF DIFFERENCED ADJUSTED SERIES
        lag <- length(adjusted.spec$acov) - 1
        n <- adjusted.spec$n
        if (separate.graphics == TRUE)
            dev.new()
        par(mfrow = c(3, 1))
        plot((0:lag), adjusted.spec$acor, type = "h",
             main="Autocorrelation & Parcor of Differenced Adjusted Series",
             ylab = "Autocorrelation", xlab = "lag", ylim = c(-1, 1))
        plot(adjusted.spec$parcor, type = "h", ylab = "Parcor", xlab = "lag",
             ylim = c(-1, 1))
        abline(h = 0, lty = 1)
        abline(h = 2/sqrt(n), lty = 3)
        abline(h = -2/sqrt(n), lty = 3)

        it <- adjusted.spec$rpspec * 10
        ymin <- (min(adjusted.spec$rpspec) - 1) * 10
        ymax <- (max(adjusted.spec$rpspec) + 1) * 10
        plot(it, type = 'l',
             main = "High Order AR-Spectrum as an approximation
 to periodgram( Order is fixed at 30 )",
             xlab = "", ylab = "Rational Spectrum", ylim = c(ymin, ymax))
        par(mfrow = c(1, 1))

#  PARCOR OF TREND.ORDER TIME(S) DIFFERENCED TREND SERIES
        lag <- length(differenced.trend$acov) - 1
        n <- differenced.trend$n
        if (separate.graphics == TRUE)
            dev.new()
        par(mfrow = c(2, 1))
        plot((0:lag), differenced.trend$acor, type = "h",
             main=paste(trend.order, "Time(s) Differenced Trend Series"),
             ylab = "Autocorrelation", cex.main = 0.9, xlab = "lag",
             ylim = c(-1, 1))
        plot(differenced.trend$parcor, type = "h", ylab = "Parcor", xlab = "lag",
             ylim = c(-1, 1))
        abline(h = 0, lty = 1)
        abline(h = 2/sqrt(n), lty = 3)
        abline(h = -2/sqrt(n), lty = 3)
        par(mfrow = c(1, 1))

#  PARCOR OF SEASONAL.ORDER TIME(S) DIFFERENCED SEASONAL SERIES
        lag <- length(differenced.season$acov) - 1
        n <- differenced.season$n
        if (separate.graphics == TRUE)
            dev.new()
        par(mfrow = c(2, 1))
        plot((0:lag), differenced.season$acor, type = "h",
             main=paste(seasonal.order, "time(s) Differenced Seasonal Series"),
             cex.main = 0.9, ylab = "Autocorrelation", xlab = "lag",
             ylim = c(-1, 1))
        plot(differenced.season$parcor, type = "h", ylab = "Parcor",
             xlab = "lag", ylim = c(-1, 1))
        abline(h = 0, lty = 1)
        abline(h = 2/sqrt(n), lty = 3)
        abline(h = -2/sqrt(n), lty = 3)
    }

    par(oldpar)
}

decomp <- function(y, trend.order = 2, ar.order = 2, seasonal.order = 1,
                   period = 1, log = FALSE, trade = FALSE, diff = 1, miss = 0, 
                   omax = 99999.9, plot = TRUE)
{
    m1 <- trend.order
    m2 <- ar.order
    ilog <- 0
    if (log == TRUE)
        ilog <- 1
    itrade <- 0
    if (trade == TRUE)
        itrade <- 1

	if (is.null(tsp(y)) == TRUE) {
	  year <- 1
	  month <- 1
	} else if (is.null(tsp(y)) == FALSE) {
	  year <- start(y)[1]
	  month <- start(y)[2]
	  period <- tsp(y)[3]
    }
	
    sorder <- seasonal.order
    if (sorder != 0)
      if (period < 2) {
        warning("period must be greater than 1")
        sorder <- 0
      }

    if (itrade != 0) 
      if(period != 4 && period !=12) {
        warning("'trade = TRUE' is available only if period is 4 or 12")
        itrade <- 0
      }

    n <- length(y)
    ipar <- rep(0, 9)
    ipar[1] <- trend.order
    ipar[2] <- ar.order
    ipar[3] <- period
    ipar[4] <- sorder
    ipar[5] <- ilog
    ipar[6] <- itrade
    ipar[7] <- diff
    ipar[8] <- year
    ipar[9] <- month

    z <- .Fortran(C_decompf,
                  as.double(y),
                  as.integer(n),
                  as.integer(ipar),
                  trend = double(n),
                  seasonal = double(n),
                  ar = double(n),
                  trad = double(n),
                  noise = double(n),
                  para = double(26),
                  as.integer(miss),
                  as.double(omax),
                  ier = integer(1))

    ier <- z$ier
    if (ier != 0)
      stop("log-transformation cannot be applied to zeros and nagative numbers")

    season = z$seasonal
    ar = z$ar
    trading = z$trad
    noise = z$noise
    aic = z$para[1]
    lkhd = z$para[2]
    sigma2 = z$para[3]
    tau1 = z$para[4]
    tau2 = z$para[5]
    tau3 = z$para[6]
    if (ar.order == 0 && sorder == 0) {
       tau2 <- NULL
       tau3 <- NULL
    } else if (ar.order == 0 && sorder != 0) {
       tau3 <- NULL
    } else if (ar.order != 0 && sorder == 0) {
       tau3 <- NULL
    }
    arcoef = z$para[7:(6+m2)]
    tdf = z$para[(7+m2):(13+m2)]

	if (sorder == 0)
       season <- NULL

	if (ar.order == 0)
       ar <- NULL

    if (itrade == 0)
       trading <- NULL
		    
    if (miss > 0) {
      for (i in 1:n)
        if (y[i] > omax) 
          noise[i] <- NA
    } else if (miss < 0) {
      for (i in 1:n)
	    if (y[i] < omax)
          noise[i] <- NA
    }

    out <- list(trend = z$trend, seasonal = season, ar = ar, trad = trading, 
                noise = noise, aic = aic, lkhd = lkhd, sigma2 = sigma2, 
				tau1 = tau1, tau2 = tau2, tau3 = tau3, arcoef = arcoef, 
				tdf = tdf)
					   
    if (plot == TRUE) 
       plot.decomp(y, log, miss, omax, out)

    return(out)
}


#======================================================

plot.decomp <- function(tsdata, log, miss, omax, x, ...)

#======================================================

{
    ts.atr <- tsp(tsdata)
    n <- length(tsdata)
    trend <- x$trend
    noise <- x$noise
    if (is.null(ts.atr) == FALSE) {
      trend <- ts(trend, start = ts.atr[1], frequency = ts.atr[3])
      noise <- ts(noise, start = ts.atr[1], frequency = ts.atr[3])
    }

    oldpar <- par(no.readonly = TRUE)
    nc <- 2
    if (is.null(x$ar)  == FALSE)
      nc <- nc + 1
    if (is.null(x$seasonal)  == FALSE)
      nc <- nc + 1
    if (is.null(x$trad) == FALSE)
      nc <-  nc + 1
    if (nc > 3)
      par(mfrow = c((nc+1)/2, 2))
    if (nc <= 3)
      par(mfrow = c(nc, 1))

    y0 <- NULL
    missing <- rep(0, n)
    if (miss == 0) {
      y0 <- tsdata
    } else {
      for (i in 1:n)
        if (miss > 0) {
          if (tsdata[i] > omax) {
            missing[i] <- 1
          } else {
            y0 <- c(y0, tsdata[i])
          }
        } else {
          if (tsdata[i] < omax) {
            missing[i] <- -1
          } else {
            y0 <- c(y0, tsdata[i])
          }
        }
    }
	
    y0max <- max(y0)
    y0min <- min(y0)
    noise.nna <- noise[!is.na(noise)]

    if(log == FALSE) {
      ymax <- max(trend, y0max)
      ymin <- min(trend, y0min)
      yy <- tsdata
    } else {
      ymax <- max(trend, log(y0max)) 
      ymin <- min(trend, log(y0min))
      yy <- log(tsdata)
    }
	yd <- (ymax - ymin) / 10
    ymax <- ymax + yd
	ymin <- ymin - yd
	
###   Original and Trend
		 
    if (miss == 0) {
      plot(trend, type = 'l', col = 2, ylim = c(ymin, ymax), 
           main = "Original and Trend", xlab = "", ylab = "")
      par(new=TRUE)
      plot(yy, type = "l", ylim = c(ymin, ymax),
           main = "", xlab = "", ylab = "", xaxt = "n", yaxt = "n")

    } else {
      yy1 <- rep(NA, n)
      ymiss <- rep(NA, n)
      if (is.null(ts.atr) == FALSE) {
        ymiss <- ts(ymiss, start = ts.atr[1], frequency = ts.atr[3])
        yy1 <- ts(yy1, start = ts.atr[1], frequency = ts.atr[3])
      }
      for (i in 1:n)
        if (missing[i] != 0) {
          yy[i] <- NA
          ymiss[i] <- ymin
        } 
      if ((missing[1] == 0) && (missing[2]) != 0)
        yy1[1] <- yy[1]
      if ((missing[n-1] != 0) && (missing[n]) == 0)
        yy1[n] <- yy[n]
      for(i in 2:(n-1)) {
        if (missing[i] == 0)
          if ((missing[i-1] != 0) && (missing[i+1]) != 0)
            yy1[i] <- yy[i]
      }
      nmiss <- sum(ymiss, na.rm = TRUE)
      if (nmiss == 0) {
        message("\n There are no missing values")
      } else if (nmiss != 0) {
        ymax <- ymax + 0.5 * yd
        ymin <- ymin - 0.5 * yd
        ymiss <- ymiss - 0.5 * yd
	  }

      plot(trend, type = 'l', col = 2, ylim = c(ymin, ymax), 
           main = "Original and Trend", xlab = "", ylab = "")
      par(new=TRUE)

      plot(yy, type = "l", ylim = c(ymin, ymax),
           main = "", xlab = "", ylab = "", xaxt = "n", yaxt = "n")

      if (nmiss != 0) {
        par(new = TRUE)
        plot(yy1, type = "p", pch = 20, ylim = c(ymin, ymax),
             main = "", xlab = "", ylab = "", xaxt = "n", yaxt = "n", cex=0.6)

        par(new = TRUE)
        plot(ymiss, type="p", pch = 20, col = 4, ylim = c(ymin, ymax),
             main = "", xlab = "", ylab = "", xaxt = "n", yaxt = "n", cex=0.6)
        legend("topright", legend = "missing", col = 4, pch=20, bty = "n")
	  }
    }

###   Seasonal

    ymax <- max(x$ar, x$seasonal, x$trad, noise.nna)
    ymin <- min(x$ar, x$seasonal, x$trad, noise.nna)

    my <- max(ymax, abs(ymin)) * 1.5
    if (is.null(x$seasonal) == FALSE) {
      season <- x$seasonal
      if (is.null(ts.atr) == FALSE)
        season <- ts(season, start = ts.atr[1], frequency = ts.atr[3])
      plot(season, type = "l", main = "Seasonal", xlab = "", ylab = "",
           ylim = c(-my, my))
    }

###   Noise
	
    if (miss != 0) {
      noise1 <- rep(NA, n)
      if (is.null(ts.atr) == FALSE)
        noise1 <- ts(noise1, start = ts.atr[1], frequency = ts.atr[3])
      if ((missing[1] == 0) && (missing[2]) != 0)
        noise1[1] <- noise[1]
      if ((missing[n-1] != 0) && (missing[n]) == 0)
        noise1[n] <- noise[n]
      for(i in 2:(n-1)) {
        if (missing[i] == 0)
          if ((missing[i-1] != 0) && (missing[i+1]) != 0)
            noise1[i] <- noise[i]
      }
      plot(noise1, type = "p", pch = 20, main = "", xlab = "", ylab = "",
           ylim = c(-my, my), cex=0.6)
      par(new = TRUE)
    }

    plot(noise, type = "l", main = "Noise", xlab = "", ylab = "",
         ylim = c(-my, my))

###   AR

    if (is.null(x$ar) == FALSE) {
      ar <- x$ar
      if (is.null(ts.atr) == FALSE)
        ar <- ts(ar, start = ts.atr[1], frequency = ts.atr[3])
      plot(ar, type = "l", main = "AR component", xlab = "", ylab = "",
           ylim = c(-my, my))
    }

### Trading Day Effect
	
    if (is.null(x$trad) == FALSE) {
      trad <- x$trad
      if (is.null(ts.atr) == FALSE)
        trad <- ts(trad, start = ts.atr[1], frequency = ts.atr[3])
      plot(trad, type = "l", main = "Trading Day Effect", xlab = "",
           ylab = "", ylim = c(-my, my))
    }

    par(mfrow=oldpar$mfrow, new = oldpar$new)
}
