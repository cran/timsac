#####   TIMSAC74   #####

autoarmafit <- function (y, max.order = NULL)
{
    n <- length(y)
    if (is.null(max.order))
        max.order <- as.integer(2*sqrt(n))

    morder <- max.order
	morder1 <- morder + 1
    autcv <- autcor(y, morder, plot = FALSE)$acov   # covariance sequence
    z1 <- canarm(y, morder, plot = FALSE)
    inc <- 1               # total number of case
    arcoef1 <- z1$arcoef   # initial estimates of AR-coefficients
    macoef1 <- z1$macoef   # initial estimates of MA-coefficients
    p1 <- length(arcoef1)  # initial AR order
    q1 <- length(macoef1)  # initial MA order

    icst <- 190
    mmax <- morder
    lmax <- morder + icst + mmax
    nmax <- 25

    z <- .Fortran(C_autarm,
                  as.integer(n),
                  as.integer(morder1),
                  as.double(autcv),
                  as.integer(inc),
                  as.integer(p1),
                  as.double(arcoef1),
                  as.integer(q1),
                  as.double(macoef1),
                  newn = integer(1),
                  p = integer(nmax),
                  a = double(mmax*nmax),
                  q = integer(nmax),
                  b = double(mmax*nmax),
                  std = double(mmax*nmax),
                  v = double(nmax),
                  gr = double(mmax*nmax),
                  aic = double(nmax),
                  aicm = double(1),
                  pbest = integer(1),
                  qbest = integer(1),
                  as.integer(lmax),
                  as.integer(mmax),
                  as.integer(nmax))

    nc <- z$newn
    p <- z$p[1:nc]
    q <- z$q[1:nc]
    a <- array(z$a, dim = c(mmax, nmax))
    b <- array(z$b, dim = c(mmax, nmax))
    std <- array(z$std, dim = c(mmax, nmax))
    v <- rep(0, nc)
    gr <- array(z$gr, dim = c(mmax, nmax))
    aaic <- z$aic[1:nc]
    arma <- vector("list", nc)
    arstd <- list()
    mastd <- list()
    grad <- list()

    aicorder <- sort(aaic, index.return = TRUE)

    for(i in 1:nc) {
        j <- aicorder$ix[i]
        if (p[j] == 0) {
            ar <- NULL
        } else {
            ar <- -a[(1:p[j]), j]
        }
        if (q[j] == 0) {
            ma <- NULL
        } else {
            ma <- -b[(1:q[j]), j]
        }
        arma[[i]] <- list(arcoef = ar, macoef = ma)
        arstd[[i]] <- std[(q[j]+1):(p[j]+q[j]), j]
        mastd[[i]] <- std[1:q[j], j]
        v[i] <- z$v[j]
        grad[[i]] <- gr[(1:(p[j]+q[j])), j]
    }

    best.model <- list(arcoef = arma[[1]]$arcoef, macoef = arma[[1]]$macoef,
                       arorder = z$pbest, maorder = z$qbest)

    model <- list()
    for(i in 1:nc)
        model[[i]] <- list(arcoef = arma[[i]]$arcoef, macoef = arma[[i]]$macoef,
                           arstd = arstd[[i]], mastd = mastd[[i]], v = v[i],
                           aic = aicorder$x[i], grad = grad[[i]])

    autoarmafit.out <- list(best.model = best.model, model = model)
    class(autoarmafit.out) <- "autoarmafit"
    return(autoarmafit.out)
}

print.autoarmafit <- function(x, ...)
{
    nc <- length(x$model)
    z <- x$model

    for(i in 1:nc) {
        iq <- length(z[[i]]$arcoef)
        ip <- length(z[[i]]$macoef)

        cat(sprintf("\n\nCase No. %i\n", i))
        if (iq != 0) {
            cat("\n AR coefficient\tStandard deviation\n")
            for(j in 1:iq)
                cat(sprintf(" %f\t%f\n", z[[i]]$arcoef[j], z[[i]]$arstd[j]))
        } 
        if (ip != 0) { 
            cat("\n MA coefficient\tStandard deviation\n")
            for(j in 1:ip)
                cat(sprintf(" %f\t%f\n", z[[i]]$macoef[j], z[[i]]$mastd[j]))
        }
        cat(sprintf("\n AIC\t%f\n", z[[i]]$aic))
        cat(sprintf(" Innovation variance\t%f\n",z[[i]]$v))
        cat(" Final gradient")
        for(j in 1:(iq+ip))
            cat(sprintf("\t%e", z[[i]]$grad[j]))
        cat("\n")
    }

    zz <- x$best.model
    cat("\n\nBest ARMA model")
    cat(sprintf("\n AR coefficient (order = %i)", zz$arorder))
    for(j in 1:zz$arorder)
        cat(sprintf("\t%f", zz$arcoef[j]))
    cat(sprintf("\n MA coefficient (order = %i)", zz$maorder))
    for(j in 1:zz$maorder)
        cat(sprintf("\t%f", zz$macoef[j]))
    cat("\n\n")
}

armafit <- function (y, model.order)
{
    n <- length(y)
    lag <- as.integer(2*sqrt(n))     # maximum lag
    lag1 <- lag + 1

    inc <- 1             # total number of case
    p1 <- rep(0,2)
    q1 <- rep(0,2)

    autcv <- autcor(y, lag, plot = FALSE)$acov		# covariance sequence
    z1 <- canarm(y, lag, plot = FALSE)
    arcoef1 <- z1$arcoef	# initial estimates of AR-coefficients
    macoef1 <- z1$macoef	# initial estimates of MA-coefficients
    p1[1] <- length(arcoef1)	# initial AR order
    q1[1] <- length(macoef1)	# initial MA order

    if (length(model.order) == 2) {
        inc <- inc + 1
        p1[2] <- model.order[1]     # AR order to be fitted successively
        q1[2] <- model.order[2]     # MA order to be fitted successively
    }

    icst <- 190
    mmax <- 50
    lmax <- lag + icst + mmax
    nmax <- 25	                # the limit of the total cases

    z <- .Fortran(C_autarm,
                  as.integer(n),
                  as.integer(lag1),
                  as.double(autcv),
                  as.integer(inc),
                  as.integer(p1),
                  as.double(arcoef1),
                  as.integer(q1),
                  as.double(macoef1),
                  newn = integer(1),
                  p = integer(nmax),
                  a = double(mmax*nmax),
                  q = integer(nmax),
                  b = double(mmax*nmax),
                  std = double(mmax*nmax),
                  v = double(nmax),
                  gr = double(mmax*nmax),
                  aic = double(nmax),
                  aicm = double(1),
                  pbest = integer(1),
                  qbest = integer(1),
                  as.integer(lmax),
                  as.integer(mmax),
                  as.integer(nmax))

    nc <- z$newn
    p <- z$p[1:nc]
    q <- z$q[1:nc]
    a <- array(z$a, dim = c(mmax, nmax))
    b <- array(z$b, dim = c(mmax, nmax))
    std <- array(z$std, dim = c(mmax, nmax))
    gr <- array(z$gr, dim = c(mmax, nmax))
    aaic <- z$aic[1:nc]
    aic <- NULL
    arcoef <- NULL
    macoef <- NULL
    arstd <- NULL
    mastd <- NULL
    grad <- NULL

    aicorder <- sort(aaic, index.return = TRUE)
    arorder <- model.order[1]
    maorder <- model.order[2]

    for(i in 1:nc) {
        if (is.null(aic) == FALSE)
            break
        j <- aicorder$ix[i]
        if (p[j]==arorder  &&  q[j]==maorder) {
            arcoef <- -a[(1:p[j]), j]
            macoef <- -b[(1:q[j]), j]
            arstd <- std[(q[j]+1):(p[j]+q[j]), j]
            mastd <- std[1:q[j], j]
            v <- z$v[j]
            aic <- z$aic[j]
            grad <- gr[(1:(p[j]+q[j])), j]
        }
    }

    armafit.out <- list(arcoef = arcoef, macoef = macoef, arstd = arstd,
                        mastd = mastd, v = v, aic = aic, grad = grad)
    return(armafit.out)
}


bispec <- function (y, lag = NULL, window = "Akaike", log = FALSE, plot = TRUE)
{
    if (window != "Akaike" && window != "Hanning")
        stop("allowed windows are 'Akaike' or 'Hanning'")

    n <- length(y)
    if (is.null(lag))
        lag <- as.integer(2*sqrt(n))      # maximum lag
    lag1 <- lag + 1

    z1 <- thirmo(y, lag, plot = FALSE)
    cv <- z1$acov			  # autocovariance

    tmnt <- array(0, dim = c(lag1, lag1))
    for(i in 1:lag1)
        tmnt[i, 1:i] <- z1$tmomnt[[i]]    # third order moments

    z <- .Fortran(C_bispecf,
                  as.integer(n),
                  as.integer(lag),
                  as.double(cv),
                  as.double(tmnt),
                  pspec1 = double(lag1),
                  pspec2 = double(lag1),
                  sig = double(lag1),
                  ch = double(lag1*lag1),
                  br = double(lag1*lag1),
                  bi = double(lag1*lag1),
                  rat = double(1))

    if (window == "Akaike")
        pspec <- z$pspec2
    if (window == "Hanning")
       pspec <- z$pspec1

    if (plot == TRUE) {
        x <- rep(0, lag1)
        for(i in 1:lag1)
            x[i] <- (i - 1) / (2 * lag)
        ylab = paste("Spectrum smoothing by ",window," window")
        if (log == TRUE)
            plot(x, pspec, type = "l", log = "y", ylab = ylab,
                 xlab = "Frequency")
        if (log == FALSE)
            plot(x, pspec, type = "l", ylab = ylab, xlab = "Frequency")
    }

    bispec.out <- list(pspec = pspec, sig = z$sig,
                       cohe = array(z$ch, dim = c(lag1, lag1)),
                       breal = array(z$br, dim = c(lag1, lag1)),
                       bimag = array(z$bi, dim = c(lag1, lag1)),
                       exval = z$rat)
    return(bispec.out)
}


canarm <- function (y, lag = NULL, max.order = NULL, plot = TRUE)
{
    n <- length(y)
    if (is.null(lag))
        lag <- as.integer(2*sqrt(n))
    morder <- lag

    z2 <- autcor(y, morder, plot = FALSE)
    autcv <- z2$acov

    if (is.null(max.order))
        max.order <- morder
    mmax <- max.order
	mmax1 <- mmax + 1
	mmax2 <- mmax * mmax
	mmax3 <- mmax2 * mmax
	
#    ifpl <- as.integer(3.0*sqrt(n))
    ifpl <- min(mmax, morder)
    l1 <- ifpl + 1			# upper limit of the model order +1
    nmax <- 2 * morder + 1

    z1 <- .Fortran(C_canarmf,
                   as.integer(n),
                   as.integer(morder+1),
                   as.double(autcv),
                   arcoef = double(nmax),
                   as.integer(l1),
                   v = double(mmax1),
                   aic = double(mmax1),
                   oaic = double(1),
                   mo = integer(1),
                   parcor = double(mmax),
                   nc = integer(1),
                   m1 = integer(mmax),
                   m2 = integer(mmax),
                   w = double(mmax3),
                   z = double(mmax2),
                   Rs = double(mmax2),
                   chi = double(mmax2),
                   ndt = double(mmax2),
                   dic = double(mmax2),
                   dicm = double(mmax),
                   po = integer(mmax),
                   k = integer(1),
                   b = double(mmax),
                   l = integer(1),
                   a = double(mmax),
                   as.integer(mmax),
                   as.integer(nmax))

    ar <- z1$arcoef
    mo <- z1$mo
    parcor <- z1$parcor[1:(l1-1)]
    nc <- z1$nc
    m1 <- z1$m1[1:nc]
    m2 <- z1$m2[1:nc]
    w <- array(z1$w, dim = c(mmax, mmax, nc))
    z <- array(z1$z, dim = c(mmax, nc))
    Rs <- array(z1$Rs, dim = c(mmax, nc))
    chi <- array(z1$chi, dim = c(mmax, nc))
    ndt <- array(z1$ndt, dim = c(mmax, nc))
    dicp <- array(z1$dic, dim = c(mmax, nc))
    cweight <- list()
    cR <- list()
    Rsquar <- list()
    chisquar <- list()
    ndf <- list()
    dic <- list()
    for(i in 1:nc) {
        cweight[[i]] <- w[(1:m1[i]), (1:m1[i]), i]
        cR[[i]] <- z[(1:m1[i]), i]
        Rsquar[[i]] <- Rs[(1:m1[i]), i]
        chisquar[[i]] <- chi[(1:m1[i]), i]
        ndf[[i]] <- ndt[(1:m1[i]), i]
        dic[[i]] <- dicp[(1:m1[i]), i]
    }
    k <- z1$k
    l <- z1$l

    if (plot == TRUE) {
        plot(parcor, type = "h", ylab = "Partial Autocorrelation", xlab = "Lag")
        abline(h = 0, lty = 1)
#        abline(h = 1/sqrt(n), lty = 3)
#        abline(h = -1/sqrt(n), lty = 3)
        abline(h = 2/sqrt(n), lty = 3)
        abline(h = -2/sqrt(n), lty = 3)
    }

    canarm.out <- list(arinit = -ar[1:mo], v = z1$v[1:l1], aic = z1$aic[1:l1],
                       aicmin = z1$oaic, order.maice = mo, parcor = parcor,
                       nc = nc, future = m1, past = m2, cweight = cweight, 
                       canocoef = cR, canocoef2 = Rsquar, chisquar = chisquar,
                       ndf = ndf, dic = dic, dicmin = z1$dicm[1:nc],
                       order.dicmin = z1$po[1:nc], arcoef = -z1$b[1:k],
                       macoef = -z1$a[1:l])
    return(canarm.out)
}


canoca <- function (y)
{
    n <- nrow(y)                        # length of data
    d <- ncol(y)                        # dimension of y(i)
    lag <- as.integer(2*sqrt(n))        # maximum lag
    lag1 <- lag + 1
# inw(k)=j means that
# the k-th component of y(i) is the j-th component of the original record z(i)
    inw <- c(1:d)

    z2 <- mulcor(y, lag, plot = FALSE)
    cov <- z2$cov                       # covariance matrix

    lmax <- 12
    mj0 <- lmax + 1
    mj1 <- mj0 * d

    l <- 0                              # upper bound of AR-order
    aic <- rep(0, mj0)                  # AIC
    oaic <- 0                           # minimum AIC
    mo <- 0                             # MAICE AR-model order                 
    v <- array(0, dim = c(d, d))        # innovation variance
    ac <- array(0, dim = c(mj0, d, d))  # autoregressive coefficients
    nc <- 0                             # number of cases
    m1 <- rep(0, mj1)                   # number of variable in the future set
    m2 <- rep(0, mj1)                   # number of variables in the past set
    w <- array(0, dim = c(mj1, mj1, mj1))  # future set canonical weight
    z <- array(0, dim = c(mj1, mj1))    # canonical R
    Rs <- array(0, dim = c(mj1, mj1))   # R-squared
    chi <- array(0, dim = c(mj1, mj1))  # chi-square
    ndt <- array(0, dim = c(mj1, mj1))  # N.D.F
    dic <- array(0, dim = c(mj1, mj1))  # DIC(=CHI**2-2*D.F)
    dicm <- rep(0, mj1)                 # minimum DIC
    po <- rep(0, mj1)                   # order of minimum DIC
    f <- array(0, dim = c(mj1, mj1))    # transition matrix F
    k <- 0                           # number of structual characteristic vector
    nh <- rep(0, mj1)                   # structual characteristic vector
    g <- array(0, dim = c(mj1, d))      # input matrix
    ivf <- 0                            # number of vector vf
    vf <- rep(0, mj1*mj1)               # F matrix in vector form

    z1 <- .Fortran(C_canocaf,
                   as.integer(d),
                   as.integer(inw),
                   as.integer(n),
                   as.integer(lag1),
                   as.integer(d),
                   as.double(cov),
                   l = integer(1),
                   aic = double(mj0),
                   oaic = double(1),
                   mo = integer(1),
                   v = double(d*d),
                   ac = double(d*d*mj0),
                   nc = integer(1),
                   m1 = integer(mj1),
                   m2 = integer(mj1),
                   w = double(mj1*mj1*mj1),
                   z = double(mj1*mj1),
                   Rs = double(mj1*mj1),
                   chi = double(mj1*mj1),
                   ndt = integer(mj1*mj1),
                   dic = double(mj1*mj1),
                   dicm = double(mj1),
                   po = integer(mj1),
                   f = double(mj1*mj1),
                   k = integer(1),
                   nh = integer(mj1),
                   g = double(d*mj1),
                   ivf = integer(1),
                   vf = double(mj1*mj1),
                   as.integer(lmax),
                   as.integer(mj0),
                   as.integer(mj1))

    l <- z1$l
    mo <- z1$mo
    ac <- array(z1$ac, dim = c(mj0, d, d))
    arcoef <- array(, dim = c(d, d, mo))
    for(i in 1:mo)
        arcoef[, , i] <- -ac[i, (1:d), (1:d)]
    nc <- z1$nc
    m1 <- z1$m1[1:nc]
    m2 <- z1$m2[1:nc]
    w <- array(z1$w, dim = c(mj1, mj1, mj1))
    z <- array(z1$z, dim = c(mj1, mj1))
    Rs <- array(z1$Rs, dim = c(mj1, mj1))
    chi <- array(z1$chi, dim = c(mj1, mj1))
    ndt <- array(z1$ndt, dim = c(mj1, mj1))
    dicp <- array(z1$dic, dim = c(mj1, mj1))
    f <- array(z1$f, dim = c(mj1, mj1))
    cweight <- list()
    cR <- list()
    Rsquar <- list()
    chisquar <- list()
    ndf <- list()
    dic <- list()
    matF <- list()
    for(i in 1:nc)  {
        cweight[[i]] <- w[(1:m1[i]), (1:m1[i]), i]
        cR[[i]] <- z[(1:m1[i]), i]
        Rsquar[[i]] <- Rs[(1:m1[i]), i]
        chisquar[[i]] <- chi[(1:m1[i]), i]
        ndf[[i]] <- ndt[(1:m1[i]), i]
        dic[[i]] <- dicp[(1:m1[i]), i]
        matF[[i]] <- f[1:(m1[i]-1), i]
    }
    k <- z1$k
    g <- array(z1$g, dim = c(mj1, d))
    ivf <- z1$ivf

    canoca.out <- list(aic = z1$aic[1:(l+1)], aicmin = z1$oaic,
                       order.maice = mo, v = array(z1$v, dim = c(d, d)),
                       arcoef = arcoef, nc = nc, future = m1, past = m2,
                       cweight = cweight, canocoef = cR, canocoef2 = Rsquar,
                       chisquar = chisquar, ndf = ndf, dic = dic,
                       dicmin = z1$dicm[1:nc], order.dicmin = z1$po[1:nc],
                       matF = matF, vectH = z1$nh[1:k],
                       matG = g[(1:k),(1:d)], vectF = z1$vf[1:ivf])
    return(canoca.out)
}


covgen <- function (lag, f, gain, plot = TRUE)
{
    k <- length(f)	# number of data points
	lag1 <- lag + 1

    z <- .Fortran(C_covgenf,
                  as.integer(lag),
                  as.integer(k),
                  as.double(f),
                  as.double(gain),
                  acov = double(lag1),
                  acor = double(lag1))

    if (plot == TRUE) {
        plot((0:lag), z$acor, type = "h", ylab = "Auto Covariance Normalized",
             xlab = "Lag")
      abline(h = 0, lty = 1)
    }

    covgen.out <- list(acov = z$acov, acor = z$acor)
    return(covgen.out)
}


markov <- function (y)
{
    n <- nrow(y)                  # length of data
    d <- ncol(y)                  # dimension of the observation vector
    lag <- as.integer(2*sqrt(n))  # maximum lag
    lag1 <- lag + 1
# output control (0:ARMA coefficients, 1:for SIMCON input, 2: for both)
    icont <- 2

    z1 <- mulcor(y, lag, plot = FALSE)
    cov <- z1$cov		  # covariance matrix

    z2 <- canoca(y)
    nh <- z2$vectH                # structural characteristic vector
    k <- length(nh)               # dimension of the state vector
# initial estimate of the vector of free parameters in F
    vectF <- z2$vectF
    nvf <- length(vectF)          # length of vector vf
# initial estimates of the free parameters in G
    matGi <- rep(z2$matG, 1)

    mj3 <- max(lag1, 100)
    mj4 <- nvf + k * d
    mj6 <- 2 * (k + d) - 1
    mj7 <- as.integer((k-1)/d+1)

    z <- .Fortran(C_markovf,
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
                  id = integer(k),
                  ir = integer(k),
                  ij = integer(d),
                  ik = integer(d),
                  ngr = integer(1),
                  gr = double(mj4),
                  a1 = double(k*k),
                  a = double(k*k),
                  b = double(k*d),
                  vd = double(mj4*mj4),
                  iqm = integer(1),
                  bm = double(d*d*mj7),
                  au = double(d*d*mj7),
                  zz = double(d*d*mj7),
                  v = double(d*d),
                  aic = double(1),
                  as.integer(mj3),
                  as.integer(mj4),
                  as.integer(mj6),
                  as.integer(mj7))

    ngr <- z$ngr
    vd <- array(z$vd, dim = c(mj4, mj4))
    iqm <- z$iqm
    bm <- array(z$bm, dim = c(d, d, mj7))
    au <- array(z$au, dim = c(d, d, mj7))
    zz <- array(z$zz, dim = c(d, d, mj7))

    arcoef <- array(-bm, dim = c(d, d, iqm))
    macoef <- array(-zz, dim = c(d, d, iqm-1))

    markov.out <- list(id = z$id, ir = z$ir, ij = z$ij, ik = z$ik,
                       grad = z$gr[1:ngr],
                       matFi = array(z$a1, dim = c(k, k)),
                       matF = array(z$a, dim = c(k, k)),
                       matG = array(z$b, dim = c(k, d)),
                       davvar = vd[(1:ngr),(1:ngr)], arcoef = arcoef,
                       impulse = array(au, dim = c(d, d, iqm-1)),
                       macoef = macoef, v = array(z$v, dim = c(d, d)),
                       aic = z$aic)
    return(markov.out)
}


nonst <- function (y, span, max.order = NULL, plot = TRUE)
{
    n <- length(y)                          # length of data
    if (span < 1)
        span <- n
    ns <- as.integer(n/span)
    if (is.null(max.order))
        max.order <- as.integer(2*sqrt(n))  # highest order of AR model
    morder <- max.order

    z <- .Fortran(C_nonstf,
                  as.integer(n),
                  as.integer(span),
                  as.double(y),
                  as.integer(ns),
                  as.integer(morder),
                  p = integer(ns),
                  coef = double(ns*morder),
                  v = double(ns),
                  aic = double(ns),
                  daic21 = double(ns),
                  daic = double(ns),
                  ks = integer(ns),
                  ke = integer(ns),
                  pspec = double(121*ns))

    coef <- array(z$coef, dim = c(morder, ns))
    arcoef <- list() 
    for(i in 1:ns)
        arcoef[[i]] <- coef[1:z$p[i], i]
    pspec <- array(z$pspec, dim = c(121, ns))

    if (plot == TRUE) {
        oldpar <- par(no.readonly = TRUE)
        x <- rep(0, 121)
        for(i in 1:121)
            x[i] <- (i - 1) / 240
        par(mfrow = c(ns, 1))
        for(i in 1:ns) {
            plot(x, pspec[,i], type = "l",
                 main = paste("y(", z$ks[i], "),...,y(", z$ke[i], ")"),
                 xlab = "Frequency", ylab = "Power Spectrum")
        }
        par(oldpar)
        nonst.out <- list(ns = ns, arcoef = arcoef, v = z$v, aic = z$aic,
                          daic21 = z$daic21, daic = z$daic, init = z$ks,
                          end = z$ke)
    } else {
        nonst.out <- list(ns = ns, arcoef = arcoef, v = z$v, aic = z$aic,
                          daic21 = z$daic21, daic = z$daic, init = z$ks,
                          end = z$ke, pspec = pspec)
    }
    class(nonst.out) <- "nonst"
    return(nonst.out)
}


print.nonst <- function(x, ...)
{
    for(i in 1:x$ns) {
        cat("\n\n..........  Current model  ..........\n\n")
        cat(sprintf(" This model was fitted to the data  y( %i ),...,y( %i )\n",
 x$init[i],x$end[i]))
        cat(sprintf(" Innovation variance = %e\n", x$v[i]))
        cat(sprintf(" aic = %f\n", x$aic[i]))
        if (i > 1) {
            cat(sprintf(" aic2-aic1 = %f\n", x$daic21[i]))
            nd <- x$end[i] - x$init[i] + 1
        cat(sprintf(" daic/%i = %f\n", nd, x$daic[i]))
        }
        mf <- length(x$arcoef[[i]])
        cat(sprintf("\n Autoregressive coefficient ( order %i )\n", mf))
        print(x$arcoef[[i]])
    }
}


prdctr <-
function (y, r, s, h, arcoef, macoef = NULL, impulse = NULL, v, plot = TRUE)
{
    if (is.array(y)) {
        n <- nrow(y)    # length of data
        d <- ncol(y)    # dimension of vector y(i)
    } else {
        n <- length(y)
        d <- 1
        y <- array(y, dim = c(n, 1))
    }

    if (is.null(arcoef))
        stop("'arcoef' must be numeric")
    if(is.array(arcoef)) {        # AR-coefficient matrices
        p <- length(arcoef)[3]
        arcoef <- -arcoef
    } else {
        p <- length(arcoef)
        arcoef <- array(-arcoef, dim = c(1, 1, p))
    }

    if (is.null(macoef) && is.null(impulse))
        stop("'macoef' or ' impulse' must be numeric")
    jsw <- 0
    if (is.null(macoef)) {
        jsw <- 1
        if (is.array(impulse)) {   # impulse response matrices
            q <- dim(impulse)[3]
        } else {
            q <- length(impulse)
            impulse <- array(impulse, dim = c(1, 1, q))
        }
        macoef <- array(0, dim = c(d, d, q))
    } else {
        if (is.array(macoef)) {   # MA-coefficient matrices
            q <- dim(macoef)[3]
            macoef <- -macoef
        } else {
            q <- length(macoef)
            macoef <- array(-macoef, dim = c(1, 1, q))
        }
        if (is.null(impulse))
            impulse <- array(0, dim = c(d, d, q))
    }
	
    k <- (s+h) * d

    z <- .Fortran(C_prdctrf,
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
                  as.double(impulse),
                  as.double(v),
                  yreal = double(k),
                  yori = double((h+1)*d),
                  ypre = double(k),
                  ys = double(n*d),
                  z1 = double(k),
                  z2 = double(k),
                  z3 = double(k),
                  zz1 = double(k),
                  zz2 = double(k),
                  zz3 = double(k))

# yreal <- array(z$yreal, dim = c(s+h, d))
# yori <- array(z$yori, dim = c(h+1, d))
# for(i in 1:(n-s+1)) yreal[s+i-1, ] <- yori[i, ]
# for(j in 1:d) yreal[, j] <- yreal[(1:n),j]

    predct <- array(z$ypre, dim = c(s+h, d))
    for(i in 1:(r-1))
        predct[i, ] <- NA

# ys <- array(z$ys, dim = c(n, d))
# for(i in 1:(r-1)) ys[i,] <- NA
# for(i in s:n) ys[i,] <- NA
    ys <- array(NA, dim = c(n, d))
    for(j in 1:d)
        for(i in r:n)
            ys[i, j] <- y[i, j] - predct[i, j]

    pstd <- array(z$z1, dim = c(s+h, d))
    pstd[1:s-1, ] <- NA
    p2std <- array(z$z2, dim = c(s+h, d))
    p2std[1:s-1, ] <- NA
    p3std <- array(z$z3, dim = c(s+h, d))
    p3std[1:s-1, ] <- NA
    mstd <- array(z$zz1, dim = c(s+h, d))
    mstd[1:s-1, ] <- NA
    m2std <- array(z$zz2, dim = c(s+h, d))
    m2std[1:s-1, ] <- NA
    m3std <- array(z$zz3, dim = c(s+h, d))
    m3std[1:s-1, ] <- NA

    if (plot == TRUE) {
        par(mfrow = c(1, d))
        for(i in 1:d) {
            ymin <- min(y[(1:n), i], predct[r:(s+h), i])
            ymax <- max(y[(1:n), i], predct[r:(s+h), i])
            plot(y[, i], type = "l", xlim = c(0, (s+h)),
 ylim = c(ymin, ymax),  xlab = "n", main = "real data & predicted values")
            par(new=TRUE)
            plot(predct[, i], type = "l", xlim = c(0, (s+h)),
 ylim = c(ymin, ymax), col = "red", xlab = "", ylab = "")
        }
        par(mfrow = c(1, 1))
    }

    prdctr.out <- list(predct = predct, ys = ys, pstd = pstd, p2std = p2std,
                       p3std = p3std, mstd = mstd, m2std = m2std, m3std = m3std)

    class(prdctr.out) <- "prdctr"
    return(prdctr.out)
}

print.prdctr <- function(x, ...)
{
    n <- dim(x$ys)[1]
    d <- dim(x$ys)[2]
    sh <- dim(x$predct)[1]
    for(i in n:1) {
        if (is.na(x$ys[i,1]) == FALSE)
            r <- i
        if (is.na(x$pstd[i,1]) == FALSE)
            s <- i
    }

    for(id in 1:d) {
        cat(sprintf("\n\n d = %i\n\n", id))
        cat(" n\tpredct\t\tys\t\tpstd\t\tp2std\t\tp3std\n")
        cat("\t\t\t\t\tmstd\t\tm2std\t\tm3std\n")
        for(i in r:(s-1))
            cat(sprintf(" %i\t%f\t%f\n", i, x$predct[i,id], x$ys[i,id]))
        for(i in s:n) {
            cat(sprintf(" %i\t%f\t%f\t%f\t%f\t%f\n",
                i, x$predct[i,id], x$ys[i,id], x$pstd[i,id], x$p2std[i,id],
                x$p3std[i,id]))
            cat(sprintf("\t\t\t\t\t%f\t%f\t%f\n", x$mstd[i,id], x$m2std[i,id],
                x$m3std[i,id]))
        }
        for(i in (n+1):sh) {
            cat(sprintf(" %i\t\t\t%f\t%f\t%f\t%f\n", i, x$pstd[i,id],
                x$predct[i,id], x$p2std[i,id], x$p3std[i,id]))
            cat(sprintf("\t\t\t\t\t%f\t%f\t%f\n", x$mstd[i,id], x$m2std[i,id],
                x$m3std[i,id]))
        }
    }
}


simcon <- function (span, len, r, arcoef, impulse, v, weight)
{
    arcoef <- -arcoef        # matrices of autoregressive coefficients
    d <- dim(arcoef)[1]	     # dimension of Y(I)
    k <- dim(arcoef)[3]	     # order of the process
    if (r > d) {
        stop("r must be less than or equal to d (dimension of a vector)")
    }

    z <- .Fortran(C_simconf,
                  as.integer(d),
                  as.integer(k),
                  as.integer(span),
                  as.integer(len),
                  as.integer(r),
                  as.double(arcoef),
                  as.double(impulse),
                  as.double(v),
                  as.double(weight),
                  bc = double(k*d*r),
                  bd = double(k*d*(d-r)),
                  g = double(k*d*r),
                  av = double(d),
                  si = double(d),
                  s2 = double(d))

    simcon.out <- list(gain = array(z$g, dim = c(r, k*d)), ave = z$av[1:d],
                       var = z$si[1:d], std = z$s2[1:d],
                       bc = array(z$bc, dim = c(k*d, r)),
                       bd = array(z$bd, dim = c(k*d, d-r)))
    return(simcon.out)
}


thirmo <- function (y, lag = NULL, plot = TRUE)
{
    n <- length(y)
    if (is.null(lag))
        lag <- as.integer(2*sqrt(n))  # maximum lag
    lag1 <- lag + 1

    z <- .Fortran(C_thirmof,
                  as.integer(n),
                  as.integer(lag),
                  as.double(y),
                  mean = double(1),
                  acov = double(lag1),
                  acor = double(lag1),
                  mnt = double(lag1*lag1))

    mnt <- array(z$mnt, dim = c(lag1, lag1))
    tmomnt <- list()
    for(i in 1:lag1)
        tmomnt[[i]] <- mnt[i, 1:i]

    if (plot == TRUE) {
        plot((0:lag), z$acor, type = "h", ylab = "Auto Covariance Normalized",
             xlab = "Lag")
        abline(h = 0, lty = 1)
    }

    thirmo.out <- list(mean = z$mean, acov = z$acov, acor = z$acor,
                       tmomnt = tmomnt)
    return(thirmo.out)
}

#--------------------------------------------------------------------------
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
        } else {
            ni <- NROW(init)
            if (ni != p) 
                stop("length of 'init' must equal length of 'filter'")
            if (NCOL(init) != 1 && NCOL(init) != nser) 
                stop(gettextf("'init'; must have 1 or %d cols", nser),
                              domain = NA)
            if (!is.matrix(init)) 
                init <- matrix(init, p, nser)
        }
        xx <- matrix(NA,nd,nser)
        for(i in 1:p) xx[i, ] <- init[i, ]
        tar <- array(0, dim = c(p,nser, nser))
        for(i in 1:p) tar[i, , ] <- t(filter[, , i])
        i <- p + 1
        while (i <= nd) {
          xx[i, ] = x[i-p, ]
          y[i, ] <- x[i-p, ]
          for(j in 1:p) y[i, ] <- y[i, ] - xx[i - j, ] %*% tar[j, , ]
             i <- i + 1
        }
    } else {
        if (missing(init)) {
            init <- matrix(0, p, nser)
        } else {
            ni <- NROW(init)
            if (ni != p) 
                stop("length of 'init' must equal length of 'filter'")
            if (NCOL(init) != 1 && NCOL(init) != nser) 
                stop(gettextf("'init'; must have 1 or %d cols", nser),
                              domain = NA)
            if (!is.matrix(init)) 
                init <- matrix(init, p, nser)
        }
        for(i in 1:p)
            y[i, ] <- init[i, ]
        tar <- array(0, dim = c(p, nser, nser))
        for(i in 1:p)
            tar[i, , ] <- t(filter[, , i])
        i <- p + 1
        while (i <= nd) {
            y[i, ] <- x[i-p, ]
            for(j in 1:p)
                y[i, ] <- y[i, ] + y[i - j, ] %*% tar[j, , ]
            i <- i + 1
        }
    }
    y <- y[(p+1):nd, ]
    y <- drop(y)
    tsp(y) <- xtsp
    class(y) <- "mts"
    return(y)
}
