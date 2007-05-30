      SUBROUTINE TVARF(Y, N, M, K, NOBS, IOPT, NOUT, LOUT, TAU20, DELTA,
     *                   TAUMAX, SIG2, FF, AIC, TAR, PAR )
C
      INCLUDE 'timsac_f.h'
C
C     PROGRAM 13.2  TVAR
C
C  ...  Time Varying AR model  ...
C
C     Inputs:
C        M:       AR Order (M =< 10)
C        K:       Trend Order
C        NOBS:    Local Stationary Span (N/NOBS =< NDIM)
C        IOPT:    Search Method
C        NOUT:    Number of Outliers
C          LOUT(I):  Position of I-th Outlier
C     Parameters:
C        NMAX:    Adjustable dimension of Y
C        MJ:      Adjustable dimension of XF, VF, etc.
C        NDIM:    Adjustable dimension OF VFS, VSS, etc.
C     @TEST.FILTER2:  SEP.08,1990
C
cc      !DEC$ ATTRIBUTES DLLEXPORT :: TVARF
C
cc      PARAMETER( NMAX=3000,MJ=20,NDIM=200 )
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION  Y(N)
cc      DIMENSION  A(2,MJ), MM(MJ), LOUT(10), LLOUT(NDIM)
cc      DIMENSION  XPS(MJ,NDIM), XFS(MJ,NDIM), XSS(MJ,NDIM)
cc      DIMENSION  VPS(MJ,MJ,NDIM), VFS(MJ,MJ,NDIM), VSS(MJ,MJ,NDIM)
cc      DIMENSION  XF(MJ), VF(MJ,MJ)
      DIMENSION  LOUT(NOUT), TAR(M,N/NOBS), PAR(M,N/NOBS)
      DIMENSION  A(K,M), MM(M), LLOUT(N/NOBS)
      DIMENSION  XPS(M*K,N/NOBS), XFS(M*K,N/NOBS), XSS(M*K,N/NOBS)
      DIMENSION  VPS(M*K,M*K,N/NOBS), VFS(M*K,M*K,N/NOBS)
      DIMENSION  XF(M*K), VF(M*K,M*K)
      DATA  OUTMIN, OUTMAX /-1.0D30, 1.0D30/
cc      DATA  LLOUT /NDIM*0/
C
      MJ = M*K
      NDIM = N/NOBS
      DO 5 I=1,N/NOBS
    5 LLOUT(I) = 0
C
cc      READ( 5,* )  M, K, NOBS, IOPT, NOUT
cc      IF( IOPT.EQ.1 )  READ(5,*)  TAU20, DELTA
      IF( NOUT.GT.0 )  THEN
cc        READ(5,*)  (LOUT(I),I=1,NOUT)
        DO 10 I=1,NOUT
        J = LOUT(I)/NOBS
        IF( J*NOBS-LOUT(I).GT.NOBS/2 )  J = J+1
   10   LLOUT(J) = 1
      END IF
      INUM = 19
      IF( IOPT.EQ.0 )  INUM = 9
C
C  ...  Read Time Series  ...
C
cc      CALL  READTS( 1,Y,N )
C
      FMAX  = -1.0D30
C
      CALL  SETCAR( M,K,A,MM )
      DO 100  II=1,INUM
      IF( IOPT.NE.0 )  TAU2 = TAU20 + DELTA*(II-9)
      IF( IOPT.EQ.0 .AND. K.EQ.1 )  TAU2 =10.0D0**(-II)
      IF( IOPT.EQ.0 .AND. K.GT.1 )  TAU2 =10.0D0**(-II-1)
      CALL  ISTCAR( M,K,MJ,XF,VF )
C
C  ...  Log-Likelihood Computation  ...
C
cc      CALL  FILTR2( Y,XF,VF,TAU2,M,K,N,NOBS,MJ,1,LLOUT,
      CALL  FILTR2K( Y,XF,VF,TAU2,M,K,N,NOBS,MJ,1,LLOUT,
     *              OUTMIN,OUTMAX,VFS,VPS,XFS,XPS,FLK,SIG2 )
      IF( FLK.GT.FMAX )  THEN
         FMAX  = FLK
         TAUMAX = TAU2
         SIG2M  = SIG2
      END IF
cc  100 WRITE(6,600)  TAU2, SIG2, FLK
  100 CONTINUE
      AIC = -2*FMAX + 2*(M+2)
cc      WRITE(6,*)  M, K
cc      WRITE(6,610)  TAUMAX, SIG2M, FMAX, AIC
cc      WRITE(6,*)  FF,SIG2
C     STOP
C
C  ... Fixed Interval Smoother  ...
C
      CALL  ISTCAR( M,K,MJ,XF,VF )
cc      CALL  FILTR2( Y,XF,VF,TAUMAX,M,K,N,NOBS,MJ,NDIM,LLOUT,
      CALL  FILTR2K( Y,XF,VF,TAUMAX,M,K,N,NOBS,MJ,NDIM,LLOUT,
     *              OUTMIN,OUTMAX,VFS,VPS,XFS,XPS,FF,SIG2 )
      NN = N/NOBS
cc      CALL  SMOTH1( A,MM,K,M,1,NN,NN,MJ,VFS,VPS,VSS,XFS,XPS,XSS )
      CALL  SMOTH1( A,MM,K,M,1,NN,NN,MJ,VFS,VPS,XFS,XPS,XSS )
C
C  ...  Plot PARCOR and print AR Coefficients  ...
C
cc      CALL  PTCAR( XSS,XPS,N,NOBS,MJ,M,K,TAUMAX,SIG2,FF,AIC )
      CALL  PTCAR( XSS,N,NOBS,MJ,M,K,TAR,PAR )
cc      CALL  PRCAR( XSS,N,NOBS,MJ,M,K,TAUMAX,SIG2,FF,AIC )
C
cc      STOP
      RETURN
  600 FORMAT( 1H ,5X,F15.10,F12.6,F13.5 )
  610 FORMAT( 1H ,'TAUMAX =',F15.10,3X,'SIG2 =',F15.10,3X,
     *            'FF =',F13.5,3X,'AIC =',F13.4 )
      E N D
      SUBROUTINE  SETCAR( M,K,A,MM )
C
C  ...  State space model for trend estimation  ...
C
C     Input:
C       M:   AR Order
C       K:   Trend Order
C     Outputs:
C       A:   Parameter of the Transition Matrix F
C       MM:  Dimension of the Component Model
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION  A(2,*), MM(*)
C
      IF( K.EQ.1 )  THEN
        DO 10 I=1,M
   10   A(1,I) = 1.0D0
      END IF
      IF( K.EQ.2 )  THEN
        DO 20 I=1,M
        A(1,I) = 2.0D0
   20   A(2,I) =-1.0D0
      END IF
      DO 30 I=1,M
   30 MM(I) = K
C
      RETURN
      E N D
cc      SUBROUTINE  FILTR2( Y,XF,VF,Q,M,K,N,NOBS,MJ,NDIM,LLOUT,
      SUBROUTINE  FILTR2K( Y,XF,VF,Q,M,K,N,NOBS,MJ,NDIM,LLOUT,
     *                    OUTMIN,OUTMAX,VFS,VPS,XFS,XPS,FF,SIG2 )
C
C  ...  Kalman Filter (for Time Varying AR model)  ...
C
C     Inputs:
C        Y:      time series
C        XF:     Initial state vector
C        VF:     Initial covariance matrix
C        Q:      K*K matrix, system noise covariance
C        M:      Dimension of the state vector
C        K:      Dimension of the system noise
C        N:      Data length
C        NOBS:   Local Stationary Span ( N/NOBS =< NDIM )
C        MJ:     Adjustable dimension of XF, VF
C        NDIM:   Adjustable dimension of VPS, VFS etc.
C        OUTMIN: Lower limit for detecting outliers
C        OUTMAX: Upper limit for detecting outliers
C     Outputs:
C        VFS:    Covariance matrices of the filter
C        VPS:    Covariance matrices of the predictor
C        XFS:    Mean vectors of the filter
C        XPS:    Mean vectors of the predictor
C        FF:     Log likelihood
C        SIG2:   Estimated observational noise variance
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION  Y(N), LLOUT(NDIM)
cc      DIMENSION  XF(MJ), VF(MJ,MJ), XP(40), VP(40,40)
      DIMENSION  XF(MJ), VF(MJ,MJ), XP(MJ), VP(MJ,MJ)
      DIMENSION  XFS(MJ,NDIM),    XPS(MJ,NDIM)
      DIMENSION  VFS(MJ,MJ,NDIM), VPS(MJ,MJ,NDIM)
cc      DIMENSION  VH(40), GAIN(40)
      DIMENSION  VH(MJ), GAIN(MJ)
      DATA   PI  /3.1415926535D0/
C
      MK = M*K
      SIG2 = 0.0D0
      SDET = 0.0D0
      NSUM = 0
      NS = M/NOBS+1
      NE = N/NOBS
C
      DO 500  II=NS,NE
C
C  ...  ONE STEP AHEAD PREDICTION  ...
C
      IF( K.EQ.1 )  THEN
         DO 10 J=1,M
         XP(J) = XF(J)
         DO 10 I=1,M
   10    VP(I,J) = VF(I,J)
      END IF
      IF( K.EQ.2 )  THEN
         DO 20 J=1,M
         J2 = 2*J
         XP(J2-1) = 2*XF(J2-1) + XF(J2)
         XP(J2)   =  -XF(J2-1)
         DO 20 I=1,M
         I2 = 2*I
         VP(I2-1,J2-1) = 4*VF(I2-1,J2-1) + 2*VF(I2,J2-1)
     *                 + 2*VF(I2-1,J2)   + VF(I2,J2)
         VP(I2-1,J2)   =-2*VF(I2-1,J2-1) - VF(I2,J2-1)
         VP(I2,J2-1)   =-2*VF(I2-1,J2-1) - VF(I2-1,J2)
   20    VP(I2,J2)     =   VF(I2-1,J2-1)
      END IF
C
      DO 30 I=1,M
   30 VP(K*(I-1)+1,K*(I-1)+1) = VP(K*(I-1)+1,K*(I-1)+1) + Q
      IF( LLOUT(II).EQ.1 )  THEN
        DO 35 I=1,MK
   35   VP(I,I) = 1.0D3
      END IF
C
C  ...  SAVE MEAN AND COVARIANCE  ...
C
      IF( NDIM.GT.1 )  THEN
      DO 40  I=1,MK
      XPS(I,II) = XP(I)
      DO 40  J=1,MK
   40 VPS(I,J,II) = VP(I,J)
      END IF
C
C
C  ...  FILTERING  ...
C
      DO 400  JJ=1,NOBS
      I1 = NOBS*(II-1) + JJ
      IF( I1.LE.M )  GO TO 400
      IF( Y(I1).GT.OUTMIN.AND.Y(I1).LT.OUTMAX.AND. I1.LE.N ) THEN
C
      DO 210  I=1,MK
      SUM = 0.0D0
      DO 200  J=1,M
  200 SUM = SUM + VP(I,K*(J-1)+1)*Y(I1-J)
  210 VH(I) = SUM
C
      PERR = Y(I1)
      PVAR = 1.0D0
      DO 220  I=1,M
      PERR = PERR - Y(I1-I)*XP(K*(I-1)+1)
  220 PVAR = PVAR + Y(I1-I)*VH(K*(I-1)+1)
C
      DO 250  I=1,MK
  250 GAIN(I) = VH(I)/PVAR
C
      DO 290  I=1,MK
  290 XF(I) = XP(I) + GAIN(I)*PERR
C
      DO 310  I=1,MK
      DO 310  J=1,MK
  310 VF(I,J) = VP(I,J) - GAIN(I)*VH(J)
C
      IF( I1.GT.MJ )  THEN
        SIG2 = SIG2 + PERR**2/PVAR
        SDET = SDET + DLOG(PVAR)
        NSUM = NSUM + 1
      END IF
C
C  ...  MISSING OBSERVATION  ...
C
      ELSE
      DO 350  I=1,MK
      XF(I) = XP(I)
      DO 350  J=1,MK
  350 VF(I,J) = VP(I,J)
      END IF
      IF( JJ.NE.NOBS )  THEN
        DO 360 J=1,MK
        XP(J) = XF(J)
        DO 360 I=1,MK
  360   VP(I,J) = VF(I,J)
      END IF
  400 CONTINUE
C
C  ...  SAVE MEAN AND COVARIANCE  ...
C
      IF( NDIM.GT.1 )  THEN
      DO 370  I=1,MK
      XFS(I,II) = XF(I)
      DO 370  J=1,MK
  370 VFS(I,J,II) = VF(I,J)
      END IF
C
  500 CONTINUE
      SIG2 = SIG2/NSUM
      FF = -0.5D0*(NSUM*(DLOG(PI*2*SIG2) + 1) + SDET)
C
      DO 510 II=1,NS-1
      DO 510 J=1,MK
      XPS(J,II) = XPS(J,NS)
      XFS(J,II) = XFS(J,NS)
      DO 510 I=1,MK
      VPS(I,J,II) = VPS(I,J,NS)
  510 VFS(I,J,II) = VFS(I,J,NS)
C
      RETURN
      E N D
      SUBROUTINE  ISTCAR( M,K,MJ,XF,VF )
C
C  ...  Initial State (for Time-Varying AR Model)  ...
C
C     Inputs:
C        M:     AR order
C        K:     Trend order
C        MJ:    Adjustable dimension of F
C     Outputs:
C         XF:   State vector, X(0|0)
C         VF:   State covarance matrix, V(0|0)
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION  XF(MJ), VF(MJ,MJ)
C
      MK = M*K
      DO 10  J=1,MK
      XF(J) = 0.0D0
      DO 10  I=1,MK
   10 VF(I,J) = 0.0D0
C
      DO 20  I=1,MK
   20 VF(I,I) = 1.0D2
C
      RETURN
      E N D
C
C
cc      SUBROUTINE  SMOTH1( A,M,MMAX,NC,NS,N,NE,MJ,VFS,VPS,VSS,
cc     *                    XFS,XPS,XSS )
      SUBROUTINE  SMOTH1( A,M,MMAX,NC,NS,N,NE,MJ,VFS,VPS,XFS,XPS,XSS )
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION  A(MMAX,NC), M(NC)
cc      DIMENSION  XS(40), VS(40,40), VP(40,40)
      DIMENSION  XS(MJ), VS(MJ,MJ), VP(MJ,MJ)
      DIMENSION  XFS(MJ,N), XPS(MJ,N), XSS(MJ,N)
      DIMENSION  VFS(MJ,MJ,N), VPS(MJ,MJ,N), VSS(MJ,MJ,N)
cc      DIMENSION  WRK(40,40), SGAIN(40,40)
cc      DIMENSION  I0(10)
      DIMENSION  WRK(MJ,MJ), SGAIN(MJ,MJ)
      DIMENSION  I0(NC)
C
      I0(1) = 0
      DO 10 I=2,NC
   10 I0(I) = I0(I-1) + M(I-1)
      MM = I0(NC) + M(NC)
C
C  ...  SMOOTHING  ...
C
      NSS = MIN0( N,NE )
      DO 20  I=1,MM
      XS(I)      = XFS(I,NSS)
      XSS(I,NSS) = XFS(I,NSS)
      DO 20  J=1,MM
      VS(I,J)      = VFS(I,J,NSS)
   20 VSS(I,J,NSS) = VFS(I,J,NSS)
      DO 30 II=NSS+1,NE
      DO 30 I=1,MM
      XSS(I,II) = XFS(I,II)
      DO 30 J=1,MM
   30 VSS(I,J,II) = VFS(I,J,II)
C
      DO 500  II=NSS-1,NS,-1
C
      NZERO = 0
      DO 100 I=1,MM
  100 IF( VFS(I,I,II).GT.1.0D-12 )  NZERO = NZERO + 1
C
      IF( NZERO.EQ.0 )  THEN
         DO 110  I=1,MM
         XS(I)     = XFS(I,II)
         XSS(I,II) = XFS(I,II)
         DO 110  J=1,MM
         VS(I,J)     = VFS(I,J,II)
  110    VSS(I,J,II) = VFS(I,J,II)
C
      ELSE
      DO 410  I=1,MM
      DO 410  J=1,MM
  410 VP(I,J) = VPS(I,J,II+1)
C
cc      CALL  GINVRS( VP,VDET,MM,40 )
      CALL  GINVRS( VP,VDET,MM )
C
      DO 420  I=1,MM
      DO 420  L=1,NC
      WRK(I,I0(L)+M(L)) = VFS(I,I0(L)+1,II)*A(M(L),L)
      DO 420  J=1,M(L)-1
  420 WRK(I,I0(L)+J) = VFS(I,I0(L)+1,II)*A(J,L) + VFS(I,I0(L)+J+1,II)
C
      DO 440  I=1,MM
      DO 440  J=1,MM
      SUM = 0.0D0
      DO 430 IJ=1,MM
  430 SUM = SUM + WRK(I,IJ)*VP(IJ,J)
  440 SGAIN(I,J) = SUM
C
      DO 450  I=1,MM
      XS(I) = XFS(I,II)
      DO 450  J=1,MM
      WRK(I,J) = 0.0D0
  450 VS (I,J) = VFS(I,J,II)
C
      DO 460  J=1,MM
      DO 460  I=1,MM
  460 XS(I) = XS(I) + SGAIN(I,J)*(XSS(J,II+1) - XPS(J,II+1))
C
      DO 470  J=1,MM
      DO 470 IJ=1,MM
      DO 470  I=1,MM
  470 WRK(I,J) = WRK(I,J) + SGAIN(I,IJ)*(VSS(IJ,J,II+1)-VPS(IJ,J,II+1))
C
      DO 480  J=1,MM
      DO 480 IJ=1,MM
      DO 480  I=1,MM
  480 VS(I,J) = VS(I,J) + WRK(I,IJ)*SGAIN(J,IJ)
      DO 485 I=1,MM
  485 IF( VS(I,I).LT.0.0D0 )  VS(I,I) = 0.0D0
C
      DO 490  I=1,MM
      XSS(I,II) = XS(I)
      DO 490  J=1,MM
  490 VSS(I,J,II) = VS(I,J)
C     WRITE(6,*) (XS(I),I=1,MM)
C     DO 88 I=1,MM
C  88 WRITE(6,89) (VS(I,J),J=1,MM)
C  89 FORMAT( 1H ,11D12.4 )
      END IF
C
  500 CONTINUE
C
      RETURN
      E N D
cc      SUBROUTINE  PTCAR( XSS,WRK,N,NOBS,MJ,M,K,TAU2,SIG2,FF,AIC )
      SUBROUTINE  PTCAR( XSS,N,NOBS,MJ,M,K,TAR,WRK )
C
C  ...  Plot Original Time Varying PARCOR  ...
C
C     Inputs:
C        XSS:    Smoothed state
C        N:      Data length
C        NOBS:   Local Stationary Span
C        MJ:     Adjustable dimension of F and G
C        M,K:    AR and Trend Orders
C        TAU2, SIG2: Variances of the State Space Model
C        FF:     Log-Likelihood of the model
C        AIC:    AIC of the Model
C     MODIFIED  2/15/93
C
      IMPLICIT REAL*8(A-H,O-Z)
cc      CHARACTER  VNAME*8
cc      DIMENSION  XSS(MJ,N), WRK(MJ,N), AR(10), PAR(10)
      DIMENSION  XSS(MJ,N), WRK(M,N/nobs), AR(M), PAR(M), TAR(M,N/NOBS)
cc      DIMENSION  DATA(400), VNAME(20), VALUE(20)
C
cc      VNAME(1) = 'M     = '
cc      VNAME(2) = 'K     = '
cc      VNAME(3) = 'NOBS  = '
cc      VNAME(4) = 'TAU2  = '
cc      VNAME(5) = 'SIG2  = '
cc      VNAME(6) = 'FF    = '
cc      VNAME(7) = 'AIC   = '
cc      VALUE(1) = M
cc      VALUE(2) = K
cc      VALUE(3) = NOBS
cc      VALUE(4) = TAU2
cc      VALUE(5) = SIG2
cc      VALUE(6) = FF
cc      VALUE(7) = AIC
cc      CALL  PLOTS
C     call  plots( 1,0,0,1,0 )
C     call  form( 1 )
C     call  factor( 10.0 )
cc      CALL  HEADER( 'PROGRAM 12.1: TVCAR MODEL ',33,7,
cc     *               VNAME,VALUE )
cc      WX = 12.0
cc      WY =  3.0
cc      DX =500.0
      FN = N
      NN = N/NOBS
cc      IY = 10
cc      YMIN1 = -1.0D0
cc      YMAX1 =  1.0D0
cc      DY = 1.0D0
      DO 40 II=1,NN
      DO 10 I=1,M
   10 AR(I) = XSS(K*(I-1)+1,II)
      CALL  PARCOR( AR,M,PAR )
      DO 20 I=1,M
      IF( PAR(I).GT. 0.95D0 )  PAR(I) = 0.95D0
   20 IF( PAR(I).LT.-0.95D0 )  PAR(I) =-0.95D0
      CALL  ARCOEF( PAR,M,AR )
      DO 30 I=1,M
      XSS(K*(I-1)+1,II) = AR(I)
c--------------------------------
      TAR(I,II) = AR(I)
c--------------------------------
   30 WRK(I,II) = PAR(I)
   40 CONTINUE
C
C  ...   PLOT  PARCOR  ...
C
cc      DO 200 II=1,M
cc      IF( II.EQ.9 )  CALL  PLOTI
C     if( ii.eq.9 )  call  plot( 0.0,0.0,777 )
cc      IF( MOD(II,8).EQ.1 )  CALL  PLOT( 2.0,17.0-SNGL(WY)-1.0,-3 )
cc      IF( II.EQ.5 )  CALL  PLOT( SNGL(WX)+2.0,3*SNGL(WY+1.0),-3 )
cc      IF( MOD(II,4).NE.1 )  CALL  PLOT( 0.0,-SNGL(WY+1.0),-3 )
C     CALL  SYMBOL( 0.25,SNGL(WY)+0.25,0.25,'TREND',0.0,5 )
cc      CALL  AXISXY( 0.0D0,0.0D0,WX,WY,0.0D0,FN,YMIN1,YMAX1,DX,DY,
cc     *              0.2D0,1,IY,2)
C
cc      DO 100 I=1,NN
cc  100 DATA(I) = WRK(II,I)
cc      CALL  NEWPEN( 1 )
cc      CALL  PLOTY( DATA,NN,YMIN1,YMAX1,WX,WY,1,1 )
cc  200 CONTINUE
C
cc      CALL  PLOTE
C     call  plot( 0.0,0.0,999 )
      RETURN
      E N D
