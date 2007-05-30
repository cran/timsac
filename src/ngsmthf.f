      SUBROUTINE  NGSMTHF( Y,N,NOISEV,TAU2,BV,NOISEW,SIG2,BW,INITD,
     *                      TREND,SS,FF,NS,NFE,NPE,K )
C
      INCLUDE 'timsac_f.h'
C
C     PROGRAM 14.1  NGSMTH
C
C  ...  NON-GAUSSIAN SMOOTHING  ...
C
C     INPUTS:
C        NOISEV:   TYPE OF SYSTEM NOISE DENSITY (0,1,2,3)
C        NOISEW:   TYPE OF OBSER. NOISE DENSITY (0,1,2,3,4)
C        TAU2:     VARIANCE OR DISPERSION OF SYSTEM NOISE
C        SIG2:     VARIANCE OR DISPERSION OF OBSERVATION NOISE
C        BV:       SHAPE PARAMETER OF SYSTEM NOISE (FOR NOISEV = 2)
C        BW:       SHAPE PARAMETER OF OBSER. NOISE (FOR NOISEW = 2)
C     PARAMETERS:
C        MJ:       ADJUSTABLE DIMENSION OF Y, PS, FS, SS (MJ >= N)
C        K:        NUMBER OF INTERVALS + 1
C     @NFILTER.F:   5/14/85, 6/11/87, 5/19/88, 2/14/91
C     MODIFIED  2/16/93
C
cc      !DEC$ ATTRIBUTES DLLEXPORT :: NGSMTHF
C
cc      PARAMETER( MJ=500, K=201 )
      IMPLICIT REAL*8(A-H,O-Z)
cc      CHARACTER  TITLE*72
cc      DIMENSION  Y(MJ), P(K), F(K), S(K), T(K), Q(-K:K)
cc      DIMENSION  YM(10), ST(7), TREND(MJ,7), LOC(MJ)
cc      REAL*4    PS(K,MJ), FS(K,MJ), SS(K,MJ)
      DIMENSION  Y(N), F(K), Q(-K:K)
      DIMENSION  YM(2), ST(7), TREND(NPE,7), LOC(NPE)
      REAL*4     SS(K,NPE)
cc      COMMON   /COMAAA/  Y
cc      COMMON   /C91214/  XMIN, XMAX, FIGMIN, FIGMAX, YM
cc      COMMON   /C91215/  NOISEV, NOISEW, INITD, ITF, ITH
cc      COMMON   /C91216/  B, OUTMIN, OUTMAX
cc      COMMON   /C91218/  SIG2, TAU2, BW, BV
cc      COMMON   /C91219/  N, KK, NS, NFE, NPE
cc      COMMON   /CMDATA/  TITLE
cc      COMMON    / DDD /  FF , AIC , SD
C     NAMELIST  /PARAM/  XMIN, XMAX, TAU20, SIG20, BC, IOPT, B,
C    *                   OUTMIN, OUTMAX, NS, NFE, NPE, FIGMIN, FIGMAX,
C    *                   NOISEV, NOISEW, ITF, ITH
C
cc      EXTERNAL  GAUSS
cc      EXTERNAL  PEARSN
cc      EXTERNAL  TWOEXP
cc      EXTERNAL  USERV
C
C  ...  READ DATA  ...
C 
cc      CALL  READTS( 1,Y,N )
cc      CALL  MOMENT( Y,N,YM(1),YM(2) )
      CALL  MOMENTK( Y,N,YM(1),YM(2) )
C
      LOC(1) = 0
cc      KK = K
cc      CALL  DEFALT
      CALL  DEFALT( Y,N,XMIN,XMAX,OUTMIN,OUTMAX )
cc      READ( 5,* )  NOISEV, TAU2, BV
cc      READ( 5,* )  NOISEW, SIG2, BW
C
      DX = (XMAX-XMIN)/(K-1)
cc      CALL  IDIST( F,K,YM(1),YM(2),XMIN,DX )
      CALL  IDIST( F,K,YM(1),YM(2),XMIN,DX,INITD )
      CALL  NORMLZ( F,K,DX,SSUM )
cc      IF( NOISEV.EQ.0 )  CALL  TRANS( USERV ,K,DX,TAU2,BV,Q )
cc      IF( NOISEV.EQ.1 )  CALL  TRANS( GAUSS ,K,DX,TAU2,BV,Q )
cc      IF( NOISEV.EQ.2 )  CALL  TRANS( PEARSN,K,DX,TAU2,BV,Q )
cc      IF( NOISEV.EQ.3 )  CALL  TRANS( TWOEXP,K,DX,TAU2,BV,Q )
      IF( NOISEV.EQ.0 )  CALL  TRANS1( K,DX,TAU2,BV,Q )
      IF( NOISEV.EQ.1 )  CALL  TRANS2( K,DX,TAU2,BV,Q )
      IF( NOISEV.EQ.2 )  CALL  TRANS3( K,DX,TAU2,BV,Q )
      IF( NOISEV.EQ.3 )  CALL  TRANS4( K,DX,TAU2,BV,Q )
cc      CALL  NGSMTH( Y,P,F,S,T,N,K,DX,XMIN,Q,FF,PS,FS,SS,LOC )
      CALL  NGSMTH( NOISEW,SIG2,BW,Y,F,N,K,DX,XMIN,Q,FF,SS,LOC,
     *              OUTMIN,OUTMAX,NS,NFE,NPE )
C
      DO 30 I=1,NPE
      DO 10 J=1,K
cc   10 S(J) = SS(J,I)
cc      CALL  PINTVL( S,K,XMIN,DX,ST )
   10 F(J) = DBLE(SS(J,I))
      CALL  PINTVL( F,K,XMIN,DX,ST )
      DO 20 J=1,7
   20 TREND(I,J) = ST(J) + DX*LOC(I)
   30 CONTINUE
C
cc      CALL  PRNGSM( NOISEV,NOISEW,TAU2,SIG2,FF,TREND,MJ,N )
C
cc      CALL  PTNGSM( N,NPE,NOISEV,NOISEW,TAU2,SIG2,FF,B,Y,TREND,MJ,
cc     *              FIGMIN,FIGMAX,SS,K,LOC )
      CALL  POST3D( SS,LOC,K,NPE )
cc      STOP
      RETURN
      E N D
cc      SUBROUTINE  DEFALT
      SUBROUTINE  DEFALT( Y,N,XMIN,XMAX,OUTMIN,OUTMAX )
C
C  ...  THIS SUBROUTINE SETS DEFAULT VALUES OF PARAMETERS  ...
C
      IMPLICIT REAL*8(A-H,O-Z)
cc      DIMENSION  Y(500)
      DIMENSION  Y(N)
cc      COMMON   /COMAAA/  Y
cc      COMMON   /C91214/  XMIN, XMAX, FIGMIN, FIGMAX, YM(10)
cc      COMMON   /C91215/  NOISEV, NOISEW, INITD, ITF, ITH
cc      COMMON   /C91216/  B, OUTMIN, OUTMAX
cc      COMMON   /C91219/  N, KK, NS, NFE, NPE
cc      NOISEV = 2
cc      NOISEW = 1
cc      INITD  = 1
cc      ITF    = 1
cc      ITH    = 1
cc      B = 1.00D0
      OUTMIN = -1.0D30
      OUTMAX =  1.0D30
cc      CALL  MAXMIN( Y,N,XMIN,XMAX,DY )
      CALL  MAXMINK( Y,N,XMIN,XMAX,DY )
cc      NS = 1
cc      NFE = N
cc      NPE = N
cc      FIGMIN = XMIN
cc      FIGMAX = XMAX
C
      RETURN
      E N D
cc      SUBROUTINE  NGSMTH( Y,P,F,S,T,N,K,DX,XMIN,Q,FF,PS,FS,SS,LOC )
      SUBROUTINE  NGSMTH( NOISEW,SIG2,BW,Y,F,N,K,DX,XMIN,Q,FF,SS,LOC,
     *                     OUTMIN,OUTMAX,NS,NFE,NPE )
C
C  ...  NON-GAUSSIAN SMOOTHER  ...
C
C     INPUTS:
C       Y(I):   TIME SERIES
C       P(I):   INITIAL DENSITY
C       N:      DATA LENGTH
C       K:      NUMBER OF INTERVALS IN STEP FUNCTION APPROXIMATION
C       DX:     WIDTH OF INTERVAL
C       XMIN:   MINIMUM OF THE INTERVAL
C       Q:      SYSTEM NOISE DENSITY
C     OUTPUTS:
C       FF:     LOG-LIKELIHOOD
C       LOC(I): LOCATION OF THE CENTER OF THE INTERVAL AT STEP I
C       SS:     SMOOTHED DENSITY
C
      IMPLICIT REAL*8(A-H,O-Z)
cc      DIMENSION  P(K), F(K), S(K), T(K), Y(N), Q(-K:K), LOC(N)
cc      REAL*4     FS(K,N), PS(K,N), SS(K,N)
      DIMENSION  P(K), F(K), S(K), T(K), Y(N), Q(-K:K), LOC(NPE)
      REAL*4     PS(K,NPE), SS(K,NPE)
cc      COMMON   /C91215/  NOISEV, NOISEW, INITD, ITF, ITH
cc      COMMON   /C91216/  B, OUTMIN, OUTMAX
cc      COMMON   /C91219/  NN, KK, NS, NFE, NPE
C
      FF = 0.0D0
C
      DO 200 II=NS,NPE
C
C  ...  CONVOLUTION (SYSTEM NOISE)  ...
C
      CALL  CONVOL( Q,F,K,P )
      CALL  NORMLZ( P,K,DX,PSUM )
C
C  ...  BAYES FORMULA  ...
C
      IF( Y(II).LE.OUTMIN .OR. Y(II).GE.OUTMAX .OR. II.GT.NFE ) THEN
      DO 110 I=1,K
  110 F(I) = P(I)
      ELSE
cc      CALL  BAYES( P,K,XMIN,DX,Y(II),F,LOC(II) )
      CALL  BAYES( NOISEW,SIG2,BW,P,K,XMIN,DX,Y(II),F,LOC(II) )
      CALL  NORMLZ( F,K,DX,FINT )
C
C  ...  LIKELIHOOD COMPUTATION  ...
C
      FF = FF + DLOG( FINT )
cc      IF( MOD(II,10).EQ.0)  WRITE(6,*) II,FF
      END IF
C
C  ...  SAVE FOR SMOOTHING  ...
C
      DO 130 I=1,K
cc      PS(I,II) = P(I)
cc  130 FS(I,II) = F(I)
      PS(I,II) = SNGL(P(I))
  130 SS(I,II) = SNGL(F(I))
C
C  ...  SHIFT ORIGIN  ...
C
cc      CALL  SHIFT( F,K,T,II,N,LOC )
      CALL  SSHIFT( F,K,T,II,N,LOC )
C
  200 CONTINUE
C
C  ...  SMOOTHING  ...
C
cc      DO 190 J=NFE,NPE
cc      DO 190 I=1,K
cc  190 SS(I,J) = FS(I,J)
      DO 195 I=1,K
cc  195 S(I) = FS(I,NFE)
  195 S(I) = DBLE(SS(I,NFE))
C
      DO 300 II=NFE-1,NS,-1
cc      IF(MOD(II,10).EQ.0)  WRITE(6,*) II
      DO 210 I=1,K
      T(I) = 0.0D0
      P(I) = 0.0D0
cc  210 F(I) = FS(I,II)
  210 F(I) = DBLE(SS(I,II))
      DO 220 I=1,K
      J = I - (LOC(II+1)-LOC(II))
cc      IF( J.GE.1.AND.J.LE.K )  P(I) = PS(J,II+1)
      IF( J.GE.1.AND.J.LE.K )  P(I) = DBLE(PS(J,II+1))
  220 IF( J.GE.1.AND.J.LE.K )  T(I) = S(J)
      DO 230 I=1,K
  230 S(I) = T(I)
C
cc      CALL  SCONVL( Q,S,P,F,K,T )
      CALL  SCONVLK( Q,S,P,F,K,T )
      CALL  NORMLZ( T,K,DX,TSUM )
C
      DO 240 I=1,K
      S(I) = T(I)
cc  240 SS(I,II) = S(I)
  240 SS(I,II) = SNGL(S(I))
  300 CONTINUE
C
      RETURN
      E N D
      SUBROUTINE  PINTVL( P,K,XMIN,DX,Y )
      IMPLICIT REAL*8(A-H,O-Z)
cc      DIMENSION  P(K), Y(7), PROB(7), P1(401)
      DIMENSION  P(K), Y(7), PROB(7), P1(K)
      DATA  PROB /0.0013D0, 0.0227D0, 0.1587D0, 0.5000D0, 0.8413D0,
     *            0.9773D0,0.9987D0/
C
      P1(1) = 0.0
      DO 10 I=2,K
   10 P1(I) = P1(I-1) + (P(I-1) + P(I))*DX/2.0
C
      DO 30 J=1,7
      PP = PROB(J)
      DO 20 I=2,K
      IF(P1(I-1).LE.PP .AND. P1(I).GT.PP)  GO TO 30
   20 CONTINUE
   30 Y(J) = XMIN + (I-2)*DX + DX*(PP - P1(I-1))/(P1(I) - P1(I-1))
C
      RETURN
      E N D
      SUBROUTINE  CONVOL( Q,S,K,P )
      IMPLICIT  REAL*8(A-H,O-Z)
      DIMENSION  S(K), P(K), Q(-K:K)
C
      DO 20 I=1,K
      J1 = 1-I
      J2 = K-I
      SUM = 0.0
      DO 10 J=J1,J2
   10 SUM = SUM + S(I+J)*Q(J)
   20 P(I) = SUM
C
      RETURN
      E N D
cc      SUBROUTINE  SCONVL( Q,P,R,S,K,T )
      SUBROUTINE  SCONVLK( Q,P,R,S,K,T )
      IMPLICIT  REAL*8(A-H,O-Z)
      DIMENSION  S(K), P(K), R(K), T(K), Q(-K:K)
C
      DO 20 I=1,K
      J1 = 1-I
      J2 = K-I
      SUM = 0.0D0
      DO 10 J=J1,J2
   10 IF(P(I+J).GT.0.0D0)  SUM = SUM + P(I+J)/R(I+J)*Q(J)
   20 T(I) = S(I)*SUM
      RETURN
      E N D
      SUBROUTINE  TRANS1(K,DX,TAU2,BV,Q )
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION  Q(-K:K), PARAM(3)
cc      EXTERNAL   FUNCT
C
      PARAM(1) = 0.0D0
      PARAM(2) = TAU2
      PARAM(3) = BV
C
      DO 20 I=1-K,K-1
      X0 = -DX*I - DX/2
      SUM = (USERV(X0,PARAM) + USERV(X0+DX,PARAM))/2
      DO 10 J=1,49
      X = X0 + (DX*J)/50
   10 SUM = SUM + USERV( X,PARAM )
   20 Q(I) = SUM*DX/50
C
      RETURN
      E N D
      SUBROUTINE  TRANS2( K,DX,TAU2,BV,Q )
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION  Q(-K:K), PARAM(3)
cc      EXTERNAL   FUNCT
C
      PARAM(1) = 0.0D0
      PARAM(2) = TAU2
      PARAM(3) = BV
C
      DO 20 I=1-K,K-1
      X0 = -DX*I - DX/2
      SUM = (GAUSS(X0,PARAM) + GAUSS(X0+DX,PARAM))/2
      DO 10 J=1,49
      X = X0 + (DX*J)/50
   10 SUM = SUM + GAUSS( X,PARAM )
   20 Q(I) = SUM*DX/50
C
      RETURN
      E N D
      SUBROUTINE  TRANS3( K,DX,TAU2,BV,Q )
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION  Q(-K:K), PARAM(3)
cc      EXTERNAL   FUNCT
C
      PARAM(1) = 0.0D0
      PARAM(2) = TAU2
      PARAM(3) = BV
C
      DO 20 I=1-K,K-1
      X0 = -DX*I - DX/2
      SUM = (PEARSN(X0,PARAM) + PEARSN(X0+DX,PARAM))/2
      DO 10 J=1,49
      X = X0 + (DX*J)/50
   10 SUM = SUM + PEARSN( X,PARAM )
   20 Q(I) = SUM*DX/50
C
      RETURN
      E N D
      SUBROUTINE  TRANS4( K,DX,TAU2,BV,Q )
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION  Q(-K:K), PARAM(3)
cc      EXTERNAL   FUNCT
C
      PARAM(1) = 0.0D0
      PARAM(2) = TAU2
      PARAM(3) = BV
C
      DO 20 I=1-K,K-1
      X0 = -DX*I - DX/2
      SUM = (TWOEXP(X0,PARAM) + TWOEXP(X0+DX,PARAM))/2
      DO 10 J=1,49
      X = X0 + (DX*J)/50
   10 SUM = SUM + TWOEXP( X,PARAM )
   20 Q(I) = SUM*DX/50
C
      RETURN
      E N D

      SUBROUTINE  NORMLZ( P,K,DX,SUM )
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION  P(K)
C
      SUM = 0.0D0
      DO 10 I=1,K
   10 SUM = SUM + P(I)
      SUM = SUM*DX
      DO 30 I=1,K
   30 P(I) = P(I)/SUM
C
      RETURN
      E N D
cc      SUBROUTINE  IDIST( P,K,P1,P2,XMIN,DX )
      SUBROUTINE  IDIST( P,K,P1,P2,XMIN,DX,INITD )
      IMPLICIT REAL*8( A-H,O-Z )
      DIMENSION  P(K), PARAM(3)
cc      COMMON   /C91215/  NOISEV, NOISEW, INITD, ITF, ITH
C
      PARAM(1) = P1
      PARAM(2) = P2
C
      DO 10 I=1,K
      X = XMIN + DX*(I-1)
      IF( INITD.EQ.0 )  P(I) = USERI( X,PARAM )
      IF( INITD.EQ.1 )  P(I) = GAUSS( X,PARAM )
      IF( INITD.EQ.2 )  P(I) = UNIF ( X,PARAM )
   10 CONTINUE
      RETURN
      END
cc      SUBROUTINE  BAYES( P,K,XMIN,DX,Y,F,LSHIFT )
      SUBROUTINE  BAYES( NOISEW,SIG2,BW,P,K,XMIN,DX,Y,F,LSHIFT )
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION  P(K), F(K), PARAM(3)
cc      COMMON   /C91215/  NOISEV, NOISEW, INITD, ITF, ITH
cc      COMMON   /C91218/  SIG2, TAU2, BW, BV
cc      EXTERNAL  USERW
cc      EXTERNAL  GAUSS
cc      EXTERNAL  PEARSN
cc      EXTERNAL  TWOEXP
cc      EXTERNAL  DBLEXP
C
      PARAM(2) = SIG2
      PARAM(3) = BW
C
      DO 10 I=1,K
      PARAM(1) = XMIN + DX*(I-1+LSHIFT)
      IF( NOISEW.EQ.0 )  F(I) = P(I)*USERW ( Y,PARAM )
      IF( NOISEW.EQ.1 )  F(I) = P(I)*GAUSS ( Y,PARAM )
      IF( NOISEW.EQ.2 )  F(I) = P(I)*PEARSN( Y,PARAM )
      IF( NOISEW.EQ.3 )  F(I) = P(I)*TWOEXP( Y,PARAM )
      IF( NOISEW.EQ.4 )  F(I) = P(I)*DBLEXP( Y,PARAM )
   10 CONTINUE
      RETURN
      E N D
      DOUBLE PRECISION FUNCTION  USERW( Y,PARAM )
      IMPLICIT  REAL*8(A-H,O-Z)
      DIMENSION  PARAM(3)
      DATA  C1  /2.506628275D0/
C
      YMEAN = 0.0D0
      VAR   = PARAM(1)
      USERW = DEXP( -(Y-YMEAN)**2/(2*VAR) )/(C1*DSQRT( VAR ))
      RETURN
      E N D
      DOUBLE PRECISION FUNCTION  USERV( X,PARAM )
      IMPLICIT  REAL*8(A-H,O-Z)
      DIMENSION  PARAM(3)
      DATA  C1  /2.506628275D0/
C
      USERV = DEXP( -(X-PARAM(1))**2/(2*PARAM(2)) )
     *                /(C1*DSQRT( PARAM(2) ))
      RETURN
      E N D
      DOUBLE PRECISION FUNCTION  TWOEXP( X,PARAM )
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION  PARAM(2)
C
      TWOEXP = DEXP( -DABS(X-PARAM(1))*PARAM(2) )*PARAM(2)/2.0D0
      RETURN
C
      E N D
cc      SUBROUTINE  SHIFT( F,K,T,II,N,LOC )
      SUBROUTINE  SSHIFT( F,K,T,II,N,LOC )
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION  F(K), T(K), LOC(*)
C
C  ...  FIND THE POSTERIOR MODE AND SHIFT ORIGIN  ...
C
      PMAX = 0.0D0
c-----------
      IMAX = 1
c-----------
      DO 10  I=1,K
      IF( F(I).LE.PMAX )  GO TO 10
         PMAX = F(I)
         IMAX = I
   10 CONTINUE
      IF(II.LT.N)  LOC(II+1) = LOC(II) + IMAX - (K+1)/2
      DO 20 I=1,K
      J = I + IMAX - (K+1)/2
      T(I) = 0.0D0
   20 IF(J.GE.1.AND.J.LE.K )  T(I) = F(J)
      DO 30 I=1,K
   30 F(I) = T(I)
C
      RETURN
      E N D
      SUBROUTINE  POST3D( F,LOC,K,N )
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*4     F(K,N)
cc      DIMENSION  YH(2000), FF(801), LOC(N)
      REAL*4     FF(-K:2*K)
      DIMENSION  LOC(N)
cc      COMMON   /C91214/  X0, X1, F0, F1, YM(10)
cc      ANGLE = 70.0
C
cc      WX = 8.0
cc      WY = 2.0
cc      WZ =12.0
cc      CALL  DELX( X0,X1,DX )
cc      Z0 =  0.0
cc      Z1 =  N
cc      KK = 1
      NDIF = 1
cc      IF( N.GE.100 )  NDIF = 2
cc      IF( N.GE.200 )  NDIF = 4
cc      IF( N.GE.300 )  NDIF = 6
cc      IF( N.GE.500 )  NDIF = N/50
cc      DZ = NDIF*10
      NN = N/NDIF
      N0 = NDIF/2 + 1
cc      DO 10 I=1,K
cc   10 FF(I) = F(I,1)
cc      CALL  MAXMIN( FF,K,Y0,Y1,DY )
C
cc      CALL  PLOT3A( K,NN,WX,WY,WZ,ANGLE,X0,X1,DX,Y0,Y1,DY,Z0,Z1,DZ,
cc     *              KN,KK,YH )
c
      DO 100 J=N0,N,NDIF
cc      DO 80 I=1,K
      DO 80 I=-K,2*K
   80 FF(I) = 0.0D0
      II = LOC(J)
      I1 = MAX0( 1,II )
      I2 = MIN0( K,K+II )
      DO 90 I=I1,I2
   90 FF(I+II) = F(I,J)
cc      CALL  PLOT3B( FF,K,NN,WX,WY,WZ,ANGLE,Y0,Y1,KK,YH,ZS,ZC )
      DO 95 I=1,K
   95 F(I,J) = FF(I)
  100 CONTINUE
cc      CALL  PLOT( -SNGL(ZC),-SNGL(ZS),-3 )
C
      RETURN
      E N D
      DOUBLE PRECISION FUNCTION  GAUSS( X,PARAM )
C
C  ...  Gaussian (normal) distribution  ...
C
C     Inputs:
C        X:
C        PARAM(1):  mean
C        PARAM(2):  variance
C     Output:
C        GAUSS:     density at X
C
      IMPLICIT  REAL*8(A-H,O-Z)
      DIMENSION  PARAM(2)
      DATA  C1  /2.506628275D0/
C
      GAUSS = DEXP( -(X-PARAM(1))**2/(2*PARAM(2)) )/(C1*DSQRT(PARAM(2)))
      RETURN
      E N D
      DOUBLE PRECISION FUNCTION  PEARSN( X,PARAM )
C
C  ...  Pearson family of  distributions  ...
C
C     Inputs:
C        X:
C        PARAM(1):  location parameter, mu
C        PARAM(2):  dispersion parameter, tau square
C        PARAM(3):  shape parameter
C     Output:
C        PEARSN:    density at X
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION  PARAM(3)
      DATA  PI/3.1415926535D0/
C
      PEARSN = DGAMMA(PARAM(3))/DGAMMA(PARAM(3)-0.5D0)
     *                  /DSQRT(PI)*PARAM(2)**(PARAM(3)-0.5D0)
     *                  /((X-PARAM(1))**2 + PARAM(2))**PARAM(3)
      RETURN
C
      END
      DOUBLE PRECISION  FUNCTION  DGAMMA( X )
C
C  ...  Gamma function  ...
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION  A(0:10)
      DATA  A /0.999999999999269D0, 0.42278433696202D0,
     *  0.41184025179616D0, 0.08157821878492D0, 0.0742379070629D0,
     * -0.00021090746731D0, 0.01097369584174D0,-0.00246674798054D0,
     *  0.00153976810472D0,-0.00034423420456D0, 0.00006771057117D0/
C
      DGAM = 1.0D0
      Y = X
      IF( X.GT.3.0D0 )  THEN
   10   Y = Y-1
        DGAM = DGAM*Y
        IF( Y.GT.3.0D0 )  GO TO 10
      END IF
      IF( X.LE.2.0D0 )  THEN
   20   DGAM = DGAM/Y
        Y = Y+1
        IF( Y.LE.2.0D0 )  GO TO 20
      END IF
C
      Z   = 1.0D0
      SUM = 0.0D0
      DO 30 I=0,10
      SUM = SUM + A(I)*Z
   30 Z = Z*(Y-2)
      DGAMMA = DGAM*SUM
      RETURN
      E N D
      DOUBLE PRECISION FUNCTION  DBLEXP( X,PARAM )
C
C  ...  double exponential distribution  f(x) = exp(x - exp(x))  ...
C
C     Inputs:
C        X:
C     Output:
C        DBLEXP:    density at X
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DBLEXP = DEXP( X-DEXP(X) )
      RETURN
C
      E N D
      DOUBLE PRECISION FUNCTION  USERI( X,PARAM )
C
C  ...  User supplied density function  ...
C       (The following is an example of two-sided exponential dist.)
C
C     Inputs:
C        X:
C        PARAM(1):  mean
C        PARAM(2):  lambda
C     Output:
C        USERI:     density at X
C
      IMPLICIT  REAL*8(A-H,O-Z)
      DIMENSION  PARAM(2)
C
      SIGMA = DSQRT( PARAM(2) )
         USERI = SIGMA*DEXP( -SIGMA*DABS(X-PARAM(1)) )/2
      RETURN
      E N D
      DOUBLE PRECISION FUNCTION  UNIF( X,PARAM )
C
C  ...  uniform distribution  f(x) = 1  ...
C
C     Output:
C        UNIF:    density at X (=1)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      UNIF = 1.0D0
      RETURN
      E N D
cc      SUBROUTINE  MAXMIN( X,N,XMIN0,XMAX0,DXL )
      SUBROUTINE  MAXMINK( X,N,XMIN0,XMAX0,DXL )
C
C  ...  This subroutine determines the minimum, the maximum and
C       the step width  ...
C
C     Inputs:
C        X(I):    data
C        N:       data length
C     Outputs:
C        XMIN0:   the minimum bound for the figure
C        XMAX0:   the maximum bound for the figure
C        DXL:     step width
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(N)
C
      XMIN = 1.0D30
      XMAX =-1.0D30
C
      DO 10 I=1,N
      IF( X(I) .LT. XMIN )   XMIN = X(I)
   10 IF( X(I) .GT. XMAX )  XMAX = X(I)
      DX = XMAX-XMIN
      IF( DLOG10(DX) .GE. 0.0D0 )  DXL = INT( DLOG10(DX) )
      IF( DLOG10(DX) .LT. 0.0D0 )  DXL = INT( DLOG10(DX) )-1.0
      DXL = 10.0**DXL
      IF( DX/DXL.GT.6.0D0 )  DXL = DXL*2.0
      DIF = INT( DX/DXL )
      XMIN0 = INT( XMIN/DXL )*DXL
      XMAX0 = XMIN0 + DIF*DXL
      IF( XMIN0 .GT. XMIN )  XMIN0 = XMIN0 - DXL
   30 IF( XMAX0 .GE. XMAX )  GO TO 40
        XMAX0 = XMAX0 + DXL
        GO TO 30
   40 CONTINUE
C
      RETURN
      E N D
