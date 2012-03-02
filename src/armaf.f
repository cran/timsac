      SUBROUTINE ARMAF(M,L,A,B,SIG2,N,K,KMAX,NF,G,COV,PAR,SP,
     * ROOTA,ROOTB,IER,JER)
C
      INCLUDE 'timsac_f.h'
C
C     PROGRAM 6.1  ARMA
C
C  ...  This program analyses time series via ARMA modeling  ...
C
C     Inputs:
C        IDEV:  Input devie (=5: CONSOLE)
C        M:     AR order
C        L:     MA order
C        SIG2:  Innovation variance
C        A(I):  AR coefficients
C        B(I):  MA coefficients
C        K:     Maximum lag of autocovariance function (K<MJ+1)
C     Parameter:
C        MJ:    Adjustable dimension of A, B, etc..
C        NF:    Number of frequencies in evaluating spectrum
C     Programmed by Y.I and G.K.
C
C     Outputs:
C         IER:  =1 : MATRIX WITH ZERO ROW IN DECOMPOSE
C               =2 : SINGULAR MATRIX IN DECOMPOSE.ZERO DIDIVIDE IN SOLVE
C               =3 : CONVERGENCE IN IMPRUV.MATRIX IS NEARLY SINGULAR
C         JER:  =1 : NON-CONVERGENCE AT POLYRT
C
cc      !DEC$ ATTRIBUTES DLLEXPORT::ARMAF
C
cc      PARAMETER( NF=200,MJ=100,IDEV=1 )
      IMPLICIT REAL*8( A-H,O-Z )
cc      DIMENSION  A(MJ), B(MJ), PAR(MJ), ROOTA(MJ,2), ROOTB(MJ,2)
cc      DIMENSION  G(0:MJ), COV(0:MJ), SP(0:NF)
cc      DIMENSION  WRK1(0:MJ), WRK2(0:MJ), WRK3(MJ)
cc      DATA  PAR/MJ*0.0D0/
      DIMENSION  A(M), B(L), PAR(K), ROOTA(M,2), ROOTB(L,2)
      DIMENSION  G(0:KMAX), COV(0:K), SP(0:NF)
      DIMENSION  WRK1(0:K), WRK2(0:K), WRK3(K,K)
C
cc      WRITE( 6,* )  'K = ?'
cc      READ( 5,* )   K
cc      OPEN( IDEV,FILE='arma.dat' )
cc      READ(IDEV,*)  M, L
cc      READ(IDEV,*)  SIG2
cc      READ(IDEV,*)  (A(I),I=1,M)
cc      READ(IDEV,*)  (B(I),I=1,L)
cc      CLOSE( IDEV )
C
cc      KMAX = MAX(M,L,K)
      CALL  IMPULS( M,L,A,B,K,G )
cc      CALL  ARMCOV( M,L,A,B,SIG2,K,COV )
      CALL  ARMCOV( M,L,A,B,SIG2,K,COV,KMAX,IER )
c------------
      PAR = 0
c------------
      CALL  PARCOR( A,M,PAR )
      CALL  ARCOEF( PAR,M,A )
cc      IF( L.GT.0 )  CALL  ARYULE( COV,1000,K,WRK1,WRK2,PAR,WRK3,MAR )
      IF( L.GT.0 )  CALL  ARYULE( COV,N,K,WRK1,WRK2,PAR,WRK3,MAR )
      CALL  ARMASP( A,M,B,L,SIG2,NF,SP )
cc      CALL  CHROOT( A,M,ROOTA,MJ )
cc      CALL  CHROOT( B,L,ROOTB,MJ )
      CALL  CHROOT( A,M,ROOTA,M,JER )
      CALL  CHROOT( B,L,ROOTB,L,JER )
cc      CALL  PRARMA( M,L,A,B,G,K,COV,K,PAR,SP,NF,ROOTA,ROOTB,MJ )
C      CALL  PTARMA( G,K,COV,K,PAR,SP,NF,ROOTA,M,ROOTB,L,MJ )
cc      STOP
      RETURN
      E N D
cc      SUBROUTINE  CHROOT( A,M,ROOT,MJ )
      SUBROUTINE  CHROOT( A,M,ROOT,MJ,IER )
C
C  ...  Characteristic roots of the AR or MA operator  ...
C
C     Inputs:
C        A:     Vector of AR or MA coefficients
C        M:     Order of the model
C        MJ:    Adjustable dimension of ROOT
C     Output:
C        ROOT:  Characteristic roots (real part,imaginary part)
C
      IMPLICIT  REAL*8(A-H,O-Z)
      DIMENSION  A(M), ROOT(MJ,2)
cc      DIMENSION  C(50), CW(50)
      DIMENSION  C(M+1)
C
      IF( M.EQ.0 )  RETURN
      DO 10  I=1,M
   10 C(I) = -A(M-I+1)
      C(M+1) = 1.0D0
      MMAX = M
C
C  ... Characteristic roots of operator 1-A(1)*S- ... -A(M)*S**M=0  ...
C
cc      CALL  POLYRT( C,CW,MMAX,ROOT(1,1),ROOT(1,2),IER )
      CALL  POLYRT( C,MMAX,ROOT(1,1),ROOT(1,2),IER )
cc      IF( IER.NE.0 )   WRITE(6,600)
C
      RETURN
  600 FORMAT( 1H0,'*****  NON-CONVERGENCE AT POLYRT  *****' )
      E N D
cc      SUBROUTINE  POLYRT( A,B,M,ROOTR,ROOTI,IER )
      SUBROUTINE  POLYRT( A,M,ROOTR,ROOTI,IER )
C
C  ...  This subroutine finds the roots of the equation
C            A(1) + A(2)*S + ... + A(M)*(S**M) = 0
C       by Newton-Raphson method
C
C     Inputs:
C        A:     Coefficients of the equation
C        B:     Working area
C        M:     Degree of the polynomial
C     Outputs:
C        ROOTR:   Real parts of the roots
C        ROOTI:   Imaginary parts of the roots
C        IER:     Error code to indicate non-convergence
C
      IMPLICIT  REAL*8(A-H,O-Z)
cc      DIMENSION A(1), B(1), ROOTR(1), ROOTI(1)
      DIMENSION A(1), B(M+1), ROOTR(1), ROOTI(1)
C
      IFIT = 0
      ISW = 0
      JSW = 0
      K = M
      IER = 0
      KX = K
      KXX = K+1
      KJI = KXX
      K2 = 1
c------------
      XPR = 0.0D0
      YPR = 0.0D0
c------------
      DO 10  I=1,KXX
      A(I) = A(I)/A(K+1)
      J = KXX - I+1
   10 B(J) = A(I)
   20 XO = 0.5D-02
      YO = 0.1D-01
      IN = 0
   30 X = XO
C
      XO = -10.D0*YO
      YO = -10.D0*X
      X = XO
      Y = YO
      ICT = 0
      IN = IN+1
      GO TO 50
C
   40 IFIT = 1
      XPR = X
      YPR = Y
C
   50 UX = 0.D0
      UY = 0.D0
      V  = 0.D0
      YT = 0.D0
      XT = 1.D0
      U = B(K+1)
      IF( U .EQ. 0.D0 )  GO TO 140
      DO 60  I=1,K
      L = K-I+1
      XT2 = X*XT - Y*YT
      YT2 = X*YT+Y*XT
      U = U + B(L)*XT2
      V = V + B(L)*YT2
      FI = I
      UX = UX + FI*XT*B(L)
      UY = UY - FI*YT*B(L)
      XT = XT2
      YT = YT2
   60 CONTINUE
      SUM = UX**2 + UY**2
      IF( SUM .EQ. 0.D0 )  GO TO 100
      DX = (V*UY - U*UX)/SUM
      X = X + DX
      DY = -(U*UY + V*UX)/SUM
      Y = Y + DY
      IF( DABS(DY)+DABS(DX) .LT. 1.0D-10 )  GO TO 80
C
      ICT = ICT+1
      IF( ICT .LT. 500 )  GO TO 50
      ISW = 1
      IF( IN .GE. 5 )  GO TO 70
      IF( IFIT .NE. 0 )  GO TO 80
   65 ISW = 0
      IFIT = 0
      GO TO 30
C
   70 IF( IFIT .EQ. 0 )  GO TO 300
      JSW = 1
   80 DO 90  L=1,KXX
      J = KJI-L+1
      TEM = A(J)
      A(J) = B(L)
   90 B(L) = TEM
      ITEM = K
      K = KX
      KX = ITEM
  100 IF( IFIT .EQ. 0 )  GO TO 40
      IF( JSW .EQ. 1 )   GO TO 110
cc      IF( ISW-1 )  120,65,120
      IF( ISW-1 .LT. 0 )  GO TO 120
      IF( ISW-1 .EQ. 0 )  GO TO 65
      IF( ISW-1 .GT. 0 )  GO TO 120
  110 X = XPR
      Y = YPR
      ISW = 0
      JSW = 0
  120 IFIT = 0
      IF( X .EQ. 0.D0 )  GO TO 130
      IF( DABS( Y/X ) .LT. 1.0D-08 )  GO TO 150
  130 ALPH = 2*X
      SUM = X*X + Y*Y
      K = K-2
      GO TO 160
C
  140 X = 0.D0
      KX = KX-1
      KXX = KXX-1
  150 Y = 0.D0
      SUM = 0.D0
      ALPH = X
      K = K-1
  160 B(2) = B(2) + ALPH*B(1)
  170 DO 180  L=2,K
  180 B(L+1) = B(L+1) + ALPH*B(L) - SUM*B(L-1)
  190 ROOTI(K2) = Y
      ROOTR(K2) = X
      K2 = K2+1
      IF( SUM .EQ. 0.D0 )  GO TO 200
      Y = -Y
      SUM = 0.D0
      GO TO 190
  200 IF (K .LE. 0 )  RETURN
      GO TO 20
C
  300 IER = 1
      RETURN
      E N D
      SUBROUTINE ARYULE( C,N,MAXM,SIG2,AIC,PARCOR,A,MAR )
C
C  ...  Yule-Walker method  ...
C
C     Inputs:
C        C(I):    Autocovariance function
C        N:       Data length
C        MAXM:    Highest AR order
C     Outputs:
C        SIG2(I): Innovation variance
C        AIC(I):  AIC
C        PARCOR(I):  PARCOR
C        AMIN:     AR coefficients of the best model
C        MAR:      Selected order of the model
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION  C(0:MAXM), SIG2(0:MAXM), AIC(0:MAXM)
      DIMENSION  PARCOR(MAXM), A(MAXM,MAXM)
      CONST = N*(DLOG(2*3.1415926535D0) + 1)
C
      SIG2(0) = C(0)
      AIC(0) = CONST + N*DLOG(SIG2(0)) + 2
      AICMIN = AIC(0)
      MAR = 0
C
      DO 50 M=1,MAXM
      SUM  = C(M)
      DO 10 I=1,M-1
   10 SUM  = SUM - A(I,M-1)*C(M-I)
      A(M,M) = SUM /SIG2(M-1)
      DO 20 J=1,M-1
   20 A(J,M) = A(J,M-1)-A(M,M)*A(M-J,M-1)
      SIG2(M) = SIG2(M-1)*(1.0D0-A(M,M)**2)
      AIC(M) = CONST + N*DLOG(SIG2(M)) + 2*(M+1)
      PARCOR(M) = A(M,M)
      IF( AIC(M).LT.AICMIN )  THEN
         AICMIN = AIC(M)
         MAR = M
      END IF
   50 CONTINUE
      RETURN
      E N D
cc      SUBROUTINE ARMCOV( M,L,A,B,SIG2,K,COV )
      SUBROUTINE ARMCOV( M,L,A,B,SIG2,K,COV,KMAX,IER )
C
C ...  Autocovariance Function of ARMA model  ...
C
C     Inputs:
C        M:     AR order
C        L:     MA order
C        A(I):  AR coefficient
C        B(I):  MA coefficient
C        SIG2:  innovation variance
C        K:     Required maximum lag of autocovariance
C     Output:
C        COV(I):  Autocovariance function
C     Y.I.
      IMPLICIT REAL*8( A-H,O-Z )
cc      DIMENSION  A(*), B(*), COV(0:K), G(0:100), X(30,30)
cc      DIMENSION  Z(100), UL(30,30), IPS(100)
cxxxxxxxxx      DIMENSION  A(*), B(*), COV(0:K), G(0:KMAX), X(M+1,M+1)
      DIMENSION  A(1), B(1), COV(0:K), G(0:KMAX), X(M+1,M+1)
      DIMENSION  Z(M+1), UL(M+1,M+1), IPS(M+1)
C
cc      KMAX = MAX(M,L,K)
      CALL  IMPULS( M,L,A,B,KMAX,G )
C
      DO 10 I=1,M+1
      DO 10 J=1,M+1
   10 X(I,J) = 0.0D0
      DO 20 I=1,M+1
   20 X(I,I) = 1.0D0
      DO 30 I=1,M
      DO 30 J=2,M-I+2
   30 X(I,J) = X(I,J) - A(I+J-2)
      DO 40 I=2,M+1
      DO 40 J=1,I-1
   40 X(I,J) = X(I,J) - A(I-J)
C
cc      CALL  DECOM( M+1,X,30,UL,IPS )
      CALL  DECOM( M+1,X,UL,IPS,IER )
C
      SUM = 1.0D0
      DO 50 J=1,L
   50 SUM = SUM - B(J)*G(J)
      Z(1)= SIG2*SUM
      DO 70 I=2,M+1
      SUM = 0.0D0
      DO 60 J=I-1,L
   60 SUM = SUM - B(J)*G(J-I+1)
   70 Z(I) = SIG2*SUM
C
cc      CALL  SOLVE( M+1,UL,30,Z,COV,IPS)
      CALL  SOLVE1( M+1,UL,Z,COV,IPS)
C
      DO 100 J=M+1,K
      SUM = 0.0D0
      DO 80 I=1,M
   80 SUM = SUM + A(I)*COV(J-I)
      DO 90 I=J,L
   90 SUM = SUM - B(I)*G(I-J)*SIG2
  100 COV(J) = SUM
C
      RETURN
      E N D
C
cc      SUBROUTINE DECOM( N,A,MJ,UL,IPS )
      SUBROUTINE DECOM( N,A,UL,IPS,IER )
C
C  ...  UL decomposition:  A = L*U  ...
C
C     Inputs:
C        N:      Dimension of the matrix A
C        A(I,J): N*N positive definite matrix
C        MJ:     Adjustable dimension of A and UL
C     Outputs:
C        UL(I,J):  L-I and U
C        IPS:    Index vector
C     Y.I.
C        IER:    Error code
C
      IMPLICIT REAL*8(A-H,O-Z )
cc      DIMENSION A(MJ,*),UL(MJ,*),SCALES(100),IPS(100)
      DIMENSION A(N,N),UL(N,N),SCALES(N),IPS(N)
C
      IER = 0
C
      DO 20 I=1,N
      IPS(I) = I
      RNORM = 0.0D0
      DO 10 J=1,N
      UL(I,J) = A(I,J)
      IF(RNORM.LT.ABS(UL(I,J)) ) RNORM = ABS( UL(I,J) )
   10 CONTINUE
      IF( RNORM .NE. 0.0D0 )  THEN
          SCALES(I) = 1/RNORM
      ELSE
          SCALES(I) = 0.0D0
cc          CALL  SING(0)         
          IER = 1
      END IF
   20 CONTINUE
C
cc-------------
      INDEX = 0
cc-------------
      DO 60 K=1,N-1
      BIG = 0.0D0
      DO 30 I=K,N
      SIZE = ABS( UL(IPS(I),K) )*SCALES( IPS(I) )
      IF( BIG.LT.SIZE ) THEN
          BIG = SIZE
          INDEX = I
      END IF
   30 CONTINUE
      IF( BIG.EQ. 0.0D0 )  THEN
cc          CALL  SING(1)
          IER = 2
          GO TO 60
      END IF
      IF( INDEX.NE.K ) THEN
      J = IPS(K)
      IPS(K) = IPS(INDEX)
      IPS(INDEX) = J
      END IF
C
      PIVOT = UL(IPS(K),K)
      DO 50 I=K+1,N
      TM = UL( IPS(I),K)/PIVOT
      UL( IPS(I),K) = TM
      IF( TM.NE. 0.0D0 )  THEN
      DO 40 J = K+1,N
   40 UL( IPS(I),J ) = UL( IPS(I),J)-TM*UL( IPS(K),J)
C     WRITE(6,*) (UL(IPS(I),J),J=1,N)
      END IF
   50 CONTINUE
   60 CONTINUE
C
cc      IF( UL(IPS(N),N) .EQ. 0.0D0 )   CALL  SING(2)
      IF( UL(IPS(N),N) .EQ. 0.0D0 )   IER = 3
      RETURN
      E N D
C
      SUBROUTINE IMPULS( M,L,A,B,K,G )
C
C ...  Impulse Response Function  ...
C
C     Inputs:
C        M:     AR order
C        L:     MA order
C        A(I):  AR coefficient
C        B(I):  MA coefficient
C        K:     Required maximum lag of impulse respose
C     Output:
C        G(I):  Impulse response function
C     Y.I.
      IMPLICIT REAL*8( A-H,O-Z )
cxxxxx      DIMENSION A(*), B(*), G(0:K)
      DIMENSION A(1), B(1), G(0:K)
C
      G(0) = 1.0
      DO  20 I=1,K
      SUM = 0.0D0
      IF(I.LE.L) SUM = -B(I)
      DO  10 J=1,I
   10 IF(J.LE.M) SUM = SUM + A(J)*G(I-J)
   20 G(I) = SUM
C
      RETURN
      E N D
C
cc      SUBROUTINE  SOLVE( N,UL,MJ,B,X,IPS )
      SUBROUTINE  SOLVE1( N,UL,B,X,IPS )
C
C  ...  Solve Ax=b using UL obtained by DECOM  ...
C
C     Inputs:
C        N:     Dimension of UL and B
C        UL:    LU decomposition of A
C        MJ:    Adjustable dimension of A
C        B:
C        IPS:   index vector
C     Output:
C        X:     Solution
C     Y.I.
      IMPLICIT REAL*8( A-H,O-Z )
cc      DIMENSION UL(MJ,*),B(*),X(*),IPS(100)
      DIMENSION UL(N,N),B(N),X(N),IPS(N)
C
      DO 20 I=1,N
      SUM = 0.0D0
      DO 10 J=1,I-1
   10 SUM = SUM + UL(IPS(I),J)*X(J)
   20 X(I) = B(IPS(I)) - SUM
C
      DO 40 I=N,1,-1
      SUM = 0.0D0
      DO 30 J=I+1,N
   30 SUM = SUM + UL(IPS(I),J)*X(J)
   40 X(I) = ( X(I)-SUM )/UL(IPS(I),I)
      RETURN
  600 FORMAT(1H ,'N=',I10,/,(5X,'IPS=',I10 ) )
      E N D
