      SUBROUTINE FFTCORF(LD,LAGH1,N,N2P,ISW,X1,Y1,XA,X,Y,
     1 CNA1,CN1,CN2,AMEAN)
C
      INCLUDE 'timsac.h'
C
C     PROGRAM 5.1.3   AUTO AND/OR CROSS CORRELATIONS VIA FFT.
C-----------------------------------------------------------------------
C     ** THIS PROGRAM IS AN ADAPTED VERSION OF THE ALGOL PROCEDURE
C         //CROSS CORRELATE// PREPARED BY DR GORDON SANDE AT PRINCETON
C         UNIVERSITY IN 1966-7.
C     ** DESIGNED BY H. AKAIKE, THE INSTITUTE OF STATISTICAL MATHEMATICS
C     ** PROGRAMMED BY E. ARAHATA, THE INSTITUTE OF STATISTICAL MATHEMAT
C         TOKYO
C     ** DATE OF THE LATEST REVISION: FEB. 16, 1978
C     ** THIS PROGRAM WAS ORIGINALLY PUBLISHED IN
C         "DAINAMIKKU SISTEMU NO TOKEI-TEKI KAISEKI TO SEIGYO (STATISTICA
C         ANALYSIS AND CONTROL OF DYNAMIC SYSTEMS)" BY H. AKAIKE AND
C         T. NAKAGAWA, SAIENSU-SHA, TOKYO, 1972 (IN JAPANESE)
C-----------------------------------------------------------------------
C     THIS PROGRAM COMPUTES AUTO AND/OR CROSS
C     COVARIANCES AND CORRELATIONS VIA FFT.
C     IT REQUIRES FOLLOWING INPUTS:
C     ISW: ISW=1...AUTO CORRELATION OF X (ONE-CHANNEL)
C           ISW=2...AUTO CORRELATIONS OF X AND Y (TWO-CHANNEL)
C           ISW=4...AUTO,CROSS CORRELATIONS OF X AND Y (TWO-CHANNEL)
C     LD: LENGTH OF DATA
C     LAGH: MAXIMUM LAG
C     DFORM: INPUT FORMAT SPECIFICATION STATEMENT IN ONE CARD,
C     FOR EXAMPLE
C     (8F10.4)
C     (X(I); I=1,LD): DATA OF CHANNEL X
C     (Y(I); I=1,LD): DATA OF CHANNEL Y (FOR ISW=2 OR 4 ONLY)
C
cxx      IMPLICIT REAL*8(A-H,O-Y)
cxx      IMPLICIT COMPLEX*16(Z)
c      REAL*4 DFORM
c      DIMENSION X(2048),Y(2048),Z(2048),ZS(1025)
c      DIMENSION CN1(501),CN2(501)
c      DIMENSION DFORM(20)
c      REAL*4 XS,YS
c      DIMENSION XS(2048),YS(2048)
cxx      DIMENSION X1(LD),Y1(LD),XA(N,2),X(N),Y(N)
cxx      DIMENSION Z(N),ZS(N/2+1)
cxx      DIMENSION CNA1(LAGH1,2),CN1(LAGH1),CN2(LAGH1),AMEAN(2)
      INTEGER LD, LAGH1, N, N2P, ISW
      DOUBLE PRECISION X1(LD), Y1(LD), XA(N,2), X(N), Y(N),
     1                 CNA1(LAGH1,2), CN1(LAGH1), CN2(LAGH1), AMEAN(2)
c local
      DOUBLE PRECISION CST0, CST1, CST2, AN, ALD, ALD1, RF, SF, RG,
     1                 SG, XI, XNI, YI, YNI, X0, XMEAN, YMEAN, CX0, Y0
      COMPLEX(kind(0d0)) :: Z(N), ZS(N/2+1), ZI, ZNI
C     INPUT / OUTPUT DATA FILE OPEN
c	CHARACTER(100) DFNAM
c	DFNAM='fftcor.out'
c	CALL SETWND
c	CALL FLOPN3(DFNAM,NFL)
c	IF (NFL.EQ.0) GO TO 999
C     INITIAL CONDITION INPUT AND PUNCH OUT
c      READ(5,1) ISW,LD,LAGH
c      WRITE(6,50)
c      WRITE(6,51)
c      WRITE(6,52) ISW,LD,LAGH
c      WRITE(7,1) LD,LAGH
c      LAGH1=LAGH+1
      LAGH=LAGH1-1
      ND=LD+LAGH1
C     N2P, N: DEFINITION
c      I0=1
c   10 IR1=2**I0
c      IF(IR1-ND) 11,12,12
c   11 I0=I0+1
c      GO TO 10
c   12 N2P=I0
c      N=2**N2P
      NP1=N+1
      NP2=N+2
      M=N/2
      M1=M+1
      CST0=0.0D-00
      CST1=1.0D-00
      CST2=0.25D-00
      AN=N
      ALD=LD
      ALD1=CST1/(AN*ALD)
C     INPUT FORMAT SPECIFICATION
c      READ(5,4) (DFORM(I),I=1,20)
c    4 FORMAT(20A4)
C     ORIGINAL DATA INPUT AND OUTPUT
      DO 20 I=1,N
      X(I)=CST0
      Y(I)=CST0
   20 CONTINUE
c      READ(5,DFORM) (XS(I),I=1,LD)
      DO 1200 I=1,LD
c      X(I)=DBLE(XS(I))
      X(I)=X1(I)
 1200 CONTINUE
      IF(ISW.EQ.1) GO TO 200
c      READ(5,DFORM) (YS(I),I=1,LD)
      DO 1201 I=1,LD
c      Y(I)=DBLE(YS(I))
      Y(I)=Y1(I)
 1201 CONTINUE
c  200 WRITE(6,53)
  200 CONTINUE
c      IF(ISW.NE.1) GO TO 201
c      WRITE(6,54)
c      CALL PRCOL1(X,1,LD,0)
c      GO TO 202
c  201 WRITE(6,55)
c      CALL PRCOL2(X,Y,1,LD,0)
C     MEAN DELETION
cxx  202 CALL DMEADL(X,LD,XMEAN)
      CALL DMEADL(X,LD,XMEAN)
      IF(ISW.EQ.1) GO TO 203
      CALL DMEADL(Y,LD,YMEAN)
C     DOUBLE PRECISION COMPLEX REPRESENTATION
  203 DO 31 I=1,N
cxx   31 Z(I)=DCMPLX(X(I),Y(I))
cxx      Z(I)=DCMPLX(X(I),Y(I))
      Z(I)=CMPLX(X(I),Y(I),KIND=8)
   31 CONTINUE
C     FOURIER TRANSFORM OF Z
      ISG=-1
      CALL MIXRAD(Z,N,N2P,ISG)
      IF(ISW.NE.1) GO TO 204
C     RAW SPECTRUM COMPUTATION
      DO 32 I=2,M
cxx      X(I)=DREAL(Z(I))**2+DIMAG(Z(I))**2
      X(I)=REAL(Z(I))**2+AIMAG(Z(I))**2
      NI=NP2-I
cxx   32 X(NI)=X(I)
      X(NI)=X(I)
   32 CONTINUE
cxx      X(1)=DREAL(Z(1))**2
cxx      X(M1)=DREAL(Z(M1))**2
      X(1)=REAL(Z(1))**2
      X(M1)=REAL(Z(M1))**2
      GO TO 205
C     DECOMPOSITION AND RAW SPECTRUM COMPUTATION
  204 DO 125 I=2,M
      NI=NP2-I
      ZI=Z(I)
      ZNI=Z(NI)
cxx      RF=DREAL(ZI)
cxx      SF=DIMAG(ZI)
cxx      RG=DREAL(ZNI)
cxx      SG=DIMAG(ZNI)
      RF=REAL(ZI)
      SF=AIMAG(ZI)
      RG=REAL(ZNI)
      SG=AIMAG(ZNI)
      XI=RF+RG
      XNI=SF-SG
cxx      Z(I)=DCMPLX(XI,XNI)
      Z(I)=CMPLX(XI,XNI,KIND=8)
      X(I)=CST2*(XI**2+XNI**2)
      X(NI)=X(I)
      YI=SF+SG
      YNI=RF-RG
cxx      Z(NI)=DCMPLX(YI,YNI)
      Z(NI)=CMPLX(YI,YNI,KIND=8)
      Y(I)=CST2*(YI**2+YNI**2)
      Y(NI)=Y(I)
  125 CONTINUE
cxx      X(1)=DREAL(Z(1))**2
cxx      Y(1)=DIMAG(Z(1))**2
cxx      X(M1)=DREAL(Z(M1))**2
cxx      Y(M1)=DIMAG(Z(M1))**2
      X(1)=REAL(Z(1))**2
      Y(1)=AIMAG(Z(1))**2
      X(M1)=REAL(Z(M1))**2
      Y(M1)=AIMAG(Z(M1))**2
      IF(ISW.NE.4) GO TO 205
C     RAW CROSS SPECTRUM COMPUTATION
      DO 126 I=2,M
      NI=NP2-I
cxx  126 ZS(I)=CST2*Z(I)*Z(NI)
      ZS(I)=CST2*Z(I)*Z(NI)
  126 CONTINUE
cxx      ZS(1)=DREAL(Z(1))*DIMAG(Z(1))
cxx      ZS(M1)=DREAL(Z(M1))*DIMAG(Z(M1))
      ZS(1)=REAL(Z(1))*AIMAG(Z(1))
      ZS(M1)=REAL(Z(M1))*AIMAG(Z(M1))
C     AUTO COVARIANCE COMPUTATION
  205 DO 33 I=1,N
cxx   33 Z(I)=DCMPLX(X(I),Y(I))
cxx      Z(I)=DCMPLX(X(I),Y(I))
      Z(I)=CMPLX(X(I),Y(I),KIND=8)
   33 CONTINUE
C     FOURIER TRANSFORM
cxx  215 CALL MIXRAD(Z,N,N2P,ISG)
      CALL MIXRAD(Z,N,N2P,ISG)
      II=1
      DO 34 I=1,LAGH1
cxx      X(I)=DREAL(Z(I))*ALD1
      X(I)=REAL(Z(I))*ALD1
cxx   34 XA(I,II)=X(I)
      XA(I,II)=X(I)
   34 CONTINUE
      X0=X(1)
c      AMEAN=XMEAN
      AMEAN(II)=XMEAN
C     NORMALIZATION
   36 CX0=X(1)
c      CALL CORNOM(X,CN1,LAGH1,CX0,CX0)
      CALL CORNOM(X,CNA1(1,II),LAGH1,CX0,CX0)
C     AUTO COVARIANCE PRINT OUT
c      WRITE(6,162) II,II,AMEAN
c      WRITE(6,163)
c      CALL PRCOL2(X,CN1,1,LAGH1,1)
C     AUTO COVARIANCE PUNCH OUT
c      WRITE(7,1) II,II
c      WRITE(7,2) (X(I),I=1,LAGH1)
      IF(ISW.EQ.1) GO TO 300
      IF(II.EQ.2) GO TO 216
      II=2
      DO 35 I=1,LAGH1
cxx      X(I)=DIMAG(Z(I))*ALD1
      X(I)=AIMAG(Z(I))*ALD1
cxx   35 XA(I,II)=X(I)
      XA(I,II)=X(I)
   35 CONTINUE
      Y0=X(1)
c      AMEAN=YMEAN
      AMEAN(II)=YMEAN
      GO TO 36
  216 IF(ISW.NE.4) GO TO 300
C     CROSS COVARIANCE COMPUTATION
      DO 127 I=2,M
      NI=NP2-I
      Z(I)=ZS(I)
cxx  127 Z(NI)=DCONJG(ZS(I))
cxx      Z(NI)=DCONJG(ZS(I))
      Z(NI)=CONJG(ZS(I))
  127 CONTINUE
      Z(1)=ZS(1)
      Z(M1)=ZS(M1)
C     FOURIER TRANSFORM
      CALL MIXRAD(Z,N,N2P,ISG)
      DO 41 I=1,LAGH
      I1=I+1
      J1=NP1-I
cxx      X(I1)=DREAL(Z(I1))*ALD1
      X(I1)=REAL(Z(I1))*ALD1
cxx   41 Y(I1)=DREAL(Z(J1))*ALD1
cxx      Y(I1)=DREAL(Z(J1))*ALD1
      Y(I1)=REAL(Z(J1))*ALD1
   41 CONTINUE
cxx      X(1)=DREAL(Z(1))*ALD1
      X(1)=REAL(Z(1))*ALD1
      Y(1)=X(1)
C     NORMALIZATION
      CALL CORNOM(X,CN1,LAGH1,X0,Y0)
      CALL CORNOM(Y,CN2,LAGH1,X0,Y0)
C     CROSS COVARIANCE PRINT OUT
c      JJ=1
c      WRITE(6,165) II,JJ
c      WRITE(6,166)
c      CALL PRCOL4(X,CN1,Y,CN2,1,LAGH1,1)
C     CROSS COVARIANCE PUNCH OUT
c      WRITE(7,1) II,JJ
c      WRITE(7,2) (X(I),I=1,LAGH1)
c      WRITE(7,1) JJ,II
c      WRITE(7,2) (Y(I),I=1,LAGH1)
  300 CONTINUE
c	CALL FLCLS3(NFL)
c  999 CONTINUE
      RETURN
c    1 FORMAT(10I5)
c    2 FORMAT(4D20.10)
c   50 FORMAT(1H ,71HPROGRAM 5.1.3   AUTO AND/OR CROSS COVARIANCES AND CO
c     ARRELATIONS VIA FFT.)
c   51 FORMAT(1H ,17HINITIAL CONDITION)
c   52 FORMAT(1H ,4HISW=,I5,5X,3HLD=,I5,5X,5HLAGH=,I5)
c   53 FORMAT(1H ,13HORIGIANL DATA)
c   54 FORMAT(1H ,4X,1HI,12X,4HX(I))
c   55 FORMAT(1H ,4X,1HI,12X,4HX(I),10X,4HY(I))
c  162 FORMAT(//1H ,14HAUTOCOVARIANCE,5X,6HCIJ(L),5X,2HI=,I5,5X,2HJ=,I5,5
c     AX,5HMEAN=,D15.5)
c  163 FORMAT(1H ,4X,1HL,5X,6HCIJ(L),8X,10HNORMALIZED)
c  165 FORMAT(//1H ,16HCROSS COVARIANCE,5X,6HCIJ(L),5X,2HI=,I5,5X,2HJ=,I5
c     A)
c  166 FORMAT(1H ,4X,1HL,5X,6HCIJ(L),8X,10HNORMALIZED,4X,6HCJI(L),8X,10HN
c     AORMALIZED)
      END SUBROUTINE
