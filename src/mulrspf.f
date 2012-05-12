      SUBROUTINE MULRSPF(H,L,IP,K,SD,A,B,Y,CH)
C
      INCLUDE 'timsac_f.h'
C
C     PROGRAM 5.4.2   MULTIPLE RATIONAL SPECTRUM
C-----------------------------------------------------------------------
C      SUBROUTINE XYCTRX(X,Y,Z,MM,NN)
C-----------------------------------------------------------------------
C     ** DESIGNED BY H. AKAIKE, THE INSTITUTE OF STATISTICAL MATHEMATICS
C     ** PROGRAMMED BY E. ARAHATA, THE INSTITUTE OF STATISTICAL MATHEMAT
C	  TOKYO
C     ** DATE OF THE LATEST REVISION: MARCH 25, 1977
C     ** THIS PROGRAM WAS ORIGINALLY PUBLISHED IN
C	 "DAINAMIKKU SISTEMU NO TOKEI-TEKI KAISEKI TO SEIGYO (STATISTICA
C	 ANALYSIS AND CONTROL OF DYNAMIC SYSTEMS)" BY H. AKAIKE AND
C	 T. NAKAGAWA, SAIENSU-SHA, TOKYO, 1972 (IN JAPANESE)
C-----------------------------------------------------------------------
C     THIS PROGRAM COMPUTES RATIONAL SPECTRUM FOR IP-DIMENSIONAL
C     AR-MA PROCESS
C     X(N)=A(1)X(N-1)+...+A(L)X(N-L)+E(N)+B(1)E(N-1)+...+B(K)E(N-K),
C     WHERE E(N) IS A WHITE NOISE WITH ZERO MEAN VECTOR AND COVARIANCE
C     MATRIX SD.
C     OUTPUTS ARE SPECTRUM MATRIX P(I) AT FREQUENCIES I/(2*H)
C     (I=0,1,...,H).
C
cc      !DEC$ ATTRIBUTES DLLEXPORT :: MULRSPF
c      USE DFLIB
C
	IMPLICIT REAL*8(A-H,O-W)
      IMPLICIT COMPLEX*16(X-Z)
      INTEGER H,H1
c      COMMON G,GR,GI,LG,H,JJF
c      DIMENSION SD(10,10),A(30,10,10),B(30,10,10),X(10,10),Y(10,10)
c      DIMENSION Z(10,10),G(31),CH(10,10)
      DIMENSION SD(IP,IP),A(L,IP,IP),B(K,IP,IP),X(IP,IP),Y(IP,IP,H+1)
      DIMENSION Z(IP,IP),G(L+K+1),CH(IP,IP,H+1)
C     INPUT / OUTPUT DATA FILE OPEN
c	CALL SETWND
c	CALL FLOPN2(NFL)
c	IF (NFL.EQ.0) GO TO 999
C     ABSOLUTE DIMENSIONS USED FOR SUBROUTINE CALL
c      MJ0=30
c      MJ1=10
      CST0=0.0D-00
      CST1=1.0D-00
C     H SPECIFICATION
c      READ(5,1) H
C     SD AND A INPUT
C     THE OUTPUTS OF PROGRAM 5.3.2 FPEC(WITH IL=0) CAN BE USED AS THE
C     FOLLWOING INPUTS WITH K=0.
c      READ(5,1) N,L,IP
C     SD INPUT
c      CALL REMATX(SD,IP,IP,1,MJ1,MJ1)
c      IF(L.LE.0) GO TO 300
C     A INPUT
c      CALL REMAT3(A,L,IP,IP,1,MJ0,MJ1,MJ1)
C     K INPUT
c  300 READ(5,1) K
c      IF(K.LE.0) GO TO 310
C     B INPUT
c      CALL REMAT3(B,K,IP,IP,1,MJ0,MJ1,MJ1)
  310 H1=H+1
C     INITIAL CONDITION PRINT OUT
c      WRITE(6,59)
c      WRITE(6,60)
c      WRITE(6,61) H,N,L,IP,K
c      WRITE(6,161)
c      CALL SUBMPR(SD,IP,IP,MJ1,MJ1)
c      IF(L.LE.0) GO TO 400
C     A PRINT OUT
c      WRITE(6,420)
c      CALL PRMAT3(A,L,IP,IP,0,MJ0,MJ1,MJ1)
c  400 IF(K.LE.0) GO TO 410
C     B PRINT OUT
c      WRITE(6,430)
c      CALL PRMAT3(B,K,IP,IP,0,MJ0,MJ1,MJ1)
C     SPECTRUM COMPUTATION
  410 DO 10 JF=1,H1
      JJF=JF
C     SD STORE
      DO 631 II=1,IP
      DO 631 JJ=1,IP
c  631 Y(II,JJ)=SD(II,JJ)
  631 Y(II,JJ,JF)=SD(II,JJ)
      IF(K.GT.0) GO TO 100
      DO 110 II=1,IP
      DO 110 JJ=1,IP
  110 Z(II,JJ)=SD(II,JJ)
      GO TO 224
C     BF COMPUTATION
  100 DO 20 II=1,IP
      DO 21 JJ=1,IP
      IF(II.NE.JJ) GO TO 22
      G(1)=CST1
      GO TO 23
   22 G(1)=CST0
   23 DO 25 I=1,K
      I1=I+1
   25 G(I1)=B(I,II,JJ)
   24 LG=K
c      CALL FGER1
      CALL FGER1(G,GR,GI,LG,H,JJF)
      X(II,JJ)=DCMPLX(GR,GI)
   21 CONTINUE
   20 CONTINUE
C     BF*SD*CONJG(BF') COMPUTATION
c      CALL XYCTRX(X,Y,Z,IP,IP,MJ1,MJ1)
      CALL XYCTRX(X,Y(1,1,JF),Z,IP,IP)
  224 IF(L.GT.0) GO TO 120
      DO 130 II=1,IP
      DO 130 JJ=1,IP
c  130 Y(II,JJ)=Z(II,JJ)
  130 Y(II,JJ,JF)=Z(II,JJ)
      GO TO 244
C     AF COMPUTATION
  120 DO 40 II=1,IP
      DO 41 JJ=1,IP
      IF(II.NE.JJ) GO TO 42
      G(1)=CST1
      GO TO 43
   42 G(1)=CST0
   43 DO 45 I=1,L
      I1=I+1
   45 G(I1)=-A(I,II,JJ)
   44 LG=L
c      CALL FGER1
      CALL FGER1(G,GR,GI,LG,H,JJF)
      X(II,JJ)=DCMPLX(GR,GI)
   41 CONTINUE
   40 CONTINUE
C     INVERSE OF AF (COMPLEX) COMPUTATION
c      CALL INVDET(X,XDET,IP,MJ1)
      CALL INVDETC(X,XDET,IP)
C     (INVERSE OF AF)*(BF*SD*CONJG(BF'))*CONJG((INVERSE OF AF)')
C     COMPUTATION
c      CALL XYCTRX(X,Z,Y,IP,IP,MJ1,MJ1)
      CALL XYCTRX(X,Z,Y(1,1,JF),IP,IP)
C     SIMPLE COHERENCE COMPUTATION
c  244 CH(1,1)=CST1
  244 CH(1,1,JF)=CST1
cc      IF(IP.EQ.1) GO TO 260
      IF(IP.EQ.1) GO TO 10
      DO 50 II=2,IP
      IM1=II-1
c      RYI=DREAL(Y(II,II))
      RYI=DREAL(Y(II,II,JF))
      DO 51 JJ=1,IM1
c      RYJ=DREAL(Y(JJ,JJ))
c      RRYIJ=DREAL(Y(II,JJ))
c      RIYIJ=DIMAG(Y(II,JJ))
c      CH(II,JJ)=(RRYIJ**2+RIYIJ**2)/(RYI*RYJ)
c   51 CH(JJ,II)=CH(II,JJ)
c   50 CH(II,II)=CST1
      RYJ=DREAL(Y(JJ,JJ,JF))
      RRYIJ=DREAL(Y(II,JJ,JF))
      RIYIJ=DIMAG(Y(II,JJ,JF))
      CH(II,JJ,JF)=(RRYIJ**2+RIYIJ**2)/(RYI*RYJ)
   51 CH(JJ,II,JF)=CH(II,JJ,JF)
   50 CH(II,II,JF)=CST1
C     RATIONAL SPECTRUM AND SIMPLE COHERENCE PRINT OUT
cc  260 JFM1=JF-1
c      WRITE(6,65) JFM1
c      WRITE(6,66)
c      CALL PRCPMA(Y,IP,IP,MJ1,MJ1)
c      WRITE(6,67)
c      CALL SUBMPR(CH,IP,IP,MJ1,MJ1)
C
   10 CONTINUE
c	CALL FLCLS2(NFL)
c  999 CONTINUE
c    1 FORMAT(10I5)
c   59 FORMAT(1H ,42HPROGRAM 5.4.2   MULTIPLE RATIONAL SPECTRUM)
c   60 FORMAT(1H ,17HINITIAL CONDITION)
c   61 FORMAT(1H ,2HH=,I5,5X,2HN=,I5,5X,2HL=,I5,5X,3HIP=,I5,5X,2HK=,I5)
c   62 FORMAT(1H ,2HI=,I5)
   65 FORMAT(///1H ,2HF=,I5)
   66 FORMAT(1H ,5X,17HRATIONAL SPECTRUM)
   67 FORMAT(/1H ,5X,16HSIMPLE COHERENCE)
c  161 FORMAT(//1H ,7HSD(I,J))
c  420 FORMAT(//1H ,6HA(I,J))
c  430 FORMAT(//1H ,6HB(I,J))
c    2 FORMAT(4D20.10)
      RETURN
      END
C
c      SUBROUTINE XYCTRX(X,Y,Z,MM,NN,MJ1,MJ2)
      SUBROUTINE XYCTRX(X,Y,Z,MM,NN)
C     Z=X*Y*CONJG(X')
C     Y,Z: HERMITIAN
C     (UPPER LEFT MM X MM OF Z)=(UPPER LEFT MM X NN OF X)*(UPPER LEFT
C     NN X NN OF Y)*CONJG((UPPER LEFT MM X NN OF X)')
C     (MJ1,MJ2): ABSOLUTE DIMENSION OF X IN THE MAIN ROUTINE
C     (MJ2,MJ2): ABSOLUTE DIMENSION OF Y IN THE MAIN ROUTINE
C     (MJ1,MJ1): ABSOLUTE DIMENSION OF Z IN THE MAIN ROUTINE
C     MM,NN: SHOULD BE LESS THAN 11.
      IMPLICIT COMPLEX*16(X-Z)
c      DIMENSION X(MJ1,MJ2),Y(MJ2,MJ2),Z(MJ1,MJ1)
c      DIMENSION Y1(10,10)
      DIMENSION X(MM,NN),Y(NN,NN),Z(MM,MM)
      DIMENSION Y1(MM,NN)
	DOUBLE PRECISION CST0
      CST0=0.0D-00
      DO 10 I=1,MM
      DO 10 J=1,NN
      XSUM=CST0
      DO 12 K=1,NN
   12 XSUM=XSUM+X(I,K)*Y(K,J)
   10 Y1(I,J)=XSUM
      DO 110 I=1,MM
      DO 110 J=1,I
      XSUM=CST0
      DO 112 K=1,NN
  112 XSUM=XSUM+Y1(I,K)*DCONJG(X(J,K))
      Z(I,J)=XSUM
  110 Z(J,I)=DCONJG(Z(I,J))
      RETURN
      END
