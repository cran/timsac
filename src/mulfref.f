      SUBROUTINE MULFRFF(K,INW,N,LAGH1,IP0,P,X,C,S,G,PH,PCH,R,CHM)
C
      INCLUDE 'timsac_f.h'
C
C     PROGRAM 5.2.4   FREQUENCY RESPONSE FUNCTION (MULTIPLE CHANNEL)
C-----------------------------------------------------------------------
C      SUBROUTINE FQCPIV(X,XDET,MM,MJ)
C      SUBROUTINE MPHASE(C,S,OARC,PH,K,JJF)
C      SUBROUTINE MULARC(C,S,ARC,K)
C      SUBROUTINE MULERR(PCH,R,N,LAGH1,K,JJF,D1,D2)
C      SUBROUTINE MULPAC(ARC,OARC,PH,K,JJF)
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
C     THIS PROGRAM COMPUTES MULTIPLE FREQUENCY RESPONSE FUNCTION, GAIN,
C     PHASE, MULTIPLE COHERENCY, PARTIAL COHERENCY AND RELATIVE ERROR
C     STATISTICS.
C     A CARD WITH THE TOATL NUMBER(K) OF INPUT VARIABLES AND ANOTHER
C     WITH SPECIFICATION OF INPUT VARIABLES(INW(I),I=1,K) AND OUTPUT
C     VARIABLE(INW(K+1)) SHOULD BE ADDED ON TOP OF THE OUTPUT OF
C     PROGRAM 5.2.2 MULSPE TO FORM THE INPUT TO THIS PROGRAM.
C     WITHIN IP0 VARIABLES OF MULSPE OUTPUT, ONLY THOSE K+1 INW(I)-TH
C     VARIABLES ARE TAKEN INTO COMPUTATION.
C
cc      !DEC$ ATTRIBUTES DLLEXPORT :: MULFRFF
c      USE DFLIB
C
	IMPLICIT REAL*8(A-H,O-W)
      IMPLICIT COMPLEX*16(X-Z)
c      DIMENSION P(10,10),X(10,10),C(10),S(10),G(10)
c      DIMENSION OARC(10),PH(10),PCH(10),R(10),INW(10)
      DIMENSION P(LAGH1,IP0,IP0),X(IP0,IP0,LAGH1)
	DIMENSION C(K,LAGH1),S(K,LAGH1),G(K,LAGH1)
      DIMENSION OARC(K),PH(K,LAGH1),PCH(K,LAGH1),R(K,LAGH1),INW(K+1)
      DIMENSION CHM(LAGH1)
	DIMENSION XFR(IP0,IP0,LAGH1)
C     INPUT / OUTPUT DATA FILE OPEN
c	CALL SETWND
c	CALL FLOPN2(NFL)
c	IF (NFL.EQ.0) GO TO 999
C     ABSOLUTE DIMENSION USED FOR SUBROUTINE CALL
c      MJ=10
C     INPUT OUTPUT VARIABLE SPECIFICATION
c      READ(5,1) K
      K1=K+1
c      READ(5,1) (INW(I),I=1,K1)
C     FOLLOWING INPUTS ARE OUTPUTS OF PROGRAM 5.2.2 MULSPE.
c      READ(5,1) N,LAGH,IP0
c      LAGH1=LAGH+1
C     INITIAL CONDITION PRINT OUT
c      WRITE(6,55)
c      WRITE(6,56)
c      WRITE(6,57) N,LAGH,K
c      WRITE(6,259) (INW(I),I=1,K1)
C     COMPUTATION START
      DO 10 JF=1,LAGH1
      JJF=JF
c      JFM1=JF-1
c      WRITE(6,58) JFM1
C     SPECTRUM INPUT
c      CALL REMATX(P,IP0,IP0,1,MJ,MJ)
C     REAL TO COMPLEX TRANSFORMATION
      DO 401 I=1,IP0
c      X(I,I)=P(I,I)
      X(I,I,JF)=P(JF,I,I)
      IF(I.EQ.1) GO TO 401
      IM1=I-1
      DO 402 J=1,IM1
c      X(I,J)=DCMPLX(P(I,J),P(J,I))
c  402 X(J,I)=DCONJG(X(I,J))
      X(I,J,JF)=DCMPLX(P(JF,I,J),P(JF,J,I))
  402 X(J,I,JF)=DCONJG(X(I,J,JF))
  401 CONTINUE
C     MATRIX REARRANGEMENT AND PRINT OUT (COMPLEX)
c      CALL REARRA(X,INW,IP0,K1,MJ)
      CALL REARRAC(X(1,1,JF),INW,IP0,K1)
c      WRITE(6,159)
c      CALL PRCPMA(X,K1,K1,MJ,MJ)
C     FREQUENCY RESPONSE FUNCTION COMPUTATION
c      P00=DREAL(X(K1,K1))
c      CALL FQCPIV(X,XDET,K,MJ)
      P00=DREAL(X(K1,K1,JF))
c      CALL FQCPIV(X(1,1,JF),XDET,K,MJ)
      DO 31 I=1,IP0
	DO 30 II=1,IP0
      XFR(I,II,JF)=X(I,II,JF)
   30 CONTINUE
   31 CONTINUE
      CALL FQCPIV(XFR(1,1,JF),XDET,K,IP0)
      DO 20 I=1,K
c      C(I)=DREAL(X(I,K1))
c   20 S(I)=-DIMAG(X(I,K1))
      C(I,JF)=DREAL(XFR(I,K1,JF))
   20 S(I,JF)=-DIMAG(XFR(I,K1,JF))
C     GAIN COMPUTATION
      DO 21 I=1,K
c   21 G(I)=DSQRT(C(I)**2+S(I)**2)
   21 G(I,JF)=DSQRT(C(I,JF)**2+S(I,JF)**2)
C     PHASE COMPUTATION
c      CALL MPHASE(C,S,OARC,PH,K,JJF)
      IF(JJF.NE.1) THEN
	 DO 24 I=1,K
   24	 PH(I,JF)=PH(I,JF-1)
	END IF
      CALL MPHASE(C(1,JF),S(1,JF),OARC,PH(1,JF),K,JJF)
C     PARTIAL COHERENCY AND MULTIPLE COHERENCY COMPUTATION
c      EP=DREAL(X(K1,K1))
      EP=DREAL(XFR(K1,K1,JF))
      DO 22 I=1,K
c      G2=G(I)**2
c      G3=G2+EP*X(I,I)
      G2=G(I,JF)**2
      G3=G2+EP*XFR(I,I,JF)
      IF(G3.NE.0.0) GO TO 23
c      PCH(I)=100.0D-00
      PCH(I,JF)=100.0D-00
      GO TO 22
c   23 PCH(I)=G2/G3
   23 PCH(I,JF)=G2/G3
   22 CONTINUE
c      CHM=1.0D-00-EP/P00
      CHM(JF)=1.0D-00-EP/P00
C     RELATIVE ERROR STATISTICS COMPUTATION
c      CALL MULERR(PCH,R,N,LAGH1,K,JJF,D1,D2)
      CALL MULERR(PCH(1,JF),R(1,JF),N,LAGH1,K,JJF,D1,D2)
C     FREQUENCY RESPONSE FUNCTION, GAIN, PHASE, PARTIAL COHERENCY,
C     MULTIPLE COHERENCY, RELATIVE ERROR STATISTICS PRINT OUT
c      WRITE(6,60)
c      WRITE(6,61)
c      CALL PRCOL6(C,S,G,PH,PCH,R,1,K,0)
c      WRITE(6,65) CHM
c      WRITE(6,65) CHM(JF)
   10 CONTINUE
c      CALL FLCLS2(NFL)
c  999 CONTINUE
c    1 FORMAT(10I5)
c   55 FORMAT(1H ,62HPROGRAM 5.2.4   FREQUENCY RESPONSE FUNCTION (MULTIPL
c     AE CHANNEL))
c   56 FORMAT(1H ,17HINITIAL CONDITION)
c   57 FORMAT(1H ,2HN=,I5,5X,5HLAGH=,I5,5X,2HK=,I5)
c   58 FORMAT(//1H ,2HF=,I5)
c   60 FORMAT(//1H ,4X,1HI,3X,27HFREQUENCY RESPONSE FUNCTION,10X,4HGAIN,9
c     AX,5HPHASE,7X,7HPARTIAL,6X,8HRELATIVE)
c   61 FORMAT(1H ,12X,9HREAL PART,4X,10HIMAG. PART,33X,9HCOHERENCY,9X,5HE
c     ARROR)
c   65 FORMAT(1H ,69X,8HMULTIPLE/1H ,68X,9HCOHERENCY/1H ,63X,D14.5)
c  159 FORMAT(1H ,28HSPECTRUM MATRIX (REARRANGED))
c  259 FORMAT(/1H ,6HINW(I),5X,10I5)
      RETURN
      END SUBROUTINE
C
      SUBROUTINE FQCPIV(X,XDET,MM,MJ)
C     THIS SUBROUTINE COMPUTES MULTIPLE FREQUENCY RESPONSE FUNCTION.
C     MM: THE TOTAL NUMBER OF INPUTS (LESS THAN 10)
C     MJ: ABSOLUTE DIMENSION OF X IN THE MAIN ROUTINE
      IMPLICIT COMPLEX*16(X)
      DIMENSION X(MJ,MJ)
      DIMENSION IDS(10)
      XDET=1.0D-00
      MP1=MM+1
      DO 10 L=1,MM
C     PIVOTING AT L-TH STAGE
      XMAXP=0.10000D-10
      MAXI=0
      DO 110 I=L,MM
      IF(CDABS(XMAXP).GE.CDABS(X(I,L))) GO TO 110
      XMAXP=X(I,L)
      MAXI=I
  110 CONTINUE
      IDS(L)=MAXI
      IF(MAXI.EQ.L) GO TO 120
      IF(MAXI.GT.0) GO TO 121
      XDET=0.0D-00
      GO TO 140
C     ROW INTERCHANGE
  121 DO 14 J=1,MP1
      XC=X(MAXI,J)
      X(MAXI,J)=X(L,J)
   14 X(L,J)=XC
      XDET=-XDET
  120 XDET=XDET*XMAXP
      XC=1.0D-00/XMAXP
      X(L,L)=1.0D-00
      DO 11 J=1,MP1
   11 X(L,J)=X(L,J)*XC
      DO 12 I=1,MP1
      IF(I.EQ.L) GO TO 12
      XC=X(I,L)
      X(I,L)=0.0D-00
      DO 13 J=1,MP1
   13 X(I,J)=X(I,J)-XC*X(L,J)
   12 CONTINUE
   10 CONTINUE
      IF(MM.GT.1) GO TO 123
      GO TO 140
C     COLUMN INTERCHANGE
  123 MM1=MM-1
      DO 130 J=1,MM1
      MMJ=MM-J
      JJ=IDS(MMJ)
      IF(JJ.EQ.MMJ) GO TO 130
      DO 131 I=1,MP1
      XC=X(I,JJ)
      X(I,JJ)=X(I,MMJ)
  131 X(I,MMJ)=XC
  130 CONTINUE
  140 RETURN
      END SUBROUTINE
C
      SUBROUTINE MPHASE(C,S,OARC,PH,K,JJF)
C     THIS SUBROUTINE COMPUTES PHASE.
C     (MULTIPLE CHANNEL)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION C(K),S(K),OARC(K),PH(K)
c      DIMENSION ARC(10)
      DIMENSION ARC(K)
C     ARCTANGENT COMPUTATION
      CALL MULARC(C,S,ARC,K)
C     PHASE COMPUTATION
      CALL MULPAC(ARC,OARC,PH,K,JJF)
      RETURN
      END SUBROUTINE
C
      SUBROUTINE MULARC(C,S,ARC,K)
C     THIS SUBROUTINE COMPUTES RAW PHASE.
C     (MULTIPLE CHANNEL)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION C(K),S(K),ARC(K)
      PI=3.1415926536
      CST5=0.5D-00
      DO 10 I=1,K
cc      IF(C(I)) 11,12,13
cc   11 IF(S(I)) 14,15,16
cc   12 IF(S(I)) 17,18,19
      IF(C(I).EQ.0) GO TO 12
      IF(C(I).GT.0) GO TO 13
   11 IF(S(I).LT.0) GO TO 14
      IF(S(I).EQ.0) GO TO 15
      IF(S(I).GT.0) GO TO 16
   12 IF(S(I).LT.0) GO TO 17
      IF(S(I).EQ.0) GO TO 18
      IF(S(I).GT.0) GO TO 19
   13 ARC(I)=DATAN(S(I)/C(I))
      GO TO 10
   14 ARC(I)=DATAN(S(I)/C(I))-PI
      GO TO 10
   15 ARC(I)=-PI
      GO TO 10
   16 ARC(I)=DATAN(S(I)/C(I))+PI
      GO TO 10
   17 ARC(I)=-PI*CST5
      GO TO 10
   18 ARC(I)=0.0D-00
      GO TO 10
   19 ARC(I)=PI*CST5
   10 CONTINUE
      RETURN
      END SUBROUTINE
C
      SUBROUTINE MULERR(PCH,R,N,LAGH1,K,JJF,D1,D2)
C     THIS SUBROUTINE COMPUTES RELATIVE ERROR STATISTICS.
C     (MULTIPLE CHANNEL)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION PCH(K),R(K)
      CST0=0.0D-00
      CST1=1.0D-00
      CST100=100.0D-00
      IF(JJF.NE.1) GO TO 30
C     CONSTANTS D1,D2 COMPUTATION
      LAGH=LAGH1-1
      CALL SUBD12(N,LAGH,K,D1,D2)
C     RELATIVE ERROR STATISTICS COMPUTATION
   30 DO 20 I=1,K
      IF(PCH(I).LE.CST0) GO TO 22
      IF(PCH(I).GT.CST1) GO TO 22
      E1=CST1/PCH(I)-CST1
      ER=DSQRT(E1)
      IF(JJF.EQ.1) GO TO 23
      IF(JJF.EQ.LAGH1) GO TO 23
      R(I)=D2*ER
      GO TO 20
   23 R(I)=D1*ER
      GO TO 20
   22 R(I)=CST100
   20 CONTINUE
      RETURN
      END SUBROUTINE
C
      SUBROUTINE MULPAC(ARC,OARC,PH,K,JJF)
C     THIS SUBROUTINE MAKES PHASE CURVE CONTINUOUS.
C     (MULTIPLE CHANNEL)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION ARC(K),OARC(K),PH(K)
      PI=3.1415926536
      PI2=PI+PI
      IF(JJF.NE.1) GO TO 20
      DO 9 I=1,K
      PH(I)=ARC(I)
    9 OARC(I)=ARC(I)
      GO TO 30
   20 DO 10 I=1,K
      DK=ARC(I)-OARC(I)
      IF(DK.GT.PI) GO TO 11
      IF(DK.LT.-PI) GO TO 12
      PH(I)=PH(I)+DK
      GO TO 10
   11 PH(I)=PH(I)+DK-PI2
      GO TO 10
   12 PH(I)=PH(I)+DK+PI2
   10 OARC(I)=ARC(I)
   30 CONTINUE
      RETURN
      END SUBROUTINE
