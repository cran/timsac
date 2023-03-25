      SUBROUTINE MULCORF(X1,N,K,LAGH1,SM,C,CN)
C
      INCLUDE 'timsac.h'
C
C     PROGRAM 5.1.2   MULTIPLE CORRELATION
C-----------------------------------------------------------------------
C     ** DESIGNED BY H. AKAIKE, THE INSTITUTE OF STATISTICAL MATHEMATICS
C     ** PROGRAMMED BY E. ARAHATA, THE INSTITUTE OF STATISTICAL MATHEMAT
C         TOKYO
C     ** DATE OF THE LATEST REVISION: MARCH 25, 1977
C     ** THIS PROGRAM WAS ORIGINALLY PUBLISHED IN
C         "DAINAMIKKU SISTEMU NO TOKEI-TEKI KAISEKI TO SEIGYO (STATISTICA
C         ANALYSIS AND CONTROL OF DYNAMIC SYSTEMS)" BY H. AKAIKE AND
C         T. NAKAGAWA, SAIENSU-SHA, TOKYO, 1972 (IN JAPANESE)
C-----------------------------------------------------------------------
C     THIS PROGRAM REQUIRES FOLLOWING INPUTS:
C     N: LENGTH OF DATA
C     K: DIMENSION OF THE OBSERVATION VECTOR
C     LAGH: MAXIMUM LAG
C     ISW: ISW=1...ROWWISE DATA INPUT
C           ISW=2...COLUMNWISE DATA INPUT
C     DFORM: INPUT FORMAT SPECIFICATION STATEMENT IN ONE CARD,
C     FOR EXAMPLE
C     (8F10.4)
C     (X1(S,I); S=1,...,N, I=1,...,K): ORIGINAL DATA MATRIX.
C     THE OUTPUTS ARE (CIJ(L): L=0,1,...,LAGH) (I=1,...,K; J=1,...,K),
C     WHERE CIJ(L)=COVARIANCE(XI(S+L),XJ(S)),
C     AND THEIR NORMALIZED (CORRELATION) VALUES.
C
cxx      IMPLICIT REAL*8 (A-H,O-Z)
c      DIMENSION X1(2000,10)
c      DIMENSION X(2000),Y(2000)
c      DIMENSION C1(501),C2(501),CN1(501),CN2(501)
c      DIMENSION SM(10),C0(10)
c      REAL*4 DFORM
c      DIMENSION DFORM(20)
cxx      DIMENSION X1(N,K),X2(N,K)
cxx      DIMENSION X(N),Y(N)
cxx      DIMENSION C(LAGH1,K,K),CN(LAGH1,K,K)
cxx      DIMENSION C1(LAGH1),C2(LAGH1),CN1(LAGH1),CN2(LAGH1)
cxx      DIMENSION SM(K),C0(K)
      INTEGER N, K, LAGH1
      DOUBLE PRECISION X1(N,K), SM(K), C(LAGH1,K,K), CN(LAGH1,K,K)
c local
      DOUBLE PRECISION X2(N,K), X(N), Y(N), C1(LAGH1), C2(LAGH1),
     1                 CN1(LAGH1), CN2(LAGH1), C0(K), CX0, CY0, XMEAN
C
C     INPUT / OUTPUT DATA FILE OPEN
c      CHARACTER(100) DFNAM
c      DFNAM='mulcor.out'
c      CALL SETWND
c      CALL FLOPN3(DFNAM,NFL)
c      IF (NFL.EQ.0) GO TO 999
C     INITIAL CONDITION INPUT AND OUTPUT
c      READ(5,1) N,LAGH,K,ISW
c      LAGH1=LAGH+1
c      WRITE(6,50)
c      WRITE(6,51)
c      WRITE(6,52) N,LAGH,K
C     INITIAL CONDITION PUNCH OUT
c      WRITE(7,1) N,LAGH,K
C     INPUT FORMAT SPECIFICATION
c      READ(5,4) (DFORM(I),I=1,20)
c    4 FORMAT(20A4)
C     ORIGINAL DATA INPUT AND OUTPUT
c      GO TO(8,9),ISW
c    8 DO 208 I=1,N
c  208 READ(5,DFORM) (X(I,II),II=1,K)
c      GO TO 400
c    9 DO 209 II=1,K
c  209 READ(5,DFORM) (X(I,II),I=1,N)
c  400 WRITE(6,53)
c      WRITE(6,54)
c      WRITE(6,154)
c      DO 220 I=1,N
c  220 WRITE(6,55) I,(X(I,II),II=1,K)
C
C     MEAN DELETION
      DO 300 II=1,K
         DO 310 I=1,N
cxx  310  X(I)=X1(I,II)
         X(I)=X1(I,II)
  310    CONTINUE
         CALL DMEADL(X,N,XMEAN)
         SM(II)=XMEAN
         DO 320 I=1,N
cxx  320   X2(I,II)=X(I)
         X2(I,II)=X(I)
  320    CONTINUE
  300 CONTINUE
C
C     COVARIANCE COMPUTATION
      DO 10 II=1,K
         DO 110 I=1,N
c  110 X(I)=X1(I,II)
cxx  110 X(I)=X2(I,II)
         X(I)=X2(I,II)
  110    CONTINUE
C     AUTO COVARIANCE COMPUTATION
      CALL CROSCO(X,X,N,C1,LAGH1)
C     NORMALIZATION
      C0(II)=C1(1)
      CX0=C0(II)
      CALL CORNOM(C1,CN1,LAGH1,CX0,CX0)
C     AUTO COVARIANCE PRINT OUT
c      WRITE(6,162) II,II,SM(II)
c      WRITE(6,163)
c      CALL PRCOL2(C1,CN1,1,LAGH1,1)
C     AUTO COVARIANCE PUNCH OUT
c      WRITE(7,1) II,II
c      WRITE(7,2) (C1(I),I=1,LAGH1)
      DO 115 I=1,LAGH1
         C(I,II,II)=C1(I)
         CN(I,II,II)=CN1(I)
  115 CONTINUE
      IF(II.EQ.1) GO TO 10
      IM1=II-1
      DO 11 JJ=1,IM1
         DO 120 I=1,N
c  120 Y(I)=X1(I,JJ)
cxx  120 Y(I)=X2(I,JJ)
         Y(I)=X2(I,JJ)
  120    CONTINUE
C     CROSS COVARIANCE COMPUTATION
      CALL CROSCO(X,Y,N,C1,LAGH1)
      CALL CROSCO(Y,X,N,C2,LAGH1)
C �@�@NORMALIZATION
      CX0=C0(II)
      CY0=C0(JJ)
      CALL CORNOM(C1,CN1,LAGH1,CX0,CY0)
      CALL CORNOM(C2,CN2,LAGH1,CX0,CY0)
C     CROSS COVARIANCE PRINT OUT
c      WRITE(6,165) II,JJ
c      WRITE(6,166)
c      CALL PRCOL4(C1,CN1,C2,CN2,1,LAGH1,1)
C     CROSS COVARIANCE PUNCH OUT
c      WRITE(7,1) II,JJ
c      WRITE(7,2) (C1(I),I=1,LAGH1)
c      WRITE(7,1) JJ,II
c      WRITE(7,2) (C2(I),I=1,LAGH1)
      DO 125 I=1,LAGH1
         C(I,II,JJ)=C1(I)
         C(I,JJ,II)=C2(I)
         CN(I,II,JJ)=CN1(I)
         CN(I,JJ,II)=CN2(I)
  125 CONTINUE
   11 CONTINUE
   10 CONTINUE
c      CALL FLCLS3(NFL)
c  999 CONTINUE
c    1 FORMAT(10I5)
c    2 FORMAT(4D20.10)
c   50 FORMAT(1H ,36HPROGRAM 5.1.2   MULTIPLE CORRELATION)
c   51 FORMAT(1H ,17HINITIAL CONDITION)
c   52 FORMAT(1H ,2HN=,I5,5X,5HLAGH=,I5,5X,2HK=,I5)
c   53 FORMAT(1H ,13HORIGINAL DATA)
c   54 FORMAT(1H ,4X,1HI,2X,8HX1(I,II))
c  154 FORMAT(1H ,16X,91H1	   2	     3	       4	 5
c     A	 6	   7	     8	       9	10)
c   55 FORMAT(1H ,I5,2X,10F10.4)
c  162 FORMAT(//1H ,14HAUTOCOVARIANCE,5X,6HCIJ(L),5X,2HI=,I5,5X,2HJ=,I5,5
c     AX,5HMEAN=,D15.5)
c  163 FORMAT(1H ,4X,1HL,5X,6HCIJ(L),8X,10HNORMALIZED)
c  165 FORMAT(//1H ,16HCROSS COVARIANCE,5X,6HCIJ(L),5X,2HI=,I5,5X,2HJ=,I5
c     A)
c  166 FORMAT(1H ,4X,1HL,5X,6HCIJ(L),8X,10HNORMALIZED,4X,6HCJI(L),8X,10HN
c     AORMALIZED)
      RETURN
      END SUBROUTINE
