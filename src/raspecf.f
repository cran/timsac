      SUBROUTINE RASPECF(H,L,K,SGME2,A,B,PXX)
C
      INCLUDE 'timsac.h'
C
C     PROGRAM 5.4.1   RATIONAL SPECTRUM
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
C     THIS PROGRAM COMPUTES POWER SPECTRUM OF AR-MA PROCESS
C     X(N)=A(1)X(N-1)+...+A(L)X(N-L)+E(N)+B(1)E(N-1)+...+B(K)E(N-K)
C     WHERE E(N) IS A WHITE NOISE WITH ZERO MEAN AND VARIANCE EQUAL TO
C     SGME2.  OUTPUTS PXX(I) ARE GIVEN AT FREQUENCIES I/(2*H)
C     I=0,1,...,H.
C     REQUIRED INPUTS ARE:
C     L,K,H,SGME2,(A(I),I=1,L), AND (B(I),I=1,K).
C     0 IS ALLOWABLE AS L AND/OR K.
C
cxx      IMPLICIT REAL*8(A-H,O-Z)
cxx      INTEGER H,H1
CC      DIMENSION A(501),B(501)
CC      DIMENSION G(501),GR1(501),GI1(501),GR2(501),GI2(501)
CC      DIMENSION PXX(510)
cxx      DIMENSION A(L),B(K)
CC      DIMENSION G(MAX(L,K)+1),GR1(H+1),GI1(H+1),GR2(H+1),GI2(H+1)
cxx      DIMENSION G(L+K+1),GR1(H+1),GI1(H+1),GR2(H+1),GI2(H+1)
cxx      DIMENSION PXX(H+1)
      INTEGER H, L, K
      DOUBLE PRECISION SGME2, A(L), B(K), PXX(H+1)
c local
      INTEGER I, I1, K1, L1, H1
      DOUBLE PRECISION G(L+K+1), GR1(H+1), GI1(H+1), GR2(H+1), GI2(H+1)
C
C     INPUT / OUTPUT DATA FILE OPEN
CC      CALL SETWND
CC      CALL FLOPN2(NFL)
CC      IF (NFL.EQ.0) GO TO 999
C     H SPECIFICATION
CC      READ(5,1) H
C     SGME2 AND A INPUT
C     THE OUTPUTS OF PROGRAM 5.3.1 FPE AUTO CAN BE USED AS THE FOLLOWING
C     INPUTS WITH K=0.
CC      READ(5,1) N,L
CC      READ(5,2) SGME2
CC      IF(L.LE.0) GO TO 300
CC      READ(5,2) (A(I),I=1,L)
C     K INPUT
CC  300 READ(5,1) K
CC      IF(K.LE.0) GO TO 310
CC      READ(5,2) (B(I),I=1,K)
CC  310 H1=H+1
      H1=H+1
      L1=L+1
      K1=K+1
      G(1)=1.0D-00
      IF(L.LE.0) GO TO 400
      DO 10 I=1,L
      I1=I+1
cxx   10 G(I1)=-A(I)
      G(I1)=-A(I)
   10 CONTINUE
  400 CALL FOUGER(G,L1,GR1,GI1,H1)
      G(1)=1.0D-00
      IF(K.LE.0) GO TO 410
      DO 20 I=1,K
      I1=I+1
cxx   20 G(I1)=B(I)
      G(I1)=B(I)
   20 CONTINUE
  410 CALL FOUGER(G,K1,GR2,GI2,H1)
      DO 30 I=1,H1
cxx   30 PXX(I)=(GR2(I)**2+GI2(I)**2)/(GR1(I)**2+GI1(I)**2)*SGME2
      PXX(I)=(GR2(I)**2+GI2(I)**2)/(GR1(I)**2+GI1(I)**2)*SGME2
   30 CONTINUE
CC      WRITE(6,60)
CC      WRITE(6,160)
CC      WRITE(6,61) L,K,H
CC      WRITE(6,164) SGME2
CC      IF(L.LE.0) GO TO 500
CC      WRITE(6,62)
CC      CALL PRCOL1(A,1,L,0)
CC  500 IF(K.LE.0) GO TO 510
CC      WRITE(6,63)
CC      CALL PRCOL1(B,1,K,0)
CC  510 WRITE(6,64)
CC      CALL PRCOL1(PXX,1,H1,1)
CC      CALL FLCLS2(NFL)
CC  999 CONTINUE
CC    1 FORMAT(10I5)
CC    2 FORMAT(4D20.10)
CC   60 FORMAT(1H ,33HPROGRAM 5.4.1   RATIONAL SPECTRUM)
CC   61 FORMAT(1H ,2HL=,I5,2X,2HK=,I5,2X,2HH=,I5)
CC   62 FORMAT(1H ,4X,1HI,12X,4HA(I))
CC   63 FORMAT(1H ,4X,1HI,12X,4HB(I))
CC   64 FORMAT(1H ,4X,1HI,10X,6HPXX(I))
CC  160 FORMAT(1H ,17HINITIAL CONDITION)
CC  164 FORMAT(1H ,6HSGME2=,D12.5)
      RETURN
      END
