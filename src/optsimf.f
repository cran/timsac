      SUBROUTINE OPTSIMF(NS,M,IR,L,A,B,G,W,X,Y,XMEAN,YMEAN,XS2,YS2,
     * XS2MEA,YS2MEA,XVAR,YVAR)
C
      INCLUDE 'timsac.h'
C
cc      PROGRAM OPTSIM
C     PROGRAM 5.5.2   OPTIMAL CONTROL SIMULATION
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
C     THIS PROGRAM PERFORMS OPTIMAL CONTROL SIMULATION FOR THE
C     CONTROLLER DESIGNED BY PROGRAM 5.5.1 AND EVALUATES THE MEANS AND
C     VARIANCES OF THE CONTROLLED AND MANIPULATED VARIABLES X AND Y.
C     FOLLOWING CONSTANTS SHOULD BE PROVIDED BESIDES THE OUTPUTS OF
C     PROGRAM 5.5.1   OPTIMAL CONTROLLER DESIGN TO START THIS PROGRAM.
C     NS: NUMBER OF STEPS OF SIMULATION
C     INTP=1: TO SUPPRESS HISTORY OUTPUT
C     INTP=2: TO PRINT OUT THE HISTORY
C     THE SEQUENCE OF NS IR-DIMENSIONAL VECTORS W, REPRESENTING WHITE
C     NOISE (OR IMPULSE) IS ALSO REQUIRED AS INPUT.
C
cxx      IMPLICIT REAL*8(A-H,O-Z)
cc	DIMENSION A(75,5),B(75,5),Q1(5,5),R(5,5),G(5,75)
cc	DIMENSION X(5),XS(5),XS2(5)
cc	DIMENSION Y(5),YS(5),YS2(5)
cc	DIMENSION W(5),Z(75),C(75)
cc	DIMENSION XMEAN(5),XS2MEA(5),XVAR(5)
cc	DIMENSION YMEAN(5),YS2MEA(5),YVAR(5)
cxx      DIMENSION A(IR*M,IR),B(IR*M,L),G(L,IR*M)
cxx      DIMENSION X(IR,NS),XS(IR),XS2(IR)
cxx      DIMENSION Y(L,NS),YS(L),YS2(L)
cxx      DIMENSION W(IR,NS),Z(IR*M),C(IR*M)
cxx      DIMENSION XMEAN(IR),XS2MEA(IR),XVAR(IR)
cxx      DIMENSION YMEAN(L),YS2MEA(L),YVAR(L)
      INTEGER NS, M, IR, L
      DOUBLE PRECISION A(IR*M,IR), B(IR*M,L), G(L,IR*M), W(IR,NS),
     1                 X(IR,NS), Y(L,NS), XMEAN(IR), YMEAN(L), XS2(IR),
     2                 YS2(L), XS2MEA(IR), YS2MEA(L), XVAR(IR), YVAR(L)
c local
      INTEGER I, INS, IPR, MR, MR1
      DOUBLE PRECISION XS(IR), YS(L), Z(IR*M), C(IR*M), CST0, CST1,
     1                 ANS, BNS
C     INPUT / OUTPUT DATA FILE OPEN
cc	CALL SETWND
cc	CALL FLOPN2(NFL)
cc	IF (NFL.EQ.0) GO TO 999
C     ABSOLUTE DIMENSIONS USED FOR SUBROUTINE CALL
cc	MJ1=5
cc	MJ2=5
cc	MJ3=75
      CST0=0.0D-00
C     INITIAL CONDITION INPUT AND OUTPUT
cc	READ(5,1) NS,INTP
C     READING THE OUTPUTS OF PROGRAM 5.5.1 OPTDES
cc	READ(5,1) N,M,IR,L
      MR=M*IR
cc	CALL REMATX(Q1,IR,IR,1,MJ1,MJ1)
cc	CALL REMATX(R,L,L,1,MJ2,MJ2)
cc	CALL REMATX(A,MR,IR,1,MJ3,MJ1)
cc	CALL REMATX(B,MR,L,1,MJ3,MJ2)
cc	CALL REMATX(G,L,MR,1,MJ2,MJ3)
cc	WRITE(6,60)
cc	WRITE(6,61)
cc	WRITE(6,62) N,M,IR,L,NS
cc	WRITE(6,65)
cc	CALL SUBMPR(Q1,IR,IR,MJ1,MJ1)
cc	WRITE(6,66)
cc	CALL SUBMPR(R,L,L,MJ2,MJ2)
cc	WRITE(6,63)
cc	CALL SUBMPR(A,MR,IR,MJ3,MJ1)
cc	WRITE(6,64)
cc	CALL SUBMPR(B,MR,L,MJ3,MJ2)
cc	WRITE(6,67)
cc	CALL SUBMPR(G,L,MR,MJ2,MJ3)
C     INITIAL CONDITIONING
cxx      DO 6 I=1,IR
cc	X(I)=CST0
cxx      DO 5 J=1,NS
cxx    5 X(I,J)=CST0
cxx      XS(I)=CST0
cxx      XS2(I)=CST0
cxx    6 CONTINUE
cxx      DO 7 I=1,L
cc	Y(I)=CST0
cxx      DO 77 J=1,NS
cxx   77 Y(I,J)=CST0
cxx      YS(I)=CST0
cxx      YS2(I)=CST0
cxx    7 CONTINUE
cxx      DO 8 I=1,MR
cxx      C(I)=CST0
cxx    8 CONTINUE
      X(1:IR,1:NS)=CST0
      XS(1:IR)=CST0
      XS2(1:IR)=CST0
      Y(1:L,1:NS)=CST0
      YS(1:L)=CST0
      YS2(1:L)=CST0
      C(1:MR)=CST0
      MR1=MR-IR
C     START OF SIMULATION
C     NOISE INPUT
      DO 10 INS=1,NS
cc	READ(5,3) (W(I),I=1,IR)
C     X COMPUTATION
cc	CALL VECADL(C,W,IR)
      CALL VECADL(C,W(1,INS),IR)
      DO 9 I=1,IR
cc    9 X(I)=C(I)
cxx    9 X(I,INS)=C(I)
      X(I,INS)=C(I)
    9 CONTINUE
C     Y COMPUTATION
cc	CALL MULVER(G,C,Y,L,MR,MJ2,MJ3)
      CALL MULVER(G,C,Y(1,INS),L,MR)
      IF(INS.EQ.NS) GO TO 101
cc	CALL MULVER(A,X,Z,MR,IR,MJ3,MJ1)
      CALL MULVER(A,X(1,INS),Z,MR,IR)
      IF(M.EQ.1) GO TO 360
      DO 20 I=1,MR1
      IPR=I+IR
cxx   20 Z(I)=Z(I)+C(IPR)
      Z(I)=Z(I)+C(IPR)
   20 CONTINUE
cc  360 CALL MULVER(B,Y,C,MR,L,MJ3,MJ2)
  360 CALL MULVER(B,Y(1,INS),C,MR,L)
      CALL VECADL(C,Z,MR)
C     SUM AND SUM OF SQUARES COMPUTATION
cc  101 CALL VECADL(XS,X,IR)
cc	CALL VECADL(YS,Y,L)
  101 CALL VECADL(XS,X(1,INS),IR)
      CALL VECADL(YS,Y(1,INS),L)
      DO 30 I=1,IR
cc   30 XS2(I)=XS2(I)+X(I)**2
cxx   30 XS2(I)=XS2(I)+X(I,INS)**2
      XS2(I)=XS2(I)+X(I,INS)**2
   30 CONTINUE
      DO 31 I=1,L
cc   31 YS2(I)=YS2(I)+Y(I)**2
cxx   31 YS2(I)=YS2(I)+Y(I,INS)**2
      YS2(I)=YS2(I)+Y(I,INS)**2
   31 CONTINUE
cc	IF(INTP.EQ.1) GO TO 10
C     X,Y,W PRINT OUT
cc	WRITE(6,260) INS
cc	WRITE(6,261)
cc	CALL PRCOL2(X,W,1,IR,0)
cc	WRITE(6,262)
cc	CALL PRCOL1(Y,1,L,0)
   10 CONTINUE
C     MEAN, MEAN SQUARE AND VARIANCE COMPUTATION
      ANS=NS
      CST1=1.0D-00
      BNS=CST1/ANS
      DO 40 I=1,IR
      XMEAN(I)=BNS*XS(I)
      XS2MEA(I)=BNS*XS2(I)
cxx   40 XVAR(I)=XS2MEA(I)-XMEAN(I)**2
      XVAR(I)=XS2MEA(I)-XMEAN(I)**2
   40 CONTINUE
      DO 41 I=1,L
      YMEAN(I)=BNS*YS(I)
      YS2MEA(I)=BNS*YS2(I)
cxx   41 YVAR(I)=YS2MEA(I)-YMEAN(I)**2
      YVAR(I)=YS2MEA(I)-YMEAN(I)**2
   41 CONTINUE
cc	WRITE(6,160)
cc	WRITE(6,161)
cc	CALL PRCOL4(XMEAN,XS2,XS2MEA,XVAR,1,IR,0)
cc	WRITE(6,163)
cc	WRITE(6,164)
cc	CALL PRCOL4(YMEAN,YS2,YS2MEA,YVAR,1,L,0)
cc	CALL FLCLS2(NFL)
cc  999 CONTINUE
      RETURN
cxx    1 FORMAT(10I5)
cxx    3 FORMAT(6D12.3)
cxx   60 FORMAT(1H ,42HPROGRAM 5.5.2   OPTIMAL CONTROL SIMULATION)
cxx   61 FORMAT(1H ,17HINITIAL CONDITION)
cxx   62 FORMAT(1H ,2HN=,I5,5X,2HM=,I5,5X,3HIR=,I5,5X,2HL=,I5,5X,3HNS=,I5)
cxx   63 FORMAT(//1H ,44HFIRST IR COLUMNS OF TRANSITION MATRIX (AI'S))
cxx   64 FORMAT(//1H ,19HGAMMA MATRIX (BI'S))
cxx   65 FORMAT(//1H ,7HQ1(I,J))
cxx   66 FORMAT(//1H ,6HR(I,J))
cxx   67 FORMAT(//1H ,13HGAIN MATRIX G)
cxx  160 FORMAT(//1H ,29X,4HX(I))
cxx  161 FORMAT(1H ,4X,1HI,5X,9HMEAN OF X,5X,11HSUM OF X**2,3X,12HMEAN OF X
cxx     A**2,2X,13HVARIANCE OF X)
cxx  163 FORMAT(//1H ,29X,4HY(I))
cxx  164 FORMAT(1H ,4X,1HI,5X,9HMEAN OF Y,5X,11HSUM OF Y**2,3X,12HMEAN OF Y
cxx     A**2,2X,13HVARIANCE OF Y)
cxx  260 FORMAT(/1H ,4HINS=,I5)
cxx  261 FORMAT(1H ,4X,1HI,5X,4HX(I),10X,4HW(I))
cxx  262 FORMAT(1H ,4X,1HI,5X,4HY(I))
      END
C
      SUBROUTINE VECADL(X,Y,MM)
C     X=X+Y (X,Y: VECTORS)
cxx      IMPLICIT REAL*8(A-H,O-Z)
cxx      DIMENSION X(MM),Y(MM)
      INTEGER MM
      DOUBLE PRECISION X(MM), Y(MM)
c local
      INTEGER I
      DO 10 I=1,MM
cxx   10 X(I)=X(I)+Y(I)
      X(I)=X(I)+Y(I)
   10 CONTINUE
      RETURN
      END
