      SUBROUTINE FPEC7F(N,L,IR,IP,IP0,INW,R1,R2,FPEC,RFPEC,AIC,
     & IFPEC,OFPEC,ORFPEC,OAIC,OSD,AO)
C
      INCLUDE 'timsac_f.h'
C
C     PROGRAM 5.3.2   FPEC(AR-MODEL FITTING FOR CONTROL)
C-----------------------------------------------------------------------
C      SUBROUTINE FPEC7F(N,L,IR,IP,IP0,LAGH1,INW,MJ,MJ0,
C      SUBROUTINE RECOVA(X,LAGH1,L1,IP0)
C      SUBROUTINE SFPEC(SD,N,K,IR,MS,Z,RZ,OOZ,AIC)
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
C     THIS PROGRAM PERFORMS FPEC(AR-MODEL FITTING FOR CONTROL)
C     COMPUTATION.
C     BESIDES THE OUTPUTS OF PROGRAM 5.1.2   MULCOR, THE FOLLOWING
C     INPUTS ARE REQUIRED:
C     L: UPPER LIMIT OF MODEL ORDER M (LESS THAN 30)
C     IR: NUMBER OF CONTROLLED VARIABLES
C     IL: NUMBER OF MANINPULATED VARIABLES, IL=0 FOR MFPE COMPUTATION
C     INW(I): INDICATOR; FIRST IR INDICATE THE CONTROLLED VARIABLES
C     AND THE REST THE MANIPULATE VARIABLES WITHIN THE IP0 VARIABLES
C     IN THE OUTPUT OF PROGRAM 5.1.2   MULCOR.
C     THE OUTPUTS ARE THE PREDICTION ERROR COVARIANCE MATRIX OSD AND
C     THE SET OF COEFFICIENT MATRICES A AND B TO BE USED IN
C     PROGRAM 5.5.1   OPTIMAL CONTROLLER DESIGN.
C
cc      !DEC$ ATTRIBUTES DLLEXPORT :: FPEC7F
c      USE DFLIB
C
      IMPLICIT REAL*8(A-H,O-Z)
c      DIMENSION R1(30,10,10),A1(30,10,10),B1(30,10,10),AO(30,10,10)
c      DIMENSION SD(10,10),SE(10,10),SF(10,10),OSD(10,10)
c      DIMENSION XSD(10,10),XSF(10,10),D(10,10),E(10,10),Z1(10,10)
c      DIMENSION INW(10),C1(10,10)
      DIMENSION INW(IP),R1(L+1,IP0,IP0),R2(L+1,IP,IP)
      DIMENSION FPEC(0:L),RFPEC(0:L),AIC(0:L)
      DIMENSION OSD(IR,IR),AO(L,IR,IP)
      DIMENSION A1(L,IP,IP),B1(L,IP,IP),C1(IP0,IP0)
      DIMENSION SD(IP,IP),SE(IP,IP),SF(IP,IP)
      DIMENSION XSD(IP,IP),XSF(IP,IP),D(IP,IP),E(IP,IP),Z1(IP,IP)
C
C     INPUT / OUTPUT DATA FILE OPEN
c	CHARACTER(100) DFNAM
c	DFNAM='fpec.out'
c	CALL SETWND
c	CALL FLOPN3(DFNAM,NFL)
c	IF (NFL.EQ.0) GO TO 999
C     INITIAL CONDITION INPUT
c      READ(5,1) L,IR,IL
c      IP=IR+IL
c      READ(5,1) (INW(I),I=1,IP)
C     READING THE OUTPUTS OF PROGRAM 5.1.2 MULCOR
c      READ(5,1) N,LAGH,IP0
      L1=L+1
c      LAGH1=LAGH+1
c      CALL RECOVA(R1,LAGH1,L1,IP0,MJ0,MJ)
      DO 10 II=1,L1
      DO 20 I=1,IP0
      DO 20 J=1,IP0
   20 C1(I,J)=R1(II,I,J)
C     MATRIX REARRANGEMENT BY INW
c	  CALL REARRA(C1,INW,IP0,IP,MJ)
      CALL REARRA(C1,INW,IP0,IP)
      DO 21 I=1,IP
      DO 21 J=1,IP
   21 R2(II,I,J)=C1(I,J)
   10 CONTINUE
C     INITIAL CONDITION AND COVARIANCE PRINT OUT
c      WRITE(6,39)
c      WRITE(6,40)
c      WRITE(6,41) N,L,IR,IL
c      WRITE(6,259) (INW(I),I=1,1,IP)
c      WRITE(6,42)
c      CALL PRMAT3(R1,L1,IP,IP,1,MJ0,MJ,MJ)
C     INITIAL SD, SF, SE COMPUTATION
      DO 330 II=1,IP
      DO 330 JJ=1,IP
c      SD(II,JJ)=R1(1,II,JJ)
      SD(II,JJ)=R2(1,II,JJ)
      SF(II,JJ)=SD(II,JJ)
c      SE(II,JJ)=R1(2,II,JJ)
      SE(II,JJ)=R2(2,II,JJ)
      XSD(II,JJ)=SD(II,JJ)
  330 XSF(II,JJ)=SF(II,JJ)
C     0-TH STEP COMPUTATION
      IFPEC=0
      MS=0
C     OFPEC, ORFPEC COMPUTATION
c      CALL SFPEC(SD,N,IP,IR,MS,OFPEC,ORFPEC,OOFPEC,MJ)
	CALL SFPEC(SD,N,IP,IR,MS,FPEC(0),RFPEC(0),OOFPEC,AIC(0))
C     OFPEC, ORFPEC PRINT OUT
c      WRITE(6,600)
c      WRITE(6,264) MS,OFPEC,ORFPEC,AIC
c      OAIC=AIC
      OAIC=AIC(0)
	OFPEC=FPEC(0)
	ORFPEC=RFPEC(0)
C     ITERATION M=1 TO L
      DO 400 M=1,L
C     INVERSE OF SD, SF COMPUTATION
c      CALL INVDET(XSD,SDDET,IP,MJ)
c      CALL INVDET(XSF,SFDET,IP,MJ)
      CALL INVDET(XSD,SDDET,IP,IP)
      CALL INVDET(XSF,SFDET,IP,IP)
C     D, E, SD, SF COMPUTATION
c      CALL MULPLY(SE,XSF,D,IP,IP,IP,MJ,MJ,MJ)
c      CALL TRAMDL(SE,XSD,E,IP,IP,IP,MJ,MJ,MJ)
c      CALL TRAMDR(D,SE,Z1,IP,IP,IP,MJ,MJ,MJ)
c      CALL SUBTAL(SD,Z1,IP,IP,MJ,MJ)
c      CALL MULPLY(E,SE,Z1,IP,IP,IP,MJ,MJ,MJ)
c      CALL SUBTAL(SF,Z1,IP,IP,MJ,MJ)
      CALL MULPLY(SE,XSF,D,IP,IP,IP)
      CALL TRAMDL(SE,XSD,E,IP,IP,IP)
      CALL TRAMDR(D,SE,Z1,IP,IP,IP)
      CALL SUBTAL(SD,Z1,IP,IP)
      CALL MULPLY(E,SE,Z1,IP,IP,IP)
      CALL SUBTAL(SF,Z1,IP,IP)
      MS=M
	DO 410 II=1,IP
      DO 410 JJ=1,IP
      XSD(II,JJ)=SD(II,JJ)
  410 XSF(II,JJ)=SF(II,JJ)
C     FPEC,RFPEC COMPUTATION
c      CALL SFPEC(SD,N,IP,IR,MS,FPEC,RFPEC,OOFPEC,MJ)
      CALL SFPEC(SD,N,IP,IR,MS,FPEC(M),RFPEC(M),OOFPEC,AIC(M))
C     FPEC,RFPEC PRINT OUT
c      WRITE(6,264) MS,FPEC,RFPEC,AIC
C     FORWARD AND BACKWARD PREDICTOR COMPUTATION
c      CALL COEFAB(A1,B1,D,E,MS(M),IP,MJ0,MJ)
      CALL COEFAB(A1,B1,D,E,MS,L,IP)
C     MIN.FPEC, MIN.RFPEC COMPUTATION
c      IF(OFPEC.LE.FPEC) GO TO 440
c      OAIC=AIC
c      OFPEC=FPEC
c      ORFPEC=RFPEC
      IF(OFPEC.LE.FPEC(M)) GO TO 440
      OAIC=AIC(M)
      OFPEC=FPEC(M)
      ORFPEC=RFPEC(M)
      IFPEC=M
      DO 560 II=1,IR
      DO 560 JJ=1,IR
  560 OSD(II,JJ)=SD(II,JJ)
      DO 561 I=1,M
      DO 562 II=1,IR
      DO 562 JJ=1,IP
  562 AO(I,II,JJ)=A1(I,II,JJ)
  561 CONTINUE
  440 IF(M.EQ.L) GO TO 400
C     SE COMPUTATION
c      CALL NEWSE(A1,R1,SE,MS(M),IP,MJ0,MJ)
      CALL NEWSE(A1,R2,SE,MS,L,IP,L+1)
  400 CONTINUE
C     MIN.FPEC, MIN.RFPEC PRINT OUT
c      WRITE(6,607) OFPEC,ORFPEC,IFPEC
c      WRITE(6,1607) OAIC
c 1607 FORMAT(1H ,'MINIMUM AIC=',D12.5)
C     OSD, AO PRINT AND PUNCH OUT
c      WRITE(6,608)
c      CALL SUBMPR(OSD,IR,IR,MJ,MJ)
c  690 WRITE(7,1) N,IFPEC,IR,IL
c      DO 680 II=1,IR
c  680 WRITE(7,2) (OSD(II,JJ),JJ=1,IR)
c      IF(IFPEC.LE.0) GO TO 699
c      WRITE(6,609)
c      CALL PRMAT3(AO,IFPEC,IR,IP,0,MJ0,MJ,MJ)
c      DO 581 I=1,IFPEC
c      DO 582 II=1,IR
c  582 WRITE(7,2) (AO(I,II,JJ),JJ=1,IP)
c  581 CONTINUE
c  699 CONTINUE
c	CALL FLCLS3(NFL)
c  999 CONTINUE
c	STOP
c    1 FORMAT(10I5)
c    2 FORMAT(4D20.10)
c   39 FORMAT(1H ,50HPROGRAM 5.3.2   FPEC(AR-MODEL FITTING FOR CONTROL))
c   40 FORMAT(1H ,17HINITIAL CONDITION)
c   41 FORMAT(1H ,2HN=,I5,5X,2HL=,I5,5X,3HIR=,I5,5X,3HIL=,I5)
c   42 FORMAT(//1H ,17HCOVARIANCE MATRIX)
c  264 FORMAT(1H ,I5,2X,3D14.5)
c  259 FORMAT(/1H ,6HINW(I),5X,10I5)
c  600 FORMAT(///1H ,4X,1HI,12X,4HFPEC,9X,5HRFPEC,11X,3HAIC)
c  607 FORMAT(1H ,13HMINIMUM FPEC=,D12.5,2X,14HMINIMUM RFPEC=,D12.5,2X,14
c     AHATTAINED AT M=,I5)
c  608 FORMAT(//1H ,10X,10HOSD(II,JJ))
c  609 FORMAT(//1H ,10X,10H(A(I)B(I)))
	RETURN
      END SUBROUTINE
C
c      SUBROUTINE SFPEC(SD,N,K,IR,MS,Z,RZ,OOZ,MJ,AIC)
      SUBROUTINE SFPEC(SD,N,K,IR,MS,Z,RZ,OOZ,AIC)
C     FPEC COMPUTATION
      IMPLICIT REAL*8(A-H,O-Z)
c      COMMON /COMA/AIC
      DIMENSION SD(K,K)
      DIMENSION SD1(IR,IR)
      AN=N
      KM=K*MS
      ANP=N+1+KM
      ANM=N-1-KM
      AP=ANP/ANM
      APR=AP**IR
      CST1=1.0D-00
      DO 9 I=1,IR
      DO 9 J=1,IR
    9 SD1(I,J)=SD(I,J)
c      CALL SUBDET(SD1,SDRM,IR,MJ)
      CALL SUBDET(SD1,SDRM,IR,IR)
      Z=APR*SDRM
      ARM2=2*MS*K*IR
      AIC=AN*DLOG(SDRM)+ARM2
      IF(MS.NE.0) GO TO 10
      OOZ=CST1/Z
   10 RZ=Z*OOZ
      RETURN
      END SUBROUTINE