	SUBROUTINE FPEAUTF(L,N,SD,CXX,SSD,FPE,RFPE,D,CHI2,
     & OFPE1,OFPE2,ORFPE,MO,OSD,A,AO)
C
      INCLUDE 'timsac_f.h'
C
C     PROGRAM 5.3.1   FPE AUTO
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
C     THIS PROGRAM PERFORMS FPE(FINAL PREDICTION ERROR) COMPUTATION FOR
C     ONE-DIMENSIONAL AR-MODEL. A CARD CONTAINING THE FOLLOWING
C     INFORMATION OF L, UPPER LIMIT OF MODEL ORDER, SHOULD BE ADDED ON
C     TOP OF THE OUTPUT OF PROGRAM 5.1.1 AUTCOR TO FORM THE INPUT TO
C     THIS PROGRAM.
C     CXX(0) IS READ AS INITIAL SD.
C     THE OUTPUTS ARE THE COEFFICIENTS A(I) OF AR-PROCESS
C     X(N)=A(1)X(N-1)+...+A(M)X(N-M)+E(N)
C     AND THE VARIANCE SIGMA**2 OF E(N).
C     CHI**2 SHOWS THE SIGNIFICANCE OF PARCOR=A(M) AS A CHI-SQUARED
C     VARIABLE WITH D.F.=1.
C
cc      !DEC$ ATTRIBUTES DLLEXPORT :: FPEAUTF
c      USE DFLIB
C
      IMPLICIT REAL*8(A-H,O-Z)
c      DIMENSION CXX(501),A(501),B(501),AO(501)
      DIMENSION CXX(L),A(L,L),B(L),AO(L)
	DIMENSION SSD(L),FPE(L),RFPE(L),D(L),CHI2(L)
C
C     INPUT / OUTPUT DATA FILE OPEN
c	CHARACTER(100) DFNAM
c	CALL SETWND
c	DFNAM='fpeaut.out'
c	CALL FLOPN3(DFNAM,NFL)
c	IF (NFL.EQ.0) GO TO 999
C     L SPECIFICATION
c      READ(5,1) L
C     READING THE OUTPUT OF PROGRAM 5.1.1 AUTCOR
c      READ(5,1) N,LAGH
c      READ(5,2) SD,(CXX(I),I=1,LAGH)
C
C     COMPUTATION START
      AN=N
      NP1=N+1
      NM1=N-1
      ANP1=NP1
      ANM1=NM1
c      OFPE=(ANP1/ANM1)*SD
      OFPE1=(ANP1/ANM1)*SD
      CST1=1.0D-00
c      OOFPE=CST1/OFPE
      OOFPE=CST1/OFPE1
      ORFPE=CST1
      OSD=SD
      MO=0
	OFPE2=OFPE1
c      WRITE(6,155)
c      WRITE(6,156)
c      WRITE(6,57) N,L
c      WRITE(6,140)
c      WRITE(6,141) SD
c      CALL PRCOL1(CXX,1,L,0)
c      WRITE(6,157)
c      WRITE(6,58) OFPE
      SE=CXX(1)
	SD0=SD
C
      DO 400 M=1,L
      MP1=M+1
c      D=SE/SD
c      A(M)=D
c      D2=D*D
c      SD=(CST1-D2)*SD
      D(M)=SE/SD0
      A(M,M)=D(M)
      D2=D(M)*D(M)
      SSD(M)=(CST1-D2)*SD0
	SD0=SSD(M)
      ANP1=NP1+M
      ANM1=NM1-M
c      FPE=(ANP1/ANM1)*SD
c      RFPE=FPE*OOFPE
c      CHI2=D2*ANM1
      FPE(M)=(ANP1/ANM1)*SSD(M)
      RFPE(M)=FPE(M)*OOFPE
      CHI2(M)=D2*ANM1
      IF(M.EQ.1) GO TO 410
C     A(I) COMPUTATION
      LM=M-1
      DO 420 I=1,LM
c  420 A(I)=A(I)-D*B(I)
	A(I,M)=A(I,M-1)-D(M)*B(I)
  420 CONTINUE
c  410 DO 421 I=1,M
  410 CONTINUE
	DO 421 I=1,M
      IM=MP1-I
c  421 B(I)=A(IM)
	B(I)=A(IM,M)
  421 CONTINUE
c      WRITE(6,60) M
c      WRITE(6,61) SD,FPE,RFPE
c      WRITE(6,62) D,CHI2
c      WRITE(6,160)
c      CALL PRCOL1(A,1,M,0)
C
c      IF(OFPE.LT.FPE) GO TO 440
c      OFPE=FPE
c      ORFPE=RFPE
c      OSD=SD
      IF(OFPE2.LT.FPE(M)) GO TO 440
      OFPE2=FPE(M)
      ORFPE=RFPE(M)
      OSD=SSD(M)
      MO=M
      DO 430 I=1,M
c  430 AO(I)=A(I)
  430 AO(I)=A(I,M)
  440 IF(M.EQ.L) GO TO 400
      SE=CXX(MP1)
      DO 441 I=1,M
  441 SE=SE-B(I)*CXX(I)
  400 CONTINUE
C
c      WRITE(6,63) OFPE,MO
c      WRITE(6,64) ORFPE
c      WRITE(7,1) N,MO
c      WRITE(7,2) OSD
c      IF(MO.LE.0) GO TO 699
c      WRITE(7,2) (AO(I),I=1,MO)
c  699 CONTINUE
c	CALL FLCLS3(NFL)
c    1 FORMAT(10I5)
c    2 FORMAT(4D20.10)
c   57 FORMAT(1H ,2HN=,I5,5X,2HL=,I5)
c   58 FORMAT(1H ,5HOFPE=,D12.5)
c   60 FORMAT(1H ,2HM=,I5)
c   61 FORMAT(1H ,9HSIGMA**2=,D12.5,2X,4HFPE=,D12.5,2X,5HRFPE=,D12.5)
c   62 FORMAT(1H ,7HPARCOR=,D14.5,2X,15HCIH**2(D.F.=1)=,D12.5)
c  160 FORMAT(1H ,4X,1HI,12X,4HA(I))
c   63 FORMAT(1H ,13HMINIMUM FPE =,D12.5,2X,14HATTAINED AT M=,I5)
c   64 FORMAT(1H ,13HMINIMUM RFPE=,D12.5)
c  140 FORMAT(1H ,4X,1HI,5X,15HAUTO COVARIANCE)
c  141 FORMAT(1H ,4X,1H0,D16.5)
c  155 FORMAT(1H ,24HPROGRAM 5.3.1   FPE AUTO)
c  156 FORMAT(1H ,17HINITIAL CONDITION)
c  157 FORMAT(1H ,2HM=,4X,1H0)
	RETURN
      END SUBROUTINE
