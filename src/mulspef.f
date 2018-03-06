      SUBROUTINE MULSPEF(N,K,LAGH1,LAGH3,CV,P1,P2,PS,PCH1,PCH2)
C
      INCLUDE 'timsac_f.h'
C
C     PROGRAM 5.2.2   MULTIPLE SPECTRUM
C-----------------------------------------------------------------------
C      SUBROUTINE MULSPEF(N,K,LAGH1,LAGH3,IR0,IR1,IR2,IC0,IC1,IC2,
C      SUBROUTINE CROSSP(FC,FS,P1,P2,LAGH1,A,LA1)
C      SUBROUTINE ECORSI(FS,LAGH1,FS1,LAGSHF,LA1)
C      SUBROUTINE FGERSI(G,LGP1,FS,LF1)
C      SUBROUTINE SIMCOH(P1,P2,C,S,P3,LAGH1)
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
C     THIS PROGRAM COMPUTES MULTIPLE SPECTRUM ESTIMATES FROM THE OUTPUT
C     OF PROGRAM 5.1.2 MULCOR, USING WINDOWS W1 AND W2.
C     ONLY ONE CARD OF LAGH(MAXIMUM LAG OF COVARIANCES TO BE USED FOR
C     SPECTRUM COMPUTATION) SHOULD BE ADDED ON TOP OF THE OUTPUT OF
C     PROGRAM 5.1.2 MULCOR TO FORM THE INPUT TO THIS PROGRAM.
C     IN THE CARD OUTPUT OF SPECTRUM MATRIX ON AND LOWER DIAGONAL ARE
C     REAL PARTS AND UPPER DIAGONAL ARE IMAGINARY PARTS OF ON AND LOWER
C     DIAGONAL SPECTRAL ELEMENTS.
C
cc      !DEC$ ATTRIBUTES DLLEXPORT :: MULSPEF
c      USE DFLIB
C
cxx      IMPLICIT REAL*8(A-H,O-Z)
c      DIMENSION C(501),S(501),G(501)                                    
c      DIMENSION FC(501),FS(501),P1(501),P2(501),P3(501)                 
c      DIMENSION A1(10),A2(10)                                           
c      DIMENSION P(10020)                                                
c      DIMENSION X(101,10,10)                                            
cxx      DIMENSION C(LAGH1),S(LAGH1),G(LAGH1),CV(LAGH3,K,K)
cxx      DIMENSION FC(LAGH1),FS(LAGH1)
cxx      DIMENSION P1(LAGH1,K,K),P2(LAGH1,K,K),PS(LAGH1,K)
cxx      DIMENSION PCH1(LAGH1,K,K),PCH2(LAGH1,K,K)
cxx      DIMENSION A1(2),A2(3)
cxx      DIMENSION P(LAGH1*2*K)
      INTEGER :: N, K, LAGH1, LAGH3
      REAL(8) :: CV(LAGH3,K,K), P1(LAGH1,K,K), P2(LAGH1,K,K),
     1           PS(LAGH1,K), PCH1(LAGH1,K,K), PCH2(LAGH1,K,K)
      REAL(8) :: C(LAGH1), S(LAGH1), G(LAGH1), FC(LAGH1), FS(LAGH1),
     1           A1(2), A2(3), P(LAGH1*2*K), CT5
C
C     INPUT / OUTPUT DATA FILE OPEN
c      CHARACTER(100) DFNAM
c      DFNAM='mulspe.out'
c      CALL SETWND
c      CALL FLOPN3(DFNAM,NFL)
c      IF (NFL.EQ.0) GO TO 999
C     WINDOW W1 DEFINITION
      MLA1=2
      A1(1)=0.5D-00
      A1(2)=0.25D-00
C     WINDOW W2 DEFINITION
      MLA2=3
      A2(1)=0.625D-00
      A2(2)=0.25D-00
      A2(3)=-0.0625D-00
C     LAGH SPECIFICATION
c      READ(5,1) LAGH
c      LAGH1=LAGH+1
C     READING THE OUTPUTS OF PROGRAM 5.1.2 MULCOR
c      READ(5,1) N,LAGH0,K
c      LAGH3=LAGH0+1
C     INITIAL CONDITION PRINT AND PUNCH OUT
c      WRITE(6,60)
c      WRITE(6,61)
c      WRITE(6,62) N,LAGH,K
c      WRITE(6,63)
c      CALL PRCOL1(A1,1,MLA1,1)
c      WRITE(6,64)
c      CALL PRCOL1(A2,1,MLA2,1)
c      WRITE(7,1) N,LAGH,K
C     COMPUTATION STARTS HERE WITH DISK AS EXTERNAL MEMORY.
CC      REWIND 1
c      IP=0
      CT5=0.5D-00
      LAG2=LAGH1+LAGH1
      DO 10 II=1,K
C     AUTO COVARIANCE INPUT
c      READ(5,1) IR0,IC0
c      READ(5,2) (C(I),I=1,LAGH3)
      DO 15 I=1,LAGH1
cxx   15 C(I)=CV(I,II,II)
      C(I)=CV(I,II,II)
   15 CONTINUE
      DO 20 I=1,LAGH1
cxx   20 G(I)=C(I)+C(I)
      G(I)=C(I)+C(I)
   20 CONTINUE
      G(1)=CT5*G(1)
      G(LAGH1)=CT5*G(LAGH1)
C     F-COS TRANSFORMATION
      CALL FGERCO(G,LAGH1,FC,LAGH1)
C     SPECTRUM SMOOTHING BY WINDOW W1
c      CALL AUSP(FC,P1,LAGH1,A1,MLA1)
      CALL AUSP(FC,P1(1,II,II),LAGH1,A1,MLA1)
C     SPECTRUM SMOOTHING BY WINDOW W2
c      CALL AUSP(FC,P2,LAGH1,A2,MLA2)
      CALL AUSP(FC,P2(1,II,II),LAGH1,A2,MLA2)
C     TEST STATISTICS COMPUTATION
c      CALL SIGNIF(P1,P2,P3,LAGH1,N)
      CALL SIGNIF(P1(1,II,II),P2(1,II,II),PS(1,II),LAGH1,N)
C     AUTO SPECTRUM AND TEST STATISTICS PRINT OUT
c      WRITE(6,65) IR0,IC0
c      WRITE(6,65) IR0(II),IC0(II)
c      WRITE(6,66)
c      WRITE(6,67)
c      CALL PRCOL3(P1,P2,P3,1,LAGH1,1)
C      AUTO SPECTRUM STORE
cc	LAG2=LAGH1+LAGH1
      II1=(II-1)*LAG2
      DO 21 I=1,LAGH1
      I1=II1+I
      I2=I1+LAGH1
c      P(I1)=P1(I)
c   21 P(I2)=P2(I)
      P(I1)=P1(I,II,II)
cxx   21 P(I2)=P2(I,II,II)
      P(I2)=P2(I,II,II)
   21 CONTINUE
      IF(II.EQ.1) GO TO 10
C     CROSS COVARIANCE INPUT
      IM1=II-1
      DO 11 JJ=1,IM1
c      READ(5,1) IR1,IC1
c      READ(5,2) (C(I),I=1,LAGH3)
c      READ(5,1) IR2,IC2
c      READ(5,2) (S(I),I=1,LAGH3)
      DO 25 I=1,LAGH1
      C(I)=CV(I,II,JJ)
cxx   25 S(I)=CV(I,JJ,II)
      S(I)=CV(I,JJ,II)
   25 CONTINUE
C     F-COS TRANSFORMATION
      DO 30 I=1,LAGH1
cxx   30 G(I)=C(I)+S(I)
      G(I)=C(I)+S(I)
   30 CONTINUE
      G(1)=CT5*G(1)
      G(LAGH1)=CT5*G(LAGH1)
      CALL FGERCO(G,LAGH1,FC,LAGH1)
C     F-SIN TRANSFORMATION
      DO 31 I=1,LAGH1
cxx   31 G(I)=S(I)-C(I)
      G(I)=S(I)-C(I)
   31 CONTINUE
      G(1)=CT5*G(1)
      G(LAGH1)=CT5*G(LAGH1)
      CALL FGERSI(G,LAGH1,FS,LAGH1)
C     SMOOTHING BY WINDOW W1
      ISW=1
      IWD=0
c      CALL CROSSP(FC,FS,P1,P2,LAGH1,A1,MLA1)
      CALL CROSSP(FC,FS,P1(1,II,JJ),P1(1,JJ,II),LAGH1,A1,MLA1)
C     SIMPLE COHERENCE COMPUTATION
   33 II1=(II-1)*LAG2+IWD
      JJ1=(JJ-1)*LAG2+IWD
      DO 32 I=1,LAGH1
      I1=II1+I
      I2=JJ1+I
      C(I)=P(I1)
cxx   32 S(I)=P(I2)
      S(I)=P(I2)
   32 CONTINUE
c      CALL SIMCOH(P1,P2,C,S,P3,LAGH1)
      IF (ISW.EQ.1) THEN
         CALL SIMCOH(P1(1,II,JJ),P1(1,JJ,II),C,S,PCH1(1,II,JJ),LAGH1)
      ELSE
         CALL SIMCOH(P2(1,II,JJ),P2(1,JJ,II),C,S,PCH2(1,II,JJ),LAGH1)
      END IF
C     CROSS SPECTRUM AND SIMPLE COHERENCE   PRINT OUT
c      WRITE(6,65) IR1,IC1
c      WRITE(6,65) IR1(JJ,II),IC1(JJ,II)
c      IF(ISW.NE.1) GO TO 260
c      WRITE(6,166)
c      GO TO 268
c  260 WRITE(6,266)
c  268 WRITE(6,167)
c      CALL PRCOL3(P1,P2,P3,1,LAGH1,1)
      IF(ISW.LT.0) GO TO 11
C     CROSS SPECTRUM STORE (DISK)
CC     WRITE(1) (P1(I),I=1,LAGH1),(P2(I),I=1,LAGH1)
C     SMOOTHING BY WINDOW W2
      ISW=-1
      IWD=LAGH1
c      CALL CROSSP(FC,FS,P1,P2,LAGH1,A2,MLA2)
      CALL CROSSP(FC,FS,P2(1,II,JJ),P2(1,JJ,II),LAGH1,A2,MLA2)
      GO TO 33
   11 CONTINUE
   10 CONTINUE
CC      END FILE 1
C     SPECTRUM (SMOOTHED BY WINDOW W1) PUNCH OUT
c      LC=LAGH1
c      IL=101
c      ILM1=IL-1
c      J1=0
c      J2=0
c      IB=0
c  416 J1=J2+1
c      J2=J1+ILM1
c      IF(J2.LE.LC) GO TO 413
c  412 J2=LC
c  413 CONTINUE
CC      REWIND 1
c      IP=0
c      DO 500 II=1,K
c      II1=(II-1)*LAG2
c      DO 520 I=J1,J2
c      I0=I-IB
c      I1=II1+I
c  520 X(I0,II,II)=P(I1)
c      IF(II.EQ.1) GO TO 500
c      IM1=II-1
c      DO 501 JJ=1,IM1
CC      READ(1) (P1(I),I=1,LAGH1),(P2(I),I=1,LAGH1)
c      DO 530 I=J1,J2
c      I0=I-IB
c      X(I0,II,JJ)=P1(I)
c  530 X(I0,JJ,II)=P2(I)
c  501 CONTINUE
c  500 CONTINUE
c      DO 610 I=J1,J2
c      I0=I-IB
c      DO 611 II=1,K
c  611 WRITE(7,2) (X(I0,II,JJ),JJ=1,K)
c  610 CONTINUE
c  417 IB=IB+IL
c      IF(J2.LT.LC) GO TO 416
c      CALL FLCLS3(NFL)
c  999 CONTINUE
c    1 FORMAT(10I5)
c    2 FORMAT(4D20.10)
c   60 FORMAT(1H ,13HPROGRAM 5.2.2,3X,17HMULTIPLE SPECTRUM)
c   61 FORMAT(1H ,17HINITIAL CONDITION)
c   62 FORMAT(1H ,2HN=,I5,5X,5HLAGH=,I5,5X,2HK=,I5)
c   63 FORMAT(1H ,12X,9HWINDOW W1/1H ,4X,1HI,11X,5HA1(I))
c   64 FORMAT(1H ,12X,9HWINDOW W2/1H ,4X,1HI,11X,5HA2(I))
c   65 FORMAT(//1H ,8HP(II,JJ),5X,3HII=,I5,3X,3HJJ=,I5)
c   66 FORMAT(1H ,7X,14HPOWER SPECTRUM)
c   67 FORMAT(1H ,4X,1HI,8X,8HPOWER W1,6X,8HPOWER W2,2X,12HSIGNIFICANCE)
c  166 FORMAT(1H ,7X,14HCROSS SPECTRUM,8X,2HW1)
c  266 FORMAT(1H ,7X,14HCROSS SPECTRUM,8X,2HW2)
c  167 FORMAT(1H ,4X,1HI,5X,11HCO-SPECTRUM,1X,13HQUAD-SPECTRUM,2X,16HSIMP
c     ALE COHERENCE)
      RETURN
      END SUBROUTINE
C
      SUBROUTINE CROSSP(FC,FS,P1,P2,LAGH1,A,LA1)
C     THIS SUBROUTINE COMPUTES SMOOTHED CROSS SPECTRUM.
C     FC,FS: OUTPUTS OF FGERCO AND FGERSI
C     P1,P2: REAL AND IMAGINARY PART OF SMOOTHED CROSS SPECTRUM
C     LAGH1: DIMENSION OF FC, FS AND PI (I=1,2)
C     A: SMOOTHING COEFFICIENTS
C     LA1: DIMENSION OF A (LESS THAN 11)
cxx      IMPLICIT REAL*8(A-H,O-Z)
cxx      DIMENSION FC(LAGH1),FS(LAGH1),P1(LAGH1),P2(LAGH1)
cxx      DIMENSION A(LA1)
cxx      DIMENSION FC1(521),FS1(521)
      INTEGER :: LAGH1, LA1
      REAL(8) :: FC(LAGH1), FS(LAGH1), P1(LAGH1), P2(LAGH1), A(LA1)
      REAL(8) :: FC1(521), FS1(521)
      LA=LA1-1
      LAGSHF=LAGH1+2*LA
C     FC SHIFT-RIGHT BY LA FOR END CORRECTION
      CALL ECORCO(FC,LAGH1,FC1,LAGSHF,LA1)
C     REAL PART SMOOTHING
      CALL SMOSPE(FC1,LAGSHF,A,LA1,P1,LAGH1)
C     FS SHIFT-RIGHT BY LA FOR END CORRECTION
      CALL ECORSI(FS,LAGH1,FS1,LAGSHF,LA1)
C     IMAGINARY PART SMOOTHING
      CALL SMOSPE(FS1,LAGSHF,A,LA1,P2,LAGH1)
      RETURN
      END SUBROUTINE
C
      SUBROUTINE ECORSI(FS,LAGH1,FS1,LAGSHF,LA1)
C     FS SHIFT-RIGHT BY LA FOR IMAGINARY PART END CORRECTION
cxx      IMPLICIT REAL*8(A-H,O-Z)
cxx      DIMENSION FS(LAGH1),FS1(LAGSHF)
      INTEGER :: LAGH1, LAGSHF, LA1
      REAL(8) :: FS(LAGH1), FS1(LAGSHF)
      LAGH2=LAGH1+1
      LA=LA1-1
      DO 100 I=1,LAGH1
      I1=LAGH2-I
      I2=I1+LA
cxx  100 FS1(I2)=FS(I1)
      FS1(I2)=FS(I1)
  100 CONTINUE
      LA2=LAGH1+LA
      DO 110 I=1,LA
      I1=LA1-I
      I2=LA1+I
      I3=LA2-I
      I4=LA2+I
      FS1(I1)=-FS1(I2)
cxx  110 FS1(I4)=-FS1(I3)
      FS1(I4)=-FS1(I3)
  110 CONTINUE
      RETURN
      END SUBROUTINE
C
      SUBROUTINE FGERSI(G,LGP1,FS,LF1)
C     FOURIER TRANSFORM (GOERTZEL METHOD)
C     THIS SUBROUTINE COMPUTES FOURIER TRANSFORM OF G(I),I=0,1,...,LG AT
C     FREQUENCIES K/(2*LF),K=0,1,...,LF AND RETURNS SIN TRANSFORM IN
C     FS(K).
cxx      IMPLICIT REAL*8(A-H,O-Z)
cxx      DIMENSION G(LGP1),FS(LF1)
      INTEGER :: LGP1 ,LF1
      REAL(8) :: G(LGP1), FS(LF1)
      REAL(8) :: T, PI, ALF, AK, TK, CK, SK, CK2, UM0, UM1, UM2
      LG=LGP1-1
      LF=LF1-1
C     REVERSAL OF G(I),I=1,...,LGP1 INTO G(LG3-I)   LG3=LGP1+1
      IF(LGP1.LE.1) GO TO 110
      LG3=LGP1+1
      LG4=LGP1/2
      DO 100 I=1,LG4
      I2=LG3-I
      T=G(I)
      G(I)=G(I2)
cxx  100 G(I2)=T
      G(I2)=T
  100 CONTINUE
  110 PI=3.1415926536
      ALF=LF
      T=PI/ALF
      DO 10 K=1,LF1
      AK=K-1
      TK=T*AK
      CK=DCOS(TK)
      SK=DSIN(TK)
      CK2=CK+CK
      UM2=0.0D-00
      UM1=0.0D-00
      IF(LG.EQ.0) GO TO 12
      DO 11 I=1,LG
      UM0=CK2*UM1-UM2+G(I)
      UM2=UM1
cxx   11 UM1=UM0
      UM1=UM0
   11 CONTINUE
   12 FS(K)=SK*UM1
   10 CONTINUE
      RETURN
      END SUBROUTINE
C
      SUBROUTINE SIMCOH(P1,P2,C,S,P3,LAGH1)
C     THIS SUBROUTINE COMPUTES SIMPLE COHERENCE.
cxx      IMPLICIT REAL*8(A-H,O-Z)
cxx      DIMENSION P1(LAGH1),P2(LAGH1),C(LAGH1),S(LAGH1),P3(LAGH1)
      INTEGER :: LAGH1
      REAL(8) :: P1(LAGH1), P2(LAGH1), C(LAGH1), S(LAGH1), P3(LAGH1)
      DO 10 I=1,LAGH1
cxx   10 P3(I)=(P1(I)**2+P2(I)**2)/(C(I)*S(I))
      P3(I)=(P1(I)**2+P2(I)**2)/(C(I)*S(I))
   10 CONTINUE
      RETURN
      END SUBROUTINE
