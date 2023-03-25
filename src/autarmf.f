      SUBROUTINE AUTARM(N,LAGH01,CYY1,NEWL1,IQI1,B1,IPI1,A1,
cx     *	          NEWN,IQ,B2,IP,A2,STD,CXX2,G,SAIC,AICM,KQ,KP,TMP,IER,
     *          NEWN,IQ,B2,IP,A2,STD,CXX2,G,SAIC,AICM,KQ,KP,
     *          LMAX,MMAX,NMAX)
C
      INCLUDE 'timsac.h'
C
cc      PROGRAM AUTARM
C     PROGRAM 74.1.2. AUTOMATIC AR-MA MODEL FITTING; SCALAR CASE.
C-----------------------------------------------------------------------
C     ** DESIGNED BY H. AKAIKE, THE INSTITUTE OF STATISTICAL MATHEMATICS
C     ** PROGRAMMED BY E. ARAHATA, THE INSTITUTE OF STATISTICAL MATHEMAT
C        TOKYO
C     ** DATE OF THE LATEST REVISION: MARCH 25, 1977
C     ** THIS PROGRAM WAS ORIGINALLY PUBLISHED IN
C        "TIMSAC-74 A TIME SERIES ANALYSIS AND CONTROL PROGRAM PACKAGE(1
C        BY H. AKAIKE, E. ARAHATA AND T. OZAKI, COMPUTER SCIENCE MONOGRA
C        NO.5, MARCH 1975, THE INSTITUTE OF STATISTICAL MATHEMATICS
C     ** FOR THE BASIC THEORY SEE "CANONICAL CORRELATION ANALYSIS OF TIM
C        AND THE USE OF AN INFORMATION CRITERION" BY H. AKAIKE, IN
C        "SYSTEM IDENTIFICATION: ADVANCES AND CASE STUDIES" R. K. MEHRA
C        D. G. LAINIOTIS EDS. ACADEMIC PRESS, NEW YORK, 1976
C-----------------------------------------------------------------------
C     THIS PROGRAM PROVIDES AN AUTOMATIC AR-MA MODEL FITTING PROCEDURE.
C     MODELS WITH VARIOUS ORDERS ARE FITTED AND THE BEST CHOICE IS DETER
C     WITH THE AID OF THE STATISTICS AIC.
C     THE MAXIMUM LIKELIHOOD ESTIMATES OF THE COEFFICIENTS OF A SCALAR
C     AUTOREGRESSIVE MOVING AVERAGE MODEL Y(I)+B(1)Y(I-1)+...+B(IQ)Y(I-I
C     =X(I)+A(1)X(I-1)+...+A(IP)X(I-IP) OF A TIME SERIES Y(I)
C     ARE OBTAINED BY USING DAVIDON'S VARIANCE ALGORITHM.
C     PURE AUTOREGRESSION IS NOT ALLOWED.
C     FOR AR-MODELS USE THE INTERMEDIATE OUTPUTS OF CANARM.
C
C     THIS PROGRAM REQUIRES THE FOLLWING INPUTS:
C     (N,LAGH): N, LENGTH OF THE ORIGINAL DATA
C                  LAGH, MAXIMUM LAG OF COVARIANCE, NOT GREATER THAN 500
C     CYY(I) (I=0,LAGH): COVARIANCE SEQUENCE ... LAGH SHOULD BE LARGE EN
C                             TO KEEP THE INNOVATION VARIANCE AND
C                             GRADIENT COMPUTATION MEANINGFUL.
C     ***** WHEN N IS NOT GREATER THAN 501, PUT LAGH EQUAL TO N-1.
C     NEWL: TOTAL NUMBER OF CASES, NOT GREATER THAN 25
C     IQ: INITIAL AR ORDER
C     B(I)(I=1,IQ): INITIAL ESTIMATES OF AR-COEFFICIENTS
C     IP: INITIAL MA ORDER
C     A(I)(I=1,IP): INITIAL ESTIMATES OF MA-COEFFICIENTS
C     (IQI(I),IPI(I))(I=1,NEWL-1): (AR,MA) ORDERS TO BE FITTED SUCCESSIV
C                                       UNNECESSARY WHEN NEWL=1
C     ***** WHEN THE BEST CHOICE IS ON THE BORDER, SOME (AR,MA) ORDERS A
C              AUTOMATICALLY WITHIN THE LIMIT OF THE TOTAL NUMBER 25.
C
C     OUTPUTS: FOR EACH PAIR OF AR-MA ORDERS
C     ONE CARD WITH THE STATEMENTS OF THE PROBLEM
C     IQ: AR ORDER
C     B(I)(I=1,IQ): MAXIMUM LIKELIHOOD ESTIMATES OF AR COEFFICIENTS
C     IP: MA ORDER
C     A(I)(I=1,IP): MAXIMUM LIKELIHOOD ESTIMATES OF MA COEFFICIENTS
C     CXX0: INNOVATION VARIANCE
C
cc      PARAMETER (LMAX=500)
cc      PARAMETER (MMAX=50)
cc      PARAMETER (NMAX=25)
      PARAMETER (ICST=190)
C
cxx      IMPLICIT REAL*8(A-H,O-Z)
cc      COMMON /COM50/VD
      COMMON /COM70/ISWRO
      COMMON /COM71/IDOS
      COMMON /COM72/ISIK
cc      DIMENSION CYY(1001)
cc      DIMENSION A(50),B(50)
cc      DIMENSION OA(50),OB(50)
cc      DIMENSION X(50)
cc      DIMENSION G(50),STD(50)
cc      DIMENSION C(50)
cc      DIMENSION CC(190)
cc      DIMENSION VD(50,50)
cc      DIMENSION CN(51)
cc      DIMENSION IQI(25),IPI(25),SAIC(25)
cxx      DIMENSION CYY(LMAX*2+1)
cxx      DIMENSION A(MMAX),B(MMAX)
cxx      DIMENSION OA(MMAX),OB(MMAX)
cxx      DIMENSION X(MMAX)
cxx      DIMENSION G(MMAX,NMAX),STD(MMAX,NMAX)
cxx      DIMENSION C(MMAX)
cxx      DIMENSION CC(ICST)
cxx      DIMENSION VD(MMAX,MMAX)
cxx      DIMENSION CN(MMAX+1)
cxx      DIMENSION IQI(NMAX),IPI(NMAX),SAIC(NMAX)
c
cxx      DIMENSION CYY1(LAGH01)
cxx      DIMENSION IQI1(NEWL1),IPI1(NEWL1)
cxx      DIMENSION A1(IPI1(1)),B1(IQI1(1))
cxx      DIMENSION IP(NMAX),IQ(NMAX),IPO(NMAX),IQO(NMAX)
cxx      DIMENSION A2(MMAX,NMAX),B2(MMAX,NMAX)
cxx      DIMENSION SMAIC2(NMAX),CXX2(NMAX)
      INTEGER N, LAGH01, NEWL1, IQI1(NEWL1), IPI1(NEWL1), NEWN,
     1        IQ(NMAX), IP(NMAX), KQ, KP, LMAX, MMAX, NMAX
      DOUBLE PRECISION CYY1(LAGH01), A1(IPI1(1)), B1(IQI1(1)),
     1                 A2(MMAX,NMAX), B2(MMAX,NMAX), STD(MMAX,NMAX),
     2                 CXX2(NMAX), G(MMAX,NMAX), SAIC(NMAX), AICM
c local
      INTEGER IQI(NMAX), IPI(NMAX), IPO(NMAX), IQO(NMAX)
      DOUBLE PRECISION CYY(LMAX*2+1), A(MMAX), B(MMAX), OA(MMAX),
     1                 OB(MMAX), X(MMAX), C(MMAX), CC(ICST),
     2                 VD(MMAX,MMAX), CN(MMAX+1), SMAIC2(NMAX), CST0,
     3                 CST1, CST2, CST05, SMAIC, AN, CXX0, SUM, AIPQ,
     4                 DMAIC, CONST1, SAN
C
cx      INTEGER*1  TMP(1)
cx      CHARACTER  CNAME*80
C
C     INPUT / OUTPUT DATA FILE OPEN
cc	CHARACTER(100) DFNAM
cc	CALL SETWND
cc	DFNAM='autarm.out'
cc	CALL FLOPN3(DFNAM,NFL)
cc	IF (NFL.EQ.0) GO TO 999
c
cx      IER=0
cx      LU=3
cx      DO 7 I = 1,80
cx    7 CNAME(I:I) = ' '
cx      I = 1
cx      IFG = 1
cx      DO WHILE( (IFG.EQ.1) .AND. (I.LE.80) )
cx	   IF ( TMP(I).NE.ICHAR(' ') ) THEN
cx            CNAME(I:I) = CHAR(TMP(I))
cx            I = I+1
cx         ELSE
cx            IFG = 0
cx         END IF
cx      END DO
cx      IF ( I.GT.1 ) THEN
cx         IFG = 1
cx         OPEN (LU,FILE=CNAME,IOSTAT=IVAR)
cx         IF (IVAR .NE. 0) THEN
cxcx            WRITE(*,*) ' ***  autarm temp FILE OPEN ERROR :',CNAME,IVAR
cx            IER=IVAR
cx            IFG=0
cx         END IF
cx      ELSE
cx         IFG = 0
cx      END IF
C
      CST0=0.0D-00
      CST1=1.0D-00
      CST2=2.0D-00
      CST05=0.00005D-00
cc      DO 8 I=1,1001
      DO 8 I=1,2*LMAX+1
cxx    8 CYY(I)=CST0
      CYY(I)=CST0
    8 CONTINUE
cc	LAGH4=501
cc	LAGH1=500
      LAGH4=LMAX+1
      LAGH1=LMAX
C     AUTOCOVARIANCE LOADING
cc	READ(5,1) N,LAGH
      LAGH=LAGH01-1
      LAGH2=LAGH4+LAGH
cc	READ(5,2) (CYY(I),I=LAGH4,LAGH2)
      DO 88 I=1,LAGH01
      CYY(LMAX+I)=CYY1(I)
   88 CONTINUE
      DO 9 I=1,LAGH1
      IIM=LAGH4-I
      IIP=LAGH4+I
cxx    9 CYY(IIM)=CYY(IIP)
      CYY(IIM)=CYY(IIP)
    9 CONTINUE
      AN=N
      SMAIC=AN*DLOG(CYY(LAGH4))
cc      IQO=0
cc      IPO=0
      JQO=0
      JPO=0
cc 3200 READ(5,1) NEWL
      NEWL=NEWL1
      NEWN=1
C     INITIAL CONDITION LOADING FOR AR-MA (IQ,IP)
cc 2000 READ(5,1) IQ
      JQ=IQI1(1)
      JP=IPI1(1)
cc      IF(IQ.LE.0) GO TO 4205
      IF(JQ.LE.0) GO TO 4205
cc      READ(5,2) (B(I),I=1,IQ)
      DO 10 I=1,JQ
         B(I)=B1(I)
   10 CONTINUE
cc 4205 READ(5,1) IP
 4205 CONTINUE
cc      IF(IP.LE.0) GO TO 4204
      IF(JP.LE.0) GO TO 4204
cc      READ(5,2) (A(I),I=1,IP)
      DO 20 I = 1,JP
         A(I)=A1(I)
   20 CONTINUE
c
 4204 NEWLM1=NEWL-1
cc      IQI(1)=IQ
cc      IPI(1)=IP
      IQI(1)=JQ
      IPI(1)=JP
      IF(NEWLM1.EQ.0) GO TO 4208
      DO 4206 I=2,NEWL
cc 4206 READ(5,1) IQI(I),IPI(I)
      IQI(I)=IQI1(I)
cxx 4206 IPI(I)=IPI1(I)
      IPI(I)=IPI1(I)
 4206 CONTINUE
 4208 IDOS=3
      ISIK=0
 4207 CONTINUE
      ISWRO=0
C     INITIAL PRINT OUT
cc      WRITE (6,11111)
cx      IF (IFG.NE.0) WRITE(LU,11111)
cc	WRITE(6,1600)
cc	WRITE(6,1601)
cc	WRITE(6,1610) N,LAGH,IQ,IP
cc 4210 IPQ=IP+IQ
 4210 IPQ=JP+JQ
      AIPQ=IPQ
cc	WRITE(6,3112) ISWRO
cc 3112 FORMAT(/1H ,'ISWRO=',I5)
cc	IF(IQ.LE.0) GO TO 4215
cc	WRITE(6,1622)
cc	DO 1632 I=1,IQ
cc 1632 WRITE(6,1611) I,B(I)
cc 4215 IF(IP.LE.0) GO TO 4216
cc	WRITE(6,1621)
cc	DO 1631 I=1,IP
cc 1631 WRITE(6,1611) I,A(I)
cc 4216 CONTINUE
cx      IF (IFG.NE.0) THEN
cx	 WRITE(LU,3112) ISWRO
cx 3112	 FORMAT(/' ISWRO=',I5)
cx	 IF (JQ.GT.0) THEN
cx	      WRITE(LU,1622)
cx	    DO 1632 I=1,JQ
cx 1632	    WRITE(LU,1611) I,B(I)
cx      END IF
cx 4215	 IF (JP.GT.0) THEN
cx         WRITE(LU,1621)
cx         DO 1631 I=1,JP
cx 1631	    WRITE(LU,1611) I,A(I)
cx      END IF
cx      END IF
cc      DO 100 I=1,IP
      DO 100 I=1,JP
cxx  100 X(I)=A(I)
      X(I)=A(I)
  100 CONTINUE
cc      IF(IQ.LE.0) GO TO 420
cc	DO 110 I=1,IQ
cc	II=IP+I
      IF(JQ.LE.0) GO TO 420
      DO 110 I=1,JQ
      II=JP+I
cxx  110 X(II)=B(I)
      X(II)=B(I)
  110 CONTINUE
  420 CONTINUE
C     INNOVATION VARIANCE, GRADIENT AND HESSIAN COMPUTATION
cc	CALL SC0GRH(X,CYY,G,CN,CXX0,IP,IQ)
cxx      CALL SC0GRH(X,CYY,G(1,NEWN),CN,CXX0,JP,JQ,VD,MMAX,LMAX,ICST,
cxx     *                   IFG,LU)
      CALL SC0GRH(X,CYY,G(1,NEWN),CN,CXX0,JP,JQ,VD,MMAX,LMAX,ICST)
C     INVERSE OF HESSIAN COMPUTATION
C     COMMON SUBROUTINE CALL
cc	CALL MATINV(HESDET,IPQ)
cxx      CALL MATINV(IPQ,VD,MMAX,0,LU)
      CALL MATINV(IPQ,VD,MMAX)
C     CORRECTION TERM C(X)=V*G(X) COMPUTATION
      DO 900 I=1,IPQ
      SUM=CST0
      DO 910 J=1,IPQ
cc  910 SUM=SUM+VD(I,J)*G(J)
cxx  910 SUM=SUM+VD(I,J)*G(J,NEWN)
      SUM=SUM+VD(I,J)*G(J,NEWN)
  910 CONTINUE
cxx  900 C(I)=SUM
      C(I)=SUM
  900 CONTINUE
C     DAVODON'S PROCEDURE; MINIMIZATION OF INNOVATION VARIANCE
cc	CALL SDAV1(X,CYY,CXX0,G,C,IP,IQ,N)
cxx      CALL SDAV1(X,CYY,CXX0,G(1,NEWN),C,JP,JQ,N,VD,MMAX,LMAX,ICST,
cxx     *                IFG,LU)
      CALL SDAV1(X,CYY,CXX0,G(1,NEWN),C,JP,JQ,N,VD,MMAX,LMAX,ICST)
      IF(ISWRO.LE.0) GO TO 940
      ISFIN=0
      IF(ISWRO.GE.10) GO TO 940
cc	DO 902 I=1,IP
cc	IF(DABS(A(I)-X(I)).GE.CST05) GO TO 904
      DO 902 I=1,JP
      IF(DABS(A(I)-X(I)).GE.CST05) GO TO 904
  902 CONTINUE
cc	IF(IQ.LE.0) GO TO 919
cc	DO 903 I=1,IQ
cc	II=IP+I
cc	IF(DABS(B(I)-X(II)).GE.CST05) GO TO 904
      IF(JQ.LE.0) GO TO 919
      DO 903 I=1,JQ
      II=JP+I
      IF(DABS(B(I)-X(II)).GE.CST05) GO TO 904
  903 CONTINUE
      GO TO 919
  904 ISFIN=1
cc  919 DO 920 I=1,IP
cc  920 A(I)=X(I)
cc	IF(IQ.LE.0) GO TO 925
cc	DO 930 I=1,IQ
cc	II=IP+I
cc  930 B(I)=X(II)
  919 DO 920 I=1,JP
cxx  920 A(I)=X(I)
      A(I)=X(I)
  920 CONTINUE
      IF(JQ.LE.0) GO TO 925
      DO 930 I=1,JQ
      II=JP+I
cxx  930 B(I)=X(II)
      B(I)=X(II)
  930 CONTINUE
  925 CONTINUE
      IF(ISFIN.EQ.0) GO TO 940
cc	WRITE(6,926)
cx      IF (IFG.NE.0) WRITE(LU,926)
cxx  926 FORMAT(/1H ,'HESSIAN RESET')
      GO TO 4210
  940 CONTINUE
cc	WRITE(6,1008) CXX0
cx      IF (IFG.NE.0) WRITE(LU,1008) CXX0
C     HESSIAN COMPUTATION
      ISWRO=1
cc	CALL SC0GRH(X,CYY,G,CN,CXX0,IP,IQ)
cxx      CALL SC0GRH(X,CYY,G(1,NEWN),CN,CXX0,JP,JQ,VD,MMAX,LMAX,ICST,
cxx     *                   IFG,LU)
      CALL SC0GRH(X,CYY,G(1,NEWN),CN,CXX0,JP,JQ,VD,MMAX,LMAX,ICST)
cc	WRITE(6,1008) CXX0
cx      IF (IFG.NE.0) WRITE(LU,1008) CXX0
C     INVERSE OF HESSIAN COMPUTATION
C     COMMON SUBROUTINE CALL
cc	CALL MATINV(HESD2,IPQ)
cxx      CALL MATINV(IPQ,VD,MMAX,0,LU)
      CALL MATINV(IPQ,VD,MMAX)
C     INVERSE OF HESSIAN PRINT OUT
cc	WRITE(6,3000)
cc	DO 3100 I=1,IPQ
cc 3100 WRITE(6,3110) I,(VD(I,J),J=1,IPQ)
cx      IF (IFG.NE.0) THEN
cx	 WRITE(LU,3000)
cx	 DO 3100 I=1,IPQ
cx 3100	 WRITE(LU,3110) I,(VD(I,J),J=1,IPQ)
cx      END IF
C     PARAMETER VARIANCE MATRIX COMPUTATION
      AN=N
      CONST1=CXX0/AN
cxx      DO 6000 I=1,IPQ
      DO 6001 I=1,IPQ
      DO 6000 J=1,IPQ
cxx 6000 VD(I,J)=CONST1*VD(I,J)
      VD(I,J)=CONST1*VD(I,J)
 6000 CONTINUE
 6001 CONTINUE
cc	WRITE(6,6100)
cc	DO 6200 I=1,IPQ
cc 6200 WRITE(6,3110) I,(VD(I,J),J=1,IPQ)
cx      IF (IFG.NE.0) THEN
cx	 WRITE(LU,6100)
cx	 DO 6200 I=1,IPQ
cx 6200	 WRITE(LU,3110) I,(VD(I,J),J=1,IPQ)
cx      END IF
      DO 6400 I=1,IPQ
      IF(VD(I,I).LT.CST0) VD(I,I)=CST0
cc 6400 STD(I)=DSQRT(VD(I,I))
cxx 6400 STD(I,NEWN)=DSQRT(VD(I,I))
      STD(I,NEWN)=DSQRT(VD(I,I))
 6400 CONTINUE
C     CN(I)=CXX(I)/CXX(0) I=1,50
cc	WRITE(6,8000)
      SAN=CST2/DSQRT(AN)
cc	WRITE(6,7999) SAN
cc	WRITE(6,7998) (CN(I),I=2,51)
cx      IF (IFG.NE.0) THEN
cx	 WRITE(LU,8000) SAN
cx	 WRITE(LU,7998) (CN(I),I=2,MMAX+1)
cx      END IF
cc	DO 800 I=1,IP
cc  800 A(I)=X(I)
cc	IF(IQ.LE.0) GO TO 820
cc	DO 810 I=1,IQ
cc	II=IP+I
cc  810 B(I)=X(II)
      DO 800 I=1,JP
cxx  800 A(I)=X(I)
      A(I)=X(I)
  800 CONTINUE
      IF(JQ.LE.0) GO TO 820
      DO 810 I=1,JQ
      II=JP+I
cxx  810 B(I)=X(II)
      B(I)=X(II)
  810 CONTINUE
  820 CONTINUE
cc	WRITE(6,7910)
cc	WRITE(6,1014) NEWN
cx      IF (IFG.NE.0) WRITE(LU,1014) NEWN
cx 1014 FORMAT(//1H ,'CASE NO.',I5)
cc	IF(IQ.LE.0) GO TO 4290
cc	IF(IQ(NEWN).LE.0) GO TO 4290
      IF(JQ.LE.0) GO TO 4291
cc	WRITE(6,862)
cc	DO 863 I=1,IQ
cc	II=IP+I
cc  863 WRITE(6,864) I,B(I),STD(II)
      DO 863 I=1,JQ
      II=JP+I
  863 CONTINUE
C     INVERSE OF AR(B) COMPUTATION
      IG=0
cc	CALL INVERS(B,IQ,A,0,CC,IB,IG)
cx      CALL INVERS(B,JQ,A,0,CC,IB,IG,IFG,LU)
cxx      CALL INVERS(B,JQ,A,0,CC,IB,ICST,IG,IFG,LU)
      CALL INVERS(B,JQ,A,0,CC,IB,ICST,IG)
cc	WRITE(6,4289) IB
cx      IF (IFG.NE.0) WRITE(LU,4289) IB
cc 4290 IF(IP.LE.0) GO TO 4291
cc	WRITE(6,865)
cc	DO 866 I=1,IP
cc  866 WRITE(6,864) I,A(I),STD(I)
 4291 CONTINUE
cc	WRITE(6,1008) CXX0
      SAIC(NEWN)=AN*DLOG(CXX0)+CST2*AIPQ
cc	WRITE(6,1001) SAIC(NEWN)
cc 1001 FORMAT(/1H ,'AIC=N*LOG(CXX0)+2.0*(IQ+IP)=',D12.5)
cc	WRITE(7,4) NEWN,IQ,IP,SAIC(NEWN)
cc    4 FORMAT(/'CASE NO.',I2,1X,'AR',I2,1X,'MA',I2,2X,'AIC=',D12.5)
cc	WRITE(7,1) IQ
cc	IF(IQ.LE.0) GO TO 4292
cc	WRITE(7,2) (B(I),I=1,IQ)
      IF(JQ.LE.0) GO TO 4292
      IQ(NEWN)=JQ
      DO 4200 I=1,JQ
         B2(I,NEWN)=B(I)
 4200 CONTINUE
cc 4292 WRITE(7,1) IP
 4292 CONTINUE
cc	IF(IP.LE.0) GO TO 1000
cc	WRITE(7,2) (A(I),I=1,IP)
      IF(JP.LE.0) GO TO 1000
      IP (NEWN)=JP
      DO 4201 I=1,JP
         A2(I,NEWN)=A(I)
 4201 CONTINUE
cc 1000 WRITE(7,2) CXX0
 1000 CXX2(NEWN)=CXX0
C     FINAL GRADIENT PRINT OUT
cc	WRITE(6,1010)
cc 1010 FORMAT(/1H ,'FINAL GRADIENT')
cc	WRITE(6,7998) (G(I),I=1,IPQ)
      DMAIC=SMAIC-SAIC(NEWN)
      IF(DMAIC.LT.CST0) GO TO 1013
      SMAIC=SAIC(NEWN)
cc	IQO=IQ
cc	IPO=IP
cc      DO 1011 I=1,IQ
cc 1011 OB(I)=B(I)
cc      DO 1012 I=1,IP
cc 1012 OA(I)=A(I)
      JQO=JQ
      JPO=JP
      DO 1011 I=1,JQ
cxx 1011 OB(I)=B(I)
      OB(I)=B(I)
 1011 CONTINUE
      DO 1012 I=1,JP
cxx 1012 OA(I)=A(I)
      OA(I)=A(I)
 1012 CONTINUE
 1013 CONTINUE
      IF(NEWN.GE.NEWL)  GO TO 2100
cc 2085 IQ1=IQ+1
cc      IP1=IP+1
 2085 JQ1=JQ+1
      JP1=JP+1
cc      DO 2090 I=IQ1,50
cxx      DO 2090 I=JQ1,MMAX
cxx 2090 B(I)=CST0
      B(JQ1:MMAX)=CST0
cc      DO 2095 I=IP1,50
cxx      DO 2095 I=JP1,MMAX
cxx 2095 A(I)=CST0
      A(JP1:MMAX)=CST0
      NEWN=NEWN+1
      IDOS=0
cc 2096 IQ=IQI(NEWN)
cc      IF(IP.LE.IPI(NEWN)) GO TO 2976
      JQ=IQI(NEWN)
      IF(JP.LE.IPI(NEWN)) GO TO 2976
      IDOS=3
cc 2976 IP=IPI(NEWN)
 2976 JP=IPI(NEWN)
      GO TO 4207
cc 2100 WRITE(6,3) SMAIC,IQO,IPO
cc    3 FORMAT(/1H ,'MINUMUM AIC =',D12.5,' ATTAINED AT THE BEST CHOICE
cc     AAR=',I5,'  MA=',I5)
 2100 CONTINUE
      SMAIC2(NEWN)=SMAIC
      IQO(NEWN)=JQO
      IPO(NEWN)=JPO
C     BORDER CHECK
cc 2111 IF(NEWL.GE.25) GO TO 2120
cc 2112 IQM1=IQO-1
cc      IPM1=IPO-1
      IF(NEWL.GE.NMAX) GO TO 2120
      IQM1=JQO-1
      IPM1=JPO-1
      IDOS=0
      IDO=-1
cc      IP=IPO+1
cc      IQ=IQO+1
      JP=JPO+1
      JQ=JQO+1
      GO TO 2150
 2109 IDO=0
cc      IP=MAX0(IPM1,1)
cc      IQ=MAX0(IQM1,0)
cc      IF(IPO.LE.IP) GO TO 2150
      JP=MAX0(IPM1,1)
      JQ=MAX0(IQM1,0)
      IF(JPO.LE.JP) GO TO 2150
      IDOS=3
      GO TO 2150
 2110 IDO=1
cc      IP=IPO
cc      IQ=MAX0(IQM1,0)
      JP=JPO
      JQ=MAX0(IQM1,0)
      GO TO 2150
 2113 IDO=2
cc      IQ=IQO+1
      JQ=JQO+1
      GO TO 2150
 2114 IDO=3
cc      IQ=IQO
cc      IP=MAX0(IPM1,1)
cc      IF(IPO.LE.IP) GO TO 2150
      JQ=JQO
      JP=MAX0(IPM1,1)
      IF(JPO.LE.JP) GO TO 2150
      IDOS=3
      GO TO 2150
 2115 IDO=4
cc      IP=IPO+1
      JP=JPO+1
      GO TO 2150
 2116 IDO=5
cc      IQ=MAX0(IQM1,0)
      JQ=MAX0(IQM1,0)
      GO TO 2150
 2117 IDO=6
cc      IQ=IQO+1
cc      IP=MAX0(IPM1,1)
cc      IF(IPO.LE.IP) GO TO 2150
      JQ=JQO+1
      JP=MAX0(IPM1,1)
      IF(JPO.LE.JP) GO TO 2150
      IDOS=3
 2150 DO 2151 I=1,NEWL
cc      IDE=IABS(IQI(I)-IQ)+IABS(IPI(I)-IP)
      IDE=IABS(IQI(I)-JQ)+IABS(IPI(I)-JP)
      IF(IDE.EQ.0) GO TO 2152
 2151 CONTINUE
      GO TO 2154
 2152 IDOS=0
      IF(IDO.EQ.-1) GO TO 2109
      IF(IDO.EQ.0) GO TO 2110
      IF(IDO.EQ.1) GO TO 2113
      IF(IDO.EQ.2) GO TO 2114
      IF(IDO.EQ.3) GO TO 2115
      IF(IDO.EQ.4) GO TO 2116
      IF(IDO.EQ.5) GO TO 2117
cc      WRITE(6,2153)
cx      IF (IFG.NE.0) WRITE(LU,2153)
cxx 2153 FORMAT(//1H ,'BORDER CHECK COMPLETED')
      GO TO 2120
 2154 NEWL=NEWL+1
cc      IQI(NEWL)=IQ
cc      IPI(NEWL)=IP
cc      IQ=IQO
cc      IP=IPO
      IQI(NEWL)=JQ
      IPI(NEWL)=JP
      JQ=JQO
      JP=JPO
cc      DO 2155 I=1,IQO
      DO 2155 I=1,JQO
cxx 2155 B(I)=OB(I)
      B(I)=OB(I)
 2155 CONTINUE
cc      DO 2156 I=1,IPO
      DO 2156 I=1,JPO
cxx 2156 A(I)=OA(I)
      A(I)=OA(I)
 2156 CONTINUE
      GO TO 2085
cc 2120 CALL FLCLS3(NFL)
 2120 CONTINUE
cxx  999 CONTINUE
      CONTINUE
      AICM=SMAIC2(NEWN)
      KQ=IQO(NEWN)
      KP=IPO(NEWN)
cx      IF (IFG.NE.0) CLOSE(LU)
      RETURN
cxx    1 FORMAT(16I5)
cxx    2 FORMAT(4D20.10)
cxx  862 FORMAT(/1H ,4X,1HI,13X,5HAR(I),1X,'STANDARD DEVIATION')
cxx  864 FORMAT(1H ,I5,2D17.5)
cxx  865 FORMAT(/1H ,4X,1HI,13X,5HMA(I),1X,'STANDARD DEVIATION')
cxx 1600 FORMAT(/1H ,'AUTOMATIC AR-MA MODEL FITTING; SCALAR CASE')
cxx 1601 FORMAT(/1H ,'DAVIDON''S (MINIMIZATION) PROCEDURE')
cxx 1610 FORMAT(/1H ,'INITIAL CONDITION / N=',I5,',LAGH=',I5,',AR-ORDER=',
cxx     A I5,',MA-ORDER=',I5)
cxx 1622 FORMAT(1H ,4X,1HI,12X,5HAR(I))
cxx 1611 FORMAT(1H ,I5,D17.5)
cxx 1621 FORMAT(/1H ,4X,1HI,12X,5HMA(I))
cxx 1008 FORMAT(/1H ,'CXX0=',D12.5)
cxx 3000 FORMAT(/1H ,'INVERSE OF HESSIAN')
cxx 3110 FORMAT(/1H ,I5,4X,10D12.5,/(1H ,9X,10D12.5))
cxx 6100 FORMAT(/1H ,'PARAMETER VARIANCE MATRIX ESTIMATE')
cc 8000 FORMAT(/1H ,'NORMALIZED AUTOCOVARIANCE OF INNOVATION')
cc 7999 FORMAT(1H ,43X,'(2*(INVERSE OF SQUARE ROOT OF N)=',D12.5,')')
cxx 8000 FORMAT(/1H ,'NORMALIZED AUTOCOVARIANCE OF INNOVATION',
cxx     *'  (2*(INVERSE OF SQUARE ROOT OF N)=',D12.5,')')
cxx 7998 FORMAT(1H ,9X,10D12.5/(1H ,9X,10D12.5))
cxx 7910 FORMAT(/1H )
cxx 4289 FORMAT(1H ,7X,'ORDER OF THE INVERSE OF AR=',I5)
cxx11111 FORMAT(//1H ,'PROGRAM 74.1.2. AUTARM')
      END
C
C
cc	SUBROUTINE SC0GRH(X,CYY,G,CN,CXX0,IP,IQ)
cxx      SUBROUTINE SC0GRH(X,CYY,G,CN,CXX0,IP,IQ,AL,MM,LL,ICST,IFG,LU)
      SUBROUTINE SC0GRH(X,CYY,G,CN,CXX0,IP,IQ,AL,MM,LL,ICST)
C     THIS SUBROUTINE COMPUTES CXX0,GRADIENT AND HESSIAN.
cxx      IMPLICIT REAL*8(A-H,O-Z)
cc      COMMON /COM50/AL
      COMMON /COM70/ISWRO
      COMMON /COM71/IDOS
      COMMON /COM72/ISIK
cc      DIMENSION X(50),A(50),B(50),AI(190),G(50)
cc      DIMENSION CYY(1001)
cc      DIMENSION CN(51)
cxx      DIMENSION X(IP+IQ),A(IP),B(IQ),AI(ICST),G(IP+IQ)
cxx      DIMENSION CYY(LL*2+1),CN(MM+1)
cxx      DIMENSION AL(MM,MM)
c
cc      DIMENSION A2(100),A2B(190),AIB(190),Y(1001)
cc      DIMENSION CXX(1001),CXY(1001),CUU(1001)
cc      DIMENSION CYX(1001),CUX(1001)
cc      DIMENSION CUZ(1001),CUY(1001)
cc      DIMENSION CYU(1001),CYZ(1001),CZX(1001)
cc      DIMENSION CZY(1001),CZZ(1001)
cc      DIMENSION AL(50,50)
cxx      DIMENSION A2(IP*2),A2B(ICST),AIB(ICST),Y(LL*2+1)
cxx      DIMENSION CXX(LL*2+1)
cxx      DIMENSION CXY(LL*2+1)
cxx      DIMENSION CYX(LL*2+1)
cxx      DIMENSION CUZ(LL*2+1)
cxx      DIMENSION CYU(LL*2+1)
cxx      DIMENSION CZY(LL*2+1)
cc      EQUIVALENCE (CXX(1),CXY(1),CUU(1))
cc      EQUIVALENCE (CYX(1),CUX(1))
cc      EQUIVALENCE (CUZ(1),CUY(1))
cc      EQUIVALENCE (CZY(1),CZZ(1))
cc      EQUIVALENCE (CYU(1),CYZ(1),CZX(1))
c
      INTEGER IP, IQ, MM, LL, ICST
      DOUBLE PRECISION X(IP+IQ), CYY(LL*2+1), G(IP+I Q), CN(MM+1),
     1                 CXX0, AL(MM,MM)
c local
      DOUBLE PRECISION A(IP), B(IQ), AI(ICST), A2(IP*2), A2B(ICST),
     1                 AIB(ICST), Y(LL*2+1), CXX(LL*2+1), CXY(LL*2+1),
     2                 CYX(LL*2+1), CUZ(LL*2+1), CYU(LL*2+1),
     3                 CZY(LL*2+1), CST0, CST1, DSR2, CAI1
c
      CST0=0.0D-00
      CST1=1.0D-00
cc      IORIG=501
      IORIG=LL+1
      DSR2=0.95D-00
      DO 100 I=1,IP
cxx  100 A(I)=X(I)
      A(I)=X(I)
  100 CONTINUE
      IF(IQ.LE.0) GO TO 420
      DO 110 I=1,IQ
      II=IP+I
cxx  110 B(I)=X(II)
      B(I)=X(II)
  110 CONTINUE
  420 CONTINUE
      IG=1
      IF(ISIK.NE.0) GO TO 25
cxx   24 IG=0
      IG=0
C     ADJUSTMENT FOR FEASIBLE INITIAL
C     INVERSE OF A(I) COMPUTATION
cc   25 CALL INVERS(A,IP,B,0,AI,IA,IG)
cx   25 CALL INVERS(A,IP,B,0,AI,IA,IG,IFG,LU)
cxx   25 CALL INVERS(A,IP,B,0,AI,IA,ICST,IG,IFG,LU)
   25 CALL INVERS(A,IP,B,0,AI,IA,ICST,IG)
      IF(ISIK.EQ.0) GO TO 26
      IF(ISWRO.NE.0) GO TO 1900
      IF(IDOS.NE.3) GO TO 1900
   26 IF(IG.EQ.0) GO TO 1900
      CAI1=CST1
      DO 1941 I=1,IP
      CAI1=CAI1*DSR2
cxx 1941 A(I)=A(I)*CAI1
      A(I)=A(I)*CAI1
 1941 CONTINUE
cc      WRITE(6,1940)
cx      IF (IFG.NE.0) WRITE(LU,1940)
cxx 1940 FORMAT(1H ,'NON-INVERTIBLE MA PART')
      IF(ISIK.NE.0) GO TO 1899
      IG=0
 1899 GO TO 25
 1900 IF(IA.NE.0) GO TO 1901
      IA=1
      AI(1)=CST0
 1901 CONTINUE
C
      DO 2100 I=1,IP
cxx 2100 X(I)=A(I)
      X(I)=A(I)
 2100 CONTINUE
      IF(IQ.LE.0) GO TO 2420
      DO 2110 I=1,IQ
      II=IP+I
cxx 2110 X(II)=B(I)
      X(II)=B(I)
 2110 CONTINUE
 2420 CONTINUE
      ISIK=1
      IPM1=IP-1
      IQM1=IQ-1
C     AIB=(INVERSE OF A)*B
      IG=0
cc      CALL INVERS(A,IP,B,IQ,AIB,IAIB,IG)
cx      CALL INVERS(A,IP,B,IQ,AIB,IAIB,IG,IFG,LU)
cxx      CALL INVERS(A,IP,B,IQ,AIB,IAIB,ICST,IG,IFG,LU)
      CALL INVERS(A,IP,B,IQ,AIB,IAIB,ICST,IG)
C     A2B=(INVERSE OF A*A)*B
C     A2=A*A
      Y(IORIG)=CST1
      DO 502 I=1,IP
      IJ=IORIG+I
      IK=IORIG-I
      A2B(I)=A(I)
      Y(IK)=A(I)
cxx  502 Y(IJ)=CST0
      Y(IJ)=CST0
  502 CONTINUE
      IK=IORIG-IP
      DO 503 I=1,IP
      IK=IK-1
cxx  503 Y(IK)=CST0
      Y(IK)=CST0
  503 CONTINUE
      L2=-IP-IP
      M2=-1
cc      CALL SCONVL(Y,A2B,Y,IP,L2,M2)
      CALL SCONVL(Y,A2B,Y,IP,L2,M2,LL)
      I2P=IP+IP
      DO 504 I=1,I2P
      IJ=IORIG-I
cxx  504 A2(I)=Y(IJ)
      A2(I)=Y(IJ)
  504 CONTINUE
      IG=1
cc      CALL INVERS(A2,I2P,B,IQ,A2B,IA2B,IG)
cx      CALL INVERS(A2,I2P,B,IQ,A2B,IA2B,IG,IFG,LU)
cxx      CALL INVERS(A2,I2P,B,IQ,A2B,IA2B,ICST,IG,IFG,LU)
      CALL INVERS(A2,I2P,B,IQ,A2B,IA2B,ICST,IG)
      LXX=0
cc      MXX=50
      MXX=MM
      LXY=0
      MXY=MXX+IAIB
      LYX=-MXY
      MYX=-LXY
      LZX=-IQ
      MZX=-1
      LZY=LZX
      MZY=MZX+IAIB
      LZZ=0
      MZZ=IQM1
      LZYZ=LZZ
      MZYZ=MZZ+IA
      LZY=MIN0(LZY,LZYZ)
      MZY=MAX0(MZY,MZYZ)
      LYZ=-MZY
      MYZ=-LZY
      LUX=-IP
      MUX=-1
      LUY=-IP
      MUY=MUX+IAIB
      LUU=0
      MUU=IPM1
      LUYU=0
      MUYU=MUU+IA2B
      LUY=MIN0(LUY,LUYU)
      MUY=MAX0(MUY,MUYU)
      LUZ=-IPM1
      MUZ=IQM1
      LUYZ=LUZ
      MUYZ=MUZ+IA
      LUY=MIN0(LUY,LUYZ)
      MUY=MAX0(MUY,MUYZ)
      LYU=-MUY
      MYU=-LUY
cc      WRITE(6,2502)
cc 2500 WRITE(6,3000) LXX,LXY,LZX,LZY,LUX,LUU,LUZ,LUY,LUYU,LZYZ,IA,IAIB,
cc     AIA2B
cc 2501 WRITE(6,3000) MXX,MXY,MZX,MZY,MUX,MUU,MUZ,MUY,MUYU,MZYZ
cx      IF (IFG.NE.0) THEN
cx      WRITE(LU,2502)
cx 2500	 WRITE(LU,3000) LXX,LXY,LZX,LZY,LUX,LUU,LUZ,LUY,LUYU,LZYZ,
cx     *	 IA,IAIB,IA2B
cx 2501	 WRITE(LU,3000) MXX,MXY,MZX,MZY,MUX,MUU,MUZ,MUY,MUYU,MZYZ
cx      END IF
C     CXX0 COMPUTATION
C     CYX=CYY*AIB'
cc	CALL SCONVL(CYY,AIB,CYX,IAIB,LYX,MYX)
      CALL SCONVL(CYY,AIB,CYX,IAIB,LYX,MYX,LL)
C
cc      CALL TURN(CYX,CXY,LYX,MYX)
      DO 505 I=1,LL*2+1
         CXY(I)=CXX(I)
  505 CONTINUE
      CALL TURN(CYX,CXY,LYX,MYX,LL)
C     CXX=CXY*AIB'
cc      CALL SCONVL(CXY,AIB,CXX,IAIB,LXX,MXX)
      CALL SCONVL(CXY,AIB,CXX,IAIB,LXX,MXX,LL)
      DO 506 I=1,LL*2+1
         CXX(I)=CXY(I)
  506 CONTINUE
      CXX0=CXX(IORIG)
      IST=IORIG+LXX
      IEN=IORIG+MXX
      IJ=0
      DO 510 I=IST,IEN
      IJ=IJ+1
cxx  510 CN(IJ)=CXX(I)/CXX0
      CN(IJ)=CXX(I)/CXX0
  510 CONTINUE
C     GA COMPUTATION
C     CYU=CYY*A2B'
cc      CALL SCONVL(CYY,A2B,CYU,IA2B,LYU,MYU)
cc      CALL TURN(CYU,CUY,LYU,MYU)
cc      CALL SCONVL(CUY,AIB,CUX,IAIB,LUX,MUX)
      CALL SCONVL(CYY,A2B,CYU,IA2B,LYU,MYU,LL)
      CALL TURN(CYU,CUZ,LYU,MYU,LL)
      CALL SCONVL(CUZ,AIB,CYX,IAIB,LUX,MUX,LL)
C     HAA
cc      CALL SCONVL(CUY,A2B,CUU,IA2B,LUU,MUU)
      CALL SCONVL(CUZ,A2B,CXX,IA2B,LUU,MUU,LL)
      IF(IQ.EQ.0) GO TO 550
C     HAB
cc      CALL SCONVL(CUY,AI,CUZ,IA,LUZ,MUZ)
      CALL SCONVL(CUZ,AI,CUZ,IA,LUZ,MUZ,LL)
C     GB COMPUTATION
C     CYZ=CYY*AI'
C     IF(IQ.EQ.0) GO TO 550
cc      CALL SCONVL(CYY,AI,CYZ,IA,LYZ,MYZ)
cc      CALL TURN(CYZ,CZY,LYZ,MYZ)
cc      CALL SCONVL(CZY,AIB,CZX,IAIB,LZX,MZX)
      CALL SCONVL(CYY,AI,CYU,IA,LYZ,MYZ,LL)
      CALL TURN(CYU,CZY,LYZ,MYZ,LL)
      CALL SCONVL(CZY,AIB,CYU,IAIB,LZX,MZX,LL)
C     HBB
cc      CALL SCONVL(CZY,AI,CZZ,IA,LZZ,MZZ)
      CALL SCONVL(CZY,AI,CZY,IA,LZZ,MZZ,LL)
C     HESSIAN ARRANGEMENT FOR U
  550 CONTINUE
      DO 211 I=1,IP
      DO 212 J=1,I
      IJ1=IORIG+I-J
cc      AL(I,J)=CUU(IJ1)
      AL(I,J)=CXX(IJ1)
cxx  212 AL(J,I)=AL(I,J)
      AL(J,I)=AL(I,J)
  212 CONTINUE
  211 CONTINUE
      IF(IQ.LE.0) GO TO 4220
C     HESSIAN ARRANGEMENT FOR V
      DO 231 I=1,IQ
      II=IP+I
      DO 232 J=1,I
      JJ=IP+J
      IJ1=IORIG+I-J
cc      AL(II,JJ)=CZZ(IJ1)
      AL(II,JJ)=CZY(IJ1)
cxx  232 AL(JJ,II)=AL(II,JJ)
      AL(JJ,II)=AL(II,JJ)
  232 CONTINUE
  231 CONTINUE
C     HESSIAN ARRANGEMENT FOR -W AND -W'
      DO 251 I=1,IQ
      II=IP+I
      DO 252 J=1,IP
      IJ1=IORIG+I-J
      AL(II,J)=-CUZ(IJ1)
cxx  252 AL(J,II)=AL(II,J)
      AL(J,II)=AL(II,J)
  252 CONTINUE
  251 CONTINUE
 4220 CONTINUE
C     GRADIENT ARRANGEMENT
      DO 280 I=1,IP
      I1=IORIG-I
cc  280 G(I)=-CUX(I1)
cxx  280 G(I)=-CYX(I1)
      G(I)=-CYX(I1)
  280 CONTINUE
      IF(IQ.LE.0) GO TO 4230
      DO 281 I=1,IQ
      II=IP+I
      I1=IORIG-I
cc  281 G(II)=CZX(I1)
cxx  281 G(II)=CYU(I1)
      G(II)=CYU(I1)
  281 CONTINUE
      IDOS=0
 4230 RETURN
cxx 2502 FORMAT(/1H ,'PARAMETER PRINT OUT AT THE STATEMENT NUMBER',
cxx     A' 2500-2501 OF SC0GRH')
cxx 3000 FORMAT(1H ,16I5)
      END
C
cxx      SUBROUTINE SC0GR1(X,CYY,G,CXX0,IP,IQ,IG,LL,ICST,IFG,LU)
      SUBROUTINE SC0GR1(X,CYY,G,CXX0,IP,IQ,IG,LL,ICST)
C     THIS SUBROUTINE COMPUTES CXX0 AND GRADIENT.
cxx      IMPLICIT REAL*8(A-H,O-Z)
cc      DIMENSION X(50),A(50),B(50),AI(190),G(50)
cc      DIMENSION A2(100),A2B(190),AIB(190),Y(1001)
cc      DIMENSION CYY(1001)
cc      DIMENSION CXX(1001),CXY(1001)
cc      DIMENSION CYX(1001),CUX(1001)
cc      DIMENSION CUY(1001)
cc      DIMENSION CYU(1001),CYZ(1001),CZX(1001)
cc      DIMENSION CZY(1001)
cxx      DIMENSION X(IP+IQ),A(IP),B(IQ),AI(ICST),G(IP+IQ)
cxx      DIMENSION A2(IP*2),A2B(ICST),AIB(ICST),Y(LL*2+1)
cxx      DIMENSION CYY(LL*2+1),CXX(LL*2+1),CYX(LL*2+1)
cxx      DIMENSION CUY(LL*2+1),CYU(LL*2+1),CZY(LL*2+1)
cc      EQUIVALENCE (CXX(1),CXY(1))
cc      EQUIVALENCE (CYX(1),CUX(1))
cc      EQUIVALENCE (CYU(1),CYZ(1),CZX(1))
      INTEGER IP, IQ, IG, LL, ICST
      DOUBLE PRECISION X(IP+IQ), CYY(LL*2+1), G(IP+IQ), CXX0
c local
      DOUBLE PRECISION A(IP), B(IQ), AI(ICST), A2(IP*2), A2B(ICST),
     1                 AIB(ICST), Y(LL*2+1), CXX(LL*2+1), CYX(LL*2+1),
     2                 CUY(LL*2+1), CYU(LL*2+1), CZY(LL*2+1), CST0, CST1
c
      CST0=0.0D-00
      CST1=1.0D-00
cc      IORIG=501
      IORIG=LL+1
      IGA2B=IG
      DO 100 I=1,IP
cxx  100 A(I)=X(I)
      A(I)=X(I)
  100 CONTINUE
      IF(IQ.LE.0) GO TO 420
      DO 110 I=1,IQ
      II=IP+I
cxx  110 B(I)=X(II)
      B(I)=X(II)
  110 CONTINUE
  420 CONTINUE
      IB=IQ
C     INVERSE OF A(I) COMPUTATION
cc   24 CALL INVERS(A,IP,B,0,AI,IA,IG)
cx   24 CALL INVERS(A,IP,B,0,AI,IA,IG,IFG,LU)
cxx   24 CALL INVERS(A,IP,B,0,AI,IA,ICST,IG,IFG,LU)
      CALL INVERS(A,IP,B,0,AI,IA,ICST,IG)
      IF(IG.NE.1) GO TO 1900
      GO TO 1000
 1900 IF(IA.NE.0) GO TO 1901
      IA=1
      AI(1)=CST0
 1901 CONTINUE
      DO 2100 I=1,IP
cxx 2100 X(I)=A(I)
      X(I)=A(I)
 2100 CONTINUE
      IF(IQ.LE.0) GO TO 2420
      DO 2110 I=1,IQ
      II=IP+I
cxx 2110 X(II)=B(I)
      X(II)=B(I)
 2110 CONTINUE
 2420 CONTINUE
      IPM1=IP-1
      IQM1=IQ-1
C     AIB=(INVERSE OF A)*B
      IGAIB=0
cc   25 CALL INVERS(A,IP,B,IQ,AIB,IAIB,IGAIB)
cx   25 CALL INVERS(A,IP,B,IQ,AIB,IAIB,IGAIB,IFG,LU)
cxx   25 CALL INVERS(A,IP,B,IQ,AIB,IAIB,ICST,IGAIB,IFG,LU)
      CALL INVERS(A,IP,B,IQ,AIB,IAIB,ICST,IGAIB)
C     A2B=(INVERSE OF A*A)*B
C     A2=A*A
      Y(IORIG)=CST1
      DO 502 I=1,IP
      IJ=IORIG+I
      IK=IORIG-I
      A2B(I)=A(I)
      Y(IK)=A(I)
cxx  502 Y(IJ)=CST0
      Y(IJ)=CST0
  502 CONTINUE
      IK=IORIG-IP
      DO 503 I=1,IP
      IK=IK-1
cxx  503 Y(IK)=CST0
      Y(IK)=CST0
  503 CONTINUE
      L2=-IP-IP
      M2=-1
cc      CALL SCONVL(Y,A2B,Y,IP,L2,M2)
      CALL SCONVL(Y,A2B,Y,IP,L2,M2,LL)
      I2P=IP+IP
      DO 504 I=1,I2P
      IJ=IORIG-I
cxx  504 A2(I)=Y(IJ)
      A2(I)=Y(IJ)
  504 CONTINUE
cc      CALL INVERS(A2,I2P,B,IQ,A2B,IA2B,IGA2B)
cx      CALL INVERS(A2,I2P,B,IQ,A2B,IA2B,IGA2B,IFG,LU)
cxx      CALL INVERS(A2,I2P,B,IQ,A2B,IA2B,ICST,IGA2B,IFG,LU)
      CALL INVERS(A2,I2P,B,IQ,A2B,IA2B,ICST,IGA2B)
      LXX=0
      MXX=0
      LXY=0
      MXY=MXX+IAIB
      LYX=-MXY
      MYX=-LXY
      LZX=-IQ
      MZX=-1
      LZY=LZX
      MZY=MZX+IAIB
      LYZ=-MZY
      MYZ=-LZY
      LUX=-IP
      MUX=-1
      LUY=-IP
      MUY=MUX+IAIB
      LYU=-MUY
      MYU=-LUY
C     CXX0 COMPUTATION
C     CYX=CYY*AIB'
cc      CALL SCONVL(CYY,AIB,CYX,IAIB,LYX,MYX)
cc      CALL TURN(CYX,CXY,LYX,MYX)
      CALL SCONVL(CYY,AIB,CYX,IAIB,LYX,MYX,LL)
      CALL TURN(CYX,CXX,LYX,MYX,LL)
C     CXX=CXY*AIB'
cc      CALL SCONVL(CXY,AIB,CXX,IAIB,LXX,MXX)
      CALL SCONVL(CXX,AIB,CXX,IAIB,LXX,MXX,LL)
      CXX0=CXX(IORIG)
C     GA COMPUTATION
C     CYU=CYY*A2B'
cc      CALL SCONVL(CYY,A2B,CYU,IA2B,LYU,MYU)
cc      CALL TURN(CYU,CUY,LYU,MYU)
cc      CALL SCONVL(CUY,AIB,CUX,IAIB,LUX,MUX)
      CALL SCONVL(CYY,A2B,CYU,IA2B,LYU,MYU,LL)
      CALL TURN(CYU,CUY,LYU,MYU,LL)
      CALL SCONVL(CUY,AIB,CYX,IAIB,LUX,MUX,LL)
C     GB COMPUTATION
C     CYZ=CYY*AI'
      IF(IQ.EQ.0) GO TO 5279
cc      CALL SCONVL(CYY,AI,CYZ,IA,LYZ,MYZ)
cc      CALL TURN(CYZ,CZY,LYZ,MYZ)
cc      CALL SCONVL(CZY,AIB,CZX,IAIB,LZX,MZX)
      CALL SCONVL(CYY,AI,CYU,IA,LYZ,MYZ,LL)
      CALL TURN(CYU,CZY,LYZ,MYZ,LL)
      CALL SCONVL(CZY,AIB,CYU,IAIB,LZX,MZX,LL)
 5279 CONTINUE
C     GRADIENT ARRANGEMENT
      DO 5280 I=1,IP
      I1=IORIG-I
cc 5280 G(I)=-CUX(I1)
cxx 5280 G(I)=-CYX(I1)
      G(I)=-CYX(I1)
 5280 CONTINUE
      IF(IQ.LE.0) GO TO 5290
      DO 5281 I=1,IQ
      II=IP+I
      I1=IORIG-I
cc 5281 G(II)=CZX(I1)
cxx 5281 G(II)=CYU(I1)
      G(II)=CYU(I1)
 5281 CONTINUE
 5290 CONTINUE
C     CXX0, GRADIENT PRINT OUT
      IPQ=IP+IQ
 1000 RETURN
      END
C
cc      SUBROUTINE SDAV1(X,CYY,CXX0,G,C,IP,IQ,N)
cxx      SUBROUTINE SDAV1(X,CYY,CXX0,G,C,IP,IQ,N,VD,NN,LL,ICST,IFG,LU)
      SUBROUTINE SDAV1(X,CYY,CXX0,G,C,IP,IQ,N,VD,NN,LL,ICST)
C      DADIDON'S (MINIMIZATION) PROCEDURE
cxx      IMPLICIT REAL*8(A-H,O-Z)
cc      COMMON /COM50/VD
      COMMON /COM70/ISWRO
cc      DIMENSION VD(50,50)
cc      DIMENSION X(50),G(50),SX(50),SG(50),SR(50)
cc      DIMENSION C(50)
cxx      DIMENSION X(IP+IQ),CYY(LL*2+1)
cxx      DIMENSION G(IP+IQ),C(IP+IQ),VD(NN,NN)
cxx      DIMENSION SX(IP+IQ),SG(IP+IQ),SR(IP+IQ)
      INTEGER IP, IQ, N, NN, LL, ICST
      DOUBLE PRECISION X(IP+IQ), CYY(LL*2+1), CXX0, G(IP+IQ), C(IP+IQ),
     1                 VD(NN,NN)
c local
      DOUBLE PRECISION SX(IP+IQ), SG(IP+IQ), SR(IP+IQ), CST0, CST1,
     1                 CST2, CST05, CONSTA, CONSTB, EPS1, EPS3, EPS4,
     2                 AIPQ, AN, PHAI, EPHAI1, T1, RO, RAM, RAMRO,
     3                 RAMROT, SUM, SRO, SROD, DGAM, DGAM1, GSR, RAMT,
     4                 RAMSRO, RAM1, CONSDR, SPHAI, OAIC, OPHAI, AIC,
     5                 DAIC
C
C     CONSTANT
      CST0=0.0D-00
      CST1=1.0D-00
      CST2=2.0D-00
      CST05=0.5D-00
      CONSTA=0.5D-00
      CONSTB=2.0D-00
      EPS1=0.01D-00
      EPS3=0.000001D-00
      EPS4=0.1D-10
      ISPHAI=0
      ITN=1
      IPHAI=1
      IPQ=IP+IQ
      AIPQ=IPQ
      AN=N
      PHAI=CXX0
  150 CONTINUE
C     RO=G'*C COMPUTATION
      ITNS=0
C     COMMON SUBROUTINE CALL
   40 CALL INNERP(G,C,RO,IPQ)
      IF(IPHAI.EQ.0) GO TO 101
      PHAI=CXX0
  101 OPHAI=PHAI
      EPHAI1=EPS1*PHAI
      T1=RO-CST2*PHAI
      IF(T1.LE.EPHAI1) GO TO 140
      RAM=CST2*PHAI/RO
C     V=V+((RAM-1.0)/RO)*(C*C')
      RAMRO=(RAM-CST1)/RO
cxx      DO 110 I=1,IPQ
      DO 111 I=1,IPQ
      RAMROT=RAMRO*C(I)
      DO 110 J=1,IPQ
cxx  110 VD(I,J)=VD(I,J)+RAMROT*C(J)
      VD(I,J)=VD(I,J)+RAMROT*C(J)
  110 CONTINUE
  111 CONTINUE
C     C=RAM*C
      DO 120 I=1,IPQ
cxx  120 C(I)=RAM*C(I)
      C(I)=RAM*C(I)
  120 CONTINUE
      IF(ITNS.GE.10) GO TO 140
      ITNS=ITNS+1
      GO TO 40
C     SX=X-R
  140 CONTINUE
      IG=0
 1210 CONTINUE
      DO 210 I=1,IPQ
cxx  210 SX(I)=X(I)-C(I)
      SX(I)=X(I)-C(I)
  210 CONTINUE
C     SPHAI=CXX0, SG=GRADIENT COMPUTATION
cc      CALL SC0GR1(SX,CYY,SG,SPHAI,IP,IQ,IG,LL,ICST)
cxx      CALL SC0GR1(SX,CYY,SG,SPHAI,IP,IQ,IG,LL,ICST,IFG,LU)
      CALL SC0GR1(SX,CYY,SG,SPHAI,IP,IQ,IG,LL,ICST)
      IF(IG.NE.1) GO TO 309
cxx      DO 303 I=1,IPQ
      DO 304 I=1,IPQ
      C(I)=CST05*C(I)
      DO 303 J=1,IPQ
cxx  303 VD(I,J)=CST05*VD(I,J)
      VD(I,J)=CST05*VD(I,J)
  303 CONTINUE
  304 CONTINUE
      GO TO 1210
  309 CONTINUE
C     SR=V*SG
      DO 310 I=1,IPQ
      SUM=CST0
      DO 311 J=1,IPQ
cxx  311 SUM=SUM+VD(I,J)*SG(J)
      SUM=SUM+VD(I,J)*SG(J)
  311 CONTINUE
cxx  310 SR(I)=SUM
      SR(I)=SUM
  310 CONTINUE
C     SRO=(SG)'*(SR)
C     COMMON SUBROUTINE CALL
      CALL INNERP(SG,SR,SRO,IPQ)
      SROD=SRO/PHAI
C     DGAM=-G'*(SR)/SRO
C     COMMON SUBROUTINE CALL
      CALL INNERP(G,SR,GSR,IPQ)
      DGAM=-GSR/SRO
      DGAM1=DGAM+CST1
      DGAM1=DABS(DGAM1)+0.1D-70
      RAM=DABS(DGAM)/DGAM1
C     IF RAM . LE. CONSTA THEN RAM=CONSTA
      IF(RAM.GT.CONSTA) GO TO 430
      RAM=CONSTA
      IRAM=1
      GO TO 470
C     IF RAM . GE. CONSTB THEN RAM=CONSTB
  430 IF(RAM.LT.CONSTB) GO TO 450
      RAM=CONSTB
      IRAM=-1
      GO TO 470
C     RAM=RAM
  450 CONTINUE
      IRAM=0
C     V=V+((RAM-1.0)/SRO)*(SR)*(SR)'
  470 RAMSRO=(RAM-CST1)/SRO
cxx      DO 480 I=1,IPQ
      DO 481 I=1,IPQ
      RAMT=RAMSRO*SR(I)
      DO 480 J=1,IPQ
cxx  480 VD(I,J)=VD(I,J)+RAMT*SR(J)
      VD(I,J)=VD(I,J)+RAMT*SR(J)
  480 CONTINUE
  481 CONTINUE
      IF(PHAI.GE.SPHAI) GO TO 540
C     SPHAI.GT.PHAI: TEST OF CORRECTION
      RAM1=RAM-CST1
      IF(DABS(RAM1).LT.EPS3) GO TO 555
      CONSDR=DGAM*RAM1
      DO 550 I=1,IPQ
cxx  550 C(I)=C(I)-CONSDR*SR(I)
      C(I)=C(I)-CONSDR*SR(I)
  550 CONTINUE
      IPHAI=0
      IF(SROD.GT.EPS4) GO TO 900
C     END OF ITERATION
  555 ISWRO=ISWRO+1
      GO TO 1000
C     SPHAI LE. PHAI: SUCCESSFUL REDUCTION
  540 DO 560 I=1,IPQ
      X(I)=SX(I)
      G(I)=SG(I)
cxx  560 C(I)=RAM*SR(I)
      C(I)=RAM*SR(I)
  560 CONTINUE
      CXX0=SPHAI
      PHAI=SPHAI
      IPHAI=1
cxx  800 CONTINUE
      OAIC=AN*DLOG(OPHAI)+CST2*AIPQ
      AIC=AN*DLOG(PHAI)+CST2*AIPQ
      DAIC=OAIC-AIC
      IF(IRAM.NE.0) GO TO 901
      IF(SROD.LT.EPS4) GO TO 555
C     ITERATION CHECK
  900 IPQ2=IPQ+IPQ
      IF(ITN.GE.IPQ2) GO TO 555
      ISPHAI=(ISPHAI+(1-IPHAI))*(1-IPHAI)
      IF(ISPHAI.GT.10) GO TO 555
      ITN=ITN+1
      GO TO 150
  901 IF(SROD.LT.EPS4) GO TO 555
      GO TO 900
C     END OF MINIMIZATION
cxx  999 ISWRO=0
      ISWRO=0
 1000 CONTINUE
cxx 1001 RETURN
      RETURN
      END
C
cc      SUBROUTINE SCONVL(Y,A,Z,K,L,M)
      SUBROUTINE SCONVL(Y,A,Z,K,L,M,LL)
C     Y(I), Z(I) CENTERED AT I=IORIG
C     Z(I)=Y(I)+Y(I+1)A(1)+...+Y(I+K)A(K) (I=L,M)
cxx      IMPLICIT REAL*8(A-H,O-Z)
cc      DIMENSION Y(1001),A(190),Z(1001)
cxx      DIMENSION
      INTEGER K, L, M, LL
      DOUBLE PRECISION Y(LL*2+1), A(K), Z(LL*2+1)
c local
      DOUBLE PRECISION SUM
cc      IORIG=501
      IORIG=LL+1
      IST=IORIG+L
      IEN=IORIG+M
      DO 3 I=IST,IEN
      SUM=Y(I)
      DO 2 J=1,K
      IJ=I+J
cxx    2 SUM=SUM+Y(IJ)*A(J)
      SUM=SUM+Y(IJ)*A(J)
    2 CONTINUE
cxx    3 Z(I)=SUM
      Z(I)=SUM
    3 CONTINUE
      RETURN
      END
C
cc      SUBROUTINE TURN(Y,Z,L,M)
      SUBROUTINE TURN(Y,Z,L,M,LL)
C     Z(IORIG+I)=Y(IORIG-I) (I=1,M)
cxx      IMPLICIT REAL*8(A-H,O-Z)
cc      DIMENSION Z(1001),Y(1001)
cxx      DIMENSION Z(LL-L+1),Y(LL+M+1)
      INTEGER L, M, LL
      DOUBLE PRECISION Y(LL+M+1), Z(LL-L+1)
cc      IORIG=501
      IORIG=LL+1
      IST=IORIG+L
      IEN=IORIG+M
      DO 1 I=IST,IEN
      IJ=IORIG-(I-IORIG)
cxx    1 Z(IJ)=Y(I)
      Z(IJ)=Y(I)
    1 CONTINUE
      RETURN
      END
C
cc      SUBROUTINE INVERS(A,IP,B,IQ,X,IX,IG)
cx      SUBROUTINE INVERS(A,IP,B,IQ,X,IX,IG,IFG,LU)
cxx      SUBROUTINE INVERS(A,IP,B,IQ,X,IX,ICST,IG,IFG,LU)
      SUBROUTINE INVERS(A,IP,B,IQ,X,IX,ICST,IG)
C     X=(INVERSE OF B )*A
C     W(I)+B(1)W(I-1)+...B(IQ)W(I-IQ)=X(I)+A(1)X(I-1)+...+A(IP)X(I-IP)
C     INPUT W(0)=1, W(I)=0 FOR I. NE. 0
C     OUTPUT X(I) (I=1,IX)
cxx      IMPLICIT REAL*8(A-H,O-Z)
cx      DIMENSION A(1),B(1),X(1)
cxx      DIMENSION A(IP),B(IQ),X(ICST)
      INTEGER IP, IQ, IX, ICST, IG
      DOUBLE PRECISION A(IP), B(IQ), X(ICST)
c local
      DOUBLE PRECISION CST0, GCONST, GAMMAX, SUM, GAM2
      CST0=0.0D-00
      IPQ=IP+IQ
      IF(IPQ.LE.0) GO TO 999
      GCONST=0.0005D-00
      GAMMAX=1.0D+10
      K=0
      LH=6
      IH=0
      IF(IG.EQ.0) GO TO 13
      GCONST=0.01D-00
      IG=0
cx   13 DO 10 I=1,190
   13 DO 10 I=1,ICST
      IX=I
      SUM=CST0
      IF(I.GT.IQ) GO TO 2
      SUM=B(I)
    2 IF(I.GT.IP) GO TO 3
      SUM=SUM-A(I)
    3 IM1=I-1
      JM=MIN0(IM1,IP)
      IF(JM.LE.0) GO TO 5
      DO 4 J=1,JM
      IMJ=I-J
cxx    4 SUM=SUM-X(IMJ)*A(J)
      SUM=SUM-X(IMJ)*A(J)
    4 CONTINUE
    5 X(I)=SUM
      GAM2=DABS(SUM)
      IF(GAM2.GE.GCONST) GO TO 24
      IH=IH+1
      IF(IH.LT.LH) GO TO 10
      GO TO 1000
   24 IF(GAM2.LE.GAMMAX) GO TO 26
      IG=1
cc	WRITE(6,60)
cx      IF (IFG.NE.0) WRITE(LU,60)
      GO TO 1000
   26 IH=0
   10 CONTINUE
      IF(IH.GE.LH) GO TO 1000
      IG=1
cc      WRITE(6,59)
cx      IF (IFG.NE.0) WRITE(LU,59)
      GO TO 1000
  999 IX=0
 1000 RETURN
cxx   59 FORMAT(1H ,'INCOMPLETE CONVERGENCE OF INVERSE')
cxx   60 FORMAT(1H ,'DIVERGENT INVERSE')
      END
