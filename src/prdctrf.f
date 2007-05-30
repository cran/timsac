      SUBROUTINE PRDCTRF(N,P,Q,H,D,K,L,JSW,YY,B,A,WW,S,Y,YORI,YD,X,Z1,
     *Z2,Z3,ZZ1,ZZ2,ZZ3,TMP)
C
      INCLUDE 'timsac_f.h'
C
cc	PROGRAM PRDCTR
C     PROGRAM 74.3.1.  PREDICTION PROGRAM
C-----------------------------------------------------------------------
C     ** DESIGNED BY H. AKAIKE, THE INSTITUTE OF STATISTICAL MATHEMATICS
C     ** PROGRAMMED BY E. ARAHATA, THE INSTITUTE OF STATISTICAL MATHEMAT
C	  TOKYO
C     ** DATE OF THE LATEST REVISION: MARCH 25, 1977
C     ** THIS PROGRAM WAS ORIGINALLY PUBLISHED IN
C	 "TIMSAC-74 A TIME SERIES ANALYSIS AND CONTROL PROGRAM PACKAGE(2
C	 BY H. AKAIKE, E. ARAHATA AND T. OZAKI, COMPUTER SCIENCE MONOGRA
C	 NO.6 MARCH 1976, THE INSTITUTE OF STATISTICAL MATHEMATICS
C-----------------------------------------------------------------------
C     THIS PROGRAM OPERATES ON A REAL RECORD OF A VECTOR PROCESS
C     Y(I) (I=1,N) AND COMPUTES PREDICTED VALUES. ONE STEP AHEAD
C     PREDICTION STARTS AT TIME P AND ENDS AT TIME Q. PREDICTION IS
C     CONTINUED WITHOUT NEW OBSERVATIONS UNTIL TIME Q+H.
C     BASIC MODEL IS THE AUTOREGRESSIVE MOVING AVERAGE
C     MODEL OF Y(I) WHICH IS GIVEN BY
C     Y(I)+B(1)Y(I-1)+...+B(K)Y(I-K) = X(I)+A(1)X(I-1)+...+A(L)X(I-L).
C
C    THE FOLLOWING INPUTS ARE REQUIRED:
C     (N,P,Q,H):
C		 N, LENGTH OF DATA
C		 P, ONE STEP AHEAD PREDICTION STARTING POSITION
C		 Q, LONG RANGE FORECAST STARTING POSITION
C		 H, MAXIMUM SPAN OF FORECAST (LESS THAN OR EQUAL TO 100)
C		 (Q+H MUST BE LESS THAN 1001)
C     JSW: JSW=0 FOR DIRECT LOADING OF AR-MA COEFFICIENTS,
C		 THE OUTPUTS OF PROGRAM MARKOV WITH ICONT=0.
C	   JSW=1 FOR LOADING OF THE OUTPUTS OF PROGRAM MARKOV,
C		 THE OUTPUTS OF PROGRAM MARKOV WITH ICONT=1.
C     (D,K,L):
C	       D, DIMENSION OF THE VECTOR Y(I)
C	       K, AR-ORDER (LESS THAN OR EQUAL TO 10)
C	       L, MA-ORDER (LESS THAN OR EQUAL TO 10)
C     N,L,K,H,P,Q,D,JSW,ARE ALL INTEGERS
C     (DFORM(I),I=1,20): INPUT FORMAT STATEMENT IN ONE CARD,
C			 FOR EXAMPLE, (8F10.4)
C     (NAME(I,J),I=1,20,J=1,D): NAME OF THE I-TH COMPONENT
C     (Y(I,J),I=1,N;J=1,D): ORIGINAL DATA
C     (B(I1,I2,J),I1=1,D,I2=1,D,J=1,K): AR-COEFFICIENT MATRICES.
C     FOR JSW=0,
C	 (A(I1,I2,J),I1=1,D,I2=1,D,J=1,L): MA-COEFFICIENT MATRICES.
C     FOR JSW=1,
C	 (W(I1,I2,J),I1=1,D,I2=1,D,J=1,L): IMPULSE RESPONSE MATRICES.
C     (S(I,J),I=1,D,J=1,D): INNOVATION VARIANCE MATRIX
C
C     THE OUTPUTS OF THIS PROGRAM ARE THE REAL
C     AND PREDICTED VALUES OF Y(I).
C
cc      !DEC$ ATTRIBUTES DLLEXPORT :: PRDCTRF
C
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER H,H1,P,Q,D
cc	INTEGER DFORM(20),XX(121)
cc	INTEGER IKA/1HA/, IKB/1HB/, IKW/1HW/
cc	REAL*8 A,B,W,S,Z,C,CX
cc	REAL*8 CST0,CST1
cc	DIMENSION A(10,10,10),B(10,10,10),S(10,10),W(10,10,101)
cc	DIMENSION X(1000,10),C(10,10),SD(10,101),CX(10),EY(10),NAME(20,10)
cc	DIMENSION Y(1000,10),YORI(500,10)
cc	DIMENSION DMXT(10),DMIT(10),AV(10)
      DIMENSION A(D,D,L),B(D,D,K),S(D,D),WW(D,D,L),W(D,D,H+1)
      DIMENSION X(N,D),C(D,D),SD(D,H+1),CX(D),EY(D)
      DIMENSION YY(N,D),YORI(H+1,D)
      DIMENSION DMXT(D),DMIT(D),AV(D)
cc	EQUIVALENCE(W(1,1,1),X(1,1))
      DIMENSION Y(Q+H,D),YD(Q+H,D)
      DIMENSION Z1(Q+H,D),Z2(Q+H,D),Z3(Q+H,D)
      DIMENSION ZZ1(Q+H,D),ZZ2(Q+H,D),ZZ3(Q+H,D)
C
      INTEGER*1  TMP(1)
      CHARACTER  CNAME*80
C
C      DATA A/1000*0.0D-00/, B/1000*0.0D-00/, W/10000*0.0D-00/
cc	DATA A/1000*0.0D-00/, B/1000*0.0D-00/, W/10100*0.0D-00/
cc	DATA S/100*0.0D-00/, C/100*0.0D-00/
cc	DATA CX/10*0.0D-00/
cc	DATA X/10000*0.0/,SD/1010*0.0/,EY/10*0.0/,Y/10000*0.0/
cc	DATA YORI/5000*0.0/,DMXT/10*0.0/,DMIT/10*0.0/,AV/10*0.0/
cc	DATA X/10000*0.0D-00/,SD/1010*0.0D-00/,EY/10*0.0D-00/
cc	DATA Y/10000*0.0D-00/,YORI/5000*0.0D-00/
cc	DATA DMXT/10*0.0/,DMIT/10*0.0/,AV/10*0.0/
cc      INTEGER XX(121)
cc      DATA XX/121*1H /
cc      DATA K1 / 1H* /, K2 / 1HX /, K3 / 1H+ /, K4 / 1HY /, K5 / 1H  /
cc      DATA K6/1H!/
      CHARACTER XX(121)
      CHARACTER K1, K2, K3, K4, K5, K6
      CHARACTER KSTOR
      DATA XX/121*' '/
      DATA K1,K2,K3,K4,K5,K6 / '*', 'X', '+', 'Y', ' ', '!' /
C
C     INPUT / OUTPUT DATA FILE OPEN
cc	CALL SETWND
cc	CALL FLOPN2(NFL)
cc	IF (NFL.EQ.0) GO TO 999
      LU=3
      DO 100 I = 1,80
  100 CNAME(I:I) = ' '
      I = 1
      IFG = 1
      DO WHILE( (IFG.EQ.1) .AND. (I.LE.80) )
	   IF ( TMP(I).NE.ICHAR(' ') ) THEN
            CNAME(I:I) = CHAR(TMP(I))
            I = I+1
         ELSE
            IFG = 0
         END IF
      END DO
      IF ( I.GT.1 ) THEN
         IFG = 1
         OPEN (LU,FILE=CNAME,IOSTAT=IVAR)
         IF (IVAR .NE. 0) THEN
            WRITE(*,*) ' ***  prdctr temp FILE OPEN ERROR :',CNAME,IVAR
            IFG=0
         END IF
      END IF
C
cc	MJ1=10
cc	MJ2=10
      CST0=0.0D-00
      CST1=1.0D-00
      ISW=1
C
C      DATA INITIALIZE
C
      IF (JSW.NE.0) CALL DINIT(A,D*D*L,CST0)
      IF (JSW.EQ.0) CALL DINIT(WW,D*D*L,CST0)
      CALL DINIT(W,D*D*(H+1),CST0)
      CALL DINIT(C,D*D,CST0)
      CALL DINIT(CX,D,CST0)
      CALL DINIT(X,N*D,CST0)
      CALL DINIT(SD,D*(H+1),CST0)
      CALL DINIT(EY,D,CST0)
      CALL DINIT(Y,(Q+H)*D,CST0)
      CALL DINIT(YORI,(H+1)*D,CST0)
      CALL DINIT(DMXT,D,CST0)
      CALL DINIT(DMIT,D,CST0)
      CALL DINIT(AV,D,SCT0)
C
C     INITIAL CONDITION INPUT
cc	READ(5,800) N,D,P,Q,H
cc	READ(5,800) JSW
cc	DO 199 J=1,D
cc	READ(5,801) (NAME(I,J),I=1,20)
cc  199 CONTINUE
cc	READ(5,801) (DFORM(I),I=1,20)
C     ORIGINAL DATA VECTOR (Y(I),I=1,N) INPUT
cc	DO 100 J=1,D
cc	READ(5,DFORM) (Y(I,J),I=1,N)
cc  100 CONTINUE
      DO 70 I = 1,N
      DO 70 J = 1,D
	 Y(I,J) = YY(I,J)
   70 CONTINUE
C
cc	READ(5,800) D,K,L
cc	WRITE(6,900)
cc	WRITE(6,901) N,D,K,L,P,Q,H,JSW
cc	WRITE(6,902) (DFORM(I),I=1,20)
cc	DO 1101 J=1,D
cc	WRITE(6,903) J,(NAME(I,J),I=1,20)
cc	WRITE(6,904) (Y(I,J),I=1,N)
cc 1101 CONTINUE
C
C     AR-COEFFICIENTS INPUT
cc	CALL REMT3X(B,D,D,K,ISW,MJ1,MJ1,MJ2,IKB)
cc	IF  (JSW.NE.0) GO TO 80
C     MA-COEFFICIENTS INPUT
cc	CALL REMT3X(A,D,D,L,ISW,MJ1,MJ1,MJ2,IKA)
cc	GO TO 81
C     IMPULSE RESPONSE INPUT
cc   80 CALL REMT3X(W,D,D,L,ISW,MJ1,MJ1,MJ2,IKW)
      DO 80 II = 1,D
      DO 80 JJ = 1,D
      DO 80 KK = 1,L
	 W(II,JJ,KK) = WW(II,JJ,KK)
   80 CONTINUE
C     INNOVATION VARIANCE INPUT
C     COMMON SUBROUTINE CALL
cc   81 CALL REMATX(S,D,D,ISW,MJ1,MJ1)
cc	WRITE(6,1357)
cc	CALL SUBMPR(S,D,D,MJ1,MJ1)
      IF  (JSW.EQ.0) GO TO 99
C
C     MA-COEFFICIENTS COMPUTATION
C     A(M)=W(M)+B(1)W(M-1)+...+B(M-1)W(1)+B(M)
C     M=1,L
      DO  84 I=1,L
      DO  83 JX=1,D
      DO  82 JY=1,D
   82 A(JX,JY,I)=W(JX,JY,I)+B(JX,JY,I)
   83 CONTINUE
   84 CONTINUE
      DO  89 I=2,L
      IM1=I-1
      DO  88 IX=1,IM1
      IY=I-IX
      DO  87 JX=1,D
      DO  86 JY=1,D
      Z=CST0
      DO  85 JZ=1,D
   85 Z=Z+B(JX,JZ,IX)*W(JZ,JY,IY)
   86 A(JX,JY,I)=A(JX,JY,I)+Z
   87 CONTINUE
   88 CONTINUE
   89 CONTINUE
      DO 8210 I=1,L
      DO 8220 IX=1,D
      DO 8220 JX=1,D
 8220 C(IX,JX)=A(IX,JX,I)
cc	WRITE(6,8121) I
cc	CALL SUBMPR(C,D,D,MJ1,MJ1)
 8210 CONTINUE
C
C     W(I) CLEAR
      DO  98 I=1,L
      DO  97 JX=1,D
      DO  96 JY=1,D
   96 W(JX,JY,I)=CST0
   97 CONTINUE
   98 CONTINUE
C
C     IMPULSE RESPONSE W(M),M=0,H COMPUTATION
C     W(0) = UNIT MATRIX
C     W(M) = A(M)-(B(1)W(M-1)+...+B(K)W(M-K)) (FOR M LESS THAN OR EQUAL
C     W(M) =-B(1)W(M-1)-...-B(K)W(M-K) (FOR M GREATER THAN L)
   99 DO  101 I=1,D
  101 W(I,I,1)=CST1
      H1=H+1
      DO 120 I=2,H1
      DO  115 J=1,K
      IJ=I-J
cc      IF  (IJ) 115,115,102
      IF  (IJ.LE.0) GO TO 115
  102 DO  105 JX=1,D
      DO  104 IX=1,D
      Z=CST0
      DO  103 I2=1,D
  103 Z=Z+(B(IX,I2,J)*W(I2,JX,IJ))
      W(IX,JX,I) = W(IX,JX,I) - Z
  104 CONTINUE
  105 CONTINUE
  115 CONTINUE
cc      IF  (I-L-1) 116,116,120
      IF  (I-L-1.GT.0) GO TO 120
  116 DO  118 JX=1,D
      DO  117 IX=1,D
  117 W(IX,JX,I) = W(IX,JX,I)+A(IX,JX,I-1)
  118 CONTINUE
  120 CONTINUE
C
C     PREDICTION ERROR VARIANCE C(I) AND STANDARD DEVIATION SD(IX,I)
C     COMPUTATION (I=1,H1)
C     C(I)=W(0)*S*W(0)'+...+W(I)*S*W(I)' AND SD(IX,I) IS THE VECTOR OF
C     THE POSITIVE SQUARE ROOTS OF THE DIAGONAL ELEMENTS OF THE PREDICTI
C     ERROR VARIANCE MATRIX C(I).  ONLY THE DIAGONAL ELEMENTS OF C(I)
C     ARE COMPUTED.
      DO 140 I=1,H1
      DO  125 IX=1,D
      DO  124 JX=1,D
      Z=CST0
      DO 123 I2=1,D
  123 Z=Z+W(IX,I2,I)*S(I2,JX)
      C(IX,JX)=Z
  124 CONTINUE
  125 CONTINUE
      DO  139 IX=1,D
      Z=CST0
      DO  129 I2=1,D
  129 Z=Z+C(IX,I2)*W(IX,I2,I)
      CX(IX)=CX(IX)+Z
      CCC=CX(IX)
cc      IF  (CCC) 130,130,131
cc  130 SD(IX,I)=0.0
      IF  (CCC.GT.0) GO TO 131
  130 SD(IX,I)=0.0D-00
      GO TO 139
cc  131 SD(IX,I)=SQRT(CCC)
  131 SD(IX,I)=DSQRT(CCC)
  139 CONTINUE
  140 CONTINUE
C
C     SUBTRACTION OF THE MEAN VALUES FROM THE ORIGINAL DATA Y(I)
      AN=N
cc	AN=1.0/AN
      AN=1.0D-00/AN
      DO  146  I=1,D
      Z=CST0
      DO  141  J=1,N
  141 Z=Z+Y(J,I)
      AVE=Z*AN
      AV(I)=AVE
      DO  145  J=1,N
      Y(J,I)=Y(J,I)-AVE
  145 CONTINUE
  146 CONTINUE
C
C
C     PREDICTIONS AND INNOVATIONS (X(I)) COMPUTATION OF Y(I) (I=1,Q+H)
C     FOR I GREATER THAN OR EQUAL TO Q, X(I) IS SET EQUAL TO 0.
C
      IQH=Q+H
C
      DO 6300 J=1,D
      DMXT(J)=Y(1,J)
 6300 DMIT(J)=Y(1,J)
      ISR=0
      DO 300 I=1,IQH
C     EY, PREDICTED VALUE OF Y, COMPUTATION
      DO  153 J=1,D
cc  153 EY(J)=0.0
  153 EY(J)=0.0D-00
C     B(1)Y(I-1)+...+B(K)Y(I-K) COMPUTATION
      DO  160 J=1,K
      IJ=I-J
cc      IF  (IJ) 160,160,154
      IF  (IJ.LE.0) GO TO 160
  154 DO  159 IX=1,D
      Z=CST0
      DO  158 I2=1,D
  158 Z=Z+B(IX,I2,J)*Y(IJ,I2)
      EY(IX)=EY(IX)-Z
  159 CONTINUE
  160 CONTINUE
C     A(1)X(I-1)+...+A(L)X(I-L) COMPUTATION
      DO  170 J=1,L
      IJ=I-J
cc      IF  (IJ) 170,170,161
cc  161 IF  (IJ-Q) 162,170,170
      IF  (IJ.LE.0) GO TO 170
  161 IF  (IJ-Q.GE.0) GO TO 170
  162 DO  169 IX=1,D
      Z=CST0
      DO  168 I2=1,D
  168 Z=Z+(A(IX,I2,J)*X(IJ,I2))
      EY(IX)=EY(IX)+Z
  169 CONTINUE
  170 CONTINUE
C     MAXIMUM, MINIMUM SEARCH
      DO  249 J=1,D
      DMAX=DMXT(J)
      DMIN=DMIT(J)
cc	YD=Y(I,J)
cc	CALL MAXMIN(DMAX,DMIN,YD)
cc	YD=EY(J)
cc	CALL MAXMIN(DMAX,DMIN,YD)
      YYD=Y(I,J)
      CALL MAXMIN(DMAX,DMIN,YYD)
      YYD=EY(J)
      CALL MAXMIN(DMAX,DMIN,YYD)
      DMXT(J)=DMAX
      DMIT(J)=DMIN
  249 CONTINUE
C     INNOVATION X(I) IS GIVEN BY (Y-EY).
C     IF I IS GREATER THAN OR EQUAL TO Q, THEN X(I)=0.
cc      IF(I-Q) 171,200,200
      IF(I-Q.GE.0) GO TO 200
  171 DO  172  J=1,D
      X(I,J)=Y(I,J)-EY(J)
  172 CONTINUE
      GO TO 300
  200 ISR=ISR+1
      DO 250 J=1,D
      YORI(ISR,J)=Y(I,J)+AV(J)
      Y(I,J)=EY(J)
  250 CONTINUE
  300 CONTINUE
C
      DO 281 I=1,IQH
      DO  280 J=1,D
  280 Y(I,J)=Y(I,J)+AV(J)
  281 CONTINUE
C
C
C     ********************
C     PRINT OUT
C     ********************
      I2=0
cc	CALL HEADPR(1,P,Q,D,N)
      DO 500 II=1,IQH
      I=II
cc	IF  (I-Q)   301,449,450
cc  301 IF  (I-P)   302,399,400
cc      IF  (I-Q)	  301,450,450
cc  301 IF  (I-P)	  500,400,400
      IF  (I-Q.GE.0) GO TO 450
  301 IF  (I-P.LT.0) GO TO 500
cc  302 DO  310 J=1,D
cc	IF  (J-1) 303,303,304
cc  303 WRITE(6,910) I,Y(I,J)
cc	GO TO 310
cc  304 WRITE(6,911) Y(I,J)
cc  310 CONTINUE
cc	WRITE(6,920)
cc	GO TO 500
cc  399 CALL HEADPR(I,P,Q,D,N)
  400 DO  410 J=1,D
cc	YD=Y(I,J)-X(I,J)
      YD(I,J)=Y(I,J)-X(I,J)
cc	IF  (J-1) 403,403,404
cc  403 WRITE(6,910) I,Y(I,J),YD,X(I,J)
cc	GO TO 410
cc  404 WRITE(6,911) Y(I,J),YD,X(I,J)
  410 CONTINUE
cc	WRITE(6,920)
      GO TO 500
cc  449 CALL HEADPR(I,P,Q,D,N)
  450 I2=I2+1
      DO  480 J=1,D
cc	YD=Y(I,J)
cc	Z1=YD+SD(J,I2)
cc	Z2=Z1+SD(J,I2)
cc	Z3=Z2+SD(J,I2)
      YD(I,J)=Y(I,J)
      Z1(I,J)=YD(I,J)+SD(J,I2)
      Z2(I,J)=Z1(I,J)+SD(J,I2)
      Z3(I,J)=Z2(I,J)+SD(J,I2)
c      IF  (J-1) 453,453,454
c  453 WRITE(6,912) I,YD,Z1,Z2,Z3
c      IF(I.LE.N) WRITE(6,937) YORI(I2,J)
c      GO TO 455
c  454 WRITE(6,913) YD,Z1,Z2,Z3
c      IF(I.LE.N) WRITE(6,937) YORI(I2,J)
cc	IF  ((J-1).LE.0) THEN
cc	   IF (I.LE.N) THEN
cc	      WRITE(6,9370) I,YORI(I2,J),YD,Z1,Z2,Z3
cc	   ELSE
cc	      WRITE(6,912) I,YD,Z1,Z2,Z3
cc	   END IF
cc	ELSE
cc	   IF (I.LE.N) THEN
cc	      WRITE(6,9371) YORI(I2,J),YD,Z1,Z2,Z3
cc	   ELSE
cc	      WRITE(6,913) YD,Z1,Z2,Z3
cc	   END IF
cc	END IF
cc  455 Z1=YD-SD(J,I2)
cc	Z2=Z1-SD(J,I2)
cc	Z3=Z2-SD(J,I2)
  455 ZZ1(I,J)=YD(I,J)-SD(J,I2)
      ZZ2(I,J)=ZZ1(I,J)-SD(J,I2)
      ZZ3(I,J)=ZZ2(I,J)-SD(J,I2)
cc	WRITE(6,914) Z1,Z2,Z3
  480 CONTINUE
cc	WRITE(6,920)
  500 CONTINUE
C
C
C     ********************
C     GRAPHIC PRINT OUT
C     ********************
      IF (IFG .NE. 0) THEN
      KSTOR=K2
      DO  600  J=1,D
      DMAX=DMXT(J)
      DMIN=DMIT(J)
cc	FMAX=ABS(DMAX)
cc	FMIN=ABS(DMIN)
cc	FMAX=AMAX1(FMAX,FMIN)
      FMAX=DABS(DMAX)
      FMIN=DABS(DMIN)
      FMAX=DMAX1(FMAX,FMIN)
      FMIN=-FMAX
cc	YST=FMAX/60.0
      YST=FMAX/60.0D-00
      TFMIN=FMIN+AV(J)
      TFMID=AV(J)
      TFMAX=FMAX+AV(J)
cc	WRITE(6,931) J,(NAME(I,J),I=1,20)
cc	WRITE(6,932)
cc	WRITE(6,933)
cc	WRITE(6,934) TFMIN,TFMID,TFMAX
cc	WRITE(6,935)
      WRITE(LU,931) J
      WRITE(LU,932)
      WRITE(LU,933)
      WRITE(LU,934) TFMIN,TFMID,TFMAX
      WRITE(LU,935)
C
      DO 599 I=1,IQH
      YA=Y(I,J)-AV(J)
      XX(61)=K6
cc      IF  (I-Q)	 501,569,570
cc  501 IF  (I-P)	 502,550,550
      IF  (I-Q.EQ.0) GO TO 569
      IF  (I-Q.GT.0) GO TO 570
  501 IF  (I-P.GE.0) GO TO 550
C
C     REAL DATA
  502 CALL SBSCAL(YA,YST,IX)
      XX(IX)=K1
cc	WRITE(6,936) I,(XX(I2),I2=1,121)
      WRITE(LU,936) I,(XX(I2),I2=1,121)
      XX(IX)=K5
      GO TO 599
C
C     REAL DATA
  550 CALL SBSCAL(YA,YST,IX)
C     ONE-STEP PREDICTION
      YB=YA-X(I,J)
      CALL SBSCAL(YB,YST,JX)
C
cc  566 IF  (IX-JX) 568,567,568
  566 IF  (IX-JX.NE.0) GO TO 568
  567 XX(IX)=K3
cc	WRITE(6,936) I,(XX(I2),I2=1,121)
      WRITE(LU,936) I,(XX(I2),I2=1,121)
      XX(IX)=K5
      GO TO 599
C
  568 XX(IX)=K1
      XX(JX)=K2
cc	WRITE(6,936) I,(XX(I2),I2=1,121)
      WRITE(LU,936) I,(XX(I2),I2=1,121)
      XX(IX)=K5
      XX(JX)=K5
      GO TO 599
C
  569 KSTOR=K2
      K2=K4
      ISR=0
C
C     LONG RANGE PREDICTION
  570 CALL SBSCAL(YA,YST,JX)
cc  576 IF  (I-N) 577,577,598
  576 IF  (I-N.GT.0) GO TO 598
C     REAL DATA
  577 ISR=ISR+1
      YC=YORI(ISR,J)-AV(J)
      CALL SBSCAL(YC,YST,IX)
      GO TO 566
C
  598 XX(JX)=K4
cc	WRITE(6,936) I,(XX(I2),I2=1,121)
      WRITE(LU,936) I,(XX(I2),I2=1,121)
      XX(JX)=K5
  599 CONTINUE
      K2=KSTOR
  600 CONTINUE
      END IF
C
cc	CALL FLCLS2(NFL)
cc  999 CONTINUE
cc	STOP
      IF (IFG.NE.0) CLOSE(LU)
  800 FORMAT(8I5)
  801 FORMAT(20A4)
  900 FORMAT(1H ,' PROGRAM 74.3.1. PREDICTION')
  901 FORMAT(/1H ,' INITIAL CONDITION: N=',I4,', D=',I2,', K=',I2,
     A', L=',I4,', P=',I4,', Q=',I4,', H=',I4,', JSW=',I1)
  902 FORMAT(/1H ,' ORIGINAL DATA (',20A4,'):')
  903 FORMAT(1H ,' ** Y ',I2,2X,20A4)
  904 FORMAT(1H ,10E12.5)
cc  910 FORMAT(1H ,'   N=',I5,4X,E20.5,4(4X,E20.5))
cc  911 FORMAT(1H ,14X,E20.5,4(4X,E20.5))
cc  912 FORMAT(1H ,'   N=',I5,24X,4(4X,E20.5))
cc  913 FORMAT(1H ,34X,4(4X,E20.5))
cc  914 FORMAT(1H ,58X,3(4X,E20.5))
  910 FORMAT(1H ,'   N=',I5,4X,D20.5,4(4X,D20.5))
  911 FORMAT(1H ,14X,D20.5,4(4X,D20.5))
  912 FORMAT(1H ,'   N=',I5,24X,4(4X,D20.5))
  913 FORMAT(1H ,34X,4(4X,D20.5))
  914 FORMAT(1H ,58X,3(4X,D20.5))
  920 FORMAT(1H ,4X)
  931 FORMAT(/1H ,10X,'J=',I2,2X,20A4)
  932 FORMAT(/1H ,50X,'(*)=OBSERVED, (X)=PREDICTED, (+)= IF * AND X OR'
     A,' Y COINSIDE.')
  933 FORMAT(1H ,50X,'(Y)=LONG RANGE FORECASTING')
cc  934 FORMAT(1H ,3X,E20.5,2(34X,E20.5))
  934 FORMAT(/1H ,3X,D20.5,2(34X,D20.5))
  935 FORMAT(1H ,10X,2H++,2(59(1H-),1H+))
  936 FORMAT(1H ,' N =',I4,2X,1HI,121A1)
cc  937 FORMAT(1H ,14X,E20.5)
  937 FORMAT(1H ,14X,D20.5)
 9370 FORMAT(1H ,'   N=',I5,5(4X,D20.5))
 9371 FORMAT(1H ,10X,5(4X,D20.5))
 1357 FORMAT(/1H ,'MATRIX S')
 8121 FORMAT(/1H ,'MATRIX  A(',I3,')')
      RETURN
      END
C
C
      SUBROUTINE MAXMIN(DMAX,DMIN,YD)
	REAL*8 DMAX,DMIN,YD
cc 1014 IF(DMAX-YD) 1001,1002,1002
 1014 IF(DMAX-YD.GE.0) GO TO 1002
 1001 DMAX=YD
cc 1002 IF(DMIN-YD) 1004,1004,1003
 1002 IF(DMIN-YD.LE.0) GO TO 1004
 1003 DMIN=YD
 1004 RETURN
      END
C
C
      SUBROUTINE SBSCAL(YA,YST,IX)
	REAL*8 YA,YST
C     YA SCALING
      IF(YA.LT.0) GO TO 11
   10 IFC=YA/YST+0.5
      GO TO 12
   11 IFC=YA/YST-0.5
   12 IX=IFC+61
      IF(IX.GT.121) IX=121
      IF(IX.LE.0) IX=1
      RETURN
      END