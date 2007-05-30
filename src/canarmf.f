      SUBROUTINE CANARMF(N,LAGH3,CYY,COEF,IFPL1,SD,AIC,OAIC,MO,A,
     *NC,MM1,MM2,V,Z,Y,XX,NDT,X3,X3MIN,MIN3,M1M,BETA,M1N,ALPHA,TMP,
     *MJ1,MJ2)
C
      INCLUDE 'timsac_f.h'
C
cc	PROGRAM CANARM
C     PROGRAM 74.1.1. CANONICAL CORRELATION ANALYSIS OF SCALAR TIME SERI
C-----------------------------------------------------------------------
C     ** DESIGNED BY H. AKAIKE, THE INSTITUTE OF STATISTICAL MATHEMATICS
C     ** PROGRAMMED BY E. ARAHATA, THE INSTITUTE OF STATISTICAL MATHEMAT
C	  TOKYO
C     ** DATE OF THE LATEST REVISION: MARCH 25, 1977
C     ** THIS PROGRAM WAS ORIGINALLY PUBLISHED IN
C	 "TIMSAC-74 A TIME SERIES ANALYSIS AND CONTROL PROGRAM PACKAGE(1
C	 BY H. AKAIKE, E. ARAHATA AND T. OZAKI, COMPUTER SCIENCE MONOGRA
C	 NO.5, MARCH 1975, THE INSTITUTE OF STATISTICAL MATHEMATICS
C     ** FOR THE BASIC THEORY SEE "CANONICAL CORRELATION ANALYSIS OF TIM
C	 AND THE USE OF AN INFORMATION CRITERION" BY H. AKAIKE, IN
C	 "SYSTEM IDENTIFICATION: ADVANCES AND CASE STUDIES" R. K. MEHRA
C	 D. G. LAINIOTIS EDS. ACADEMIC PRESS, NEW YORK, 1976
C-----------------------------------------------------------------------
C     THIS PROGRAM FITS AN AR-MA MODEL TO STATIONARY SCALAR TIME SERIES
C     THROUGH THE ANALYSIS OF CANONICAL CORRELATIONS
C     BETWEEN THE FUTURE AND PAST SETS OF OBSERVATIONS.
C     THE OUTPUTS OF THIS PROGRAM SHOULD BE ADDED TO THE INPUTS
C     TO THIS PROGRAM TO FORM AN INPUT TO THE PROGRAM AUTARM.
C
C     INPUTS REQUIRED:
C     (N,LAGH0): N, LENGTH OF ORIGINAL DATA Y(I) (I=1,N)
C		 LAGH0, MAXIMUM LAG OF COVARIANCE
C     CYY(I),I=0,LAGH0: AUTOCOVARIANCE SEQUENCE OF Y(I)
C
C     OUTPUTS:
C     NEWL: NEWL=1, FOR DIRECT INPUT TO PROGRAM AUTARM
C     M1M: ORDER OF AR
C     BETA(I)(I=1,M1M): AR-COEFFICIENTS
C     M1N: ORDER OF MA (=M1M-1)
C     ALPHA(I)(I=1,M1N): MA-COEFFICIENTS
C
C     THE AR-MA MODEL IS GIVEN BY
C     Y(N)+BETA(1)Y(N-1)+...+BETA(M1M)Y(N-M1M) = X(N)+ALPHA(1)X(N-1)+...
C						 ...+ALPHA(M1N)X(N-M1N)
C
cc      !DEC$ ATTRIBUTES DLLEXPORT :: CANARMF
C
cc      PARAMETER (MJ1=50,MJ2=101)
      IMPLICIT REAL*8(A-H,O-Z)
cc	DIMENSION X(50),Y(50),Z(50)
cc	DIMENSION WL(50)
cc	DIMENSION XX(50),X3(50),BETA(50)
cc	DIMENSION CYY(1001)
cc	DIMENSION A(101),ALPHA(50)
cc	DIMENSION VC(101),VT(101)
cc	DIMENSION ST(101,50),T(101,50),V(50,50)
cc	COMMON /COM9/AST1
cc	DIMENSION AST1(5200)
cc	DIMENSION VV(50,50)
c
      DIMENSION CYY(LAGH3)
      DIMENSION COEF(MJ2)
      DIMENSION SD(0:MJ1),AIC(0:MJ1),A(MJ1)
      DIMENSION MM1(MJ1),MM2(MJ1),V(MJ1,MJ1,MJ1)
      DIMENSION Z(MJ1,MJ1),Y(MJ1,MJ1),XX(MJ1,MJ1)
      DIMENSION NDT(MJ1,MJ1),X3(MJ1,MJ1)
      DIMENSION X3MIN(MJ1),MIN3(MJ1)
      DIMENSION BETA(MJ1),ALPHA(MJ1)
c
      DIMENSION WL(MJ1)
      DIMENSION VC(MJ2),VT(MJ2)
      DIMENSION ST(MJ2,MJ1),T(MJ2,MJ1)
      DIMENSION AST1((MJ2-1)*MJ2/2)
      DIMENSION VV(MJ1,MJ1)
c
      INTEGER*1 TMP(1)
      CHARACTER CNAME*80
c
C     INITIAL CLEARING
cc	DATA BETA /50*0.0D-00/, ALPHA/50*0.0D-00/
cc	DO 100 I=1,50
cc	   BETA(I)=0.0D0
cc	   ALPHA(I)=0.0D0
cc  100 CONTINUE
C
C     INPUT / OUTPUT DATA FILE OPEN
cc	CHARACTER(100) DFNAM
cc	CALL SETWND
cc	DFNAM='canarm.out'
cc	CALL FLOPN3(DFNAM,NFL)
cc	IF (NFL.EQ.0) GO TO 999
cc	IF (NFL.EQ.0) GO TO 999
c
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
            WRITE(*,*) ' ***  canarm temp FILE OPEN ERROR :',CNAME,IVAR
            IFG=0
         END IF
      END IF
C
C     ABSOLUTE DIMENSIONS USED FOR SUBROUTINE CALL
cc	MJ1=50
cc	MJ2=101
      CST0=0.0D-00
      CST1=1.0D-00
      CST2=2.0D-00
      CST9=9999.0D-00
cc	DO 18 I=1,1001
cc   18 CYY(I)=CST0
C     INITIAL CONDITION INPUT AND PRINT OUT
cc	READ(5,1) N,LAGH0
cc	LAGH3=LAGH0+1
      LAGH0=LAGH3-1
cc	READ(5,2) (CYY(I),I=1,LAGH3)
cc	WRITE(6,11130)
cc	WRITE(6,11111) N,LAGH0
      AN=N
      IFPL=3.0D-00*DSQRT(AN)
cc	IFPL=MIN0(IFPL,49,LAGH0)
      IFPL=MIN0(IFPL,(MJ1-1),LAGH0)
      IFPL1=IFPL+1
C     INITIAL AUTOREGRESSIVE MODEL FITTING
cc	CALL NSICP(CYY,IFPL1,N,A,MO,SD,AIC)
      NA=(IFPL1*(IFPL1+1)/2)
      CALL NSICP(CYY,LAGH3,IFPL1,N,AST1,NA,COEF,
     *		 SD,AIC,A,MO,OAIC,IFG,LU)
      NINEW=MO
      M9=MO+1
C     NINEW0=0
      INDX=1
C     MATRIX MULTIPLICATION: ORTHO-NORMALIZATION OF VARIABLES BY USING
C     AR-MODELS OF SUCCESSIVELY INCREASING ORDER
C     VC=AST1*CYY COMPUTATION
      DO 210 I=1,M9
  210 VC(I)=AST1(1)*CYY(I)
C     VT=VC*AST1' COMPUTATION
C     AST1 OBTAINED BY SUBROUTINE NSICP
cc	CALL SVCMAT(VC,VT,M9)
      CALL SVCMAT(VC,VT,M9,AST1,NA)
      DO 220 I=1,M9
  220 ST(I,1)=VT(I)
C     NINEW0: NUMBER OF ITERATIONS
      NINEW0=1
      NC=0
  500 M1=NINEW0+1
      NC=NC+1
      M2=M9
      M=M1+M2
C     VC=(M1-TH ROW OF AST1)*CYY
cc	CALL SVECT(CYY,VC,M9,M1,INDX)
cc	CALL SVCMAT(VC,VT,M9)
      CALL SVECT(CYY,LAGH3,AST1,NA,VC,M9,M1,INDX)
      CALL SVCMAT(VC,VT,M9,AST1,NA)
      DO 240 I=1,M9
  240 ST(I,M1)=VT(I)
      DO 250 I=1,M9
      DO 250 J=1,M1
  250 T(I,J)=ST(I,J)
C     SVD OF T
C     T IS THE COVARIANCE MATRIX BETWEEN THE SETS OF THE
C     ORTHO-NORMALIZED FUTURE AND PAST VARIABLES
C     SINGULAR VALUES (Z) ARE THE CANONICAL CORRELATION COEFFICIENTS.
C     COMMON SUBROUTINE CALL
cc	CALL MSVD(T,VV,Z,M2,M1,MJ2,MJ1)
cc	CALL SVTR(VV,V,M1,MJ1)
      CALL MSVD(T,VV,Z(1,NC),M2,M1,MJ2,MJ1)
      CALL SVTR(VV,V(1,1,NC),AST1,NA,M1,MJ1)
      DO 260 J=1,M1
cc  260 Y(J)=Z(J)*Z(J)
  260 Y(J,NC)=Z(J,NC)*Z(J,NC)
C     FUTURE CANONICAL WEIGHTS (V) PRINT OUT
cc	WRITE(6,6)
cc	WRITE(6,7) M1
cc	WRITE(6,8) M2
cc	WRITE(6,9) N
cc	WRITE(6,35)
C     COMMON SUBROUTINE CALL
cc	CALL SUBMPR(V,M1,M1,MJ1,MJ1)
      MM1(NC)=M1
      MM2(NC)=M2
C     TEST OF DEPENDENCE OF THE LAST PREDICTOR BY DIC (DIFFERENCE OF AIC
C     DIC(J) = AIC(J) - AIC(MAXIMUM J)
      EM=M
      EN=N
      J=M1
      WL(J+1)=CST1
cc   42 WL(J)=WL(J+1)*(CST1-Y(J))
   42 WL(J)=WL(J+1)*(CST1-Y(J,NC))
      J=J-1
      IF(J.GT.0) GO TO 42
      DO 45 J=1,M1
      IF(WL(J).GT.CST0) GO TO 145
cc	X(J)=CST9
cc	XX(J)=CST9
      XX(J,NC)=CST9
      GO TO 45
cc  145 X(J)=-EN*DLOG(WL(J))
cc	XX(J)=-EN*DLOG(WL(J))
  145 XX(J,NC)=-EN*DLOG(WL(J))
   45 CONTINUE
cc	NDT=M1*M2
cc	ANDT=NDT
      NDT(1,NC)=M1*M2
      ANDT=NDT(1,NC)
C     DIC(J)=X3(J)
cc	X3(1)=XX(1)-CST2*ANDT
      X3(1,NC)=XX(1,NC)-CST2*ANDT
cc	WRITE(6,49)
      J=0
cc	WRITE(6,50) J,Z(1),Y(1),XX(1),NDT,X3(1)
cc	X3MIN=X3(1)
cc	MIN3=0
      X3MIN(NC)=X3(1,NC)
      MIN3(NC)=0
      IF(M1.LT.2) GO TO 4110
      DO 51 J=2,M1
      J1=J-1
cc	NDT=(M1-J1)*(M2-J1)
cc	ANDT=NDT
cc	X3(J)=XX(J)-CST2*ANDT
      NDT(J,NC)=(M1-J1)*(M2-J1)
      ANDT=NDT(J,NC)
      X3(J,NC)=XX(J,NC)-CST2*ANDT
cc   51 WRITE(6,50) J1,Z(J),Y(J),XX(J),NDT,X3(J)
   51 CONTINUE
      DO 4300 J=2,M1
C     MINIMUM OF DIC SERCH
cc	IF(X3(J).GE.X3MIN) GO TO 4300
cc	X3MIN=X3(J)
cc	MIN3=J-1
      IF(X3(J,NC).GE.X3MIN(NC)) GO TO 4300
      X3MIN(NC)=X3(J,NC)
      MIN3(NC)=J-1
 4300 CONTINUE
cc 4110 WRITE(6,4410) X3MIN,MIN3
 4110 CONTINUE
cc	WRITE(6,11112)
      IF(NINEW0.EQ.NINEW) GO TO 6999
C     DEPENDENCE ACCEPTED WHEN M1N-DIC IS NEGATIVE
cc	IF(X3(M1).GT.CST0) GO TO 110
      IF(X3(M1,NC).GT.CST0) GO TO 110
 6999 M1M=M1-1
      IF(M1M.LE.0) GO TO 110
C     BETA(AR-COEFF) COMPUTATION
cc	AII=CST1/V(M1,M1)
      AII=CST1/V(M1,M1,NC)
      DO 5100 I=1,M1M
      II=M1-I
cc	BETA(II)=V(M1,I)*AII
      BETA(II)=V(M1,I,NC)*AII
 5100 CONTINUE
C     ALPHA(MA-COEFF) COMPUTATION
cc	CALL ALPHAS(A,M1M,BETA,ALPHA)
      CALL ALPHAS(COEF,M1M,BETA,ALPHA)
C     THE INPUTS TO THE PROGRAM AUTAMA PUNCH OUT
 5400 NEWL=1
cc	WRITE(7,1) NEWL
C     BETA(AR-COEFF), ALPHA(MA-COEFF) PUNCH OUT
cc	WRITE(7,1) M1M
      IF  (M1M.LE.0) GO TO 5200
cc	WRITE(7,2) (BETA(I),I=1,M1M)
 5200 M1N=M1M-1
cc	WRITE(7,1) M1N
      IF  (M1N.LE.0) GO TO 1100
cc	WRITE(7,2) (ALPHA(I),I=1,M1N)
      GO TO 1100
  110 IF(NINEW0.LT.NINEW) GO TO 5300
      M1M=0
      GO TO 5400
 5300 NINEW0=NINEW0+1
      GO TO 500
cc 1100 CALL FLCLS3(NFL)
cc  999 STOP
 1100 CONTINUE
      IF (IFG.NE.0) CLOSE(LU)
      RETURN
    1 FORMAT(16I5)
    2 FORMAT(4D20.10)
    6 FORMAT(//1H ,21HCANONICAL CORRELATION)
    7 FORMAT(/1H ,'NUMBER OF PRESENT AND FUTURE VARIABLES',2X,
     A	     'M1=',I5)
    8 FORMAT(/1H ,'NUMBER OF PRESENT AND PAST VARIABLES',4X,
     A	     'M2=',I5)
    9 FORMAT(/1H ,'DATA LENGTH=N=',I6)
   35 FORMAT(/1H ,'FUTURE SET CANONICAL WEIGHTS, ROWWISE')
   49 FORMAT(/1H ,5X,8HORDER(P),4X,11HCANONICAL R,3X,9HR-SQUARED,3X,
     A 10HCHI-SQUARE,3X,6HN.D.F.,2X,23HDIC (P)(=CHI**2-2*D.F.))
   50 FORMAT(1H ,10X,I3,4X,F8.4,4X,F8.3,7X,F8.2,4X,I6,3X,F10.4)
 4410 FORMAT(/1H ,'MINIMUM DIC(P) =',F12.2,1X,'ATTAINED AT P=',I5)
 1130 FORMAT(/1H ,'STRUCTURAL CHARACTERISTIC VECTOR',
     A	     ' (H(I),I=1,P)')
11111 FORMAT(/1H ,'INITIAL AUTO REGRESSIVE MODEL FITTING ',
     A	     'BY THE MINIMUM AIC PROCEDURE. / N=',I5,',LAGH0=',I5)
11112 FORMAT(1H ,5X,'THE VALUES OF CHI-SQUARE AND DIC (P) ',
     A	     'CORRESPONDING TO CANONICAL R=','1.000 ',
     A	     'SHOULD BE IGNORED')
11130 FORMAT(/1H ,'PROGRAM 74.1.1. CANARM')
      END
C
      SUBROUTINE ALPHAS(A,M1M,BETA,ALPHA)
C     THIS SUBROUTINE COMPUTES ALPHA(MA-COEFFICIENTS).
      IMPLICIT REAL*8(A-H,O-Z)
cc	DIMENSION A(101)
      DIMENSION A(M1M)
      DIMENSION BETA(M1M),ALPHA(M1M)
      ALPHA(M1M)=0.0D-00
      IF  (M1M.LE.1) GO TO 20
      ALPHA(1)=BETA(1)-A(1)
      IF(M1M.LE.2) GO TO 20
      IPM=M1M-1
      DO 10 K=2,IPM
      KM1=K-1
      SUM=0.0
      DO 11 I=1,KM1
      KMI=K-I
   11 SUM=SUM-ALPHA(I)*A(KMI)
      ALPHA(K)=BETA(K)-A(K)+SUM
   10 CONTINUE
C     ALPHA,BETA PRINT OUT
cc   20 WRITE(6,60) M1M
   20 CONTINUE
cc	WRITE(6,61)
cc	DO  21 I=1,M1M
cc   21 WRITE(6,62) I,BETA(I),ALPHA(I)
      RETURN
   60 FORMAT(//1H ,'M1M=',I5)
   61 FORMAT(/1H ,4X,1HI,6X,17HBETA(I) (AR-COEF),4X,
     A 18HALPHA(I) (MA-COEF))
   62 FORMAT(1H ,I5,6X,F17.5,5X,F17.5)
      END
C
C
cc	SUBROUTINE NSICP(CYY,L1,N,COEF,MO,OSD,OAIC)
      SUBROUTINE NSICP(CYY,L3,L1,N,AST1,NA,COEF,SD,AIC,AA,
     *MO,OAIC,IFG,LU)
C     COMMON SUBROUTINE
C     THIS SUBROUTINE FITS AUTOREGRESSIVE MODELS OF SUCCESSIVELY
C     INCREASING ORDER UP TO L(=L1-1).
C     INPUT:
C     CYY(I),I=0,L1; AUTOCOVARIANCE SEQUENCE
C     L1: L1=L+1, L IS THE UPPER LIMIT OF THE MODEL ORDER
C     N; LENGTH OF ORIGINAL DATA
C     OUT PUT:
C     COEF; AR-COEFFICIENTS
C     MO: ORDER OF AR
C     OSD: INNOVATION VARIANCE
C     OAIC: VALUE OF AIC
C     AST1: MATRIX OF AR-COEFFICIENTS (IN VECTOR FORM)
      IMPLICIT REAL*8(A-H,O-Z)
cc	COMMON /COM9/AST1
cc	DIMENSION AST1(5200)
cc	DIMENSION A(101),B(101)
      DIMENSION CYY(L3),COEF(L1)
      DIMENSION AST1(NA)
      DIMENSION A(L1),B(L1)
      DIMENSION AIC(0:L1),SD(0:L1),AA(L1)
cc      REAL*4 AX,BL,STA,DASH,PLUS
cc      REAL*4 FFFF
cc      REAL*4  F(41) / 41*1H  /, AMES(41) / 41*1H- /
      CHARACTER AX,BL,STA,DASH,PLUS
      CHARACTER FFFF
      CHARACTER F(41),AMES(41)
cc      DATA AX,BL,STA,DASH,PLUS/1H!,1H ,1H*,1H-,1H+/
      DATA AX,BL,STA,DASH,PLUS/'!',' ','*','-','+'/
      DATA F,AMES/ 41*' ', 41*'-' /
      CST0=0.0D-00
      CST1=1.0D-00
      CST2=2.0D-00
      CST20=20.0D-00
      CST05=0.05D-00
      CST01=0.00001D-00
      L=L1-1
cc	SD=CYY(1)
      SD(1)=CYY(1)
      INX=1
cc	AST1(1)=CST1/DSQRT(SD)
      AST1(1)=CST1/DSQRT(SD(1))
      AN=N
cc	OAIC=AN*DLOG(SD)
cc	OSD=SD
      OAIC=AN*DLOG(SD(1))
      MO=0
cc	SD0=OSD
cc	AIC0=OAIC
      SD(0)=SD(1)
	AIC(0)=OAIC
C     INITIAL CONDITION PRINT OUT
cc  991 WRITE(6,1100)
  991 CONTINUE
      RAN=CST1/DSQRT(AN)
      SCALH=CST20
      JJ0=SCALH+CST1
      JJL=SCALH*CST2+CST1
      JJL1=JJL-1
      AMES(1)=PLUS
      AMES(11)=PLUS
      AMES(JJ0)=PLUS
      AMES(JJ0+10)=PLUS
      AMES(JJL)=PLUS
      IAN=SCALH*(RAN+CST05)
      IAN1=IAN+JJ0
      IAN2=2*IAN+JJ0
      LAN1=-IAN+JJ0
      LAN2=-2*IAN+JJ0
cc	WRITE(6,26100)
cc	WRITE(6,26101)
cc	WRITE(6,261)
cc	WRITE(6,26102)
cc	WRITE(6,262)
cc	WRITE(6,264) (AMES(J),J=1,JJL)
      IF (IFG.NE.0) THEN
	 WRITE(LU,261)
	 WRITE(LU,262)
	 WRITE(LU,264) (AMES(J),J=1,JJL)
      END IF
cc	WRITE(6,859) MO,OSD,OAIC
      F(JJ0)=AX
      F(IAN1)=AX
      F(IAN2)=AX
      F(LAN1)=AX
      F(LAN2)=AX
cc	WRITE(6,861) (F(J),J=1,JJL)
      IF (IFG.NE.0) WRITE(LU,264) (F(J),J=1,JJL)
      SE=CYY(2)
C     ITERATION START
      DO 400 M=1,L
cc	SDR=SD/CYY(1)
      SDR=SD(M)/CYY(1)
      IF(SDR.GE.CST01) GO TO 399
cc	WRITE(6,2600)
      GO TO 402
  399 MP1=M+1
cc	D=SE/SD
      D=SE/SD(M)
      A(M)=D
      D2=D*D
cc	SD=(CST1-D2)*SD
cc	CONST=CST1/DSQRT(SD)
      SD(M)=(CST1-D2)*SD(M)
      CONST=CST1/DSQRT(SD(M))
      AM=M
cc	AIC=AN*DLOG(SD)+CST2*AM
      AIC(M)=AN*DLOG(SD(M))+CST2*AM
C
C
cc	DLSD=DLOG(SD)
      DLSD=DLOG(SD(M))
      IF(M.EQ.1) GO TO 410
C     A(I) COMPUTATION
      LM=M-1
      DO 420 I=1,LM
      A(I)=A(I)-D*B(I)
  420 CONTINUE
  410 DO 460 I=1,M
      II=M+1-I
      INX=INX+1
  460 AST1(INX)=-A(II)*CONST
      INX=INX+1
      AST1(INX)=CONST
      DO 421 I=1,M
      IM=MP1-I
  421 B(I)=A(IM)
C     M,SD,AIC	PRINT OUT
      IF(A(M).LT.CST0) GO TO 300
      NFC=SCALH*(A(M)+CST05)
      GO TO 310
  300 NFC=SCALH*(A(M)-CST05)
  310 ANFC=NFC
      JJ=ANFC+SCALH+CST1
      FFFF=F(JJ)
      F(JJ)=STA
cc	WRITE(6,860) M,SD,AIC,A(M)
      AA(M)=A(M)
cc	WRITE(6,861) (F(J),J=1,JJL)
      IF (IFG.NE.0) WRITE(LU,264) (F(J),J=1,JJL)
      F(JJ)=FFFF
cc  990 IF(OAIC.LT.AIC) GO TO 440
cc	OAIC=AIC
  990 IF(OAIC.LT.AIC(M)) GO TO 440
      OAIC=AIC(M)
cc	OSD=SD
      MO=M
      DO 430 I=1,M
  430 COEF(I)=-A(I)
  440 IF(M.EQ.L) GO TO 400
      SE=CYY(M+2)
      DO 441 I=1,M
  441 SE=SE-B(I)*CYY(I+1)
      SD(M+1)=SD(M)
  400 CONTINUE
  402 CONTINUE
C     MO, COEF(I) OUT PUT
cc	WRITE(6,870) OAIC,MO
cc	WRITE(6,1871)
cc	WRITE(6,871)
cc	CALL SUBVCP(COEF,MO)
cc  699 F(JJ0)=BL
  699 CONTINUE
      F(IAN1)=BL
      F(IAN2)=BL
      F(LAN1)=BL
      F(LAN2)=BL
      AMES(JJ0)=DASH
      AMES(JJ0+10)=DASH
      AMES(JJL)=DASH
      RETURN
26100 FORMAT(/1H ,16X,'SD(M)',15X,'AIC(M)',13X,'A(M)')
26101 FORMAT(1H ,4X,'M',11X,'INNOVATION',10X,'AIC(M)=    ',8X,
     A	    'PARTIAL AUTO-')
cc  261 FORMAT(1H ,16X,'   VARIANCE',9X,'N*DLOG(SD(M))+2*M',
cc     A      '  CORRELATION    ',
cc     A      'PARTIAL CORRELATION (LINES SHOW +SD AND +2SD)')
  261 FORMAT(' PARTIAL CORRELATION (LINES SHOW +/-SD AND +/-2SD)')
26102 FORMAT(1H ,102X,'_',7X,'_')
cc  262 FORMAT(1H ,69X,'-1',19X,'0',19X,'1')
cc  264 FORMAT(1H ,70X,41A1)
  262 FORMAT(' -1',19X,'0',19X,'1')
  264 FORMAT(2X,41A1)
  859 FORMAT(1H ,I5,2X,2D20.5)
  860 FORMAT(1H ,I5,2X,3D20.5)
  861 FORMAT(1H ,70X,41A1)
  960 FORMAT(/1H ,5X,'A(I)')
  870 FORMAT(/1H ,'MINIMUM AIC(M)=',D12.5,2X,'ATTAINED AT M=',I5)
  980 FORMAT(/1H ,5X,'COEF(I)')
 1100 FORMAT(/1H ,'AIC(M)=N*DLOG(SD)+2*M')
  871 FORMAT(/1H ,'AR-COEFFICIENTS')
 1871 FORMAT(/1H ,'AR MODEL: Y(N)+AR(1)Y(N-1)+...+AR(M)Y(N-M)=X(N)')
 2600 FORMAT(/1H ,'ACCURACY OF COMPUTATION LOST')
      END
C
cc	SUBROUTINE SVCMAT(VC,VT,M9)
      SUBROUTINE SVCMAT(VC,VT,M9,AST1,NA)
C     THIS SUBROUTINE COMPUTES VT=VC*AST1'.
C     AST1 IS AN OUTPUT OF NSICP.
      IMPLICIT REAL*8(A-H,O-Z)
cc	COMMON /COM9/AST1
cc	DIMENSION AST1(5200)
      DIMENSION AST1(NA)
      DIMENSION VC(M9),VT(M9)
      CST0=0.0D-00
      INX=0
      DO  10 I=1,M9
      SUM=CST0
      DO 11 K=1,I
      INX=INX+1
   11 SUM=SUM+VC(K)*AST1(INX)
      VT(I)=SUM
   10 CONTINUE
      RETURN
      END
C
      SUBROUTINE SVTR(VV,V,AST1,NA,M1,MJ1)
C     THIS SUBROUTINE COMPUTES FUTURE SET CANONICAL WEIGHTS DEFINED BY
C     V=VV'*AST1.
C     AST1 IS AN OUTPUT OF NSICP.
      IMPLICIT REAL*8(A-H,O-Z)
cc	COMMON /COM9/AST1
cc	DIMENSION AST1(5200),ISUM1(50)
      DIMENSION VV(MJ1,MJ1)
      DIMENSION AST1(NA),ISUM1(M1)
      DIMENSION V(MJ1,MJ1)
      CST0=0.0D-00
      ISUM=0
      DO 15 I=1,M1
      ISUM=ISUM+I
   15 ISUM1(I)=ISUM
      DO 10 I=1,M1
      DO 11 J=1,M1
      SUM=CST0
      LL=ISUM1(J)
      INX=0
      DO 12 K=J,M1
      IJK=LL+INX
      SUM=SUM+VV(K,I)*AST1(IJK)
      INX=INX+K
   12 CONTINUE
      V(I,J)=SUM
   11 CONTINUE
   10 CONTINUE
      RETURN
      END
C
      SUBROUTINE SVECT(CYY,L3,AST1,NA,VC,M9,M1,INDX)
C     THIS SUBROUTINE COMPUTES VC=(M1-TH ROW OF AST1)*(CYY MATRIX)
C     AST1 IS AN OUTPUT OF NSICP.
      IMPLICIT REAL*8(A-H,O-Z)
cc	COMMON /COM9/AST1
cc	DIMENSION AST1(5200)
cc	DIMENSION CYY(1001),VC(M9)
      DIMENSION AST1(NA)
      DIMENSION CYY(L3),VC(M9)
      CST0=0.0D-00
      DO  10 I=1,M9
   10 VC(I)=CST0
      DO 20 IS=1,M1
      INDX=INDX+1
      ISM=IS-1
      DO  30 I=1,M9
      II=ISM+I
      VC(I)=VC(I)+AST1(INDX)*CYY(II)
   30 CONTINUE
   20 CONTINUE
      RETURN
      END
C