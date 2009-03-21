      SUBROUTINE CANOCAF(IR,INW,N,LAGH1,IP0,CCV,L,AIC,OAIC,MO,OSD,AAO,
     *NC,N1,N2,VV,Z,Y,XX,NDT,X3,X3MIN,MIN3,F,M1NH,NH,G,IAW,VF,
     *LMAX,MJ0,MJ1)
C
      INCLUDE 'timsac_f.h'
C
cc	PROGRAM CANOCA
C     PROGRAM 74.2.1. CANONICAL CORRELATION ANALYSIS OF VECTOR TIME SERI
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
C     THIS PROGRAM DOES CANONICAL CORRELATION ANALYSIS OF AN IR-DIMENSIO
C     MULTIVARIATE TIME SERIES Y(I) (I=1,N).
C
C     FIRST AR-MODEL IS FITTED BY THE MINIMUM  A I C  PROCEDURE.
C     THE RESULTS ARE USED TO ORTHO-NORMALIZE THE PRESENT AND PAST VARIA
C     THE PRESENT AND FUTURE VARIABLES ARE TESTED SUCCESSIVELY TO DECIDE
C     ON THE DEPENDENCE OF THEIR PREDICTORS. WHEN THE LAST DIC (AN INFOR
C     CRITERION) IS NEGATIVE THE PREDICTOR OF THE VARIABLE IS DECIDED
C     TO BE LINEARLY DEPENDENT ON THE ANTECEDENTS.
C
C     THE STRUCTURAL CHARACTERISTIC VECTOR H OF THE CANONICAL MARKOVIAN
C     REPRESENTATION AND THE ESTIMATE OF THE TRANSITION MATRIX F, IN
C     VECTOR FORM, ARE PUNCHED OUT. THE ESTIMATE OF THE INPUT MATRIX G A
C     THE COVARIANCE MATRIX C OF THE INNOVATION, OBTAINED BY USING
C     THE F-MATRIX AND THE AR-MODEL, ARE ALSO PUNCHED OUT.
C
C     INPUTS REQUIRED:
C     IR:		 DIMENSION OF Y(I)
C     INW(K)(K=1,IP):	 INW(K)=J MEANS THAT THE K-TH COMPONENT OF Y(I)
C	  IS THE J-TH COMPONENT OF THE ORIGINAL RECORD Z(I) USED FOR
C	  THE COMPUTATION OF THE COVARIANCE SEQUENCE CZZ(I).
C
C     ALSO THE FOLLOWING INPUTS ARE REQUESTED BY THE SUBROUTINE	 R E C O
C     THESE INPUTS CAN BE OBTAINED AS THE OUTPUTS OF THE TIMSAC
C     PROGRAM 5.1.2 ( MULCOR ).
C     N:		 DATA LENGTH
C     LAGH:		 MAXIMUM LAG OF COVARIANCE
C     IP0:		 DIMENSION OF THE ORIGINAL RECORD
C     CZZ(I)(I=0,LAGH):	 COVARIANCE MATRIX SEQUENCE OF Z(I) GIVEN IN THE
C	  OF THE SUCCESSION OF THE SEQUENCES OF COVARIANCES BETWEEN THE
C	  C-TH COMPONENTS OF Z(I), ARRANGED IN THE ORDER (R=1,C=1),(R=2,
C	  (R=2,C=1),(R=1,C=2),(R=3,C=3),(R=3,C=1),(R=1,C=3),(R=3,C=2),(R
C	  ...... EACH SEQUENCE HAS THE HEADING (R,C).
C	  THE (J+1)ST ELEMENT OF THE COVARIANCE SEQUENCE WITH R=K AND C=
C	  IS AN ESTIMATE OF E(Z(I+J,K)*Z(I,L)), WHERE Z(I,K) DENOTES
C	  THE K-TH COMPONENT OF Z(I).
C
cc      !DEC$ ATTRIBUTES DLLEXPORT :: CANOCAF
C
      IMPLICIT REAL*8(A-H,O-Z)
cc	DIMENSION Y(91),Z(91),WL(91)
cc	DIMENSION CV(25,7,7)
cc	DIMENSION XX(91),X3(91),F(91),VF(637)
cc	DIMENSION U(91,46),V(46,46),VV(46,46)
cc	DIMENSION AST1(91,7,7)
cc	DIMENSION NH(91),IH(91)
cc	DIMENSION FL(46,46)
cc	DIMENSION RGT(91,91),FRG(46,91)
cc	DIMENSION VTG(91),SR(46),ZZ(46),SF(46),SFRG(91)
cc	DIMENSION INW(7),C1(7,7)
cc	DIMENSION OSD(7,7),AO(13,7,7)
cc	DIMENSION G(46,7)
      DIMENSION Y(MJ1,MJ1),Z(MJ1,MJ1),WL(MJ1)
      DIMENSION CV(LMAX*2+1,IR,IR)
      DIMENSION XX(MJ1,MJ1),X3(MJ1,MJ1),F(MJ1,MJ1),VF(MJ1*MJ1)
      DIMENSION U(MJ1,MJ1),V(MJ1,MJ1),VV(MJ1,MJ1,MJ1)
      DIMENSION AST1((LMAX+1)*LMAX,IR,IR)
      DIMENSION NH(MJ1),IH(MJ1)
      DIMENSION FL(MJ1,MJ1)
      DIMENSION RGT(MJ1,MJ1),FRG(MJ1,MJ1)
      DIMENSION VTG(MJ1),SR(MJ1),ZZ(MJ1),SF(MJ1),SFRG(MJ1)
      DIMENSION INW(IR),C1(IP0,IP0)
      DIMENSION OSD(IR,IR),AO(MJ0,IR,IR),AAO(MJ0,IR,IR)
      DIMENSION G(MJ1,IR)
C
      DIMENSION CCV(LAGH1,IP0,IP0)
      DIMENSION AIC(MJ0)
      DIMENSION N1(MJ1),N2(MJ1),NDT(MJ1,MJ1)
      DIMENSION X3MIN(MJ0*IR),MIN3(MJ0*IR)
C
cc	COMMON /COM9/AST1
cc	COMMON /COM10/CV
cc	COMMON /COM11/RGT
cc	COMMON /COM19/U
cc	COMMON /COM20/V
cc	COMMON /COM21/VV
cc	COMMON /COM25/FL
cc      EQUIVALENCE (AST1(1,1,1),VF(1))
C
C     INPUT / OUTPUT DATA FILE OPEN
cc	CHARACTER(100) DFNAM
cc	CALL SETWND
cc	DFNAM='canoca.out'
cc	CALL FLOPN3(DFNAM,NFL)
cc	IF (NFL.EQ.0) GO TO 999
C
cc	MJ1=46
cc	MJ2=91
cc	MJ=7
cc	MJ0=13
      NC = 0
      CST0=0.0D-00
      CST1=1.0D-00
      CST2=2.0D-00
cc	DO 3111 I=1,46
cc	DO 3112 J=1,7
      DO 3111 I=1,MJ1
      DO 3112 J=1,IR
      G(I,J)=CST0
 3112 CONTINUE
 3111 CONTINUE
C     INITIAL CONDITION INPUT
cc	READ(5,1) IR
      IL=0
      IP=IR
cc	READ(5,1) (INW(I),I=1,IP)
cc	WRITE(6,11114)
cc	WRITE(6,11111)
cc	WRITE(6,11113) IR
cc	LCV=24
      LCV=LMAX*2
      LCV1=LCV+1
C     AUTO COVARIANCE INPUT
cc	CALL RECOVA(LCV1,IP0,N)
      DO 90 I=1,MIN(LCV1,LAGH1)
	DO 90 J=1,IP
	DO 90 K=1,IP
	 CV(I,J,K) = CCV(I,J,K)
   90 CONTINUE
C     L (UPPER BOUND OF AR-ORDER) COMPUTATION
cc	LMAX=12
      XL=2*N
      XL=DSQRT(XL)
      YL=IR
      J=XL/YL
      L= MIN0(LMAX,J)
      L1=L+1
C     MATRIX ARRANGEMENT
      DO 10 II=1,LCV1
      DO 20 I=1,IP0
      DO 20 J=1,IP0
   20 C1(I,J)=CV(II,I,J)
C     MATRIX REARRANGEMENT BY INW
C     COMMON SUBROUTINE CALL
cc	CALL REARRA(C1,INW,IP0,IP,MJ)
      CALL REARRA(C1,INW,IP0,IP)
      DO 21 I=1,IP
      DO 21 J=1,IP
   21 CV(II,I,J)=C1(I,J)
   10 CONTINUE
cc	WRITE(6,259) (INW(I),I=1,IP)
      IAW=0
      ID=IR
      LMAX2=(LMAX+1)*LMAX
C     AR-MODEL FITTING BY THE MINIMUM AIC PROCEDURE
C     COMMON SUBROUTINE CALL
cc	CALL NWFPEC(OSD,AO,L,IR,IL,N,MO,MJ0,MJ)
      CALL NWFPEC(AIC,OAIC,CV,AST1,OSD,AO,AAO,L,IR,IL,N,MO,MJ0,
     * LMAX2,LCV1)
C     MO: MAICE AR-MODEL ORDER DETERMINED BY NWFPEC
C     RGT=R12*G' COMPUTATION
cc	CALL SBRUGT(MO,ID)
      CALL SBRUGT(MO,ID,AST1,CV,RGT,MJ1,IR,LMAX2,LCV1)
      MMMH=(MO+1)*ID
      NINEW=MMMH
C     FRG' COMPUTATION
C     M1: NUMBER OF VARIABLE IN THE FUTURE SET
C     M2: NUMBER OF VARIABLES IN THE PAST SET =MMMH
      M1=1
      R11C=CV(1,1,1)
c-----------
      DO 22 I=1,MJ1
         DO 22 J=1,MJ1
            FL(I,J)=CST0
   22 CONTINUE
c------------
      FL(1,1)=CST1/DSQRT(R11C)
      FCT=FL(1,1)
      M2=(MO+1)*ID
      DO 100 J=1,M2
      FRG(1,J)=FCT*RGT(1,J)
  100 CONTINUE
      ISW=0
      NINEW0=2
      DO 191 I=1,IR
  191 NH(I)=I
      DO 634 I=1,MMMH
  634 IH(I)=0
  500 M1=M1+1
      M2=MMMH
      M=M1+M2
      MP1=M1+1
C     COVARIANCE MATRIX ARRANGEMENT
      NINEWI=NINEW0
      DO 636 I=NINEWI,MMMH
      IF(IH(I).EQ.0) GO TO 637
C     TEST OF EXHAUSTION OF FUTURE VARIABLES
      IF(I.GE.MMMH) GO TO 1140
      NINEW0=NINEW0+1
      GO TO 636
  637 NH(M1)=I
      GO TO 638
  636 CONTINUE
  638 CONTINUE
C     THE AUGMENTED FRG' MATRIX (F+)(G+)G' IS OBTAINED BY THE
C     BORDERING TECHNIQUE:
C			+-	 -+
C		 (F+) = I F    0  I
C			I	  I
C			I	  I
C			I SF'  SG I
C			+-	 -+,
C			+-	 -+
C	       (RVV+) = I RVV  SR I
C			I	  I
C			I	  I
C			I SR'  SS I
C			+-	 -+,
C	       (F+)(RVV+)(F+)'=I,
C			+-	 -+
C		 (R+) = I    R	  I
C			I	  I
C			I	  I
C			I    ST'  I
C			+-	 -+,
C	       AND
C			+-	 -+
C	   (F+)(R+)G' = I   FRG'  I
C			I	  I
C			I	  I
C			I  LRFRG  I
C			+-	 -+,
C	   FRG' IS STORED AS FRG.
C     SR,SS,VTG=(ST)'G' ARRANGEMENT
      II=NH(M1)
      IDB=II/ID
      IRES=II-IDB*ID
      IF(IRES.NE.0) IDB=IDB+1
      IF(IRES.EQ.0) IRES=ID
      DO 640 J=1,M2
  640 VTG(J)=RGT(II,J)
      DO 642 J=1,M1
      JJ=NH(J)
      JDB=JJ/ID
      JRES=JJ-JDB*ID
      IF(JRES.NE.0) JDB=JDB+1
      IF(JRES.EQ.0) JRES=ID
      IJDB=IDB-JDB
      IF(IJDB.LT.0) GO TO 646
      SR(J)=CV(IJDB+1,IRES,JRES)
      GO TO 642
  646 IJDB1=-IJDB+1
      SR(J)=CV(IJDB1,JRES,IRES)
  642 CONTINUE
      SS=SR(M1)
C     ZZ= F*SR COMPUTATION (F STORED AS FL IN COMMON AREA)
      MM1=M1-1
cc	CALL BLMLVC(SR,ZZ,MM1)
      CALL BLMLVC(SR,ZZ,MM1,FL,MJ1)
C     INNER PRODUCT OF FL*SR COMPUTATION
C     COMMON SUBROUTINE CALL
      CALL INNERP(ZZ,ZZ,FSRINP,MM1)
C     SG COMPUTATION
      SSFR2=SS-FSRINP
cc      IF  (SSFR2) 731,731,732
      IF  (SSFR2 .LT. 0) GO TO 731
      IF  (SSFR2 .EQ. 0) GO TO 731
      IF  (SSFR2 .GT. 0) GO TO 732
  731 SG=CST0
      GO TO 733
  732 SG=CST1/DSQRT(SSFR2)
C     SF=-(F'*(F*SR))*SG COMPUTATION
cc  733 CALL AVMLVC(ZZ,SF,MM1)
  733 CALL AVMLVC(ZZ,SF,MM1,FL,MJ1)
      DO 350 I=1,MM1
  350 SF(I)=-SF(I)*SG
C     F AUGMENTATION
      DO 360 J=1,MM1
  360 FL(M1,J)=SF(J)
      FL(M1,M1)=SG
C     SFRG=SF'*R12*G' COMPUTATION
cc	CALL VECMTX(SF,SFRG,NH,MM1,M2)
      CALL VECMTX(SF,SFRG,NH,RGT,MM1,M2,MJ1)
C      LRFRG=SFRG+SG*VTG (=THE LAST ROW OF (F+)(R+)G')
      DO 370 I=1,M2
  370 FRG(M1,I)=SFRG(I)+SG*VTG(I)
      DO 230 I=1,M2
      DO 230 J=1,M1
      U(I,J)=FRG(J,I)
  230 CONTINUE
C     CANONICAL WEIGHTS FOR THE SET OF ORTHO-NORMALIZED FUTURE VARIABLES
C     ARE OBTAINED IN V(COMMON AREA).
C     CANONICAL CORRELATION COEFFICIENTS ARE RETURNED IN Z.
C     COMMON SUBROUTINE CALL
cc	CALL MSVD(Z,M2,M1)
	NC = NC+1
      CALL MSVD(U,V,Z(1,NC),M2,M1,MJ1,MJ1)
C     W=V'*FL : CANONICAL WEIGHTS COMPUTATION
cc	CALL MWTFL(M1)
      CALL MWTFL(V,VV(1,1,NC),M1,FL,MJ1)
      DO 260 J=1,M1
cc  260 Y(J)=Z(J)*Z(J)
  260 Y(J,NC)=Z(J,NC)*Z(J,NC)
cc	WRITE(6,6)
cc	WRITE(6,7) M1
cc	WRITE(6,8) M2
      N1(NC)=M1
      N2(NC)=M2
cc	WRITE(6,9) N
cc	WRITE(6,35)
cc	DO 8100 I=1,M1
cc 8100 WRITE(6,8200) I,(VV(I,J),J=1,M1)
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
cc	XX(J)=9999.0D-00
      XX(J,NC)=9999.0D-00
      GO TO 45
cc  145 XX(J)=-EN*DLOG(WL(J))
  145 XX(J,NC)=-EN*DLOG(WL(J))
   45 CONTINUE
cc	NDT=M1*M2
cc	ANDT=NDT
      NDT(1,NC)=M1*M2
      ANDT=NDT(1,NC)
C     DIC(=CHI-SQUARE-2.0*D.F.) COMPUTATION
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
cc   51 WRITE(6,50) J1,Z(J),Y(J),XX(J),NDT,X3(J)
      NDT(J,NC)=(M1-J1)*(M2-J1)
      ANDT=NDT(J,NC)
      X3(J,NC)=XX(J,NC)-CST2*ANDT
   51 CONTINUE
      DO 4300 J=2,M1
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
      IF(X3(M1,NC).GT.CST0) GO TO 110
 6999 M1M=M1-1
      IF(M1M.LE.0) GO TO 110
C     TRANSITION MATRIX (F) COMPUTATION
cc	AII=CST1/VV(M1,M1)
      AII=CST1/VV(M1,M1,NC)
      DO 5100 I=1,M1M
      IAW=IAW+1
cc	F(I)=-VV(M1,I)*AII
cc	VF(IAW)=F(I)
      F(I,NC)=-VV(M1,I,NC)*AII
      VF(IAW)=F(I,NC)
 5100 CONTINUE
cc	WRITE(6,5200)
C     COMMON SUBROUTINE CALL
cc	CALL SUBVCP   (F,M1M)
      M1=M1-1
      ISW=ISW+1
      DO 1120 I=NINEW0,MMMH,ID
 1120 IH(I)=1
      IF(ISW.GE.ID) GO TO 1100
  110 IF(NINEW0.GE.NINEW) GO TO 1100
      NINEW0=NINEW0+1
      GO TO 500
 1140 M1=M1-1
cc 1100 WRITE(6,1130)
 1100 CONTINUE
C     STRUCTURAL CHARACTERISTIC VECTOR PRINT AND PUNCH OUT
      M1NH=M1
cc	WRITE(6,1131) (NH(I),I=1,M1)
cc	WRITE(7,1) IR,M1
cc	WRITE(7,1) (NH(I),I=1,M1)
C     F MATRIX (IN VECTOR FORM) PUNCH OUT
cc	WRITE(7,2) (VF(I),I=1,IAW)
C     INPUT MATRIX (G) COMPUTATION
cc 1011 CALL SUBBMA(AO,G,NH,M1,ID,MO,MJ0,MJ)
 1011 CALL SUBBMA(AO,G,NH,M1,ID,MO,MJ1,MJ0)
C     INPUT MATRIX (G) PRINT AND PUNCH OUT
cc      WRITE(6,61)
C     COMMON SUBROUTINE CALL
cc      CALL SUBMPR(G,M1,ID,MJ1,MJ)
cc      IRP1=IR+1
cc      IF(IRP1.GT.M1) GO TO 1400
cc      DO 2010 I=IRP1,M1
cc 2010 WRITE(7,2) (G(I,J),J=1,ID)
cc 1400 CALL FLCLS3(NFL)
cc  999 STOP
      CONTINUE
      RETURN
    1 FORMAT(16I5)
    2 FORMAT(4D20.10)
    6 FORMAT(//1H ,21HCANONICAL CORRELATION)
    7 FORMAT(/1H ,'NUMBER OF PRESENT AND FUTURE VARIABLES',2X,
     A'M1=',I5)
    8 FORMAT(/1H ,'NUMBER OF PRESENT AND PAST VARIABLES',4X,
     A'M2=',I5)
    9 FORMAT(/1H ,'DATA LENGTH=N=',I6)
   35 FORMAT(/1H ,'FUTURE SET CANONICAL WEIGHTS, ROWWISE')
   49 FORMAT(/1H ,5X,8HORDER(P),4X,11HCANONICAL R,3X,9HR-SQUARED,3X,
     A 10HCHI-SQUARE,3X,6HN.D.F.,2X,23HDIC (P)(=CHI**2-2*D.F.))
   50 FORMAT(/1H ,10X,I3,4X,F8.4,4X,F8.3,7X,F8.2,4X,I6,3X,F10.4)
 4410 FORMAT(/1H ,'MINIMUM DIC(P) =',F12.2,1X,'ATTAINED AT P='
     A,I5)
 5200 FORMAT(/1H ,5X,'F(I)')
 1130 FORMAT(/1H ,'STRUCTURAL CHARACTERISTIC VECTOR',
     A' (H(I),I=1,P)')
 1131 FORMAT(1H ,4X,10I12)
   61 FORMAT(/1H ,'G-MATRIX')
  259 FORMAT(1H ,6HINW(I),5X,10I5)
 8200 FORMAT(/1H ,I5,10D12.5,/(1H ,5X,10D12.5))
11111 FORMAT(/1H ,'INITIAL AUTO REGRESSIVE MODEL FITTING ',
     A'BY THE MINIMUM AIC PROCEDURE.')
11112 FORMAT(1H ,8X,'        THE VALUE OF CHI-SQUARE AND DIC (P) ',
     A'CORRESPONDING TO CANONICAL R=1.000 ',
     A'SHOULD BE IGNORED')
11114 FORMAT(1H ,'PROGRAM 74.2.1. CANOCA')
11113 FORMAT(/1H ,'IR=',I5)
      END
C
cc	SUBROUTINE AVMLVC(Y,Z,MM)
      SUBROUTINE AVMLVC(Y,Z,MM,FL,MJ1)
C     Z=FL'*Y
C     FL: LOWER TRIANGLE
C     Y: VECTOR
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION Y(MM),Z(MM)
cc	DIMENSION FL(46,46)
cc	COMMON /COM25/FL
      DIMENSION FL(MJ1,MJ1)
      CST0=0.0D-00
      DO 10 I=1,MM
      SUM=CST0
      DO 11 J=I,MM
   11 SUM=SUM+FL(J,I)*Y(J)
      Z(I)=SUM
   10 CONTINUE
      RETURN
      END
C
cc	SUBROUTINE BLMLVC(Y,Z,MM)
      SUBROUTINE BLMLVC(Y,Z,MM,FL,MJ1)
C     Z=FL*Y
C     FL: LOWER TRIANGLE
C     Y: VECTOR
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION Y(MM),Z(MM)
cc	DIMENSION FL(46,46)
cc	COMMON /COM25/FL
      DIMENSION	 FL(MJ1,MJ1)
      CST0=0.0D-00
      DO 10 I=1,MM
      SUM=CST0
      DO 11 J=1,I
   11 SUM=SUM+FL(I,J)*Y(J)
      Z(I)=SUM
   10 CONTINUE
      RETURN
      END
C
cc	SUBROUTINE MWTFL(V,VV,MM)
      SUBROUTINE MWTFL(V,VV,MM,FL,MJ1)
C     VV=V'*FL
      IMPLICIT REAL *8(A-H,O-Z)
cc	COMMON /COM20/V
cc	COMMON /COM21/VV
cc	COMMON /COM25/FL
cc	DIMENSION V(46,46),VV(46,46),FL(46,46)
      DIMENSION V(MJ1,MJ1),VV(MJ1,MJ1),FL(MJ1,MJ1)
      CST0=0.0D-00
      DO 10 I=1,MM
      DO 11 J=1,MM
      SUM=CST0
      DO 12 K=1,MM
   12 SUM=SUM+V(K,I)*FL(K,J)
      VV(I,J)=SUM
   11 CONTINUE
   10 CONTINUE
      RETURN
      END
C
cc	SUBROUTINE NLTIV(R,RIN,DET,K,MJ)
      SUBROUTINE NLTIV(R,RIN,DET,K)
C     INVERSE OF R IS FACTORI
C     INVERSE OF R IS FACTORED INTO THE FORM L'*L IS
C     INVERSE OF R IS FACTORED INTO THE FORM L'*L. L,LOWER
C     INVERSE OF R IS FACTORED INTO THE FORM L'*L. L, LOWER TRIANGULAR,
C     IS RETURNED IN R. OTHER OUTPUTS ARE
C     RIN: INVERSE OF DIAGONAL L(I,J)
C     DET: DETERMINANT OF R(I,J)
C     MJ IS THE ABSOLUTE DIMENSION OF R IN THE MAIN ROUTINE
      IMPLICIT REAL*8(A-H,O-Z)
cc	DIMENSION R(MJ,MJ),RIN(K)
      DIMENSION R(K,K),RIN(K)
      CST0=0.0D-00
      CST1=1.0D-00
      KM1=K-1
      ABC0=CST0
      DET=CST1
      DO 10 L=1,K
      RR=R(L,L)
      DET=DET*RR
      RPIVOT=CST1/DSQRT(RR)
      R(L,L)=RPIVOT
      RIN(L)=CST1/RPIVOT
      DO 12 I=1,K
      IF(I.EQ.L) GO TO 12
      R(L,I)=RPIVOT*R(L,I)
   12 CONTINUE
      IF(L.EQ.K) GO TO 11
      L1=L+1
      DO 13 I=L1,K
      RIL=-RPIVOT*R(I,L)
      R(I,L)=RIL*RPIVOT
      DO 14 M=1,K
      IF(M.EQ.L) GO TO 14
      R(I,M)=R(I,M)+RIL*R(L,M)
   14 CONTINUE
   13 CONTINUE
   10 CONTINUE
   11 RETURN
      END
C
cc	SUBROUTINE NWFPEC(OSD,AO,L,IR,IL,N,IFPEC,MJ0,MJ)
      SUBROUTINE NWFPEC(AAIC,OAIC,CV,AST1,OSD,AO,AAO,L,IR,IL,N,IFPEC,
     *                  MJ0,LMAX2,LCV1)
C     AR-FITTING
C     AUTOREGRESSIVE MODEL FITTING BY THE MINIMUM AIC PROCEDURE.
      IMPLICIT REAL*8(A-H,O-Z)
cc	DIMENSION CV(25,7,7)
cc	DIMENSION A1(13,7,7),B1(13,7,7)
cc	DIMENSION SD(7,7),SE(7,7),SF(7,7)
cc	DIMENSION XSD(7,7),XSF(7,7),D(7,7),E(7,7),Z1(7,7)
cc	DIMENSION SFL(7,7),RIN(7),B(7,7)
cc	DIMENSION AST1(91,7,7)
cc	DIMENSION OSD(7,7),AO(13,7,7)
      DIMENSION CV(LCV1,IR+IL,IR+IL)
      DIMENSION A1(L,IR+IL,IR+IL),B1(L,IR+IL,IR+IL)
      DIMENSION SD(IR+IL,IR+IL),SE(IR+IL,IR+IL),SF(IR+IL,IR+IL)
      DIMENSION XSD(IR+IL,IR+IL),XSF(IR+IL,IR+IL)
      DIMENSION D(IR+IL,IR+IL),E(IR+IL,IR+IL),Z1(IR+IL,IR+IL)
      DIMENSION SFL(IR+IL,IR+IL),RIN(IR+IL),B(IR+IL,IR+IL)
      DIMENSION AST1(LMAX2,IR+IL,IR+IL)
      DIMENSION OSD(IR,IR),AO(MJ0,IR,IR+IL),AAO(MJ0,IR,IR+IL)
      DIMENSION AAIC(0:MJ0-1)
cc	COMMON /COM9/AST1
cc	COMMON /COM10/CV
c
	IP=IR+IL
      L1=L+1
      CST0=0.0D-00
C     INITIAL CONDITION AND COVARIANCE PRINT OUT
C     INITIAL SD, SF, SE COMPUTATION
      DO 330 II=1,IP
      DO 330 JJ=1,IP
      SD(II,JJ)=CV(1,II,JJ)
      SF(II,JJ)=SD(II,JJ)
      SFL(II,JJ)=SF(II,JJ)
      SE(II,JJ)=CV(2,II,JJ)
      XSD(II,JJ)=SD(II,JJ)
  330 XSF(II,JJ)=SF(II,JJ)
C     0-TH STEP COMPUTATION
      IFPEC=0
      MS=0
C     AIC COMPUTATION
cc	CALL SAIC (SD,N,IP,MS,AIC,MJ)
      CALL SAIC (SD,N,IP,MS,AIC)
      OAIC=AIC
C     AIC PRINT OUT
cc	WRITE(6,600)
cc	WRITE(6,264) MS,AIC
      AAIC(0)=AIC
cc	CALL NLTIV(SFL,RIN,SFDT,IP,MJ)
      CALL NLTIV(SFL,RIN,SFDT,IP)
C     SFL IS FACTORED INTO L'*L AND L IS RETURNED IN SFL
      INX=1
      DO 8 II=1,IP
      DO 8 JJ=1,II
      AST1(1,JJ,II)=CST0
    8 AST1(1,II,JJ)=SFL(II,JJ)
C     ITERATION M=1 TO L
      DO 400 M=1,L
C     INVERSE OF SD, SF COMPUTATION
C     COMMON SUBROUTINE CALL
cc	CALL INVDET(XSD,SDDET,IP,MJ)
      CALL INVDET(XSD,SDDET,IP,IP)
C     COMMON SUBROUTINE CALL
cc	CALL INVDET(XSF,SFDET,IP,MJ)
      CALL INVDET(XSF,SFDET,IP,IP)
C     D, E, SD, SF COMPUTATION
C     COMMON SUBROUTINE CALL
cc	CALL MULPLY(SE,XSF,D,IP,IP,IP,MJ,MJ,MJ)
      CALL MULPLY(SE,XSF,D,IP,IP,IP)
C     COMMON SUBROUTINE CALL
cc	CALL TRAMDL(SE,XSD,E,IP,IP,IP,MJ,MJ,MJ)
      CALL TRAMDL(SE,XSD,E,IP,IP,IP)
C     COMMON SUBROUTINE CALL
cc	CALL TRAMDR(D,SE,Z1,IP,IP,IP,MJ,MJ,MJ)
      CALL TRAMDR(D,SE,Z1,IP,IP,IP)
C     COMMON SUBROUTINE CALL
cc	CALL SUBTAL(SD,Z1,IP,IP,MJ,MJ)
      CALL SUBTAL(SD,Z1,IP,IP)
C     COMMON SUBROUTINE CALL
cc	CALL MULPLY(E,SE,Z1,IP,IP,IP,MJ,MJ,MJ)
      CALL MULPLY(E,SE,Z1,IP,IP,IP)
C     COMMON SUBROUTINE CALL
cc	CALL SUBTAL(SF,Z1,IP,IP,MJ,MJ)
      CALL SUBTAL(SF,Z1,IP,IP)
      MS=M
      DO 410 II=1,IP
      DO 410 JJ=1,IP
      XSD(II,JJ)=SD(II,JJ)
      XSF(II,JJ)=SF(II,JJ)
      SFL(II,JJ)=SF(II,JJ)
  410 CONTINUE
C     AIC COMPUTATION
cc	CALL SAIC (SD,N,IP,MS,AIC,MJ)
      CALL SAIC (SD,N,IP,MS,AIC)
C     AIC PRINT OUT
cc	WRITE(6,264) MS,AIC
      AAIC(MS) = AIC
cc	CALL NLTIV(SFL,RIN,SFDT,IP,MJ)
      CALL NLTIV(SFL,RIN,SFDT,IP)
C     SFL IS FACTORED INTO L'*L AND L IS RETURNED IN SFL
C     FORWARD AND BACKWARD PREDICTOR COMPUTATION
cc	CALL COEFAB(A1,B1,D,E,MS,IP,MJ0,MJ)
      CALL COEFAB(A1,B1,D,E,MS,L,IP)
      DO 14 I=1,M
      IMI=M+1-I
      DO 15 II=1,IP
      DO 15 JJ=1,IP
   15 B(II,JJ)=B1(IMI,II,JJ)
cc	CALL BLMULP(SFL,B,Z1,IP,IP,MJ,MJ)
      CALL BLMULP(SFL,B,Z1,IP,IP)
      INX=INX+1
      DO 16 II=1,IP
      DO 16 JJ=1,IP
   16 AST1(INX,II,JJ)=-Z1(II,JJ)
   14 CONTINUE
      INX=INX+1
      DO 17 II=1,IP
      DO 17 JJ=1,II
      AST1(INX,JJ,II)=CST0
   17 AST1(INX,II,JJ)=SFL(II,JJ)
C     MINIMUM AIC SEARCH
      IF(OAIC.LE.AIC) GO TO 440
      OAIC=AIC
      IFPEC=M
      DO 560 II=1,IR
      DO 560 JJ=1,IR
  560 OSD(II,JJ)=SD(II,JJ)
      DO 561 I=1,M
      DO 562 II=1,IR
      DO 562 JJ=1,IP
      AAO(I,II,JJ)=-A1(I,II,JJ)
  562 AO(I,II,JJ)=-A1(I,II,JJ)
  561 CONTINUE
  440 IF(M.EQ.L) GO TO 400
C     SE COMPUTATION
C     COMMON SUBROUTINE CALL
cc	CALL NEWSE(A1,SE,MS,IP,MJ0,MJ)
      CALL NEWSE(A1,CV,SE,MS,L,IP,LCV1)
  400 CONTINUE
C     MIN.AIC PRINT OUT
cc	WRITE(6,607) OAIC,IFPEC
C     OSD, AO PRINT OUT
cc	WRITE(6,608)
C     COMMON SUBROUTINE CALL
cc	CALL SUBMPR(OSD,IR,IR,MJ,MJ)
      IF(IFPEC.LE.0) GO TO 699
cc	WRITE(6,1601)
cc	WRITE(6,1602)
cc	WRITE(6,609)
C     COMMON SUBROUTINE CALL
cc	CALL PRMAT3(AO,IFPEC,IR,IP,0,MJ0,MJ,MJ)
      DO 1611 I=1,M
      DO 1612 II=1,IR
      DO 1612 JJ=1,IP
 1612 AO(I,II,JJ)=-AO(I,II,JJ)
 1611 CONTINUE
  699 RETURN
    1 FORMAT(10I5)
    2 FORMAT(4D20.10)
   42 FORMAT(//1H ,17HCOVARIANCE MATRIX)
  264 FORMAT(1H ,I5,2X,D16.7)
  600 FORMAT(///1H ,4X,1HI,15X,3HAIC)
  607 FORMAT(/1H ,'MINIMUM AIC ',D12.5,2X,'ATTAINED AT M=',I5)
  608 FORMAT(//1H ,10X,10HOSD(II,JJ),': INNOVATION VARIANCE')
  609 FORMAT(//1H ,10X,'A(I): AUTOREGRESSIVE COEFFICIENTS')
 1601 FORMAT(/1H ,'AR-MODEL:')
 1602 FORMAT(1H ,'Y(N)+A(1)Y(N-1)+...+A(L)Y(N-L)=X(N)')
      END
C
cc	SUBROUTINE SAIC (SD,N,K,MS,AIC,MJ)
      SUBROUTINE SAIC (SD,N,K,MS,AIC)
C     AIC COMPUTATION.
C     SD: COVARIANCE MATRIX OF INNOVATION
      IMPLICIT REAL*8(A-H,O-Z)
cc	DIMENSION SD(MJ,MJ)
cc	DIMENSION SD1(7,7)
      DIMENSION SD(K,K)
      DIMENSION SD1(K,K)
      AN=N
      DO 9 I=1,K
      DO 9 J=1,K
    9 SD1(I,J)=SD(I,J)
C     COMMON SUBROUTINE CALL
cc	CALL SUBDET(SD1,SDRM,K,MJ)
      CALL SUBDETC(SD1,SDRM,K)
      ARM2=2*MS*K*K
      AIC=AN*DLOG(SDRM)+ARM2
      RETURN
      END
C
cc	SUBROUTINE SBRUGT(MO,ID)
      SUBROUTINE SBRUGT(MO,ID,AST1,CV,RGT,MJ1,MJ,LMAX2,LCV1)
C     THIS SUBROUTINE COMPUTES MATRIX R12*G'.
C     INPUTS REQUIRED:
C     ID: DIMENSION OF Y(I) =IR=IP
C     MO: MAICE ORDER OF AR-MODEL
C     CV(I,II,JJ): COVARIANCE MATRIX
C     AST1(I,II,JJ): MATRIX G FOR ORTHO-NORMALIZATION OF THE PRESENT AND
C     PAST VARIABLES. G IS A LOWER TRIANGULAR BLOCK MATRIX WITH THE
C     STRUCTURE
C		 +-			   -+
C	     G = I  L0	   0	  0   ..... I
C		 I -L1B11  L1	  0   ..... I
C		 I -L2B22 -L2B21  L2   .... I
C		 I   .	    .	   ..	... I
C		 I   .	    .	   . .	 .. I
C		 I   .	    .	   .  .	  . I
C		 I   .	    .	   .   .    I
C		 I   .	    .	   .	.   I
C		 I   .	    .	   .	 .  I
C		 I   .	    .	   .	  . I
C		 +-			   -+
C     WHERE BMI'S ARE THE COEFFICIENTS OF THE M-TH ORDER BACKWARD AUTO-
C     REGRESSION,
C     Y(I-M)-BM1Y(I-M+1)-..-BMMY(I)=ZM(I-M)
C     AND LM(LM)'=INVERSE OF SD(M),WHERE SD(M) IS THE COVARIANCE MATRIX
C     RESIDUAL ZM(I).
C     OUTPUT:
C     RGT(I,J): MATRIX OF R12*G'
      IMPLICIT REAL*8(A-H,O-Z)
cc	DIMENSION CV(25,7,7),AST1(91,7,7)
cc	DIMENSION RGT(91,91)
cc	DIMENSION X(7,7),Y(7,7)
      DIMENSION CV(LCV1,ID,ID),AST1(LMAX2,MJ,MJ)
      DIMENSION RGT(MJ1,MJ1)
      DIMENSION X(ID,ID),Y(ID,ID)
cc	COMMON /COM9/AST1
cc	COMMON /COM10/CV
cc	COMMON /COM11/RGT
      CST0=0.0D-00
      MP1=MO+1
cc	IRG=91
      IRG=MJ1
      DO 9 I=1,IRG
      DO 9 J=1,IRG
      RGT(I,J)=CST0
    9 CONTINUE
      INC=0
      DO 10 I=1,MP1
      IM1=I-1
      INX=0
      JNC=0
      DO 11 J=1,MP1
      DO 12 K=1,J
      IB=IM1+K
      JB=INX+K
      DO 18 II=1,ID
      DO 18 JJ=1,ID
      X(II,JJ)=CV(IB,II,JJ)
      Y(II,JJ)=AST1(JB,II,JJ)
   18 CONTINUE
      MJ7=7
      DO 15 II=1,ID
      DO 16 JJ=1,ID
      SUM=CST0
      DO 17 KK=1,ID
   17 SUM=SUM+X(II,KK)*Y(JJ,KK)
      IN=INC+II
      JN=JNC+JJ
      RGT(IN,JN)=RGT(IN,JN)+SUM
   16 CONTINUE
   15 CONTINUE
   12 CONTINUE
      JNC=JNC+ID
      INX=INX+J
   11 CONTINUE
      INC=INC+ID
   10 CONTINUE
      RETURN
      END
C
cc	SUBROUTINE SUBBMA(AO,B,NH,M1,ID,IQ,MJ0,MJ)
      SUBROUTINE SUBBMA(AO,B,NH,M1,ID,IQ,MJ1,MJ0)
C     AR-FITTING
C     B-MATRIX COMPUTATION
C     (USE AR-COEFFICIENTS)
C     M1 IS LESS THAN 47.
      IMPLICIT REAL*8(A-H,O-Z)
cc      DIMENSION AO(MJ0,MJ,MJ)
cc	DIMENSION X(13,7,7)
cc	DIMENSION XX(7,7),AA(7,7),C(7,7)
cc	DIMENSION W(100,7),B(46,7)
cc	DIMENSION NH(91)
      DIMENSION AO(MJ0,ID,ID)
      DIMENSION X(IQ,ID,ID)
      DIMENSION XX(ID,ID),AA(ID,ID),C(ID,ID)
      DIMENSION W(100,ID),B(MJ1,ID)
      DIMENSION NH(M1)
      CST0=0.0D-00
      CST1=1.0D-00
      DO 10 I=1,ID
      DO 11 J=1,ID
   11 W(I,J)=CST0
   10 W(I,I)=CST1
      IF(IQ.LE.1) GO TO 110
      MQ1=IQ-1
      DO 200 II=1,MQ1
      DO 210 I=1,ID
      DO 210 J=1,ID
  210 X(II,I,J)=AO(II,I,J)
      IF(II.LE.1) GO TO 201
      IIM1=II-1
      DO 220 JJ=1,IIM1
      DO 230 I=1,ID
      DO 230 J=1,ID
  230 AA(I,J)=AO(JJ,I,J)
      IJ=II-JJ
      DO 240 I=1,ID
      DO 240 J=1,ID
  240 XX(I,J)=X(IJ,I,J)
cc	CALL MULPLY(AA,XX,C,ID,ID,ID,MJ,MJ,MJ)
      CALL MULPLY(AA,XX,C,ID,ID,ID)
      DO 250 I=1,ID
      DO 250 J=1,ID
  250 X(II,I,J)=X(II,I,J)+C(I,J)
  220 CONTINUE
  201 DO 260 I=1,ID
      DO 260 J=1,ID
      IC=ID+(II-1)*ID+I
  260 W(IC,J)=X(II,I,J)
  200 CONTINUE
  110 DO 300 I=1,M1
      NHC=NH(I)
      DO 310 J=1,ID
  310 B(I,J)=W(NHC,J)
  300 CONTINUE
  100 RETURN
      END
C
cc	SUBROUTINE VECMTX(X,Z,NH,MM,NN)
      SUBROUTINE VECMTX(X,Z,NH,RGT,MM,NN,MJ1)
C     Z=X*Y (X,Z: VECTORS, Y: SUBMATRIX)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(MM),Z(NN)
cc	DIMENSION RGT(91,91)
cc	DIMENSION NH(91)
      DIMENSION RGT(MJ1,MJ1)
      DIMENSION NH(MM)
cc	COMMON /COM11/RGT
      CST0=0.0D-00
      DO 10 I=1,NN
      SUM=CST0
      DO 11 J=1,MM
      JJ=NH(J)
   11 SUM=SUM+X(J)*RGT(JJ,I)
   10 Z(I)=SUM
      RETURN
      END
C
cc	SUBROUTINE BLMULP(X,Y,Z,MM,NN,MJ1,MJ2)
      SUBROUTINE BLMULP(X,Y,Z,MM,NN)
C     COMMON SUBROUTINE
C     Z=X*Y (X: LOWER TRIANGLE)
      IMPLICIT REAL*8(A-H,O-Z)
cc	DIMENSION X(MJ1,MJ1),Y(MJ1,MJ2),Z(MJ1,MJ2)
      DIMENSION X(MM,MM),Y(MM,NN),Z(MM,NN)
      CST0=0.0D-00
      DO 10 I=1,MM
      DO 11 J=1,NN
      SUM=CST0
      DO 12 K=1,I
   12 SUM=SUM+X(I,K)*Y(K,J)
      Z(I,J)=SUM
   11 CONTINUE
   10 CONTINUE
      RETURN
      END
C
C
cc	SUBROUTINE SUBDET(X,XDETMI,MM,MJ)
      SUBROUTINE SUBDETC(X,XDETMI,MM)
C     COMMON SUBROUTINE
C     THIS SUBROUTINE COMPUTES THE DETERMINANT OF UPPER LEFT MM X MM
C     OF X.  FOR GENERAL USE STATEMENTS 20-21 SHOULD BE RESTORED.
C     X: ORIGINAL MATRIX
C     XDETMI: DETERMINANT OF UPPER LEFT MM X MM OF X
C     MJ: ABSOLUTE DIMENSION OF X IN THE MAIN ROUTINE
      IMPLICIT REAL*8(X)
cc	DIMENSION X(MJ,MJ)
      DIMENSION X(MM,MM)
      CST0=0.0D-00
      CST1=1.0D-00
      XDETMI=CST1
      IF(MM.EQ.1) GO TO 18
      MM1=MM-1
      DO 10 I=1,MM1
C   20 IF(X(I,I).NE.CST0) GO TO 11
C     DO 12 J=I,MM
C     IF(X(I,J).EQ.CST0) GO TO 12
C     JJ=J
C     GO TO 13
C  12 CONTINUE
C     XDETMI=CST0
C     GO TO 17
C  13 DO 14 K=I,MM
C     XXC=X(K,JJ)
C     X(K,JJ)=X(K,I)
C  14 X(K,I)=XXC
C  21 XDETMI=-XDETMI
   11 XDETMI=XDETMI*X(I,I)
      XC=CST1/X(I,I)
      I1=I+1
      DO 15 J=I1,MM
      XXC=X(J,I)*XC
      DO 16 K=I1,MM
   16 X(J,K)=X(J,K)-X(I,K)*XXC
   15 CONTINUE
   10 CONTINUE
   18 XDETMI=XDETMI*X(MM,MM)
   17 RETURN
      END
