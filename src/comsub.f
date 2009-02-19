      INCLUDE 'timsac_f.h'
C------------------------------------------------- 
C     AUSP	---  (72) auspec, mulspe
C     COEFAB  ---  (72) fpec7 (74) canoca
C     CORNOM  ---  (72) autcor, mulcor, fftcor
C     CROSCO  ---  (72) autcor, mulcor
C     DMEADL  ---  (72) autcor, mulcor, fftcor
C     ECORCO  ---  (72) mulspe
C     FGER1	---  (72) mulnos, mulrsp
C     FGERCO  ---  (72) auspec, mulspe
C     INVDET  ---  (72) fpec7, optdes (74) canoca (78) mulmar, perars, mlomar, exsar
C     INVDETC ---	 (72) mulnos, mulrsp
C     LTINV   ---  (72) wnoise (74) simcon
C     LTRVEC  ---  (72) wnoise (74) simcon
C     MATADL  ---  (72) fpec7, optdes
C     MIXRAD  ---  (72) fftcor (74) covgen
C     MULPLY  ---  (72) fpec7, optdes
C     MULVER  ---  (72) optsim (74) simcon
C     NEWSE   ---  (72) fpec7, (74) canoca
C     REARRA  ---  (72) fpec7
C     REARRAC ---  (72) mulfre
C     SIGNIF  ---  (72) auspec, mulspe
C     SMOSPE  ---  (72) mulspe
C     SUBD12  ---  (72) sglfre, mulfre
C     SUBDET  ---  (72) fpec7, (78) exsar
C     SUBNOS  ---  (72) mulnos
C     SUBTAC  ---  (72) fpec7, coefab (74) canoca
C     TRAMDL  ---  (72) fpec7, optdes
C     TRAMDR  ---  (72) fpec7, (74) canoca
C
C     DSUMF   ---  (72) autcor, mulcor, fftcor
C     RANDM   ---  (72) wnoise
C
C------------------------------------------------- 
C     DINT    ---  (74) prdctr, simcon
C     INNERP  ---  (74) autarm, canoca, markov
C     MATINV  ---  (74) autarm, markov
C     MSVD    ---  (74) canarm, canoca
C     SUBTAC  ---  (74) canoca, simcon
C
C------------------------------------------------- 
C     ADDVAR  ---  (78) mulmar, perars, mlomar
C     AICCOM  ---  (78) mulmar, perars, mlomar
C     ARBAYS  ---  (78) unibar, blocar
C     ARCOEF  ---  (78) unibar, blocar, exsar
C     ARMFIT  ---  (78) unimar, bsubst, mlocar, exsar
C     BAYSPC  ---  (78) unibar, blocar
C     BAYSWT  ---  (78) unibar, mulbar, blocar, blomar
C     COEF2   ---  (78) mulmar, perars, mlomar
C     COMAIC  ---  (78) unimar, unibar, bsubst, mlocar, blocar, exsar
C     COPY    ---  (78) mulmar, mulbar, perars, mlocar, mlomar, blomar
C     DELETE  ---  (78) mulmar, perars, mlomar
C     FOUGER  ---  (78) unibar, mlocar, blocar, (72) raspec (74) nonst
C     HUSHL1  ---  (78) mulmar, mulmar, perars, mlomar, blomar
C     HUSHLD  ---  (78) unimar, unibar, bsubst, mulabr, perars, mlocar, blocar, mlomar, blomar, easer
C     MAICE   ---  (78) mulmar, perars, mlocar, blocar, mlomar, blomar, exsar
C     MARCOF  ---  (78) mulbar, blomar
C     MARFIT  ---  (78) mulmar, perars, mlomar
C     MBYSAR  ---  (78) mulbar, blomar
C     MBYSPC  ---  (78) mulbar, blomar
C     MCOEF   ---  (78) mulmar, perars, mlomar
C     MPARCO  ---  (78) mulbar, blomar
C     MRDATA  ---  (78) mulmar, mulbar, mlomar, blomar
C     MREDCT  ---  (78) mulmar, mulbar, perars, mlomar, blomar
C     MSDCOM  ---  (78) mulbar, blomar
C     MSETX1  ---  (78) mulmar, mulbar, perars, mlomar, blomar
C     NRASPE  ---  (78) unibar, bsubst, mlocar, blocar
C     PARCOR  ---  (78) exsar (84) decomp
C     RECOEF  ---  (78) unimar, bsubst, mlocar, exsar (84) decomp
C     REDATA  ---  (78) unimar, unibar, bsubst, perars, mlocar, blocar, exsar
C     REDUCT  ---  (78) unimar, unibar, bsubst, mlocar, blocar, exsar
C     SDCOMP  ---  (78) unibar, bsubst, blocar
C     SETX1   ---  (78) unimar, unibar, bsubst, mlocar, blocar, exsar
C     SOLVE   ---  (78) mulbar, blomar
C     SRCOEF  ---  (78) mulmar, perars, mlomar
C     TRIINV  ---  (78) mulmar, perars
C
C     DMIN*8  ---  (78) blocar, blomar
C
C-------------------------------------------------
C     ARMASP  ---  (Iwanami) arma, tvspc
C     FOURIE  ---  (Iwanami) arma, tvspc
C     MOMENTK ---  (Iwanami) ngsmth, tvvar
C     SMOOTH  ---  (Iwanami) tvvar, smooth
C     GINVRS  ---  (Iwanami) tvvar, smooth
C
C-------------------------------------------------
C
      SUBROUTINE AUSP(FC,P1,LAGH1,A,LA1)
C     THIS SUBROUTINE COMPUTES SMOOTHED AUTO SPECTRUM.
C     FC: OUTPUT OF FGERCO
C     P1: SMOOTHED SPECTRUM
C     LAGH1: DIMENSION OF FC AND P1
C     A: SMOOTHING COEFFICIENTS
C     LA1: DIMENSION OF A (LESS THAN 11)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION FC(LAGH1),P1(LAGH1),A(LA1)
c      DIMENSION FC1(521)
      DIMENSION FC1(LAGH1+2*(LA1-1))
      LA=LA1-1
      LAGSHF=LAGH1+2*LA
C     FC SHIFT-RIGHT BY LA FOR END CORRECTION
      CALL ECORCO(FC,LAGH1,FC1,LAGSHF,LA1)
C     SMOOTHING
      CALL SMOSPE(FC1,LAGSHF,A,LA1,P1,LAGH1)
      RETURN
      END
C
C
cc	SUBROUTINE COEFAB(A1,B1,D,E,MS,K,MJ0,MJ)
      SUBROUTINE COEFAB(A1,B1,D,E,MS,L,K)
C     AR-FITTING
C     THIS SUBROUTINE COMPUTES FORWARD(A) AND BACKWARD(B) PREDICTOR
C     COEFFICIENTS.
      IMPLICIT REAL*8(A-H,O-Z)
cc	DIMENSION A1(MJ0,MJ,MJ),B1(MJ0,MJ,MJ)
cc	DIMENSION D(MJ,MJ),E(MJ,MJ)
cc	DIMENSION A(7,7),B(7,7),Z1(7,7),Z2(7,7)
      DIMENSION A1(L,K,K),B1(L,K,K)
      DIMENSION D(K,K),E(K,K)
      DIMENSION A(K,K),B(K,K),Z1(K,K),Z2(K,K)
      IF(MS.EQ.1) GO TO 40
      MSM1=MS-1
      DO 10 I=1,MSM1
      MMI=MS-I
      DO 20 II=1,K
      DO 20 JJ=1,K
      A(II,JJ)=A1(I,II,JJ)
   20 B(II,JJ)=B1(MMI,II,JJ)
cc	CALL MULPLY(D,B,Z1,K,K,K,MJ,MJ,MJ)
cc	CALL MULPLY(E,A,Z2,K,K,K,MJ,MJ,MJ)
cc	CALL SUBTAL(A,Z1,K,K,MJ,MJ)
cc	CALL SUBTAL(B,Z2,K,K,MJ,MJ)
      CALL MULPLY(D,B,Z1,K,K,K)
      CALL MULPLY(E,A,Z2,K,K,K)
      CALL SUBTAL(A,Z1,K,K)
      CALL SUBTAL(B,Z2,K,K)
      DO 21 II=1,K
      DO 21 JJ=1,K
      A1(I,II,JJ)=A(II,JJ)
   21 B1(MMI,II,JJ)=B(II,JJ)
   10 CONTINUE
   40 DO 30 II=1,K
      DO 30 JJ=1,K
      A1(MS,II,JJ)=D(II,JJ)
   30 B1(MS,II,JJ)=E(II,JJ)
      RETURN
      END
C
      SUBROUTINE CORNOM(C,CN,LAGH1,CX0,CY0)
C     NORMALIZATION OF COVARIANCE
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION C(LAGH1),CN(LAGH1)
      CST1=1.0D-00
      DS=CST1/DSQRT(CX0*CY0)
      DO 10 I=1,LAGH1
   10 CN(I)=C(I)*DS
      RETURN
      END
C
      SUBROUTINE CROSCO(X,Y,N,C,LAGH1)
C     THIS SUBROUTINE COMPUTES C(L)=COVARIANCE(X(S+L),Y(S))
C     (L=0,1,...,LAGH1-1).
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(N),Y(N),C(LAGH1)
      AN=N
      BN1=1.0D-00
      BN=BN1/AN
      CT0=0.0D-00
      DO 10 II=1,LAGH1
      I=II-1
      T=CT0
      IL=N-I
      DO 20 J=1,IL
      J1=J+I
   20 T=T+X(J1)*Y(J)
   10 C(II)=T*BN
      RETURN
      END
C
      SUBROUTINE DMEADL(X,N,XMEAN)
C     DOUBLE PRECISION MEAN DELETION
cc      DOUBLE PRECISION X,XMEAN,AN
      DOUBLE PRECISION X,XMEAN,AN,DSUMF
      DIMENSION X(N)
      AN=N
      XMEAN=DSUMF(X,N)/AN
      DO 10 I=1,N
   10 X(I)=X(I)-XMEAN
      RETURN
      END
C
      SUBROUTINE ECORCO(FC,LAGH1,FC1,LAGSHF,LA1)
C     FC SHIFT-RIGHT BY LA FOR REAL PART END CORRECTION
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION FC(LAGH1),FC1(LAGSHF)
      LAGH2=LAGH1+1
      LA=LA1-1
      DO 100 I=1,LAGH1
      I1=LAGH2-I
      I2=I1+LA
  100 FC1(I2)=FC(I1)
      LA2=LAGH1+LA
      DO 110 I=1,LA
      I1=LA1-I
      I2=LA1+I
      I3=LA2-I
      I4=LA2+I
      FC1(I1)=FC1(I2)
  110 FC1(I4)=FC1(I3)
      RETURN
      END
C
c      SUBROUTINE FGER1
      SUBROUTINE FGER1(G,GR,GI,LG,H,JJF)
C     FOURIER TRANSFORM(GOERTZEL METHOD)
C     THIS SUBROUTINE COMPUTES ONE VALUE OF THE FOURIER TRANSFORM BY
C     GOERTZEL METHOD.
      IMPLICIT REAL*8(A-H,O-W)
      INTEGER H
c      COMMON G,GR,GI,LG,H,JJF
c      DIMENSION G(31)
      DIMENSION G(LG+1)
      CST0=0.0D-00
      LGP1=LG+1
C     REVERSAL OF G(I),I=1,...,LGP1 INTO G(LG3-I)   LG3=LGP1+1
      IF(LGP1.LE.1) GO TO 110
      LG3=LGP1+1
      LG4=LG3/2
      DO 100 I=1,LG4
      I2=LG3-I
      T=G(I)
      G(I)=G(I2)
  100 G(I2)=T
  110 PI=3.1415926536
      AH=H
      T=PI/AH
      AK=JJF-1
      TK=T*AK
      CK=DCOS(TK)
      SK=DSIN(TK)
      CK2=CK+CK
      UM2=CST0
      UM1=CST0
      IF(LG.EQ.0) GO TO 12
      DO 11 I=1,LG
      UM0=CK2*UM1-UM2+G(I)
      UM2=UM1
   11 UM1=UM0
   12 GR=CK*UM1-UM2+G(LGP1)
      GI=-SK*UM1
      RETURN
      END
C
C
      SUBROUTINE FGERCO(G,LGP1,FC,LF1)
C     FOURIER TRANSFORM (GOERTZEL METHOD)
C     THIS SUBROUTINE COMPUTES FOURIER TRANSFORM OF G(I),I=0,1,...,LG AT
C     FREQUENCIES K/(2*LF),K=0,1,...,LF AND RETURNS COSIN TRANSFORM IN
C     FC(K).
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION G(LGP1),FC(LF1)
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
  100 G(I2)=T
  110 PI=3.1415926536
      ALF=LF
      T=PI/ALF
      DO 10 K=1,LF1
      AK=K-1
      TK=T*AK
      CK=DCOS(TK)
      CK2=CK+CK
      UM2=0.0D-00
      UM1=0.0D-00
      IF(LG.EQ.0) GO TO 12
      DO 11 I=1,LG
      UM0=CK2*UM1-UM2+G(I)
      UM2=UM1
   11 UM1=UM0
   12 FC(K)=CK*UM1-UM2+G(LGP1)
   10 CONTINUE
      RETURN
      END
C
      SUBROUTINE INVDET(X,XDET,MM,MJ)                                   
C                                                                       
C       THE INVERSE AND DETERMINANT OF X COMPUTATION                    
C                                                                       
C       INPUTS:                                                         
C          X:     MM*MM SQUARE MATRIX                                   
C          MM:    DIMENSION OF X                                        
C          MJ:    ABSOLUTE DIMENSION OF X IN THE MAIN PROGRAM           
C                                                                       
C       OUTPUTS:                                                        
C          X:     INVERSE OF X                                          
C          XDET:  DETERMINANT OF X                                      
C                                                                       
      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
      DIMENSION X(MJ,MJ)                                                
cc      DIMENSION  IDS(100)                                               
      DIMENSION  IDS(MM)
      XDET = 1.0D00                                                     
      DO 10 L=1,MM                                                      
C     PIVOTING AT L-TH STAGE                                            
      XMAXP=0.10000D-10                                                 
      MAXI=0                                                            
      DO 110 I=L,MM                                                     
    1 IF( DABS(XMAXP) .GE. DABS(X(I,L)) )     GO TO 110                 
      XMAXP=X(I,L)                                                      
      MAXI=I                                                            
  110 CONTINUE                                                          
      IDS(L)=MAXI                                                       
      IF(MAXI.EQ.L) GO TO 120                                           
      IF(MAXI.GT.0) GO TO 121                                           
      XDET = 0.0D00                                                     
      GO TO 140                                                         
C     ROW INTERCHANGE                                                   
  121 DO 14 J=1,MM                                                      
      XC=X(MAXI,J)                                                      
      X(MAXI,J)=X(L,J)                                                  
   14 X(L,J)=XC                                                         
      XDET=-XDET                                                        
  120 XDET=XDET*XMAXP                                                   
      XC = 1.0D00 / XMAXP                                               
      X(L,L)=1.0D00                                                     
      DO 11 J=1,MM                                                      
   11 X(L,J)=X(L,J)*XC                                                  
      DO 12 I=1,MM                                                      
      IF(I.EQ.L) GO TO 12                                               
      XC=X(I,L)                                                         
      X(I,L) = 0.0D00                                                   
      DO 13 J=1,MM                                                      
   13 X(I,J)=X(I,J)-XC*X(L,J)                                           
   12 CONTINUE                                                          
   10 CONTINUE                                                          
      IF(MM.GT.1) GO TO 123                                             
      GO TO 140                                                         
C     COLUMN INTERCHANGE                                                
  123 MM1=MM-1                                                          
      DO 130 J=1,MM1                                                    
      MMJ=MM-J                                                          
      JJ=IDS(MMJ)                                                       
      IF(JJ.EQ.MMJ) GO TO 130                                           
      DO 131 I=1,MM                                                     
      XC=X(I,JJ)                                                        
      X(I,JJ)=X(I,MMJ)                                                  
  131 X(I,MMJ)=XC                                                       
  130 CONTINUE                                                          
  140 RETURN                                                            
      END                                                               
C
c      SUBROUTINE INVDET(X,XDET,MM,MJ)
      SUBROUTINE INVDETC(X,XDET,MM)
C     THIS SUBROUTINE COMPUTES THE INVERSE AND DETERMINANT OF
C     UPPER LEFT MM X MM OF COMPLEX MATRIX X.
C     X: ORIGINAL MATRIX
C     MM: DIMENSION OF UPPER LEFT OF X (SHOULD BE LESS THAN 11)
C     XDET: DETERMINANT OF UPPER LEFT MM X MM OF X
C     MJ: ABSOLUTE DIMENSION OF X IN THE MAIN ROUTINE
C     THE INVERSE MATRIX IS OVERWRITTEN ON THE ORIGINAL.
c      DIMENSION X(MJ,MJ)
c      DIMENSION IDS(10)
      IMPLICIT COMPLEX*16(X)
c      DIMENSION X(MJ,MJ)
c      DIMENSION IDS(10)
      DIMENSION X(MM,MM)
      DIMENSION IDS(MM)
	DOUBLE PRECISION CST0, CST1
      CST0=0.0D-00
      CST1=1.0D-00
      XDET=CST1
      DO 10 L=1,MM
C     PIVOTING AT L-TH STAGE
      XMAXP=0.10000D-10
      MAXI=0
      DO 110 I=L,MM
      IF(CDABS(XMAXP).GE.CDABS(X(I,L))) GO TO 110
      XMAXP=X(I,L)
      MAXI=I
  110 CONTINUE
      IDS(L)=MAXI
      IF(MAXI.EQ.L) GO TO 120
      IF(MAXI.GT.0) GO TO 121
      XDET=CST0
      GO TO 140
C     ROW INTERCHANGE
  121 DO 14 J=1,MM
      XC=X(MAXI,J)
      X(MAXI,J)=X(L,J)
   14 X(L,J)=XC
      XDET=-XDET
  120 XDET=XDET*XMAXP
      XC=CST1/XMAXP
      X(L,L)=CST1
      DO 11 J=1,MM
   11 X(L,J)=X(L,J)*XC
      DO 12 I=1,MM
      IF(I.EQ.L) GO TO 12
      XC=X(I,L)
      X(I,L)=CST0
      DO 13 J=1,MM
   13 X(I,J)=X(I,J)-XC*X(L,J)
   12 CONTINUE
   10 CONTINUE
      IF(MM.GT.1) GO TO 123
      GO TO 140
C     COLUMN INTERCHANGE
  123 MM1=MM-1
      DO 130 J=1,MM1
      MMJ=MM-J
      JJ=IDS(MMJ)
      IF(JJ.EQ.MMJ) GO TO 130
      DO 131 I=1,MM
      XC=X(I,JJ)
      X(I,JJ)=X(I,MMJ)
  131 X(I,MMJ)=XC
  130 CONTINUE
  140 RETURN
      END
C
cc      SUBROUTINE LTINV(R,K,MJ)                                          
      SUBROUTINE LTINV(R,K)
C     COMMON SUBROUTINE                                                 
C     THIS SUBROUTINE FACTORIZES (R(I,J): I,J=1,K) INTO R=L*L',         
C     WITH L LOWER TRIANGLE, AND GIVES L' ON AND ABOVE THE DIAGONAL OF R
C     MJ: ABSOLUTE DIMENSION OF R IN THE MAIN ROUTINE                   
      IMPLICIT REAL*8(A-H,O-Z)                                          
cc      DIMENSION R(MJ,MJ)                                                
      DIMENSION R(K,K)
      CST1=1.0D-00                                                      
      DO 10 L=1,K                                                       
      RPIVOT=CST1/DSQRT(R(L,L))                                         
      R(L,L)=CST1/RPIVOT                                                
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
cc      SUBROUTINE LTRVEC(X,Y,Z,MM,NN,MJ1,MJ2)                            
      SUBROUTINE LTRVEC(X,Y,Z,MM,NN)
C     Z=X*Y                                                             
C     (VECTOR Z)=(LOWER TRIANGLE OF UPPER LEFT MM X NN OF X)*(VECTOR Y) 
C     (MJ1,MJ2): ABSOLUTE DIMENSION OF X IN THE MAIN ROUTINE            
      IMPLICIT REAL*8(A-H,O-Z)                                          
cc      DIMENSION X(MJ1,MJ2),Y(NN),Z(MM)                                  
      DIMENSION X(MM,NN),Y(NN),Z(MM)                                  
      CST0=0.0D-00                                                      
      DO 10 I=1,MM                                                      
      SUM=CST0                                                          
      DO 11 J=1,I                                                       
   11 SUM=SUM+X(I,J)*Y(J)                                               
   10 Z(I)=SUM                                                          
      RETURN                                                            
      END                                                               
C
c      SUBROUTINE MATADL(X,Y,MM,NN,MJ1,MJ2)
      SUBROUTINE MATADL(X,Y,MM,NN)
C     MATRIX ADDITION
C     X=X+Y
C     (UPPER LEFT MM X NN OF X)=(UPPER LEFT MM X NN OF X)+(UPPER LEFT
C     MM X NN OF Y).
C     (MJ1,MJ2): ABSOLUTE DIMENSION OF X AND Y IN THE MAIN ROUTINE
      IMPLICIT REAL*8(A-H,O-Z)
c      DIMENSION X(MJ1,MJ2),Y(MJ1,MJ2)
      DIMENSION X(MM,NN),Y(MM,NN)
      DO 10 I=1,MM
      DO 10 J=1,NN
   10 X(I,J)=X(I,J)+Y(I,J)
      RETURN
      END
C
      SUBROUTINE MIXRAD(Z,N,N2P,ISG)
C     COMMON SUBROUTINE
C     MIXED RADIX FAST FOURIER TRANSFORM
C     ISG=-1...FOURIER TRANSFORM
C     ISG=1...INVERSE FOURIER TRANSFORM
      IMPLICIT REAL*8(A-H,O-Y)
      IMPLICIT COMPLEX*16(Z)
      DIMENSION Z(N)
cc	DIMENSION MS(11)
      DIMENSION MS(N2P)
      CST0=0.0D-00
      CST1=1.0D-00
      AN=N
      PI=3.1415926536
      PI2=PI+PI
      SG=ISG
      ZCI=SG*DCMPLX(CST0,CST1)
      DO 10 I=1,N2P
   10 MS(I)=2**(N2P-I)
      N3=N2P/2
      M=N
      DO 11 L=1,N3
      M=M/4
      M4=M*4
      LM4=N-M4+1
      AM4=M4
      AM5=SG*PI2/AM4
      DO 12 J=1,M
      JM1=J-1
      AJM1=JM1
      ARG=AJM1*AM5
      C1=DCOS(ARG)
      S1=DSIN(ARG)
      C2=C1*C1-S1*S1
      S2=C1*S1+C1*S1
      C3=C1*C2-S1*S2
      S3=C1*S2+C2*S1
      ZW1=DCMPLX(C1,S1)
      ZW2=DCMPLX(C2,S2)
      ZW3=DCMPLX(C3,S3)
      DO 13 I=1,LM4,M4
      J1=I+JM1
      J2=J1+M
      J3=J2+M
      J4=J3+M
      ZC1=Z(J1)+Z(J3)
      ZC2=Z(J1)-Z(J3)
      ZC3=Z(J2)+Z(J4)
      ZC4=Z(J2)-Z(J4)
      Z(J1)=ZC1+ZC3
      Z(J2)=(ZC1-ZC3)*ZW2
      ZC4=ZCI*ZC4
      Z(J3)=(ZC2+ZC4)*ZW1
      Z(J4)=(ZC2-ZC4)*ZW3
   13 CONTINUE
   12 CONTINUE
   11 CONTINUE
      N5=N2P-2*N3
      IF(N5.NE.1) GO TO 120
      NM1=N-1
      DO 110 I=1,NM1,2
      I1=I+1
      ZC=Z(I)+Z(I1)
      Z(I1)=Z(I)-Z(I1)
      Z(I)=ZC
  110 CONTINUE
C     UNSCRAMBLING
  120 JF=0
      DO 16 I=1,N
      IF(JF.LT.I) GO TO 17
      ZC=Z(I)
      Z(I)=Z(JF+1)
      Z(JF+1)=ZC
   17 DO 18 L=1,N2P
      LL=L
      IF(JF.LT.MS(L)) GO TO 19
   18 JF=JF-MS(L)
      LL=N2P
   19 JF=JF+MS(LL)
   16 CONTINUE
      IF(ISG.LT.0) GO TO 30
      DO 20 I=1,N
   20 Z(I)=Z(I)/AN
   30 RETURN
      END
C
c      SUBROUTINE MULPLY(X,Y,Z,MM,NN,NC,MJ1,MJ2,MJ3)
      SUBROUTINE MULPLY(X,Y,Z,MM,NN,NC)
C     MATRIX MULTIPLICATION
C     Z=X*Y
C     (UPPER LEFT MM X NC OF Z)=(UPPER LEFT MM X NN OF X)*(UPPER LEFT
C     NN X NC OF Y).
C     (MJ1,MJ2): ABSOLUTE DIMENSION OF X IN THE MAIN ROUTINE
C     (MJ2,MJ3): ABSOLUTE DIMENSION OF Y IN THE MAIN ROUTINE
C     (MJ1,MJ3): ABSOLUTE DIMENSION OF Z IN THE MAIN ROUTINE
      IMPLICIT REAL*8(A-H,O-Z)
c      DIMENSION X(MJ1,MJ2),Y(MJ2,MJ3),Z(MJ1,MJ3)
      DIMENSION X(MM,NN),Y(NN,NC),Z(MM,NC)
      CST0=0.0D-00
      DO 10 I=1,MM
      DO 11 J=1,NC
      SUM=CST0
      DO 12 K=1,NN
   12 SUM=SUM+X(I,K)*Y(K,J)
      Z(I,J)=SUM
   11 CONTINUE
   10 CONTINUE
      RETURN
      END
C
cc      SUBROUTINE MULVER(X,Y,Z,MM,NN,MJ1,MJ2)                            
      SUBROUTINE MULVER(X,Y,Z,MM,NN)
C     COMMON SUBROUTINE                                                 
C     Z=X*Y (X: MATRIX  Y,Z: VECTORS)                                   
C     (MJ1,MJ2): ABSOLUTE DIMENSION OF X IN THE MAIN ROUTINE            
      IMPLICIT REAL*8(A-H,O-Z)                                          
cc      DIMENSION X(MJ1,MJ2),Y(NN),Z(MM)                                  
      DIMENSION X(MM,NN),Y(NN),Z(MM)
      CST0=0.0D-00                                                      
      DO 10 I=1,MM                                                      
      SUM=CST0                                                          
      DO 11 J=1,NN                                                      
   11 SUM=SUM+X(I,J)*Y(J)                                               
   10 Z(I)=SUM                                                          
      RETURN                                                            
      END                                                               
C
cc	SUBROUTINE NEWSE(A1,SE,MS,K,MJ0,MJ)
      SUBROUTINE NEWSE(A1,CV,SE,MS,L,K,LCV1)
C     SE COMPUTATION
      IMPLICIT REAL*8(A-H,O-Z)
cc	DIMENSION A1(MJ0,MJ,MJ)
cc	DIMENSION SE(MJ,MJ)
cc	DIMENSION CV(25,7,7)
cc	DIMENSION A(7,7),R(7,7),Z(7,7)
      DIMENSION A1(L,K,K), CV(LCV1,K,K)
      DIMENSION SE(K,K)
      DIMENSION A(K,K),R(K,K),Z(K,K)
cc	COMMON /COM10/CV
      CST0=0.0D-00
      DO 10 II=1,K
      DO 10 JJ=1,K
   10 Z(II,JJ)=CST0
C
      MSP2=MS+2
      DO 11 I=1,MS
      MMI=MSP2-I
      DO 12 II=1,K
      DO 12 JJ=1,K
      A(II,JJ)=A1(I,II,JJ)
   12 R(II,JJ)=CV(MMI,II,JJ)
cc	CALL MULPLY(A,R,SE,K,K,K,MJ,MJ,MJ)
cc   11 CALL MATADL(Z,SE,K,K,MJ,MJ)
      CALL MULPLY(A,R,SE,K,K,K)
   11 CALL MATADL(Z,SE,K,K)
C
      DO 14 II=1,K
      DO 14 JJ=1,K
   14 R(II,JJ)=CV(MSP2,II,JJ)
cc	CALL SUBTAC(R,Z,SE,K,K,MJ,MJ)
      CALL SUBTAC(R,Z,SE,K,K)
      RETURN
      END
C
c      SUBROUTINE REARRA(X,INW,IP0,IP,MJ)
      SUBROUTINE REARRA(X,INW,IP0,IP)
C     SUBMATRIX REARRANGEMENT
C     X: ORIGINAL MATRIX
C     INW: INDICATOR OF ADOPTED ROWS
C     IP0: DIMENSION OF ORIGINAL MATRIX, SHOULD BE LESS THAN 11
C     IP: DIMENSION OF REARRANGED SUBMATRIX
C     MJ: ABSOLUTE DIMENSION OF X IN THE MAIN ROUTINE
C     THE REARRANGED SUBMATRIX IS OVERWRITTEN ON THE ORIGINAL.
C     NEXT STATEMENT SHOULD BE REPLACED BY
C     IMPLICIT COMPLEX*16(X)
C     FOR COMPLEX VERSION.
      IMPLICIT REAL*8(X)
c      DIMENSION X(MJ,MJ),INW(IP)
c      DIMENSION IOD(10)
      DIMENSION X(IP0,IP0),INW(IP)
      DIMENSION IOD(IP0)
C
      DO 300 I=1,IP0
  300 IOD(I)=I
      DO 301 I=1,IP
      I1=INW(I)
      I2=IOD(I1)
      IF(I.EQ.I2) GO TO 301
C     ROW INTERCHANGE
      DO 312 JJ=1,IP0
      XC=X(I,JJ)
      X(I,JJ)=X(I2,JJ)
  312 X(I2,JJ)=XC
C     COLUMN INTERCHANGE
      DO 314 II=1,IP0
      XC=X(II,I)
      X(II,I)=X(II,I2)
  314 X(II,I2)=XC
      ID=IOD(I)
      IOD(I2)=ID
      IOD(ID)=I2
  301 CONTINUE
      RETURN
      END
C
      SUBROUTINE REARRAC(X,INW,IP0,IP)
C     SUBMATRIX REARRANGEMENT
C     X: ORIGINAL MATRIX
C     INW: INDICATOR OF ADOPTED ROWS
C     IP0: DIMENSION OF ORIGINAL MATRIX, SHOULD BE LESS THAN 11
C     IP: DIMENSION OF REARRANGED SUBMATRIX
C     MJ: ABSOLUTE DIMENSION OF X IN THE MAIN ROUTINE
C     THE REARRANGED SUBMATRIX IS OVERWRITTEN ON THE ORIGINAL.
C     NEXT STATEMENT SHOULD BE REPLACED BY
      IMPLICIT COMPLEX*16(X)
c      DIMENSION X(MJ,MJ),INW(IP)
c      DIMENSION IOD(10)
      DIMENSION X(IP0,IP0),INW(IP)
      DIMENSION IOD(IP0)
C
      DO 300 I=1,IP0
  300 IOD(I)=I
      DO 301 I=1,IP
      I1=INW(I)
      I2=IOD(I1)
      IF(I.EQ.I2) GO TO 301
C     ROW INTERCHANGE
      DO 312 JJ=1,IP0
      XC=X(I,JJ)
      X(I,JJ)=X(I2,JJ)
  312 X(I2,JJ)=XC
C     COLUMN INTERCHANGE
      DO 314 II=1,IP0
      XC=X(II,I)
      X(II,I)=X(II,I2)
  314 X(II,I2)=XC
      ID=IOD(I)
      IOD(I2)=ID
      IOD(ID)=I2
  301 CONTINUE
      RETURN
      END
C
      SUBROUTINE SIGNIF(P1,P2,P3,LAGH1,N)
C     SIGNIFICANCE TEST
C     P1: SPECTRUM SMOOTHED BY WINDOW W1
C     P2: SPECTRUM SMOOTHED BY WINDOW W2
C     P3: TEST STATISTICS
C     LAGH1: DIMENSION OF PI (I=1,2,3)
C     N: LENGTH OF THE ORIGINAL DATA
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION P1(LAGH1),P2(LAGH1),P3(LAGH1)
      LAGH=LAGH1-1
      H=LAGH
      AN=N
      HAN=H/AN
      SD2=0.43D-00*DSQRT(HAN)
      SD3=1.0D-00/SD2
      DO 10 I=1,LAGH1
      T=P2(I)/P1(I)-1.0D-00
   10 P3(I)=DABS(T)*SD3
      RETURN
      END
C
      SUBROUTINE SMOSPE(X,LAGSHF,A,LA1,Z,LAGH1)
C     SPECTRUM SMOOTHING BY THE FORMULA
C     Z(I)=A(0)X(I)+A(1)(X(I+1)+X(I-1))+...+A(LA)(X(I+LA)+X(I-LA))
C     I=0,1,...,LAGH.
C     ACTUAL X(I) IS SHIFTED TO THE RIGHT BY LA FOR END CORRECTION.
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(LAGSHF),A(LA1),Z(LAGH1)
      LA=LA1-1
      DO 10 I=1,LAGH1
      I0=I+LA
      SUM1=0.0D-00
      DO 11 J=1,LA
      J1=I0-J
      J2=I0+J
   11 SUM1=SUM1+A(J+1)*(X(J1)+X(J2))
   10 Z(I)=A(1)*X(I0)+SUM1
      RETURN
      END
C
      SUBROUTINE SUBD12(N,LAGH,K,D1,D2)
C     CONSTANTS D1,D2 COMPUTATION
      IMPLICIT REAL*8(A-H,O-Z)
C     L1: NUMBER OF A(I)S (LESS THAN 5)
      DIMENSION A(4)
      CST0=0.0D-00
      L1=2
      A(1)=0.5D-00
      A(2)=0.25D-00
      AN=N
      H=LAGH
      SUM=0.0D-00
      DO 20 I=2,L1
   20 SUM=SUM+A(I)**2
      SUM=SUM+SUM+A(1)**2
      SUM=SUM+SUM
      NF=AN/(H*SUM)+0.5D-00
      FK=NF-K
      IF(FK.EQ.CST0) GO TO 100
      C1=FK-1.40D-00
      IF(C1.EQ.CST0) GO TO 100
      D1=(3.84D-00+10.0D-00/C1)/FK
      IF(D1.LT.CST0) GO TO 100
      D1=DSQRT(D1)
      GO TO 110
  100 D1=100.0D-00
  110 C2=FK+FK-1.40D-00
      IF(C2.EQ.CST0) GO TO 120
      D2=(3.0D-00+10.0D-00/C2)/FK
      IF(D2.LT.CST0) GO TO 120
      D2=DSQRT(D2)
      GO TO 130
  120 D2=100.0D-00
  130 RETURN
      END
C
      SUBROUTINE SUBDET(X,XDETMI,MM,MJ)                                 
C                                                                       
C       DETERMINANT OF X COMPUTATION                                    
C                                                                       
C     THIS SUBROUTINE COMPUTES THE DETERMINANT OF UPPER LEFT MM X MM    
C     OF X.  FOR GENERAL USE STATEMENTS 20-21 SHOULD BE RESTORED.       
C     X: ORIGINAL MATRIX                                                
C     XDETMI: DETERMINANT OF UPPER LEFT MM X MM OF X                    
C     MJ: ABSOLUTE DIMENSION OF X IN THE MAIN ROUTINE                   
      IMPLICIT REAL*8(X)                                                
      DIMENSION X(MJ,MJ)                                                
      XDETMI=1.0D0                                                      
      IF(MM.EQ.1) GO TO 18                                              
      MM1=MM-1                                                          
      DO 10 I=1,MM1                                                     
   20 IF(X(I,I).NE.0.0D0) GO TO 11                                      
      DO 12 J=I,MM                                                      
      IF(X(I,J).EQ.0.0D0) GO TO 12                                      
      JJ=J                                                              
      GO TO 13                                                          
   12 CONTINUE                                                          
      XDETMI=0.0D0                                                      
      GO TO 17                                                          
   13 DO 14 K=I,MM                                                      
      XXC=X(K,JJ)                                                       
      X(K,JJ)=X(K,I)                                                    
   14 X(K,I)=XXC                                                        
   21 XDETMI=-XDETMI                                                    
   11 XDETMI=XDETMI*X(I,I)                                              
      XC=1.0D0/X(I,I)                                                   
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
C
      SUBROUTINE SUBNOS(X,SD,IP,RS,R,MJ)
C     THIS SUBROUTINE COMPUTES RELATIVE POWER CONTRIBUTIONS.
C     MJ: ABSOLUTE DIMENSION OF X IN THE MAIN ROUTINE
C     IP: DIMENSION OF RS1 OR RL (LESS THAN 11)
      IMPLICIT REAL*8(A-H,O-W)
      IMPLICIT COMPLEX*16(X-Z)
      DIMENSION X(MJ,MJ)
      DIMENSION SD(MJ,MJ),RS(MJ,MJ),R(MJ,MJ)
c      DIMENSION RS1(10),RL(10)
      DIMENSION RS1(MJ),RL(MJ)
      CST0=0.0D-00
      CST1=1.0D-00
      DO 10 II=1,IP
      SUM=CST0
      DO 11 JJ=1,IP
      RX=DREAL(X(II,JJ))
      RIX=DIMAG(X(II,JJ))
      RS1(JJ)=(RX**2+RIX**2)*SD(JJ,JJ)
      SUM=SUM+RS1(JJ)
   11 RL(JJ)=SUM
      RCONST=CST1/RL(IP)
      DO 14 JJ=1,IP
   14 RS(II,JJ)=RS1(JJ)*RCONST
      DO 12 LL=1,IP
   12 R(II,LL)=RL(LL)*RCONST
   10 CONTINUE
      RETURN
      END
C
cc	SUBROUTINE SUBTAL(X,Y,MM,NN,MJ1,MJ2)
      SUBROUTINE SUBTAL(X,Y,MM,NN)
C     COMMON SUBROUTINE
C     MATRIX SUBTRACTION
C     X=X-Y
C     (UPPER LEFT MM X NN OF X)=(UPPER LEFT MM X NN OF X)-(UPPER LEFT
C     MM X NN OF Y).
C     (MJ1,MJ2): ABSOLUTE DIMENSION OF X AND Y IN THE MAIN ROUTINE
      IMPLICIT REAL*8(A-H,O-Z)
cc	DIMENSION X(MJ1,MJ2),Y(MJ1,MJ2)
      DIMENSION X(MM,NN),Y(MM,NN)
      DO 10 I=1,MM
      DO 10 J=1,NN
   10 X(I,J)=X(I,J)-Y(I,J)
      RETURN
      END
C
c      SUBROUTINE TRAMDL(X,Y,Z,MM,NN,NC,MJ1,MJ2,MJ3)
      SUBROUTINE TRAMDL(X,Y,Z,MM,NN,NC)
C     TRANSPOSE MULTIPLY (LEFT)
C     Z=X'*Y
C     (UPPER LEFT NN X NC OF Z)=(UPPER LEFT MM X NN OF X)'*(UPPER LEFT
C     MM X NC OF Y).
C     (MJ1,MJ2): ABSOLUTE DIMENSION OF X IN THE MAIN ROUTINE
C     (MJ1,MJ3): ABSOLUTE DIMENSION OF Y IN THE MAIN ROUTINE
C     (MJ2,MJ3): ABSOLUTE DIMENSION OF Z IN THE MAIN ROUTINE
      IMPLICIT REAL*8(A-H,O-Z)
c      DIMENSION X(MJ1,MJ2),Y(MJ1,MJ3),Z(MJ2,MJ3)
      DIMENSION X(MM,NN),Y(MM,NC),Z(NN,NC)
      CST0=0.0D-00
      DO 10 I=1,NN
      DO 11 J=1,NC
      SUM=CST0
      DO 12 K=1,MM
   12 SUM=SUM+X(K,I)*Y(K,J)
      Z(I,J)=SUM
   11 CONTINUE
   10 CONTINUE
      RETURN
      END
C
cc	SUBROUTINE TRAMDR(X,Y,Z,MM,NN,NC,MJ1,MJ2,MJ3)
      SUBROUTINE TRAMDR(X,Y,Z,MM,NN,NC)
C     COMMON SUBROUTINE
C     TRANSPOSE MULTIPLY (RIGHT)
C     Z=X*Y'
C     (UPPER LEFT MM X NC OF Z)=(UPPER LEFT MM X NN OF X)*(UPPER LEFT
C     NC X NN OF Y)'.
C     (MJ1,MJ2): ABSOLUTE DIMENSION OF X IN THE MAIN ROUTINE
C     (MJ3,MJ2): ABSOLUTE DIMENSION OF Y IN THE MAIN ROUTINE
C     (MJ1,MJ3): ABSOLUTE DIMENSION OF Z IN THE MAIN ROUTINE
      IMPLICIT REAL*8(A-H,O-Z)
cc	DIMENSION X(MJ1,MJ2),Y(MJ3,MJ2),Z(MJ1,MJ3)
      DIMENSION X(MM,NN),Y(NC,NN),Z(MM,NC)
      CST0=0.0D-00
      DO 10 I=1,MM
      DO 11 J=1,NC
      SUM=CST0
      DO 12 K=1,NN
   12 SUM=SUM+X(I,K)*Y(J,K)
      Z(I,J)=SUM
   11 CONTINUE
   10 CONTINUE
      RETURN
      END
C
C
      DOUBLE PRECISION FUNCTION DSUMF(X,N)
C     DOUBLE PRECISION SUMMATION
      DOUBLE PRECISION X
      DIMENSION X(N)
      DSUMF=0.0D-00
      DO 10 I=1,N
   10 DSUMF=DSUMF+X(I)
      RETURN
      END
C
cc	FUNCTION RANDOM(K)
ccc	 FUNCTION RANDM(K)
      DOUBLE PRECISION FUNCTION RANDM(K,K1,K2,K3,K4)
C     RANDOM NUMBER GENERATOR
      MCST11=11
      MC100=100
cc      IF(K) 1,2,1
      IF(K .LT. 0) GO TO 1
      IF(K .EQ. 0) GO TO 2
      IF(K .GT. 0) GO TO 1
C     STARTING NUMBER FOR GENERATOR
    1 K1=53
      K2=95
      K3=27
      K4=04
cc	WRITE(6,4) K1,K2,K3,K4
    2 M1=MCST11*K4
      M2=MCST11*K3
      M3=MCST11*K2+K4
      M4=MCST11*K1+K3
      J=M1/MC100
      K4=M1-MC100*J
      M2=M2+J
      J=M2/MC100
      K3=M2-MC100*J
      M3=M3+J
      J=M3/MC100
      K2=M3-MC100*J
      M4=M4+J
      J=M4/MC100
      K1=M4-MC100*J
      X1=K1
      X2=K2
cc	RANDOM=X1*1.E-2+X2*1.E-4
ccc      RANDM=X1*1.E-2+X2*1.E-4
      RANDM=X1*1.D-2+X2*1.D-4
      RETURN
    4 FORMAT(1H ,29HSTARTING NUMBER FOR GENERATOR,5X,4I3)
      END
C
C
C------------------------------------------------- TIMSAC74
      SUBROUTINE DINIT(A,N,DD)
      REAL*8 A(N),DD
      DO 100 I=1,N
	 A(I)=DD
  100 CONTINUE
      RETURN
      END
C
      SUBROUTINE INNERP(DD1,DD2,DINP12,INP)
C     COMMON SUBROUTINE
C     INNER-PRODUCT OF DD1 AND DD2.
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION DD1(INP),DD2(INP)
      CST0=0.0D-00
      SUM=CST0
      DO 100 I=1,INP
  100 SUM=SUM+DD1(I)*DD2(I)
      DINP12=SUM
      RETURN
      END
C
C
cc	SUBROUTINE MATINV(DET,M)
      SUBROUTINE MATINV(M,HS,NN,IFG,LU)
C      HS IS AN M*M MATRIX (IN COMMON AREA).
C     HS-INVERSE IS RETURNED IN HS.
C     DETERMINANT IS RETURNED IN DET.
      IMPLICIT REAL*8(A-H,O-Z)
cc	COMMON /COM50/HS
cc	DIMENSION HS(50,50)
      DIMENSION HS(NN,NN)
C
      CALL INVDET(HS,XDET,M,NN)
C
cc      CST0=0.0D-00
cc      CST1=1.0D-00
cc      DET=CST1
cc      DO 1 J=1,M
cc      PVT=HS(J,J)
cc      DET=DET*PVT
cc	WRITE(6,60) J,PVT,DET
cc      IF (IFG.NE.0) WRITE(LU,60) J,PVT,DET
cc      IF (PVT.GT.0) GO TO 5
cc	WRITE(6,61)
cc      IF (IFG.NE.0) WRITE(LU,61)
cc    5 HS(J,J)=CST1
cc      DO 2 K=1,M
cc    2 HS(J,K)=HS(J,K)/PVT
cc      DO 1 K=1,M
cc      IF(K-J) 3,1,3
cc    3 T=HS(K,J)
cc      HS(K,J)=CST0
cc      DO 4 L=1,M
cc    4 HS(K,L)=HS(K,L)-HS(J,L)*T
cc    1 CONTINUE
cc   60 FORMAT(1H ,'J=',I5,5X,'PVT=',D12.5,5X,'DET=',D12.5)
cc   61 FORMAT(1H ,'WARNING: NON POSITIVE DEFINITE HESSIAN')
      RETURN
      END
C
      SUBROUTINE MSVD(U,V,Q,M,N,MJ2,MJ1)
C     COMMON SUBROUTINE
C     THIS SUBROUTINE COMPLETES THE SINGULAR VALUE
C     DECOMPOSITION OF A REAL RECTANGULAR MATRIX A INTO THE FORM
C     A=U*DIAG(Q)*V' WITH ORTHGONAL MATRICES U AND V.
C     INPUTS:
C     U: ORIGINAL MATRIX A
C     M: NUMBER OF ROWS OF A, NOT LESS THAN N.
C     N: NUMBER OF COLUMNS OF A, NOT GREATER THAN 104.
C     OUTPUTS:
C     V: ORTHOGONAL MATRIX V
C     Q: SINGULAR VALUES IN DECREASING ORDER
      IMPLICIT REAL*8(A-H,O-Z)
cc	DIMENSION E(50)
      DIMENSION U(MJ2,MJ1),V(MJ1,MJ1)
      DIMENSION E(M)
      DIMENSION Q(N)
C     ESP: SMALL POSITIVE CONSTANT TO BE USED FOR THE DECISION OF CONVER
C     TOL: (SMALLEST POSITIVE NUMBER REPRESENTABLE IN THE COMPUTER)/EPS
      EPS=1.0D-15
      TOL=1.0D-60
C     HOUSEHOLDER'S REDUCTION TO BIDIAGONAL FORM
      CST0=0.0D-00
      CST1=1.0D-00
      CST2=2.0D-00
      G=CST0
      X=CST0
      DO 10 I=1,N
      E(I)=G
      S=CST0
      L=I+1
      DO 11 J=I,M
   11 S=S+U(J,I)**2
      IF(S.GE.TOL) GO TO 12
      G=CST0
      GO TO 19
   12 F=U(I,I)
   13 IF(F.LT.CST0) GO TO 14
      G=-DSQRT(S)
      GO TO 15
   14 G=DSQRT(S)
   15 H=F*G-S
      U(I,I)=F-G
      IF(L.GT.N) GO TO 19
      DO 16 J=L,N
      S=CST0
      DO 17 K=I,M
   17 S=S+U(K,I)*U(K,J)
      F=S/H
      DO 18 K=I,M
   18 U(K,J)=U(K,J)+F*U(K,I)
   16 CONTINUE
   19 Q(I)=G
      S=CST0
      IF(L.GT.N) GO TO 201
      DO 20 J=L,N
   20 S=S+U(I,J)**2
  201 IF(S.GE.TOL) GO TO 22
      G=CST0
      GO TO 30
   22 F=U(I,I+1)
   23 IF(F.LT.CST0) GO TO 24
      G=-DSQRT(S)
      GO TO 25
   24 G=DSQRT(S)
   25 H=F*G-S
      U(I,I+1)=F-G
      IF(L.GT.N) GO TO 30
      DO 26 J=L,N
   26 E(J)=U(I,J)/H
      DO 27 J=L,M
      S=CST0
      DO 28 K=L,N
   28 S=S+U(J,K)*U(I,K)
      DO 29 K=L,N
cc   29 U(J,K)=U(J,K)+S*E(K)
      U(J,K)=U(J,K)+S*E(K)
   29 CONTINUE
   27 CONTINUE
   30 Y=DABS(Q(I))+DABS(E(I))
C     X=MAX.Y
      IF(Y.GT.X) X=Y
   10 CONTINUE
      NP1=N+1
	L=N+1
C     ACCUMULATION OF RIGHT-HAND TRANSFORMATIONS
      DO 110 II=1,N
      I=NP1-II
      IF(L.GT.N) GO TO 202
      IF(G.EQ.CST0) GO TO 115
      H=U(I,I+1)*G
      DO 111 J=L,N
  111 V(J,I)=U(I,J)/H
      DO 112 J=L,N
      S=CST0
      DO 113 K=L,N
  113 S=S+U(I,K)*V(K,J)
      DO 114 K=L,N
  114 V(K,J)=V(K,J)+S*V(K,I)
  112 CONTINUE
  115 DO 116 J=L,N
      V(I,J)=CST0
  116 V(J,I)=CST0
  202 V(I,I)=CST1
      G=E(I)
      L=I
  110 CONTINUE
C     DIAGONALIZATION OF THE BIDIAGONAL FORM
  129 EPS=EPS*X
      DO 130 KK=1,N
      K=NP1-KK
      KP1=K+1
C     TEST F SPLITTING
  145 DO 131 LL=1,K
      L=KP1-LL
      IF(DABS(E(L)).LE.EPS) GO TO 134
      IF(DABS(Q(L-1)).LE.EPS) GO TO 132
  131 CONTINUE
C     CANCELLATION OF E(L) IF L>1.
  132 C=CST0
      S=CST1
      L1=L-1
      DO 133 I=L,K
      F=S*E(I)
      E(I)=C*E(I)
      IF(DABS(F).LE.EPS) GO TO 134
      G=Q(I)
      H=DSQRT(F*F+G*G)
      Q(I)=H
      C=G/H
      S=-F/H
  133 CONTINUE
C     TEST F CONVERGENCE
  134 Z=Q(K)
      IF(L.EQ.K) GO TO 142
C     SHIFT FROM BOTTOM 2X2 MINOR
      X=Q(L)
      Y=Q(K-1)
      G=E(K-1)
      H=E(K)
      F=((Y-Z)*(Y+Z)+(G-H)*(G+H))/(CST2*H*Y)
      T=F*F+CST1
      G=DSQRT(T)
      IF(F.GE.CST0) GO TO 135
      FG=F-G
      GO TO 136
  135 FG=F+G
  136 F=((X-Z)*(X+Z)+H*(Y/FG-H))/X
C     QR TRANSFORMATION
C     NEXT QR TRANSFORMATION
      C=CST1
      S=CST1
      LP1=L+1
      DO 137 I=LP1,K
      G=E(I)
      Y=Q(I)
      H=S*G
      G=C*G
      Z=DSQRT(F*F+H*H)
      E(I-1)=Z
      C=F/Z
      S=H/Z
      F=X*C+G*S
      G=-X*S+G*C
      H=Y*S
      Y=Y*C
      DO 138 J=1,N
      X=V(J,I-1)
      Z=V(J,I)
      V(J,I-1)=X*C+Z*S
  138 V(J,I)=-X*S+Z*C
  139 Z=DSQRT(F*F+H*H)
      Q(I-1)=Z
      C=F/Z
      S=H/Z
      F=C*G+S*Y
      X=-S*G+C*Y
  137 CONTINUE
      E(L)=CST0
      E(K)=F
      Q(K)=X
      GO TO 145
C     CONVERGENCE
  142 IF(Z.GE.CST0) GO TO 130
C     Q(K) IS MADE NON-NEGATIVE.
      Q(K)=-Z
      DO 144 J=1,N
  144 V(J,K)=-V(J,K)
  130 CONTINUE
C     Q AND V ARE ORDERED IN DECREASING ORDER OF THE SINGULAR VALUES.
      IF(N.LE.1) GO TO 340
      NM1=N-1
      DO 311 I=1,NM1
      IP1=I+1
      DO 312 J=IP1,N
      IF(Q(I).GE.Q(J)) GO TO 312
      T=Q(I)
      Q(I)=Q(J)
      Q(J)=T
      DO 321 L=1,N
      T=V(L,I)
      V(L,I)=V(L,J)
      V(L,J)=T
  321 CONTINUE
  312 CONTINUE
  311 CONTINUE
  340 RETURN
      END
C

C
cc	SUBROUTINE SUBTAC(X,Y,Z,MM,NN,MJ1,MJ2)
      SUBROUTINE SUBTAC(X,Y,Z,MM,NN)
C     COMMON SUBROUTINE
C     MATRIX SUBTRACTION
C     Z=X-Y
C     (UPPER LEFT MM X NN OF Z)=(UPPER LEFT MM X NN OF X)-(UPPER LEFT
C     MM X NN OF Y).
C     (MJ1,MJ2): ABSOLUTE DIMENSION OF X, Y AND Z IN THE MAIN ROUTINE
      IMPLICIT REAL*8(A-H,O-Z)
cc	DIMENSION X(MJ1,MJ2),Y(MJ1,MJ2),Z(MJ1,MJ2)
      DIMENSION X(MM,NN),Y(MM,NN),Z(MM,NN)
      DO 10 I=1,MM
      DO 10 J=1,NN
   10 Z(I,J)=X(I,J)-Y(I,J)
      RETURN
      END
C
C-------------------------------------------------------------------
C
cc      SUBROUTINE  ADDVAR( X,D,IND,JND,K,L,M,MJ )                        
      SUBROUTINE  ADDVAR( X,IND,JND,K,L,M,MJ )                        
C                                                                       
C         +-----------------------------------------------------------+ 
C         ! ADDITION OF THE VARIABLE M AS THE JJ-TH REGRESSOR (JJ<=L) ! 
C         +-----------------------------------------------------------+ 
C                                                                       
C       ----------------------------------------------------------------
C       THE FOLLOWING SUBROUTINE IS DIRECTLY CALLED BY THIS SUBROUTINE: 
C             HUSHL1                                                    
C       ----------------------------------------------------------------
C       INPUTS:                                                         
C          X:            (K+1)*(K+1)  MATRIX                            
C          D:            WORKING AREA                                   
C          IND:          WORKING AREA                                   
C          JND(I)=J:     PRESENT STATUS, I-TH REGRESSOR IS VARIABLE J.  
C          K:            NUMBER OF VARIABLES                            
C          L:            NUMBER OF REGRESSORS IN THE PRESENT MODEL      
C          M:            INDICATION OF THE VARIABLE TO BE ADDED         
C          MJ:           ABSOLUTE DIMENSION OF X                        
C                                                                       
C       OUTPUTS:                                                        
C          X:            (K+1)*(K+1)  MATRIX                            
C          JND(I)=J:     UPDATED STATUS, I-TH REGRESSOR IS VARIABLE J.  
C                                                                       
      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
cc      DIMENSION  X(MJ,1) , D(1) , IND(1) , JND(1)                       
      DIMENSION  X(MJ,1) , IND(1) , JND(1)                       
C                                                                       
      K1 = K + 1                                                        
      DO  60     I=1,K1                                                 
      J = JND(I)                                                        
   60 IND(J) = I                                                        
      JJ = IND(M)                                                       
cc      IF( JJ-L )     40,40,10                                           
      IF( JJ-L .LT. 0 )  GO TO 40
      IF( JJ-L .EQ. 0 )  GO TO 40
      IF( JJ-L .GT. 0 )  GO TO 10
   10 II1 = L + 1                                                       
      DO  20     I=II1,JJ                                               
      I1 = JJ+L-I                                                       
   20 JND(I1+1) = JND(I1)                                               
      JND(L) = M                                                        
C                                                                       
cc      CALL  HUSHL1( X,D,MJ,K1,JJ,L,IND,JND )                            
      CALL  HUSHL1( X,MJ,K1,JJ,L,IND,JND )                            
C                                                                       
   30 L = L + 1                                                         
   40 RETURN                                                            
C                                                                       
      END                                                               
c
c
      SUBROUTINE  AICCOM( X,N,M,K,MJ,SD,AIC )                           
C                                                                       
C     THIS SUBROUTINE COMPUTES INNOVATION VARIANCE AND AIC OF THE MODEL 
C     WITH M REGRESSORS.                                                
C                                                                       
      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
      DIMENSION  X(MJ,1)                                                
C                                                                       
C       INPUTS:                                                         
C          X:     (K+1)*(K+1) REDUCED MATRIX                            
C          N:     DATA LENGTH                                           
C          M:     NUMBER OF REGRESSORS OF THE MODEL                     
C          K:     POSSIBLE LARGEST NUMBER OF REGRESSORS                 
C          MJ:    ABSOLUTE DIMENSION OF X                               
C                                                                       
C       OUTPUTS:                                                        
C          SD:    INNOVATION VARIANCE OF THE MODEL                      
C          AIC:   AIC OF THE MODEL                                      
C                                                                       
      FN = N                                                            
      M1 = M + 1                                                        
      K1 = K + 1                                                        
      SUM = 0.D0                                                        
      DO 10  I=M1,K1                                                    
   10 SUM = SUM + X(I,K1)*X(I,K1)                                       
      SD = SUM / FN                                                     
   20 AIC = FN*DLOG( SD ) + 2.D0*M                                      
      RETURN                                                            
      END                                                               
c
c
cc      SUBROUTINE  ARBAYS( X,D,K,LAG,N,ISW,TITLE,MJ1,A,B,SDB,AICB )      
      SUBROUTINE  ARBAYS( X,D,K,LAG,N,ISW,MJ1,SD,AIC,DIC,AICM,SDMIN,
     *IMIN,A,B1,B,C,SDB,PN,AICB )      
C                                                                       
C         +-----------------------------------------------+             
C         ! AUTOREGRESSIVE MODEL FITTING (BAYESIAN MODEL) !             
C         +-----------------------------------------------+             
C                                                                       
C     THIS SUBROUTINE PRODUCES AN AUTOREGRESSIVE MODEL BY A BAYESIAN    
C     PROCEDURE USING THE OUTPUT OF SUBROUTINE REDUCT.                  
C                                                                       
C       ----------------------------------------------------------------
C       THE FOLLOWING SUBROUTINES ARE DIRECTLY CALLED BY THIS SUBROUTINE
C             MAICE                                                     
C             ARCOEF                                                    
C             BAYSPC                                                    
C             BAYSWT                                                    
C             COMAIC                                                    
C             PRINTA                                                    
C             SDCOMP                                                    
C       ----------------------------------------------------------------
C       INPUTS:                                                         
C         X:    (K+1)*(K+1) TRIANGULAR MATRIX, OUTPUT OF SUBROUTINE REDU
C         D:    WORKING AREA                                            
C         K:    HIGHEST ORDER OF THE MODELS                             
C         N:    DATA LENGTH                                             
C         ISW:  =0     TO SUPPRESS THE OUTPUTS                          
C               =1     TO PRINT OUT THE OUTPUTS                         
C         TITLE:     TITLE OF DATA                                      
C         MJ1:  ABSOLUTE DIMENSION OF X                                 
C                                                                       
C       OUTPUTS:                                                        
C         A(I) (I=1,K):     AR-COEFFICIENTS                             
C         B(I) (I=1,K):     PARTIAL AUTOCORRELATION COEFFICIENTS        
C         SDB:  INNOVATION VARIANCE OF BAYESIAN MODEL                   
C         AICB: AIC OF BAYESIAN MODEL                                   
C                                                                       
      IMPLICIT REAL * 8( A-H,O-Z )                                      
cc      REAL * 4   TITLE(20) , TTL(6)                                     
cc      DIMENSION  X(MJ1,1) , D(1) , A(1) , B(1)                          
cc      DIMENSION  SD(101) , AIC(101) , C(101)                            
      DIMENSION  X(MJ1,1), D(1), A(1), B(1), B1(1)
      DIMENSION  SD(K+1), AIC(K+1), C(K+1), DIC(K+1)
cc      DATA       TTL / 4H..  ,4H BAY,4HESIA,4HN MO,4HDEL ,4H  .. /      
C                                                                       
cc      IF(ISW.GE.1)  WRITE( 6,2 )                                        
      FN = N                                                            
C                                                                       
C         +-----------------------------------------+                   
C         ! INNOVATION VARIANCE AND AIC COMPUTATION !                   
C         +-----------------------------------------+                   
C                                                                       
      CALL  COMAIC( X,N,K,MJ1,SD,AIC )                                  
C                                                                       
C         +-------------+                                               
C         ! AIC DISPLAY !                                               
C         +-------------+                                               
C                                                                       
cc      CALL  MAICE( AIC,SD,K,ISW,AICM,SDMIN,IMIN )                       
      CALL MAICE( AIC,SD,K,ISW,AICM,SDMIN,IMIN,DIC )                   
C                                                                       
C         +-----------------------------+                               
C         ! BAYSIAN WEIGHTS COMPUTATION !                               
C         +-----------------------------+                               
C                                                                       
      CALL  BAYSWT( AIC,AICM,K,0,C )                                    
C                                                                       
C         +-----------------------------------------------------+       
C         ! PARTIAL AUTOCORRELATION ESTIMATION (BAYESIAN MODEL) !       
C         +-----------------------------------------------------+       
C                                                                       
cc      CALL  BAYSPC( X,C,N,K,ISW,MJ1,B,D )                               
      CALL  BAYSPC( X,C,N,K,ISW,MJ1,B1,B,D )
C                                                                       
C         +-----------------------------------------+                   
C         ! AUTOREGRESSIVE COEFFICIENTS COMPUTATION !                   
C         +-----------------------------------------+                   
C                                                                       
      CALL  ARCOEF( B,K,A )                                             
C                                                                       
C                                                                       
C          EQUIVALENT NUMBER OF PARAMETERS                              
C                                                                       
      PN = 1.D0                                                         
      DO   10    I=1,K                                                  
   10 PN = PN + D(I)*D(I)                                               
C                                                                       
C         +-------------------------------------------+                 
C         ! INNOVATION VARIANCE OF THE BAYESIAN MODEL !                 
C         +-------------------------------------------+                 
C                                                                       
cc      CALL  SDCOMP( X,A,D,N,K,MJ1,SDB )                                 
      CALL  SDCOMP( X,A,N,K,MJ1,SDB )                                 
C                                                                       
C          -------------------------                                    
C          AIC OF THE BAYESIAN MODEL                                    
C          -------------------------                                    
C                                                                       
      AICB = FN * DLOG( SDB ) + 2.D0*PN                                 
C                                                                       
cc      IF( ISW .EQ. 0 )     RETURN                                       
C                                                                       
cc      WRITE( 6,12 )                                                     
cc      WRITE( 6,9 )     (B(I),I=1,K)                                     
cc      WRITE( 6,15 )     SDB , PN , AICB                                 
cc      L1 = LAG + 1                                                      
cc      NL = N + LAG                                                      
cc      IF( ISW .GE. 1 )   CALL  PRINTA( A,SDB,K,TTL,6,TITLE,L1,NL )      
C                                                                       
      RETURN                                                            
C                                                                       
    2 FORMAT( //1H ,18(1H-),/,' BAYESIAN PROCEDURE',/,1H ,18(1H-) )     
    5 FORMAT( 1H ,4X,'I',10X,'SD',19X,'AIC',11X,'WEIGHTS' )             
    6 FORMAT( 1H ,I5,D20.10,F16.3,F16.5 )                               
    7 FORMAT( 1H ,'*****  MINIMUM AIC =',D17.10,3X,'ATTAINED AT M =',I3,
     1'  *****' )                                                       
    8 FORMAT( 1H ,'M =',I3,5X,'WEIGHT =',F8.5 )                         
    9 FORMAT( 1H ,10F13.8 )                                             
   11 FORMAT( //,1H ,132(1H-),/,' AR-COEFFICIENTS  (UP TO WEIGHT = 0.01)
     1' )                                                               
   12 FORMAT( ///' < BAYESIAN MODEL >',/,9H PARCOR'S )                  
   13 FORMAT( ' AR-COEFFICIENTS' )                                      
   15 FORMAT( 1H ,'INNOVATION VARIANCE',13X,'=',D24.8,/,1H ,'EQUIVALENT 
     1NUMBER OF PARAMETERS =',F15.3,/,' EQUIVALENT AIC',18X,'=',F15.3 ) 
   16 FORMAT( 1H+,32X,'***  MAICE MODEL  ***' )                         
C                                                                       
      E N D                                                             
c
c
      SUBROUTINE  ARCOEF( B,K,A )                                       
C                                                                       
C     THIS SUBROUTINE PRODUCES THE AUTOREGRESSIVE COEFFICIENTS FROM A SE
C     OF PARTIAL AUTOCORRELATION COEFFICIENTS.                          
C                                                                       
C       INPUTS:                                                         
C         B:    VECTOR OF PARTIAL AUTOCORRELATIONS                      
C         K:    ORDER OF THE MODEL                                      
C                                                                       
C       OUTPUTS:                                                        
C         A:    VECTOR OF AR-COEFFICIENTS                               
C                                                                       
C                                                                       
      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
cc      DIMENSION  B(1) , A(1) , AA(100)                                  
      DIMENSION  B(K) , A(K) , AA(K)
C                                                                       
      DO  30     II=1,K                                                 
      A(II) = B(II)                                                     
      AA(II) = B(II)                                                    
      IM1 = II - 1                                                      
      IF( IM1 .LE. 0 )     GO TO 30                                     
      DO  10     J=1,IM1                                                
      JJ = II - J                                                       
   10 A(J) = AA(J) - B(II)*AA(JJ)                                       
      IF( II .EQ. K )     GO TO 40                                      
      DO  20     J=1,IM1                                                
   20 AA(J) = A(J)                                                      
   30 CONTINUE                                                          
   40 CONTINUE                                                          
      RETURN                                                            
      END                                                               
c
c
cc      SUBROUTINE  ARMFIT( X,K,LAG,N,ISW,TITLE,MJ1,A,SDMIN,IMIN )        
      SUBROUTINE  ARMFIT( X,K,LAG,N,ISW,MJ1,A,IMIN,SD,AIC,DIC,SDMIN,
     *                    AICM,IFG,LU )
C                                                                       
C          +------------------------------+                             
C          ! AUTOREGRESSIVE MODEL FITTING !                             
C          +------------------------------+                             
C                                                                       
C       THIS SUBROUTINE PRODUCES PARAMETERS OF AR-MODELS USING THE OUTPU
C       SUBROUTINE REDUCT.                                              
C                                                                       
C       ----------------------------------------------------------------
C       THE FOLLOWING SUBROUTINES ARE DIRECTLY CALLED BY THIS SUBROUTINE
C             MAICE                                                     
C             COMAIC                                                    
C             PRINTA                                                    
C             RECOEF                                                    
C       ----------------------------------------------------------------
C       INPUTS:                                                         
C         X:    (K+1)*(K+1) TRIANGULAR MATRIX, OUTPUT OF SUBROUTINE REDU
C         K:    NUMBER OF REGRESSORS OF THE MODEL                       
C         LAG:  MAXIMUM TIME LAG OF THE MODEL                           
C         N:    DATA LENGTH                                             
C         ISW:  =0   TO PRODUCE THE MAICE MODEL ONLY (OUTPUTS SUPPRESSED
C               =1   TO PRODUCE THE MAICE MODEL ONLY                    
C               =2   TO PRODUCE ALL AR-MODELS (UP TO THE ORDER K)       
C         TITLE:     TITLE OF DATA                                      
C         MJ1:  ABSOLUTE DIMENSION OF X                                 
C                                                                       
C       OUTPUTS:                                                        
C         A:      VECTOR OF MAICE AR-COEFFICIENTS                       
C         SDMIN:  MAICE INNOVATION VARIANCE                             
C         IMIN:   MAICE ORDER                                           
C                                                                       
      IMPLICIT REAL * 8( A-H,O-Z )                                      
cc      REAL * 4   TITLE(20) , TTL(4)                                     
cc      DIMENSION  X(MJ1,1) , A(1) , SD(101) , AIC(101)                   
      DIMENSION  X(MJ1,1) , A(K)
      DIMENSION  SD(K+1), AIC(K+1), DIC(K+1)
cc      DATA     TTL / 4H   M,4H A I ,4H C E,4H   . /                    
c
cc      IF(ISW.GE.1)  WRITE( 6,3 )                                        
      IF( (ISW.GE.2) .AND. (IFG.NE.0) )  WRITE( LU,3 )
C                                                                       
C          +-----------------------------------------+                +-
C          ! INNOVATION VARIANCE AND AIC COMPUTATION !                ! 
C          +-----------------------------------------+                +-
C                                                                       
      CALL  COMAIC( X,N,K,MJ1,SD,AIC )                                  
C                                                                       
C          +-------------+                                            +-
C          ! AIC DISPLAY !                                            ! 
C          +-------------+                                            +-
C                                                                       
cc      CALL  MAICE( AIC,SD,K,ISW,AICM,SDMIN,IMIN )                       
      CALL  MAICE( AIC,SD,K,ISW,AICM,SDMIN,IMIN,DIC )
C                                                                       
C          +-----------------------------------------+                +-
C          ! AUTOREGRESSIVE COEFFICIENTS COMPUTATION !                ! 
C          +-----------------------------------------+                +-
C                                                                       
      IF( ISW .LT. 2 )     GO TO  20                                    
C                                                                       
cc           WRITE( 6,5 )                                                 
      IF ( IFG.NE.0 )  WRITE( LU,5 )
      M = 0                                                             
cc      WRITE( 6,7 )     M , SD(1)                                        
      IF ( IFG.NE.0 )  WRITE( LU,7 )  M , SD(1)
      DO  10     M=1,K                                                  
      M1 = M + 1                                                        
      CALL  RECOEF( X,M,K,MJ1,A )                                       
cc   10 WRITE( 6,6 )     M , SD(M1) , (A(I),I=1,M)                        
   10 IF ( IFG.NE.0 )  WRITE( LU,6 )  M, SD(M1), (A(I),I=1,M)
C                                                                       
   20 IF( IMIN .GE. 1 )  CALL  RECOEF( X,IMIN,K,MJ1,A )                 
      L1 = LAG + 1                                                      
      NL = N + LAG                                                      
cc      IF( ISW .GE. 1 )   CALL  PRINTA( A,SDMIN,IMIN,TTL,4,TITLE,L1,NL ) 
C                                                                       
      RETURN                                                            
C                                                                       
    3 FORMAT( // 1H ,15(1H-),/,' MAICE PROCEDURE',/,1H ,15(1H-))        
    5 FORMAT( //' AR-COEFFICIENTS' )                                    
    6 FORMAT( ' M =',I3,9X,'(SD(M) =',D15.8,1X,')',/,(1X,10F13.8) )     
    7 FORMAT( ' M =',I3,9X,'(SD(M) =',D15.8,1X,')' )                    
C                                                                       
      E N D                                                             
c
cc      SUBROUTINE  BAYSPC( X,C,N,K,ISW,MJ1,B,D )                         
      SUBROUTINE  BAYSPC( X,C,N,K,ISW,MJ1,B,B1,D )
C                                                                       
C     THIS SUBROUTINE PRODUCES PARTIAL AUTOCORRELATION COEFFICIENTS B(I)
C     (I=1,K) OF THE BAYESIAN MODEL.                                    
C                                                                       
C         +-----------------------------------------------------+       
C         ! PARTIAL AUTOCORRELATION ESTIMATION (BAYESIAN MODEL) !       
C         +-----------------------------------------------------+       
C                                                                       
C       INPUTS:                                                         
C         X:    N*(K+1) TRIANGULAR MATRIX, OUTPUT OF SUBROUTINE REDUCT  
C         C:    BAYESIAN WEIGHTS, OUTPUT OF SUBROUTINE BAYSWT           
C         N:    DATA LENGTH                                             
C         K:    HIGHEST ORDER OF THE MODELS                             
C         ISW:  =0     OUTPUTS ARE SUPPRESSED                           
C         ISW:  >0     OUTPUTS ARE PRINTED OUT                          
C         MJ1:  ABSOLUTE DIMENSION OF X                                 
C                                                                       
C       OUTPUTS:                                                        
C         B:    VECTOR OF PARTIAL AUTOCORRELATIONS                      
C         D:    INTEGRATED BAYESIAN WEIGHTS                             
C                                                                       
      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
cc      REAL * 4  TABLE(41) , AST , AXIS , TERM                           
      CHARACTER TABLE(41) , AST , AXIS , TERM                           
      DIMENSION  X(MJ1,1) , B(1) , C(1) , D(1)                          
      DIMENSION  B1(1)
cc      DATA  AST , AXIS / 1H* , 1H! /                                    
cc      DATA  TABLE / 41*1H /                                             
      DATA  AST , AXIS / '*' , '!' /                                    
      DATA  TABLE/ 41*' ' /                                             
      K1 = K + 1                                                        
      FN = N                                                            
C                                                                       
C          PARTIAL AUTOCORRELATION (LEAST SQUARES ESTIMATE)             
C                                                                       
      SUM = X(K1,K1)**2                                                 
      DO 10  I=1,K                                                      
      J = K1 - I                                                        
      SUM = SUM + X(J,K1)**2                                            
      G = DSQRT( SUM )                                                  
   10 B(J) = X(J,K1)*X(J,J) / (G*DABS(X(J,J)))                          
C                                                                       
C          INTEGRATED BAYESIAN WEIGHT                                   
C                                                                       
      D(K) = C(K1)                                                      
      DO  80     I=2,K                                                  
      J = K1 - I                                                        
   80 D(J) = D(J+1) + C(J+1)                                            
      IF( ISW .EQ. 0 )   GO TO 100                                      
C                                                                       
C          PARCOR AND BAYESIAN WEIGHTS DISPLAY                          
C                                                                       
      SC = 20.D0                                                        
      SIG = 1.D0/DSQRT(FN)                                              
      ISIG = SIG*SC + 0.5D0                                             
      J0 = SC + 1                                                       
      ISIGP1 = J0 + ISIG                                                
      ISIGM1 = J0 - ISIG                                                
      ISIGP2 = J0 + ISIG*2                                              
      ISIGM2 = J0 - ISIG*2                                              
      TABLE(ISIGP1) = AXIS                                              
      TABLE(ISIGM1) = AXIS                                              
      TABLE(ISIGP2) = AXIS                                              
      TABLE(ISIGM2) = AXIS                                              
      TABLE(J0) = AXIS                                                  
C                                                                       
cc      WRITE( 6,11 )                                                     
cc      WRITE( 6,12 )                                                     
      DO  90     I=1,K                                                  
      BB = B(I)*SC                                                      
      IB = BB
      IF( BB .GT. 0.D0 )     IB = BB + 0.5D0                            
      IF( BB .LT. 0.D0 )     IB = BB - 0.5D0                            
      IB = IB + J0                                                      
      TERM = TABLE(IB)                                                  
      TABLE(IB) = AST                                                   
cc      WRITE( 6,13 )     I , B(I) , C(I+1) , D(I)                        
cc      WRITE( 6,14 )     (TABLE(J),J=1,41)                               
      TABLE(IB) = TERM                                                  
   90 CONTINUE                                                          
C                                                                       
C          PARCOR OF BAYESIAN MODEL                                     
C                                                                       
  100 CONTINUE                                                          
      DO  110     I=1,K                                                 
cc  110 B(I) = B(I) * D(I)                                                
  110 B1(I) = B(I) * D(I)
C                                                                       
      RETURN                                                            
C                                                                       
cc   11 FORMAT( ///,73X,'PARCOR (LINES SHOW +SD AND +2SD), SD=SQRT(1/N)',/
cc     1,1H+,91X,'_',7X,'_' )                                             
cc   12 FORMAT( 1H ,12X,'PARCOR',9X,'WEIGHT',7X,'INTEGRATED',17X,'-1',19X,
cc     1'0',19X,'1',/,1H ,4X,'I',8X,'B(I)',11X,'W(I)',10X,'WEIGHT',19X,'+'
cc     2,4(10H---------+) )                                               
cc   13 FORMAT( 1H ,I5,3F15.8 )                                           
cc   14 FORMAT( 1H ,67X,41A1 )                                            
   11 FORMAT(//65X,'PARCOR (LINES SHOW +/-SD AND +/-2SD), SD=SQRT(1/N)')
   12 FORMAT( 1H ,12X,'PARCOR',9X,'WEIGHT',7X,'INTEGRATED',16X,'-1',19X,
     1'0',19X,'1',/,1H ,4X,'I',8X,'B(I)',11X,'W(I)',10X,'WEIGHT',19X,'+'
     2,4(10H---------+) )                                               
   13 FORMAT( 1H ,I5,3F15.8,17X,41A1 )
C                                                                       
      END                                                               
c
      SUBROUTINE  BAYSWT( AIC,AICM,K,ISW,C )                            
C                                                                       
C     THIS SUBROUTINE COMPUTES BAYESIAN WEIGHT OF AR-MODEL OF EACH ORDER
C                                                                       
C         +-----------------------------+                               
C         ! BAYESIAN WEIGHT COMPUTATION !                               
C         +-----------------------------+                               
C                                                                       
C     BAYESIAN WEIGHT OF THE M-TH ORDER MODEL IS DEFINED BY             
C            W(M)  =  CONST * C(M) / ( M+1 ),                           
C     WHERE                                                             
C            CONST =  NORMALIZING FACTOR                                
C            C(M)  =  EXP( -0.5*AIC(M) )                                
C                                                                       
C       INPUTS:                                                         
C         AIC:  VECTOR OF AIC'S                                         
C         AICM: MINIMUM AIC                                             
C         K:    HIGHEST ORDER OF THE MODELS                             
C         ISW:  =0   NON-ADAPTIVE DAMPER                                
C               =1   DATA ADAPTIVE DAMPER                               
C               =2   DAMPER DOES NOT USED                               
C                                                                       
C       OUTPUT:                                                         
C         C(I) (I=1,K+1):   VECTOR OF BAYESIAN WEIGHTS                  
C                                                                       
      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
      DIMENSION  AIC(1) , C(1)                                          
C                                                                       
C          C(I)  =  EXP( -AIC(I)/2 )                                    
C                                                                       
      K1 = K + 1                                                        
      SUM = 0.D0                                                        
      EK = 0.D0                                                         
      DO  10     I=1,K1                                                 
      DIC = -0.5D0 * (AIC(I) - AICM)                                    
      C(I) = 0.D0                                                       
      IF( DIC .LT. -40.D0 )     GO TO 10                                
      C(I) = DEXP(DIC)                                                  
      EK = EK + (I-1) * C(I)                                            
   10 SUM = SUM + C(I)                                                  
C                                                                       
C          DAMPING OF C(I)                                              
C                                                                       
      IF( ISW .EQ. 1 )     GO TO 30                                     
      IF( ISW .EQ. 2 )     GO TO 50                                     
      SUM = 0.D0                                                        
      DO  20     I=1,K1                                                 
      AI = I                                                            
      C(I) = C(I) / AI                                                  
   20 SUM = SUM + C(I)                                                  
      GO TO 50                                                          
C                                                                       
C          DATA ADAPTIVE DAMPER                                         
C                                                                       
   30 EK = EK / ( SUM+EK )                                              
      SUM = 0.D0                                                        
      DO  40     I=1,K1                                                 
      C(I) = C(I) * EK**(I-1)                                           
   40 SUM = SUM + C(I)                                                  
C                                                                       
C          NORMALIZATION OF C                                           
C                                                                       
   50 DO  60     I=1,K1                                                 
   60 C(I) = C(I) / SUM                                                 
C                                                                       
      RETURN                                                            
C                                                                       
      END                                                               
c
c
      SUBROUTINE  COEF2( A,M,ID,II,JND,LMAX,MM,KSW,MSW,MJ1,B,C,E )      
C                                                                       
C     COMPOSITION OF AR-COEFFICIENT MATRICES WITH INSTANTANEOUS RESPONSE
C                                                                       
C       INPUTS:                                                         
C          A(I)  (I=1,M):     REGRESSION COEFFICIENTS                   
C          JND(I) (I=1,M):    SPECIFICATION OF REGRESSORS               
C          M:     NUMBER OF REGRESSORS                                  
C          ID:    DIMENSION OF VECTOR OF OBSERVATIONS                   
C          II:    REGRESSAND                                            
C          LMAX:  HIGHEST ORDER OF THE MODEL                            
C          MM:    ORDER                                                 
C          KSW:   =0  CONSTANT VECTOR IS NOT INCLUDED AS A REGRESSOR    
C                 =1  CONSTANT VECTOR IS INCLUDED AS THE FIRST REGRESSOR
C          MSW:   =0  CONSTANT VECTOR IS NOT ADOPTED AS A REGRESSOR     
C                 =1  CONSTANT VECTOR IS ACTUALLY ADOPTED AS REGRESSOR  
C          MJ1:   ABSOLUTE DIMENSION OF B                               
C       OUTPUTS:                                                        
C          B:     AR-COEFFICIENT MATRICES                               
C          C:     CONSTANT VECTOR                                       
C          E:     COEFFICIENT FOR INSTANTANEOUS RESPONSE                
C                                                                       
      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
      DIMENSION  A(1) , JND(1) , C(1)                                   
      DIMENSION  B(MJ1,MJ1,1) , E(MJ1,1)                                
C                                                                       
      M0 = MSW + 1                                                      
      C(II) = 0.D0                                                      
      IF( MSW .EQ. 1 )     C(II) = A(1)                                 
      DO  100   JJ=M0,M                                                 
      I = JND(JJ) - KSW                                                 
      L = I / ID                                                        
      J = I - L*ID                                                      
      IF( J .NE. 0 )     GO TO 10                                       
      J = ID                                                            
      L = L - 1                                                         
   10 L1 = L + 1                                                        
      IF( I .LE. MM*ID )     GO TO  20                                  
      E(II,J) = -A(JJ)                                                  
      GO TO 30                                                          
   20 B(II,J,L1) = A(JJ)                                                
      IF( LMAX .LT. L1 )     LMAX = L1                                  
   30 CONTINUE                                                          
  100 CONTINUE                                                          
      DO  40     I=1,ID                                                 
   40 E(I,I) = 1.D0                                                     
      RETURN                                                            
      E N D                                                             
c
c
      SUBROUTINE  COMAIC( X,N,K,MJ1,SD,AIC )                            
C                                                                       
C          +-----------------------------------------+                  
C          ! INNOVATION VARIANCE AND AIC COMPUTATION !                  
C          +-----------------------------------------+                  
C                                                                       
C       INPUTS:                                                         
C          X:     (K+1)*(K+1) TRIANGULAR MATRIX, OUTPUT OF SUBROUTINE  R
C          N:     DATA LENGTH                                           
C          K:     NUMBER OF REGRESSORS                                  
C          MJ1:   ABSOLUTE DIMENSION OF X                               
C                                                                       
C       OUTPUTS:                                                        
C          SD(M)  (M=1,K+1):   VECTOR OF INNOVATION VARIANCES           
C          AIC(M) (M=1,K+1):   VECTOR OF AIC'S.                         
C                                                                       
      IMPLICIT  REAL * 8( A-H,O-Z )                                     
      DIMENSION  X(MJ1,1) , AIC(1) , SD(1)                              
      FN = N                                                            
      K1 = K + 1                                                        
C                                                                       
      OSD = 0.0D00                                                      
      DO 10   I = 1,K1                                                  
      M = K1 - I + 1                                                    
      OSD = OSD + X(M,K1)**2                                            
      SD(M)  = OSD / FN                                                 
   10 AIC(M) = FN*DLOG( SD(M) ) + 2*M                                   
C                                                                       
      RETURN                                                            
C                                                                       
      E N D                                                             
c
c
      SUBROUTINE  COPY( X,K,II,JJ,MJ1,MJ2,Y )                           
C                                                                       
C         +-----------------------+                                     
C         ! MAKE A COPY OF X ON Y !                                     
C         +-----------------------+                                     
C                                                                       
C        INPUTS:                                                        
C             X(I+II,J):     K*K MATRIX (I,J=1,...,K)                   
C             II:            INDICATES ORIGIN OF X                      
C             JJ:            INDICATES ORIGIN OF Y                      
C             MJ1:           ABSOLUTE DIMENSION OF X                    
C             MJ2:           ABSOLUTE DIMENSION OF Y                    
C                                                                       
C        OUTPUT:                                                        
C             Y(I+JJ,J):     COPY OF X (I,J=1,...,K)                    
C                                                                       
      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
      DIMENSION  X(MJ1,1) , Y(MJ2,1)                                    
C                                                                       
      DO  10     I=1,K                                                  
      I1 = I + II                                                       
      I2 = I + JJ                                                       
      DO  10     J=1,K                                                  
   10 Y(I2,J) = X(I1,J)                                                 
C                                                                       
      RETURN                                                            
C                                                                       
      END                                                               
c
c
cc      SUBROUTINE  DELETE( X,D,IND,JND,K,L,M,MJ )                        
      SUBROUTINE  DELETE( X,IND,JND,K,L,M,MJ )                        
C                                                                       
C         +------------------------------------------------+            
C         ! DELETION OF THE VARIABLE M FROM THE REGRESSORS !            
C         +------------------------------------------------+            
C                                                                       
C       ----------------------------------------------------------------
C       THE FOLLOWING SUBROUTINE IS DIRECTLY CALLED BY THIS SUBROUTINE: 
C             HUSHL1                                                    
C       ----------------------------------------------------------------
C       INPUTS:                                                         
C          IND(J)=I:     PRESENT STATUS,   VARIABLE J IS THE I-TH REGRES
C          JND(I)=J:     REQUIRED STATUS,  I-TH REGRESSOR IS VARIABLE J.
C          X:            (K+1)*(K+1)  MATRIX                            
C          D:            WORKING AREA                                   
C          K:            NUMBER OF VARIABLES                            
C          L:            NUMSER OF REGRESSORS IN THE PRESENT MODEL      
C          M:            INDICATION OF THE VARIABLE TO BE DELETED       
C          MJ:           ABSOLUTE DIMENSION OF X                        
C       OUTPUTS:                                                        
C          X:            (K+1)*(K+1)  MATRIX                            
C                                                                       
      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
cc      DIMENSION  X(MJ,1) , D(1) , IND(1) , JND(1)                       
      DIMENSION  X(MJ,1) , IND(1) , JND(1)                       
C                                                                       
      K1 = K + 1                                                        
      DO  60     I=1,K1                                                 
      J = JND(I)                                                        
   60 IND(J) = I                                                        
      II = IND(M)                                                       
cc      IF( II-L )     10,30,40                                           
      IF( II-L .LT. 0 )  GO TO 10
      IF( II-L .EQ. 0 )  GO TO 30
      IF( II-L .GT. 0 )  GO TO 40
   10 I1 = II + 1                                                       
      DO  20     I=I1,L                                                 
   20 JND(I-1) = JND(I)                                                 
      JND(L) = M                                                        
      LM1 = L - 1                                                       
C                                                                       
cc      CALL  HUSHL1( X,D,MJ,K1,LM1,II,IND,JND )                          
      CALL  HUSHL1( X,MJ,K1,LM1,II,IND,JND )                          
C                                                                       
   30 L = L - 1                                                         
C                                                                       
   40 RETURN                                                            
      END                                                               
C
C
      SUBROUTINE FOUGER(G,LGP1,FC,FS,LF1)                               
C                                                                       
C     FOURIER TRANSFORM (GOERTZEL METHOD)                               
C     THIS SUBROUTINE COMPUTES FOURIER TRANSFORM OF G(I),I=0,1,...,LG AT
C     FREQUENCIES K/(2*LF), K=0,1,...,LF AND RETURNS COSINE TRANSFORM IN
C     FC(K) AND SINE TRANSFORM IN FS(K).                                
C                                                                       
C       INPUTS:                                                         
C          G(I):   ORIGINAL DATA (I=0,1,...,LG)                         
C          LG1:    = LG+1                                               
C          LF1:    = LF+1                                               
C                                                                       
C       OUTPUTS:                                                        
C          FC(I):  COSINE TRANSFORM OF G  (I=0,1,...,LF)                
C          FB(I):  SINE TRANSFORM OF G  (I=0,1,...,LF)                  
C                                                                       
      IMPLICIT REAL*8(A-H,O-Z)                                          
      DIMENSION G(LGP1),FC(LF1),FS(LF1)                                 
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
  100 G(I2)=T                                                           
  110 PI=3.1415926536D0                                                 
      ALF=LF                                                            
      T=PI/ALF                                                          
      DO 10 K=1,LF1                                                     
      AK=K-1                                                            
      TK=T*AK                                                           
      CK=DCOS(TK)                                                       
      SK=DSIN(TK)                                                       
      CK2=CK+CK                                                         
      UM1=0.0D0                                                         
      UM2=0.0D0                                                         
      IF(LG.EQ.0) GO TO 12                                              
      DO 11 I=1,LG                                                      
      UM0=CK2*UM1-UM2+G(I)                                              
      UM2=UM1                                                           
   11 UM1=UM0                                                           
   12 FC(K)=CK*UM1-UM2+G(LGP1)                                          
      FS(K)=SK*UM1                                                      
   10 CONTINUE                                                          
      RETURN                                                            
      END                                                               
C
C
cc      SUBROUTINE  HUSHL1( X,D,MJ1,K,L,M,IND,JND )                       
      SUBROUTINE  HUSHL1( X,MJ1,K,L,M,IND,JND )                       
C                                                                       
C     THIS SUBROUTINE PERFORMS THE HOUSEHOLDER TRANSFORMATION OF THE MAT
C                                                                       
C       INPUTS:                                                         
C         X:   ORIGINAL (K+1)*(K+1) MATRIX                              
C         N:   NUMBER OF ROWS OF X,  NOT GREATER THAN MJ1               
C         K:   NUMBER OF COLUMNS OF X                                   
C         L:   END POSITION OF THE HOUSEHOLDER TRANSFORMATION           
C         M:   STARTING POSITION OF THE HOUSEHOLDER TRANSFORMATION      
C         IND:   SPECIFICATION OF THE PRESENT FORM OF X                 
C                IND(J) = I;     VARIABLE I IS THE J-TH REGRESSOR       
C         JND:   SPECIFICATION OF THE REQUIRED FORM OF X                
C                JND(I) = J;     THE I-TH REGRESSOR IS VARIABLE J       
C                                                                       
C       OUTPUTS:                                                        
C         X:   TRANSFORMED MATRIX                                       
C                                                                       
C                                                                       
      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
      DIMENSION  X(MJ1,1) , D(MJ1) , IND(1) , JND(1)                    
C                                                                       
      TOL = 1.0D-60                                                     
C                                                                       
      NN = 0                                                            
      DO  100     II=M,L                                                
      JJ = JND(II)                                                      
      NN = MAX0( NN,IND(JJ) )                                           
      H = 0.0D00                                                        
      DO  10     I=II,NN                                                
      D(I) = X(I,JJ)                                                    
   10 H = H + D(I)*D(I)                                                 
      IF( H .GT. TOL )     GO TO 20                                     
      G = 0.0D00                                                        
      GO TO 100                                                         
   20 G = DSQRT( H )                                                    
      F = X(II,JJ)                                                      
      IF( F .GE. 0.0D00 )     G = -G                                    
      D(II) = F - G                                                     
      H = H - F * G                                                     
C                                                                       
C     ( I - D*D'/H ) * X                                                
C                                                                       
      II1 = II + 1                                                      
      IF( II .EQ. NN )   GO TO 35                                       
      DO  30     I=II1,NN                                               
   30 X(I,JJ) = 0.D0                                                    
   35 CONTINUE                                                          
      IF( II .EQ. K )     GO TO 100                                     
      DO  60     J1=II1,K                                               
      J = JND(J1)                                                       
      S = 0.0D00                                                        
      DO  40     I=II,NN                                                
   40 S = S + D(I)*X(I,J)                                               
      S = S / H                                                         
      DO  50     I=II,NN                                                
   50 X(I,J) = X(I,J) - D(I)*S                                          
   60 CONTINUE                                                          
  100 X(II,JJ) = G                                                      
      RETURN                                                            
      E N D                                                             
C
C
cc      SUBROUTINE  HUSHLD( X,D,MJ1,N,K )                                 
      SUBROUTINE  HUSHLD( X,MJ1,N,K )                                 
C                                                                       
C          +----------------------------+                               
C          ! HOUSEHOLDER TRANSFORMATION !                               
C          +----------------------------+                               
C                                                                       
C     THIS SUBROUTINE  TRANSFORMS MATRIX X INTO AN UPPER TRIANGULAR FORM
C     BY HOUSEHOLDER TRANSFORMATION.                                    
C                                                                       
C       INPUTS:                                                         
C          X:     ORIGINAL N*K DATA MATRIX                              
C          D:     WORKING AREA                                          
C          MJ1:   ABSOLUTE DIMENSION OF X                               
C          N:     NUMBER OF ROWS OF X, NOT GREATER THAN MJ1             
C          K:     NUMBER OF COLUMNS OF X                                
C                                                                       
C       OUTPUT:                                                         
C          X:     SQUARE ROOT OF DATA COVARIANCE MATRIX (UPPER TRIANGULA
C                                                                       
C                                                                       
      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
      DIMENSION  X(MJ1,1) , D(MJ1)                                      
C                                                                       
           TOL = 1.0D-60                                                
C                                                                       
      DO 100  II=1,K                                                    
         H = 0.0D00                                                     
         DO 10  I=II,N                                                  
            D(I) = X(I,II)                                              
   10       H = H + D(I)*D(I)                                           
         IF( H .GT. TOL )  GO TO 20                                     
         G = 0.0D00                                                     
         GO TO 100                                                      
   20    G = DSQRT( H )                                                 
         F = X(II,II)                                                   
         IF( F .GE. 0.0D00 )   G = -G                                   
         D(II) = F - G                                                  
         H = H - F*G                                                    
C                                                                       
C          FORM  (I - D*D'/H) * X, WHERE H = D'D/2                      
C                                                                       
         II1 = II+1                                                     
         DO 30  I=II1,N                                                 
   30    X(I,II) = 0.0D00                                               
         IF( II .EQ. K )  GO TO 100                                     
         DO 60  J=II1,K                                                 
            S = 0.0D00                                                  
            DO 40  I=II,N                                               
   40       S = S + D(I)*X(I,J)                                         
            S = S/H                                                     
            DO 50  I=II,N                                               
   50      X(I,J) = X(I,J) - D(I)*S                                     
   60    CONTINUE                                                       
  100 X(II,II) = G                                                      
C                                                                       
      RETURN                                                            
C                                                                       
      E N D                                                             
C

C
cc      SUBROUTINE  MAICE( AIC,SD,K,ISW,AICM,SDM,IMIN )                   
      SUBROUTINE  MAICE( AIC,SD,K,ISW,AICM,SDM,IMIN,DIC )
C                                                                       
C             +-------------+                                           
C             ! AIC DISPLAY !                                           
C             +-------------+                                           
C                                                                       
C        THIS SUBROUTINE PRODUCES NUMERICAL AND GRAPHICAL DISPLAYS OF AI
C                                                                       
C       INPUTS:                                                         
C          AIC:   VECTOR OF AIC'S                                       
C          SD:    VECTOR OF INNOVATION VARIANCES                        
C          K:     UPPER LIMIT OF THE ORDER                              
C          ISW:   =0   OUTPUTS ARE SUPPRESSED                           
C                 >0   AIC'S ARE DIPLAIED                               
C                                                                       
C       OUTPUTS:                                                        
C          AICM:  MINIMUM AIC                                           
C          SDM:   MAICE INNOVATION VARIANCE                             
C          IMIN:  MAICE ORDER                                           
C                                                                       
      IMPLICIT  REAL * 8( A-H,O-Z )                                     
cc      REAL * 4  TBL(41) , AST , BLNK , AXIS , PERI                      
      DIMENSION  AIC(1) , SD(1)                                         
      DIMENSION  DIC(K+1)
cc      DATA  AST / 1H* / , BLNK / 1H  / , AXIS / 1H! / , PERI / 1H. /    
C                                                                       
cc      DO  10     I=1,41                                                 
cc   10 TBL(I) = BLNK                                                     
C                                                                       
C       SEARCH FOR THE MINIMUM OF AIC(I)                                
C                                                                       
      K1 = K + 1                                                        
      IMIN = 0                                                          
      SDM  = SD(1)                                                      
      AICM = AIC(1)                                                     
      DO  20   I = 1,K                                                  
      IF( AIC(I+1) .GE. AICM )  GO TO 20                                
      IMIN = I                                                          
      SDM = SD(I+1)                                                     
      AICM = AIC(I+1)                                                   
   20 CONTINUE                                                          
      DO 25    I = 1,K1
   25 DIC(I) = AIC(I) - AICM
cc      IF( ISW .LE. 0 )     RETURN                                       
C                                                                       
C       DISPLAY OF AIC'S                                                
C                                                                       
cc      WRITE( 6,5 )                                                      
cc      DO  30   I = 1,K1                                                 
cc      II = I - 1                                                        
cc      DIC = AIC(I) - AICM                                               
cc      TBL(1) = AXIS                                                     
cc      ID = DIC + 1.5D0                                                  
cc      ID = DIC(I) + 1.5D0                                                  
cc      IF( ID .LE. 41 )     TBL(ID) = AST                                
cc      IF( ID .GT. 41 )     TBL(41) = PERI                               
cc      WRITE( 6,6 )     II , SD(I) , AIC(I) , DIC , (TBL(J),J=1,41)      
cc      IF( ID .LE. 41 )     TBL(ID) = BLNK                               
cc      IF( ID .GT. 41 )     TBL(41) = BLNK                               
cc   30 CONTINUE                                                          
cc      WRITE( 6,7 )     AICM , IMIN , SDM                                
C                                                                       
      RETURN                                                            
C                                                                       
    5 FORMAT( /1H ,70X,'AIC(M)-AICMIN (TRUNCATED AT 40.0)',/,2X,'ORDER',
     1 2X,'INNOVATION VARIANCE',40X,'0',8X,'10',8X,'20',8X,'30',8X,'40',
     2/,5X,'M',10X,'SD(M)',14X,'AIC(M)',7X,'AIC(M)-AICMIN',7X,'+',4(10H-
     3--------+) )                                                      
    6 FORMAT( 1H ,I5,D20.10,2F16.3,10X,41A1 )                           
    7 FORMAT( /1H ,'*****  MINIMUM AIC =',D17.10,3X,'ATTAINED AT M =',
     1I3,5X,'SD(M) =',D17.10,'  *****' )
C                                                                       
      E N D                                                             
C
C
      SUBROUTINE  MARCOF( D,E,ID,M,MJ3,A,B )                            
C                                                                       
C     THIS SUBROUTINE COMPUTES COEFFICIENT MATRICES OF MULTI-VARIATE AUT
C     REGRESSIVE MODEL FROM PARTIAL AUTOREGRESSION COEFFICIENT MATRICES 
C     FORWARD AND BACKWARD MODELS.                                      
C                                                                       
C       INPUTS:                                                         
C          D:     PARTIAL AUTOCORRELATIONS OF FORWARD MODEL             
C          E:     PARTIAL AUTOCORRELATIONS OF BACKWARD MODEL            
C          ID:    DIMENSION OF THE PROCESS                              
C          M:     ORDER OF THE MODELS                                   
C          MJ3:   ABSOLUTE DIMENSION OF A,B,D AND E IN THE MAIN PROGRAM 
C                                                                       
C       OUTPUTS:                                                        
C          A:     AR-COEFFICIENT MATRICES OF FORWARD MODEL              
C          B:     AR-COEFFICIENT MATRICES OF BACKWARD MODEL             
      IMPLICIT  REAL * 8 ( A-H , O-Z )                                  
      DIMENSION  A(MJ3,MJ3,1) , B(MJ3,MJ3,1)                            
      DIMENSION  D(MJ3,MJ3,1) , E(MJ3,MJ3,1)                            
cc      DIMENSION  F(10,10) , G(10,10)                                    
      DIMENSION  F(ID,ID) , G(ID,ID)                                    
C                                                                       
      DO 10  II=1,M                                                     
      DO 10  I=1,ID                                                     
      DO 10  J=1,ID                                                     
      A(I,J,II) = D(I,J,II)                                             
   10 B(I,J,II) = E(I,J,II)                                             
      IF(  M .EQ. 1 )  RETURN                                           
C                                                                       
      DO 60  II=2,M                                                     
C                                                                       
      II1 = II - 1                                                      
      DO 50  JJ=1,II1                                                   
      IMJ = II - JJ                                                     
      DO 20  I=1,ID                                                     
      DO 20  J=1,ID                                                     
      F(I,J) = A(I,J,IMJ)                                               
   20 G(I,J) = B(I,J,JJ)                                                
      DO 40  I=1,ID                                                     
      DO 40  J=1,ID                                                     
      SUMA = F(I,J)                                                     
      SUMB = G(I,J)                                                     
      DO 30  L=1,ID                                                     
      SUMA = SUMA - A(I,L,II)*G(L,J)                                    
   30 SUMB = SUMB - B(I,L,II)*F(L,J)                                    
      A(I,J,IMJ) = SUMA                                                 
   40 B(I,J,JJ) = SUMB                                                  
   50 CONTINUE                                                          
C                                                                       
   60 CONTINUE                                                          
      RETURN                                                            
C                                                                       
      E N D                                                             
C
C
cc      SUBROUTINE  MARFIT( X,Y,D,N,ID,M,KSW,MJ1,MJ2,MJ3,MJ4,ISW,IPR,B,E, 
cc     *                    EX,C,LMAX,AICS )                              
      SUBROUTINE  MARFIT( X,N,ID,M,KSW,MJ1,MJ2,MJ3,MJ4,ISW,IPR,AIC,SD,
     *DIC,AICM,SDM,IM,BI,EI,B,E,EX,C,LMAX,AICS,JNDF,AF,NPR,AAIC,IFG,LU )
C                                                                       
C         MULTI-VARIATE AUTOREGRESSIVE MODEL FITTING                    
C       ----------------------------------------------------------------
C       THE FOLLOWING SUBROUTINES ARE DIRECTLY CALLED BY THIS SUBROUTINE
C             COPY                                                      
C             COEF2                                                     
C             MAICE                                                     
C             MCOEF                                                     
C             ADDVAR                                                    
C             AICCOM                                                    
C             DELETE                                                    
C             HUSHL1                                                    
C             SRCOEF                                                    
C       ----------------------------------------------------------------
C                                                                       
C       INPUTS:                                                         
C          X:     ((M+1)*ID)*((M+1)*ID)  UPPER TRIANGULAR MATRIX,       
C                 OUTPUT OF SUBROUTINE MREDCT                           
C          Y:     WORKING AREA (MATRIX)                                 
C          D:     WORKING AREA                                          
C          N:     DATA LENGTH                                           
C          ID:    DIMENSION OF DATA                                     
C          M:     HIGHEST ORDER OF THE MODELS                           
C          KSW:   =0  CONSTANT VECTOR IS NOT INCLUDED AS A REGRESSOR    
C                 =1  CONSTANT VECTOR IS INCLUDED AS THE FIRST REGRESSOR
C          MJ1:   ABSOLUTE DIMENSION OF X IN THE MAIN PROGRAM           
C          MJ2:   ABSOLUTE DIMENSION OF E IN THE MAIN PROGRAM           
C          MJ3:   ABSOLUTE DIMENSION OF B IN THE MAIN PROGRAM           
C          MJ4:   ABSOLUTE DIMENSION OF Y IN THE MAIN PROGRAM           
C          ISW:   =0  MULTI-VARIATE AUTOREGRESSIVE MODEL IS REQUESTED   
C                 =1  INSTANTANEOUS RESPONSE MODEL IS REQUESTED         
C          IPR:   PRINT OUT CONTROL                                     
C                                                                       
C       OUTPUTS:                                                        
C          B:     AR-COEFFICIENT MATRICES                               
C          E:     INNOVATION VARIANCE MATRIX                            
C          EX:    RESIDUAL VARIANCES OF INSTANTENEOUS RESPONSE MODELS   
C          C:     CONSTANT VECTOR                                       
C          LMAX:  ORDER OF THE MAICE MODEL                              
C                                                                       
      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
cc      DIMENSION  X(MJ1,1) , Y(MJ4,1) , D(1) , E(MJ2,1) , B(MJ2,MJ2,MJ3) 
cc      DIMENSION  C(1) , A(100) , AIC(51) , SD(51) , EX(1)               
cc      DIMENSION  IND(100) , JND(100) , KND(100)                         
      DIMENSION  X(MJ1,1), Y(MJ4,MJ4), E(ID,ID), B(ID,ID,M) 
      DIMENSION  C(1), A(MJ4)
      DIMENSION  AIC(M+1,ID), SD(M+1,ID), DIC(M+1,ID), EX(ID)
      DIMENSION  AICM(ID), SDM(ID), IM(ID)
      DIMENSION  IND(MJ4) , JND(MJ4) , KND(MJ4)                         
      DIMENSION  EI(ID,ID) , BI(ID,ID,M)
      DIMENSION  JNDF(MJ4,ID), AF(MJ4,ID)
      DIMENSION  NPR(ID), AAIC(ID)
C                                                                       
C                                                                       
C         INITIAL SETTING                                               
C                                                                       
      M1 = M + 1                                                        
      MD0 = M * ID                                                      
      MD = M*ID + KSW                                                   
      MD2 = M1*ID + KSW                                                 
      AICSUM = 0.D0                                                     
      LMAX = 0                                                          
      NSW = 0                                                           
      DO 20  I=1,MJ2                                                    
      DO 20  J=1,MJ2                                                    
      DO 10  K=1,MJ3                                                    
   10 B(I,J,K) = 0.D0                                                   
   20 E(I,J) = 0.D0                                                     
      DO 30  I=1,MD2                                                    
   30 JND(I) = I                                                        
C                                                                       
      CALL  COPY( X,MD2,0,MD2,MJ1,MJ1,X )                               
C                                                                       
C                                                                       
      DO 500     II=1,ID                                                
      MSW = KSW                                                         
C                                                                       
cc      IF(IPR.GE.3)  WRITE( 6,3 )                                        
cc      IF( IPR.GE.2)   WRITE( 6,645 )   II                               
cc      IF(IPR.GE.3)  WRITE( 6,642 )                                      
      IF( (IPR.GE.2) .AND. (IFG.NE.0) )  WRITE( LU,645 )   II
      IF( (IPR.GE.3) .AND. (IFG.NE.0) )  WRITE( LU,642 )
      JJ = II - 1                                                       
      KK = MD + JJ                                                      
      KK1 = KK + 1                                                      
C                                                                       
C         ADDITION OF REGRESSOR (INSTANTANEOUS RESPONSE FROM II-TH VARIA
C                                                                       
      CALL  COPY( X,MD2,MD2,0,MJ1,MJ4,Y )                               
      CALL  COPY( X,MD2,MD2,0,MJ1,MJ1,X )                               
C                                                                       
      IF( II .LE. 1 )  GO TO 40                                         
      J1 = JJ + KSW                                                     
      DO 70  I=1,MD2                                                    
      J = JND(I)                                                        
   70 IND(J) = I                                                        
      JND(1) = 1                                                        
      DO 75  I=1,JJ                                                     
      J = MD + I                                                        
      I1 = I + KSW                                                      
   75 JND(I1) = J                                                       
      DO 80  I=1,MD0                                                    
      J = I + KSW                                                       
      I1 = I + KSW + JJ                                                 
   80 JND(I1) = J                                                       
      I1 = JJ + 1                                                       
      DO 85  I=I1,ID                                                    
      J = MD + I                                                        
   85 JND(J) = J                                                        
C                                                                       
cc      CALL  HUSHL1( X,D,MJ1,MD2,MD2,1,IND,JND )                         
      CALL  HUSHL1( X,MJ1,MD2,MD2,1,IND,JND )                         
C                                                                       
      CALL  COPY( X,MD2,0,MD2,MJ1,MJ1,X )                               
C                                                                       
C--------------------   FIRST STEP OF AIC MINIMIZATION   ---------------
C                                                                       
C          AIC'S OF INITIAL MODELS COMPUTATION                          
C                                                                       
   40 DO 50  I=1,M1                                                     
      K = (I-1)*ID + JJ + KSW                                           
      CALL  AICCOM( X,N,K,KK,MJ1,OSD,OAIC )                             
cc      SD(I) = OSD                                                       
cc   50 AIC(I) = OAIC                                                     
      SD(I,II) = OSD                                                    
   50 AIC(I,II) = OAIC                                                  
      DO 60  I=1,MD2                                                    
   60 KND(I) = JND(I)                                                   
C                                                                       
C         ORDER DETERMINATION BY AIC ( INITIAL ESTMATE )                
C                                                                       
      IPR2 = IPR - 2                                                    
cc      CALL  MAICE( AIC,SD,M,IPR2,AICMIN,SDMIN,IMIN )                    
      CALL MAICE( AIC(1,II),SD(1,II),M,IPR2,AICMIN,SDMIN,IMIN,
     *            DIC(1,II) )
      AICM(II) = AICMIN
      SDM(II) = SDMIN
      IM(II) = IMIN
C                                                                       
      K0 = IMIN*ID + JJ + KSW                                           
C                                                                       
C                                                                       
C         REGRESSION COEFFICIENTS COMPUTATION ( INITIAL ESTIMATE )      
C                                                                       
      IF( IPR .LT. 3 )  GO TO 90                                        
cc      WRITE( 6,4 )                                                      
cc      CALL  SRCOEF( X,K0,KK,N,MJ1,JND,IPR,A,SDD )                       
cc      WRITE( 6,643 )                                                    
      IF ( IFG.NE.0 )  WRITE( LU,4 )
      CALL  SRCOEF( X,K0,KK,N,MJ1,JND,IPR,A,SDD,AAIC(II),IFG,LU )
      IF ( IFG.NE.0 )  WRITE( LU,643 )
   90 CONTINUE                                                          
C                                                                       
C                                                                       
C--------------------   SECOND STEP OF AIC MINIMIZATION   --------------
C                                                                       
      IMP1 = IMIN + 1                                                   
      DO 230  I1=1,ID                                                   
      DO 100  I=1,KK                                                    
  100 KND(I) = JND(I)                                                   
C                                                                       
      CALL  COPY( X,KK1,0,0,MJ1,MJ4,Y )                                 
cc      IF(IPR.GE.3)  WRITE( 6,11 )     I1                                
      IF( (IPR.GE.3) .AND. (IFG.NE.0) )  WRITE( LU,11 )     I1
      AIC1 = AICMIN                                                     
      AIC2 = AICMIN                                                     
      K01 = K0                                                          
cc      IMP1 = IMIN + 1                                                   
cc      IF(IPR.GE.3)  WRITE( 6,13 )     AIC1                              
      IF( (IPR.GE.3) .AND. (IFG.NE.0) )  WRITE( LU,13 )     AIC1
      IF( IMIN .GE. M )  GO TO 140                                      
C                                                                       
C          CHECK REGRESSOR MADD  < ADD? >                               
C                                                                       
      DO 110  J1=IMP1,M                                                 
      MADD = (J1-1)*ID + I1 + KSW                                       
      IF( IND(MADD) .EQ. 1 )  GO TO 110                                 
      K = K01 + 1                                                       
      KX1 = K                                                           
cc      CALL  ADDVAR( X,D,IND,JND,KK,KX1,MADD,MJ1 )                       
      CALL  ADDVAR( X,IND,JND,KK,KX1,MADD,MJ1 )                       
      CALL  AICCOM( X,N,K,KK,MJ1,OSD,OAIC )                             
C                                                                       
C         DECISION BY AIC                                               
C                                                                       
cc      IF(IPR.GE.3)  WRITE( 6,5 )   MADD,OAIC,MADD,J1,I1                 
      IF( OAIC .GT. AIC1 )  GO TO 120                                   
cc      IF( IPR .GE. 3 )     WRITE( 6,6 )                                 
      IF((IPR.GE.3).AND.(IFG.NE.0)) WRITE( LU,56 ) MADD,OAIC,MADD,J1,I1
      AIC1 = OAIC                                                       
      K01 = K01 + 1                                                     
  110 CONTINUE                                                          
      GO TO 140                                                         
  120 K = K - 1                                                         
cc      IF( IPR .GE. 3 )     WRITE( 6,9 )                                 
      IF((IPR.GE.3).AND.(IFG.NE.0)) WRITE( LU,59 ) MADD,OAIC,MADD,J1,I1
      IF( J1 .NE. IMP1 )  GO TO 140                                     
      CALL  COPY( Y,KK1,0,0,MJ4,MJ1,X )                                 
      DO 130  J=1,KK                                                    
      I = IND(J)                                                        
  130 JND(I) = J                                                        
  140 CONTINUE                                                          
C                                                                       
C         CHECK REGRESSOR MDEL  < DELETE ? >                            
C                                                                       
      IF( IMIN .LT. 1 )  GO TO 200                                      
      K02 = K0                                                          
      DO 170  JM1=1,IMIN                                                
      J1 = IMIN - JM1 + 1                                               
      MDEL = (J1-1)*ID + I1 + KSW                                       
      K = K02 - 1                                                       
      DO 150  I0=1,MD2                                                  
      IF( JND(I0) .EQ. MDEL )  GO TO 160                                
  150 CONTINUE                                                          
      GO TO 170                                                         
  160 CONTINUE                                                          
cc      CALL  DELETE( Y,D,IND,KND,KK,K02,MDEL,MJ4 )                       
      CALL  DELETE( Y,IND,KND,KK,K02,MDEL,MJ4 )                       
      CALL  AICCOM( Y,N,K,KK,MJ4,OSD,OAIC )                             
C                                                                       
C         DECISION BY AIC                                               
C                                                                       
cc      IF(IPR.GE.3)  WRITE( 6,7 )   MDEL,OAIC,MDEL,J1,I1                 
      IF( OAIC .GT. AIC2 )  GO TO 180                                   
cc      IF( IPR .GE. 3 )     WRITE( 6,8 )                                 
      IF((IPR.GE.3).AND.(IFG.NE.0)) WRITE( LU,78 ) MDEL,OAIC,MDEL,J1,I1
      AIC2 = OAIC                                                       
      IF( AIC2 .GE. AIC1 )  GO TO 170                                   
      DO 165  I=1,KK                                                    
  165 JND(I) = KND(I)                                                   
      CALL  COPY( Y,KK1,0,0,MJ4,MJ1,X )                                 
  170 CONTINUE                                                          
      GO TO  200                                                        
  180 K02 = K02 + 1                                                     
cc      IF( IPR .GE. 3 )     WRITE( 6,18 )                                
      IF((IPR.GE.3).AND.(IFG.NE.0)) WRITE( LU,718 ) MDEL,OAIC,MDEL,J1,I1
  200 CONTINUE                                                          
      AICMIN = DMIN1( AIC1,AIC2 )                                       
C                                                                       
C         COMPARISON OF AIC1 AND AIC2                                   
C                                                                       
      K0 = K01                                                          
      IF( AIC1 .LE. AIC2 )  GO TO 220                                   
      K0 = K02                                                          
  220 CONTINUE                                                          
cc      IF(IPR.GE.3)  WRITE( 6,12 )     AICMIN ,K0                        
      IF( (IPR.GE.3).AND.(IFG.NE.0) )  WRITE( LU,12 )  AICMIN ,K0
  230 CONTINUE                                                          
C                                                                       
C                                                                       
C       SECOND ESTIMATE                                                 
C                                                                       
      IF( IPR .LT. 3 )  GO TO 240                                       
cc      WRITE( 6,14 )                                                     
cc      CALL  SRCOEF( X,K0,KK,N,MJ1,JND,IPR,A,SDD )                       
cc      WRITE( 6,644 )                                                    
      IF( IFG.NE.0 )  WRITE( LU,14 )
      CALL  SRCOEF( X,K0,KK,N,MJ1,JND,IPR,A,SDD,AAIC(II),IFG,LU )
      IF( IFG.NE.0 )  WRITE( LU,644 )                                   
  240 CONTINUE                                                          
C                                                                       
C--------------------   FINAL STEP OF AIC MINIMIZATION   ---------------
C                                                                       
      IF(KSW.NE.1)  GO TO 280                                           
C                                                                       
C          CHECK CONSTANT VECTOR  < DELETE ? >                          
C                                                                       
      CALL  COPY( X,KK1,0,0,MJ1,MJ4,Y )                                 
cc      IF(IPR.GE.3)  WRITE( 6,19 )                                       
cc      IF(IPR.GE.3)  WRITE( 6,13 )  AICMIN                               
      IF( (IPR.GE.3).AND.(IFG.NE.0) )  WRITE( LU,19 )
      IF( (IPR.GE.3).AND.(IFG.NE.0) )  WRITE( LU,13 )  AICMIN
      MDEL = 1                                                          
      K = K0-1                                                          
cc      CALL  DELETE( Y,D,IND,JND,KK,K0,MDEL,MJ4 )                        
      CALL  DELETE( Y,IND,JND,KK,K0,MDEL,MJ4 )                        
C                                                                       
      CALL  AICCOM( Y,N,K,KK,MJ4,OSP,OAIC )                             
cc      IF(IPR.GE.3)  WRITE( 6,21)  MDEL,OAIC,MDEL                        
      IF(OAIC.GE.AICMIN)  GO TO 250                                     
C                                                                       
cc      IF(IPR.GE.3)  WRITE( 6,8 )                                        
      IF( (IPR.GE.3).AND.(IFG.NE.0) )  WRITE( LU,218 )  MDEL,OAIC,MDEL
      AICMIN = OAIC                                                     
      MSW = 0                                                           
      CALL  COPY( Y,KK1,0,0,MJ4,MJ1,X )                                 
      GO TO 270                                                         
C                                                                       
  250 K0 = K0+1                                                         
cc      IF(IPR.GE.3)  WRITE( 6,18 )                                       
      IF( (IPR.GE.3).AND.(IFG.NE.0) )  WRITE( LU,2118 )  MDEL,OAIC,MDEL
      DO 260  I=1,KK                                                    
      J = IND(I)                                                        
  260 JND(J) = I                                                        
cc  270 IF(IPR.GE.3)  WRITE( 6,16 )  AICMIN,K0                            
  270 IF( (IPR.GE.3).AND.(IFG.NE.0) )  WRITE( LU,16 )  AICMIN,K0
  280 CONTINUE                                                          
C                                                                       
C         CHECK REGRESSOR MDEL  < DELETE ? >                            
C                                                                       
      DO 400     I1=1,ID                                                
C                                                                       
      CALL  COPY( X,KK1,0,0,MJ1,MJ4,Y )                                 
C                                                                       
cc      IF(IPR.GE.3)  WRITE( 6,11 )     I1                                
cc      IF(IPR.GE.3)  WRITE( 6,13 )     AICMIN                            
      IF( (IPR.GE.3).AND.(IFG.NE.0) )  WRITE( LU,11 )     I1
      IF( (IPR.GE.3).AND.(IFG.NE.0) )  WRITE( LU,13 )     AICMIN
      DO 330  J0=1,IMP1                                                 
      J1 = J0 - 1                                                       
      IF( J1 .EQ. 0 .AND. I1.GE. II )  GO TO 330                        
      MDEL = (J1-1)*ID + I1 + KSW                                       
      IF( J1 .EQ. 0 )     MDEL = MD + I1                                
      DO 310  I=1,MD2                                                   
      IF( JND(I) .EQ. MDEL )  GO TO 320                                 
  310 CONTINUE                                                          
      GO TO 330                                                         
  320 K = K0 - 1                                                        
      IF( I .GT. K0 )   GO TO 400                                       
cc      CALL  DELETE( Y,D,IND,JND,KK,K0,MDEL,MJ4 )                        
      CALL  DELETE( Y,IND,JND,KK,K0,MDEL,MJ4 )                        
C                                                                       
      CALL  AICCOM( Y,N,K,KK,MJ4,OSD,OAIC )                             
C                                                                       
cc      IF(IPR.GE.3)  WRITE( 6,7 )   MDEL,OAIC,MDEL,J1,I1                 
      IF( OAIC .GE. AICMIN )  GO TO 340                                 
C                                                                       
cc      IF( IPR .GE. 3 )     WRITE( 6,8 )                                 
      IF((IPR.GE.3).AND.(IFG.NE.0)) WRITE( LU,78 ) MDEL,OAIC,MDEL,J1,I1
      AICMIN = OAIC                                                     
      CALL  COPY( Y,KK1,0,0,MJ4,MJ1,X )                                 
C                                                                       
  330 CONTINUE                                                          
      GO TO 400                                                         
cc  340 K0 = K0 + 1                                                       
  340 CONTINUE
      K0 = K0 + 1                                                       
cc      IF( IPR .GE. 3 )     WRITE( 6,18 )                                
      IF((IPR.GE.3).AND.(IFG.NE.0)) WRITE( LU,718 ) MDEL,OAIC,MDEL,J1,I1
      DO 350  I=1,KK                                                    
      J = IND(I)                                                        
  350 JND(J) = I                                                        
  400 CONTINUE                                                          
cc      IF( IPR .GE. 3 )     WRITE( 6,16 )   AICMIN , K0                  
      IF( (IPR.GE.3).AND.(IFG.NE.0) )   WRITE( LU,16 )   AICMIN , K0
C                                                                       
C                                                                       
C                                                                       
cc      IF(IPR.GE.3)  WRITE( 6,15 )                                       
cc      CALL  SRCOEF( X,K0,KK,N,MJ1,JND,IPR,A,SDD )                       
      IF( (IPR.GE.3).AND.(IFG.NE.0) )  WRITE( LU,15 )
      CALL  SRCOEF( X,K0,KK,N,MJ1,JND,IPR,A,SDD,AAIC(II),IFG,LU )
      NPR(II) = K0
      DO 450  KK = 1,K0
         JNDF(KK,II) = JND(KK)
         AF(KK,II) = A(KK)
  450 CONTINUE
cc      EX(II) = SDD / DFLOAT(N)                                          
      EX(II) = SDD / DBLE(N)                                          
      CALL  COEF2( A,K0,ID,II,JND,LMAX,M,KSW,MSW,MJ2,B,C,E )            
      AICSUM = AICSUM + AICMIN                                          
      NSW = MAX0( NSW,MSW )                                             
C                                                                       
C                                                                       
  500 CONTINUE                                                          
C                                                                       
      IF( ISW .EQ. 1 )   RETURN                                         
cc      IF(IPR.GE.1)  WRITE( 6,3 )                                        
cc      CALL  MCOEF( B,C,E,EX,ID,LMAX,NSW,IPR,MJ2,MJ3 )                   
      CALL  MCOEF( BI,B,C,EI,E,EX,ID,LMAX,NSW,IPR,MJ2,MJ3 )
cc      IF(IPR.GE.1)  WRITE( 6,611 )  AICSUM                              
C                                                                       
      AICS = AICSUM                                                     
C                                                                       
C                                                                       
      RETURN                                                            
    3 FORMAT( 1H ,132(1H-) )                                            
cc    4 FORMAT( ' *****  INITIAL ESTIMATE *****' )                        
    4 FORMAT( /,' *****  INITIAL ESTIMATE *****' )
    5 FORMAT( ' ADD',4X ,I3,5X,'AIC =',F15.3,17X,'REGRESSOR',I4,
     *' ( LAG =',I2,' , I=',I2,' ) ' )
    6 FORMAT( 1H+,45X,'-----',36X,'ADDED    -----' )                    
   56 FORMAT( ' ADD',4X ,I3,5X,'AIC =',F15.3,10X,'-----  REGRESSOR',I4,
     *' ( LAG =',I2,' , I=',I2,' )  ADDED    -----' )
    7 FORMAT( ' DELETE ',I3,5X,'AIC =',F15.3,17X,'REGRESSOR',I4,
     *' ( LAG =',I2,' , I=',I2,' ) ' )
    8 FORMAT( 1H+,45X,'-----',36X,'DELETED  -----' )                    
   78 FORMAT( ' DELETE ',I3,5X,'AIC =',F15.3,10X,'-----  REGRESSOR',I4,
     *' ( LAG =',I2,' , I=',I2,' )  DELETED  -----' )
    9 FORMAT( 1H+,86X,'NOT ADDED' )                                     
   59 FORMAT( ' ADD',4X ,I3,5X,'AIC =',F15.3,17X,'REGRESSOR',I4,
     *' ( LAG =',I2,' , I=',I2,' )  NOT ADDED' )
cc   11 FORMAT( 1H ,10X,'--------  CHECK',I3,'-TH  VARIABLE  --------' )  
cc   12 FORMAT( 1H ,10X,'<<<  AIC =',F15.3,'  >>>   .....TEMPORARY MINIMUM
   11 FORMAT( /10X,' --------  CHECK',I3,'-TH  VARIABLE  --------' )  
   12 FORMAT( /10X,' <<<  AIC =',F15.3,'  >>>   .....TEMPORARY MINIMUM',
     1' AIC.....   ( NUMBER OF PARAMETERS=',I3,' )' )
cc   13 FORMAT( 1H ,15X,'AIC =',F15.3 )                                   
cc   14 FORMAT( 1H ,'*****  SECONDARY ESTIMATE  *****' )                  
cc   15 FORMAT( 1H ,'*****  FINAL ESTIMATE  *****' )                      
cc   16 FORMAT( 1H ,11X,'<<< MINIMUM AIC =',F15.3,' >>>    ( NUMBER OF PAR
   13 FORMAT( /15X,' AIC =',F15.3 )                                   
   14 FORMAT( /' *****  SECONDARY ESTIMATE  *****' )                  
   15 FORMAT( /' *****  FINAL ESTIMATE  *****' )                      
   16 FORMAT( /11X,' <<< MINIMUM AIC =',F15.3,' >>>    ( NUMBER OF ',
     *'PARAMETERS =',I3,' )' )
   17 FORMAT( 1H+,22X,'<<< MAICE >>>' )                                 
   18 FORMAT( 1H+,86X,'NOT DELETED' )                                   
  718 FORMAT( ' DELETE ',I3,5X,'AIC =',F15.3,17X,'REGRESSOR',I4,
     *' ( LAG =',I2,' , I=',I2,' )  NOT DELETED' )
cc   19 FORMAT( 1H ,10X,'--------  CHECK CONSTANT VECTOR  --------' )     
   19 FORMAT( /10X,' --------  CHECK CONSTANT VECTOR  --------' )     
   21 FORMAT( ' DELETE ',I3,5X,'AIC =',F15.3,17X,'REGRESSOR',I4,        
     * ' (CONSTANT  VECTOR)')                                           
  218 FORMAT( ' DELETE ',I3,5X,'AIC =',F15.3,10X,'-----  REGRESSOR',I4,
     * ' (CONSTANT  VECTOR)  DELETED  -----')
 2118 FORMAT( ' DELETE ',I3,5X,'AIC =',F15.3,10X,'-----  REGRESSOR',I4,
     * ' (CONSTANT  VECTOR)  DELETED  -----')
cc  600 FORMAT( 1H ,'-----  X  -----' )                                   
cc  602 FORMAT( 1H ,'AICMIN =',D13.5,5X,'AIC1 =',D13.5,5X,'AIC2 =',D13.5, 
  600 FORMAT( /' -----  X  -----' )                                   
  602 FORMAT( /' AICMIN =',D13.5,5X,'AIC1 =',D13.5,5X,'AIC2 =',D13.5, 
     1 5X,'K0 =',I5,5X,'K01 =',I5,5X,'K02 =',I5 )                       
cc  603 FORMAT( 1H ,'-----  Y  -----' )                                   
cc  605 FORMAT( 1H ,4X,'I',4X,'JND(I)',10X,'A(I)' )                       
cc  606 FORMAT( 1H ,'LMAX =',I5 )                                         
  603 FORMAT( /' -----  Y  -----' )                                   
  605 FORMAT( /5X,'I',4X,'JND(I)',10X,'A(I)' )                       
  606 FORMAT( /' LMAX =',I5 )                                         
  607 FORMAT( /I5,I10,D20.5 )                                        
cc  611 FORMAT( 1H ,'AIC =',F15.3 )                                       
cc  614 FORMAT( 1H ,'--  ADDVAR VARIABLES     MA =',I5,5X,'J1 =',I5,5X,   
  611 FORMAT( /' AIC =',F15.3 )                                       
  614 FORMAT( /' --  ADDVAR VARIABLES     MA =',I5,5X,'J1 =',I5,5X,   
     1  'I1 =',I5 )                                                     
cc  615 FORMAT( 1H ,'--  DELETE VARIABLES     MDEL =',I5,5X,'J1 =',I5,5X, 
  615 FORMAT( /' --  DELETE VARIABLES     MDEL =',I5,5X,'J1 =',I5,5X, 
     1  'I1 =',I5 )                                                     
  616 FORMAT( 1H ,I5,F15.3 )                                            
  630 FORMAT( 1H ,40I3 )                                                
  635 FORMAT( 1H ,'-----  JND  -----' )                                 
  637 FORMAT( 1H ,'-----  IND  -----' )                                 
  638 FORMAT( 1H ,'OAIC =',F15.3,5X,'AIC1 =',F15.3 )                    
  639 FORMAT( 1H ,'OAIC =',F15.3,5X,'AIC2 =',F15.3 )                    
  640 FORMAT( 1H ,4X,'I',10X,'AIC ' )                                   
  641 FORMAT( 1H ,'OAIC =',F15.3,5X,'AICMIN =',F15.3 )                  
cc  642 FORMAT( 1H ,30(1H-),/,1H ,'FIRST STEP OF AIC MINIMIZATION',/,1H , 
  642 FORMAT( /1H ,30(1H-),/,1H ,'FIRST STEP OF AIC MINIMIZATION',/,1H ,
     1 30(1H-) )                                                        
cc  643 FORMAT( 1H ,31(1H-),/,1H ,'SECOND STEP OF AIC MINIMIZATION',/,1H ,
  643 FORMAT(/1H ,31(1H-),/,1H ,'SECOND STEP OF AIC MINIMIZATION',/,1H ,
     1  31(1H-) )                                                       
cc  644 FORMAT( 1H ,30(1H-),/,1H ,'FINAL STEP OF AIC MINIMIZATION',/,1H , 
  644 FORMAT( /1H ,30(1H-),/,1H ,'FINAL STEP OF AIC MINIMIZATION',/,1H ,
     1 30(1H-) )                                                        
cc  645 FORMAT( 1H ,23X,10(1H.),'REGRESSION MODEL FOR THE REGRESSAND  II =
  645 FORMAT( /24X,10(1H.),'REGRESSION MODEL FOR THE REGRESSAND  II =
     1',I2,2X,10(1H.) )                                                 
      END                                                               
C
C
cc      SUBROUTINE  MBYSAR( X,D,N,M,ID,KSW,IPR,MJ1,MJ2,A,B,G,H,E,AICB,EK )
      SUBROUTINE  MBYSAR( X,N,M,ID,KSW,IPR,MJ1,MJ2,SD1,AIC1,DIC1,
     * AICM1,SDMIN1,IMIN1,C,D,A,B,G,H,E,AICB,EK )
C                                                                       
C     THIS SUBROUTINE PRODUCES MULTI-VARIATE AUTOREGRESSIVE MODELS BY A 
C     BAYESIAN PROCEDURE USING THE OUTPUT OF SUBROUTINE MREDCT.         
C       ----------------------------------------------------------------
C       THE FOLLOWING SUBROUTINES ARE DIRECTLY CALLED BY THIS SUBROUTINE
C             COPY                                                      
C             MAICE                                                     
C             BAYSWT                                                    
C             HUSHLD                                                    
C             HUSHL1                                                    
C             MARCOF                                                    
C             MBYSPC                                                    
C             MPARCO                                                    
C             MSDCOM                                                    
C             PRINT3                                                    
C       ----------------------------------------------------------------
C                                                                       
C       INPUTS:                                                         
C         X:      ((M+1)*ID+KSW)*((M+1)*ID+KSW) TRIANGULAR MATRIX, OUTPU
C                 SUBROUTINE MREDCT                                     
C         D:      WORKING AREA                                          
C         N:      DATA LENGTH                                           
C         M:      MAXIMUM TIME LAG OF THE MODEL                         
C         ID:     DIMENSION OF THE OBSERVATION                          
C         KSW:    =0   CONSTANT TERM IS NOT INCLUDED AS A REGRESSOR     
C                 =1   CONSTANT TERM IS INCLUDED AS THE FIRST REGRESSOR 
C         IPR:    PRINT OUT CONTROL                                     
C         MJ1:    ABSOLUTE DIMENSION OF X IN THE MAIN PROGRAM           
C         MJ2:    ABSOLUTE DIMENSION OF A,B,G,H AND E IN THE MAIN PROGRA
C                                                                       
C       OUTPUTS:                                                        
C         A:      AR-COEFFICIENT MATRICES OF FORWARD MODEL              
C         B:      AR-COEFFICIENT MATRICES OF BACKWARD MODEL             
C         G:      PARTIAL AUTOCORRELATION COEFFICIENT MATRICES OF FORWAR
C         H:      PARTIAL AUTOCORRELATION COEFFICIENT MATRICES OF BACKWA
C         E:      INNOVATION VARIANCE-COVARIANCE MATRIX                 
C         AICB:   EQUIVALENT AIC OF THE BAYESIAN (FORWARD) MODEL        
C         EK:     EQUIVALENT NUMBER OF AUTOREGRESSIVE COEFFICIENTS      
C                                                                       
      IMPLICIT REAL * 8  ( A-H , O-Z )                                  
cc      DIMENSION  X(MJ1,1) , D(1) , A(MJ2,MJ2,1) , B(MJ2,MJ2,1)          
      DIMENSION  X(MJ1,1), D(M), A(MJ2,MJ2,1), B(MJ2,MJ2,1)          
      DIMENSION  G(MJ2,MJ2,1) , H(MJ2,MJ2,1)                            
      DIMENSION  E(MJ2,1)                                               
cc      DIMENSION  AIC(51) , SD(51)                                       
      DIMENSION  AIC1(M+1) , SD1(M+1), DIC1(M+1), DIC2(M+1)
      DIMENSION  AIC(M+1) , SD(M+1)
cc      DIMENSION  IND(100) , JND(100) , C(100)                           
cc      DIMENSION  Y(100,10)                                              
      DIMENSION  IND((M+1)*ID+KSW) , JND((M+1)*ID+KSW)
      DIMENSION  C(M+1)
      DIMENSION  Y((M+1)*ID+KSW,ID)                                     
C
C                                                                       
C          ---------------                                              
C          INITIAL SETTING                                              
C          ---------------                                              
C                                                                       
cc      MJ4 = 100                                                         
      M1 = M + 1                                                        
      MD = M * ID + KSW                                                 
      MD2= M1*ID + KSW                                                  
	MJ4 = MD2
C                                                                       
      DO 30  I=1,M1                                                     
      SD(I)  = 1.D0                                                     
   30 AIC(I) = 0.D0                                                     
      CALL  COPY( X,MD2,0,MD2,MJ1,MJ1,X )                               
C                                                                       
C                                                                       
cc      IF(IPR.GE.3)  WRITE( 6,3 )                                        
C                                                                       
C          ------------------------------------------------------------ 
C          PARTIAL AUTOREGRESSION COEFFICIENT-MATRICES OF FORWARD MODEL 
C          ------------------------------------------------------------ 
      CALL  MPARCO( X,ID,M,KSW,0,MJ1,MJ2,G,H )                          
cc      IF( IPR .GE. 3 )     WRITE( 6,4 )                                 
cc      IF( IPR .GE. 3 )     CALL PRINT3( G,ID,ID,M,MJ2,MJ2,1 )           
C                                                                       
C          -----------------------------------------------              
C          AIC COMPUTATION ( FORWARD AND BACKWARD MODELS )              
C          -----------------------------------------------              
      J0 = MD                                                           
cc      ASSIGN 70 TO ISUB                                                 
      ISUB =  70
C                                                                       
   40 DO 60  II=1,M1                                                    
      K = (II-1)*ID + KSW                                               
      K1= K+1                                                           
      DO 50  J=1,ID                                                     
      JJ = J + J0                                                       
      DO 50  I=K1,MD2                                                   
      I1 = I - K                                                        
   50 Y(I1,J) = X(I,JJ)                                                 
C                                                                       
      MDMK = MD2 - K                                                    
cc      CALL  HUSHLD( Y,D,MJ4,MDMK,ID )                                   
      CALL  HUSHLD( Y,MJ4,MDMK,ID )                                   
C                                                                       
      DO 60  I=1,ID                                                     
      YY = Y(I,I)**2 / N                                                
      SD(II)  = SD(II) * YY                                             
   60 AIC(II) = AIC(II) + N*DLOG( YY ) + 2.D0*(K+1)                     
cc      IF( IPR .GE. 2 )     WRITE( 6,615 )                               
cc      GO TO ISUB, ( 70,100 )                                            
      IF( ISUB .EQ. 70 ) GO TO 70
      IF( ISUB .EQ. 100 ) GO TO 100
C                                                                       
cc   70 CALL  MAICE( AIC,SD,M,IPR-1,AICM,SDMIN,IMIN )                     
   70 CALL  MAICE( AIC,SD,M,IPR-1,AICM,SDMIN,IMIN,DIC1 )
      DO 71  I=1,M1
         AIC1(I) = AIC(I)
         SD1(I) = SD(I)
   71 CONTINUE
      AICM1 = AICM
      SDMIN1 = SDMIN
      IMIN1 = IMIN
C          -------------------------------------------------------------
C          PARTIAL AUTOREGRESSION COEFFICIENT-MATRICES OF BACKWARD MODEL
C          -------------------------------------------------------------
      DO 80  I=1,MD2                                                    
   80 IND(I) = I                                                        
      JND(1) = 1                                                        
      J2 = KSW                                                          
      DO 85  JJ=2,M                                                     
      J = (M-JJ)*ID                                                     
      DO 85  I=1,ID                                                     
      J2 = J2 + 1                                                       
      J1 = J +I                                                         
   85 JND(J2) = J1                                                      
      DO 90  I=1,ID                                                     
      J2 = J2 + 1                                                       
      J1 = MD + I                                                       
      JND(J2) = J1                                                      
      J3 = J2 + ID                                                      
      J1 = J1 - ID                                                      
   90 JND(J3) = J1                                                      
C                                                                       
cc      CALL  HUSHL1( X,D,MJ1,MD2,MD2,1,IND,JND )                         
      CALL  HUSHL1( X,MJ1,MD2,MD2,1,IND,JND )                         
      CALL  MPARCO( X,ID,M,KSW,1,MJ1,MJ2,G,H )                          
cc      IF( IPR .GE. 3 )      WRITE( 6,5 )                                
cc      IF( IPR .GE. 3 )      WRITE( 6,9 )                                
cc      IF( IPR .GE. 3 )      CALL PRINT3( H,ID,ID,M,MJ2,MJ2,1 )          
C                                                                       
C          ---------------                                              
C          AIC COMPUTATION                                              
C          ---------------                                              
      J0 = MD - ID                                                      
cc      ASSIGN 100 TO ISUB                                                
      ISUB = 100
      GO TO 40                                                          
C                                                                       
  100 DO 110  I=1,M1                                                    
      SD(I) = DSQRT( SD(I) )                                            
  110 AIC(I) = AIC(I) / 2                                               
C                                                                       
cc      IF( IPR .GE. 3 )     WRITE( 6,616 )                               
cc      CALL  MAICE( AIC,SD,M,IPR-2,AICM,SDMIN,IMIN )                     
      CALL  MAICE( AIC,SD,M,IPR-2,AICM,SDMIN,IMIN,DIC2 )
C                                                                       
C          ----------------------------                                 
C          BAYESIAN WEIGHTS COMPUTATION                                 
C          ----------------------------                                 
C                                                                       
      CALL  BAYSWT( AIC,AICM,M,0,C )                                    
C                                                                       
C          -------------------------------------------------            
C          PARTIAL AUTOREGRESSION MATRICES OF BAYESIAN MODEL            
C          -------------------------------------------------            
C                                                                       
      CALL  MBYSPC( G,H,C,D,M,ID,IPR,MJ2 )                              
cc      IF( IPR .GE. 1 )     WRITE( 6,11 )                                
cc      IF( IPR .LT. 2 )     GO TO 115                                    
cc      WRITE( 6,4 )                                                      
cc      CALL  PRINT3( G,ID,ID,M,MJ2,MJ2,1 )                               
cc      WRITE( 6,5 )                                                      
cc      CALL  PRINT3( H,ID,ID,M,MJ2,MJ2,1 )                               
C                                                                       
C          -----------------------------------------------------        
C          AUTOREGRESSION COEFFICIENT MATRICES OF BAYESIAN MODEL        
C          -----------------------------------------------------        
C                                                                       
  115 CALL  MARCOF( G,H,ID,M,MJ2,A,B )                                  
cc      IF( IPR .GE. 1 )   WRITE( 6,6 )                                   
cc      IF( IPR .GE. 1 )   CALL  PRINT3( A,ID,ID,M,MJ2,MJ2,1 )            
cc      IF( IPR .GE. 3 )   WRITE( 6,7 )                                   
cc      IF( IPR .GE. 3 )   CALL  PRINT3( B,ID,ID,M,MJ2,MJ2,1 )            
C          -------------------------                                    
C          AIC OF THE BAYESIAN MODEL                                    
C          -------------------------                                    
      EK = 0.D0                                                         
      DO 120  I=1,M                                                     
cc  120 EK = EK + C(I)**2                                                 
  120 EK = EK + D(I)**2                                                 
      EK = EK*(ID**2)                                                   
      CALL  COPY( X,MD2,MD2,0,MJ1,MJ1,X )                               
cc      CALL  MSDCOM( X,A,Y,D,N,M,ID,KSW,IPR,MJ1,MJ2,MJ4,E,OSD )          
      CALL  MSDCOM( X,A,N,M,ID,KSW,IPR,MJ1,E,OSD )          
C                                                                       
      AICB = N*DLOG( OSD ) + 2.D0*EK + 2.D0*KSW*ID + ID*(ID+1)          
cc      IF( IPR .GE. 1 )     WRITE( 6,614 )   AICB                        
C                                                                       
      RETURN                                                            
    3 FORMAT( 1H ,132(1H-) )                                            
    4 FORMAT( /,1H ,'PARTIAL AUTOREGRESSION COEFFICIENTS  ( FORWARD MODEL
     1)' )                                                              
    5 FORMAT( /,1H ,'PARTIAL AUTOREGRESSION COEFFICIENTS ( BACKWARD MODEL
     1)' )                                                              
    6 FORMAT( /,1H ,'AR-COEFFICIENT MATRICES  ( FORWARD MODEL )' )      
    7 FORMAT( /,1H ,'AR-COEFFICIENT MATRICES  ( BACKWARD MODEL )' )     
    8 FORMAT( 1H ,10X,'SD  = RESIDUAL VARIANCE',15X,'=',D23.12,/,11X,'EK
     1  = EQUIVALENT NUMBER OF PARAMETERS =',F15.3,/,11X,'AIC =  N*LAG(S
     2D) + 2*EK',15X,'=',F15.3 )                                        
    9 FORMAT( /,1H ,23(1H-),/,' LEAST SQUARES ESTIMATES',/,1H ,23(1H-) )
   11 FORMAT( /,1H ,18(1H-),/,' BAYESIAN ESTIMATES',/,1H ,18(1H-) )     
  610 FORMAT( 1H ,10D13.5 )                                             
  612 FORMAT( 1H ,'-----  LEAST SQUARE ESTIMATES OF FULL ORDER MODEL  --
     1---' )                                                            
  613 FORMAT( 1H ,'-----  BAYESIAN ESTIMATES OF REGRESSION COEFFICIENTS 
     * -----' )                                                         
  614 FORMAT( 1H ,10X,'AIC = ',F15.3)                                   
  615 FORMAT( /,1H ,'-----  FORWARD MODELS  -----' )                    
  616 FORMAT( /,1H ,'-----  AVERAGE OF FORWARD AND BACKWARD MODELS  
     :-----' )                                                          
  642 FORMAT( 1H ,30(1H-),/,1H ,'FIRST STEP OF AIC MINIMIZATION',/,1H , 
     1 30(1H-) )                                                        
  645 FORMAT( 1H ,23X,10(1H.),2X,'REGRESSION MODEL FOR THE REGRESSAND  I
     1I =',I2,2X,10(1H.) )                                              
      END                                                               
C
C
      SUBROUTINE  MBYSPC( G,H,C,D,M,ID,IPR,MJ2 )                        
C                                                                       
C     THIS SUBROUTINE PRODUCES PARTIAL AUTOREGRESSION COEFFICIENTS G(I),
C     (I=1,K) OF THE MULTI-VARIATE AUTOREGRESSIVE MODEL.                
C                                                                       
C       INPUTS:                                                         
C          G:      LEAST SQUARES ESTIMATES OF "PARCOR'S" (FORWARD MODEL)
C          H:      LEAST SQUARES ESTIMATES OF "PARCOR'S" (BACKWARD MODEL
C          C(I+1): BAYESIAN WEIGHT OF EACH ORDER  (I=0,...,M)           
C          M:      MAXIMUM TIME LAG OF THE MODEL                        
C          ID:     DIMENSION OF THE PROCESS                             
C          IPR:    PRINT OUT CONTROL                                    
C          MJ2:    ABSOLUTE DIMENSION OF G AND H                        
C                                                                       
C       OUTPUTS:                                                        
C          G:      BAYESIAN ESTIMATES OF "PARCOR'S" (FORWARD MODEL)     
C          H:      BAYESIAN ESTIMATES OF "PARCOR'S" (BACKWARD MODEL)    
C          D:      INTEGRATED BAYESIAN WEIGHT EACH ORDER (I=1,...,M)    
C                                                                       
      IMPLICIT  REAL * 8 ( A-H , O-Z )                                  
      DIMENSION  G(MJ2,MJ2,1) , H(MJ2,MJ2,1) , C(1) , D(1)              
C                                                                       
C          INTEGRATED BAYESIAN WEIGHT                                   
C                                                                       
      M1 = M+1                                                          
      D(M) = C(M1)                                                      
      DO 10 I=2,M                                                       
      J = M1 - I                                                        
   10 D(J) = D(J+1) + C(J+1)                                            
cc      IF( IPR .LE. 1 )   GO TO 20                                       
cc      WRITE( 6,6 )                                                      
cc      DO 15  I=1,M1                                                     
cc       IM1 = I - 1                                                      
cc      IF( I .EQ. 1 )   WRITE( 6,7 )   IM1 , C(I)                        
cc      IF( I .NE. 1 )   WRITE( 6,7 )   IM1 , C(I) , D(IM1)               
cc   15 CONTINUE                                                          
C                                                                       
cc   20 DO 30  I=1,M                                                      
cc   30 C(I) = D(I)                                                       
C                                                                       
C          PARTIAL CORRELATION                                          
C                                                                       
      DO 40  II=1,M                                                     
      DO 40  J=1,ID                                                     
      DO 40  I=1,ID                                                     
      G(I,J,II) = G(I,J,II)*D(II)                                       
   40 H(I,J,II) = H(I,J,II)*D(II)                                       
C                                                                       
C                                                                       
      RETURN                                                            
    6 FORMAT( 1H ,4X,'M',14X,'BAYESIAN WEIGHTS',8X,'INTEGRATES BAYESIAN 
     1WEIGHTS' )                                                        
    7 FORMAT( 1H ,I5,2D30.7 )                                           
      END                                                               
C
C
cc      SUBROUTINE  MCOEF( B,C,E,EX,ID,LMAX,KSW,IPR,MJ2,MJ3 )             
      SUBROUTINE  MCOEF( BI,B,C,EI,E,EX,ID,LMAX,KSW,IPR,MJ2,MJ3 )
C                                                                       
C     THIS SUBROUTINE COMPUTES AND PRINTS OUT THE COEFFICIENT MATRICES O
C     MULTI-VARIATE AUTOREGRESSIVE MODEL FROM THE COEFFICIENT MATRICES O
C     THE MODEL WITH INSTANTANEOUS RESPONSE                             
C       ----------------------------------------------------------------
C       THE FOLLOWING SUBROUTINES ARE DIRECTLY CALLED BY THIS SUBROUTINE
C             INVDET                                                    
C             PRINT3                                                    
C             TRIINV                                                    
C       ----------------------------------------------------------------
C                                                                       
C       INPUTS:                                                         
C          B:     REGRESSION COEFFICIENTS MATRIX                        
C          C:     CONSTANT VECTOR                                       
C          E:     COEFFICIENT OF QUICK RESPONSE                         
C          EX:    RESIDUAL VARIANCES OF ORTHOGONALIZED MODEL            
C          ID:    DIMENSION OF THE OBSERVATION                          
C          LMAX:  HIGHEST ORDER OF THE AR-MODEL                         
C          KSW:   =0  CONSTANT VECTOR IS NOT INCLUDED AS A REGRESSOR    
C                 =1  CONSTANT VECTOR IS INCLUDED AS THE FIRST REGRESSOR
C          IPR:   PRINT OUT CONTROL                                     
C          MJ2:   ABSOLUTE DIMENSION OF B AND E                         
C          MJ3:   ABSOLUTE DIMENSION OF B                               
C       OUTPUTS:                                                        
C          B:     AR-COEFFICIENT MATRIX                                 
C          C:     CONSTANT VECTOR                                       
C          E:     INNOVATION COVARIANCE MATRIX                          
C                                                                       
      IMPLICIT  REAL*8 ( A-H,O-Z )                                      
      DIMENSION  B( MJ2,MJ2,MJ3 ) , E( MJ2,1 ) , EX( 1 ) , C( 1 )       
cc      DIMENSION  C1(24) , EE(24,24)                                     
      DIMENSION  C1(ID) , EE(ID,ID)
C
C    INPUT  E ---> EI
C    INPUT  B ---> BI
      DIMENSION  EI(ID,ID), BI(ID,ID,LMAX)
C
cc      MJ5 = 24                                                          
      MJ5 = ID
C                                                                       
      IF(IPR.LE.1)  GO TO 300                                           
cc      WRITE( 6,603 )                                                    
      DO 310  I=1,ID                                                    
      DO 310  J=1,ID
  310 EI(I,J) = E(I,J)
cc  310 WRITE( 6,610 )  (E(I,J),J=1,ID)                                   
cc      WRITE( 6,604 )                                                    
      DO  330     II=1,LMAX                                             
      DO  320     I=1,ID                                                
      DO  320     J=1,ID
  320 BI(I,J,II) = B(I,J,II)
cc  320 WRITE( 6,610 )     (B(I,J,II),J=1,ID)                             
cc  330 WRITE( 6,609 )                                                    
  330 CONTINUE
      IF( KSW .EQ. 0 )     GO TO 300                                    
cc      WRITE( 6,601 )                                                    
cc      WRITE( 6,610 )     (C(I),I=1,ID)                                  
  300 CONTINUE                                                          
C                                                                       
      IF( KSW .NE. 1 )     GO TO 375                                    
      DO  335     I=1,ID                                                
      DO  335     J=1,ID                                                
      SUM = E(I,J)                                                      
      DO  325     II=1,LMAX                                             
  325 SUM = SUM - B(I,J,II)                                             
  335 EE(I,J) = SUM                                                     
      CALL  INVDET( EE,EDET,ID,MJ5 )                                    
      DO  355     I=1,ID                                                
      SUM = 0.D0                                                        
      DO  345     J=1,ID                                                
  345 SUM = SUM + EE(I,J)*C(J)                                          
  355 C1(I) = SUM                                                       
      DO  365     I=1,ID                                                
  365 C(I) = C1(I)                                                      
C                                                                       
  375 CONTINUE                                                          
      CALL  TRIINV( E,ID,MJ2,MJ5,EE )                                   
C                                                                       
      DO  360     II=1,LMAX                                             
      DO  350     I=1,ID                                                
      DO  350     J=1,ID                                                
      SUM = 0.D0                                                        
      DO  340     JJ=1,I                                                
  340 SUM = SUM + EE(I,JJ)*B(JJ,J,II)                                   
  350 E(I,J) = SUM                                                      
      DO  360     I=1,ID                                                
      DO  360     J=1,ID                                                
  360 B(I,J,II) = E(I,J)                                                
      DO  370     I=1,ID                                                
      SUM = 0.D0                                                        
      DO  371     J=1,ID                                                
  371 SUM = SUM + EE(I,J)*C(J)                                          
  370 C1(I) = SUM                                                       
      DO  372     I=1,ID                                                
  372 C(I) = C1(I)                                                      
      DO  380     I=1,ID                                                
      DO  380     J=1,I                                                 
      SUM = 0.D0                                                        
      DO  385     II=1,J                                                
  385 SUM = SUM + EE(I,II) * EE(J,II)*EX(II)                            
      E(I,J) = SUM                                                      
  380 E(J,I) = SUM                                                      
cc      IF(IPR.EQ.0)   RETURN                                             
cc      WRITE( 6,612 )                                                    
cc      DO  390     I=1,ID                                                
cc  390 WRITE( 6,611 )     (E(I,J),J=1,ID)                                
cc      WRITE( 6,613 )                                                    
cc      CALL  PRINT3( B,ID,ID,LMAX,MJ2,MJ2,0 )                            
cc      IF(KSW.NE.1)  RETURN                                              
cc      WRITE( 6,602 )                                                    
cc      WRITE( 6,610 )     (C(I),I=1,ID)                                  
C                                                                       
      RETURN                                                            
  601 FORMAT( 1H ,'-----  C  -----' )                                   
  602 FORMAT( 1H ,'-----  CONSTANT VECTOR  -----' )                     
  603 FORMAT( 1H ,'-----  INSTANTANEOUS RESPONSE  -----' )              
  604 FORMAT( 1H ,'-----  B  -----' )                                   
  609 FORMAT( 1H  )                                                     
  610 FORMAT( 1H ,8F15.8 )                                              
  611 FORMAT( 1H ,8D15.7 )                                              
  612 FORMAT( 1H ,'*  INNOVATION VARIANCE MATRIX  *' )                  
  613 FORMAT( 1H ,'*  AR-COEFFICIENTS  *' )                             
      END                                                               
C
C
      SUBROUTINE  MPARCO( X,ID,M,KSW,IFG,MJ1,MJ3,G,H )                  
C                                                                       
C     THIS SUBROUTINE PRODUCES LEAST SQUARES ESTIMATES OF PARTIAL       
C     AUTOREGRESSION COEFFICIENT MATRICES OF FORWARD AND BACKWARD       
C     MULTI-DIMENSIONAL AUTOREGRESSIVE MODEL                            
C       ----------------------------------------------------------------
C       THE FOLLOWING SUBROUTINE IS DIRECTLY CALLED BY THIS SUBROUTINE: 
C             SOLVE                                                     
C       ----------------------------------------------------------------
C                                                                       
C       INPUTS:                                                         
C          X:      ((M+1)*ID+KSW)*((M+1)*ID+KSW) MATRIX                 
C          ID:     DIMENSION OF OBSERVATION                             
C          M:      MAXIMUM TIME LAG OF THE MODEL                        
C          KSW:    =0  CONSTANT TERM IS NOT INCLUDED AS A REGRESSOR     
C                  =1  CONSTANT TERM IS INCLUDED AS THE FIRST REGRESSOR 
C          IFG:    =0   TO REQUEST THE COEFFICIENTS OF FORWARD MODEL    
C                  =1   TO REQUEST THE COEFFICIENTS OF BACKWARD MODEL   
C          MJ1:    ABSOLUTE DIMENSION OF X                              
C          MJ3:    ABSOLUTE DIMENSION OF G AND H                        
C                                                                       
C       OUTPUTS:                                                        
C          G:      AR-COEFFICIENT MATRICES OF FORWARD MODEL             
C          H:      AR-COEFFICIENT MATRICES OF BACKWARD MODEL            
C                                                                       
      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
      DIMENSION  X(MJ1,1) , G(MJ3,MJ3,1) , H(MJ3,MJ3,1)                 
cc      DIMENSION  C(10,10) , R(10,10)                                    
      DIMENSION  C(ID,ID) , R(ID,ID)                                    
C                                                                       
cc      MJ2 = 10                                                          
      MJ2 = ID
      IF( IFG .NE. 0 )     GO TO 30                                     
      MD = ID*M + KSW                                                   
      DO  20     II=1,M                                                 
      I0 = (II-1)*ID + KSW                                              
      DO  10     J=1,ID                                                 
      J1 = J + I0                                                       
      J2 = J + MD                                                       
      DO  10     I=1,ID                                                 
      I1 = I + I0                                                       
      C(I,J) = X(I1,J1)                                                 
   10 R(I,J) = X(I1,J2)                                                 
   20 CALL  SOLVE( C,R,ID,II,MJ2,MJ3,G )                                
      GO TO 60                                                          
C                                                                       
   30 MM1 = M - 1                                                       
      JS0 = MM1*ID + KSW                                                
      DO  50     II=1,M                                                 
      I0 = (II-1)*ID + KSW                                              
      J0 = (MM1-II)*ID + KSW                                            
      IF( II .EQ. M )     J0 = M*ID + KSW                               
      DO  40     J=1,ID                                                 
      J1 = J0 + J                                                       
      J2 = JS0 + J                                                      
      DO  40     I=1,ID                                                 
      I1 = I0 + I                                                       
      C(I,J) = X(I1,J1)                                                 
   40 R(I,J) = X(I1,J2)                                                 
   50 CALL  SOLVE( C,R,ID,II,MJ2,MJ3,H )                                
C                                                                       
   60 RETURN                                                            
      END                                                               
C
C
cc      SUBROUTINE  MRDATA( MT,MJ,Z,N,ID )                                
      SUBROUTINE MRDATA( ZS,Z,N,ID,C,ZMEAN,ZVARI )
C                                                                       
C         +-----------------------------------------+                   
C         ! ORIGINAL DATA LOADING AND MEAN DELETION !                   
C         +-----------------------------------------+                   
C                                                                       
C     THIS SUBROUTINE IS USED FOR THE LOADING OF ORIGINAL DATA AND      
C     DELETION OF THE MEAN VALUES.  THE DATA IS LOADED THROUGH THE DEVIC
C     SPECIFIED BY MT.  EACH DATA SET IS COMPOSED OF TITLE, DATA LENGTH,
C     CALIBRATIONS OF DATA, DATA FORMAT AND THE ORIGINAL DATA.          
C                                                                       
C       INPUTS:                                                         
C          MT:      INPUT DEVICE SPECIFICATION                          
C          MJ:      ABSOLUTE DIMENSION OF Z                             
C                                                                       
C          TITLE:   SPECIFICATION OF DATA                               
C          N:       DATA LENGTH                                         
C          ID:      DIMENSION OF OBSERVATION                            
C          IFM:     CONTROL FOR INPUT                                   
C          FORM:    INPUT DATA FORMAT SPECIFICATION                     
C          C(I):    CALIBRATION OF CHANNEL I (I=1,ID)                   
C          Z:       ORIGINAL DATA; Z(K,I) (K=1,N) REPRESENTS THE I-TH CH
C                   RECORD                                              
C                                                                       
C       OUTPUTS:                                                        
C          Z:       ORIGINAL DATA ( MEAN DELETED )                      
C          N:       DATA LENGTH                                         
C          ID:      DIMENSION OF OBSERVATION                            
C                                                                       
      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
CC      REAL * 4  Z(MJ,1)                                                 
cc      DIMENSION  Z(MJ,1)
cc      DIMENSION  C(24)                                                  
      DIMENSION  ZS(N,ID), Z(N,ID), C(ID)
      DIMENSION  ZMEAN(ID), ZVARI(ID)
cc      REAL * 4  FORM(20) , TITLE(20)                                    
C                                                                       
C                                                                       
cc      READ( MT,2 )     TITLE                                            
cc      READ( MT,1 )     N , ID , IFM                                     
cc      READ( MT,2 )     FORM                                             
cc      READ( MT,3 )     (C(I),I=1,ID)                                    
C                                                                       
cc      GO TO ( 10,20,30,40,50,60 ) , IFM                                 
C                                                                       
cc   10 DO  15    I=1,N                                                   
cc   15 READ( MT,FORM )     (Z(I,J),J=1,ID)                               
cc      GO TO 70                                                          
cc   20 DO  25     J=1,ID                                                 
cc   25 READ( MT,FORM )     (Z(I,J),I=1,N)                                
cc      GO TO 70                                                          
cc   30 READ( MT,FORM )     ((Z(I,J),J=1,ID),I=1,N)                       
cc      GO TO 70                                                          
cc   40 READ( MT,FORM )     ((Z(I,J),I=1,N),J=1,ID)                       
cc      GO TO 70                                                          
cc   50 READ( MT )     ((Z(I,J),J=1,ID),I=1,N)                            
cc      GO TO 70                                                          
cc   60 READ( MT )     ((Z(I,J),I=1,N),J=1,ID)                            
cc   70 CONTINUE                                                          
C                                                                       
C                                                                       
C       ORIGINAL DATA PRINT OUT                                         
C                                                                       
cc      WRITE( 6,601 )     TITLE                                          
cc      WRITE( 6,602 )     N , ID , FORM                                  
cc      WRITE( 6,603 )                                                    
      DO  75     J=1,ID                                                 
cc      WRITE( 6,604 )     J                                              
cc   75 WRITE( 6,610 )     (Z(I,J),I=1,N)                                 
      DO  75     I=1,N
   75 Z(I,J) = ZS(I,J)
C                                                                       
      DO  85     J=1,ID                                                 
      CC = C(J)                                                         
      DO  85     I=1,N                                                  
   85 Z(I,J) = Z(I,J) * CC                                              
C                                                                       
C                                                                       
cc      WRITE( 6,606 )                                                    
      DO  120     J=1,ID                                                
C                                                                       
C     MEAN DELETION   AND  VARIANCE COMPUTATION                         
C                                                                       
      SUM = 0.D0                                                        
      DO  90     I=1,N                                                  
   90 SUM = SUM + Z(I,J)                                                
cc      ZMEAN = SUM / DFLOAT(N)                                           
      ZMEAN(J) = SUM / DFLOAT(N)
      DO  100     I=1,N                                                 
cc  100 Z(I,J) = Z(I,J) - ZMEAN                                           
  100 Z(I,J) = Z(I,J) - ZMEAN(J)
C                                                                       
      SUM = 0.D0                                                        
      DO  110     I=1,N                                                 
  110 SUM = SUM + Z(I,J)*Z(I,J)                                         
cc      ZVARI = SUM / DFLOAT(N)                                           
      ZVARI(J) = SUM / DFLOAT(N)
C                                                                       
C         MEAN AND VARIANCE PRINT OUT                                   
C                                                                       
cc  120 WRITE( 6,605 )     J , ZMEAN , ZVARI                              
  120 CONTINUE
C                                                                       
C                                                                       
      RETURN                                                            
C                                                                       
    1 FORMAT( 16I5 )                                                    
    2 FORMAT( 20A4 )                                                    
    3 FORMAT( 8F10.0 )                                                  
  601 FORMAT( 1H ,'TITLE  ---  ',20A4 )                                 
  602 FORMAT( 1H ,'N =',I5,5X,'ID =',I5,5X,'FORMAT =',20A4 )            
  603 FORMAT( 1H ,'** ORIGINAL DATA **' )                               
  604 FORMAT( 1H ,I5,'-CHANNEL' )                                       
  605 FORMAT( 1H ,I6,2D15.8 )                                           
  606 FORMAT( 1H ,5X,'I',10X,'MEAN',7X,'VARIANCE' )                     
  610 FORMAT( 1H ,10D13.5 )                                             
      END                                                               
C
C
cc      SUBROUTINE  MREDCT( Z,D,NMK,N0,LAG,ID,MJ,MJ1,KSW,X )              
      SUBROUTINE  MREDCT( Z,NMK,N0,LAG,ID,MJ,MJ1,KSW,X )              
C                                                                       
C         +-------------------------+                                   
C         ! HOUSEHOLDER'S REDUCTION !                                   
C         +-------------------------+                                   
C                                                                       
C     THIS SUBROUTINE FIRST SETS UP DATA MATRIX X BY AUGMENTING         
C     SUCCESSIVELY SHIFTED ORIGINAL DATA MATRIX Z AND THEN TRANSFORMS X 
C     INTO TRIANGULAR FORM BY HOUSEHOLDER TRANSFORMATION.               
C                                                                       
C       ----------------------------------------------------------------
C       THE FOLLOWING SUBROUTINES ARE DIRECTLY CALLED BY THIS SUBROUTINE
C             HUSHLD                                                    
C             MSETX1                                                    
C       ----------------------------------------------------------------
C                                                                       
C       INPUTS:                                                         
C          Z:      ORIGINAL DATA MATRIX                                 
C          D:      WORKING AREA                                         
C          NMK:    DIMENSION OF THE VECTOR OF REGRESSAND (Z(N0+K-KSW+1),
C                  ...,Z(N0+K-KSW+NMK) (J=1,ID)                         
C          N0:     INDEX OF THE END POINT OF DISCARDED FORMER OBSERVATIO
C          LAG:    HIGHEST TIME LAG OF THE MODEL                        
C          ID:     DIMENSION OF OBSERVATIONS                            
C          MJ:     ABSOLUTE DIMENSION OF Z                              
C          MJ1:    ABSOLUTE DIMENSION OF X                              
C          KSW:    =0   CONSTANT VECTOR IS NOT INCLUDED AS A REGRESSOR  
C                  =1   CONSTANT VECTOR IS INCLUDED AS THE FIRST REGRESS
C                                                                       
C         OUTPUT:                                                       
C          X:      REDUCED MATRIX ( UPPER TRIANGULAR FORM )             
C                                                                       
      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
cc      DIMENSION  X(MJ1,1) , D(1)                                        
CC      REAL * 4  Z(MJ,1)                                                 
      DIMENSION  Z(MJ,1), X(MJ1,1)
C                                                                       
      L = MIN0( NMK,MJ1 )                                               
      K1 = LAG + 1                                                      
      KD1 = K1*ID + KSW                                                 
      N1 = L                                                            
C                                                                       
C          +-----------------+                                        +-
C          ! MATRIX X SET UP !                                        ! 
C          +-----------------+                                        +-
C                                                                       
      CALL  MSETX1( Z,N0,L,LAG,ID,MJ,MJ1,0,KSW,X )                      
C                                                                       
C          +----------------------------+                             +-
C          ! HOUSEHOLDER TRANSFORMATION !                             ! 
C          +----------------------------+                             +-
C                                                                       
cc      CALL  HUSHLD( X,D,MJ1,L,KD1 )                                     
      CALL  HUSHLD( X,MJ1,L,KD1 )                                     
C                                                                       
      IF( N1 .GE. NMK )     RETURN                                      
C                                                           +------->>--
C                                                           !           
  100 L = MIN0( NMK-N1,MJ1-KD1 )                                        
C                                                           !           
      LK = L + KD1                                                      
      N2 = N0 + N1                                                      
C                                                           !           
C          +-----------------+                              !         +-
C          ! MATRIX X SET UP !                              !         ! 
C          +-----------------+                              !         +-
C                                                           !           
      CALL  MSETX1( Z,N2,L,LAG,ID,MJ,MJ1,1,KSW,X )                      
C                                                           !           
C          +----------------------------+                   !         +-
C          ! HOUSEHOLDER TRANSFORMATION !                   !         ! 
C          +----------------------------+                   !         +-
C                                                           !           
cc      CALL  HUSHLD( X,D,MJ1,LK,KD1 )                                    
      CALL  HUSHLD( X,MJ1,LK,KD1 )                                    
C                                                           !         +-
C                                                           +<-- NO --!N
C                                                                     +-
      N1 = N1 + L                                                       
      IF( N1 .LT. NMK )     GO TO 100                                   
C                                                                       
C                                                                     +-
C                                                                     ! 
C                                                                     +-
      RETURN                                                            
C                                                                       
      END                                                               
C
C
cc      SUBROUTINE  MSDCOM( X,A,Y,D,N,M,ID,KSW,IPR,MJ,MJ2,MJ4,E,SD )      
      SUBROUTINE  MSDCOM( X,A,N,M,ID,KSW,IPR,MJ,E,SD )      
C                                                                       
C     THIS SUBROUTINE PRODUCES THE ONE-STEP AHEAD PREDICTION ERROR VARIA
C     MATRIX AND ITS DETERMINANT FOR THE MULTI-VARIATE AUTOREGRESSIVE MO
C     DEFINED BY THE AR-COEFFICIENT MATRICES GIVEN BY A.                
C       ----------------------------------------------------------------
C       THE FOLLOWING SUBROUTINE IS DIRECTLY CALLED BY THIS SUBROUTINE: 
C             HUSHLD                                                    
C       ----------------------------------------------------------------
C                                                                       
C       INPUTS:                                                         
C          X:    ((M+1)*ID+KSW)*((M+1)*ID+KSW) TRIANGULAR MATRIX        
C          A:    AR-COEFFICIENT MATRICES                                
C          Y:    WORKING AREA                                           
C          D:    WORKING AREA                                           
C          N:    DATA LENGTH                                            
C          M:    ORDER OF THE AR-MODEL                                  
C          ID:   DIMENSION OF THE PROCESS                               
C          KSW:  =0   CONSTANT TERM IS NOT INCLUDED AS A REGRESSOR      
C                =1   CONSTANT TERM IS INCLUDED AS THE FIRST REGRESSOR  
C          IPR:  PRINT OUT CONTROL                                      
C          MJ:   ABSOLUTE DIMENSION OF X IN THE MAIN PROGRAM            
C          MJ2:  ABSOLUTE DIMENSION OF A AND E IN THE MAIN PROGRAM      
C          MJ4:  ABSOLUTE DIMENSION OF Y IN THE MAIN PROGRAM            
C                                                                       
C       OUTPUTS:                                                        
C          E:    ONE STEP AHEAD PREDICTION ERROR VARIANCE MATRIX        
C          SD:   DETERMINANT OF SD                                      
C                                                                       
      IMPLICIT  REAL * 8 ( A-H,O-Z )                                    
cc      DIMENSION  X(MJ,1) , A(MJ2,MJ2,1) , E(MJ2,1) , Y(MJ4,1)           
cc      DIMENSION  D(1)                                                   
      DIMENSION  X(MJ,1) , A(ID,ID,M) , E(ID,ID) , Y((M+1)*ID,ID)
C                                                                       
      MD = M * ID                                                       
      M1D= MD + ID                                                      
      DO 30  JJ=1,ID                                                    
      DO 20  II=1,MD                                                    
      I0 = II + KSW                                                     
      SUM = 0.D0                                                        
      DO 10  J=II,MD                                                    
      J0 = J + KSW                                                      
      M0 = (J-1)/ID + 1                                                 
      J1 = J - (M0-1)*ID                                                
   10 SUM = SUM + X(I0,J0)*A(JJ,J1,M0)                                  
      J0 = MD + KSW + JJ                                                
   20 Y(II,JJ) = X(I0,J0) - SUM                                         
   30 CONTINUE                                                          
      DO 40  J=1,ID                                                     
      J0 = MD + KSW + J                                                 
      DO 40  I=1,ID                                                     
      I0 = MD + KSW + I                                                 
      I1 = MD + I                                                       
   40 Y(I1,J) = X(I0,J0)                                                
C                                                                       
cc      CALL  HUSHLD( Y,D,MJ4,M1D,ID )                                    
      CALL  HUSHLD( Y,M1D,M1D,ID )                                    
      SD = 1.D0                                                         
      DO 50  I=1,ID                                                     
   50 SD = SD * (Y(I,I)**2)/N                                           
      DO 70  I=1,ID                                                     
      DO 70  J=1,ID                                                     
      SUM = 0.D0                                                        
      DO 60  II=1,ID                                                    
   60 SUM = SUM + Y(II,I)*Y(II,J)                                       
   70 E(I,J) = SUM / N                                                  
cc      IF( IPR .LT. 1 )     RETURN                                       
cc      WRITE( 6,5 )                                                      
cc      DO 180  I=1,ID                                                    
cc  180 WRITE( 6,6 )   (E(I,J),J=1,ID)                                    
      RETURN                                                            
    5 FORMAT( /,1H ,'*  INNOVATION VARIANCE MATRIX  *' )
    6 FORMAT( 1H ,5D20.10 )                                             
      E N D                                                             
C
C
      SUBROUTINE  MSETX1( Z,N0,L,LAG,ID,MJ,MJ1,JSW,KSW,X )              
C                                                                       
C          +-----------------+                                          
C          ! MATRIX X SET UP !                                          
C          +-----------------+                                          
C                                                                       
C     THIS SUBROUTINE PREPARES DATA MATRIX X FROM DATA MATRIX Z(I,J) (I=
C     ,N0+K-KSW+L;J=1,...,ID) FOR AUTOREGRESSIVE MODEL FITTING.  X IS TH
C     AS INPUT TO SUBROUTINE HUSHLD.                                    
C       INPUTS:                                                         
C          Z:      ORIGINAL DATA MATRIX                                 
C          N0:     INDEX OF THE END POINT OF DISCARDED FORMER OBSERVATIO
C          L:      DIMENSION OF THE VECTOR OF NEW OBSERVATIONS          
C          LAG:    MAXIMUM TIME LAG OF THE MODEL                        
C          ID:     DIMENSION OF OBSERVATION                             
C          MJ:     ABSOLUTE DIMENSION OF Z                              
C          MJ1:    ABSOLUTE DIMENSION OF X                              
C          JSW:   =0   TO CONSTRUCT INITIAL L*(LAG+1) DATA MATRIX       
C                 =1   TO AUGMENT ORIGINAL (LAG+1)*(LAG+1) MATRIX X BY A
C                      L*(LAG+1) DATA MATRIX OF ADDTIONAL OBSERVATIONS  
C          KSW:   =0   CONSTANT VECTOR IS NOT INCLUDED AS A REGRESSOR   
C                 =1   CONSTANT VECTOR IS INCLUDED AS THE FIRST REGRESSO
C                                                                       
C       OUTPUT:                                                         
C          X:      L*(LAG+1) MATRIX            IF  JSW = 0              
C                  (LAG+1+L)*(LAG+1) MATRIX    IF  JSW = 1              
C                                                                       
      REAL * 8  X(MJ1,1) , Z(MJ,1)
CC      DIMENSION  Z(MJ,1)                                                
C                                                                       
      KD = LAG*ID + KSW                                                 
      KD1 = (LAG+1)*ID + KSW                                            
      I0 = 0                                                            
      IF( JSW .EQ. 1 )     I0 = KD1                                     
C                                                                       
      DO  30     II=1,L                                                 
      I1 = N0 + LAG + II                                                
      I2 = I0 + II                                                      
      DO  10     J=1,ID                                                 
      J2 = KD + J                                                       
   10 X(I2,J2) = Z(I1,J)                                                
      DO  30     JJ=1,LAG                                               
      I1 = I1 - 1                                                       
      J1 = (JJ-1)*ID + KSW                                              
      DO  20     J=1,ID                                                 
      J2 = J1 + J                                                       
   20 X(I2,J2) = Z(I1,J)                                                
   30 CONTINUE                                                          
C                                                                       
      IF( KSW .NE. 1 )     RETURN                                       
      DO  40     II=1,L                                                 
      I = II + I0                                                       
   40 X(I,1) = 1.D0                                                     
C                                                                       
      RETURN                                                            
      END                                                               
C
C
cc      SUBROUTINE  NRASPE( SGME2,A,B,L,K,H,TITLE )                       
      SUBROUTINE  NRASPE( SGME2,A,B,L,K,H,SXX )                       
C     THIS SUBROUTINE COMPUTES POWER SPECTRUM OF AN AR-MA PROCESS       
C     X(N)=A(1)X(N-1)+...+A(L)X(N-L)+E(N)+B(1)E(N-1)+...+B(K)E(N-K),    
C     WHERE E(N) IS A WHITE NOISE WITH ZERO MEAN AND VARIANCE EQUAL TO  
C     SGME2.  OUTPUTS PXX(I) ARE GIVEN AT FREQUENCIES I/(2*H)           
C     I=0,1,...,H.                                                      
C     REQUIRED INPUTS ARE;                                              
C     L,K,H,SGME2,(A(I),I=1,L), AND (B(I),I=1,K).                       
C        SGME2: NOISE VARIANCE                                          
C        A: AR-COEFFICIENTS                                             
C        B: MA-COEFFICIENTS                                             
C        L: ORDER OF AR                                                 
C        K: ORDER OF MA                                                 
C        N: LENGTH OF DATA                                              
C        H: NUMBER OF SEGMENTS OF FREQUENCY AXIS                        
C     0 IS ALLOWABLE AS L AND/OR K.                                     
C       ----------------------------------------------------------------
C       THE FOLLOWING SUBROUTINES ARE DIRECTLY CALLED BY THIS SUBROUTINE
C             FOUGER                                                    
C             SPEGRH                                                    
C       ----------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)                                          
CC      REAL*4  SXX(501) , TITLE(1)                                       
cc      REAL*4  TITLE(1)
      INTEGER H,H1                                                      
cc      DIMENSION A(501),B(501)                                           
cc      DIMENSION G(501),GR1(501),GI1(501),GR2(501),GI2(501)              
cc      DIMENSION SXX(501) , PXX(510)                                                
      DIMENSION A(L),B(K)
      DIMENSION G(L+K+1),GR1(H+1),GI1(H+1),GR2(H+1),GI2(H+1)
      DIMENSION SXX(H+1), PXX(H+1)
  310 H1=H+1                                                            
      L1=L+1                                                            
      K1=K+1                                                            
      G(1)=1.0                                                          
      IF(L.LE.0) GO TO 400                                              
      DO 10 I=1,L                                                       
      I1=I+1                                                            
   10 G(I1)=-A(I)                                                       
  400 CALL FOUGER(G,L1,GR1,GI1,H1)                                      
      G(1)=1.0                                                          
      IF(K.LE.0) GO TO 410                                              
      DO 20 I=1,K                                                       
      I1=I+1                                                            
   20 G(I1)=B(I)                                                        
  410 CALL FOUGER(G,K1,GR2,GI2,H1)                                      
      DO 30 I=1,H1                                                      
   30 PXX(I)=(GR2(I)**2+GI2(I)**2)/(GR1(I)**2+GI1(I)**2)*SGME2          
  510 DO 520 I=1,H1                                                     
  520 SXX(I)=DLOG10(PXX(I))                                             
 1000 CONTINUE                                                          
cc      CALL  SPEGRH( SXX,TITLE )                                         
C                                                                       
      RETURN                                                            
      END                                                               
C
C
      SUBROUTINE  PARCOR( AR,K,PAC )                                    
C                                                                       
C  ...  TRANSFORMATION FROM AR COEFFICIENTS TO PARCOR                   
C                                                                       
C       INPUTS:                                                         
C          AR:   VECTOR OF AR-COEFFICIENTS                              
C          K:    ORDER OF THE MODEL                                     
C       OUTPUT:                                                         
C          PAC:  VECTOR OF PARTIAL AUTOCORRELATIONS                     
C                                                                       
      IMPLICIT  REAL * 8( A-H,O-Z )                                     
cc      DIMENSION  AR(K) , PAC(K) , W(20)                                 
      DIMENSION  AR(K) , PAC(K) , W(K)
      DO 10  I=1,K                                                      
   10 PAC(I) = AR(I)                                                    
      IF( K .EQ. 1 )   RETURN                                           
      DO 20  JJ=1,K-1                                                   
      II = K - JJ                                                       
      S = 1.D0 - PAC(II+1)**2                                           
      DO 30  I=1,II                                                     
      J = II - I + 1                                                    
   30 W(I) = (PAC(I) + PAC(II+1)*PAC(J))/S                              
      II1 = II + 1                                                      
      I2 = II1 / 2                                                      
      IF( MOD( II,2 ) .EQ. 1 )  W(I2) = PAC(I2)/(1.D0 - PAC(II1))       
      DO 40  I=1,II                                                     
   40 PAC(I) = W(I)                                                     
   20 CONTINUE                                                          
C                                                                       
      RETURN                                                            
      E N D                                                             
C
      SUBROUTINE  RECOEF( X,M,K,MJ,A )                                  
C          +-------------------------------------+                      
C          ! REGRESSION COEFFICIENTS COMPUTATION !                      
C          +-------------------------------------+                      
C                                                                       
C     THIS SUBROUTINE PRODUCES REGRESSION COEFFICIENTS OF REGRESSION MOD
C     WITH M REGRESSORS FROM THE TRIANGULAR MATRIX X PREPARED BY SUBROUT
C     REDUCT.                                                           
C                                                                       
C       INPUTS:                                                         
C          X:     (K+1)*(K+1) TRIANGULAR MATRIX, WITH J-TH COLUMN REPRES
C                 J-TH REGRESSOR (J=1,K) AND (K+1)-TH REGRESSAND        
C          K:     NUMBER OF REGRESSORS                                  
C          M:     REGRESSION ON THE FIRST M REGRESSORS IS REQUESTED     
C          MJ:    ABSOLUTE DIMENSION OF X                               
C                                                                       
C       OUTPUT:                                                         
C          A(I) (I=1,M):   VECTOR OF REGRESSION COEFFICIENTS            
C                                                                       
      IMPLICIT REAL * 8 (A-H,O-Z )                                      
      DIMENSION X(MJ,1) , A(1)                                          
C                                                                       
      K1 = K + 1                                                        
      A(M) = X(M, K1) / X(M,M)                                          
      IF( M .EQ. 1 )     RETURN                                         
      MM1 = M - 1                                                       
      DO  10   II = 1,MM1                                               
      I = M - II                                                        
      SUM = X(I,K1)                                                     
      I1 = I + 1                                                        
      DO  20   J = I1,M                                                 
   20 SUM = SUM - A(J) * X(I,J)                                         
   10 A(I) = SUM / X(I,I)                                               
C                                                                       
      RETURN                                                            
C                                                                       
      E N D                                                             
C
C
cc      SUBROUTINE  REDATA( X,N,MT,TITLE )                                
      SUBROUTINE  REDATA( XS,X,N,XMEAN,SUM )
C                                                                       
C          +---------------------------------------+                    
C          ! ORIGINAL DATA INPUT AND MEAN DELETION !                    
C          +---------------------------------------+                    
C                                                                       
C     THIS SUBROUTINE IS USED FOR THE LOADING OF ORIGINAL DATA AND DELET
C     THE MEAN VALUE.  THE DATA IS LOADED THROUGH THE DEVICE SPECIFIED B
C     EACH DATA SET IS COMPOSED OF TITLE, DATA LENGTH, DATA FORMAT AND  
C     ORIGINAL DATA.                                                    
C                                                                       
C       INPUTS:                                                         
C         MT:     INPUT DEVICE SPECIFICATION                            
C         TITLE:  TITLE OF DATA                                         
C         N:      DATA LENGTH                                           
C         DFORM:  INPUT DATA FORMAT SPECIFICATION                       
C         X(I) (I=1,N):  ORIGINAL DATA                                  
C                                                                       
C       OUTPUTS:                                                        
C         X:   ORIGINAL DATA ( MEAN DELETED )                           
C         N:   DATA LENGTH                                              
C         TITLE:  TITLE OF DATA                                         
C                                                                       
      IMPLICIT REAL*8( A-H,O-Z )                                        
CC      REAL * 4     DFORM(20) , TITLE(20) , X(1)                         
      DIMENSION     XS(N), X(N)
cc	REAL * 4     DFORM(20) , TITLE(20)
C                                                                       
C       LOADING OF TITLE, DATA LENGTH, FORMAT SPECIFICATION AND DATA    
C                                                                       
cc      READ( MT,5 )     TITLE                                            
cc      READ( MT,1 )      N                                               
cc      READ( MT,5 )     DFORM                                            
cc      READ( MT,DFORM )     (X(I),I=1,N)                                 
C                                                                       
C       ORIGINAL DATA PRINT OUT                                         
C                                                                       
cc      WRITE( 6,6 )                                                      
cc      WRITE( 6,8 )     TITLE                                            
cc      WRITE( 6,9 )     N , (DFORM(I),I=1,20)                            
cc      WRITE( 6,4 )     (X(I),I=1,N)                                     
      DO 100 I=1,N
         X(I)=XS(I)
  100 CONTINUE
C                                                                  
C          MEAN DELETION                                                
C                                                                       
      SUM = 0.0D00                                                      
      DO 10     I=1,N                                                   
   10 SUM = SUM + X(I)                                                  
      XMEAN = SUM / DFLOAT(N)                                           
      DO 20   I=1,N                                                     
   20 X(I) = X(I) - XMEAN                                               
C                                                                       
C          VARIANCE COMPUTATION                                         
C                                                                       
      SUM = 0.D0                                                        
      DO  30     I=1,N                                                  
   30 SUM = SUM + X(I) * X(I)                                           
      SUM = SUM / DFLOAT(N)                                             
C                                                                       
C          MEAN AND VARIANCE PRINT OUT                                  
C                                                                       
cc      WRITE( 6,7 )   XMEAN, SUM                                         
C
      RETURN                                                            
C                                                                       
    1 FORMAT( 16I5 )                                                    
    4 FORMAT ( 1H ,10D13.5 )                                            
    5 FORMAT ( 20A4 )                                                   
    6 FORMAT( //1H ,'** ORIGINAL DATA **' )                             
    7 FORMAT ( 1H ,'MEAN =',D15.8,5X,'VARIANCE =',D15.8 )               
    8 FORMAT( 1H ,'TITLE  --  ',1H",20A4,1H" )                          
    9 FORMAT( 1H ,'N =',I5,5X,'FORMAT =',20A4 )                         
C                                                                       
      E N D                                                             
C
C
cc       SUBROUTINE  REDUCT( SETX,Z,D,NMK,N0,K,MJ1,LAG,X )                 
       SUBROUTINE  REDUCT( SETX,Z,NMK,N0,K,MJ1,LAG,X )                 
C                                                                       
C          +-----------------------+                                    
C          ! HOUSEHOLDER REDUCTION !                                    
C          +-----------------------+                                    
C                                                                       
C     THIS SUBROUTINE FIRST SETS UP DATA MATRIX X BY AUGMENTING         
C     SUCCESSIVELY SHIFTED ORIGINAL DATA VECTOR Z AND THEN TRANSFORMS   
C     X INTO TRIANGULAR FORM BY HOUSEHOLDER TRANSFORMATION.             
C       ----------------------------------------------------------------
C       THE FOLLOWING SUBROUTINES ARE DIRECTLY CALLED BY THIS SUBROUTINE
C             HUSHLD                                                    
C             (SETX)                                                    
C       ----------------------------------------------------------------
C                                                                       
C       INPUTS:                                                         
C          SETX:  EXTERNAL SUBROUTINE DESIGNATION                       
C          Z:     ORIGINAL DATA VECTOR                                  
C          D:     WORKING AREA                                          
C          NMK:   DIMENSION OF THE VECTOR OF REGRESSAND (Z(N0+LAG+1),...
C                 Z(N0+LAG+NMK))                                        
C          N0:    INDEX OF THE END POINT OF DISCARDED FORMER OBSERVATION
C          K:     NUMBER OF REGRESSORS                                  
C          MJ1:   ABSOLUTE DIMENSION OF X                               
C          LAG:   MAXIMUM TIME LAG OF THE MODEL                         
C                                                                       
C       OUTPUT:                                                         
C          X:     REDUCED MATRIX ( UPPER TRIANGULAR FORM )              
C                                                                       
      IMPLICIT  REAL*8( A-H,O-Z )                                       
CC      DIMENSION  X(MJ1,1) , D(1)                                        
CC      REAL * 4   Z(1)                                                   
      DIMENSION  X(MJ1,1) , Z(1)
C                                                                       
      L = MIN0( NMK,MJ1 )                                               
      K1 = K + 1                                                        
      N1 = L                                                            
C                                                                       
C          +-----------------+                                        +-
C          ! MATRIX X SET UP !                                        ! 
C          +-----------------+                                        +-
C                                                                       
      CALL  SETX( Z,N0,L,K,MJ1,0,LAG,X )                                
C                                                                       
C          +----------------------------+                             +-
C          ! HOUSEHOLDER TRANSFORMATION !                             ! 
C          +----------------------------+                             +-
C                                                                       
cc      CALL  HUSHLD( X,D,MJ1,L,K1 )                                      
      CALL  HUSHLD( X,MJ1,L,K1 )                                      
C                                                                       
      IF( N1 .GE. NMK )     RETURN                                      
C                                                           +------->>--
C                                                           !           
  100 L = MIN0( NMK-N1,MJ1-K1 )                                         
C                                                           !           
      LK = L + K1                                                       
      N2 = N0 + N1                                                      
C                                                           !           
C          +-----------------+                              !         +-
C          ! MATRIX X SET UP !                              !         ! 
C          +-----------------+                              !         +-
C                                                           !           
      CALL  SETX( Z,N2,L,K,MJ1,1,LAG,X )                                
C                                                           !           
C          +----------------------------+                   !         +-
C          ! HOUSEHOLDER TRANSFORMATION !                   !         ! 
C          +----------------------------+                   !         +-
C                                                           !           
cc      CALL  HUSHLD( X,D,MJ1,LK,K1 )                                     
      CALL  HUSHLD( X,MJ1,LK,K1 )                                     
C                                                           !         +-
C                                                           +<-- NO --!N
C                                                                     +-
      N1 = N1 + L                                                       
      IF( N1 .LT. NMK )     GO TO 100                                   
C                                                                       
C                                                                     +-
C                                                                     ! 
C                                                                     +-
      RETURN                                                            
C                                                                       
      END                                                               
C
C
cc      SUBROUTINE  SDCOMP( X,A,Y,N,K,MJ,SD )                             
      SUBROUTINE  SDCOMP( X,A,N,K,MJ,SD )                             
C                                                                       
C     THIS SUBROUTINE COMPUTES THE RESIDUAL VARIANCE OF THE REGRESSION M
C     WITH THE REGRESSION COEFFICIENTS A(I) (I=1,K).                    
C                                                                       
C       INPUTS:                                                         
C         X:   N*(K+1) TRIANGULAR MATRIX, OUTPUT OF SUBROUTINE REDUCT   
C         A:   VECTOR OF REGRESSION COEFFICIENTS                        
C         Y:   WORKING AREA                                             
C         N:   DATA LENGTH                                              
C         K:   HIGHEST ORDER OF THE MODEL                               
C         MJ:  ABSOLUTE DIMENSION OF X                                  
C                                                                       
C       OUTPUT:                                                         
C         SD:  RESIDUAL VARIANCE                                        
C                                                                       
      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
cc      DIMENSION  X(MJ,1) , Y(1) , A(1)                                  
      DIMENSION  X(MJ,1) , Y(K+1) , A(K)
C                                                                       
      K1 = K + 1                                                        
C                                                                       
      DO  20     I=1,K                                                  
      SUM = 0.D0                                                        
      DO  10     J=I,K                                                  
   10 SUM = SUM + X(I,J)*A(J)                                           
   20 Y(I) = SUM                                                        
      Y(K1) = 0.D0                                                      
C                                                                       
      SUM = 0.D0                                                        
      DO  30     I=1,K1                                                 
   30 SUM = SUM + (Y(I)-X(I,K1))**2                                     
      SD = SUM / N                                                      
C                                                                       
      RETURN                                                            
C                                                                       
      END                                                               
C
C
      SUBROUTINE  SETX1( Z,N0,L,K,MJ1,JSW,LAG,X )                       
C                                                                       
C          +-----------------+                                          
C          ! MATRIX X SET UP !                                          
C          +-----------------+                                          
C                                                                       
C     THIS SUBROUTINE PREPARES DATA MATRIX X FROM DATA VECTOR Z(I) (I=N0
C     N0+LAG+L) FOR AUTOREGRESSIVE MODEL FITTING.  X IS THEN USED AS INP
C     SUBROUTINE HUSHLD.                                                
C                                                                       
C       INPUTS:                                                         
C          Z:     ORIGINAL DATA VECTOR                                  
C          N0:    INDEX OF THE END POINT OF DISCARDED FORMER OBSERVATION
C                 (NEW OBSERVATION STARTS AT N0+LAG+1 AND ENDS AT N0+LAG
C          L:     DIMENSION OF THE VECTOR OF NEW OBSERVATIONS           
C          K:     NUMBER OF REGRESSORS (=LAG OR LAG+1)                  
C          MJ1:   ABSOLUTE DIMENSION OF X                               
C          JSW:   =0   TO CONSTRUCT INITIAL L*(K+1) DATA MATRIX         
C                 =1   TO AUGUMENT ORIGINAL (K+1)*(K+1) MATRIX X BY AN  
C                      L*(K+1) DATA MATRIX OF ADDITIONAL OBSERVATIONS   
C          LAG:   MAXIMUM TIME LAG OF THE MODEL                         
C                 =K   CONSTANT VECTOR IS NOT INCLUDED AS A REGRESSOR   
C                 <K   CONSTANT VECTOR IS INCLUDED AS THE FIRST REGRESSO
C                                                                       
C       OUTPUT:                                                         
C          X:      L*(K+1) MATRIX          IF  JSW = 0                  
C                 (K+1+L)*(K+1) MATRIX     IF  JSW = 1                  
C                                                                       
CC      REAL * 8  X(MJ1,1)                                                
CC      DIMENSION  Z(1)                                                   
      REAL * 8  X(MJ1,1) , Z(1)
C                                                                       
      KSW = 0                                                           
      IF( K .NE. LAG )     KSW = 1                                      
      K1 = K + 1                                                        
      I0 = 0                                                            
      IF( JSW .EQ. 1 )     I0 = K1                                      
      DO  10     I=1,L                                                  
      II = I + I0                                                       
      JJ = N0 + LAG + I                                                 
      X(II,K1) = Z(JJ)                                                  
      DO  10   J=1,LAG                                                  
      JJ = JJ - 1                                                       
      JKSW = J + KSW                                                    
   10 X(II,JKSW) = Z(JJ)                                                
C                                                                       
      IF( KSW .NE. 1 )     RETURN                                       
C                                                                       
      DO  20     I=1,L                                                  
   20 X(I,1) = 1.D0                                                     
      RETURN                                                            
C                                                                       
      E N D                                                             
C
C
      SUBROUTINE  SOLVE( C,R,ID,II,MJ2,MJ3,G )                          
C                                                                       
C     THIS SUBROUTINE SOLVES THE MATRIX EQUATION  C*G=R,  WHERE THE MATR
C     IS UPPER TRIANGULAR.                                              
C                                                                       
C       INPUTS:                                                         
C          C:      ID*ID COEFFICIENT MATRIX                             
C          R:      ID*ID MATRIX                                         
C          ID:     NUMBER OF ROWS AND COLUMNS OF C AND R                
C          II:     DESIGNATION OF POSITION WHERE THE SOLUSION TO BE STOR
C          MJ2:    ABSOLUTE DIMENSION OF C AND R                        
C          MJ3:    ABSOLUTE DIMENSION OF G AND H                        
C                                                                       
C       OUTPUT:                                                         
C          G:      SOLUTION OF THE MATRIX EQUATION                      
C                                                                       
      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
      DIMENSION  C(MJ2,1) , R(MJ2,1) , G(MJ3,MJ3,1)                     
C                                                                       
      IDM1 = ID - 1                                                     
      DO  10   J=1,ID                                                   
   10 G(J,ID,II) = R(ID,J)/C(ID,ID)                                     
C                                                                       
      DO  30   JJ=1,IDM1                                                
      I = ID - JJ                                                       
      IP1 = I + 1                                                       
      DO  30   J=1,ID                                                   
      SUM = 0.D0                                                        
      DO  20   L=IP1,ID                                                 
   20 SUM = SUM + G(J,L,II)*C(I,L)                                      
   30 G(J,I,II) = (R(I,J)-SUM) / C(I,I)                                 
C                                                                       
      RETURN                                                            
      END                                                               
C
C
cc      SUBROUTINE  SRCOEF( X,M,K,N,MJ,JND,IPR,A,SD )                     
      SUBROUTINE  SRCOEF( X,M,K,N,MJ,JND,IPR,A,SD,AIC,IFG,LU )          
C                                                                       
C     SUBSET REGRESSION COEFFICIENTS AND RESIDUAL VARIANCE COMPUTATION. 
C                                                                       
C                                                                       
C       INPUTS:                                                         
C         X:     TRIANGULAR MATRIX                                      
C         M:     NUMBER OF REGRESSORS                                   
C         K:     HEIGHEST ORDER OF THE MODELS                           
C         N:     DATA LENGTH                                            
C        JND(I):   (I=1,...,M)  SPECIFICATION OF I-TH REGRESSOR         
C       OUTPUTS:                                                        
C         A:     REGRESSION COEFFICIENTS                                
C         SD:    INNOVATION VARIANCE                                    
C                                                                       
      IMPLICIT  REAL * 8(A-H,O-Z)                                       
      DIMENSION  X(MJ,1) , A(1) , JND(1)                                
      K1 = K + 1                                                        
      M1 = M + 1                                                        
C                                                                       
C                                                                       
C          REGRESSION COEFFICIENTS COMPUTATION                          
C                                                                       
      L = JND(M)                                                        
      A(M) = X(M,K1) / X(M,L)                                           
      MM1 = M - 1                                                       
      IF( MM1 .EQ. 0 )     GO TO  60                                    
      DO  10     II=1,MM1                                               
      I = M - II                                                        
      SUM = X(I,K1)                                                     
      I1 = I + 1                                                        
      DO  20     J=I1,M                                                 
      L = JND(J)                                                        
   20 SUM = SUM - A(J) * X(I,L)                                         
      L = JND(I)                                                        
   10 A(I) = SUM / X(I,L)                                               
C                                                                       
C                                                                       
C          RESIDUAL VARIANCE AND AIC COMPUTATION                        
C                                                                       
   60 CONTINUE                                                          
      SD = 0.0D00                                                       
      DO  30     I=M1,K1                                                
   30 SD = SD + X(I,K1) * X(I,K1)                                       
      OSD = SD / N                                                      
      AIC = N * DLOG( OSD ) + 2.0D00 * M                                
C                                                                       
C                                                                       
C          REGRESSION COEFFICIENTS AND RESIDUAL VARIANCE PRINT OUT      
C                                                                       
      IF( IPR .LT. 2 )  RETURN                                          
      IF( IFG .EQ. 0 )  RETURN
cc      WRITE( 6,5 )                                                      
cc      WRITE( 6,6 )                                                      
      IF( IFG.NE.0 )  WRITE( LU,5 )
      IF( IFG.NE.0 )  WRITE( LU,6 )
      DO  40     I=1,M                                                  
      L = JND(I)                                                        
cc   40 WRITE( 6,7 )     L , A(I)                                         
cc      WRITE( 6,8 )     OSD , M , AIC                                    
   40 IF( IFG.NE.0 )  WRITE( LU,7 )     L , A(I)
      IF( IFG.NE.0 )  WRITE( LU,8 )     OSD , M , AIC
C                                                                       
      RETURN                                                            
cc    5 FORMAT( 1H0,10X,'SUBSET REGRESSION COEFFICIENTS' )                
    5 FORMAT( /,11X,'SUBSET REGRESSION COEFFICIENTS' )                
cc    6 FORMAT( 1H ,14X,1HI,12X,'A(I)' )                                  
    6 FORMAT( /,15X,1HI,12X,'A(I)' )                                  
    7 FORMAT( 1H ,10X,I5,5F20.10 )                                      
cc    8 FORMAT( 1H0,10X,'SD  = RESIDUAL VARIANCE    =',D19.12,/,11X,'M   =
    8 FORMAT( /,11X,'SD  = RESIDUAL VARIANCE    =',D19.12,/,11X,'M   =
     1 NUMBER OF PARAMETERS =',I3,/,11X,'AIC =  N*LOG(SD) + 2*M     =', 
     2F15.3 )                                                           
      E N D                                                             
C
C
      SUBROUTINE  TRIINV( X,M,MJ,MJ1,Y )                                
C                                                                       
C       LOWER TRIANGULAR MATRIX INVERSION                               
C                                                                       
C       INPUTS:                                                         
C          X:    TRIANGULAR MATRIX, DIAGONAL ELEMENTS ARE ASSUMED TO BE 
C          M:    DIMENSION OF MATRIX X                                  
C          MJ:   ABSOLUTE DIMENSION OF X                                
C          MJ1:  ABSOLUTE DIMENSION OF Y                                
C       OUTPUT:                                                         
C          Y:    INVERSE OF X                                           
C                                                                       
      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
      DIMENSION  X(MJ,1) , Y(MJ1,1)                                     
      MM1 = M - 1                                                       
      DO  10     I=1,MM1                                                
      DO  10     J=I,M                                                  
   10 Y(I,J) = 0.D0                                                     
      DO  20     I=1,M                                                  
   20 Y(I,I) = 1.D0                                                     
      DO  40     J=1,MM1                                                
      J1 = J + 1                                                        
      DO  40     I=J1,M                                                 
      SUM = 0.D0                                                        
      IJ = I - J                                                        
      DO  30     II=1,IJ                                                
      JJ = II + J - 1                                                   
   30 SUM = SUM + X(I,JJ) * Y(JJ,J)                                     
   40 Y(I,J) = -SUM                                                     
      RETURN                                                            
      END
C
cc      REAL FUNCTION  DMIN*8( X,N )                                      
      REAL*8 FUNCTION  DMIN( X,N )                                      
C                                                                       
C       THIS FUNCTION RETURNS THE MINIMUM VALUE AMONG X(I) (I=1,N).     
C                                                                       
C       INPUTS:                                                         
C          X:     VECTOR                                                
C          N:     DIMENSION OF VECTOR X                                 
C                                                                       
C       OUTPUT:                                                         
C          DMIN:  MINIMUM OF X(I) (I=1,N)                               
C                                                                       
      REAL * 8  X(1)                                                    
C                                                                       
      DMIN = X(1)                                                       
      DO  10   I=2,N                                                    
      IF( DMIN .GT. X(I) )         DMIN = X(I)                          
   10 CONTINUE                                                          
      RETURN                                                            
      END                                                               
C
C
      SUBROUTINE ARMASP( A,M,B,L,SIG2,NF,SP )
C
C  ...  Logarithm of the power spectrum of the ARMA model  ...
C
C     Inputs:
C        M:     AR order
C        L:     MA order
C        A(I):  AR coefficient
C        B(I):  MA coefficient
C        SIG2:  Innovation variance
C        NF:    Number of frequencies
C     Output:
C        SP(I): Power spectrum (in log scale)
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(M), B(L)
cc      DIMENSION SP(0:NF), H(0:500), FR(0:500), FI(0:500)
      DIMENSION SP(0:NF), H(0:M+L), FR(0:NF), FI(0:NF)
C
      H(0) = 1.0D0
      DO 10 I=1,M
   10 H(I) = -A(I)
C
      CALL  FOURIE( H,M+1,NF+1,FR,FI )
C
      DO 20 I=0,NF
   20 SP(I) = SIG2/( FR(I)**2 + FI(I)**2 )
C
      H(0) = 1.0D0
      DO 30 I=1,L
   30 H(I) = -B(I)
      CALL  FOURIE( H,L+1,NF+1,FR,FI )
      DO 40 I=0,NF
   40 SP(I) = SP(I)*( FR(I)**2 + FI(I)**2 )
C
      DO 50 I=0,NF
   50 SP(I) = DLOG10( SP(I) )
C
      RETURN
      E N D
C
C

      SUBROUTINE FOURIE( X,N,M,FC,FS )
C
C  ...  Discrete Fourier transformation by Goertzel method  ...
C
C     Inputs:
C        X(I):   data (I=1,N)
C        N:      data length
C        M:      number of Fourier components
C        FC(J):  Fourier cosine transform (J=1,M)
C        FS(J):  Fourier sine transform   (J=1,M)
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION  X(N), FC(M), FS(M)
      DATA  PI/3.14159265358979D0/
C
      W = PI/(M-1)
      DO 20 I=1,M
      CI = DCOS(W*(I-1))
      SI = DSIN(W*(I-1))
      T1 = 0.0
      T2 = 0.0
      DO 10 J=N,2,-1
        T0 = 2*CI*T1 - T2 + X(J)
        T2 = T1
   10   T1 = T0
      FC(I) = CI*T1 - T2 + X(1)
   20 FS(I) = SI*T1
C
      RETURN
      E N D
cc      SUBROUTINE  MOMENT( Y,N,YMEAN,VAR )
      SUBROUTINE  MOMENTK( Y,N,YMEAN,VAR )      
C
C  ...  Mean and variance of the data  ...
C
C     Inputs:
C       Y(I):   data
C       N:      data length
C     Outputs:
C       YMEAN:  mean
C       VAR:    variance
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION  Y(N)
C
      SUM = 0.0D0
      DO 10 I=1,N
   10 SUM = SUM + Y(I)
      YMEAN = SUM/N
      SUM = 0.0D0
      DO 20 I=1,N
   20 SUM = SUM + (Y(I)-YMEAN)**2
      VAR = SUM/N
C
      RETURN
      E N D
cc      SUBROUTINE  SMOOTH( F,M,MJ,NDIM,NS,NFE,NPE,VFS,VPS,XFS,XPS,
cc     *                    VSS,XSS )
      SUBROUTINE  SMOOTH( F,M,NDIM,NS,NFE,NPE,VFS,VPS,XFS,XPS,VSS,XSS )
C
C  ...  Fixed Interval Smoother (General Form)  ...
C
C     Inputs:
C        NS:     Start position of filtering
C        NFE:    End position of filtering
C        NPE:    End position of prediction
C        M:      Dimension of the state vector
C        F:      M*M matrix
C        MJ:     Adjustable dimension of XF, VF
C        NDIM:   Adjustable dimension of XFS, XPS, VFS, VPS
C        NMAX    Adjustable dimension of Y
C        VFS:    Covariance matrices of the filter
C        VPS:    Covariance matrices of the predictor
C        XFS:    Mean vectors of the filter
C        XPS:    Mean vectors of the predictor
C     Outputs:
C        VSS:    Covariance matrices of the smoother
C        XSS:    Mean vectors of the smoother
C
      IMPLICIT REAL*8(A-H,O-Z)
cc      DIMENSION  F(MJ,MJ)
cc      DIMENSION  XS(40), VS(40,40), VP(40,40)
cc      DIMENSION  XFS(MJ,NDIM), XPS(MJ,NDIM), XSS(MJ,NDIM)
cc      DIMENSION  VFS(MJ,MJ,NDIM), VPS(MJ,MJ,NDIM), VSS(MJ,MJ,NDIM)
cc      DIMENSION  WRK(40,40), SGAIN(40,40)
      DIMENSION  F(M,M)
      DIMENSION  XS(M), VS(M,M), VP(M,M)
      DIMENSION  XFS(M,NDIM), XPS(M,NDIM), XSS(M,NDIM)
      DIMENSION  VFS(M,M,NDIM), VPS(M,M,NDIM), VSS(M,M,NDIM)
      DIMENSION  WRK(M,M), SGAIN(M,M)
C
C
C  ...  SMOOTHING  ...
C
      DO 10 II=NFE,NPE
      DO 10  I=1,M
      XSS(I,II)   = XFS(I,II)
      DO 10  J=1,M
   10 VSS(I,J,II) = VFS(I,J,II)
      DO 20  I=1,M
      XS(I)   = XFS(I,NFE)
      DO 20  J=1,M
   20 VS(I,J) = VFS(I,J,NFE)
C
      DO 500  II=NFE-1,NS,-1
C
      NZERO = 0
      DO 100 I=1,M
  100 IF( VFS(I,I,II).GT.1.0D-12 )  NZERO = NZERO + 1
C
      IF( NZERO.EQ.0 )  THEN
         DO 110  I=1,M
         XS(I)     = XFS(I,II)
         XSS(I,II) = XFS(I,II)
         DO 110  J=1,M
         VS(I,J)     = VFS(I,J,II)
  110    VSS(I,J,II) = VFS(I,J,II)
C
      ELSE
      DO 410  I=1,M
      DO 410  J=1,M
  410 VP(I,J) = VPS(I,J,II+1)
C
cc      CALL  GINVRS( VP,VDET,M,40 )
      CALL  GINVRS( VP,VDET,M )
C
      DO 425  I=1,M
      DO 425  J=1,M
      SUM = 0.0D0
      DO 420  IJ=1,M
  420 SUM = SUM + VFS(I,IJ,II)*F(J,IJ)
  425 WRK(I,J) = SUM
C
      DO 440  I=1,M
      DO 440  J=1,M
      SUM = 0.0D0
      DO 430 IJ=1,M
  430 SUM = SUM + WRK(I,IJ)*VP(IJ,J)
  440 SGAIN(I,J) = SUM
C
      DO 450  I=1,M
      XS(I) = XFS(I,II)
      DO 450  J=1,M
      WRK(I,J) = 0.0D0
  450 VS (I,J) = VFS(I,J,II)
C
      DO 460  J=1,M
      DO 460  I=1,M
  460 XS(I) = XS(I) + SGAIN(I,J)*(XSS(J,II+1) - XPS(J,II+1))
C
      DO 470  J=1,M
      DO 470 IJ=1,M
      DO 470  I=1,M
  470 WRK(I,J)=WRK(I,J) + SGAIN(I,IJ)*(VSS(IJ,J,II+1)-VPS(IJ,J,II+1))
C
      DO 480  J=1,M
      DO 480 IJ=1,M
      DO 480  I=1,M
  480 VS(I,J) = VS(I,J) + WRK(I,IJ)*SGAIN(J,IJ)
      DO 485 I=1,M
  485 IF( VS(I,I).LT.0.0D0 )  VS(I,I) = 0.0D0
C
      DO 490  I=1,M
      XSS(I,II) = XS(I)
      DO 490  J=1,M
  490 VSS(I,J,II) = VS(I,J)
      END IF
C
  500 CONTINUE
C
      RETURN
      E N D
cc      SUBROUTINE  GINVRS( A,DET,M,MJ )
      SUBROUTINE  GINVRS( A,DET,M )
C
C  ...  Generalized inverse of a square matrix A  ...
C
C     Inputs:
C        A:     M*M matrix
C        M:     Dimension of A
C        MJ:    Adjustable dimension of A
C     Outputs:
C        A:     Generalize inverse of A
C        DET:   Determinant of A
C
      IMPLICIT REAL*8(A-H,O-Z)
cc      DIMENSION  A(MJ,MJ), IND(50)
      DIMENSION  A(M,M), IND(M+1)
      EPS = 1.0D-10
C
      DO 10  I=1,M
   10 IND(I) = I
cc------------
      I0 = 0
      LMAX = 0
cc------------
      DO 60  L=1,M
      AMAX = 0.0D0
      DO 20  I=L,M
      IF( A(IND(I),IND(I)).GT.AMAX )  THEN
        AMAX = A(IND(I),IND(I))
        I0 = I
      END IF
   20 CONTINUE
      IF( AMAX.GT.EPS*A(IND(1),IND(1)) )  THEN
         IMAX = IND(I0)
         DO 30  I=I0,L+1,-1
   30    IND(I) = IND(I-1)
         IND(L) = IMAX
         LMAX   = L
         DO 40  I=L+1,M
         A(IND(I),IMAX) = -A(IND(I),IMAX)/A(IMAX,IMAX)
         DO 40  J=L+1,M
   40    A(IND(I),IND(J)) = A(IND(I),IND(J))
     *                    + A(IND(I),IMAX)*A(IMAX,IND(J))
      ELSE
         DO 50  I=L,M
         DO 50  J=L,M
   50    A(IND(I),IND(J)) = 0.0D0
         GO TO 70
      END IF
   60 CONTINUE
C
   70 DET = 1.0D0
      DO 80  I=1,M
C     DET = DET*A(IND(I),IND(I))
C     IF( A(IND(I),IND(I)).GT.EPS*AMAX )  THEN
      IF( A(IND(I),IND(I)).GT.0.0D0 )  THEN
        A(IND(I),IND(I)) = 1.0D0/A(IND(I),IND(I))
      ELSE
        A(IND(I),IND(I)) = 0.0D0
      END IF
   80 CONTINUE
C
      MS = MIN0( M-1,LMAX )
      DO 200  L=MS,1,-1
      DO 100 J=L+1,M
      SUM = 0.0D0
      DO 90  I=L+1,M
   90 SUM = SUM + A(IND(I),IND(L))*A(IND(I),IND(J))
  100 A(IND(L),IND(J)) = SUM
      SUM = A(IND(L),IND(L))
      DO 110  I=L+1,M
  110 SUM = SUM + A(IND(L),IND(I))*A(IND(I),IND(L))
      A(IND(L),IND(L)) = SUM
      DO 120  I=L+1,M
  120 A(IND(I),IND(L)) = A(IND(L),IND(I))
      IMAX = IND(L)
      DO 130 I=L+1,M
      IF( IMAX.GT.IND(I) )  THEN
         IND(I-1) = IND(I)
         IND(I)   = IMAX
      END IF
  130 CONTINUE
  200 CONTINUE
C
      RETURN
      E N D
      
