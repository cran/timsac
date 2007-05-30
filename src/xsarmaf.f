      SUBROUTINE XSARMAF( YS,N,IQ,IP,P01,G1,TL1,P02,G2,ALPHB,ALPHA,TL2,
     *                    SIGMA2 )
C
	INCLUDE 'timsac_f.h'
C
cc      PROGRAM XSARMA                                                    
C.......................................................................
C.....PLANNED BY H.AKAIKE...............................................
C.....DESIGNED BY H.AKAIKE..............................................
C.....PROGRAMMED BY E.ARAHATA...........................................
C.....ADDRESS: THE INSTITUTE OF STATISTICAL MATHEMATICS, 4-6-7 MINAMI-AZ
C..............MINATO-KU, TOKYO 106, JAPAN..............................
C.....DATE OF THE LATEST REVISION:  MAR. 6,1979.........................
C.......................................................................
C.....THIS PROGRAM WAS ORIGINALLY PUBLISHED IN "TIMSAC-78", BY H.AKAIKE,
C.....G.KITAGAWA, E.ARAHATA AND F.TADA, COMPUTER SCIENCE MONOGRAPHS, NO.
C.....THE INSTITUTE OF STATISTICAL MATHEMATICS, TOKYO, 1979.............
C.......................................................................
C     TIMSAC 78.5.2                                                     
C     __                                 _      __ __                   
C     EXACT MAXIMUM LIKELIHOOD METHOD OF SCALAR AR-MA MODEL FITTING     
C                                                                       
C-----------------------------------------------------------------------
C     THIS PROGRAM PRODUCES EXACT MAXIMUM LIKELIHOOD ESTIMATES OF THE   
C     PARAMETERS OF A SCALAR AR-MA MODEL.                               
C                                                                       
C     THE AR-MA MODEL IS GIVEN BY                                       
C                                                                       
C     Y(I)+B(1)Y(I-1)+...+B(IQ)Y(I-IQ)=X(I)+A(1)X(I-1)+...+A(IP)X(I-IP),
C                                                                       
C     WHERE X(I) IS A ZERO MEAN WHITE NOISE.                            
C                                                                       
C     ----------------------------------------------------------------- 
C       REFERENCE:                                                      
C          H.AKAIKE(1978), "COVARIANCE MATRIX COMPUTATION OF THE STATE  
C          VARIABLE OF A STATIONARY GAUSSIAN PROCESS,",  RESEARCH MEMO. 
C          NO.139, THE INSTITUTE OF STATISTICAL MATHEMATICS; TOKYO.     
C          TO BE PUBLISHED IN ANN. INST. STATIST. MATH..                
C     ----------------------------------------------------------------- 
C     THE FOLLOWING SUBROUTINES ARE DIRECTLY CALLED BY THIS PROGRAM:    
C             SDATPR                                                    
C             SMINOP                                                    
C             SUBRST                                                    
C     ----------------------------------------------------------------- 
C     THE FOLLOWING INPUTS ARE REQUESTED BY SUBROUTINE SDATPR:          
C          IQ:  AR-ORDER                                                
C          (B(I),I=1,IQ):  INITIAL ESTIMATES OF AR COEFFICIENTS         
C          IP:  MA-ORDER                                                
C          (A(I),I=1,IP):  INITIAL ESTIMATES OF MA COEFFICIENTS         
C          TITL:  TITLE OF THE DATA                                     
C          N:  DATA LENGTH                                              
C          (DFORM(I),I=1,20):  INPUT FORMAT SPECIFICATION IN ONE CARD,  
C                              EXAMPLE,                                 
C                              (8F10.4)                                 
C          (Y(I),I=1,N):  ORIGINAL DATA                                 
C                                                                       
C     OUTPUTS:                                                          
C          TITL:  TITLE OF THE DATA                                     
C          IQ:  AR-ORDER                                                
C          (B(I),I=1,IQ):  AR COEFFICIENTS                              
C          IP:  MA-ORDER                                                
C          (A(I),I=1,IP):  MA COEFFICIENTS                              
C          SIGMA2:  WHITE NOISE VARIANCE                                
C-----------------------------------------------------------------------
C                                                                       
cc      !DEC$ ATTRIBUTES DLLEXPORT :: XSARMAF
C
      IMPLICIT REAL*8(A-H,O-Z)                                          
cc      REAL*4 TITL(20)                                                   
cc      DIMENSION Y(2000),P0(50)                                          
      DIMENSION YS(N), Y(N), P01(IP+IQ)
      DIMENSION P02(IP+IQ), G1(IP+IQ), G2(IP+IQ)
      DIMENSION ALPHB(IQ), ALPHA(IP)
C
cc      CHARACTER(100)  IFLNAM,OFLNAM,MFLNAM
cc      CALL FLNAM3( IFLNAM,OFLNAM,MFLNAM,NFL )
cc      IF ( NFL.EQ.0 ) GO TO 999
cc      IF ( (NFL.EQ.2).OR.(NFL.EQ.4) ) THEN
cc         OPEN (6,FILE=OFLNAM,ERR=900,IOSTAT=IVAR)
cc      ELSE
cc         CALL SETWND
cc      END IF
cc      OPEN( 5,FILE=IFLNAM,ERR=910,IOSTAT=IVAR,STATUS='OLD' )
cc      IF ((NFL.EQ.3) .OR. (NFL.EQ.4)) THEN
cc         OPEN (7,FILE=MFLNAM,ERR=920,IOSTAT=IVAR)
cc      ELSE
cc         OPEN (7,FILE='xsarma.out',ERR=930,IOSTAT=IVAR)
cc      END IF
C
C                                                                       
C-----  READ IN AND PRINT OUT OF INITIAL CONDITION  -----               
cc      CALL SDATPR(TITL,Y,N,P0,IQ,IP)                                    
      CALL SDATPR(YS,Y,N,P01,IQ,IP)
cc      CLOSE( 5 )
cc#ifdef __linux__
ccC     reopen #5 as stdin
cc      OPEN(5, FILE='/dev/fd/0')
cc#endif
ccC /* __linux */
C                                                                       
C-----  MINIMIZATION OF TL(=-2)LOG LIKELIHOOD)                          
C       BY DAVIDON'S VARIANCE ALGORITHM  -----                          
cc      CALL SMINOP(TL,SIGMA2,Y,N,P0,IQ,IP)                               
      IPRNT = 0
C     IPRNT=0:  NOT TO PRINT OUT INTERMEDIATE RESULTS
C     IPRNT=1:  TO PRINT OUT INTERMEDIATE RESULTS ( LU: UNIT NUMBER)
      CALL SMINOP( TL1,TL2,SIGMA2,Y,N,P01,G1,P02,G2,ALPHB,ALPHA,IQ,IP,
     *             IPRNT,LU )
C                                                                       
C-----  PRINT AND PUNCH OUT OF FINAL RESULT  -----                      
cc      CALL SUBRST(TITL,TL,SIGMA2,P0,IQ,IP)                              
cc      GO TO 990
C                                                                       
C
cc  900 CONTINUE
cc      WRITE(6,690) IVAR,OFLNAM
cc  690 FORMAT(' !!! Output_Data_File OPEN ERROR ',I8,//,5X,100A)
cc      GO TO 999
C
cc  910 CONTINUE
cc      IF ((NFL.EQ.2) .OR. (NFL.EQ.4)) CLOSE(6)
cc#ifdef __linux__
ccC     reopen #6 as stdout
cc      IF ((NFL.EQ.2) .OR. (NFL.EQ.4)) OPEN(6, FILE='/dev/fd/1')
cc#endif
ccC /* __linux__ */
cc      WRITE(6,691) IVAR,IFLNAM
cc  691 FORMAT(' !!! Input_Data_File OPEN ERROR ',I8,//,5X,100A)
cc      GO TO 999
C
cc  920 CONTINUE
cc      CLOSE(5)
cc      IF ((NFL.EQ.2) .OR. (NFL.EQ.4)) CLOSE(6)
cc#ifdef __linux__
ccC     reopen #5 as stdin, #6 as stdout
cc      OPEN(5, FILE='/dev/fd/0')
cc      IF ((NFL.EQ.2) .OR. (NFL.EQ.4)) OPEN(6, FILE='/dev/fd/1')
cc#endif
ccC /* __linux__ */
cc      WRITE(6,692) IVAR,MFLNAM
cc  692 FORMAT(' !!! Intermediate_Data_File OPEN ERROR ',I8,//,5X,100A)
cc      GO TO 999
C
cc  930 CONTINUE
cc      CLOSE(5)
cc      IF ((NFL.EQ.2) .OR. (NFL.EQ.4)) CLOSE(6)
cc#ifdef __linux__
ccC     reopen #5 as stdin, #6 as stdout
cc      OPEN(5, FILE='/dev/fd/0')
cc      IF ((NFL.EQ.2) .OR. (NFL.EQ.4)) OPEN(6, FILE='/dev/fd/1')
cc#endif
ccC /* __linux__ */
cc      WRITE(6,693) IVAR
cc  693 FORMAT(' !!! xsarma.out  OPEN ERROR ',I8)
cc      GO TO 999
C
cc  990 CONTINUE
cc      CLOSE( 7 )
cc      IF ((NFL.EQ.2) .OR. (NFL.EQ.4)) CLOSE(6)
cc#ifdef __linux__
ccC     reopen #6 as stdout
cc      IF ((NFL.EQ.2) .OR. (NFL.EQ.4)) OPEN(6, FILE='/dev/fd/1')
cc#endif
ccC /* __linux__ */
cc  999 CONTINUE
cc      STOP                                                              
      RETURN
      END                                                               
cc      SUBROUTINE ARCHCK(A,M,ICOND)                                      
      SUBROUTINE ARCHCK(A,ALPH,M,ICOND)                                      
C                                                                       
C-----------------------------------------------------------------------
C     THIS SUBROUTINE CHECKS STABILITY OF AR OR MA PART.                
C                                                                       
C     INPUTS:                                                           
C          (A(I),I=1,M):  AR OR MA COEFFICIENTS                         
C          M:  AR- OR MA-ORDER                                          
C                                                                       
C     OUTPUTS:                                                          
C          (A(I),I=1,M):  AR OR MA COEFFICIENTS                         
C          ICOND:  ICOND=0, WHEN STABLE                                 
C                  ICOND=1, WHEN NOT STABLE                             
C-----------------------------------------------------------------------
C                                                                       
      IMPLICIT REAL*8(A-H,O-Z)                                          
cc      DIMENSION A(50),B(50),ALPH(50)                                    
cc      COMMON /CMALPH/ALPH                                               
      DIMENSION A(M),B(M),ALPH(M)                                    
      DATA CST0,CST1,CST099/0.0D-00,1.0D-00,0.99999D-00/                
C                                                                       
C                                                                       
C                                                                       
C-----  STABILITY CHECK  -----                                          
cc      DO 10 I=1,50                                                      
      DO 10 I=1,M
   10 B(I)=CST0                                                         
      DO 20 I=1,M                                                       
      MI=M-I                                                            
      MIP1=MI+1                                                         
      AT=A(MIP1)                                                        
      IF(DABS(AT).LT.CST099) GO TO 210                                  
      ICOND=1                                                           
      AT=CST099*AT/DABS(AT)                                             
  210 ALPH(MIP1)=AT                                                     
      IF(MI.EQ.0) GO TO 20                                              
      CT=CST1/(CST1-AT**2)                                              
      DO 25 J=1,MI                                                      
      KI=MIP1-J                                                         
      B(J)=A(KI)                                                        
   25 CONTINUE                                                          
      DO 30 K=1,MI                                                      
      A(K)=(A(K)-AT*B(K))*CT                                            
   30 CONTINUE                                                          
   20 CONTINUE                                                          
C                                                                       
C-----  RECOVERY OF COEFFICIENTS  -----                                 
      DO 400 I=1,M                                                      
      D=ALPH(I)                                                         
      A(I)=D                                                            
      IF(I.EQ.1) GO TO 410                                              
      IM1=I-1                                                           
      DO 420 J=1,IM1                                                    
      A(J)=A(J)+D*B(J)                                                  
  420 CONTINUE                                                          
  410 CONTINUE                                                          
      IP1=I+1                                                           
      DO 421 K=1,I                                                      
      IK=IP1-K                                                          
  421 B(K)=A(IK)                                                        
  400 CONTINUE                                                          
C                                                                       
      RETURN                                                            
      END                                                               
cc      SUBROUTINE  FUNCT2( F,SD,Y,N,P0,IQ,IP )                           
      SUBROUTINE  FUNCT2( F,SD,Y,N,P0,IQ,IP,IR )
C                                                                       
C-----------------------------------------------------------------------
C     THIS SUBROUTINE COMPUTES SD AND F(=(-2)LOG LIKELIHOOD) BY         
C     A PROCEDURE OF MORF, SIDHU AND KAILATH (IEEE TRANS. AUTOMAT.      
C     CONTR., AC19, 315-323,1974).                                      
C       ----------------------------------------------------------------
C       THE FOLLOWING SUBROUTINE IS DIRECTLY CALLED BY THIS SUBROUTINE: 
C             SUBPM                                                     
C       ----------------------------------------------------------------
C                                                                       
C     INPUTS:                                                           
C          N:  DATA LENGTH                                              
C          (Y(I),I=1,N):  ORIGINAL DATA, MEAN DELETED                   
C          IQ:  AR-ORDER                                                
C          IP:  MA-ORDER                                                
C          (P0(I),I=1,IPQ):  THE VECTOR OF AR AND MA COEFFICIENTS (IPQ=I
C                                                                       
C     OUTPUTS:                                                          
C          F:  (-2)LOG LIKELIHOOD                                       
C          SD:  WHITE NOISE VARIANCE                                    
C-----------------------------------------------------------------------
C                                                                       
      IMPLICIT REAL*8(A-H,O-Z)                                          
cc      DIMENSION Y(2000)                                                 
cc      DIMENSION A(50),B(50),P(51,51),P0(50)                             
cc      DIMENSION P1(51),AKI(51),YY(51),PY(51),Z(51),PZ(51)               
      DIMENSION Y(N)                                                 
      DIMENSION A(IR),B(IR),P(IR,IR),P0(IP+IQ)                             
      DIMENSION P1(IR),AKI(IR),YY(IR),PY(IR),Z(IR),PZ(IR)
cc      DATA A,B,AKI/151*0.0D-00/                                         
      DATA CST0,CST1,CSTA/0.0D-00,1.0D-00,0.1D-05/                      
C                                                                       
      DO 600 I=1,IR
         A(I)=CST0
         B(I)=CST0
         AKI(I)=CST0
  600 CONTINUE
C                                                                       
      IF(IQ.EQ.0) GO TO 620                                             
      DO 610 I=1,IQ                                                     
      B(I)=P0(I)                                                        
  610 CONTINUE                                                          
  620 IF(IP.EQ.0) GO TO 640                                             
      DO 630 I=1,IP                                                     
      II=IQ+I                                                           
      A(I)=P0(II)                                                       
  630 CONTINUE                                                          
  640 CONTINUE                                                          
C                                                                       
cc      IP1=IP+1                                                          
cc      IR=MAX0(IQ,IP1)                                                   
C                                                                       
C-----  P MATRIX COMPUTATION  -----                                     
      CALL SUBPM(P,B,A,IQ,IP,IR)                                        
C                                                                       
C-----  INITIAL CONDITION  -----                                        
      REM=P(1,1)                                                        
      DO 20 I=1,IR                                                      
      P1(I)=P(I,1)                                                      
   20 CONTINUE                                                          
      IRM1=IR-1                                                         
      IRP1=IR+1                                                         
      IF(IRM1.LE.0) GO TO 510                                           
      DO 50 I=1,IRM1                                                    
      AKI(I)=P1(I+1)                                                    
   50 CONTINUE                                                          
  510 CONTINUE                                                          
      SUM=CST0                                                          
      DO 60 I=1,IR                                                      
      I1=IRP1-I                                                         
      SUM=SUM+B(I)*P1(I1)                                               
   60 CONTINUE                                                          
      AKI(IR)=-SUM                                                      
      AMI=-CST1/REM                                                     
      DO 100 I=1,IR                                                     
      YY(I)=AKI(I)                                                      
  100 CONTINUE                                                          
      EM=Y(1)                                                           
      DO 110 I=1,IR                                                     
      Z(I)=CST0                                                         
  110 CONTINUE                                                          
C                                                                       
      SERI=EM**2/REM                                                    
      SLR=DLOG(REM)                                                     
C******************************                                         
      INDX = 1
      DO 1000 NS=2,N                                                    
      Y2=YY(1)*YY(1)                                                    
C-----  RE COMPUTATION  -----                                           
      REW=REM+AMI*Y2                                                    
      INDX=NS                                                           
C                                                                       
C-----  PHI*Z COMPUTATION  -----                                        
      IF(IRM1.LE.0) GO TO 520                                           
      DO 210 I=1,IRM1                                                   
      PZ(I)=Z(I+1)                                                      
  210 CONTINUE                                                          
  520 CONTINUE                                                          
      SUM=CST0                                                          
      DO 220 I=1,IR                                                     
      I1=IRP1-I                                                         
      SUM=SUM+B(I)*Z(I1)                                                
  220 CONTINUE                                                          
      PZ(IR)=-SUM                                                       
C----- Z COMPUTATION  -----                                             
      ERI=EM/REM                                                        
      DO 230 I=1,IR                                                     
      Z(I)=PZ(I)+ERI*AKI(I)                                             
  230 CONTINUE                                                          
C                                                                       
C-----  PHI*YY COMPUTATION  -----                                       
      IF(IRM1.LE.0) GO TO 530                                           
      DO 310 I=1,IRM1                                                   
      PY(I)=YY(I+1)                                                     
  310 CONTINUE                                                          
  530 CONTINUE                                                          
      SUM=CST0                                                          
      DO 320 I=1,IR                                                     
      I1=IRP1-I                                                         
      SUM=SUM+B(I)*YY(I1)                                               
  320 CONTINUE                                                          
      PY(IR)=-SUM                                                       
C                                                                       
C-----  K COMPUTATION  -----                                            
      AMY=AMI*YY(1)                                                     
      DO 330 I=1,IR                                                     
      AKI(I)=AKI(I)+AMY*PY(I)                                           
  330 CONTINUE                                                          
C                                                                       
C-----  M COMPUTATION  -----                                            
      AMRI=AMI/REM                                                      
      AMI=AMI*(CST1+AMRI*Y2)                                            
C                                                                       
C-----  YY COMPUTATION  -----                                           
      YRI=YY(1)/REW                                                     
      DO 400 I=1,IR                                                     
      YY(I)=PY(I)-YRI*AKI(I)                                            
  400 CONTINUE                                                          
C                                                                       
C-----  E COMPUTATION  -----                                            
      EW=Y(NS)-Z(1)                                                     
C                                                                       
C-----  SERI, SLR COMPUTATION  -----                                    
      SERI=SERI+EW**2/REW                                               
      SLR=SLR+DLOG(REW)                                                 
      REM=REW                                                           
      EM=EW                                                             
      IF(DABS(REW-CST1).LT.CSTA) GO TO 1100                             
 1000 CONTINUE                                                          
C                                                                       
 1100 CONTINUE                                                          
      IF(INDX.GE.N) GO TO 1500                                          
      INDX1=INDX+1                                                      
      DO 1110 NS=INDX1,N                                                
C-----  PHI*Z  -----                                                    
      IF(IRM1.LE.0) GO TO 540                                           
      DO 1210 I=1,IRM1                                                  
      PZ(I)=Z(I+1)                                                      
 1210 CONTINUE                                                          
  540 CONTINUE                                                          
      SUM=CST0                                                          
      DO 1220 I=1,IR                                                    
      I1=IRP1-I                                                         
      SUM=SUM+B(I)*Z(I1)                                                
 1220 CONTINUE                                                          
      PZ(IR)=-SUM                                                       
C                                                                       
C-----  Z COMPUTATION  -----                                            
      DO 1230 I=1,IR                                                    
      Z(I)=PZ(I)+EM*AKI(I)                                              
 1230 CONTINUE                                                          
C                                                                       
C-----  E COMPUTATION  -----                                            
      EW=Y(NS)-Z(1)                                                     
C                                                                       
C-----  SERI COMPUTATION  -----                                         
      SERI=SERI+EW**2                                                   
      EM=EW                                                             
 1110 CONTINUE                                                          
 1500 CONTINUE                                                          
C******************************                                         
C                                                                       
C                                                                       
      AN=N                                                              
      SD=SERI/AN                                                        
      F=SLR+AN*DLOG(SD)                                                 
C                                                                       
      IF(IQ.EQ.0) GO TO 1620                                            
      DO 1610 I=1,IQ                                                    
      P0(I)=B(I)                                                        
 1610 CONTINUE                                                          
 1620 IF(IP.EQ.0) GO TO 1640                                            
      DO 1630 I=1,IP                                                    
      II=IQ+I                                                           
      P0(II)=A(I)                                                       
 1630 CONTINUE                                                          
 1640 CONTINUE                                                          
 2100 CONTINUE                                                          
C                                                                       
      RETURN                                                            
      END                                                               
cc      SUBROUTINE MSDAV2(PHAI,SIGMA2,G,C,Y,N,X,IQ,IP,ISWRO,IPRNT)        
      SUBROUTINE MSDAV2(PHAI,SIGMA2,G,C,Y,N,X,IQ,IP,ISWRO,VD,IPRNT,LU)
C                                                                       
C-----------------------------------------------------------------------
C     DAVIDON'S (MINIMIZATION) PROCEDURE                                
C       THE FOLLOWING SUBROUTINES ARE DIRECTLY CALLED BY THIS SUBROUTINE
C             SGRAD                                                     
C             ARCHCK                                                    
C                                                                       
C     INPUTS:                                                           
C          N:  DATA LENGTH                                              
C          (Y(I),I=1,N):  ORIGINAL DATA, MEAN DELETED                   
C          IQ:  AR-ORDER                                                
C          IP:  MA-ORDER                                                
C          (X(I),I=1,IPQ):  THE VECTOR OF AR AND MA COEFFICIENTS (IPQ=IP
C          PHAI:  (-2)LOG LIKELIHOOD                                    
C          SIGMA2:  WHITE NOISE VARIANCE                                
C          (G(I),I=1,IPQ):  GRADIENT                                    
C          (C(I),I=1,IPQ):  CORRECTION TERM                             
C          ISWRO:  ITERATION COUNT OF SUBROUTINE MSDAV2                 
C          ((VD(I,J),I=1,IPQ),J=1,IPQ):  INVERSE OF HESSIAN             
C                                                                       
C     OUTPUTS:                                                          
C          PHAI:  NEW PHAI                                              
C          SIGMA2:  NEW SIGMA2                                          
C          (G(I),I=1,IPQ):  NEW GRADIENT                                
C          (X(I),I=1,IPQ): NEW X(I)                                     
C          ISWRO:  NEW ISWRO                                            
C          ((VD(I,J),I=1,IPQ),J=1,IPQ):  NEW INVERSE OF HESSIAN         
C                                                                       
C     PARAMETERS:                                                       
C          IPRNT=0:  NOT TO PRINT OUT INTERMEDIATE RESULTS              
C          IPRNT=1:  TO PRINT OUT INTERMEDIATE RESULTS                  
C-----------------------------------------------------------------------
C                                                                       
      IMPLICIT REAL*8(A-H,O-Z)                                          
      REAL*8 MAXVD                                                      
cc      DIMENSION VD(50,50),X(50),G(50),SX(50),SG(50),SR(50),C(50)        
cc      DIMENSION Y(2000),SSX(50)                                         
cc      COMMON /COM50/VD                                                  
      DIMENSION G(IP+IQ),C(IP+IQ),Y(N),X(IP+IQ),VD(IP+IQ,IP+IQ)
      DIMENSION SX(IP+IQ),SG(IP+IQ),SR(IP+IQ),SSX(IP+IQ)                                         
      DIMENSION ALPH(IP+IQ)
      DATA CST0,CST1,CST4,CST10/0.0D-00,1.0D-00,4.0D-00,10.0D-00/       
      DATA CONSTA,CONSTB,EPS3,EPS4/0.5D-00,2.0D-00,0.1D-05,0.1D-05/     
cc      DATA SSX/50*0.0D-00/                                              
C                                                                       
C                                                                       
      ISPHAI=0                                                          
      IPHAI=1                                                           
      IPQ=IQ+IP                                                         
      IPQ2=IPQ+IPQ                                                      
C                                                                       
      ITNS=0                                                            
  150 CONTINUE                                                          
      ITNS=ITNS+1                                                       
      ITN=0                                                             
C                                                                       
 1210 CONTINUE                                                          
      ICOND=0                                                           
C-----  SX=X-C  -----                                                   
      DO 210 I=1,IPQ                                                    
  210 SX(I)=X(I)-C(I)                                                   
C                                                                       
      IF(IPRNT.EQ.0) GO TO 3200                                         
cc      WRITE(6,3900)                                                     
cc      WRITE(6,3910) (C(I),I=1,IPQ)                                      
      WRITE(LU,3900)                                                     
      WRITE(LU,3910) (C(I),I=1,IPQ)                                      
 3200 CONTINUE                                                          
C                                                                       
C                                                                       
C-----  GRADIENT COMPUTATION  -----                                     
C                                                                       
 4000 CONTINUE                                                          
      ICOND=0                                                           
      IF(IQ.LE.0) GO TO 4510                                            
      DO 4500 I=1,IQ                                                    
 4500 SSX(I)=SX(I)                                                      
cc      CALL ARCHCK(SSX,IQ,ICOND)                                         
      CALL ARCHCK(SSX,ALPH,IQ,ICOND)
 4510 IF(IP.LE.0) GO TO 4600                                            
      DO 4520 I=1,IP                                                    
      II=IQ+I                                                           
 4520 SSX(I)=SX(II)                                                     
cc      CALL ARCHCK(SSX,IP,ICOND)                                         
      CALL ARCHCK(SSX,ALPH,IP,ICOND)
 4600 CONTINUE                                                          
C                                                                       
      IF(ICOND.EQ.0) GO TO 309                                          
C                                                                       
      IF(IPRNT.EQ.0) GO TO 1220                                         
cc      WRITE(6,4700)                                                     
      WRITE(LU,4700)                                                     
 1220 CONTINUE                                                          
C                                                                       
      ITN=ITN+1                                                         
C                                                                       
C------------------------------                                         
      MAXVD=CST0                                                        
      DO 4900 I=1,IPQ                                                   
      IF(VD(I,I).GT.MAXVD)MAXVD=VD(I,I)                                 
 4900 CONTINUE                                                          
      VDN=MAXVD/CST4                                                    
      DO 304 I=1,IPQ                                                    
      DO 303 J=1,IPQ                                                    
  303 VD(I,J)=VD(I,J)/CST10                                             
  304 VD(I,I)=VD(I,I)+VDN                                               
      DO 306 I=1,IPQ                                                    
      SCI=CST0                                                          
      DO 305 J=1,IPQ                                                    
  305 SCI=SCI+VD(I,J)*G(J)                                              
  306 C(I)=SCI                                                          
C------------------------------                                         
C                                                                       
      GO TO 1210                                                        
  309 CONTINUE                                                          
cc      CALL SGRAD(SPHAI,SD,SG,Y,N,SX,IQ,IP,IPRNT)                        
      CALL SGRAD(SPHAI,SD,SG,Y,N,SX,IQ,IP,IPRNT,LU)                        
      IF(ICOND.EQ.1) GO TO 1220                                         
      IF(ITN.GE.10) GO TO 312                                           
  312 CONTINUE                                                          
C                                                                       
C-----  SR=V*SG  -----                                                  
      DO 310 I=1,IPQ                                                    
      SUM=CST0                                                          
      DO 311 J=1,IPQ                                                    
  311 SUM=SUM+VD(I,J)*SG(J)                                             
  310 SR(I)=SUM                                                         
C                                                                       
C-----  SRO=(SG)'*(SR)  -----                                           
      SRO=0.0D-00                                                       
      DO 1050 I=1,IPQ                                                   
 1050 SRO=SRO+SG(I)*SR(I)                                               
C                                                                       
C-----  DGAM=-G'*(SR)/SRO  -----                                        
      GSR=0.0D-00                                                       
      DO 1060 I=1,IPQ                                                   
 1060 GSR=GSR+G(I)*SR(I)                                                
      DGAM=-GSR/SRO                                                     
      DGAM1=DGAM+CST1                                                   
      DGAM1=DABS(DGAM1)+0.1D-70                                         
      RAM=DABS(DGAM)/DGAM1                                              
C                                                                       
      IF(RAM.GT.CONSTA) GO TO 430                                       
      RAM=CONSTA                                                        
      IRAM=1                                                            
      GO TO 470                                                         
C                                                                       
  430 IF(RAM.LT.CONSTB) GO TO 450                                       
      RAM=CONSTB                                                        
      IRAM=-1                                                           
      GO TO 470                                                         
C                                                                       
  450 CONTINUE                                                          
      IRAM=0                                                            
  470 RAMSRO=(RAM-CST1)/SRO                                             
      DO 480 I=1,IPQ                                                    
      RAMT=RAMSRO*SR(I)                                                 
      DO 480 J=1,IPQ                                                    
  480 VD(I,J)=VD(I,J)+RAMT*SR(J)                                        
      IF(PHAI.GE.SPHAI) GO TO 540                                       
C                                                                       
C-----  SPHAI.GT.PHAI: TEST OF CORRECTION  -----                        
      RAM1=RAM-CST1                                                     
      IF(DABS(RAM1).LT.EPS3) GO TO 555                                  
      CONSDR=DGAM*RAM1                                                  
      DO 550 I=1,IPQ                                                    
  550 C(I)=C(I)-CONSDR*SR(I)                                            
      IPHAI=0                                                           
      IF(SRO.GT.EPS4) GO TO 900                                         
C     END OF ITERATION                                                  
  555 ISWRO=ISWRO+1                                                     
      GO TO 1000                                                        
C                                                                       
C-----  SHPAI.LE.PHAI: SUCCESSFUL REDUCTION  -----                      
  540 DO 560 I=1,IPQ                                                    
      X(I)=SX(I)                                                        
      G(I)=SG(I)                                                        
  560 C(I)=RAM*SR(I)                                                    
      PHAI=SPHAI                                                        
      SIGMA2=SD                                                         
C                                                                       
      IF(IPRNT.EQ.0) GO TO 571                                          
cc      WRITE(6,570) PHAI                                                 
      WRITE(LU,570) PHAI                                                 
  571 CONTINUE                                                          
C                                                                       
      IPHAI=1                                                           
  800 CONTINUE                                                          
      IF(IRAM.NE.0) GO TO 901                                           
      IF(SRO.LT.EPS4) GO TO 555                                         
C     ITERATION CHECK                                                   
  900 CONTINUE                                                          
      ISPHAI=(ISPHAI+(1-IPHAI))*(1-IPHAI)                               
C                                                                       
      IF(ISPHAI.GT.IPQ2) GO TO 555                                      
C                                                                       
      GO TO 150                                                         
  901 IF(SRO.LT.EPS4) GO TO 555                                         
      GO TO 900                                                         
C     END OF MINIMIZATION                                               
 1000 CONTINUE                                                          
C                                                                       
 1001 RETURN                                                            
  570 FORMAT(1H ,'NEW PHAI=',D20.10)                                    
 3900 FORMAT(//1H ,'C(I)')                                              
 3910 FORMAT(1H ,5D20.10)                                               
 4700 FORMAT(1H ,'ON THE BOUNDARY')                                     
      END                                                               
cc      SUBROUTINE SDATPR(TITL,Y,N,P0,IQ,IP)                              
      SUBROUTINE SDATPR(YS,Y,N,P0,IQ,IP)                              
C                                                                       
C-----------------------------------------------------------------------
C     THIS SUBROUTINE READS IN AND PRINTS OUT INITIAL CONDITION AND DELE
C     THE MEAN OF THE DATA.                                             
C                                                                       
C     THE FOLLOWING INPUTS ARE REQUIRED:                                
C          IQ:  AR-ORDER                                                
C          (B(I),I=1,IQ):  INITIAL ESTIMATES OF AR COEFFICIENTS         
C          IP:  MA-ORDER                                                
C          (A(I),I=1,IP):  INITIAL ESTIMATES OF MA COEFFICIENTS         
C          TITL:  TITLE OF THE DATA                                     
C          N:  DATA LENGTH                                              
C          (DFORM(I),I=1,20):  INPUT FORMAT SPECIFICATION IN ONE CARD,  
C                              EXAMPLE,                                 
C                              (8F10.4)                                 
C          (Y(I),I=1,N):  ORIGINAL DATA                                 
C                                                                       
C     THE AR-MA MODEL IS GIVEN BY                                       
C     Y(I)+B(1)Y(I-1)+...+B(IQ)Y(I-IQ)=X(I)+A(1)X(I-1)+...+A(IP)X(I-IP),
C     WHERE X(I) IS A ZERO MEAN WHITE NOISE.                            
C-----------------------------------------------------------------------
C                                                                       
      IMPLICIT REAL*8(A-H,O-Z)                                          
cc      REAL*4 DFORM(20),TITL(20)                                         
cc      DIMENSION Y(2000),A(50),B(50),P0(50)                              
cc      DATA A,B/100*0.0D-00/                                             
      DIMENSION YS(N), Y(N), P0(IP+IQ)
      DATA CST0/0.0D-00/                                                
C                                                                       
C-----  INITIAL CONDITION LOADING FOR AR-MA(IQ,IP)  -----               
cc      READ(5,1) IQ                                                      
cc      IF(IQ.LE.0) GO TO 4205                                            
cc      READ(5,2) (B(I),I=1,IQ)                                           
cc 4205 READ(5,1) IP                                                      
cc      IF(IP.LE.0) GO TO 4206                                            
cc      READ(5,2) (A(I),I=1,IP)                                           
cc 4206 CONTINUE                                                          
C                                                                       
C-----  P0 ARRANGEMENT  -----                                           
cc      IF(IQ.LE.0) GO TO 300                                             
cc      DO 200 I=1,IQ                                                     
cc      P0(I)=B(I)                                                        
cc  200 CONTINUE                                                          
cc  300 IF(IP.LE.0) GO TO 310                                             
cc      DO 210 I=1,IP                                                     
cc      II=IQ+I                                                           
cc      P0(II)=A(I)                                                       
cc  210 CONTINUE                                                          
cc  310 CONTINUE                                                          
C                                                                       
C-----  DATA INPUT  -----                                               
cc      READ(5,4) (TITL(I),I=1,20)                                        
cc      READ(5,1) N                                                       
cc      READ(5,4) (DFORM(I),I=1,20)                                       
C     ORIGINAL DATA INPUT AND PRINT OUT                                 
cc      READ(5,DFORM) (Y(I),I=1,N)                                        
      DO 320 I=1,N
  320 Y(I)=YS(I)
C                                                                       
cc      WRITE(6,59)                                                       
cc      WRITE(6,162) (TITL(I),I=1,20)                                     
cc      WRITE(6,62) IQ                                                    
cc      IF(IQ.LE.0) GO TO 4215                                            
cc      WRITE(6,65) (B(I),I=1,IQ)                                         
cc 4215 CONTINUE                                                          
cc      WRITE(6,63) IP                                                    
cc      IF(IP.LE.0) GO TO 4216                                            
cc      WRITE(6,65) (A(I),I=1,IP)                                         
cc 4216 CONTINUE                                                          
C                                                                       
cc      WRITE(6,285)                                                      
cc      WRITE(6,164) N                                                    
cc      WRITE(6,610) (Y(I),I=1,N)                                         
C                                                                       
C-----  MEAN DELETION  -----                                            
      AN=N                                                              
      SUM=CST0                                                          
      DO 9 I=1,N                                                        
    9 SUM=SUM+Y(I)                                                      
      YMEAN=SUM/AN                                                      
      DO 10 I=1,N                                                       
   10 Y(I)=Y(I)-YMEAN                                                   
C                                                                       
      RETURN                                                            
    1 FORMAT(16I5)                                                      
    2 FORMAT(4D20.10)                                                   
    4 FORMAT(20A4)                                                      
   59 FORMAT( //,' PROGRAM TIMSAC 78.5.2',/,
     *'   EXACT MAXIMUM LIKELIHOOD METHOD OF AR-MA MODEL FITTING;',
     *'  SCALAR CASE',/,'   < AR-MA MODEL >',
     *//,11X,'Y(I) + B(1)*Y(I-1) + ... + B(IQ)*Y(I-IQ)  =  ',
     *'X(I) + A(1)*X(I-1) + ... + A(IP)*X(I-IP)',/,' WHERE',/,11X,
     *'IQ:    AR-ORDER',/,11X,'IP:    MA-ORDER',/,11X,
     *'X(I): ZERO MEAN WHITE NOISE' ) 
   62 FORMAT(//1H ,5('-'),2X,'INITIAL AR(I)',2X,'IQ=',I3,2X,5('-'))     
   63 FORMAT(//1H ,5('-'),2X,'INITIAL MA(I)',2X,'IP=',I3,2X,5('-'))     
   65 FORMAT(1H ,5D20.10,/(1H ,5D20.10))                                
  162 FORMAT( 1H ,' TITLE:',/,2X,20A4 )                                 
  164 FORMAT(1H ,'N=',I5)                                               
  285 FORMAT(///1H ,5('-'),2X,'ORIGINAL DATA',2X,5('-'))                
  610 FORMAT(1H ,10F10.5,/(1H ,10F10.5))                                
      END                                                               
cc      SUBROUTINE SGRAD(F0,SD,G,Y,N,P0,IQ,IP,IPRNT)                      
      SUBROUTINE SGRAD(F0,SD,G,Y,N,P0,IQ,IP,IPRNT,LU)                      
C                                                                       
C-----------------------------------------------------------------------
C     THIS SUBROUTINE COMPUTES AN APPROXIMATION TO GRADIENT BY DIFFERENC
C     THIS SUBROUTINE SHOULD EVENTUALLY BE REPLACED BY AN ANALYTIC EVALU
C     PROCEDURE OF GRADIENTS.                                           
C       ----------------------------------------------------------------
C       THE FOLLOWING SUBROUTINES ARE DIRECTLY CALLED BY THIS SUBROUTINE
C             FUNCT2                                                    
C             ARCHCK                                                    
C       ----------------------------------------------------------------
C                                                                       
C     INPUTS:                                                           
C          N:  DATA LENGTH                                              
C          (Y(I),I=1,N):  ORIGINAL DATA, MEAN DELETED                   
C          IQ:  AR-ORDER                                                
C          IP:  MA-ORDER                                                
C          (P0(I),I=1,IPQ):  THE VECTOR OF AR AND MA COEFFICIENTS (IPQ=I
C                                                                       
C     OUTPUTS:                                                          
C          F0:  (-2)LOG LIKELIHOOD                                      
C          SD:  WHITE NOISE VARIANCE                                    
C          (G(I),I=1,IPQ):  GRADIENT                                    
C                                                                       
C     PARAMETERS:                                                       
C          EPSA:  ORDINATE DIFFERENCE FOR GRADIENT COMPUTATION BY DIFFER
C          IPRNT=0:  NOT TO PRINT OUT INTERMEDIATE RESULTS              
C          IPRNT=1:  TO PRINT OUT INTERMEDIATE RESULTS                  
C-----------------------------------------------------------------------
C                                                                       
      IMPLICIT REAL*8(A-H,O-Z)                                          
cc      DIMENSION P0(50),P1(50),Y(2000),G(50),PP0(50)                     
      DIMENSION P0(IP+IQ),P1(IP+IQ),Y(N),G(IP+IQ),PP0(IP+IQ)                     
      DIMENSION ALPH(IP+IQ)
      DATA EPSA,CST07/0.1D-03,0.7D-00/                                  
C                                                                       
cc      CALL FUNCT2(F0,SD,Y,N,P0,IQ,IP)                                   
      IP1=IP+1                                                          
      IR=MAX0(IQ,IP1)                                                   
      CALL FUNCT2(F0,SD,Y,N,P0,IQ,IP,IR)                                   
      IPQ=IP+IQ                                                         
      DO 9 J=1,IPQ                                                      
      P1(J)=P0(J)                                                       
    9 CONTINUE                                                          
C                                                                       
C-----  GRADIENT COMPUTATION  -----                                     
      DO 10 I=1,IPQ                                                     
      EPSAS=EPSA                                                        
      ITR=1                                                             
 4000 CONTINUE                                                          
      ICOND=0                                                           
      P1(I)=P0(I)+EPSAS                                                 
      IF(I.GT.IQ) GO TO 4500                                            
      DO 4010 II=1,IQ                                                   
      PP0(II)=P1(II)                                                    
 4010 CONTINUE                                                          
cc      CALL ARCHCK(PP0,IQ,ICOND)                                         
      CALL ARCHCK(PP0,ALPH,IQ,ICOND)
      GO TO 5000                                                        
 4500 CONTINUE                                                          
      DO 5010 II=1,IP                                                   
      III=IQ+II                                                         
 5010 PP0(II)=P1(III)                                                   
cc      CALL ARCHCK(PP0,IP,ICOND)                                         
      CALL ARCHCK(PP0,ALPH,IP,ICOND)
 5000 CONTINUE                                                          
      IF(ICOND.EQ.0) GO TO 4100                                         
      IF(ITR.LT.10) GO TO 4110                                          
C     WRITE(6,4200)                                                     
      RETURN                                                            
 4110 EPSAS=EPSAS*CST07                                                 
      EPSAS=-EPSAS                                                      
      ITR=ITR+1                                                         
      GO TO 4000                                                        
 4100 CONTINUE                                                          
cc      CALL FUNCT2(F1,SDN,Y,N,P1,IQ,IP)                                  
      CALL FUNCT2(F1,SDN,Y,N,P1,IQ,IP,IR)
      G(I)=(F1-F0)/EPSAS                                                
      P1(I)=P0(I)                                                       
   10 CONTINUE                                                          
C--------------------                                                   
C                                                                       
      IF(IPRNT.EQ.0) RETURN                                             
cc      WRITE(6,450)                                                      
cc      WRITE(6,520) (P0(I),I=1,IPQ)                                      
cc      WRITE(6,500)                                                      
cc      WRITE(6,520) (G(I),I=1,IPQ)                                       
      WRITE(LU,450)                                                      
      WRITE(LU,520) (P0(I),I=1,IPQ)                                      
      WRITE(LU,500)                                                      
      WRITE(LU,520) (G(I),I=1,IPQ)                                       
C                                                                       
      RETURN                                                            
  450 FORMAT(1H ,'P0(I)')                                               
  500 FORMAT(1H ,'GRADIENT G(I)')                                       
  520 FORMAT(1H ,5D20.10)                                               
 4200 FORMAT(1H ,'ICOND=1; GRADIENT UNOBTAINED')                        
      END                                                               
cc      SUBROUTINE SMINOP(TL,SIGMA2,Y,N,P0,IQ,IP)                         
      SUBROUTINE SMINOP( TL,TL2,SIGMA2,Y,N,P0,G,P02,G2,ALPHB,ALPHA,IQ,
     *                   IP,IPRNT,LU )
C                                                                       
C-----------------------------------------------------------------------
C     THIS SUBROUTINE CONTROLS THE MAXIMUM LIKELIHOOD COMPUTATION.      
C                                                                       
C       THE FOLLOWING SUBROUTINES ARE DIRECTLY CALLED BY THIS SUBROUTINE
C             SGRAD                                                     
C             ARCHCK                                                    
C                                                                       
C     INPUTS:                                                           
C          N:  DATA LENGTH                                              
C          (Y(I),I=1,N):  ORIGINAL DATA, MEAN DELETED                   
C          IQ:  AR-ORDER                                                
C          IP:  MA-ORDER                                                
C          (P0(I),I=1,IPQ):  THE VECTOR OF AR AND MA COEFFICIENTS (IPQ=I
C                                                                       
C     OUTPUTS:                                                          
C          TL:  (-2)LOG LIKELIHOOD                                      
C          SIGMA2:  WHITE NOISE VARIANCE                                
C                                                                       
C     PARAMETERS:                                                       
C          IPRNT=0:  NOT TO PRINT OUT INTERMEDIATE RESULTS              
C          IPRNT=1:  TO PRINT OUT INTERMEDIATE RESULTS                  
C-----------------------------------------------------------------------
C                                                                       
      IMPLICIT REAL*8(A-H,O-Z)                                          
      REAL*8 MAXAB                                                      
cc      DIMENSION Y(2000),P0(50),G(50),HS(50,50),CR(50),PP0(50),ALPH(50)  
      DIMENSION Y(N),P0(IP+IQ),P02(IP+IQ)
      DIMENSION G(IP+IQ),G2(IP+IQ),HS(IP+IQ,IP+IQ),CR(IP+IQ)
      DIMENSION PP0(IP+IQ),ALPH(IP+IQ),ALPHB(IQ),ALPHA(IP)
cc      COMMON /COM50/HS /CMALPH/ALPH                                     
cc      DATA PP0/50*0.0D-00/                                              
      DATA CST0,CST10,CST05,CST005/0.0D-00,10.0D-00,0.1D-03,0.00005D-00/
C                                                                       
C                                                                       
cc      IPRNT=0                                                           
C                                                                       
      IPQ=IP+IQ                                                         
      DO 310 I=1,IPQ                                                    
      G(I)=CST0                                                         
      PP0(I)=CST0
      DO 310 J=1,IPQ                                                    
  310 HS(I,J)=CST0                                                      
C                                                                       
C                                                                       
C-----  INITIAL GRADIENT COMPUTATION  -----                             
 4000 CONTINUE                                                          
      ICOND=0                                                           
      IF(IQ.LE.0) GO TO 4510                                            
      DO 4500 I=1,IQ                                                    
 4500 PP0(I)=P0(I)                                                      
cc      CALL ARCHCK(PP0,IQ,ICOND)                                         
      CALL ARCHCK(PP0,ALPH,IQ,ICOND)
      DO 5000 I=1,IQ                                                    
 5000 P0(I)=PP0(I)                                                      
 4510 IF(IP.LE.0) GO TO 4800                                            
      DO 4700 I=1,IP                                                    
      II=IQ+I                                                           
 4700 PP0(I)=P0(II)                                                     
cc      CALL ARCHCK(PP0,IP,ICOND)                                         
      CALL ARCHCK(PP0,ALPH,IP,ICOND)
      DO 5100 I=1,IP                                                    
      II=IQ+I                                                           
      P0(II)=PP0(I)                                                     
 5100 CONTINUE                                                          
 4800 CONTINUE                                                          
      ISWRO=0                                                           
cc      CALL SGRAD(TL,SIGMA2,G,Y,N,P0,IQ,IP,IPRNT)                        
      CALL SGRAD(TL,SIGMA2,G,Y,N,P0,IQ,IP,IPRNT,LU)                        
C                                                                       
cc      WRITE(6,450)                                                      
cc      WRITE(6,520) (P0(I),I=1,IPQ)                                      
cc      WRITE(6,500)                                                      
cc      WRITE(6,520) (G(I),I=1,IPQ)                                       
cc      WRITE(6,550) TL                                                   
      DO 4850 I=1,IPQ
         P02(I)=P0(I)
         G2(I)=G(I)
 4850 CONTINUE
      TL2=TL
C                                                                       
 4890 CONTINUE                                                          
      MAXAB=CST0                                                        
      DO 4900 I=1,IPQ                                                   
cc      PAB=DABS(G(I))                                                    
      PAB=DABS(G2(I))
      IF(PAB.GT.MAXAB) MAXAB=PAB                                        
 4900 CONTINUE                                                          
C                                                                       
C                                                                       
C-----  INVERSE OF HESSIAN COMPUTATION  -----                           
      BN=CST05/MAXAB                                                    
      DO 3010 I=1,IPQ                                                   
      DO 3009 J=1,IPQ                                                   
 3009 HS(I,J)=HS(I,J)/CST10                                             
      HS(I,I)=BN+HS(I,I)                                                
 3010 CONTINUE                                                          
C                                                                       
      IF(IPRNT.EQ.0) GO TO 3120                                         
cc      WRITE(6,3000)                                                     
      WRITE(LU,3000)                                                     
      DO 3100 I=1,IPQ                                                   
cc 3100 WRITE(6,3110) I,(HS(I,J),J=1,IPQ)                                 
 3100 WRITE(LU,3110) I,(HS(I,J),J=1,IPQ)                                 
 3120 CONTINUE                                                          
C                                                                       
C-----  CORRECTION TERM CR(X)=HS*G(X) COMPUTATION  -----                
      DO 900 I=1,IPQ                                                    
      SUM=CST0                                                          
      DO 910 J=1,IPQ                                                    
cc  910 SUM=SUM+HS(I,J)*G(J)                                              
  910 SUM=SUM+HS(I,J)*G2(J)
  900 CR(I)=SUM                                                         
C                                                                       
      IF(IPRNT.EQ.0) GO TO 3920                                         
cc      WRITE(6,3900)                                                     
cc      WRITE(6,3910) (CR(I),I=1,IPQ)                                     
      WRITE(LU,3900)                                                     
      WRITE(LU,3910) (CR(I),I=1,IPQ)                                     
 3920 CONTINUE                                                          
C                                                                       
C-----  DAVIDON'S PROCEDURE  -----                                      
C     MINIMIZATION OF INNOVATION VARIANCE                               
cc      CALL MSDAV2(TL,SIGMA2,G,CR,Y,N,P0,IQ,IP,ISWRO,IPRNT)              
      CALL MSDAV2(TL2,SIGMA2,G2,CR,Y,N,P02,IQ,IP,ISWRO,HS,IPRNT,LU)
C                                                                       
C------------------------------                                         
      IF(IPRNT.EQ.0) GO TO 3930                                         
cc      WRITE(6,1200) ISWRO                                               
      WRITE(LU,1200) ISWRO                                               
 3930 CONTINUE                                                          
C                                                                       
      IF(ISWRO.GE.IPQ) GO TO 1201                                       
      DO 1902 I=1,IPQ                                                   
cc      IF(DABS(PP0(I)-P0(I)).GE.CST005) GO TO 1919                       
      IF(DABS(PP0(I)-P02(I)).GE.CST005) GO TO 1919
 1902 CONTINUE                                                          
      GO TO 1201                                                        
 1919 CONTINUE                                                          
C                                                                       
      IF(IPRNT.EQ.0) GO TO 3950                                         
cc      WRITE(6,1926)                                                     
      WRITE(LU,1926)                                                     
 3950 CONTINUE                                                          
C                                                                       
      GO TO 4890                                                        
C                                                                       
 1201 CONTINUE                                                          
C------------------------------                                         
C                                                                       
cc      WRITE(6,650)                                                      
cc      WRITE(6,520) (P0(I),I=1,IPQ)                                      
cc      WRITE(6,660)                                                      
cc      WRITE(6,520) (G(I),I=1,IPQ)                                       
C                                                                       
C-----  ALPH(I) PRINT OUT  -----                                        
cc      WRITE(6,6100)                                                     
      ICOND=0                                                           
      IF(IQ.LE.0) GO TO 6510                                            
      DO 6500 I=1,IQ                                                    
cc 6500 PP0(I)=P0(I)                                                      
cc      CALL ARCHCK(PP0,IQ,ICOND)                                         
 6500 PP0(I)=P02(I)                                                      
      CALL ARCHCK(PP0,ALPHB,IQ,ICOND)
cc      WRITE(6,6501)                                                     
cc      WRITE(6,520) (ALPH(I),I=1,IQ)                                     
 6510 IF(IP.LE.0) GO TO 6800                                            
      DO 6700 I=1,IP                                                    
      II=IQ+I                                                           
cc 6700 PP0(I)=P0(II)                                                     
cc      CALL ARCHCK(PP0,IP,ICOND)                                         
 6700 PP0(I)=P02(II)
      CALL ARCHCK(PP0,ALPHA,IP,ICOND)
cc      WRITE(6,6502)                                                     
cc      WRITE(6,520) (ALPH(I),I=1,IP)                                     
 6800 CONTINUE                                                          
C------------------------------                                         
C                                                                       
cc      WRITE(6,670) TL                                                   
C                                                                       
C                                                                       
      RETURN                                                            
  450 FORMAT(//1H ,5('-'),2X,'INITIAL VALUES P0(I)',2X,5('-'))          
  500 FORMAT(//1H ,5('-'),2X,'INITIAL GRADIENT G(I)',2X,5('-'))         
  520 FORMAT(1H ,5D20.10,/(1H ,5D20.10))                                
  550 FORMAT(//1H ,'INITIAL (-2)LOG LIKELIHOOD=',D20.10,///)            
  650 FORMAT(////1H ,5('-'),2X,'FINAL VALUES P0(I)',2X,5('-'))          
  660 FORMAT(//1H ,5('-'),2X,'FINAL GRADIENT G(I)',2X,5('-'))           
  670 FORMAT(//1H ,'FINAL (-2)LOG LIKELIHOOD=',D20.10,///)              
 1200 FORMAT(//1H ,'ISWRO=',I5)                                         
 1926 FORMAT(1H ,'HESSIAN RESET')                                       
 3000 FORMAT(1H ,'INVERSE OF HESSIAN')                                  
 3110 FORMAT(1H ,I5,4X,10D12.5,/(1H ,8X,10D12.5))                       
 3900 FORMAT(1H ,'CR(I)')                                               
 3910 FORMAT(1H ,5D20.10)                                               
 6100 FORMAT(//1H ,5('-'),2X,'FINAL ALPH(I)''S AT SUBROUTINE ARCHCK',2X,
     A5('-'))                                                           
 6501 FORMAT(/1H ,12X,'AR-PART')                                        
 6502 FORMAT(/1H ,12X,'MA-PART')                                        
      END                                                               
      SUBROUTINE SUBPM(P,B,A,IQ,IP,IR)                                  
C                                                                       
C-----------------------------------------------------------------------
C     THIS SUBROUTINE COMPUTES THE VARIANCE MATRIX OF A STATIONARY STATE
C     VECTOR BY THE PROCEDURE OF AKAIKE (RESEARCH MEMO. 139 INST. STATIS
C     MATH. OCTOBER, 1978).                                             
C                                                                       
C     INPUTS:                                                           
C          IQ:  AR-ORDER                                                
C          (B(I),I=1,IQ):  AR-COEFFICIENTS                              
C          IP:  MA-ORDER                                                
C          (A(I),I=1,IP):  MA-COEFFICIENTS                              
C          IR:  IR=MAX(IQ,IP+1)                                         
C                                                                       
C     OUTPUTS:                                                          
C          ((P(I,J),I=1,IR),J=1,IR):  VARIANCE MATRIX OF THE STATIONARY 
C                                     STATE VECTOR                      
C-----------------------------------------------------------------------
C                                                                       
      IMPLICIT REAL*8(A-H,O-Z)                                          
cc      DIMENSION P(51,51),A(50),B(50),WS(51),R(51),DB(1300)              
      DIMENSION P(IR,IR),A(IR),B(IR),WS(IR),R(IR+1),DB(IQ*2)              
      DATA CST0,CST1/0.0D-00,1.0D-00/                                   
C                                                                       
C                                                                       
C                                                                       
cc      DO 6100 I=1,1300                                                  
      DO 6100 I=1,IQ*2
      DB(I)=CST0                                                        
 6100 CONTINUE                                                          
cc      DO 6200 I=1,51                                                    
      DO 6200 I=1,IR
      R(I)=CST0                                                         
 6200 CONTINUE                                                          
C                                                                       
C-----  IMPULSE RESPONSE COMPUTATION  -----                             
      WS(1)=CST1                                                        
      IPP1=IP+1                                                         
      IQP1=IQ+1                                                         
      IRM1=IR-1                                                         
      IF(IRM1.LE.0) GO TO 400                                           
      DO 100 I=2,IR                                                     
      IM1=I-1                                                           
      IMPR=MIN0(IQ,IM1)                                                 
      SUM=CST0                                                          
      IF(IMPR.EQ.0) GO TO 201                                           
      DO 200 J=1,IMPR                                                   
      IMJ=I-J                                                           
      SUM=SUM-B(J)*WS(IMJ)                                              
  200 CONTINUE                                                          
  201 IF(I.LE.IPP1) SUM=SUM+A(IM1)                                      
      WS(I)=SUM                                                         
  100 CONTINUE                                                          
  400 CONTINUE                                                          
C                                                                       
C-----  PREPARATION OF CONSTANT VECTOR  -----                           
      IRP1=IR+1                                                         
      R(IRP1)=CST0                                                      
      IF(IRM1.EQ.0) GO TO 7200                                          
      R(IR)=A(IRM1)                                                     
      DO 7010 I=1,IRM1                                                  
      IM1=I-1                                                           
      IPMI=IP-IM1                                                       
      SUM=CST0                                                          
      IF(IPMI.LE.0) GO TO 7021                                          
      DO 7020 J=1,IPMI                                                  
      JI=IM1+J                                                          
      SUM=SUM+A(JI)*WS(J+1)                                             
 7020 CONTINUE                                                          
 7021 IF(I.EQ.1) SUM=SUM+CST1                                           
      IF(I.GT.1) SUM=SUM+A(IM1)                                         
      R(I)=SUM                                                          
 7010 CONTINUE                                                          
      GO TO 7300                                                        
 7200 CONTINUE                                                          
      R(1)=CST1                                                         
 7300 CONTINUE                                                          
C                                                                       
C                                                                       
      IF(IQ.EQ.0) GO TO 1600                                            
C                                                                       
C-----  TRIANGULARIZATION OF THE COEFFICIENT MATRIX  -----              
      IQP1=IQ+1                                                         
      DO 2100 L=1,IQ                                                    
 2100 DB(L)=B(L)                                                        
      IT=IQ                                                             
      I=IQ                                                              
 2200 CONTINUE                                                          
      IP1=I+1                                                           
      IP2=IP1+1                                                         
      IH=IP2/2                                                          
      IHP1=IH+1                                                         
      CKB=DB(IT)                                                        
      CKI=CST1/(CST1-CKB**2)                                            
      DO 2300 K=1,IH                                                    
      IPK=IP2-K                                                         
      R(K)=(R(K)-CKB*R(IPK))*CKI                                        
 2300 CONTINUE                                                          
C                                                                       
      IF(I.LE.2) GO TO 2401                                             
      DO 2400 K=IHP1,I                                                  
      IPK=IP2-K                                                         
      R(K)=R(K)-CKB*R(IPK)                                              
 2400 CONTINUE                                                          
C                                                                       
 2401 IM1=I-1                                                           
 2410 IF(IM1.EQ.0) GO TO 2600                                           
      DO 2500 K=1,IM1                                                   
      ITPK=IT+K                                                         
      IMIK=IT-I+K                                                       
      ITMK=IT-K                                                         
      DB(ITPK)=(DB(IMIK)-CKB*DB(ITMK))*CKI                              
 2500 CONTINUE                                                          
C                                                                       
      I=IM1                                                             
      IT=IT+I                                                           
      GO TO 2200                                                        
 2600 CONTINUE                                                          
C                                                                       
C-----  SOLVING THE LINEAR EQUATION  -----                              
      IF(IQ.LE.1) GO TO 3110                                            
      ITP1=IT+1                                                         
      DO 3100 I=2,IQ                                                    
      SUM=R(I)                                                          
      IM1=I-1                                                           
      DO 3200 J=1,IM1                                                   
      ITP1=ITP1-1                                                       
      SUM=SUM-DB(ITP1)*R(J)                                             
 3200 CONTINUE                                                          
      R(I)=SUM                                                          
 3100 CONTINUE                                                          
 3110 CONTINUE                                                          
      DO 3300 I=IQP1,IRP1                                               
      SUM=R(I)                                                          
      DO 3400 J=1,IQ                                                    
      IMJ=I-J                                                           
      SUM=SUM-DB(J)*R(IMJ)                                              
 3400 CONTINUE                                                          
      R(I)=SUM                                                          
 3300 CONTINUE                                                          
C                                                                       
 1600 CONTINUE                                                          
C                                                                       
C-----  P(I,J) COMPUTATION  -----                                       
      DO 9100 I=1,IR                                                    
      DO 9200 J=1,I                                                     
      SUM=CST0                                                          
      IF(J.EQ.1) GO TO 9400                                             
      DO 9300 K=1,J                                                     
      IJK=K+I-J                                                         
      SUM=SUM+WS(IJK)*WS(K)                                             
 9300 CONTINUE                                                          
 9400 CONTINUE                                                          
      IJ1=I-J+1                                                         
      P(I,J)=R(IJ1)-SUM                                                 
      P(J,I)=P(I,J)                                                     
 9200 CONTINUE                                                          
 9100 CONTINUE                                                          
      RETURN                                                            
      END                                                               