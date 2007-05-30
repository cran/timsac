      SUBROUTINE BSUBSTF( ZS,N,IMODEL,LAG,K,IL,LG1,LG2,F,CNST,ZMEAN,SUM,
     *M,AICM,SDM,A1,SD,AIC,DIC,AICB,SDB,EK,A2,IND,C,C1,C2,B,OEIC,ESUM,
     *OMEAN,OM,E,EMEAN,VARI,SKEW,PEAK,COV,SXX )
C
      INCLUDE 'timsac_f.h'
C
cc      PROGRAM BSUBST
C.......................................................................
C.....PLANNED BY H.AKAIKE...............................................
C.....DESIGNED BY H.AKAIKE AND G.KITAGAWA...............................
C.....PROGRAMMED BY G.KITAGAWA AND F.TADA...............................
C.....ADDRESS: THE INSTITUTE OF STATISTICAL MATHEMATICS, 4-6-7 MINAMI-AZ
C..............MINATO-KU, TOKYO 106, JAPAN..............................
C.....DATE OF THE LATEST REVISION:  MAY  14, 1979.......................
C.......................................................................
C.....THIS PROGRAM WAS ORIGINALLY PUBLISHED IN "TIMSAC-78", BY H.AKAIKE,
C.....G.KITAGAWA, E.ARAHATA AND F.TADA, COMPUTER SCIENCE MONOGRAPHS, NO.
C.....THE INSTITUTE OF STATISTICAL MATHEMATICS, TOKYO, 1979.............
C.......................................................................
C     TIMSAC 78.1.3.                                                    
C     _                 ____ _                                          
C     BAYESIAN TYPE ALL SUBSET ANALYSIS OF TIME SERIES BY A MODEL LINEAR
C     IN PARAMETERS                                                     
C                                                                       
C     THIS PROGRAM PRODUCES BAYESIAN ESTIMATES OF TIME SERIES MODELS SUC
C     PURE AR MODELS, AR-MODELS WITH NON-LINEAR TERMS, AR-MODELS WITH PO
C     TYPE MEAN VALUE FUNCTIONS, ETC.  THE GOODNESS OF FIT OF A MODEL IS
C     CHECKED BY THE ANALYSIS OF SEVERAL STEPS AHEAD PREDICTION ERRORS. 
C     BY PREPARING AN EXTERNAL SUBROUTINE SETX PROPERLY, ANY TIME SERIES
C     WHICH IS LINEAR IN PARAMETERS CAN BE TREATED.                     
C       ----------------------------------------------------------------
C       THE FOLLOWING SUBROUTINES ARE DIRECTLY CALLED BY THIS PROGRAM:  
C             REDATA                                                    
C             REDLAG                                                    
C             SETLAG                                                    
C             REDREG                                                    
C             REDUCT                                                    
C             ARMFIT                                                    
C             SBBAYS                                                    
C             CHECK                                                     
C       ----------------------------------------------------------------
C                                                                       
C       INPUTS REQUIRED:                                                
C          MT:    ORIGINAL DATA INPUT DEVICE SPECIFICATION              
C          IMODEL:=1  AUTOREGRESSIVE MODEL                              
C                 =2  POLYNOMIAL TYPE NON-LINEAR MODEL (LAG'S READ IN ) 
C                 =3  POLYNOMIAL TYPE NON-LINEAR MODEL (LAG'S AUTOMATICA
C                 =4  AR-MODEL WITH POLYNOMIAL MEAN VALUE FUNCTION      
C                 =5  ANY NON-LINEAR MODEL                              
C                 =6  POLYNOMIAL TYPE EXPONENTIALLY DAMPED NON-LINEAR MO
C                 =7  THIS MODEL IS RESERVED FOR THE USER'S OPTIONAL USE
C          LAG:   MAXIMUM TIME LAG USED IN THE MODEL                    
C          K:     NUMBER OF REGRESSORS                                  
C          IL:    PREDICTION ERRORS CHECKING (UP TO IL-STEPS AHEAD) IS R
C                 N*IL SHOULD BE LESS THAN OR EQUAL TO 20000            
C                                                                       
C       --   THE FOLLOWING INPUTS ARE REQUIRED AT SUBROUTINE REDATA   --
C                                                                       
C          TITLE:   ORIGINAL DATA SPECIFICATION                         
C          N:       DATA LENGTH                                         
C          DFORM:   INPUT DATA FORMAT SPECIFICATION STATEMENT           
C                   -- EXAMPLE --  (8F10.5 )                            
C          X(I) (I=1,N):   ORIGINAL DATA                                
C       ----------------------------------------------------------------
C                                                                       
cc      !DEC$ ATTRIBUTES DLLEXPORT :: BSUBSTF
C
      PARAMETER  ( MJ2 = 101 )
C
      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
CC      REAL*4     Z(10000) , TITLE(20) , TTL(5) , E(20000)               
cc      REAL*4     TITLE(20) , TTL(5)
cc      DIMENSION  Z(10000), E(20000) 
cc      DIMENSION  X(200,101), D(200) , A(100) , B(100)                  
cc      DATA  TTL  / 4H   B,4HAYES,4HIAN ,4HMODE,4HL    /                 
      COMMON     / BBB / L1(50) , L2(50) , L3(50) , SD0 , CONST         
cc      REAL * 4   CONST , SD0                                            
      DIMENSION  ZS(N), Z(N), E(N,IL) 
      DIMENSION  X(N,K+1), A1(K), A2(K), B(K)                  
      DIMENSION  LG1(3,K), LG2(5)
      DIMENSION  SD(K+1), AIC(K+1), DIC(K+1)
      DIMENSION  IND(K), C(K), C1(K+1), C2(K), ESUM(K+1)
      DIMENSION  EMEAN(IL), VARI(IL), SKEW(IL), PEAK(IL), COV(MJ2)
      DIMENSION  SXX(121), FA(N-LAG,IL)
cc      CHARACTER*4 F(20,K+1)
      INTEGER*1  F(80,K+1)
      CHARACTER(4) G
      COMMON     / EEE /  G(20,31)                                      
      COMMON     / AAA /  NN
C                                                                       
C        EXTERNAL SUBROUTINE DECLARATION:                               
C                                                                       
      EXTERNAL  SETX1                                                   
      EXTERNAL  SETX2                                                   
      EXTERNAL  SETX4                                                   
      EXTERNAL  SETX5                                                   
      EXTERNAL  SETX6                                                   
      EXTERNAL  SETX7                                                   
      EXTERNAL  PRDCT1                                                  
      EXTERNAL  PRDCT2                                                  
      EXTERNAL  PRDCT3                                                  
      EXTERNAL  PRDCT6                                                  
C
cc	CHARACTER(100) IFLNAM,OFLNAM
cc	CALL FLNAM2( IFLNAM,OFLNAM,NFL )
cc	IF ( NFL.EQ.0 ) GO TO 999
cc	IF ( NFL.EQ.2 ) THEN
cc	   OPEN( 6,FILE=OFLNAM,ERR=900,IOSTAT=IVAR)
cc	ELSE
cc	   CALL SETWND
cc	END IF
C
C
C        PARAMETERS:                                                    
C             MJ1:  ABSOLUTE DIMENSION FOR SUBROUTINE CALL              
C                                                                       
      IF ((IMODEL.LE.0) .OR. (IMODEL.GE.8)) GO TO 150                                                         
cc      MJ = 1000                                                         
cc      MJ1 = 200                                                         
      NN = N
      MJ = N
      MJ1 = N
      ISW = 1                                                           
      IPR = 2                                                           
C                                                                       
CC      READ( 5,1 )     MT                                                
cc      MT = 5
cc      OPEN( MT,FILE=IFLNAM,ERR=910,IOSTAT=IVAR,STATUS='OLD' )
cc      READ( 5,1 )     IMODEL , LAG , K , IL                             
cc      WRITE( 6,3 )                                                      
cc      IF( IMODEL.EQ.1 )  WRITE( 6,4 )   IMODEL                          
cc      IF( IMODEL.EQ.2 )  WRITE( 6,5 )   IMODEL                          
cc      IF( IMODEL.EQ.3 )  WRITE( 6,5 )   IMODEL                          
cc      IF( IMODEL.EQ.5 )  WRITE( 6,5 )   IMODEL                          
cc      IF( IMODEL.EQ.6 )  WRITE( 6,11 )   IMODEL                         
cc      WRITE( 6,6 )                                                      
cc      IF( IMODEL.EQ.1 )  WRITE( 6,7 )  K                                
cc      IF( IMODEL.EQ.2 )  WRITE( 6,8 )   LAG , K                         
cc      WRITE( 6,2 )     MT                                               
C                                                                       
C                                                                       
C          ---------------------------------------                      
C          ORIGINAL DATA LOADING AND MEAN DELETION                      
C          ---------------------------------------                      
C                                                                       
cc      CALL  REDATA( Z,N,MT,TITLE )                                      
      CALL  REDATA( ZS,Z,N,ZMEAN,SUM )
      NMK = N - LAG                                                     
      LAG1 = LAG + 1                                                    
C                                                                       
C          ---------------------                                        
C          HOUSEHOLDER REDUCTION                                        
C          ---------------------                                        
C                                                                       
      GO TO ( 10,20,30,40,50,60,70 ), IMODEL                            
C                                                                       
   10 K = LAG                                                           
cc      CALL  REDUCT( SETX1,Z,D,NMK,0,K,MJ1,LAG,X )                       
      CALL  REDUCT( SETX1,Z,NMK,0,K,MJ1,LAG,X )                       
      GO TO 100                                                         
C                                                                       
cc   20 CALL  REDLAG( K )                                                 
   20 CONTINUE
      DO 21 I=1,K
         L1(I) = LG1(1,I)
         L2(I) = LG1(2,I)
         L3(I) = LG1(3,I)
   21 CONTINUE
cc      CALL  REDUCT( SETX2,Z,D,NMK,0,K,MJ1,LAG,X )                       
      CALL  REDUCT( SETX2,Z,NMK,0,K,MJ1,LAG,X )                       
      GO TO 100                                                         
C                                                                       
cc   30 CALL  SETLAG( K )                                                 
cc      CALL  REDUCT( SETX2,Z,D,NMK,0,K,MJ1,LAG,X )                       
   30 CALL  SETLAG( K,LG2(1),LG2(2),LG2(3),LG2(4),LG2(5) )
      DO 31 I=1,K
         LG1(1,I) = L1(I)
         LG1(2,I) = L2(I)
         LG1(3,I) = L3(I)
   31 CONTINUE
      CALL  REDUCT( SETX2,Z,NMK,0,K,MJ1,LAG,X )                       
      GO TO 100                                                         
C                                                                       
cc   40 CALL  REDUCT( SETX4,Z,D,NMK,0,K,MJ1,LAG,X )                       
   40 CALL  REDUCT( SETX4,Z,NMK,0,K,MJ1,LAG,X )                       
      GO TO 100                                                         
C                                                                       
cc   50 CALL  REDREG( K )                                                 
   50 CONTINUE
      DO 51 J=1,K+1
         DO 51 I=1,20
            II = (I-1)*4+1
            G(I,J) = CHAR(F(II,J)) // CHAR(F(II+1,J))
     *                             // CHAR(F(II+2,J)) //CHAR(F(II+3,J))
   51 CONTINUE
cc      CALL  REDUCT( SETX5,Z,D,NMK,0,K,MJ1,LAG,X )                       
      CALL  REDUCT( SETX5,Z,NMK,0,K,MJ1,LAG,X )                       
      GO TO 100                                                         
C                                                                       
   60 CONTINUE                                                          
	SD0 = SUM
      CONST = CNST
cc      SD0 = 0.D0                                                        
cc      DO 65  I=1,N                                                      
cc   65 SD0 = SD0 + Z(I)*Z(I)                                             
cc      SD0 = SD0 / N                                                     
cc      READ( 5,9 )     CONST                                             
cc      WRITE( 6,12 )   SD0 , CONST                                       
cc      CALL  REDLAG( K )                                                 
      DO 66 I=1,K
         L1(I) = LG1(1,I)
         L2(I) = LG1(2,I)
         L3(I) = LG1(3,I)
   66 CONTINUE
cc      CALL  REDUCT( SETX6,Z,D,NMK,0,K,MJ1,LAG,X )                       
      CALL  REDUCT( SETX6,Z,NMK,0,K,MJ1,LAG,X )                       
      GO TO 100                                                         
C                                                                       
cc   70 CALL  REDUCT( SETX7,Z,D,NMK,0,K,MJ1,LAG,X )                       
   70 CALL  REDUCT( SETX7,Z,NMK,0,K,MJ1,LAG,X )                       
C                                                                       
  100 CONTINUE                                                          
cc	CLOSE( MT )
C                                                                       
C          ---------------                                              
C          MAICE PROCEDURE                                              
C          ---------------                                              
C                                                                       
cc      CALL  ARMFIT( X,K,LAG,NMK,ISW,TITLE,MJ1,A,SD,M )                  
      IFG=0
      CALL ARMFIT( X,K,LAG,NMK,ISW,MJ1,A1,M,SD,AIC,DIC,SDM,AICM,
     *             IFG,LU )
C                                                                       
C          ------------------                                           
C          BAYESIAN PROCEDURE                                           
C          ------------------                                           
C                                                                       
cc      CALL  SBBAYS( X,D,K,NMK,IPR,MJ1,A,SD )                            
      CALL  SBBAYS( X,K,NMK,IPR,MJ1,A2,SDB,EK,AICB,IND,C,C1,C2,B,
     *OEIC,ESUM,OMEAN,OM  )
C                                                                       
cc      IF( IMODEL .EQ. 1 )  CALL  NRASPE( SD,A,B,K,0,121,TITLE )         
      IF( IMODEL .EQ. 1 )  CALL  NRASPE( SDB,A2,B,K,0,120,SXX )
C                                                                       
      NPS = LAG+1                                                       
C                                                                       
C          -------------------------                                    
C          PREDICTION ERROR CHECKING                                    
C          -------------------------                                    
C                                                                       
      GO TO ( 110,120,120,150,130,140 ), IMODEL                         
      GO TO 150                                                         
C                                                                       
cc  110 CALL CHECK( PRDCT1,Z,A,K,0,IL,NPS,N,0,MJ,E )                      
  110 CALL CHECK( PRDCT1,Z,A2,K,0,IL,NPS,N,0,MJ,E,FA,EMEAN,VARI,SKEW,
     *            PEAK,COV,MJ2 )
      GO TO 150                                                         
C                                                                       
cc  120 CALL  CHECK( PRDCT2,Z,A,K,0,IL,NPS,N,0,MJ,E )                     
  120 CALL  CHECK( PRDCT2,Z,A2,K,0,IL,NPS,N,0,MJ,E,FA,EMEAN,VARI,SKEW,
     *             PEAK,COV,MJ2 )
      GO TO 150                                                         
C                                                                       
cc  130 CALL CHECK( PRDCT3,Z,A,K,0,IL,NPS,N,0,MJ,E )                      
  130 CALL CHECK( PRDCT3,Z,A2,K,0,IL,NPS,N,0,MJ,E,FA,EMEAN,VARI,SKEW,
     *            PEAK,COV,MJ2 )
C                                                                       
cc  140 CALL  CHECK( PRDCT6,Z,A,K,0,IL,NPS,N,0,MJ,E )                     
  140 CALL  CHECK( PRDCT6,Z,A2,K,0,IL,NPS,N,0,MJ,E,FA,EMEAN,VARI,SKEW,
     *             PEAK,COV,MJ2 )
C                                                                       
  150 CONTINUE                                                          
      RETURN
C                                                                       
    1 FORMAT( 16I5 )                                                    
    2 FORMAT( /1H ,'ORIGINAL DATA INPUT DEVICE  MT =',I4 )              
    3 FORMAT( ' PROGRAM TIMSAC 78.1.3',/,'   SCALAR TIME SERIES MODEL FI
     *TTING;     BAYESIAN PROCEDURE ( ALL SUBSET REGRESSIN TYPE )' )    
    4 FORMAT( //1H ,'MODEL TYPE',I2,/,'   < AUTOREGRESSIVE MODEL >',/,  
     11H ,10X,'Z(I) = A(1)*Z(I-1) + A(2)*Z(I-2) + ... + A(M)*Z(I-M) + E(
     2I)' )                                                             
    5 FORMAT( //1H ,'MODEL TYPE',I2,/,'   < NON-LINEAR MODEL >',/,1H ,10
     1X,'Z(I) = A(1)*Y(I,1) + A(2)*Y(I,2) + ... + A(K)*Y(I,K) + E(I)' ) 
    6 FORMAT( 1H ,2X,'WHERE',/,11X,'M:     ORDER OF THE MODEL',/,11X,'E(
     1I):  GAUSSIAN WHITE NOISE WITH MEAN 0  AND  VARIANCE SD(M).' )    
    7 FORMAT( 1H ,'FITTING UP TO THE ORDER',I3,2X,'IS TRIED' )          
    8 FORMAT( 1H ,'MAXIMUM LAG =',I4,/,' NUMBER OF REGRESSORS =',I4 )   
    9 FORMAT( F10.0 )                                                   
   11 FORMAT( //1H ,'MODEL TYPE',I2,/,'   < EXPONENTIALLY DAMPED NON-LIN
     *EAR MODEL >',/,1H ,10X,'Z(I) = A(1)*Y(I,1) + A(2)*Y(I,2) + ... + A
     *(K)*Y(I,K) + E(I)' )                                              
   12 FORMAT( 1H ,'SIGMA2 =',D15.5,5X,'CONST =',D15.5 )                 
C                                                                       
      E N D                                                             
cc      REAL FUNCTION  BICOEF * 8( K,J )                                  
      REAL*8 FUNCTION  BICOEF( K,J )                                  
C                                                                       
C     THIS FUNCTION RETURNS BINOMIAL COEFFICIENTS                       
C                                                                       
C          F(K,J) = K]/(J]*(K-J)])                                      
C                                                                       
C       INPUTS:                                                         
C          K:     NUMBER OF OBJECTS                                     
C          J:     NUMBER OF OBJECTS TAKEN                               
C                                                                       
C       OUTPUT:                                                         
C          F:     NUMBER OF COMBINATIONS OF SELECTING J OBJECTS FROM    
C                 A SET OF K OBJECTS                                    
C                                                                       
      IMPLICIT REAL * 8 ( A-H , O-Z )                                   
C                                                                       
      KMJ = K-J                                                         
      SUM = 0.D0                                                        
      DO 10   I=1,K                                                     
      DI = I                                                            
   10 SUM = SUM + DLOG( DI )                                            
C                                                                       
      IF( J .EQ. 0 )   GO TO 30                                         
      DO 20   I=1,J                                                     
      DI = I                                                            
   20 SUM = SUM - DLOG( DI )                                            
C                                                                       
   30 IF( KMJ .EQ. 0 )   GO TO 50                                       
      DO 40   I=1,KMJ                                                   
      DI = I                                                            
   40 SUM = SUM - DLOG( DI )                                            
C                                                                       
   50 BICOEF = DEXP( SUM )                                              
      RETURN                                                            
C                                                                       
      END                                                               
cc      SUBROUTINE  CHECK( PRDCT,X,A,K,L,IL,NPS,NPE,IPR,MJ,E )            
      SUBROUTINE  CHECK( PRDCT,X,A,K,L,IL,NPS,NPE,IPR,MJ,E,F,EMEAN,VARI,
     *                   SKEW,PEAK,COV,MJ2 )            
C                                                                       
C     THIS SUBROUTINE DRAWS HISTGRAMS AND AUTOCOVARIANCE FUNCTION OF ORI
C     DATA OR PREDICTION ERRORS.                                        
C                                                                       
C       ----------------------------------------------------------------
C       THE FOLLOWING SUBROUTINES ARE DIRECTLY CALLED BY THIS SUBROUTINE
C             GRAPH1                                                    
C             GRAPH2                                                    
C             MOMENT                                                    
C             (PRDCT)                                                   
C       ----------------------------------------------------------------
C       INPUTS:                                                         
C          PRDCT:  EXTERNAL SUBROUTINE DESIGNATION                      
C          X:      ORIGINAL DATA                                        
C          A:      REGRESSION COEFFICIENTS                              
C          K:      NUMBER OF REGRESSORS                                 
C          L:      MA-ORDER ( THIS ARGUMENT IS ONLY USED FOR THE CHECKIN
C                  OF AR-MA MODEL)                                      
C          IL:     MAXIMUM SPAN OF LONG RANGE PREDICTION                
C                  =0     ANALYSIS OF ORIGINAL DATA                     
C                  >0     ANALYSIS OF MULTI-STEP (UP TO IL) PREDICTION E
C          NPS:    PREDICTION STARTING POSITION                         
C          NPE:    PREDICTION ENDING POSITION                           
C          IPR:    =0  MATRIX OF SEVERAL STEP AHEAD PREDICTION ERRORS SU
C                  =1  MATRIX OF SEVERAL STEP AHEAD PREDICTION ERRORS IS
C                      OUT                                              
C          MJ:     ABSOLUTE DIMENSION OF E IN THE MAIN PROGRAM          
C                                                                       
C       OUTPUT:                                                         
C          E:      SEVERAL-STEPS PREDICTION ERRORS                      
C                                                                       
C                                                                       
      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
CC      REAL * 4  X , F , E                                               
      DIMENSION  X(1) , E(MJ,1)                                         
cc      DIMENSION  A(1) , COV(120)                                        
      DIMENSION  A(1) , COV(MJ2)                                        
      DIMENSION  EMEAN(IL), VARI(IL), SKEW(IL), PEAK(IL)
cc      COMMON  /COMXX/ F(2000)                                           
      DIMENSION  F(NPE-NPS+1,IL)
C                                                                       
C                                                                       
      ISTEP = 1                                                         
      ISW = IL                                                          
      LAGH = 100                                                        
      N = NPE - NPS - 1                                                 
      IF( LAGH .GE. N )     LAGH = N - 1                                
      LAG1 = LAGH + 1                                                   
      NMK = N - K                                                       
      IF( ISW .GT. 0 )     GO TO 20                                     
C                                                                       
      DO 10  I=NPS,NPE                                                  
   10 E(I,1) = X(I)                                                     
      IL = 1                                                            
      GO TO 36                                                          
C                                                                       
C       ---  SEVERAL STEP AHEAD PREDICITON  ---                         
C                                                                       
   20 CONTINUE                                                          
      CALL  PRDCT( X,A,K,L,IL,NPS,NPE,MJ,E )                            
C                                                                       
C       ---  PREDICTION ERROR  ---                                      
C                                                                       
      DO 30  II=NPS,NPE                                                 
         I = NPE-II+NPS                                                 
         DO 30  J=1,IL                                                  
         JJ = I-J+1                                                     
   30 E(I,J) = X(I) - E(JJ,J)                                           
      IF( IL .EQ. 1 )     GO TO 34                                      
      DO 35  J=2,IL                                                     
      JJ = J-1                                                          
      DO 35  I=1,JJ                                                     
      II = I+NPS-1                                                      
   35 E(II,J) = 0.D0                                                    
   34 CONTINUE                                                          
cc      WRITE( 6,690 )                                                    
cc      IF( IPR .EQ. 0 )   GO TO 36                                       
cc      WRITE( 6,670 )     (I,I=1,IL)                                     
cc      DO 37  I=NPS,NPE                                                  
cc   37 WRITE( 6,640 )     I , (E(I,J),J=1,IL)                            
   36 CONTINUE                                                          
C                                                                       
C       ---  MOMENT COMPUTATION  ---                                    
C                                                                       
      DO 50  KK=1,IL                                                    
C                                                                       
      II = NPS+KK-1                                                     
      DO 40  I=II,NPE                                                   
      J = I - II + 1                                                    
cc   40 F(J) = E(I,KK)                                                    
   40 F(J,KK) = E(I,KK)
      NMK = NPE-NPS-(KK-2)                                              
C                                                                       
cc      CALL  MOMENT( F,NMK,EMEAN,VARI,SKEW,PEAK )                        
      CALL  MOMENT( F(1,KK),NMK,EMEAN(KK),VARI(KK),SKEW(KK),PEAK(KK) )
C                                                                       
cc      IF( ISW .GT. 1 )     WRITE( 6,610 )   KK                          
cc      IF( ISW .EQ. 0 )     WRITE( 6,650 )                               
cc      WRITE( 6,620 )     EMEAN , VARI , SKEW , PEAK                     
C                                                                       
C       ---  HISTOGRAM OF F(I)  ---                                     
C                                                                       
cc      SIG = DSQRT(VARI)                                                 
cc      CALL  GRAPH1( F,1,NMK,SIG )                                       
   50 CONTINUE                                                          
C                                                                       
C       ---  AUTOCORRELATION FUNCTION COMPUTATION  ---                  
C                                                                       
      DO 100  KK=1,IL                                                   
C                                                                       
      DO 70   II=1,LAG1                                                 
      JJ = NPS + KK - 1                                                 
      IE = NPE - II + 1                                                 
      SUM = 0.D0                                                        
      DO 60   I=JJ,IE                                                   
      J = I + II- 1                                                     
   60 SUM = SUM + E(I,KK)*E(J,KK)                                       
   70 COV(II) = SUM / (NPE-NPS-KK+2)                                    
C                                                                       
      COV1 = COV(1)                                                     
      DO 80   I=1,LAG1                                                  
   80 COV(I) = COV(I) / COV1                                            
C                                                                       
cc      IF( ISW .GT. 0 )   WRITE( 6,630 )     KK                          
cc      IF( ISW .EQ. 0 )   WRITE( 6,660 )                                 
cc      WRITE( 6,680 )     (COV(I),I=1,50)                                
C                                                                       
C       ---  AUTOCORRELATION FUNCTION DISPLAY  ---                      
C                                                                       
      SD = NPE - NPS + 1                                                
      SD = DSQRT( 1.D0/SD )                                             
cc      CALL  GRAPH2( COV,LAG1,SD )                                       
C                                                                       
      IF( ISTEP .EQ. KK )     GO TO 110                                 
C                                                                       
  100 CONTINUE                                                          
C                                                                       
C                                                                       
  110 CONTINUE                                                          
      RETURN                                                            
  610 FORMAT( //1H ,I3,'-STEP AHEAD PREDICTION ERROR' )                 
  620 FORMAT( 1H ,'MEAN       =',D15.8,/,' VARIANCE   =',D15.8,/,' SKEWN
     1ESS   =',D15.8,/,' PEAKEDNESS =',D15.8 )                          
  630 FORMAT( //1H ,'AUTOCORRELATION FUNCTION OF ',I3,'-STEP AHEAD PREDI
     1CTION ERROR' )                                                    
  640 FORMAT( 1H ,I5,5X,10D12.4 )                                       
  650 FORMAT( //,' ORIGINAL DATA' )                                     
  660 FORMAT( //1H ,'AUTOCORRELATION FUNCTION OF ORIGINAL DATA' )       
  670 FORMAT( 1H ,10X,10(4X,'J =',I2,3X) )                              
  680 FORMAT( 1H ,10D13.5 )                                             
  690 FORMAT( ///1H ,45(1H-),2X,'<< J-STEP AHEAD PREDICTION ERROR >>',2X
     1,45(1H-) )                                                        
      END                                                               
      SUBROUTINE  MOMENT( X,N,F1,F2,F3,F4 )                             
C                                                                       
C          +--------------------+                                       
C          ! MOMENT COMPUTATION !                                       
C          +--------------------+                                       
C                                                                       
C     THIS SUBROUTINE COMPUTES MOMENTS.                                 
C                                                                       
C       INPUTS:                                                         
C          X:     ORIGINAL DATA VECTOR                                  
C          N:     DATA LENGTH                                           
C                                                                       
C       OUTPUTS:                                                        
C          F1:    MEAN OF X                                             
C          F2:    VARIANCE OF X                                         
C          F3:    SKEWNESS OF X                                         
C          F4:    PEAKEDNESS OF X                                       
C                                                                       
      IMPLICIT  REAL * 8  ( F )                                         
CC      DIMENSION  X(1)                                                   
	REAL*8  X(1)
C                                                                       
      FN = N                                                            
      FSUM = 0.D0                                                       
      DO  10     I=1,N                                                  
   10 FSUM = FSUM + X(I)                                                
C                                                                       
      F1 = FSUM / FN                                                    
C                                                                       
      F2 = 0.D0                                                         
      F3 = 0.D0                                                         
      F4 = 0.D0                                                         
      DO  20     I=1,N                                                  
      FF = X(I) - F1                                                    
      F2 = F2 + FF*FF                                                   
      F3 = F3 + FF**3                                                   
   20 F4 = F4 + FF**4                                                   
C                                                                       
      F2 = F2 / FN                                                      
      F3 = F3 / (FN*F2*DSQRT(F2))                                       
      F4 = F4 / (FN*F2*F2)                                              
C                                                                       
      RETURN                                                            
      END                                                               
      SUBROUTINE  PRDCT1( Z,A,M,L,IL,NPS,NPE,MJ,EZ )                    
C                                                                       
C     THIS SUBROUTINE COMPUTES SEVARAL STEP AHEAD PREDICTION VALUE OF AN
C     AUTOREGRESSIVE MOVING AVERAGE MODEL.                              
C                                                                       
C       INPUTS:                                                         
C          Z:      ORIGINAL DATA VECTOR                                 
C          A:      AR-MA COEFFICIENTS                                   
C          M:      AR-ORDER                                             
C          L:      MA-ORDER                                             
C          IL:     MAXIMUM SPAN OF LONG RANGE PREDICTION                
C          NPS:    PREDICTION STARTING POSITION                         
C          NPE:    PREDICTION ENDING POSITION                           
C          MJ:     ABSOLUTE DIMENSION OF EZ                             
C                                                                       
C       OUTPUT:                                                         
C          EZ:     PREDICTION VALUE MATRIX                              
C                                                                       
      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
CC      REAL * 4  Z(1) , EZ(MJ,1)                                         
      DIMENSION  Z(1) , EZ(MJ,1)
	DIMENSION  A(1)                                                   
C                                                                       
C                                                                       
      DO  100     II=NPS,NPE                                            
C                                                                       
      DO  90     KK=1,IL                                                
      KKM1 = KK - 1                                                     
      SUM = 0.D0                                                        
      IF( KK .EQ. 1 )     GO TO 30                                      
      DO  20     I=1,KKM1                                               
      KI = KK - I                                                       
   20 SUM = SUM + A(I)*EZ(II,KI)                                        
   30 IF( KK .GT. M )     GO TO 50                                      
      DO  40     I=KK,M                                                 
      I1 = II + KKM1 - I                                                
   40 SUM = SUM + A(I)*Z(I1)                                            
C                                                                       
   50 IF( L .LE. 0 )     GO TO 90                                       
      IF( KK .GT. L )     GO TO 90                                      
      DO  60     I=KK,L                                                 
      I1 = M + I                                                        
      I2 = II + KKM1 - I                                                
      IF( I2 .GE. II )     GO TO 60                                     
      SUM = SUM + A(I1)*(Z(I2)-EZ(I2,1))                                
   60 CONTINUE                                                          
   90 EZ(II,KK) = SUM                                                   
  100 CONTINUE                                                          
C                                                                       
      RETURN                                                            
      END                                                               
      SUBROUTINE  PRDCT2( Z,A,K,L,IL,NPS,NPE,MJ1,EZ )                   
C                                                                       
C     THIS SUBROUTINE COMPUTES SEVERAL STEPS AHEAD PREDICTION VALUES OF 
C     NON-LINEAR REGRESSION MODEL.                                      
C                                                                       
C       INPUTS:                                                         
C          Z:      ORIGINAL DATA VECTOR                                 
C          A:      VECTOR OF AR-COEFFICIENTS                            
C          K:      ORDER OF THE AR MODEL                                
C          L:      THIS DUMMY VARIABLE IS NOT REFERENCED IN THIS SUBROUT
C          IL:     MAXIMUM SPAN OF LONG RANGE PREDICTION                
C          NPS:    PREDICTION STARTING POSITION                         
C          NPE:    PREDICTION ENDING POSITION                           
C          MJ1:    ABSOLUTE DIMENSION OF EZ                             
C                                                                       
C       OUTPUT:                                                         
C          EZ:     PREDICTION VALUE MATRIX                              
C                                                                       
CC      IMPLICIT  REAL*8 ( A-D,O-Y )                                      
      IMPLICIT  REAL*8 ( A-H,O-Z ) 
cc      DIMENSION  Z(1) , A(1) , EZ(MJ1,1) , Y(20)                        
      DIMENSION  Z(1) , A(1) , EZ(MJ1,1) , Y(IL)                        
cc      REAL * 4   CSTDMY, SD0DMY
      COMMON     / BBB /  LAG1(50) , LAG2(50) , LAG3(50), CSTDMY, SD0DMY
CC      INTEGER  RETURN                                                   
C                                                                       
      DO 100  II=NPS,NPE                                                
         DO 50  J1=1,IL                                                 
            JJ = J1-1                                                   
            SUM = 0.D0                                                  
            DO 40  J=1,K                                                
               XX = 1.D0                                                
               LAG = LAG1(J)                                            
CC               ASSIGN 10 TO RETURN                                      
CC               GO TO 200                                                
CC   10          XX = XX*X
               X = 1.D0
               IF ( LAG .GT. 0 ) THEN
                  I = II+JJ-LAG
                  IF ( I .GE. II ) THEN
                     I = I - II + 1
                     X = Y(I)
                  ELSE
                     X = Z(I)
                  END IF
               END IF
               XX = XX*X                                               
C                                                
               LAG = LAG2(J)                                            
CC               ASSIGN 20 TO RETURN                                      
CC               GO TO 200                                                                                               
CC   20          XX = XX*X
   	         X = 1.D0
               IF ( LAG .GT. 0 ) THEN
                  I = II+JJ-LAG
                  IF ( I .GE. II ) THEN
                     I = I - II + 1
                     X = Y(I)
                  ELSE
                     X = Z(I)
                  END IF
               END IF
               XX = XX*X
C
               LAG = LAG3(J)                                            
CC               ASSIGN 30 TO RETURN                                      
CC               GO TO 200                                                
CC   30          XX = XX*X
               X = 1.D0
               IF ( LAG .GT. 0 ) THEN
                  I = II+JJ-LAG
                  IF ( I .GE. II ) THEN
                     I = I - II + 1
                     X = Y(I)
                  ELSE
                     X = Z(I)
                  END IF
               END IF
               XX = XX*X
C                                                                       
   40       SUM = SUM + A(J)*XX                                         
   50    Y(J1) = SUM                                                    
C                                                                       
      DO  100  J=1,IL                                                   
      EZ(II,J) = Y(J)                                                   
  100 CONTINUE                                                          
CC      GO TO 300
C	-----  INTERNAL SUBROUTINE  -----
C
CC  200 X = 1.D0
CC      IF ( LAG .LE. 0 )     GO TO 220
CC      I = II+JJ-LAG
CC      IF ( I .GE. II )      GO TO 210
CC      X = Z(I)
CC      GO TO 220
CC  210 I = I - II + 1
CC      X = Y(I)
CC  220 GO TO RETURN, ( 10,20,30 )
C	----------------------------------
C                                                                                                                                                              
CC  300 RETURN                                                            
      RETURN
	END                                                               
      SUBROUTINE  PRDCT3( Z,A,K,L,IL,NPS,NPE,MJ1,EZ )                   
C                                                                       
C     THIS SUBROUTINE COMPUTES SEVERAL STEP AHEAD PREDICTION VALUE OF A 
C     NON-LINEAR REGRESSION MODEL.                                      
C                                                                       
C       INPUTS:                                                         
C          Z:      ORIGINAL DATA VECTOR                                 
C          A:      VECTOR OF AR-COEFFICIENTS                            
C          K:      ORDER OF THE AR MODEL                                
C          L:      THIS DUMMY VARIABLE IS NOT REFERENCED IN THIS SUBROUT
C          IL:     MAXIMUM SPAN OF LONG RANGE PREDICTION                
C          NPS:    PREDICTION STARTING POSITION                         
C          NPE:    PREDICTION ENDING POSITION                           
C          MJ1:    ABSOLUTE DIMENSION OF EZ                             
C                                                                       
C          G:      REGRESSOR AND REGRESSAND SPECIFICATION STATEMENT, LOA
C                  AT SUBROUTINE REDREG                                 
C                                                                       
C       OUTPUT:                                                         
C          EZ:     PREDICTION VALUE MATRIX                              
C                                                                       
CC      REAL*8  A                                                         
      CHARACTER(4) G
      COMMON     / EEE /  G(20,31)                                      
CC      DIMENSION  Z(1) , A(1) , EZ(MJ1,1)                                
      REAL*8  Z(1) , A(1) , EZ(MJ1,1)
      CHARACTER(4) F
      DIMENSION  E(20) , F(20)                                          
      CHARACTER(4) A1, A2, A4,
     *             A5, A6, A7, A8,
     *             T, X
cc      DATA      A1/4H//  / , A2/4H)   / , A4/4H+   / ,                  
cc     *          A5/4H-   / , A6/4H*   / , A7/4H/   / , A8/4H**  / ,     
cc     *          T /4HT   / , X /4HX   /                                 
      DATA      A1/'//  '/ , A2/')   '/ , A4/'+   '/ ,                  
     *          A5/'-   '/ , A6/'*   '/ , A7/'/   '/ , A8/'**  '/ ,     
     *          T /'T   '/ , X /'X   '/                                 
      CHARACTER(4) ASIN, ACOS, CLOG,
     *             BLOG, AEXP, AINV,
     *             SQ, CUBE, ROOT
cc      DATA      ASIN/4HSIN / , ACOS/4HCOS / , CLOG/4HLOG / ,            
cc     *          BLOG/4HLOG1/ , AEXP/4HEXP / , AINV/4HINV / ,            
cc     *          SQ  /4HSQ  / , CUBE/4HCUBE/ , ROOT/4HSQRT/              
      DATA      ASIN/'SIN '/ , ACOS/'COS '/ , CLOG/'LOG '/ ,            
     *          BLOG/'LOG1'/ , AEXP/'EXP '/ , AINV/'INV '/ ,            
     *          SQ  /'SQ  '/ , CUBE/'CUBE'/ , ROOT/'SQRT'/              
      CHARACTER(4) FF, GG
C                                                                       
C                                                                       
      P = 0
      Q = 0
      W = 0
      V = 0
C
      DO 3000  II=NPS,NPE                                               
      DO 2000  J1=1,IL                                                  
      JJ = J1 - 1                                                       
      SUM = 0.0                                                         
      DO 1000  J=1,K                                                    
      DO 5  I=1,20                                                      
    5 F(I) = G(I,J)                                                     
      I = 0                                                             
      Y = 1.0                                                           
  500 I = I + 1                                                         
      IF(I.GT.20)  GO TO 10                                             
      FF = F(I)                                                         
C      ENCODE( 4,1,GG )   FF                                             
      GG = FF(1:1)
C                                                                       
      IF(FF.NE.A1)  GO TO 20                                            
   10 Y = Y*U                                                           
      GO TO 1000                                                        
C                                                                       
   20 IF(FF.NE.A2)  GO TO 30                                            
      Y = Y*U                                                           
      GO TO 500                                                         
C                                                                       
   30 IF(GG.NE.T)  GO TO 50                                             
CC      ASSIGN 40 TO ISUB                                                 
CC      GO TO 300
      Q = P
      P = W                                                             
      W = V                                                             
      V = U                                                             
C
C g77 NOT supports DECODE function and uses special form.
C   40 DECODE( 4,3,FF )   S                                              
      READ(UNIT=FF, FMT=3) S
C
      U = II + S + JJ                                                   
      GO TO 500                                                         
C                                                                       
   50 IF(GG.NE.X)  GO TO 60                                             
CC      ASSIGN 55 TO ISUB                                                 
CC      GO TO 300
      Q = P
      P = W                                                             
      W = V                                                             
      V = U                                                             
C
C g77 NOT supports DECODE function and uses special form.
C   55 DECODE( 4,2,FF )   LL                                             
      READ(UNIT=FF, FMT=2) LL
C
      IT = II + LL + JJ                                                 
      I1 = LL + J1                                                      
      IF( -LL-J1 .GE. 0 )   U = Z(IT)                                   
      IF( -LL-J1 .LT. 0 )   U = E(I1)                                   
      GO TO 500                                                         
C                                                                       
   60 IF(FF.NE.A4)  GO TO 70                                            
      U = V + U                                                         
CC      GO TO 400                                                         
C-----------------------------------------------------------            
C                    STACK AREA                                         
C-----------------------------------------------------------            
C                                                                       
C          POP UP STACK                                                 
C                                                                       
      V = W                                                             
      W = P                                                             
      P = Q                                                             
      GO TO 500                                                         
C C                                                                       
   70 IF(FF.NE.A5)  GO TO 80                                            
      U = V- U                                                          
CC      GO TO 400                                                         
C-----------------------------------------------------------            
C                    STACK AREA                                         
C-----------------------------------------------------------            
C                                                                       
C          POP UP STACK                                                 
C                                                                       
      V = W                                                             
      W = P                                                             
      P = Q                                                             
      GO TO 500                                                         
C C                                                                       
   80 IF(FF.NE.A6)  GO TO 90                                            
      U = V * U                                                         
CC      GO TO 400                                                         
C-----------------------------------------------------------            
C                    STACK AREA                                         
C-----------------------------------------------------------            
C                                                                       
C          POP UP STACK                                                 
C                                                                       
      V = W                                                             
      W = P                                                             
      P = Q                                                             
      GO TO 500                                                         
C C                                                                       
   90 IF(FF.NE.A7)  GO TO 100                                           
      U = V / U                                                         
CC      GO TO 400                                                         
C-----------------------------------------------------------            
C                    STACK AREA                                         
C-----------------------------------------------------------            
C                                                                       
C          POP UP STACK                                                 
C                                                                       
      V = W                                                             
      W = P                                                             
      P = Q                                                             
      GO TO 500                                                         
C C                                                                       
  100 IF(FF.NE.A8)  GO TO 120                                           
      U = V**U                                                          
CC      GO TO 400                                                         
C-----------------------------------------------------------            
C                    STACK AREA                                         
C-----------------------------------------------------------            
C                                                                       
C          POP UP STACK                                                 
C                                                                       
      V = W                                                             
      W = P                                                             
      P = Q                                                             
      GO TO 500                                                         
C C                                                                       
  120 IF(FF.NE.ASIN)  GO TO 130                                         
      U = SIN( U )                                                      
      GO TO 500                                                         
C                                                                       
  130 IF(FF.NE.ACOS)  GO TO 140                                         
      U = COS( U )                                                      
      GO TO 500                                                         
C                                                                       
  140 IF(FF.NE.CLOG)  GO TO 150                                         
      U = ALOG( U )                                                     
      GO TO 500                                                         
C                                                                       
  150 IF(FF.NE.BLOG)  GO TO 160                                         
      U = ALOG10( U )                                                   
      GO TO 500                                                         
C                                                                       
  160 IF(FF.NE.AEXP)  GO TO 170                                         
      U = EXP( U )                                                      
      GO TO 500                                                         
C                                                                       
  170 IF(FF.NE.AINV)  GO TO  180                                        
      U = 1.D0 / U                                                      
      GO TO 500                                                         
C                                                                       
  180 IF(FF.NE.SQ)  GO TO 190                                           
      U = U**2                                                          
      GO TO 500                                                         
C                                                                       
  190 IF(FF.NE.ROOT)  GO TO 220                                         
      U = SQRT( U )                                                     
      GO TO 500                                                         
C                                                                       
  220 IF(FF.NE.CUBE)  GO TO 200                                         
      U = U**3                                                          
      GO TO 500                                                         
C                                                                       
  200 CONTINUE                                                          
CC      ASSIGN 210 TO ISUB                                                
CC      GO TO 300
      Q = P
      P = W                                                             
      W = V                                                             
      V = U                                                             
C
C g77 NOT supports DECODE function and uses special form.
C  210 DECODE( 4,4,FF )   U                                              
      READ(UNIT=FF, FMT=4) U
C
      GO TO 500                                                         
C                                                                       
 1000 SUM = SUM + A(J)*Y                                                
 2000 E(J1) = SUM                                                       
      DO 3000  J=1,IL                                                   
 3000 EZ(II,J) = E(J)                                                   
      RETURN                                                            
C                                                                       
C-----------------------------------------------------------            
C                    STACK AREA                                         
C-----------------------------------------------------------            
C                                                                       
C          POP UP STACK                                                 
C                                                                       
CC  400 V = W                                                             
CC      W = P                                                             
CC      P = Q                                                             
CC      GO TO 500                                                                      
C
C          PUSH DOWN STACK                                              
C
CC  300 Q = P                                                             
CC      P = W                                                             
CC      W = V                                                             
CC      V = U                                                             
CC      GO TO ISUB, ( 40,55,210 )        
C                                                                       
    1 FORMAT( A1,3X )                                                   
    2 FORMAT( 1X,I3 )                                                   
    3 FORMAT( 1X,F3.0 )                                                 
    4 FORMAT( F4.0 )                                                    
      E N D                                                             
      SUBROUTINE  PRDCT6( Z,A,K,L,IL,NPS,NPE,MJ1,EZ )                   
C                                                                       
C     THIS SUBROUTINE COMPUTES SEVERAL STEP AHEAD PREDICITON VALUE OF A 
C     NON-LINEAR REGRESSION MODEL.                                      
C                                                                       
C       INPUTS:                                                         
C          Z:      ORIGINAL DATA VECTOR                                 
C          A:      VECTOR OF AR-COEFFICIENTS                            
C          K:      ORDER OF THE AR MODEL                                
C          L:      THIS DUMMY VARIABLE IS NOT REFERENCED IN THIS SUBROUT
C          IL:     MAXIMUM SPAN OF LONG RANGE PREDICTION                
C          NPS:    PREDICTION STARTING POSITION                         
C          NPE:    PREDICTION ENDING POSITION                           
C          MJ1:    ABSOLUTE DIMENSION OF EZ                             
C                                                                       
C       OUTPUT:                                                         
C          EZ:     PREDICTION VALUE MATRIX                              
C                                                                       
CC      IMPLICIT  REAL*8 ( A-D,O-Y )                                      
      IMPLICIT  REAL * 8  ( A-H , O-Z )
cc      REAL * 4   SD , CONST , C                                         
cc      DIMENSION  Z(1) , A(1) , EZ(MJ1,1) , Y(20)                        
      DIMENSION  Z(1) , A(1) , EZ(MJ1,1) , Y(IL)                        
      COMMON     / BBB /  LAG1(50) , LAG2(50) , LAG3(50) , SD , CONST   
cc      INTEGER  RETURN                                                   
C                                                                       
      C = -1.0/(SD*CONST)                                               
      DO 100  II=NPS,NPE                                                
         DO 50  J1=1,IL                                                 
            JJ = J1-1                                                   
            SUM = 0.D0                                                  
            DO 40  J=1,K                                                
               XX = 1.D0                                                
               LAG = LAG1(J)                                            
CC               ASSIGN 10 TO RETURN                                      
CC               GO TO 200
CC   10          XX = XX*X
               X = 1.D0                                                          
               IF( LAG .GT. 0 )  THEN                                      
                  I = II+JJ-LAG                                                     
                  IF( I .GE. II )  THEN
                     I = I-II+1                                                    
                     ZZ = C*Y(I)**2                                                    
                     X = Y(I)*EXP( ZZ )                                    
                  ELSE
                     ZZ = C*Z(I)**2                                                    
                     X = Z(I)*EXP( ZZ )                                                                                                         
                  END IF
               END IF                                                
               XX = XX*X                                                
C                                                                       
               LAG = LAG2(J)                                            
CC               ASSIGN 20 TO RETURN                                      
CC               GO TO 200
CC   20          XX = XX*X                                                
               X = 1.D0                                                          
               IF( LAG .GT. 0 )  THEN                                      
                  I = II+JJ-LAG                                                     
                  IF( I .GE. II )  THEN
                     I = I-II+1                                                    
                     ZZ = C*Y(I)**2
                     X = Y(I)*EXP( ZZ )
                  ELSE
                     ZZ = C*Z(I)**2                                                    
                     X = Z(I)*EXP( ZZ )                                                                                                         
                  END IF
               END IF                                                
               XX = XX*X  
C                                                                       
               LAG = LAG3(J)                                            
CC               ASSIGN 30 TO RETURN                                      
CC               GO TO 200
CC   30          XX = XX*X                                                
               X = 1.D0                                                          
               IF( LAG .GT. 0 )  THEN                                      
                  I = II+JJ-LAG                                                     
                  IF( I .GE. II )  THEN
                     I = I-II+1                                                    
                     ZZ = C*Y(I)**2                                                    
                     X = Y(I)*EXP( ZZ )                                    
                  ELSE
                     ZZ = C*Z(I)**2                                                    
                     X = Z(I)*EXP( ZZ )                                                                                                         
                  END IF
               END IF                                                
               XX = XX*X
C                                                                       
   40       SUM = SUM + A(J)*XX                                         
   50    Y(J1) = SUM                                                    
C                                                                       
      DO  100  J=1,IL                                                   
      EZ(II,J) = Y(J)                                                   
  100 CONTINUE                                                          
      GO TO 300                                                         
C                                                                       
C   -----  INTERNAL SUBROUTINE  -----                                   
C                                                                       
CC  200 X = 1.D0                                                          
CC      IF( LAG .LE. 0 )   GO TO 220                                      
CC      I = II+JJ-LAG                                                     
CC      IF( I .GE. II )   GO TO 210                                       
CC      ZZ = C*Z(I)**2                                                    
CC      X = Z(I)*EXP( ZZ )                                                
CC      GO TO 220                                                         
CC  210 I = I - II + 1                                                    
CC      ZZ = C*Y(I)**2                                                    
CC      X = Y(I)*EXP( ZZ )                                                
CC  220 GO TO RETURN, ( 10,20,30 )                                      
C                                                                       
  300 RETURN                                                            
      END                                                               
cc      SUBROUTINE  SBBAYS( X,D,K,N,IPR,MJ1,A,SD )                        
      SUBROUTINE  SBBAYS( X,K,N,IPR,MJ1,A,SD,EK,AIC,IND,C,C1,C2,B,
     *OEIC,ESUM,OMEAN,OM  )                        
C                                                                       
C     THIS SUBROUTINE PRODUCES BAYESIAN MODEL BASED ON ALL SUBSET       
C     REGRESSION MODELS USING THE OUTPUT OF SUBROUTINE REDUCT.          
C       ----------------------------------------------------------------
C       THE FOLLOWING SUBROUTINES ARE DIRECTLY CALLED BY THIS SUBROUTINE
C             RECOEF                                                    
C             SDCOMP                                                    
C             SUBSPC                                                    
C       ----------------------------------------------------------------
C                                                                       
C       INPUTS:                                                         
C          X:     N*(K+1) TRIANGULAR MATRIX,OUTPUT OF SUBROUTINE REDUCT 
C          K:     NUMBER OF REGRESSORS OF THE BAYESIAN MODEL            
C          N:     DATA LENGTH                                           
C          IPR:   =0  TO SUPPRESS THE OUTPUTS                           
C                 =1  TO PRINT OUT FINAL RESULT                         
C                 =2  TO PRINT OUT INTERIM AND FINAL RESULTS            
C          MJ1:   ABSOLUTE DIMENSION OF X                               
C                                                                       
C       OUTPUTS:                                                        
C          A(I) (I=1,K):   REGRESSION COEFFICIENTS OF BAYESIAN MODEL    
C          SD:    RESIDUAL VARIANCE                                     
C                                                                       
      IMPLICIT  REAL * 8  (A-H , O-Z )                                  
cc      DIMENSION  X(MJ1,1) , A(1) , D(1)                                 
cc      DIMENSION  B(50) , G(50)                                          
cc      DIMENSION  IND(50), C(50), C1(50), C2(50), ESUM(50)
      DIMENSION  X(MJ1,1) , A(1) , D(K)                                 
      DIMENSION  B(K) , G(K)
      DIMENSION  IND(K), C(K), C1(K+1), C2(K), ESUM(K+1)
      K1 = K + 1                                                        
      FN = N                                                            
cc      IF( IPR .GE. 2 )     WRITE( 6,3 )                                 
C                                                                       
C          PARTIAL CORRELATIONS COMPUTATION                             
C                                                                       
      SUM = X(K1,K1)**2                                                 
      DO 10   I=1,K                                                     
      J = K1-I                                                          
      SUM = SUM + X(J,K1)**2                                            
      G(J) = DSQRT( SUM )                                               
   10 B(J) = X(J,K1)*X(J,J) / (G(J)*DABS(X(J,J)))                       
C                                                                       
C          PARTIAL CORRELATIONS OF BAYESIAN MODEL COMPUTATION           
C                                                                       
cc      CALL  SUBSPC( B,K,N,IPR,EK )                                      
      CALL  SUBSPC( B,K,N,IPR,EK,IND,C,C1,C2,OEIC,ESUM,OMEAN,OM )
C                                                                       
C          MODIFICATION OF CROSS-PRODUCTS  X(I,K1) (I=1,K)              
C                                                                       
      DO 30   I=1,K                                                     
cc      B(I) = B(I)*X(I,I)*G(I) / DABS(X(I,I))                            
      BB = B(I)*X(I,I)*G(I) / DABS(X(I,I))                            
      D(I) = X(I,K1)                                                    
cc   30 X(I,K1) = B(I)                                                    
   30 X(I,K1) = BB
C                                                                       
C          REGRESSION COEFFICIENTS OF BAYSIAN MODEL                     
C                                                                       
      CALL  RECOEF( X,K,K,MJ1,A )                                       
C                                                                       
      DO 40   I=1,K                                                     
   40 X(I,K1) = D(I)                                                    
C                                                                       
C          RESIDUAL VARIANCE AND AIC                                    
C                                                                       
cc      CALL  SDCOMP( X,A,D,N,K,MJ1,SD )                                  
      CALL  SDCOMP( X,A,N,K,MJ1,SD )                                  
C                                                                       
      IF( IPR .EQ. 0 )     RETURN                                       
      AIC = FN*DLOG(SD) + 2.D0*EK                                       
cc      WRITE( 6,6 )     SD , EK , AIC                                    
      RETURN                                                            
C                                                                       
    3 FORMAT( //1H ,18(1H-),/,' BAYESIAN PROCEDURE',/,1H ,18(1H-) )     
    5 FORMAT( 1H ,10D13.5 )                                             
    6 FORMAT( 1H ,'RESIDUAL VARIANCE',16X,'SD =',D19.8,/,1H ,'EQUIVALENT
     1 NUMBER OF PARAMETERS  EK =',F10.3,/,1H ,32X,'AIC =',F15.3 )      
      E N D                                                             
cc      SUBROUTINE  SETLAG( K )                                           
      SUBROUTINE  SETLAG( K,LAG1,LAG2,LAG3,LAG4,LAG5 )                                           
C                                                                       
C     THIS SUBROUTINE  PREPARES SPECIFICATION OF REGRESSORS (L1(I),L2(I)
C     (I=1,...,K) FOR THE FITTING OF (POLYNOMIAL TYPE) NON-LINEAR MODEL.
C     THE OUTPUTS ARE USED AS THE INPUTS TO SUBROUTINE SETX2.           
C                                                                       
C       INPUTS:                                                         
C          LAG1:    MAXIMUM TIME LAG OF LINEAR TERM                     
C          LAG2:    MAXIMUM TIME LAG OF SQUARED TERM                    
C          LAG3:    MAXIMUM TIME LAG OF QUADRATIC CROSS TERM            
C          LAG4:    MAXIMUM TIME LAG OF CUBIC TERM                      
C          LAG5:    MAXIMUM TIME LAG OF CUBIC CROSS TERM                
C                                                                       
C       OUTPUTS:                                                        
C          K:       NUMBER OF REGRESSORS                                
C          (L1(I),L2(I),L3(I))  (I=1,K):     SPECIFICATION OF REGRESSORS
C                                                                       
C              ......................................................   
C              I-TH REGRESSOR IS DEFINED BY                             
C                   Z(N-L1(I)) * Z(N-L2(I)) * Z(N-L3(I))                
C              WHERE  0-LAG TERM Z(N-0) IS REPLACED BY THE CONSTANT 1.  
C              ......................................................   
C                                                                       
cc      REAL * 4   CSTDMY,SD0DMY
      REAL * 8   CSTDMY,SD0DMY
      COMMON     / BBB /  L1(50),L2(50),L3(50),CSTDMY,SD0DMY
C                                                                       
cc      READ( 5,1 )     LAG1,LAG2,LAG3,LAG4,LAG5                          
cc      WRITE( 6,2 )    LAG1,LAG2,LAG3,LAG4,LAG5                          
      IF(LAG1.LE.0)  GO TO 15                                           
      DO 10  I=1,LAG1                                                   
      L1(I) = I                                                         
      L2(I) = 0                                                         
   10 L3(I) = 0                                                         
   15 K = LAG1                                                          
C                                                                       
      IF(LAG2.LE.0)   GO TO 30                                          
      DO 20  I=1,LAG2                                                   
      K = K+1                                                           
      L1(K) = I                                                         
      L2(K) = I                                                         
   20 L3(K) = 0                                                         
C                                                                       
   30 IF(LAG3.LE.1)   GO TO 50                                          
      LL = LAG3-1                                                       
      DO 40  I=1,LL                                                     
         I1 = I+1                                                       
         DO 40  J=I1,LAG3                                               
         K = K+1                                                        
         L1(K) = I                                                      
         L2(K) = J                                                      
   40    L3(K) = 0                                                      
   50 M  = K                                                            
C                                                                       
      IF(LAG4.LE.0)   GO TO 65                                          
      DO 60  I=1,LAG4                                                   
      K = K+1                                                           
      L1(K) = I                                                         
      L2(K) = I                                                         
   60 L3(K) = I                                                         
   65 CONTINUE                                                          
C                                                                       
      IF(LAG5.LE.1)   GO TO 80                                          
      DO 70  I=1,LAG5                                                   
         DO 70  J=I,LAG5                                                
            DO 70  L=J,LAG5                                             
            IF(I.EQ.J .AND. J.EQ.L)  GO TO 70                           
            K = K+1                                                     
            L1(K) = I                                                   
            L2(K) = J                                                   
            L3(K) = L                                                   
   70 CONTINUE                                                          
C                                                                       
cc   80 WRITE( 6,3 )                                                      
   80 CONTINUE
      IF( LAG1 .EQ. 0 )  GO TO 100                                      
      DO 90  I=1,LAG1                                                   
cc   90 WRITE( 6,4 )     I , L1(I)                                        
   90 CONTINUE
  100 J = LAG1+1                                                        
      IF( LAG2+LAG3 .EQ. 0 )   GO TO 120                                
      DO 110  I=J,M                                                     
cc  110 WRITE( 6,5 )     I , L1(I) , L2(I)                                
  110 CONTINUE
  120 J = M+1                                                           
      IF( LAG4+LAG5 .EQ. 0 )   GO TO 140                                
      DO 130  I=J,K                                                     
cc  130 WRITE( 6,6 )     I ,L1(I) , L2(I) ,L3(I)                          
  130 CONTINUE
  140 CONTINUE                                                          
C                                                                       
      RETURN                                                            
    1 FORMAT( 16I5 )                                                    
    2 FORMAT( /1H ,'LAG1 =',I3,5X,'LAG2 =',I3,5X,'LAG3 =',I3,5X,
     *        'LAG4 =',I3,5X,'LAG5 =',I3 )
    3 FORMAT( 1H ,4X,'M',5X,'REGRESSOR  Y(I,M)' )                       
    4 FORMAT( 1H ,I5,5X,'Z(I-',I2,')' )                                 
    5 FORMAT( 1H ,I5,5X,'Z(I-',I2,') * Z(I-',I2,')' )                   
    6 FORMAT( 1H ,I5,5X,'Z(I-',I2,') * Z(I-',I2,') * Z(I-',I2,')' )     
C                                                                       
      END                                                               
	SUBROUTINE  SETX2( Z,N0,L,K,MJ1,JSW,LAG,X )                       
C                                                                       
C          +----------------------------------------+                   
C          ! MATRIX X SET UP (FOR NON-LINEAR MODEL) !                   
C          +----------------------------------------+                   
C                                                                       
C     THIS SUBROUTINE PREPARES DATA MATRIX X FROM DATA VECTOR Z(I) (I=N0
C     N0+K+LAG) FOR THE FITTING OF NON-LINEAR AUTOREGRESSIVE MODEL.  X I
C     USED AS THE INPUT TO SUBROUTINE HUSHLD.                           
C                                                                       
C       INPUTS:                                                         
C          Z:      ORIGINAL DATA VECTOR                                 
C          N0:     INDEX OF THE END POINT OF DISCARDED FORMER OBSERVATIO
C                  (NEW OBSERVATION STARTS AT N0+LAG+1 AND ENDS AT N0+LA
C          L:      DIMENSION OF THE VECTOR OF NEW OBSERVATIONS          
C          K:      NUMBER OF REGRESSORS                                 
C          MJ1:    ABSOLUTE DIMENSION OF X                              
C          JSW:    =0   TO CONSTRUCT INITIAL L*(K+1) DATA MATRIX        
C                  =1   TO AUGMENT ORIGINAL (K+1)*(K+1) MATRIX X BY AN  
C                       L*(K+1) DATA MATRIX OF ADDITIONAL OBSERVATIONS  
C          LAG:    MAXIMUM TIME LAG                                     
C          KSW:    THIS DUMMY VARIABLE IS NOT REFERENCED IN THIS SUBROUT
C                                                                       
C--  THE FOLLOWING VARIABLE SPECIFICATION IS GIVEN EITHER BY REDLAG OR S
C         (L1(I) , L2(I) , L3(I))  (I=1,K)                              
C                                                                       
C               I-TH REGRESSOR IS DEFINED BY                            
C                    Z(N-L1(I)) * Z(N-L2(I)) * Z(N-L3(I))               
C               WHERE 0-LAG TERM Z(N-0) IS AUTOMATICALLY REPLACED BY CON
C                                                                       
C       OUTPUT:                                                         
C          X:      L*(K+1) MATRIX           IF  JSW = 0                 
C                  (K+1+L)*(K+1) MATRIX     IF  JSW = 1                 
C                                                                       
CC      REAL * 8  X(MJ1,1)                                                
CC      DIMENSION  Z(1)                                                   
cc      REAL * 4   CSTDMY, SD0DMY
      REAL * 8  X(MJ1,1) ,  Z(1) , ZTEM
      REAL * 8  CSTDMY, SD0DMY
      COMMON     / BBB /  L1(50) , L2(50) , L3(50), CSTDMY, SD0DMY  
C                                                                       
      K1 = K + 1                                                        
      I0 = K1*JSW                                                       
      DO  10     I=1,L                                                  
      I1 = I + I0                                                       
      J1 = N0 + LAG + I                                                 
   10 X(I1,K1) = Z(J1)                                                  
C                                                                       
      DO  70     II=1,K                                                 
      LL1 = L1(II)                                                      
      LL2 = L2(II)                                                      
      LL3 = L3(II)                                                      
      DO  60     I=1,L                                                  
      ZTEM = 1.D0                                                       
      I1 = I + I0                                                       
      J1 = N0 + LAG + I                                                 
      IF( LL1 .EQ. 0 )     GO TO 40                                     
      M1 = J1 - LL1                                                     
      ZTEM = ZTEM * Z(M1)                                               
   40 IF( LL2 .EQ. 0 )     GO TO 50                                     
      M2 = J1 - LL2                                                     
      ZTEM = ZTEM * Z(M2)                                               
   50 IF( LL3 .EQ. 0 )     GO TO 60                                     
      M3 = J1 - LL3                                                     
      ZTEM = ZTEM * Z(M3)                                               
   60 X(I1,II) = ZTEM                                                   
   70 CONTINUE                                                          
C                                                                       
      RETURN                                                            
C                                                                       
      E N D                                                             
      SUBROUTINE  SETX4( Z,NO,L,K,MJ1,JSW,LAG,X )                       
C                                                                       
C     THIS SUBROUTINE PREPARES DATA MATRIX X FROM DATA VECTOR Z(I) (I=NO
C     NO+K+L) FOR THE FITTING OF AUTOREGRESSIVE MODEL WITH POLYNOMIAL TY
C     VALUE FUNCTION.  X IS THEN USED AS THE INPUT TO SUBROUTINE HUSHLD.
C       INPUTS:                                                         
C          Z:      ORIGINAL DATA VECTOR                                 
C          NO:     INDEX OF THE END POINT OF DISCARDED FORMER OBSERVATIO
C          L:      DIMENSION OF THE VECTOR OF NEW OBSERVATIONS          
C          K:      NUMBER OF REGRESSORS                                 
C          MJ1:    ABSOLUTE DIMENSION OF X                              
C          JSW:    =0   TO CONSTRUCT INITIAL L*(K+1) DATA MATRIX        
C                  =1   TO AUGMENT ORIGINAL (K+1)*(K+1) MATRIX X BY AN  
C                       L*(K+1) DATA MATRIX OF ADDITIONAL OBSERVATIONS  
C          LAG:    MAXIMUM TIME LAG OF THE MODELS                       
C          N:      DATA LENGTH                                          
C                                                                       
C       OUTPUT:                                                         
C          X:      L*(K+1) MATRIX           IF   JSW = 0                
C                  (K+1+L)*(K+1) MATRIX     IF   JSW = 1                
C                                                                       
C                                                                       
      IMPLICIT  REAL  * 8 (A-H,O-Z)                                     
CC      REAL * 4 Z(1)                                                     
      DIMENSION  Z(1)
      DIMENSION  X(MJ1,1)                                               
cc      COMMON     / AAA /  N , M                                         
      COMMON     / AAA /  N
C                                                                       
C          M:      ORDER OF POLYNOMIAL OF MEAN VALUE FUNCTION           
C                                                                       
      M = K - LAG - 1                                                   
      K1 = K + 1                                                        
      I0 = JSW*K1                                                       
      M1 = M + 1                                                        
      LAG=K-M1                                                          
      BN = 2.D0/(N-LAG-1.D0)                                            
      DO 10  I=1,L                                                      
      Y= BN*(NO+I-1)-1.D0                                               
      XX= 1.D0                                                          
      DO 10 J=1,M1                                                      
      II= I+I0                                                          
      X(II,J) = XX                                                      
   10 XX = XX*Y                                                         
C                                                                       
      DO 20  I=1,L                                                      
      II = I+I0                                                         
      NN = NO+LAG+I                                                     
      X(II,K1) = Z(NN)                                                  
      DO 20   J=1,LAG                                                   
      NN = NN-1                                                         
      JJ = J+M1                                                         
   20 X(II,JJ) = Z(NN)                                                  
C                                                                       
      RETURN                                                            
C                                                                       
      END                                                               
      SUBROUTINE  SETX5( Z,N0,L,K,MJ1,JSW,LAG,X1 )                      
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
C                  =1   TO AUGMENT ORIGINAL (K+1)*(K+1) MATRIX X BY AN  
C                      L*(K+1) DATA MATRIX OF ADDITIONAL OBSERVATIONS   
C          LAG:   MAXIMUM TIME LAG OF THE MODEL                         
C                 =K   CONSTANT VECTOR IS NOT INCLUDED AS A REGRESSOR   
C                 <K   CONSTANT VECTOR IS INCLUDED AS THE FIRST REGRESSO
C                                                                       
C          G:     REGRESSOR AND REGRESSAND SPECIFICATION STATEMENT, LOAD
C                 AT SUBROUTINE REDREG                                  
C                                                                       
C       OUTPUT:                                                         
C          X:      L*(K+1) MATRIX          IF  JSW = 0                  
C                 (K+1+L)*(K+1) MATRIX     IF  JSW = 1                  
C                                                                       
      REAL*8  X1                                                        
      REAL*8  Z , P , Q , U , V , W , Y
      CHARACTER(4) G
      COMMON     / EEE /  G(20,31)                                      
      CHARACTER(4) F
      DIMENSION  Z(1) , F(20)                                           
      DIMENSION  X1(MJ1,1)                                              
      CHARACTER(4) A1, A2, A4,
     *             A5, A6, A7, A8,
     *             T, X
cc      DATA      A1/4H//  / , A2/4H)   / , A4/4H+   / ,                  
cc     *          A5/4H-   / , A6/4H*   / , A7/4H/   / , A8/4H**  / ,     
cc     *          T /4HT   / , X /4HX   /                                 
      DATA      A1/'//  '/ , A2/')   '/ , A4/'+   '/ ,                  
     *          A5/'-   '/ , A6/'*   '/ , A7/'/   '/ , A8/'**  '/ ,     
     *          T /'T   '/ , X /'X   '/                                 
      CHARACTER(4) ASIN, ACOS, CLOG,
     *             BLOG, AEXP, AINV,
     *             SQ, CUBE, ROOT
cc      DATA      ASIN/4HSIN / , ACOS/4HCOS / , CLOG/4HLOG / ,            
cc     *          BLOG/4HLOG1/ , AEXP/4HEXP / , AINV/4HINV / ,            
cc     *          SQ  /4HSQ  / , CUBE/4HCUBE/ , ROOT/4HSQRT/              
      DATA      ASIN/'SIN '/ , ACOS/'COS '/ , CLOG/'LOG '/ ,            
     *          BLOG/'LOG1'/ , AEXP/'EXP '/ , AINV/'INV '/ ,            
     *          SQ  /'SQ  '/ , CUBE/'CUBE'/ , ROOT/'SQRT'/              
      CHARACTER(4) FF, GG
C
      P = 0
      Q = 0
      W = 0
      V = 0
C                                                                       
C     ---  REGRESSOR AND REGRESSAND SPECIFICATION STATEMENT G IS LOADED 
C                                                                       
      K1 = K + 1                                                        
      I0 = 0                                                            
      IF(JSW.EQ.1)   I0 = K1                                            
C                                                                       
      DO 1000  JJ=1,K1                                                  
      DO 5  I=1,20                                                      
    5 F(I) = G(I,JJ)                                                    
      DO 1000  II=1,L                                                   
      I1 = II + I0                                                      
      I = 0                                                             
      Y = 1.D0                                                          
      U = 1.D0                                                          
  500 I = I + 1                                                         
      IF(I.GT.20)  GO TO 10                                             
      FF = F(I)                                                         
C      ENCODE( 4,1,GG )   FF                                             
      GG = FF(1:1)
C                                                                       
      IF(FF.NE.A1)  GO TO 20                                            
   10 Y = Y*U                                                           
      GO TO 1000                                                        
C                                                                       
   20 IF(FF.NE.A2)  GO TO 30                                            
      Y = Y*U                                                           
      GO TO 500                                                         
C                                                                       
   30 IF(GG.NE.T)  GO TO 50                                             
CC      ASSIGN 40 TO ISUB                                                 
CC      GO TO 300
	Q = P                                                             
      P = W                                                             
      W = V                                                             
      V = U                                                             
C
C g77 NOT supports DECODE function and uses special form.
C   40 DECODE( 4,3,FF )   S                                              
      READ(UNIT=FF, FMT=3) S
C
      U = II + N0 + LAG + S                                             
      GO TO 500                                                         
C                                                                       
   50 IF(GG.NE.X)  GO TO 60                                             
CC      ASSIGN 55 TO ISUB                                                 
CC      GO TO 300
	Q = P                                                             
      P = W                                                             
      W = V                                                             
      V = U                                                             
C
C g77 NOT supports DECODE function and uses special form.
C   55 DECODE( 4,2,FF )   LL                                             
      READ(UNIT=FF, FMT=2) LL
C
      IT = II + N0 + LAG + LL                                           
      U = Z(IT)                                                         
      GO TO 500                                                         
C                                                                       
   60 IF(FF.NE.A4)  GO TO 70                                            
      U = V + U                                                         
CC      GO TO 400                                                         
	V = W
	W = P
	P = Q
	GO TO 500
C                                                                       
   70 IF(FF.NE.A5)  GO TO 80                                            
      U = V- U                                                          
CC      GO TO 400                                                         
	V = W
	W = P
	P = Q
	GO TO 500	
C                                                                       
   80 IF(FF.NE.A6)  GO TO 90                                            
      U = V * U                                                         
CC      GO TO 400                                                         
	V = W
	W = P
	P = Q
	GO TO 500
C                                                                       
   90 IF(FF.NE.A7)  GO TO 100                                           
      U = V / U                                                         
CC      GO TO 400                                                         
	V = W
	W = P
	P = Q
	GO TO 500
C                                                                       
  100 IF(FF.NE.A8)  GO TO 120                                           
      U = V**U                                                          
CC      GO TO 400                                                         
	V = W
	W = P
	P = Q
	GO TO 500
C                                                                       
  120 IF(FF.NE.ASIN)  GO TO 130                                         
CC      U = SIN( U )                                                      
      U = DSIN( U )
      GO TO 500                                                         
C                                                                       
  130 IF(FF.NE.ACOS)  GO TO 140                                         
CC      U = COS( U )                                                      
      U = DCOS( U )
      GO TO 500                                                         
C                                                                       
  140 IF(FF.NE.CLOG)  GO TO 150                                         
CC      U = ALOG( U )                                                     
      U = DLOG( U )
      GO TO 500                                                         
C                                                                       
  150 IF(FF.NE.BLOG)  GO TO 160                                         
CC      U = ALOG10( U )                                                   
      U = DLOG10( U )
      GO TO 500                                                         
C                                                                       
  160 IF(FF.NE.AEXP)  GO TO 170                                         
CC      U = EXP( U )                                                      
      U = DEXP( U )
      GO TO 500                                                         
C                                                                       
  170 IF(FF.NE.AINV)  GO TO  180                                        
      U = 1.D0 / U                                                      
      GO TO 500                                                         
C                                                                       
  180 IF(FF.NE.SQ)  GO TO 190                                           
      U = U**2                                                          
      GO TO 500                                                         
C                                                                       
  190 IF(FF.NE.ROOT)  GO TO 220                                         
CC      U = SQRT( U )                                                     
      U = DSQRT( U )
      GO TO 500                                                         
C                                                                       
  220 IF(FF.NE.CUBE)  GO TO 200                                         
      U = U**3                                                          
      GO TO 500                                                         
C                                                                       
  200 CONTINUE                                                          
CC      ASSIGN 210 TO ISUB                                                
CC      GO TO 300
	Q = P                                                             
      P = W                                                             
      W = V                                                             
      V = U                                                             
C
C g77 NOT supports DECODE function and uses special form.
C  210 DECODE( 4,4,FF )   U                                              
      READ(UNIT=FF, FMT=4) U
C
      GO TO 500                                                         
C                                                                       
 1000 X1(I1,JJ) = Y                                                     
      RETURN                                                            
C                                                                       
C-----------------------------------------------------------            
C                        STACK AREA                                     
C-----------------------------------------------------------            
C                                                                       
C          POP UP STACK                                                 
C                                                                       
CC  400 V = W                                                             
CC      W = P                                                             
CC      P = Q                                                             
CC      GO TO 500                                                         
C                                                                       
C          PUSH DOWN STACK                                              
C                                                                       
CC  300 Q = P                                                             
CC      P = W                                                             
CC      W = V                                                             
CC      V = U                                                             
CC      GO TO ISUB, ( 40,55,210 )                                         
C                                                                       
C                                                                       
    1 FORMAT( A1,3X )                                                   
    2 FORMAT( 1X,I3 )                                                   
    3 FORMAT( 1X,F3.0 )                                                 
    4 FORMAT( F4.0 )                                                    
      END                                                               
      SUBROUTINE  SETX6( Z,N0,L,K,MJ1,JSW,LAG,X )                       
C                                                                       
C          +----------------------------------------+                   
C          ! MATRIX X SET UP (FOR NON-LINEAR MODEL) !                   
C          +----------------------------------------+                   
C                                                                       
C     THIS SUBROUTINE PREPARES DATA MATRIX X FROM DATA VECTOR Z(I) (I=N0
C     N0+K+LAG) FOR THE FITTING OF NON-LINEAR AUTOREGRESSIVE MODEL.  X I
C     USED AS THE INPUT TO SUBROUTINE HUSHLD.                           
C                                                                       
C       INPUTS:                                                         
C          Z:      ORIGINAL DATA VECTOR                                 
C          N0:     INDEX OF THE END POINT OF DISCARDED FORMER OBSERVATIO
C                  (NEW OBSERVATION STARTS AT N0+LAG+1 AND ENDS AT N0+LA
C          L:      DIMENSION OF THE VECTOR OF NEW OBSERVATIONS          
C          K:      NUMBER OF REGRESSORS                                 
C          MJ1:    ABSOLUTE DIMENSION OF X                              
C          JSW:    =0   TO CONSTRUCT INITIAL L*(K+1) DATA MATRIX        
C                  =1   TO AUGMENT ORIGINAL (K+1)*(K+1) MATRIX X BY AN  
C                       L*(K+1) DATA MATRIX OF ADDITIONAL OBSERVATIONS  
C          LAG:    MAXIMUM TIME LAG                                     
C          KSW:    THIS DUMMY VARIABLE IS NOT REFERENCED IN THIS SUBROUT
C                                                                       
C--  THE FOLLOWING VARIABLE SPECIFICATION IS GIVEN EITHER BY REDLAG OR S
C         (L1(I) , L2(I) , L3(I))  (I=1,K)                              
C                                                                       
C               I-TH REGRESSOR IS DEFINED BY                            
C                    Z(N-L1(I)) * Z(N-L2(I)) * Z(N-L3(I))               
C               WHERE 0-LAG TERM Z(N-0) IS AUTOMATICALLY REPLACED BY CON
C                                                                       
C       OUTPUT:                                                         
C          X:      L*(K+1) MATRIX           IF  JSW = 0                 
C                  (K+1+L)*(K+1) MATRIX     IF  JSW = 1                 
C                                                                       
      IMPLICIT  REAL * 8  ( Z )
      REAL * 8  X(MJ1,1)                                                
      REAL * 8  C, SD, CONST
      DIMENSION  Z(1)                                                   
      COMMON     / BBB /  L1(50) , L2(50) , L3(50) , SD , CONST         
C                                                                       
      C = -1.0/(SD*CONST)                                               
      K1 = K + 1                                                        
      I0 = K1*JSW                                                       
      DO  10     I=1,L                                                  
      I1 = I + I0                                                       
      J1 = N0 + LAG + I                                                 
   10 X(I1,K1) = Z(J1)                                                  
C                                                                       
      DO  70     II=1,K                                                 
      LL1 = L1(II)                                                      
      LL2 = L2(II)                                                      
      LL3 = L3(II)                                                      
      DO  60     I=1,L                                                  
      ZTEM = 1.0                                                        
      I1 = I + I0                                                       
      J1 = N0 + LAG + I                                                 
      IF( LL1 .EQ. 0 )     GO TO 40                                     
      M1 = J1 - LL1                                                     
      ZZ = C*Z(M1)**2                                                   
CC      ZZ = Z(M1)*EXP(ZZ)                                                
      ZZ = Z(M1)*DEXP(ZZ)
      ZTEM = ZTEM * ZZ                                                  
   40 IF( LL2 .EQ. 0 )     GO TO 50                                     
      M2 = J1 - LL2                                                     
      ZZ = C*Z(M2)**2                                                   
CC      ZZ = Z(M2)*EXP(ZZ)                                                
      ZZ = Z(M2)*DEXP(ZZ)
      ZTEM = ZTEM * ZZ                                                  
   50 IF( LL3 .EQ. 0 )     GO TO 60                                     
      M3 = J1 - LL3                                                     
      ZZ = C*Z(M3)**2                                                   
CC      ZZ = Z(M3)*EXP(ZZ)                                                
      ZZ = Z(M3)*DEXP(ZZ)
      ZTEM = ZTEM * ZZ                                                  
   60 X(I1,II) = ZTEM                                                   
   70 CONTINUE                                                          
C                                                                       
      RETURN                                                            
C                                                                       
      E N D                                                             
      SUBROUTINE  SETX7( Z,N0,L,K,MJ1,JSW,KSW,X )                       
cc      WRITE( 6,3)                                                       
      WRITE( *,3)                                                       
    3 FORMAT(///////1H ,30X,73(1H*),/,                                  
     11H ,30X,23(1H*),5X,'SUBROUTINE  SETX7',5X,23(1H*),//,             
     11H ,35X,'THIS SUBROUTINE IS RESERVED FOR THE OPTIONAL USE.',/,    
     11H ,35X,'USER SHOULD PREPARE THIS SUBROUTINE.',/,                 
     11H ,35X,'FORMAT:',/,                                              
     11H ,40X,'SUBROUTINE SETX7( Z,N0,L,K,MJ1,JSW,KSW,X )',/,           
     11H ,35X,'OTHER PARAMETERS SHOULD BE SENT BY USING COMMON AREA.',/,
     11H ,35X,'SEE THE LISTS OF SETX1,SETX2,SETX3 AND SETX4 AS THE INSTR
     1UCTION.',//,                                                      
     11H ,30X,23(1H*),5X,'SUBROUTINE  SETX6',5X,23(1H*),/,              
     11H ,30X,73(1H*))                                                  
      STOP                                                              
      END                                                               
      SUBROUTINE  SRTMIN( X,N,IX )                                      
C                                                                       
C       THIS SUBROUTINE ARRANGES X(I) (I=1,N) IN ORDER OF INCREASING    
C       MAGNITUDE OF X(I)                                               
C                                                                       
C       INPUTS:                                                         
C          X:   VECTOR                                                  
C          N:   DIMENSION OF THE VECTOR                                 
C       OUTPUTS:                                                        
C          X:   ARRANGED VECTOR                                         
C          IND: INDEX OF ARRANGED VECTOR                                
C                                                                       
      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
      DIMENSION  X(1) , IX(1)                                           
C                                                                       
      NM1 = N - 1                                                       
      DO  30     I=1,N                                                  
   30 IX(I) = I                                                         
      DO  20     II=1,NM1                                               
      XMIN = X(II)                                                      
      MIN = II                                                          
      DO  10     I=II,N                                                 
      IF( XMIN .LT. X(I) )     GO TO 10                                 
      XMIN = X(I)                                                       
      MIN = I                                                           
   10 CONTINUE                                                          
      IF( XMIN .EQ. X(II) )     GO TO 20                                
      XT = X(II)                                                        
      X(II) = X(MIN)                                                    
      X(MIN) = XT                                                       
      IT = IX(II)                                                       
      IX(II) = IX(MIN)                                                  
      IX(MIN) = IT                                                      
   20 CONTINUE                                                          
C                                                                       
      RETURN                                                            
      E N D                                                             
cc      SUBROUTINE  SUBSPC( B,K,N,IPR,EK )                                
      SUBROUTINE  SUBSPC( B,K,N,IPR,EK,IND,C,C1,C2,OEIC,ESUM1,OMEAN,OM )                                
C                                                                       
C       THIS SUBROUTINE PRODUCES BAYESIAN ESTIMATES OF PARTIAL CORRELATI
C       BY CHECKING ALL SUBSET REGRESSION MODELS.                       
C       ----------------------------------------------------------------
C       THE FOLLOWING SUBROUTINES ARE DIRECTLY CALLED BY THIS SUBROUTINE
C             BICOEF                                                    
C             SRTMIN                                                    
C       ----------------------------------------------------------------
C                                                                       
C         INPUTS:                                                       
C           B:   LEAST SQUARES ESTIMATES OF PARTIAL CORRELATIONS        
C           K:   DIMENSION OF VECTOR A                                  
C           N:   NUMBER OF OBSERVATIONS USED FOR THE ESTIMATION OF A(I) 
C           IPR: =0  TO SUPPRESS THE OUTPUTS                            
C                >0  TO PRINT OUT THE OUTPUTS                           
C                                                                       
C         OUTPUTS:                                                      
C           B(I) (I=1,K):   BAYESIAN ESTIMATES OF PARTIAL CORRELATIONS  
C           EK:   EQUIVALENT NUMBER OF FREE PARAMETERS IN THE BAYESIAN M
C                                                                       
      IMPLICIT  REAL * 8 ( A-H , O-Z )                                  
cc      DIMENSION  B(1) , C(50) , D(50,50)                                
cc      DIMENSION  IND(50) , KND(50) , ESUM(50)                           
cc      DIMENSION  C1(50), C2(50), ESUM1(50)
      DIMENSION  B(1) , C(K) , D(K+1,K+1)                                
      DIMENSION  IND(K) , KND(K+1) , ESUM(K+1)                           
      DIMENSION  C1(K+1), C2(K), ESUM1(K+1)
C                                                                       
      CC = 1.D0 + DLOG(2.D0)                                            
      K1 = K + 1                                                        
      DN = N                                                            
      DO 10   I=1,K1                                                    
      ESUM(I) = 0.D0                                                    
      DO 10   J=1,K1                                                    
   10 D(I,J) = 0.D0                                                     
C                                                                       
C          SQUARE OF PARTIAL CORRELATIONS ( NORMALISED BY MULTIPLYING N 
C                                                                       
      DO 20   I=1,K                                                     
   20 C(I) = B(I)*B(I)*DN                                               
C                                                                       
C          ARRANGEMENT OF C(I) IN ORDER OF INCREASING MAGNITUDE         
C                                                                       
      CALL  SRTMIN( C,K,IND )                                           
cc      IF( IPR .LE. 1 )     GO TO 60                                     
cc      WRITE( 6,7 )                                                      
cc      DO  50     I=1,K                                                  
cc   50 WRITE( 6,6 )     I , IND(I) , C(I)                                
cc   60 CONTINUE                                                          
C                                                                       
C          FIND THE MINIMUM OF EIC                                      
C                                                                       
      OEIC = CC*K                                                       
      SUM = 0.D0                                                        
      DO 30   I=1,K                                                     
      SUM = SUM + C(I)                                                  
      EIC = SUM + CC*(K-I)                                              
   30 IF( OEIC .GT. EIC )   OEIC = EIC                                  
cc      WRITE( 6,604 )   OEIC                                             
C                                                                       
C--------  COMPUTATION OF EIC'S OF WHOLE SUBSET REGRESSION MODELS  -----
C                                                                       
C          INITIAL SETTING                                              
C                                                                       
      DO 40   I=1,K                                                     
   40 KND(I) = 0                                                        
      KND(K1) = 1                                                       
      SUM = 0.D0                                                        
      SUMC = 0.D0
      M = K                                                             
      IP = 0                                                            
      IQ = 0                                                            
C                                                                       
  100 CONTINUE                                                          
C                                                                       
C          -----  SPECIFICATION OF NEXT SUBSET  -----                   
               DO 110   I=1,K                                           
               IF( KND(I) .EQ. 0 )   GO TO 110                          
               KND(I) = 0                                               
               GO TO 120                                                
  110          KND(I) = 1                                               
  120          CONTINUE                                                 
C          ------------------------------------------                   
C                                                                       
  130 CONTINUE                                                          
      IF( IP .GT. K )   GO TO 200                                       
      IF( KND(IP+1) .EQ. 0 )   GO TO 140                                
C                                                                       
      IF( IQ .EQ. 0 )   GO TO 165                                       
      IF( KND(IQ) .EQ. 1 )   GO TO 150                                  
      IQ = IQ-1                                                         
      SUMC = SUMC + C(IQ+1)                                             
C                                                                       
      IF( SUMC + CC*(K-IP+IQ) .GT. OEIC + 40.D0 )   GO TO 180           
      GO TO 150                                                         
C                                                                       
  140 IP = IP+1                                                         
      IQ = IP-1                                                         
      SUMC = C(IP)                                                      
      IF( SUMC + CC .GT. OEIC + 40.D0 )   GO TO 200                     
C                                                                       
  150 M = K-IP+IQ                                                       
      SUM = SUMC                                                        
      IF( IQ .EQ. 0 )   GO TO 165                                       
      DO 160   I=1,IQ                                                   
      IF( KND(I) .EQ. 1 )   GO TO 160                                   
      M = M-1                                                           
      SUM = SUM + C(I)                                                  
  160 CONTINUE                                                          
  165 CONTINUE                                                          
      EIC = SUM + CC*M - OEIC                                           
      IF( EIC .GT. 40.D0 )   GO TO 100                                  
      EXIC = DEXP( -0.5D0*EIC )                                         
      ESUM(M+1) = ESUM(M+1) + EXIC                                      
      DO 170   I=1,K                                                    
  170 IF( KND(I) .EQ. 1 )   D(I,M+1) = D(I,M+1) + EXIC                  
      GO TO 100                                                         
C         --------------------------------------------                  
  180          DO 190   I=1,IP                                          
  190          KND(I) = 1                                               
               KND(IP+1) = 0                                            
               IP = IP+1                                                
               IQ = IP-1                                                
               GO TO 130                                                
C         ---------------------------------------------                 
C                                                                       
C--------------------------  WHOLE SUBSETS CHECKED  --------------------
C                                                                       
  200 CONTINUE                                                          
cc      IF( IPR .GE. 2 )     WRITE( 6,8 )                                 
cc      IF( IPR .GE. 2 )     WRITE( 6,607 )     (ESUM(I),I=1,K1)          
      DO 201 I=1,K1
  201 ESUM1(I) = ESUM(I)
C                                                                       
C          MEAN OF NUMBER OF PARAMETERS                                 
C                                                                       
      OSUM = 0.D0                                                       
      SUM = ESUM(1)                                                     
      DO 210   I=1,K                                                    
      SUM = SUM + ESUM(I+1)                                             
  210 OSUM = OSUM + I*ESUM(I+1)                                         
      OMEAN = OSUM / SUM                                                
      OM = OMEAN / K                                                    
cc      IF( IPR .GE. 2 )     WRITE( 6,608 )     OMEAN , OM                
C                                                                       
C       --  BINOMIAL TYPE DAMPER  --                                    
C                                                                       
      DO 220   I=1,K1                                                   
      J = I-1                                                           
      KMJ = K-J                                                         
cc  220 C(I) = BICOEF(K,J)*(OM**J)*((1.D0-OM)**KMJ)                       
cc      C(1) = 1.D0 / (1.D0 + K)                                          
  220 C1(I) = BICOEF(K,J)*(OM**J)*((1.D0-OM)**KMJ)                       
      C1(1) = 1.D0 / (1.D0 + K)                                          
      DO 221  I=1,K                                                     
cc  221 C(I+1) = C(I) * I / (1.D0 + K - I)                                
  221 C1(I+1) = C1(I) * I / (1.D0 + K - I)                                
cc      IF( IPR .GE. 2 )     WRITE( 6,609 )                               
cc      IF( IPR .GE. 2 )     WRITE( 6,607 )     (C(I),I=1,K1)             
C                                                                       
      DO 230   I=1,K1                                                   
cc  230 ESUM(I) = ESUM(I)*C(I)                                            
  230 ESUM(I) = ESUM(I)*C1(I)                                            
C                                                                       
      SUM = 0.D0                                                        
      DO 240   I=1,K1                                                   
  240 SUM = SUM + ESUM(I)                                               
C                                                                       
      DO 250   J=1,K1                                                   
      DO 250   I=1,K                                                    
cc  250 D(I,J) = D(I,J)*C(J) / SUM                                        
  250 D(I,J) = D(I,J)*C1(J) / SUM                                        
C                                                                       
C          WEIGHTS OF PARTIAL CORRELATIONS                              
C                                                                       
      DO 260   I=1,K                                                    
cc  260 C(I) = 0.D0                                                       
  260 C2(I) = 0.D0                                                       
      DO 270   I=1,K                                                    
      DO 270   J=1,K1                                                   
cc  270 C(I) = C(I) + D(I,J)                                              
  270 C2(I) = C2(I) + D(I,J)                                              
cc      IF( IPR .GE. 2 )     WRITE( 6,603 )                               
cc      IF( IPR .GE. 2 )     WRITE( 6,607 )     (C(I),I=1,K)              
C                                                                       
C          AVERAGING AND REARRANGEMENT OF PARTIAL CORRELATIONS          
C                                                                       
      EK = 1.D0                                                         
      DO 280   I=1,K                                                    
      J = IND(I)                                                        
cc      B(J) = B(J)*C(I)                                                  
cc  280 EK = EK + C(I)**2                                                 
      B(J) = B(J)*C2(I)                                                  
  280 EK = EK + C2(I)**2                                                 
cc      IF( IPR .LE. 1 )     RETURN                                       
cc      WRITE( 6,602 )                                                    
cc      DO  290     I=1,K                                                 
cc  290 WRITE( 6,609 )     I , B(I)                                       
C                                                                       
      RETURN                                                            
C                                                                       
    3 FORMAT( 1H ,20I3 )                                                
    4 FORMAT( 1H ,3D20.10,2I5 )                                         
    6 FORMAT( 1H ,2I7,F13.5 )                                           
    7 FORMAT( 1H ,6X,'I IND(I)',4X,'N*B(I)**2' )                        
    8 FORMAT( 1H ,'ESUM(I) (I=1,M+1)' )                                 
    9 FORMAT( 1H ,'***  BINOMIAL TYPE  ***' )                           
  602 FORMAT( 1H ,'PARTIAL CORRELATIONS OF THE BAYESIAN MODEL',/,1H ,6X,
     1'I',9X,'A(I)' )                                                   
  603 FORMAT( 1H ,'FINAL BAYESIAN WEIGHTS OF PARTIAL CORRELATIONS' )    
  604 FORMAT( 1H ,'  OAIC =',D13.5 )                                    
  606 FORMAT( 1H ,'OSUM =',D20.10,5X,'SUM =',D20.10,5X,'COD =',D20.10 ) 
  607 FORMAT( 1H ,10D13.5 )                                             
  608 FORMAT( 1H ,'OMEAN =',D15.8,5X,'OM =',D15.8 )                     
  609 FORMAT( 1H ,I7,F13.5 )                                            
      END                                                               
