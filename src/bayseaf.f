C Reviced	M.S C85-02-19-16:46:06 BAYSEAA PAIR                                         
CC      PROGRAM BAYSEA
      SUBROUTINE BAYSEAF(Y,NDATA,FOCAST,CDATA,DMOI,TREND,SEASON,TDCMP,
     *                    IRREG,ADJUST,EST,PSDS,PSDT,AVABIC,
     *                    IPARA,PARA,ARFT,ARFS,ARFN,IART,IARS,IARN)
C                                                                      
      INCLUDE 'timsac_f.h'
C
C     ---      --       -                                               
C     BAYESIAN SEASONAL ADJUSTMENT PROCEDURE                            
C                                                                       
C     THIS IS VERSION(3/1/85) OF BAYSEA WHICH WAS ORIGINALLY            
C     PUBLISHED IN                                                      
C                                                                       
C       AKAIKE,H. AND ISHIGURO,M. (1980)                                
C         BAYSEA, A BAYESIAN SEASONAL ADJUSTMENT PROGRAM.               
C         COMPUTER SCIENCE MONOGRAPHS, NO.13,                           
C         THE INSTITUTE OF STATISTICAL MATHEMATICS, TOKYO.              
C                                                                       
C     THIS VERSION WAS DESIGNED AND PROGRAMMED BY                       
C     HIROTUGU AKAIKE AND MAKIO ISHIGURO, THE INSTITUTE OF STATISTICAL  
C     MATHEMATICS, 4-6-7 MINAMI-AZABU, MINATO-KU, TOKYO, 106, JAPAN.    
C     THE SUBROUTINES FOR OUTLIER CORRECTION WERE PREPARED BY           
C     GENSHIRO KITAGAWA.                                                
C                                                                       
C     THIS PROGRAM REALIZES A DECOMPOSITION OF TIME SERIES Y            
C     INTO THE FORM                                                     
C     Y(I) = T(I) +S(I)+I(I)+TDC(I)+OCF(I)                              
C     WHERE  T(I)=TREND  S(I)=SEASONAL  I(I)=IRREGULAR                  
C            TDC(I)=TRADING DAY COMPONENT     AND                       
C            OCF(I)=OUTLIER CORRECTION FACTOR                           
C                                                                       
C     THE PROCEDURE IS BASED ON A BAYESIAN MODEL AND ITS                
C     PERFORMANCE IS CONTROLLED BY THE SELECTION OF THE PARAMETERS OF   
C     THE PRIOR DISTRIBUTION.  THE CONSTRUCTION OF THE BASIC MODEL IS   
C     DISCUSSED IN THE FOLLOWING PAPERS:                                
C                                                                       
C       AKAIKE, H. (1980) LIKELIHOOD AND THE BAYES PROCEDURE.           
C         BAYESIAN STATISTICS, J.M.BERNARDO, M.H.DE GROOT, D.V.LINDLEY  
C         AND A.F.M.SMITH, EDS., UNIVERSITY PRESS, VALENCIA,            
C         SPAIN, 143-166.                                               
C       AKAIKE, H. (1980) SEASONAL ADJUSTMENT BY A BAYESIAN MODELING.   
C         JOURNAL OF TIME SERIES ANALYSIS, 1, 1-13.                     
C       AKAIKE, H. AND ISHIGURO, M. (1980) TREND ESTIMATION WITH        
C         MISSING OBSERVATIONS. ANNALS OF THE INSTITUTE OF STATISTICAL  
C         MATHEMATICS, 32,B, 481-488.                                   
C       ISHIGURO, M. AND AKAIKE, H. (1980) A BAYESIAN APPROACH TO       
C         THE TRADING-DAY ADJUSTMENT OF MONTHLY DATA.                   
C         TIME SERIES ANALYSIS, O.D.ANDERSON AND M.R.PERRYMAN, EDS.,    
C         NORTH-HOLLAND, AMSTERDAM, 213-226.                            
C       ISHIGURO,M. (1984) COMPUTATIONALLY EFFICIENT IMPLEMENTATION     
C         OF A BYAESIAN SEASONAL ADJUSTMENT PROCEDURE.                  
C         JOURNAL OF TIME SERIES ANALYSIS (TO APPEAR).                  
C       KITAGAWA,G. AND AKAIKE,H. (1982) A QUASI BAYESIAN APPROACH      
C         TO OUTLIER DETECTION. ANNALS OF THE INSTITUTE OF              
C         STATISTICAL MATHEMATICS, 34,B, 389-398.                       
C                                                                       
C     THE PRIOR DISTRIBUTION CONTROLS THE SMOOTHNESS OF THE TREND AND   
C     SEASONAL COMPONENTS BY ASSUMING LOW ORDER GAUSSIAN AR-MODELS      
C     FOR SOME DIFFERENCES OF THESE COMPONENTS. THE CHOICE OF           
C     THE VARIANCE OF THE GAUSSIAN DISTRIBUTION IS                      
C     REALIZED BY MAXIMIZING THE LOG LIKELIHOOD OF THE BAYESIAN MODEL.  
C                                                                       
C     FOR THE PURPOSE OF COMPARISON OF MODELS WITH DIFFERENT STRUCTURES 
C     THE CRITERION ABIC IS DEFINED BY                                  
C                                                                       
C          ABIC = (-2)LOG MAXIMUM LIKELIHOOD OF THE MODEL.              
C                                                                       
C     SMALLER VALUE OF ABIC REPRESENTS BETTER FIT.                      
C     FOR THE COMPARISON OF OVERALL PERFORMANCES OF VARIOUS             
C     MODELS THE AVERAGED ABIC (AVABIC) IS USEFUL.                      
C                                                                       
C     THIS PROGRAM REQUIRES THE FOLLOWING PARAMETERS.  WITHOUT FURTHER  
C     SPECIFICATION BY THE PROGRAM USER, THEY ARE RESPECTIVELY SET EQUAL
C     TO THE VALUES GIVEN IN THE PARENTHESES AT THE ENDS OF THEIR       
C     DESCRIPTIONS.                                                     
C                                                                       
C     PARAMETERS:                                                       
C       LOGT:   LOG-ADDITIVE-MODEL OPTION (0)                           
C        IF( LOGT.EQ.0) ADDITIVE MODEL                                  
C        IF( LOGT.EQ.1) LOG-ADDITIVE-MODEL                              
C       MT:     INPUT DEVICE SPECIFICATION (5)                          
C       RLIM:   OUTLIER LIMIT (0.0)                                     
C        IF( RLIM.GT.0.0) ANY DATA WHOSE VALUE IS GREATER               
C                         THAN RLIM IS TREATED AS A MISSING OBSERVATION 
C        IF( RLIM.LE.0.0) NO MISSING OBSERVATIONS                       
C       PERIOD: NUMBER OF SEASONALS WITHIN A PERIOD (12)                
C       SPAN:   NUMBER OF PERIODS TO BE PROCESSED AT ONE TIME (4)       
C       SHIFT:  NUMBER OF PERIODS TO BE SHIFTED                         
C               TO DEFINE THE NEW SPAN OF DATA (1)                      
C       FOCAST: LENGTH OF FORECAST AT THE END OF DATA (0)               
C       ORDER:  ORDER OF DIFFERENCING OF TREND (2)                      
C       ARFT(I): I-TH PARCOR OF DIFFERENCED TREND (I LESS THAN 4)       
C       (0.D0)                                                          
C       SORDER: ORDER OF DIFFERENCING OF SEASONAL (1)                   
C       ARFS(I): I*PERIOD-TH PARCOR OF DIFFERENCD SEASONAL              
C       (I LESS THAN 4) (0.D0)                                          
C       ARFN(I): I-TH PARCOR OF DIFFERENCD SEASONAL                     
C       (I LESS THAN 4) (0.DO)                                          
C       (FOR THE INITIAL CHOICES OF ARFT,ARFS,AND ARFN, USE             
C        PARCOR (PARTIAL AUTOCORRELATION) OUTPUTS FOR THE               
C        DEFAULT OPTION)                                                
C               (SORDER SHOULD NOT BE GREATER THAN SPAN.)               
C       RIGID:  CONTROLS THE RIGIDITY OF THE SEASONAL COMPONENT, MORE   
C               RIGID SEASONAL WITH LARGER RIGID (1.0)                  
C               (TO BE ADJUSTED ONLY AFTER THE SELECTION OF ORDER AND   
C                SORDER)                                                
C       YEAR:   TRADING-DAY ADJUSTMENT OPTION  (0)                      
C        IF(YEAR .EQ. 0) WITHOUT TRADING-DAY ADJUSTMENT                 
C        IF(YEAR .NE. 0) WITH TRADING-DAY ADJUSTMENT                    
C         NOTE: THE SERIES IS SUPPOSED TO START AT THIS 'YEAR'          
C       MONTH:  NUMBER OF THE MONTH IN WHICH THE SERIES STARTS (1)      
C        IF(YEAR .EQ. 0) THIS PARAMETER IS IGNORED                      
C       WTRD:   CONTROLS THE ADAPTIVITY OF THE TRADING-DAY-COMPONENT    
C               MORE DATA-ADAPTIVE WITH SMALLER WTRD (1.0)              
C       IOUTD: OUTLIER CORRECTION OPTION (0)                            
C        IF(IOUTD .EQ. 0) WITHOUT OUTLIER DETECTION                     
C        IF(IOUTD .EQ. 1) WITH OUTLIER DETECTION BY MARGINAL PROBABILITY
C        IF(IOUTD .EQ. 2) WITH OUTLIER DETECTION BY MODEL SELECTION     
C         (NOTE: OUTLIER DETECTION IS                                   
C                EXPENSIVE. DECLARING ABNORMAL OBSERVATIONS AS MISSING  
C                BY RLIM IS MUCH CHEAPER)                               
C       SPEC:  SPECTRUM ESTIMATION OPTION (1)                           
C        IF(SPEC .EQ. 0) NO SPECTRUM                                    
C        IF(SPEC .EQ. 1) SPECTRA OF IRREGULAR AND                       
C                         DIFFERENCED ADJUSTED                          
C         (IF THERE ARE MISSING OBSERVATIONS, SPECTRA OF THE            
C          INTERPOLATED DATA ARE GIVEN)                                 
C       PUNCH:  CARD OUTPUT CONTROL (0)                                 
C        IF(PUNCH .EQ. 0) OUTPUT SUPRESSED                              
C        IF(PUNCH .EQ. 1) CARD OUTPUT AVAILABLE                         
C                                                                       
C      --- PARAMETERS BELOW THIS LINE ARE SELDOM TO BE MODIFIED ---     
C                                                                       
C       ZERSUM: CONTROLS THE SUM OF THE SEASONALS WITHIN A PERIOD,      
C               CLOSER TO ZERO WITH LARGER ZERSUM (1.0)                 
C       DELTA:  CONTROLS THE LEAP YEAR EFFECT (7.0)                     
C       ALPHA:  CONTROLS PRIOR VARIANCE OF INITIAL TREND (0.01D0)       
C       BETA:   CONTROLS PRIOR VARIANCE OF INITIAL SEASONAL             
C               (0.01D0)                                                
C       GAMMA:  CONTROLS PRIOR VARIANCE OF INITIAL SUM OF SEASONAL      
C               (0.01D0)                                                
C                                                                       
C     THE FOLLOWING PROCEDURE OF PARAMETER MODIFICATION IS RECOMMENDED  
C     FOR USUAL APPLICATIONS:                                           
C     FIRST TRY THE COMBINATION(ORDER=2,SORDER=1). IF NECESSARY, TRY    
C     (ORDER=2,SORDER=2). THEN TRY A REDUCED VALUE OF 'RIGID' AND CHECK 
C     ABIC. WHEN 'RIGID' IS SUITABLY CHOSEN, USUALLY SHALLOW TROUGHS    
C     APPEAR AT FRIQUENCIES MARKED BY '---X' ON THE CHARTS OF SPECTRA.  
C                                                                       
C     TO MODIFY THE PARAMETERS FOLLOW THE FOLLOWING EXAMPLE:            
C       IF THE USER WANTS TO SET SORDER=2,RIGID=0.5 AND SO ON, AND IF   
C       MT IS CARD READER, PLACE THE CARDS WHICH CONTAIN THE FOLLOWING  
C       THREE TYPES OF STATEMENTS ON TOP OF THE INPUT DATA( PUNCH ONE   
C       SPACE AT THE FIRST COLUMN OF EACH CARD) :                       
C  &PARAM  (ON THE FIRST CARD)                                          
C  SORDER=2, RIGID=0.5, AND SO ON, (ON THE SECOND AND LATER CARDS)      
C  &END  (ON THE LAST CARD)                                             
C     NOTE THE COMMA AT THE END OF EACH PARAMETER SPECIFICATION.        
C                                                                       
C     INPUT DATA:                                                       
C                                                                       
C     THE FOLLOWING DATA SET SHOULD BE FED THROUGH THE INPUT DEVICE     
C     SPECIFIED BY MT IN THE FORMATS SHOWN IN THE PARENTHSES.           
C                                                                       
C     TITLE (20A4):   TITLE OF THE DATA                                 
C     NDATA (I5):     DATA LENGTH                                       
C     FORMAT SPECIFICATION OF DATA (20A4): FOR EXAMPLE, (4D20.10)       
C     DATA:       Y(I);I=1,NDATA                                        
C                                                                       
C     NUMERICAL OUTPUTS:                                                
C       OCF:    OUTLIER CORRECTION FACTOR                               
C       TREND:  ESTIMATED TREND                                         
C       SEASONAL: ESTIMATED SEASONAL                                    
C       TDCMP: ESTIMATED TRADING-DAY COMPONENT                          
C       IRREGULAR = 0.0        (IF OBSERVATION IS MISSING)              
C                 = ORIGINAL DATA - TREND - SEASONAL - TDCMP - OCF      
C                              (OTHERWISE)                              
C       ADJUSTED = TREND + IRREGULAR                                    
C       SMOOTHED = TREND + SEASONAL + TDCMP                             
C                                                                       
C     GRAPHICAL OUTPUTS:                                                
C       ORIGINAL DATA                                                   
C       OCF                                                             
C       TREND                                                           
C       SEASONAL                                                        
C       IRREGULAR                                                       
C       ADJUSTED                                                        
C       SMOOTHED                                                        
C       TDCMP                                                           
C       MISSING OBSERVATION INTERPOLATED DATA                           
C       SPECTRA OF IRREGULAR AND DIFFERENCED ADJUSTED                   
C                                                                       
C     WORKING AREA REQUIRED BY THIS PROGRAM:                            
C     VECTORS:                                                          
C       Y(IY)                                                           
C       TREND(IRSLT)                                                    
C       SEASON(IRSLT)                                                   
C       EST(IRSLT)                                                      
C       IRREG(IRSLT)                                                    
C       TDCMP(IRSLT)                                                    
C       ADJUST(IRSLT)                                                   
C       CDATA(IRSLT)                                                    
C       DMOI(IRSLT)                                                     
C       FTRN(IA)                                                        
C       FSEA(IA)                                                        
C       YS(IA)                                                          
C       YS1(IA)                                                         
C       YO(IA)                                                          
C       DC(MDC)                                                         
C       H2(8*IA)                                                        
C       WEEK(7*IA)                                                      
C   WHERE                                                               
C    IY = MAXIMUM DATA LENGTH                                           
C    IRSLT = MAXIMUM OUTPUT LENGTH                                      
C    IA = ( MAXIMUM DATA LENGTH WITHIN A SPAN)                          
C                                                                       
C   STRUCTURE OF THE PROGRAM                                            
C                                                                       
C   ****************00002150
C   *              *00002160
C   *    BAYSEA    *00002170
C   *              *00002180
C   ********I*******00002190
C           I                                                           
C           I-ARCOEF---PARTAR                                           
C           I                                                           
C           I----------------------------DSQRT                          
C           I                                                           
C           I-DATAR ------------DLOG                                    
C           I                                                           
C           I-------------------DLOG                                    
C           I                                                           
C           I-------------------COPY                                    
C           I                                                           
C           I-CALEND------------MOD                                     
C           I                                                           
C           I-SUBSEA-+-------------------DSQRT                          
C           I        I                                                  
C           I        I-SETDC -+-SETD  ---MIN0                           
C           I        I        I                                         
C           I        I        I-INIT                                    
C           I        I        I                                         
C           I        I        +-EXHSLD                                  
C           I        I                                                  
C           I        I-------------------DABS                           
C           I        I                                                  
C           I        I----------DLOG                                    
C           I        I                                                  
C           I        I-SETX  -+-HUSHLD-+-DABS                           
C           I        I        I        I                                
C           I        I        I        +-DSQRT                          
C           I        I        I                                         
C           I        I        I----------DABS                           
C           I        I        I                                         
C           I        I        I-DLOG                                    
C           I        I        I                                         
C           I        I        I-SETD  ---MIN0                           
C           I        I        I                                         
C           I        I        I-INIT                                    
C           I        I        I                                         
C           I        I        I-EXHSLD                                  
C           I        I        I                                         
C           I        I        +-MOD                                     
C           I        I                                                  
C           I        I-SOLVE                                            
C           I        I                                                  
C           I        I-DECODE-+-CLEAR                                   
C           I        I        I                                         
C           I        I        I----------DSQRT                          
C           I        I        I                                         
C           I        I        I-COPY                                    
C           I        I        I                                         
C           I        I        I-PRDCT                                   
C           I        I        I                                         
C           I        I        I-ADD                                     
C           I        I        I                                         
C           I        I        +-SBTRCT                                  
C           I        I                                                  
C           I        +-DEXP                                             
C           I                                                           
C           I-------------------SBTRCT                                  
C           I                                                           
C           I-OUTLIR-+-SRTMIN                                           
C           I        I                                                  
C           I        I----------DLOG                                    
C           I        I                                                  
C           I        I----------DFLOAT                                  
C           I        I                                                  
C           I        I----------BINARY                                  
C           I        I                                                  
C           I        I-LKOUT1-+-DFLOAT                                  
C           I        I        I                                         
C           I        I        I-DLOG                                    
C           I        I        I                                         
C           I        I        I-POOLAV                                  
C           I        I        I                                         
C           I        I        I----------DSQRT                          
C           I        I        I                                         
C           I        I        +-PERMUT---ISORT                          
C           I        I                                                  
C           I        I-DEXP                                             
C           I        I                                                  
C           I        I-PRPOST---BINARY                                  
C           I        I                                                  
C           I        +-MODIFY---BINARY                                  
C           I                                                           
C           I-------------------ADD                                     
C           I                                                           
C           I----------DEXP                                             
C           I                                                           
C           I-DFR1                                                      
C           I                                                           
C           I-GRAPH -+-LOG10                                            
C           I        I                                                  
C           I        I-SQRT                                             
C           I        I                                                  
C           I        I-------------------ABS                            
C           I        I                                                  
C           I        I----------MOD                                     
C           I        I                                                  
C           I        +-DEXP                                             
C           I                                                           
C           +-SPGRH -+-------------------MIN0                           
C                    I                                                  
C                    I-SAUTCO-+-SMEADL---SUMF                           
C                    I        I                                         
C                    I        I-CROSCO---DBLE                           
C                    I        I                                         
C                    I        +-CORNOM---DSQRT                          
C                    I                                                  
C                    I-------------------DSQRT                          
C                    I                                                  
C                    I-SICP2 -+-DLOG                                    
C                    I        I                                         
C                    I        +----------DSQRT                          
C                    I                                                  
C                    +-SNRASP-+-FOUGER-+-DCOS                           
C                             I        I                                
C                             I        +-DSIN                           
C                             I                                         
C                             I-SUBVCP                                  
C                             I                                         
C                             I-DLOG10                                  
C                             I                                         
C                             I-DSP3  -+-AMAX                           
C                             I        I                                
C                             I        I-ABS                            
C                             I        I                                
C                             I        +-AMIN                           
C                             I                                         
C                             +-MOD                                     
C                                                                       
cc      !DEC$ ATTRIBUTES DLLEXPORT::BAYSEAF
C
       IMPLICIT  REAL*8( A-H,O-Z )                                      
      INTEGER*4  ORDER, SORDER, PERIOD, SPAN, OVLAP, FOCAST, HEAD, SHIFT
cc     * ,TAIL,PUNCH,YEAR,SPEC                                            
     * ,TAIL,YEAR,SPEC                                            
      REAL*8 IRREG, IRREG0                                              
cc      DIMENSION  FTRN(200),FSEA(200),PSDS(500),PSDT(500),PSDS0(500),    
cc     *           PSDT0(500)                                             
cc      DIMENSION  SEASON(500), TREND(500), EST(500), ADJUST(500),        
cc     * IRREG(500),WEEK(2000),TDCMP(500),TDCMP0(500),DMOI(500)           
cc     *       ,SEAS0(500), TREND0(500), EST0(500), ADJ0(500), IRREG0(500)
cc      DIMENSION  DC(40000),ARFT(3),ARFS(3),ARFN(3), H2(4000)            
cc      DIMENSION  Y(500),YS(500),CDATA(500),YS1(500),YO(500)             
cc      DIMENSION   DADJ(500),H(20000),F(1000)
cc      DIMENSION  DC(40000),H2(4000),H(20000),F(1000),WEEK(2000)
cc      DIMENSION  YS(500),YS1(500),YO(500)
      DIMENSION  Y(NDATA),CDATA(NDATA),DMOI(NDATA),
     *            TREND(NDATA+FOCAST),SEASON(NDATA+FOCAST),
     *            TDCMP(NDATA+FOCAST),IRREG(NDATA),
     *            ADJUST(NDATA),EST(NDATA+FOCAST),
     *            PSDS(NDATA+FOCAST),PSDT(NDATA+FOCAST),
     *            IPARA(12),PARA(8),ARFT(3),ARFS(3),ARFN(3)
c
      DIMENSION  TREND0(NDATA+FOCAST),SEAS0(NDATA+FOCAST),
     *            TDCMP0(NDATA+FOCAST),IRREG0(NDATA+FOCAST),
     *            ADJ0(NDATA+FOCAST),EST0(NDATA+FOCAST),
     *            PSDS0(NDATA+FOCAST),PSDT0(NDATA+FOCAST)
      DIMENSION  FTRN(IPARA(4)+3),FSEA((IPARA(5)+3)*IPARA(1)+3)
      DIMENSION  F(NDATA+FOCAST+1),WEEK(7,NDATA+FOCAST)
      DIMENSION  YS(NDATA+FOCAST),YS1(NDATA+FOCAST),YO(NDATA+FOCAST)
cc      CHARACTER*80   TITLE
cc      COMMON /ILOGT/ LOGT,ISHRNK,PUNCH,IOUTD,ROUT                       
cc      COMMON /IDATA/ PERIOD,ORDER,SORDER,YEAR,NDAY,IFIX                 
cc      COMMON /RDATA/ ALPHA,BETA,GAMMA,ZER,SMTH,SMTH2,DD,WTRD,DELTA      
cc      NAMELIST /PARAM/  MT, PERIOD, RLIM, SPAN, SHIFT, FOCAST,          
cc     *                  RIGID, ORDER, SORDER, ZERSUM,  LOGT, PUNCH      
cc     * ,YEAR,MONTH,WTRD,DELTA,SPEC,IOUTD,IFIX,ARFT,ARFS,ARFN,IART       
cc     * ,IARS,IARN
C
C	FILE NAME INPUT 
C
cc 	CHARACTER(100)  IFLNAM,OFLNAM,MFLNAM
cc         CALL FLNAM3(IFLNAM,OFLNAM,MFLNAM,NFL)
cc            NFL = 4
cc 	IF (NFL.EQ.0) GO TO 999
cc 	IF ((NFL.EQ.2) .OR. (NFL.EQ.4)) THEN
cc 	   OPEN (6,FILE=OFLNAM,ERR=900,IOSTAT=IVAR)
cc 	ELSE
cc 	   CALL SETWND
cc 	END IF
C                                                                       
      NPF = NDATA+FOCAST                                                
C                                                                       
C     PRESET CONSTANTS                                                  
C                                                                       
CC 1496 REWIND 5                                                          
cc      IY = 500                                                          
cc      IA = 500                                                          
cc      IRSLT = 500                                                       
cc      MDC=40000                                                         
cc      ALPHA = 0.01D0                                                    
cc      BETA = 0.01D0                                                     
cc      GAMMA = 0.1D0                                                     
C                                                                       
C     PRESET PARAMETERS                                                 
C                                                                       
cc      RIGID = 1.D0                                                      
cc      ZERSUM = 1.D0                                                     
cc      ORDER = 2                                                         
cc      SORDER = 1                                                        
cc      PERIOD = 12                                                       
cc      SPAN = 4                                                          
cc      SHIFT = 1                                                         
cc      FOCAST = 0                                                        
      RLIM = 0.0                                                        
cc      MT = 5                                                            
cc      DD = 1.D0                                                         
cc      LOGT=0                                                            
cc      PUNCH = 0                                                         
cc      WTRD = 1.D0                                                       
cc      YEAR = 0                                                          
cc      MONTH = 1                                                         
cc      IOUTD=0                                                           
cc      SPEC = 1                                                          
cc      NDAY = 1                                                          
cc      DELTA=7.D0                                                        
      PERIOD=IPARA(1)
      SPAN=IPARA(2)
      SHIFT=IPARA(3)
      ORDER=IPARA(4)
      SORDER=IPARA(5)
      LOGT=IPARA(6)
      YEAR=IPARA(7)
      MONTH=IPARA(8)
      NDAY=IPARA(9)
      SPEC=IPARA(10)
      IOUTD=IPARA(11)
      IDC=IPARA(12)
      RIGID=PARA(1)
      WTRD=PARA(2)
      DD= PARA(3)
      ZERSUM=PARA(4)
      DELTA=PARA(5)
      ALPHA=PARA(6)
      BETA=PARA(7)
      GAMMA=PARA(8)
cc      ARFT(1) = 0.0D0                                                   
cc      ARFT(2) = 0.D0                                                    
cc      ARFT(3) = 0.D0                                                    
cc      ARFS(1) = 0.0D0                                                   
cc      ARFS(2) = 0.0D0                                                   
cc      ARFS(3) = 0.D0                                                    
cc      ARFN(1) = 0.0D0                                                   
cc      ARFN(2) = 0.0D0                                                   
cc      ARFN(3) = 0.D0                                                    
cc      WRITE(6,1)                                                        
cc    1 FORMAT(1H ,'BAYSEA',/                                             
cc     * 1H ,'BAYESIAN TREND AND SEASONAL ESTIMATION OF TIME SERIES'      
cc     * ,'  VERSION 3/1/85',/                                            
cc     * ,'   WITH TRADING-DAY AND LEAP-YEAR ADJUSTMENT AND',/            
cc     * ,'        OUTLIER DETECTION OPTION')                             
C                                                                       
C     PARAMETER MODIFICATION AND                                        
C     DATA INPUT                                                        
C                                                                       
cc      READ(5,PARAM)
      IF(   SORDER .GT. SPAN ) SORDER = SPAN                            
cc      NDAY=1                                                            
cc      WRITE(6,610) LOGT,MT,RLIM,PERIOD,SPAN,SHIFT,FOCAST,ORDER,         
cc     *             SORDER,RIGID,YEAR,MONTH,WTRD,IOUTD,SPEC,PUNCH,       
cc     *             ZERSUM,DELTA
cc      CALL ARCOEF('ARFT',ARFT,IART)                                     
cc      CALL ARCOEF('ARFS',ARFS,IARS)                                     
cc      CALL ARCOEF('ARFN',ARFN,IARN)                       
      IS = PERIOD*SORDER                                                
      AP=PERIOD                                                         
cc      IPRD=2                                                            
cc      IF(PERIOD.EQ.1)IPRD=1                                             
      LFTRN = ORDER + IART                                              
      LFSEA = (SORDER + IARS)*PERIOD + IARN                             
cc      IDC=LFTRN*IPRD+1                                                  
cc      IDCX = LFSEA * 2 + 1                                              
cc      IF(PERIOD .GT. 1 .AND. IDC .LT. IDCX)  IDC=IDCX                   
cc      IF(PERIOD .GT. 1 .AND. IDC .LT. PERIOD*2-1) IDC=PERIOD*2-1
      NH= LFSEA + 1                                                     
      N2=1                                                              
      IF(YEAR .NE. 0) N2=8                                              
C  ************                                                         
      ZER=ZERSUM/DSQRT(AP)*RIGID                                        
      SMTH = 1.D0/RIGID                                                 
      SMTH2=1.D0                                                        
C  ************                                                         
CC      IF( MT .GT. 7 ) REWIND MT
cc      OPEN (MT,FILE=IFLNAM,ERR=910,IOSTAT=IVAR,STATUS='OLD')
cc      IF (PUNCH.EQ.1) THEN
cc         IF ((NFL.EQ.3) .OR. (NFL.EQ.4)) THEN
cc            OPEN (7,FILE=MFLNAM,ERR=920,IOSTAT=IVAR)
cc         ELSE
cc            OPEN (7,FILE='baysea.out',ERR=930,IOSTAT=IVAR)
cc         END IF
cc      END IF
C
cc      CALL DATAR( TITLE,Y,IY,NDATA,MT,LOGT,RLIM )
cc      CLOSE(MT)
      IF(IOUTD .EQ. 0) GO TO 1212                                       
      ROUT = 1.D60                                                      
      RLIM = 1.D50                                                      
      IF(LOGT .EQ. 0) GO TO 1212                                        
      ROUT = DLOG(ROUT)                                                 
      RLIM = DLOG(RLIM)                                                 
c-----
      DO 20 I=1,NDATA
   20 Y(I) = DLOG(Y(I))
c-----
 1212 CONTINUE                                                          
C     -------------------                                               
C                                                                       
C                                                                       
C     WORKING AREA  CHECK                                               
C                                                                       
cc      NPAR=(SPAN*2-1)*PERIOD                                            
cc      IF(NPAR.GT.NDATA+FOCAST)NPAR=NDATA+FOCAST                         
cc      IF( NDATA+FOCAST .LE. IRSLT )   GO TO 100                         
cc      WRITE( 6,603 )                                                    
cc  603 FORMAT( 1H , 'IRSLT IS TOO SMALL OR NDATA+FOCAST IS TOO LARGE' )  
cc      STOP                                                              
cc  100 IF( NPAR.LE.IA )   GO TO 200                                      
cc      WRITE( 6,605 )                                                    
cc  604 FORMAT( 1H ,'MDC IS TOO SMALL')                                   
cc      STOP                                                              
cc  200 CONTINUE                                                          
cc      NPAR=NPAR*IPRD                                                    
cc      IF( NPAR .LE. MDC )   GO TO 300                                   
cc      WRITE( 6,604 )                                                    
cc  605 FORMAT( 1H ,'IA IS TOO SMALL' )                                   
cc      STOP                                                              
cc  300 CONTINUE                                                          
C                                                                       
C                                                                       
C     INITIALIZATION                                                    
      NF=ORDER                                                          
      IF(PERIOD.GT.1.AND.NF.LT.IS)NF=IS                                 
cc      WRITE(6,606)                                                      
cc  606 FORMAT(1H ,'INITIALIZATION'    )                                  
C     N=LENGTH OF A SPAN                                                
      OVLAP = SPAN - 1                                                  
      LIMIT = NDATA - OVLAP*PERIOD                                      
      SY = 0.D0                                                         
c-----
      IQ = 0
c-----
      DO 2468 I=1,NDATA                                                 
      YTEM = Y(I)                                                       
      IF(RLIM .LE. 0.D0) GO TO 4681                                     
      IF(YTEM .GE. RLIM) GO TO 2468                                     
 4681 IQ = IQ + 1                   
      SY = SY + YTEM                                                    
      IF(IQ .GE. PERIOD) GO TO 4680                                     
 2468 CONTINUE                                                          
 4680 CONTINUE                                                          
      YTEM = SY / AP              
      DO 8  I=1,LFTRN                                                   
    8 FTRN(I) = YTEM                                                    
      IF(LFSEA .EQ. 0) GO TO 998                                        
      DO 9  I=1,LFSEA                                                   
    9 FSEA(I) = 0.D0                                                    
  998 CONTINUE                                                          
      N = (SPAN*2-1)*PERIOD                                             
      AVABIC = 0.D0                                                     
      COUNT = 0.D0                                                      
      IEND = 0                                                          
C                                                                       
C     ***************                                                   
C     **           **                                                   
C     ** MAIN LOOP **                                                   
C     **           **                                                   
C     ***************                                                   
C     ICNT0: ITERATION CONTROL FOR THE 0-TH SPAN                        
C                                                                       
      DO 1000  ICNT1=1,1000                                             
C                                                                       
C     MANIPULATION OF THE ICNT-TH SPAN                                  
C                                                                       
      ICNT = ICNT1-1                                                    
C                                                                       
C     DATA END DETECTION                                                
C                                                                       
C     --------------------                                              
      HEAD = 1 + (ICNT1-2)*SHIFT*PERIOD+SPAN*PERIOD                     
      IF(ICNT1 .EQ. 1) HEAD=1                                           
      IF( ICNT . LE. 0 )   GO TO 2345                                   
C     --------------------                                              
C     HEAD: INITIAL POINT OF THE NEW SPAN                               
      IF(HEAD.LE.LIMIT)   GO TO 2345                                    
      GO TO 1234                                                        
C     NEW SPAN                                                          
 2345 CONTINUE                                                          
C     TAIL: END POINT OF THE NEW SPAN                                   
      TAIL = HEAD+N-1                                                   
      IF( TAIL .GT. NDATA )   N = NDATA-HEAD+1                          
      TAIL = HEAD + N-1                                                 
C     --------------------                                              
      IF( TAIL .EQ. NDATA )   IEND = 1                                  
C     IEND=1: LAST SPAN IN NORMAL ITERATION                             
C      --------------------                                             
C                                                                       
C                                                                       
cc      CALL COPY(YO     ,N,1,N,1,1,Y     ,N,1,NDATA,HEAD,1)              
      CALL BCOPY(YO     ,N,1,N,1,1,Y     ,N,1,NDATA,HEAD,1)
C                                                                       
C     INITIALIZATION FOR THE INNER LOOP SUBSEA                          
C                                                                       
 6789 CONTINUE                                                          
      IF(YEAR .NE. 0) CALL CALEND(WEEK,YEAR,MONTH+HEAD-1,N+FOCAST)      
cc      WRITE(6,2) ICNT1,HEAD,TAIL                                        
    2 FORMAT(1H ,'(',I5,' )TH SPAN  HEAD =',I6,'   TAIL =',I6 )         
cc      WRITE( 6,620 )   (FTRN(I),I=1,LFTRN)                              
  620 FORMAT(1H ,'* INITIAL TREND  *'/,(1X,12D11.3))                    
cc      IF(LFSEA .NE. 0) WRITE( 6,621 )   (FSEA(I),I=1,LFSEA)             
  621 FORMAT(1H ,'* INITIAL SEASONAL  *'/,(1X,12D11.3))                 
      NEXT = HEAD + N                                                   
      IF( IEND .EQ. 0 ) NEXT = NEXT - OVLAP*PERIOD                      
      LINKT = NEXT - LFTRN                                              
      LINKS = NEXT - LFSEA                                              
C     --------------------                                              
C     DO SEARCH CONTROL                                                 
      ITRN = ICNT1                                                      
C     --------------------                                              
      IOUT=0                                                            
cc      CALL COPY(YS,N,1,N,1,1,YO, N,1,N,1,1)                             
      CALL BCOPY(YS,N,1,N,1,1,YO, N,1,N,1,1)
 9700 CONTINUE                                                          
C     SEASONAL DECOMPOSITION OF ICN-TH LOCAL SPAN                       
      CALL SUBSEA(ABIC,SEAS0,TREND0,EST0,ADJ0,IRREG0,TDCMP0,            
cc     *  FSEA,FTRN,YS,N,FOCAST,RLIM,WEEK,DC,IDC,H,NH,F,H2,N2,ITRN,       
cc     *  IARS,ARFS,IART,ARFT,IARN,ARFN,PSDT0,PSDS0)                      
     *  FSEA,FTRN,YS,N,FOCAST,RLIM,WEEK,IDC,NH,F,N2,ITRN,       
     *  IARS,ARFS,IART,ARFT,IARN,ARFN,PSDT0,PSDS0,NPF,
     *  PERIOD,ORDER,SORDER,YEAR,NDAY,LOGT,
     *  ALPHA,BETA,GAMMA,ZER,SMTH,SMTH2,DD,WTRD,DELTA)
      IF(IOUTD .EQ. 0) GO TO 9600
      IF(IOUT .GE. 2) GO TO 9600
      CALL SBTRCT(IRREG0,N,YO,N,EST0,N)                                 
cc      CALL OUTLIR(IRREG0,N,10,2,1,YS1,RLIM) 
      CALL OUTLIR(IRREG0,N,10,2,1,YS1,RLIM,IOUTD,ROUT)
      IOUT=IOUT+1                                                       
      CALL ADD(YS,N,EST0,N,YS1,N)                                       
cc      IF(IOUT .EQ. 2) CALL COPY(YS1,N,1,N,1,1, YS,N,1,N,1,1)            
      IF(IOUT .EQ. 2) CALL BCOPY(YS1,N,1,N,1,1, YS,N,1,N,1,1)
      IF(IOUT .NE. 1) GO TO 9700                                        
      NTEM = N+1                                                        
      DO 9702 I=1,ORDER                                                 
      NTEM=NTEM-1                                                       
 9702 YS(NTEM)=ROUT                                                     
      GO TO 9700                                                        
 9600 CONTINUE                                                          
C                                                                       
C                                                                       
C     RECORDING THE BEST RESULT                                         
C                                                                       
      L=NEXT - HEAD                                                     
      LF=L+FOCAST                                                       
cc      CALL COPY(PSDT,LF,1,LF,HEAD,1,PSDT0,LF,1,LF,1,1)                  
cc      CALL COPY(PSDS,LF,1,LF,HEAD,1,PSDS0,LF,1,LF,1,1)                  
cc      CALL COPY(SEASON,LF,1,LF,HEAD,1,SEAS0,LF,1,LF,1,1)                
cc      CALL COPY(TREND,LF,1,LF,HEAD,1,TREND0,LF,1,LF,1,1)                
cc      CALL COPY(EST,LF,1,LF,HEAD,1,EST0,LF,1,LF,1,1)                    
cc      CALL COPY(ADJUST,L,1,L,HEAD,1,ADJ0,L,1,L,1,1)                     
cc      CALL COPY(IRREG,L,1,L,HEAD,1,IRREG0,L,1,L,1,1)                    
cc      CALL COPY(TDCMP,LF,1,LF,HEAD,1,TDCMP0,LF,1,LF,1,1)                
cc      CALL COPY(CDATA,L,1,L,HEAD,1,YS1,L,1,L,1,1)                       
      CALL BCOPY(PSDT,LF,1,LF,HEAD,1,PSDT0,LF,1,LF,1,1)
      CALL BCOPY(PSDS,LF,1,LF,HEAD,1,PSDS0,LF,1,LF,1,1)
      CALL BCOPY(SEASON,LF,1,LF,HEAD,1,SEAS0,LF,1,LF,1,1)
      CALL BCOPY(TREND,LF,1,LF,HEAD,1,TREND0,LF,1,LF,1,1)
      CALL BCOPY(EST,LF,1,LF,HEAD,1,EST0,LF,1,LF,1,1)
      CALL BCOPY(ADJUST,L,1,L,HEAD,1,ADJ0,L,1,L,1,1)
      CALL BCOPY(IRREG,L,1,L,HEAD,1,IRREG0,L,1,L,1,1)
      CALL BCOPY(TDCMP,LF,1,LF,HEAD,1,TDCMP0,LF,1,LF,1,1)
      CALL BCOPY(CDATA,L,1,L,HEAD,1,YS1,L,1,L,1,1)
      AN = N                                                            
      AVABIC = AVABIC + ABIC                                            
      COUNT = COUNT + AN                                                
C                                                                       
C                                                                       
C     INITIAL VALUES FOR THE NEXT SPAN                                  
C                                                                       
      IF(IEND .EQ. 1) GO TO 1234                                        
cc      CALL  COPY( FTRN,LFTRN,1,LFTRN,1,1,TREND,LFTRN,1,IOUT,            
      CALL  BCOPY( FTRN,LFTRN,1,LFTRN,1,1,TREND,LFTRN,1,IOUT,
     * LINKT,1)                                                         
      ISTEM = LFSEA                                                     
      IF(LINKS .GE. 1) GO TO 1111                                       
      LINKS=1-LINKS                                                     
      ISTEM=ISTEM-LINKS                                                 
      I1=LFSEA+1                                                        
      DO 2222 I=1,LINKS                                                 
      I1=I1-1                                                           
      I2=I1-ISTEM                                                       
 2222 FSEA(I1)=FSEA(I2)                                                 
      LINKS=1                                                           
 1111 CONTINUE                                                          
cc      CALL COPY(FSEA,ISTEM,1,ISTEM,1,1,SEASON,ISTEM,1,IOUT,LINKS,1)     
      CALL BCOPY(FSEA,ISTEM,1,ISTEM,1,1,SEASON,ISTEM,1,IOUT,LINKS,1)
      IF(ICNT1 .GT. 1) GO TO 1000                                       
      ALPHA = 1.D0                                                      
      BETA = 1.D0                                                       
      GAMMA = 1.D0                                                      
      N = SPAN*PERIOD                                                   
      IF(N .GT. NDATA) N=NDATA                                          
      OVLAP = SPAN-SHIFT                                                
      LIMIT = NDATA-OVLAP*PERIOD                                        
 1000 CONTINUE                                                          
C     ************************                                          
C     *                      *                                          
C     * END OF THE MAIN LOOP *                                          
C     *                      *                                          
C     ************************                                          
C                                                                       
 1234 CONTINUE                                                          
C                                                                       
C     NUMERICAL OUTPUTS                                                 
C                                                                       
      DO 4444 I=1,NDATA                                                 
      CDATA(I)=Y(I) - CDATA(I)                                          
      DMOI(I) = Y(I)                                                    
      IF(RLIM .LE. 0.D0) GO TO 4444                                     
      IF(Y(I) .LT. RLIM.AND.IOUTD.EQ.0) GO TO 4444                      
      IF(IOUTD .NE. 0 .AND. Y(I) .GT. RLIM) GO TO 4442                  
      IF(IOUTD.NE.0.AND.-CDATA(I).LT.RLIM)GO TO 4443                    
 4442 CONTINUE                                                          
      ADJUST(I)=TREND(I)                                                
      DMOI(I)=EST(I)                                                    
      IRREG(I)=0.D0                                                     
 4443 CDATA(I) = Y(I) - DMOI(I)                                         
 4444 CONTINUE                                                          
cc      NPF = NDATA+FOCAST                                                
      IF(LOGT .EQ. 0) GO TO 1250                                        
      DO 1240 I=1,NPF                                                   
      TREND(I) = DEXP(TREND(I))                                         
      SEASON(I) = DEXP(SEASON(I))                                       
      EST(I) = DEXP(EST(I))                                             
      TDCMP(I) = DEXP(TDCMP(I))                                         
      IF(I .GT. NDATA) GO TO 1240                                       
      IRREG(I) = DEXP(IRREG(I))                                         
      Y(I) = DEXP(Y(I))                                                 
      ADJUST(I) = DEXP(ADJUST(I))                                       
      CDATA(I) = DEXP(CDATA(I))                                         
      DMOI(I) = DEXP(DMOI(I))                                           
 1240 CONTINUE                                                          
      IF(RLIM .GT. 0.D0) RLIM = DEXP(RLIM)                              
 1250 CONTINUE                                                          
cc      DO 1324 I=1,NDATA                                                 
cc 1324 DADJ(I)=ADJUST(I)                                                 
cc      CALL DFR1(1,1,NDATA,DADJ,NDATA1)                                  
cc      IF(PUNCH .EQ. 0) GO TO 5555                                       
cc      WRITE(7,714) (TITLE(I),I=1,20)                                    
cc      WRITE(7,720) NDATA                                                
cc      WRITE(7,700)                                                      
cc      WRITE(7,702) (Y(I),I=1,NDATA)                                     
cc      WRITE(7,716)                                                      
cc      WRITE(7,720) NPF                                                  
cc      WRITE(7,700)                                                      
cc      WRITE(7,702) (TREND(I),I=1,NPF)                                   
cc      WRITE(7,717)                                                      
cc      WRITE(7,720) NPF                                                  
cc      WRITE(7,700)                                                      
cc      WRITE(7,702) (SEASON(I),I=1,NPF)                                  
cc      WRITE(7,718)                                                      
cc      WRITE(7,720) NDATA                                                
cc      WRITE(7,700)                                                      
cc      WRITE(7,702) (IRREG(I),I=1,NDATA)                                 
cc      WRITE(7,719)                                                      
cc      WRITE(7,720) NDATA                                                
cc      WRITE(7,700)                                                      
cc      WRITE(7,702) (ADJUST(I),I=1,NDATA)                                
cc      WRITE(7,703)                                                      
cc      WRITE(7,720) NDATA1                                               
cc      WRITE(7,700)                                                      
cc      WRITE(7,702) (DADJ(I),I=1,NDATA1)                                 
cc      WRITE(7,715)                                                      
cc      WRITE(7,720) NPF                                                  
cc      WRITE(7,700)                                                      
cc      WRITE(7,702) (EST(I),I=1,NPF)                                     
cc      IF(IOUTD .EQ. 0) GO TO 8500                                       
cc      WRITE(7,721)                                                      
cc      WRITE(7,720)NDATA                                                 
cc      WRITE(7,700)                                                      
cc      WRITE(7,702) (CDATA(I),I=1,NDATA)                                 
cc 8500 CONTINUE                                                          
cc      IF(YEAR .EQ. 0) GO TO 5552                                        
cc      WRITE(7,730)                                                      
cc      WRITE(7,720) NPF                                                  
cc      WRITE(7,700)                                                      
cc      WRITE(7,702) (TDCMP(I),I=1,NPF)                                   
cc 5552 IF(RLIM .LE. 0.D0) GO TO 5555                                     
cc      WRITE(7,777)                                                      
cc      WRITE(7,720) NDATA                                                
cc      WRITE(7,700)                                                      
cc      WRITE(7,702) (DMOI(I),I=1,NDATA)                                  
cc 5555 CONTINUE                                                          
      AVABIC = AVABIC/COUNT                                             
      AVABIC = AVABIC*NDATA                                             
cc      WRITE(6,600) AVABIC                                               
cc      IF(IOUTD .EQ. 0) GO TO 5553                                       
cc      WRITE(6,641)                                                      
cc      WRITE(6,602) (CDATA(I),I=1,NDATA)                                 
cc 5553 CONTINUE                                                          
cc      WRITE( 6,616 )                                                    
cc      WRITE( 6,602 )   (TREND(I),I=1,NPF)                               
cc      IF(FOCAST .NE. 0) WRITE( 6,630 ) FOCAST                           
cc      WRITE( 6,617 )                                                    
cc      WRITE( 6,602 )   (SEASON(I),I=1,NPF)                              
cc      IF(FOCAST .NE. 0) WRITE( 6,630 ) FOCAST                           
cc      IF(YEAR .EQ. 0) GO TO 5432                                        
cc      WRITE(6,640)                                                      
cc      WRITE(6,602) (TDCMP(I),I=1,NPF)                                   
cc      IF(FOCAST .NE. 0) WRITE(6,630) FOCAST                             
cc 5432 CONTINUE                                                          
cc      WRITE( 6,618 )                                                    
cc      WRITE( 6,602 )   (IRREG(I),I=1,NDATA)                             
cc      WRITE( 6,619 )                                                    
cc      WRITE( 6,602 )   (ADJUST(I),I=1,NDATA)                            
cc      WRITE( 6,615 )                                                    
cc      WRITE( 6,602 )   (EST(I),I=1,NPF)                                 
cc      IF(FOCAST .NE. 0) WRITE( 6,630 ) FOCAST                           
C                                                                       
C     GRAPHICAL OUTPUTS                                                 
C                                                                       
cc      ISD = 1                                                           
cc      IF(LOGT .NE. 0) ISD = 2                                           
cc      WRITE( 6,611 )                                                    
cc      WRITE( 6,614 )   (TITLE(I),I=1,20)                                
cc      CALL GRAPH(0,PSDT,Y,NDATA,RLIM,1)                                 
cc      IF(RLIM .LE. 0.D0) GO TO 8301                                     
cc      WRITE(6,611)                                                      
cc      WRITE(6,776)                                                      
cc      CALL GRAPH(0,PSDT,DMOI,NDATA,RLIM,0)                              
cc 8301 CONTINUE                                                          
cc      WRITE( 6,611 )                                                    
cc      WRITE( 6,616 )                                                    
cc      CALL GRAPH(ISD,PSDT,TREND,NDATA+FOCAST,RLIM,0)                    
cc      IF(FOCAST .NE. 0) WRITE( 6,630 ) FOCAST                           
cc      WRITE( 6,611 )                                                    
cc      WRITE( 6,619 )                                                    
cc      CALL GRAPH(0,PSDT,ADJUST,NDATA,RLIM,0)                            
cc      WRITE( 6,611 )                                                    
cc      WRITE( 6,615 )                                                    
cc      CALL GRAPH(0,PSDT,EST,NDATA+FOCAST,RLIM,0)                        
cc      IF(FOCAST .NE. 0) WRITE( 6,630 ) FOCAST                           
cc      ITEM = -1                                                         
cc      IF(LOGT .NE. 0) ITEM=1                                            
cc      WRITE( 6,611 )                                                    
cc      WRITE( 6,617 )                                                    
cc      CALL GRAPH(ISD,PSDS,SEASON,NDATA+FOCAST,RLIM,ITEM)                
cc      IF(FOCAST .NE. 0) WRITE( 6,630 ) FOCAST                           
cc      WRITE( 6,611 )                                                    
cc      WRITE( 6,618 )                                                    
cc      CALL GRAPH(0,PSDT,IRREG,NDATA,RLIM,0)                             
cc      IF(IOUTD .EQ. 0) GO TO 8300                                       
cc      WRITE(6,611)                                                      
cc      WRITE(6,641)                                                      
cc      CALL GRAPH(0,PSDT,CDATA,NDATA,RLIM,0)                             
cc 8300 CONTINUE                                                          
cc      IF(LOGT .EQ. 1) GO TO 3332                                        
cc      WRITE( 6,611 )                                                    
cc      WRITE( 6,622 )                                                    
cc      CALL  GRAPH(0,PSDT, IRREG,NDATA,RLIM,2 )                          
cc 3332 IF(YEAR .EQ. 0) GO TO 3333                                        
cc      WRITE(6,611)                                                      
cc      WRITE(6,640)                                                      
cc      CALL GRAPH(0,PSDT,TDCMP,NDATA+FOCAST,RLIM,0)                      
cc      IF(FOCAST .NE. 0) WRITE(6,630) FOCAST                             
cc 3333 CONTINUE                                                          
cc      IF(SPEC .EQ. 0) GO TO 9999                                        
cc      IF(PUNCH .NE. 0) WRITE(7,8888)                                    
cc 8888 FORMAT('IRREGULAR(SPECTRUM)')                                     
cc      WRITE(6,800)                                                      
cc      CALL SPGRH(IRREG,NDATA,1)                                         
cc      WRITE(6,810)                                                      
cc      IF(PUNCH .NE. 0) WRITE(7,8887)                                    
cc 8887 FORMAT('DADJ(SPECTRUM)')                                          
cc      CALL SPGRH(DADJ,NDATA1,1)                                         
cc  800 FORMAT(1H ,'SPECTRUM OF IRREGULAR')                               
cc  810 FORMAT(1H ,'SPECTRUM OF DIFFERENCED ADJUSTED SERIES')             
cc 1810 FORMAT(1H ,'PARCOR OF',I2,' TIME(S) DIFFERENCED TREND SERIES')    
cc 1811 FORMAT(1H ,'PARCOR OF',I2,' TIME(S) DIFFERENCED SEASONAL SERIES') 
cc 9999 CONTINUE                                                          
cc      WRITE(6,1810) ORDER                                               
cc      CALL DFR1(ORDER,1,NDATA,TREND,NDATA1)                             
cc      CALL SPGRH(TREND,NDATA1,0)                                        
cc      WRITE(6,1811) SORDER                                              
cc      CALL DFR1(SORDER,PERIOD,NDATA,SEASON,NDATA1)                      
cc      CALL SPGRH(SEASON,NDATA1,0)
cc      GO TO 990
C
cc  900 CONTINUE
cc      WRITE(6,690) IVAR,OFLNAM
cc      GO TO 999
cc  910 CONTINUE
cc      IF ((NFL.EQ.2) .OR. (NFL.EQ.4)) CLOSE(6)
cc      WRITE(6,691) IVAR,IFLNAM
cc      GO TO 999
cc  920 CONTINUE
cc      CLOSE(5)
cc      IF ((NFL.EQ.2) .OR. (NFL.EQ.4)) CLOSE(6)
cc      WRITE(6,692) IVAR,MFLNAM
cc      GO TO 999
cc  930 CONTINUE
cc      CLOSE(5)
cc      IF ((NFL.EQ.2) .OR. (NFL.EQ.4)) CLOSE(6)
cc      WRITE(6,693) IVAR
cc      GO TO 999
C
  690 FORMAT(1H ,' !!! Output_Data_File OPEN ERROR ',I8/1H ,100A)
  691 FORMAT(1H ,' !!! Input_Data_File OPEN ERROR ',I8/1H ,100A)
  692 FORMAT(1H ,' !!! Intermediate_Data_File OPEN ERROR ',I8/1H ,100A)
  693 FORMAT(1H ,' !!! baysea.out  OPEN ERROR ',I8)
C
cc  990 CONTINUE
cc      IF ((NFL.EQ.2) .OR. (NFL.EQ.4)) CLOSE(6)
cc      IF (PUNCH.EQ.1) CLOSE(7)
cc  999 CONTINUE
cc      STOP                                                              
      RETURN
    3 FORMAT( 1H ,'*** ABIC(',D20.10,'  ) =  ',D20.10 )                 
  600 FORMAT(1H ,'AVABIC =',F10.2)                                      
  602 FORMAT(1H ,12D11.3)                                               
  614 FORMAT( 1H ,'ORIGINAL DATA',/,1H , 20A4 )                         
  615 FORMAT(1H ,'SMOOTHED=TREND+SEASONAL+TRADING.DAY.COMP')            
  630 FORMAT(1H ,'***  LAST ',I3,' VALUES ARE FOCASTED  ***')           
  616 FORMAT( 1H ,'TREND')                                              
  617 FORMAT( 1H ,'SEASONAL')                                           
  618 FORMAT(1H ,'IRREGULAR=ORIGINAL DATA-TREND-SEASONAL-TRADING.',     
     *'DAY.COMP')                                                       
  619 FORMAT(1H ,'ADJUSTED=ORIGINAL DATA-SEASONAL-TRADING.DAY.COMP-OCF')
  622 FORMAT( 1H ,'IRREGULAR ( SCALED BY THE STANDARD DEVIATION )')     
  610 FORMAT(1H ,                                                       
     * 'LOGT  =',I10,/,                                                 
     *' MT    =',I10,/,                                                 
     *' RLIM  =',D10.3,/,                                               
     *' PERIOD=',I10,/,                                                 
     *' SPAN  =',I10,/,                                                 
     *' SHIFT =',I10,/,                                                 
     *' FOCAST=',I10,/,                                                 
     *' ORDER =',I10,/,                                                 
     *' SORDER=',I10,/,                                                 
     *' RIGID =',D10.3,/,                                               
     *' YEAR  =',I10,/,                                                 
     *' MONTH =',I10,/,                                                 
     *' WTRD  =',D10.3,/,                                               
     *' IOUTD =',I10,/,                                                 
     *' SPEC  =',I10,/,                                                 
     *' PUNCH =',I10,/,                                                 
     *' ZERSUM=',D10.3,/,                                               
     *' DELTA =',D10.3)                                                 
  611 FORMAT( 1H  )                                                     
  612 FORMAT(1H ,'D3    =',10F10.6)                                     
  640 FORMAT(1H ,'TRADING-DAY COMPONENT')                               
  641 FORMAT(1H ,'OUTLIER CORRECTION FACTOR')                           
  702 FORMAT(6D12.5)                                                    
  700 FORMAT('(6D12.5)')                                                
  714 FORMAT(20A4)                                                      
  715 FORMAT('SMOOTHED=TREND+SEASONAL+TRADING.DAY.COMP')                
  716 FORMAT('TREND')                                                   
  717 FORMAT('SEASONAL')                                                
  718 FORMAT('IRREGULAR')                                               
  719 FORMAT('ADJUSTED=DATA-SEASONAL-TRADING.DAY.COMP-OCF')             
  720 FORMAT(I5)                                                        
  721 FORMAT('OUTLIER CORRECTION FACTOR')                               
  730 FORMAT('TRADING DAY COMPONENT')                                   
  703 FORMAT('DADJ')                                                    
  777 FORMAT('MISSING OBSERVATION INTERPOLATED DATA')                   
  776 FORMAT(1H ,'MISSING OBSERVATION INTERPOLATED DATA')               
      END                                                               
      SUBROUTINE  ADD(X,MX,Y,MY,Z,MZ)                                   
C     THIS SUBROUTINE COMPUTES                                          
C          X = Y + Z.                                                   
C     INPUTS:                                                           
C       X:     MX-VECTOR                                                
C       Y:     MY-VECTOR                                                
C       Z:     MZ-VECTOR                                                
C                                                                       
      IMPLICIT REAL*8 ( A-H,O-Z )                                       
      DIMENSION X(1), Y(1), Z(1)                                        
      DO 100 I=1,MX                                                     
      TEM = 0.D0                                                        
      IF( I .LE. MY )  TEM = Y(I)                                       
      IF( I .LE. MZ )  TEM = TEM + Z(I)                                 
  100 X(I) = TEM                                                        
      RETURN                                                            
      END                                                               
      SUBROUTINE CLEAR(X,M,N,MJ,I0,J0)                                  
C                                                                       
C     THIS SUBROUTINE CLEARS MATRIX X.                                  
C     INPUTS:                                                           
C       X:     M*N MATRIX                                               
C       MJ:    ABSOLUTE DIMENSION OF X                                  
C       I0:    ABSOLUTE POSITION OF THE FIRST ROW OF X                  
C       J0:    ABSOLUTE POSITION OF THE FIRST COLUMN OF X               
C                                                                       
      IMPLICIT REAL*8 (A-H,O-Z )                                        
      DIMENSION X(MJ,1)                                                 
      I0M1 = I0 - 1                                                     
      J0M1 = J0 - 1                                                     
      DO 10 J=1,N                                                       
      DO 10 I=1,M                                                       
   10 X(I0M1+I,J0M1+J) = 0.D0                                           
      RETURN                                                            
      END                                                               
cc      SUBROUTINE COPY(X,MX,NX,MMX,IX,JX,Y,MY,NY,MMY,IY,JY)              
      SUBROUTINE BCOPY(X,MX,NX,MMX,IX,JX,Y,MY,NY,MMY,IY,JY)
C     THIS SUBROUTINE COPIES Y INTO X.                                  
C     INPUTS:                                                           
C       X:     MX*NX MATRIX                                             
C       MMX:   ABSOLUTE DIMENSION OF X                                  
C       IX:    ABSOLUTE POSITION OF THE FIRST ROW OF X                  
C       JX:    ABSOLUTE POSITION OF THE FIRST COLUMN OF X               
C       Y:     MY*NY MATRIX                                             
C       MMY:   ABSOLUTE DIMENSION OF Y                                  
C       IY:    ABSOLUTE POSITION OF THE FIRST ROW OF Y                  
C       JY:    ABSOLUTE POSITION OF THE FIRST COLUMN OF Y               
C                                                                       
      IMPLICIT  REAL*8 ( A-H,O-Z )                                      
      DIMENSION X(MMX,1), Y(MMY,1)                                      
      IXM1 = IX-1                                                       
      JXM1 = JX - 1                                                     
      IYM1 = IY - 1                                                     
      JYM1 = JY - 1                                                   
      DO 100 J=1,NX                                                     
      DO 50 I=1,MX                                                      
      TEM = 0.D0                                                        
      IF( I .GT. MY ) GO TO 50                                          
      IF( J .GT. NY ) GO TO 50                                          
      TEM = Y(IYM1+I,JYM1+J)                                            
   50 X(IXM1+I,JXM1+J) = TEM                                            
  100 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE DECODE(SEAS0,TREND0,EST0,ADJ0,IRREG0,TDC0,W,           
cc     * A,Y,NN,NF,WEEK,ERR,PSDS,PSDT,SQE)                                
     * A,Y,NN,NF,WEEK,ERR,PSDS,PSDT,SQE, IP,YEAR,NDAY)
C   THIS SUBROUTINE COMPUTES                                            
C       TREND0                                                          
C       SEAS0                                                           
C       EST0=TREND0 + SEAS                                              
C       ADJ0=Y - SEAS0                                                  
C       IRREG0=Y - EST0                                                 
C                                                                       
      IMPLICIT REAL*8 (A-H,O-Z )                                        
      INTEGER*4 YEAR                                                    
      REAL*8 IRREG0                                                     
      DIMENSION A(1),Y(1),SEAS0(1),TREND0(1),EST0(1),ADJ0(1),IRREG0(1)  
cc     * ,W(1),WEEK(7,1),TDC0(1), PSDT(500),PSDS(500),ERR(1000)           
     * ,W(1),WEEK(7,1),TDC0(1), PSDT(NN+NF),PSDS(NN+NF),ERR(2*(NN+NF))
cc      COMMON /IDATA/ IP,IDUMMY(2),YEAR,NDAY                             
C                                                                       
      N=NN+NF                                                           
      NR = 2                                                            
      IF( IP .EQ. 1 )   NR = 1                                          
      CALL  CLEAR( SEAS0,N,1,N,1,1 )                                    
      CALL  CLEAR( PSDS,N,1,N,1,1)                                      
      SD2 = DSQRT(SQE) * 2.D0                                           
      DO 10 I=1,N                                                       
      I1=NR*(I-1)+1                                                     
      I2=NR*I                                                           
      TREND0(I)=A(I1)                                                   
      PSDT(I)=DSQRT(ERR(I1))*SD2                                        
      IF(IP .LE. 1) GO TO 10                                            
      SEAS0(I)=A(I2)                                                    
      PSDS(I)=DSQRT(ERR(I2))*SD2                                        
   10 CONTINUE                                                          
      IF(YEAR .EQ. 0) GO TO 20                                          
      NTEM = N*2+1                                                      
      N7 = 6 + NDAY                                                     
cc      CALL COPY(W,N7,1,8,1,1,A,N7,1,N,NTEM,1)                           
      CALL BCOPY(W,N7,1,8,1,1,A,N7,1,N,NTEM,1)
      CALL PRDCT(TDC0,1,N,1,W,1,N7,1,WEEK,N7,N,N7)                      
   20 CONTINUE                                                          
      CALL ADD(EST0,N,TREND0,N,SEAS0,N)                                 
      IF(YEAR .NE. 0) CALL ADD(EST0,N,EST0,N,TDC0,N)                    
      CALL SBTRCT(ADJ0,N,Y,N,SEAS0,N)                                   
      IF(YEAR .NE. 0) CALL SBTRCT(ADJ0,N,ADJ0,N,TDC0,N)                 
      CALL SBTRCT(IRREG0,N,Y,N,EST0,N)                                  
      RETURN                                                            
      END                                                               
cc      SUBROUTINE  HUSHLD( X,N,K,MJ1,ICNT )                              
      SUBROUTINE  BHUSHLD( X,N,K,MJ1,ICNT )
C                                                                       
C                                                                       
C     THIS SUBROUTINE  TRANSFORMS MATRIX X INTO AN UPPER TRIANGULAR FORM
C     BY HOUSEHOLDER TRANSFORMATION.                                    
C                                                                       
C       INPUTS:                                                         
C          X:     ORIGINAL N*K DATA MATRIX                              
C          MJ1:   ABSOLUTE DIMENSION OF X                               
C          N:  NUMBER OF ROWS OF X, NOT GREATER THAN MJ1                
C          K:     NUMBER OF COLUMNS OF X. NOT GREATER THAN MJ1          
C          ICNT:  =0  WHEN X IS A FULL N*K MATRIX                       
C                 =L (L.NE.0) WHEN X IS COMPOSED OF TWO UPPER TRIANGULAR
C                             MATRICES. L SHOULD BE EQUAL TO THE NUMBER 
C                             OF NON ZERO ROWS OF THE SECOND MATRIX.    
C                             SECOND MATRIX MUST BE ROTATED AND STORED  
C                             AT THE BOTTOM LEFT PART OF THE FIRST      
C                             MATRIX.                                   
C       OUTPUT:                                                         
C          X:     IN UPPER TRIANGULAR FORM                              
C                                                                       
      IMPLICIT  REAL * 8  ( A-H , O-Z )                                 
cc      DIMENSION  X(MJ1,1) , D(1000)                                     
      DIMENSION  X(MJ1,1) , D(N)                                     
C                                                                       
      TOL=1.0D-38                                                       
      DIIO=0.0D00
C                                                                       
cc      IF( MJ1 .GT. 1000) STOP                                           
      MNK=K                                                             
      IF(N.LE.K) MNK=N-1                                                
      DO 100 II=1,MNK                                                   
         H = 0.0D00                                                     
      IIOTEM = II                                                       
      IITEM = II                                                        
      IF( ICNT .LE. 0 )   GO TO 5                                       
      H = X(II,II)*X(II,II)                                             
      IITEM = K+1-II                                                    
      IIOTEM = N+1-II                                                   
      IF( IIOTEM .LE. N-ICNT )   IIOTEM = N-ICNT+1                      
    5 CONTINUE                                                          
      DO 10  I=IIOTEM,N                                                 
      D(I) = X(I,IITEM)                                                 
      ABSLD=DABS(D(I))                                                  
      IF(ABSLD.LE.TOL) D(I)=0.0D-00                                     
   10       H = H + D(I)*D(I)                                           
         IF( H .GT. TOL )  GO TO 20                                     
         G = 0.0D00                                                     
         GO TO 100                                                      
   20    G = DSQRT( H )                                                 
      F=X(II,II)                                                        
         IF( F .GE. 0.0D00 )   G = -G                                   
      IF( ICNT .LE. 0 )   D(II) = F-G                                   
      IF( ICNT .GT. 0 )   DIIO = F-G                                    
         H = H - F*G                                                    
C                                                                       
C          FORM  (I - D*D'/H) * X, WHERE H = D'D/2                      
C                                                                       
         II1 = IITEM+1                                                  
         KTEM = K                                                       
         IF( ICNT .LE. 0 )   GO TO 25                                   
         II1 = 1                                                        
         KTEM = IITEM-1                                                 
   25    CONTINUE                                                       
      II10 = II1                                                        
      IF( ICNT .GT. 0 )   II10 = IIOTEM                                 
      DO 30 I=II10,N                                                    
   30 X(I,IITEM) = 0.D0                                                 
         IF( II .EQ. K )  GO TO 100                                     
         DO 60  J=II1,KTEM                                              
            S = 0.0D00                                                  
         JTEM = K+1-J                                                   
      IF(ICNT .GT. 0 ) S = DIIO*X(II,JTEM)                              
      DO 40  I=IIOTEM,N                                                 
   40 S = S + D(I)*X(I,J)                                               
            S = S/H                                                     
      IF(ICNT .GT. 0) X(II,JTEM) = X(II,JTEM) - DIIO*S                  
      DO 50  I=IIOTEM,N                                                 
   50 X(I,J) = X(I,J) - D(I)*S                                          
   60    CONTINUE                                                       
  100 X(II,II) = G                                                      
C                                                                       
      RETURN                                                            
C                                                                       
      E N D                                                             
      SUBROUTINE  PRDCT(X,MX,NX,MMX,Y,MY,NY,MMY,Z,MZ,NZ,MMZ)            
C     THIS SUBROUTINE COMPUTES                                          
C          X = Y * Z                                                    
C     INPUTS:                                                           
C       X:     MX*NX MATRIX                                             
C       MMX:   ABSOLUTE DIMENSION OF X                                  
C       Y:     MY*NY MATRIX                                             
C       MMY:   ABSOLUTE DIMENSION OF Y                                  
C       Z:     MZ*NZ MATRIX                                             
C       MMZ:   ABSOLUTE DIMENSION OF Z                                  
C                                                                       
      IMPLICIT REAL*8 ( A-H,O-Z )                                       
      DIMENSION X(MMX,1), Y(MMY,1), Z(MMZ,1)                            
      KK = NY                                                           
      IF( KK .GT. MZ ) KK = MZ                                          
      DO 100 J=1,NX                                                     
      DO 50 I=1,MX                                                      
      SUM = 0.D0                                                        
      IF(I .GT. MY) GO TO 50                                            
      IF( J .GT. NZ ) GO TO 50                                          
      DO 20 K=1,KK                                                      
   20 SUM = SUM + Y(I,K)*Z(K,J)                                         
   50 X(I,J) = SUM                                                      
  100 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE  SBTRCT(X,MX,Y,MY,Z,MZ)                                
C     THIS SUBROUTINE COMPUTES                                          
C          X = Y - Z                                                    
C     INPUTS:                                                           
C       X:     MX-VECTOR                                                
C       Y:     MY-VECTOR                                                
C       Z:     MZ-VECTOR                                                
C                                                                       
      IMPLICIT REAL*8 ( A-H,O-Z )                                       
      DIMENSION X(1), Y(1), Z(1)                                        
      DO 100 I=1,MX                                                     
      TEM = 0.D0                                                        
      IF( I .LE. MY )  TEM = Y(I)                                       
      IF( I .LE. MZ )  TEM = TEM - Z(I)                                 
  100 X(I) = TEM                                                        
      RETURN                                                            
      END                                                               
      SUBROUTINE SUBSEA(ABICM,SEASON,TREND,EST,ADJ,IRREG,TDC,           
cc     *   FSEA,FTRN,YS,N,NF,RLIM,WEEK,DC,IDC,H,NH,F,H2,N2,ITRN,          
cc     *  IARS,ARFS,IART,ARFT,IARN,ARFN,PSDT,PSDS)                        
     *    FSEA,FTRN,YS,N,NF,RLIM,WEEK,IDC,NH,F,N2,ITRN,IARS,ARFS,IART,
     *   ARFT,IARN,ARFN,PSDT,PSDS,NPF,PERIOD,IORD,ISOD,YEAR,NDAY,LOGT,
     *  ALPHA,BETA,GAMMA,ZER,SMTH,SMTH2,DD,WTRD,DELTA)
C     SEASONAL DECOMPOSITION PROCEDURE                                  
C     FOR THE DEFINITIONS OF THE VARIABLES APPEARING IN THE ARGUMENTS,  
C     SEE THE COMMENTS IN THE MAIN ROUTINE                              
      IMPLICIT  REAL*8  ( A-H,O-Z )                                    
cc      REAL*8  IRREG0, IRREG                                             
      REAL*8  IRREG                                             
      INTEGER*4 PERIOD,YEAR                                             
cc      DIMENSION   FSEA(1), FTRN(1), DC(IDC,1), YS(1), H2(N2,1)
cc     *   ,ERR(1000),PSDT(500),PSDS(500)                                 
cc     *  ,SEASON(1),TREND(1),EST(1),ADJ(1),IRREG(1),A(1000)
cc     * ,WEEK(7,1),TDC(1),WEEK0(8),WEEK1(8),H(NH,1),F(1)                 
cc      DIMENSION  DTRN(500),DSEAS(500), ARFS(1), ARFT(1), ARFN(1)        
      DIMENSION  SEASON(1),TREND(1),EST(1),ADJ(1),IRREG(1),TDC(1),
     *   FSEA(1),FTRN(1),YS(1),WEEK(7,1),F(1),ARFS(1),ARFT(1),ARFN(1),
     *  PSDT(NPF),PSDS(NPF)
      DIMENSION  DC(IDC,2*NPF+N2),H(NH,NPF),H2(N2,2*NPF+N2),WEEK0(8),
     *   WEEK1(8),ERR(2*(N+NF)+NDAY+7),A(2*(N+NF)+NDAY+7),DTRN(NPF),
     *  DSEAS(NPF)
cc      COMMON /IDATA/ PERIOD,IORD,ISOD,YEAR,NDAY,IFIX                    
cc      COMMON /RDATA/ ALPHA,BETA,GAMMA,ZER,SMTH,SMTH2,DD,WTRD,DELTA      
cc      COMMON /ILOGT/ LOGT,ISHRNK                                        
      IFLAG=0                                                           
      DMAX0 = 1000.D0                                                   
      DMIN0 = 1.0D0
      MODE = 0                                                          
      RO = 1.41421D0                                                    
      IF(ITRN .NE. 0) RO = DSQRT(RO)                                    
      ND=(N+NF)*2                                                       
      IF(PERIOD .EQ. 1) ND=N+NF                                         
      N7 = NDAY + 6                                                     
      IF(YEAR .NE. 0) ND=ND+N7                                          
      ABICM = 1.D50                                                     
c-----
      DMIN = DMIN0
      ANN = ND
c-----                                                     
C                                                                       
      ALNDT0=0.D0                                                       
      ITRN0 = 30                                                        
      DD=DMIN
      IF(ITRN.EQ.1) DD=5.D0                                             
      DO 9999  IIII=1,ITRN0                                             
C                                                                       
C     BASIC ROUTINE : SEASONAL ADJUSTMENT UNDER GIVEN PRIOR DISTRIBUTION
C                                                                       
C                                                                       
C      CONSTRUCTION AND HOUSEHOLDER TRANSFORMATION OF MATRIX DC         
C                                                                       
      WT=DD*SMTH                                                        
      IF(IIII.GT.1)GO TO 8888                                           
      IF(PERIOD.EQ.1)GO TO 8888                                         
cc      CALL SETDC(H,NH,F,M1,FSEA,N+NF,SMTH2,ZER,IARS,ARFS,IARN,ARFN)     
cc     *        ALPHA,BETA,GAMMA,WTRD,DELTA,PERIOD,IORD,ISOD,YEAR)
      CALL SETDC(H,NH,F,M1,FSEA,N+NF,SMTH2,ZER,IARS,ARFS,IARN,ARFN,
     *        ALPHA,BETA,GAMMA,WTRD,DELTA,PERIOD,ISOD,YEAR)
      ALNDT0=0.D0                                                       
      DO 2233 I=1,M1                                                    
      TEM=DABS(H(1,I))                                                  
 2233 ALNDT0=ALNDT0+DLOG(TEM)                                           
 8888 CONTINUE                                                          
C                                                                       
C     CALCULATION OF LOG(DET(DC'DC))*0.5                                
C                                                                       
C     --------------------                                              
                                                                        
      ALNDTD=(N+NF)*DLOG(WT)+IORD*DLOG(ALPHA)+ALNDT0                    
      IF(PERIOD.NE.1)ALNDTD=ALNDTD+(N+NF)*DLOG(DD)                      
C     --------------------                                              
C                                                                       
C     CONSTRUCTION AND HOUSEHOLDER TRANSFORMATION OF MATRIX DCX         
C                                                                       
 4567 CONTINUE                                                          
cc      CALL SETX(DC,IDC,H2,N2,M1,ICOUNT,FTRN,N+NF,H,NH,WT                
cc     *             ,YS,N,RLIM,WEEK,N7,ALNDTD,F,DD,IART,ARFT)            
      CALL SETX(DC,IDC,H2,N2,M1,ICOUNT,FTRN,N+NF,H,NH,WT,
     *   YS,N,RLIM,WEEK,N7,ALNDTD,F,DD,IART,ARFT,ALPHA,
     *  BETA,GAMMA,WTRD,DELTA,PERIOD,IORD,ISOD,YEAR)
C                                                                       
C     LEAST SQUARES COMPUTATION                                         
C                                                                       
      K=M1 + N2                                                         
      SQE = H2(N2,K)**2                                                 
C                                                                       
C     INTERPRETATION OF THE SOLUTION                                    
C                                                                       
C     ****                                                              
C                                                                       
C                                                                       
C     ABIC COMPUTATION                                                  
C                                                                       
       AN = ICOUNT                                                      
      ANN = ICOUNT + ND                                                 
      ALNDN=0.D0                                                        
      DO 3344 I=1,M1                                                    
      TEM=DABS(DC(1,I))                                                
 3344 ALNDN=ALNDN + DLOG(TEM)                                           
      IF(N2 .EQ. 1) GO TO 3346                                          
      N2M1=N2-1                                                         
      DO 3345 I=1,N2M1                                                  
      TEM=DABS(H2(I,M1+I))                                              
 3345 ALNDN=ALNDN+DLOG(TEM)                                             
 3346 CONTINUE                                                          
      ALSQE=AN*DLOG(SQE/AN)                                             
      ABIC=ALSQE + 2.D0*(ALNDN-ALNDTD)                                  
      IF(YEAR .NE. 0 .AND. WTRD .LE. 0.D0)ABIC=ABIC+N7*2.D0             
cc      WRITE( 6,3 )    DD, ABIC, ALSQE, ALNDN, ALNDTD                    
C                                                                       
C     END OF BASIC ROUTINE                                              
C                                                                       
C                                                                       
C     MINIMUM ABIC PROCEDURE                                            
C                                                                       
      IF(ABIC.GE.ABICM) GO TO 9000                                      
      IF(ABICM-ABIC .LT. 0.0001D0) GO TO 9000                           
      ABICM = ABIC                                                      
      DMIN = DD                                                         
      DD=RO*DD                                                          
      IF(IIII .EQ. 2) MODE = 1                                          
 5678 CONTINUE                                                          
      GO TO 2345                                                        
 9000 IF(MODE .EQ. 1) GO TO 6000                                        
 9001 MODE = 1                                                          
      RO = 1.D0/RO                                                      
      DD = DD*RO*RO                                                     
 2345 IF(DD .LE. DMAX0) GO TO 1234                                      
      IF(MODE .EQ. 0) GO TO 9001                                        
      DD = DMAX0                                                        
      IFLAG=1                                                           
      GO TO 6000                                                        
 1234 CONTINUE                                                          
      IF(DD .GE. DMIN0) GO TO 9999                                      
      DD = DMIN0                                                        
      IFLAG=-1                                                          
      GO TO 6000                                                        
 9999 CONTINUE                                                          
C                                                                       
 6000 CONTINUE                                                          
      DD=DMIN                                                           
      WT = DD*SMTH                                                      
cc      CALL SETX(DC,IDC,H2,N2,M1,ICOUNT,FTRN,N+NF,H,NH,WT                
cc     *             ,YS,N,RLIM,WEEK,N7,ALNDTD,F,DD,IART,ARFT)            
      CALL SETX(DC,IDC,H2,N2,M1,ICOUNT,FTRN,N+NF,H,NH,WT,
     *    YS,N,RLIM,WEEK,N7,ALNDTD,F,DD,IART,ARFT,ALPHA,
     *   BETA,GAMMA,WTRD,DELTA,PERIOD,IORD,ISOD,YEAR)
      NDTEM=ND + 1                                                      
cc      CALL SOLVE(DC,IDC,H2,N2,A,M1,SQE,NDTEM,ERR)                       
      CALL BSOLVE(DC,IDC,H2,N2,A,M1,SQE,NDTEM,ERR)
      SQE=SQE/ANN
      CALL DECODE(SEASON,TREND,EST,ADJ,IRREG,TDC,WEEK0,                 
cc     *            A,YS,N,NF,WEEK,ERR,PSDS,PSDT,SQE)                     
     *              A,YS,N,NF,WEEK,ERR,PSDS,PSDT,SQE,PERIOD,YEAR,NDAY)
      IF(LOGT .EQ. 0) GO TO 6200                                        
      AJACOB=0.D0                                                       
      DO 6100 I=1,N                                                     
 6100 IF(YS(I) .LT. RLIM .OR. RLIM .LE. 0.D0) AJACOB=AJACOB+YS(I)       
      AJACOB=AJACOB+AJACOB                                              
      ABICM=ABICM+AJACOB                                                
 6200 CONTINUE                                                          
      DO 6300 I=1,7                                                     
      WEEK1(I) = WEEK0(I)                                               
      IF(LOGT .NE. 0 .AND. YEAR .NE. 0) WEEK1(I) = DEXP(WEEK1(I))       
 6300 CONTINUE                                                          
    3 FORMAT( 1H ,'    ABIC(',D20.10,'  ) =  ',F10.2,                   
cc     *3X,'ALSQE=',D13.5,3X,'ALNDN=',D13.5,3X,'ALNDTD='D13.5)            
     *3X,'ALSQE=',D13.5,3X,'ALNDN=',D13.5,3X,'ALNDTD=',D13.5)           
      DO 3320 I=1,N                                                     
      DTRN(I)=TREND(I)                                                  
 3320 DSEAS(I)=SEASON(I)                                                
      IF(IORD .EQ. 0) GO TO 3323                                        
      DO 3321 J=1,IORD                                                  
      NMJ=N-J                                                           
      DO 3321 I=1,NMJ                                                   
 3321 DTRN(I)=DTRN(I+1)-DTRN(I)                                         
 3323 IF(ISOD .EQ. 0) GO TO 3325                                        
      DO 3324 J=1,ISOD                                                  
      NMJ=N-J*PERIOD                                                    
      DO 3324 I=1,NMJ                                                   
 3324 DSEAS(I)=DSEAS(I+PERIOD)-DSEAS(I)                                 
 3325 SSTR=0.D0                                                         
      SSEA=0.D0                                                         
      SSIR=0.D0                                                         
      SSAS=0.D0                                                         
      SAS=0.D0                                                          
      IPM1=PERIOD-1                                                     
      DO 3326 I=1,IPM1                                                  
 3326 SAS=SAS+SEASON(I)                                                 
      DO 3400 I=1,N                                                     
      IF(RLIM .GT. 0.D0 .AND. IRREG(I) .GT. RLIM) GO TO 3400            
      SSIR=SSIR+IRREG(I)**2                                             
 3400 CONTINUE                                                          
      NMJ=N-IORD                                                        
      DO 3410 I=1,NMJ                                                   
 3410 SSTR=SSTR+DTRN(I)**2                                              
      NMJ=N-ISOD*PERIOD                                                 
      DO 3420 I=1,NMJ                                                   
      SAS=SAS+SEASON(I+IPM1)                                            
      SSAS=SSAS+SAS**2                                                  
      SAS=SAS-SEASON(I)                                                 
 3420 SSEA=SSEA+DSEAS(I)**2                                             
      APRD=PERIOD                                                       
      SSAS=SSAS/APRD                                                    
cc      WRITE(6,3450) SSIR,SSTR,SSEA,SSAS                                 
 3450 FORMAT(1H ,'SS IRREGULAR =',D13.5,5X,'SS TREND =',D13.5,5X,'SS    
     * SEASONAL =',D13.5,5X,'SS AVSEAS =',D13.5,/,1H )                  
cc      IF(YEAR .NE. 0) WRITE(6,602)  (WEEK1(I),I=1,7)                    
cc      WRITE(6,606)  ABICM, DMIN                                         
  606 FORMAT(1H ,'MINIMUM ABIC =',F10.2,'  ATTAINED AT D =',D20.10 )    
cc      IF(IFLAG .EQ. 1) WRITE(6,600) DMAX0                               
cc      IF(IFLAG .EQ. -1) WRITE(6,601) DMIN0                              
 9876 RETURN                                                            
  600 FORMAT(1H ,'**** D IS HITTING THE UPPER BOUND ',D13.5,' ----- TRY'
     * ,  ' LOWER VALUES OF ORDER AND/OR SORDER')                       
  601 FORMAT(1H ,'**** D IS HITTING THE LOWER BOUND ',D13.5,' ----- TRY'
     * ,   ' HIGHER VALUES OF ORDER AND/OR SORDER')                     
  602 FORMAT(1H ,4X,'MON',9X,'TUE',9X,'WED',9X,'THU',9X,                
     *  'FRI',9X,'SAT',9X,'SUN',/,1H ,7D12.4)                           
  603 FORMAT(1H ,'SHRINKAGE FACTORS ARE  ',D12.5,' FOR TREND, ',        
     *  D12.5,' FOR SEASONAL.')                                         
      END                                                               
      SUBROUTINE CALEND(WEEK,YEAR0,MONTH0,N)                            
C      THIS SUBROUTINE COMPUTES THE DAYS-OF-WEEK DISTRIBUTION OF        
C     N SUCCESSIVE MONTHS STARTING AT MONTH0 OF YEAR0                   
C     NOTE:  THIS SUBROUTINE WORKS FOR YEARS                            
C              AD.1901 - AD.2099                                        
C                                                                       
      IMPLICIT INTEGER*4 (A-Z)                                          
      REAL*8 WEEK(7,1),W0(8)                                            
C                                                                       
      DYEAR=(MONTH0-1)/12                                               
      IF(MONTH0 .GE. 1) GO TO 20                                        
      DYEAR=-MONTH0                                                     
      DYEAR=DYEAR/12+1                                                  
      DYEAR=-DYEAR                                                      
   20 CONTINUE                                                          
      MONTH=MONTH0-12*DYEAR                                             
      YEAR=YEAR0+DYEAR                                                  
      LEAP=MOD(YEAR,4)                                                  
C  DAY-OF-WEEK OF THE FIRST DAY OF 'YEAR'                               
      Y = YEAR-1901                                                     
      L=Y/4                                                             
      D=Y+L+2                                                           
      DAY=MOD(D,7) + 1                                                  
C  DAY-OF-WEEK OF THE FIRST DAY OF 'MONTH' OF 'YEAR'                    
      GO TO (200,203,203,206,201,204,206,202,205,200,203,205), MONTH    
  201 DAY=DAY+1                                                         
      GO TO 200                                                         
  202 DAY=DAY+2                                                         
      GO TO 200                                                         
  203 DAY=DAY+3                                                         
      GO TO 200                                                         
  204 DAY=DAY+4                                                         
      GO TO 200                                                         
  205 DAY=DAY+5                                                         
      GO TO 200                                                         
  206 DAY=DAY+6                                                         
  200 IF(LEAP .EQ. 0 .AND. MONTH .GE. 3) DAY=DAY+1                      
      IF(DAY .GT. 7) DAY=DAY-7                                          
C  ITERATION                                                            
      DO 100 I=1,N                                                      
      DO 10 J=1,7                                                       
   10 W0(J) = 4.D0                                                      
      GO TO (331,328,331,330,331,330,331,331,330,331,330,331), MONTH    
  331 DIFF=3                                                            
      W0(8)=31.D0                                                       
      GO TO 300                                                         
  330 DIFF=2                                                            
      W0(8)=30.D0                                                       
      GO TO 300                                                         
  328 DIFF=0                                                            
      W0(8)=28.D0                                                       
      IF(LEAP .NE. 0) GO TO 50                                          
      DIFF=1                                                            
      W0(8)=29.D0                                                       
  300 WDAY=8-DAY                                                        
      DO 400 J=1,DIFF                                                   
      W0(WDAY) = 5.D0                                                   
      IF(J .EQ. DIFF) GO TO 50                                          
      WDAY=WDAY-1                                                       
  400 IF(WDAY .EQ. 0) WDAY=7                                            
   50 CONTINUE                                                          
      DO 410 J=1,7                                                      
  410 WEEK(J,I)=W0(J)-30.4375D0/7.D0                                    
      IF(I .EQ. N) GO TO 900                                            
      DAY = DAY + DIFF                                                  
      IF(DAY .GT. 7) DAY = DAY - 7                                      
      MONTH=MONTH+1                                                     
      IF(MONTH .LE. 12) GO TO 100                                       
      MONTH=1                                                           
      YEAR=YEAR+1                                                       
      LEAP=MOD(YEAR,4)                                                  
  100 CONTINUE                                                          
  900 RETURN                                                            
      END                                                               
C                                                                       
CC      FUNCTION AMAX(A,N)                                                
      DOUBLE PRECISION FUNCTION AMAX(A,N)
C     COMMON SUBROUTINE                                                 
C     MAXIMUM OF A(I)(I=1,N) SEARCH                                     
      DOUBLE PRECISION A
      DIMENSION A(N)                                                    
      AMAX=A(1)                                                         
      DO 10 I=2,N                                                       
      IF(AMAX.LT.A(I)) AMAX=A(I)                                        
   10 CONTINUE                                                          
      RETURN                                                            
      END                                                               
C                                                                       
CC      FUNCTION AMIN(A,N)                                                
      DOUBLE PRECISION FUNCTION AMIN(A,N)
C     COMMON SUBROUTINE                                                 
C     MINIMUM OF A(I)(I=1,N) SEARCH                                     
      DOUBLE PRECISION A
      DIMENSION A(N)                                                    
      AMIN=A(1)                                                         
      DO 10 I=2,N                                                       
      IF(AMIN.GT.A(I)) AMIN=A(I)                                        
   10 CONTINUE                                                          
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE  POOLAV( Z,K,X,SD1 )                                   
C                                                                       
C     THIS SUBROUTINE SEARCHES FOR THE MINIMUM OF                       
C          F(X;Z) = (Z(1)-X(1))**2 + ... + (Z(K)-X(K))**2               
C     SUBJECT TO X(1)<X(2)< ... < X(K), BY THE POOL-ADJACENT-VIOLALORS  
C     ALGORITHM.                                                        
C                                                                       
C     INPUTS:                                                           
C        (Z(I),I=1,K): DATA                                             
C        K:            NUMBER OF DATA                                   
C     OUTPUTS:                                                          
C        (X(I),I=1,K): VECTOR OF MINIMIZING SOLUTION                    
C        SD1:          MINIMUM OF F(X;Z)                                
C                                                                       
      IMPLICIT  REAL*8( A-H,O-Z )                                       
cc      DIMENSION  Z(K), X(K), Y(20)                                      
      DIMENSION  Z(K), X(K), Y(K)                                      
      DO 10  I=1,K                                                      
   10 X(I) = Z(I)                                                       
C                                                                       
  100 DO 20  I=2,K                                                      
   20 IF( X(I-1) .GT. X(I) )  GO TO 30                                  
      GO TO 300                                                         
C                                                                       
   30 IFG = 0                                                           
      DO 40  I=1,K                                                      
   40 Y(I) = X(I)                                                       
C
      N0=1                                                                       
      DO 200  I=1,K-1
      I0 = I
C                                                                       
      IF( X(I) .LT. X(I+1) )  GO TO 110                                 
      IF( I .EQ. K-1  .AND.  IFG .NE. 0 )  GO TO 90                     
      IF( IFG .NE. 0 )  GO TO 200                                       
      IFG = 1                                                           
      N0 = I                                                            
      IF( I .EQ. K-1 )  GO TO 90                                        
      GO TO 200                                                         
CC   90 I = K
   90 I0 = K
      GO TO 115                                                         
C                                                                       
  110 IF( IFG .EQ. 0 )  GO TO 200                                       
      IFG = 0                                                           
  115 CONTINUE                                                          
      SUM = 0.D0                                                        
CC      DO 120  J=N0,I
      DO 120 J=N0,I0
  120 SUM = SUM + Y(J)                                                  
CC      SUM = SUM/(I-N0+1)
      SUM = SUM/(I0-N0+1)
CC      DO 130  J=N0,I
      DO 130  J=N0,I0                                             
  130 Y(J) = SUM                                                        
C                                                                       
  200 CONTINUE                                                          
C                                                                       
      DO 210  I=1,K                                                     
  210 X(I) = Y(I)                                                       
      GO TO 100                                                         
  300 SD1 = 0.D0                                                        
      DO 310  I=1,K                                                     
  310 SD1 = SD1 + (X(I) - Z(I))**2                                      
    2 FORMAT( 1H ,10F13.5 )                                             
      RETURN                                                            
      E N D                                                             
      SUBROUTINE  LKOUT1( X,N,IND,JSW,F,W )                             
C                                                                       
C     THIS SUBROUTINE COMPUTES THE LIKELIHOOD OF THE MODEL THAT (X(I);  
C     IND(I)=1) ARE THE OUTLIERS.  (MEAN SLIPPAGE TYPE MODEL)           
C                                                                       
C     INPUTS:                                                           
C        (X(I),I=1,N): OBSERVATIONS                                     
C        N:            NUMBER OF OBSERVATIONS                           
C        (IND(I),I=1,N): = 0 ; IF X(I) IS A "NORMAL" OBSERVATION.       
C                        = 1 ; IF X(I) IS AN OUTLIER.                   
C        JSW:   =0;    WHOLE MODELS ARE EVALUATED                       
C               =1;    ONLY THE NATURALLY ORDERED MODEL IS EVALUATED    
C                      (SIMPLIFIED ALGORITHM)                           
C                                                                       
C     OUTPUTS:                                                          
C        F:            LIKELIHOOD OF THE MODEL                          
C        W:                                                             
C                                                                       
      IMPLICIT REAL*8  ( A-H,O-Z )                                      
cc      DIMENSION  X(N), Y(10), Z(10), ZE(10), IND(N), JND(10)            
      DIMENSION  X(N), Y(N), Z(N), ZE(N), IND(N), JND(N)            
C                                                                       
      L = 0                                                             
      SUM = 0.D0                                                        
      DO 10  I=1,N                                                      
      IF( IND(I) .EQ. 1 )   GO TO 10                                    
      L = L + 1                                                         
      SUM = SUM+X(I)                                                    
   10 CONTINUE                                                          
      XMEAN = SUM/DFLOAT(L)                                             
      K = N-L                                                           
C                                                                       
      SUM = 0.D0                                                        
      DO 20  I=1,N                                                      
   20 IF( IND(I) .EQ. 0 )   SUM = SUM+(X(I)-XMEAN)**2                   
      SIG2 = SUM/N                                                      
      W = 1.D0                                                          
      F = -0.5D0*N*DLOG(SIG2)                                           
C                                                                       
      IF( JSW .EQ.1 )   RETURN                                          
      IF( K .LE. 1 )   RETURN                                           
C                                                                       
      J = 0                                                             
      DO 30  I=1,N                                                      
      IF( IND(I) .EQ. 0 )   GO TO 30                                    
      J = J+1                                                           
      Y(J) = X(I)                                                       
   30 CONTINUE                                                          
C                                                                       
      W = 0.D0                                                          
      DO 40  I=1,K                                                      
   40 JND(I) = I                                                        
C                                                                       
   50 DO 60  I=1,K                                                      
      J = JND(I)                                                        
   60 Z(I) = Y(J)                                                       
      CALL  POOLAV( Z,K,ZE,SD )                                         
C                                                                       
      W = W + 1.D0/DSQRT(1.D0+SD/SUM)**N                                
C                                                                       
      CALL  PERMUT( JND,K,IFG )                                         
      IF( IFG .EQ. 0 )   GO TO 50                                       
C                                                                       
      RETURN                                                            
  600 FORMAT( 1H ,'IND',40I3 )                                          
  601 FORMAT( 1H ,'F =',D13.5,5X,'FSUM =',D13.5 )                       
      E N D                                                             
cc      SUBROUTINE  PRPOST( POST,X,IND,JND,KND,IC,N,L)                    
      SUBROUTINE  PRPOST( POST,X,IND,JND,KND,IC,N,L)
C                                                                       
C     THIS SUBROUTINE ARRANGES POST(I), JND(I) AND KND(I) (I=1,IC) IN   
C     DECREASING ORDER OF POST(I), AND DRAWS POSTERIOR PROBABILITY AND  
C     THE ASSUMED OUTLIERS OF THE MODEL WITH POSTERIOR PROBABILITY      
C     GREATER THAN EPS.                                                 
C                                                                       
C     INPUTS:                                                           
C        (POST(I),I=1,IC):   POSTERIOR PROBABILITIES OF THE MODELS      
C        (X(I),I=1,N):       ORIGINAL DATA                              
C        (IND(I),I=1,N):     WORK AREA                                  
C        (JND(I),I=1,IC):    SPECIFICATION OF THE OUTLIERS IN LOW SIDE  
C                            (CODED IN DECIMAL)                         
C        (KND(I),I=1,IC):    SPECIFICATION OF THE OUTLIERS IN HIGH SIDE 
C                            (CODED IN DECIMAL)                         
C        IC:                 NUMBER OF RECORDED MODELS                  
C        N:                  NUMBER OF ORIGINAL DATA                    
C        L:                  NUMBER OF POSSIBLE OUTLIERS IN BOTH SIDES  
C        EPS:                LOWEST LIMIT OF POSTERIOR PROBABILITY TO BE
C                            PRINTED                                    
C                                                                       
      REAL*8  POST, X(N)                                                
cc      DIMENSION  POST(IC), JND(IC),KND(IC), IND(N), Y(10)               
      DIMENSION  POST(IC), JND(IC),KND(IC), IND(N), Y(N)
cc      COMMON /CSPRSS/ ISPRSS                                            
C                                                                       
      DO 20  I=1,IC                                                     
      IMAX = I                                                          
      PMAX = POST(I)                                                    
      DO 10  J=I,IC                                                     
      IF( POST(J) .LE. PMAX )   GO TO 10                                
      IMAX = J                                                          
      PMAX = POST(J)                                                    
   10 CONTINUE                                                          
      IF( IMAX .EQ. I )   GO TO 20                                      
      POST(IMAX) = POST(I)                                              
      POST(I) = PMAX                                                    
      JJ = JND(I)                                                       
      KK = KND(I)                                                       
      JND(I) = JND(IMAX)                                                
      KND(I) = KND(IMAX)                                                
      JND(IMAX) = JJ                                                    
      KND(IMAX) = KK                                                    
   20 CONTINUE                                                          
   30 IC1 = IC                                                          
      NML1 = N-L+1                                                      
      DO 40  I=1,N                                                      
   40 IND(I) = 0                                                        
C                                                                       
cc      IF(ISPRSS .EQ. 0) WRITE( 6,4 )                                    
      DO 100  J=1,IC1                                                   
      CALL  BINARY( JND(J),L,IND )                                      
      CALL  BINARY( KND(J),L,IND(NML1) )                                
      ID = 0                                                            
      DO 50  I=1,N                                                      
      IF( IND(I) .EQ. 0 )    GO TO 50                                   
      ID = ID+1                                                         
      Y(ID) = X(I)                                                      
   50 CONTINUE                                                          
cc      IF(ISPRSS .NE. 0) GO TO 100                                       
cc      IF( ID .GE. 1 )   WRITE( 6,5 )   J, POST(J), (Y(I),I=1,ID)        
cc      IF( ID .EQ. 0 )   WRITE( 6,6 )   J, POST(J)                       
  100 CONTINUE                                                          
      RETURN                                                            
    4 FORMAT( 1H ,10X,'POSTERIOR',10X,'OUTLIERS' )                      
    5 FORMAT( 1H ,I5,F13.5,5X,10F10.3 )                                 
    6 FORMAT( 1H ,I5,F13.5,9X,'NONE' )                                  
      E N D                                                             
      SUBROUTINE PERMUT( IND,K,IFG )                                    
C                                                                       
C     THIS SUBROUTINE SEQUENTIALLY SPECIFIES K] CONFIGURATIONS (IND(1), 
C     ...,IND(K)) OBTAINED BY PERMUTING (1,...,K)                       
C                                                                       
C     INPUTS:                                                           
C        (IND(I),I=1,K): FORMER CONFIGURATION                           
C        K:              NUMBER OF ELEMENTS TO BE PERMUTED              
C                                                                       
C     OUTPUTS:                                                          
C        (IND(I),I=1,K): NEW CONFIGURATION                              
C        IFG:            = 0 ; IF THE NEW CONFIGURATION IS OBTAINED     
C                        = 1 ; SEARCH FOR THE CONFIGURATION COMPLETED.  
C                                                                       
      DIMENSION  IND(K)                                                 
C                                                                       
      I1 = 1                                                            
      I2 = 2                                                            
      I0 = 1                                                            
      IMAX = IND(1)                                                     
      IFG = 0                                                           
C                                                                       
   10 IF( IND(I1) .LT. IND(I2) )   GO TO 100                            
      IF( I2 .EQ. I1+1 )   GO TO 20                                     
      I1 = I1+1                                                         
      GO TO 10                                                          
C                                                                       
   20 I2 = I2+1                                                         
      IF( I2 .GT. K )   GO TO 200                                       
C                                                                       
      I2M1 = I2-1                                                       
      DO 30  I=1,I2M1                                                   
   30 IF( IND(I) .LE. IND(I2) )   GO TO 40                              
      GO TO 20                                                          
C                                                                       
   40 IMAX = 0                                                          
      DO 50  I=1,I2M1                                                   
      IF( IND(I) .GT. IND(I2) )   GO TO 50                              
      IF( IND(I) .LT. IMAX )   GO TO 50                                 
      IMAX = IND(I)                                                     
      I0 = I                                                            
   50 CONTINUE                                                          
C                                                                       
  100 IND(I0) = IND(I2)                                                 
      IND(I2) = IMAX                                                    
      I2M1 = I2-1                                                       
      IF( I2 .GT. 2 )   CALL  ISORT( IND,I2M1 )                         
      RETURN                                                            
C                                                                       
  200 IFG = 1                                                           
      RETURN                                                            
C                                                                       
      E N D                                                             
      SUBROUTINE  ISORT( IND,N )                                        
C                                                                       
C     THIS SUBROUTINE ARRANGES IND(I) (I=1,N) IN ORDER OF INCREASING    
C     MAGNITUDE                                                         
C                                                                       
C     INPUTS:                                                           
C        (IND(I),I=1,N): ORIGINAL DATA                                  
C        N:              NUMBER OF DATA                                 
C                                                                       
C     OUTPUT:                                                           
C        (IND(I),I=1,N): REORDERED DATA                                 
C                                                                       
      DIMENSION  IND(N)                                                 
C                                                                       
      NM1 = N-1                                                         
      DO 20  II=1,NM1                                                   
      MINI = IND(II)                                                    
      IMIN = II                                                         
      DO 10  I=II,N                                                     
      IF( MINI .LE. IND(I) )   GO TO 10                                 
      MINI = IND(I)                                                     
      IMIN = I                                                          
   10 CONTINUE                                                          
      IF( IMIN .EQ. II )   GO TO 20                                     
      J = IND(II)                                                       
      IND(II) = MINI                                                    
      IND(IMIN) = J                                                     
   20 CONTINUE                                                          
      RETURN                                                            
      E N D                                                             
      SUBROUTINE  BINARY( M,K,MB )                                      
C                                                                       
C       DECIMAL TO BINARY CONVERSION                                    
C                                                                       
C       INPUTS:                                                         
C          M:     NUMBER IN DECIMAL REPRESENTATION                      
C          K:     NUMBER OF BITS USED FOR THE BINARY REPRESENTATION     
C                                                                       
C       OUTPUT:                                                         
C          MB:    NUMBER IN BINARY REPRESENTATION                       
C                                                                       
      DIMENSION  MB(1)                                                  
C                                                                       
      N = M                                                             
      DO 10  I=1,K                                                      
        L = N / 2                                                       
        MB(I) = N - L*2                                                 
   10 N = L                                                             
      RETURN                                                            
C                                                                       
      E N D                                                             
cc      SUBROUTINE  SRTMIN( X,N,IX )                                      
      SUBROUTINE  BSRTMIN( X,N,IX )
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
cc      SUBROUTINE OUTLIR( Z,NN,K,ISW,JSW,Y,RLIM )                        
      SUBROUTINE OUTLIR( Z,NN,K,ISW,JSW,Y,RLIM,IOUTD,ROUT )
C     REVISED  MAY 30, 1980                                             
C                                                                       
C     INPUTS:                                                           
C        (Z(I),I=1,N):   ORIGINAL DATA                                  
C        NN:             NUMBER OF DATA                                 
C        K:              MAXIMUM OF THE NUMBER OF OUTLIERS              
C        L:              NUMBER OF POSSIBLE OUTLIERS IN BOTH SIDES.     
C        ISW   =0:       MODIFIED DATA ARE NOT GIVEN                    
C              =2:       THE OBSERVATIONS JUDGED AS OUTLIERS ARE        
C                        REPLACED BY A CONSTANT.                        
C        JSW   =0:       WHOLE MODELS ARE EVALUATED                     
C              =1:       ONLY THE NATURALLY ORDERED MODEL IS EVALUATED  
C                        (SIMPLIFIED ALGORITHM)                         
C        RLIM:           ORIGINAL DATA WHOSE VALUES ARE GREATER         
C                        THAN OR EQUAL TO RLIM ARE TREATEDAS            
C                        MISSING OBSERVATIONS.                          
C     OUTPUTS:                                                          
C        (IX(I),I=1,N):  SUBSCRIPT INDICATING THE ORDER OF THE MAGNITUDE
C                        OF ORIGINAL DATA.                              
C                        I.E.,  Z(IX(1))<Z(IX(2))< ... <Z(IX(N)).       
C                        NOTE THAT Z(IX(I))=X(I).                       
C        (PM(I),I=1,N):  MARGINAL POSTRIOR PROBABILITY THAT X(I) IS AN O
C        (POST(I),I=1,IC):   POSTERIOR PROBABILITIES OF THE MODELS      
C        (JND(I),I=1,IC):    SPECIFICATION OF THE OUTLIERS IN LOW SIDE  
C                            (CODED IN DECIMAL)                         
C        (KND(I),I=1,IC):    SPECIFICATION OF THE OUTLIERS IN HIGH SIDE 
C                            (CODED IN DECIMAL)                         
C        IC:                 NUMBER OF RECORDED MODELS                  
C                                                                       
      IMPLICIT  REAL*8  ( A-H,O-Z )                                     
cc      DIMENSION  X(500), F(501), PM(500), Z(500), C(20)                 
cc      DIMENSION  IX(500) , IND(500) , JND(1000) , KND(1000) , POST(1000)
cc      DIMENSION  Y(500)                                                 
      DIMENSION  X(NN), F(NN+1), PM(NN), Z(NN), C(K+1)                 
      DIMENSION  IX(NN) , IND(NN) , JND(2**K) , KND(2**K) , POST(2**K)
      DIMENSION  Y(NN)                                                 
cc      COMMON /CSPRSS/ ISPRSS                                            
      ISPRSS = 1                                                        
C                                                                       
      N=0                                                               
      DO 5 I=1,NN                                                       
      IF(RLIM .LE. 0.D0) GO TO 6                                        
      IF(Z(I) .GE. RLIM) GO TO 7                                        
    6 N=N+1                                                             
      X(N)=Z(I)                                                         
      IX(N)=I                                                           
    7 Y(I)=Z(I)                                                         
    5 CONTINUE                                                          
      L=K                                                               
cc      IF(ISPRSS .EQ. 0) WRITE( 6,600 )   N, K, L, ISW, JSW              
cc      IF(ISPRSS .EQ. 0) WRITE( 6,601 )   (X(I),I=1,N)                   
C                                                                       
cc      CALL  SRTMIN( X,N,IX )                                            
      CALL  BSRTMIN( X,N,IX )
      EPS = 1.0D-03                                                     
      NML1 = N-L+1                                                      
      F(1) = 0.D0                                                       
      DO 10  I=1,N                                                      
      IND(I) = 0                                                        
      PM(I) = 0.D0                                                      
      DI = I                                                            
   10 F(I+1) = F(I)+DLOG(DI)                                            
      C(1) = DFLOAT(N*2)/DFLOAT(N-3)                                    
      DO 20  I=1,K                                                      
   20 C(I+1) = DFLOAT(N*(I+2))/DFLOAT(N-I-3)+F(N+1)-F(N-I+1)            
C                                                                       
      IL = 2**L                                                         
      IF(JSW.EQ.1)  IL=K+1                                              
      IC = 0                                                            
      SUMF = 0.D0                                                       
      DLK0 = 0.D0
C                                                                       
      DO 101  II1 = 1,IL                                                
      II = II1-1                                                        
      IF(JSW.EQ.1)  II=2**II-1                                          
      CALL  BINARY( II,L,IND )                                          
      K1 = 0                                                            
      DO 30  I=1,L                                                      
   30 K1 = K1+IND(I)                                                    
      IF( K1 .GT. K )   GO TO 101                                       
C                                                                       
      DO 100  JJ1=1,IL                                                  
      JJ = JJ1-1                                                        
      IF(JSW.EQ.1)   JJ=2**K-2**JJ                                      
      CALL  BINARY( JJ,L,IND(NML1) )                                    
      K2 = K1                                                           
      DO 40  I=NML1,N                                                   
   40 K2 = K2+IND(I)                                                    
      IF( K2 .GT. K )   GO TO 100                                       
C                                                                       
      CALL  LKOUT1( X,N,IND,JSW,FF,W )                                  
      IF( IC.EQ.0 ) DLK0=FF-C(K2+1)                                     
C                                                                       
      F0 = FF-C(K2+1)-DLK0                                              
      IF( F0 .LT. -20.D0 )  GO TO 100                                   
      IF(F0.LT.20.D0)GO TO 45                                           
      DLK0=F0                                                           
      TEM=DEXP(-F0)                                                     
      SUMF=SUMF*TEM                                                     
      DO 46 I=1,N                                                       
   46 PM(I)=PM(I)*TEM                                                   
      DO 47 I=1,IC                                                      
   47 POST(I)=POST(I)*TEM                                               
      F0=0.D0                                                           
   45 CONTINUE                                                          
      EXPF = DEXP( F0 )*W                                               
      SUMF = SUMF+EXPF                                                  
      DO 50  I=1,N                                                      
   50 PM(I) = PM(I) + IND(I)*EXPF                                       
      IF( EXPF/SUMF .LT. EPS )   GO TO 100                              
C                                                                       
      IC = IC+1                                                         
      JND(IC) = II                                                      
      KND(IC) = JJ                                                      
      POST(IC) = EXPF                                                   
C                                                                       
  100 CONTINUE                                                          
C                                                                       
  101 CONTINUE                                                          
      DO 110  I=1,N                                                     
  110 PM(I) = PM(I)/SUMF                                                
C                                                                       
      DO 120  I=1,IC                                                    
  120 POST(I) = POST(I)/SUMF                                            
C                                                                       
      CALL  PRPOST( POST,X,IND,JND,KND,IC,N,L)                          
C                                                                       
cc      IF(ISPRSS .EQ. 0) WRITE( 6,602 )   (PM(I),I=1,N)                  
cc      IF( ISW .GE. 1 )   CALL  MODIFY( N,L,IX,PM,JND,KND,Y,IC )         
      IF( ISW .GE. 1 )  CALL MODIFY( N,L,IX,PM,JND,KND,Y,IC,IOUTD,ROUT )
cc      IF(ISW .GE. 1 .AND. ISPRSS .EQ. 0) WRITE(6,603)(Y(I),I=1,NN)      
C                                                                       
      RETURN                                                            
  600 FORMAT( 1H ,'N    =',I6,5X,'(NUMBER OF DATA)',/,                  
     1' K    =',I6,5X,'(MAXIMUM NUMBER OF OUTLIERS)',/,                 
     2' L    =',I6,5X,'(RANGE OF SEARCH ON BOTH SIDES)',/,              
     3' ISW  =',I6,/,' JSW  =',I6 )                                     
  601 FORMAT( 1H ,'**  DATA  **',/,(1X,12D11.3) )                       
  602 FORMAT( 1H ,'**  MARGINAL POSTERIOR PROBABILITIES  **',/,(1X,12D11
     *.3) )                                                             
  603 FORMAT(1H ,'** MODIFIED DATA **',/,(1X,12D11.3) )                 
      E N D                                                             
cc      SUBROUTINE  MODIFY( N,L,IX,POST,JND,KND,Y,IC )                    
      SUBROUTINE  MODIFY( N,L,IX,POST,JND,KND,Y,IC,IOUTD,CONST )
C     REVISED  MAY 29, 1980                                             
C                                                                       
C     THIS SUBROUTINE MODIFIES THE ORIGINAL DATA BY USING THE OUTPUTS   
C     FROM SUBROUTINE OUTLIR.                                           
C                                                                       
C     INPUTS:                                                           
C        (X(I),I=1,N):   ORIGINAL DATA                                  
C        N:              NUMBER OF ORIGINAL DATA                        
C        L:              NUMBER OF POSSIBLE OUTLIERS IN BOTH SIDES      
C        (IX(I),I=1,N):                                                 
C        (POST(I),I=1,N): MARGINAL POSTERIOR PROBABILITIES THAT X(IX(I))
C                        IS AN OUTLIER                                  
C        IOUTD:          =1  X(IX(I)) IS JUDGED AS AN OUTLIER           
C                            WHEN POST(I) IS GREATER THAN 0.01          
C                        =2  X(IX(I)) IS JUDGED AS AN OUTLIER           
C                            WHEN AT LEAST ONE MODEL, WHOSE             
C                            POSTERIOR PROBABILITY IS GREATER THAN      
C                            THAT OF NO-OUTLIER MODEL, SPECIFIES IT     
C                            AS AN OUTLIER                              
C                                                                       
C     OUTPUT:                                                           
C        (Y(I),I=1,N):   MODIFIED DATA                                  
C                                                                       
      IMPLICIT  REAL*8  ( A-H,O-Z )                                     
      DIMENSION   IX(N),  POST(1), JND(1), KND(1), Y(N)                 
cc      DIMENSION  IND(500)                                               
      DIMENSION  IND(N)                                               
cc      COMMON /CSPRSS/ ISPRSS                                            
cc      COMMON /ILOGT/ IDUMMY(3),IOUTD,CONST                              
      ICTEM = IC                                                        
      IF(IOUTD .EQ. 1) ICTEM=1                                          
      DO 110 I=1,N                                                      
  110 IND(I) = 0                                                        
      NML1 = N-L+1                                                      
      DO 200 K=1,ICTEM                                                  
      CALL  BINARY( JND(K),L,IND )                                      
      CALL  BINARY( KND(K),L,IND(NML1) )                                
C                                                                       
      ICHK = 0                                                          
      DO 120  I=1,N                                                     
      J = IX(I)                                                         
      IF(IOUTD .EQ. 1 .AND. POST(I) .LE. 0.01D0) GO TO 120              
      IF(IOUTD .EQ. 2 .AND. IND(I) .EQ. 0) GO TO 120                    
      ICHK=1                                                            
      Y(J) = CONST                                                      
  120 CONTINUE                                                          
      IF(ICHK .EQ. 0) GO TO 201                                         
  200 CONTINUE                                                          
  201 CONTINUE                                                          
      RETURN                                                            
  600 FORMAT( 1H ,'**  MODIFIED DATA (ISW=1)  **',/,(1X,10F13.5) )      
  601 FORMAT( 1H ,'**  MODIFIED DATA (ISW=2)  **',/,(1X,10F13.5) )      
  602 FORMAT( 1H ,'**  MODIFIED DATA (ISW=3)  **',/,(1X,10F13.5) )      
  605 FORMAT( 1H ,'LOCATION PARAMETER;   XM =',F13.5 )                  
      E N D                                                             
      SUBROUTINE SETD(W,IP,ID,C,IAR,AR)                                 
      IMPLICIT REAL*8 (A-H,O-Z)                                         
cc      DIMENSION W(IP,10), WW(10), AR(10)                                
      DIMENSION W(IP,ID+IAR+1), WW(ID+IAR+1), AR(1)
      IDAR = ID + IAR                                                   
      IDP1 = IDAR + 1                                                   
      W(1,IDP1) = C                                                     
      WW(IDP1) = C                                                      
      IF(IDAR .EQ. 0) GO TO 998                                         
      DO 10 J=1,IDAR                                                    
      WW(J)=0.D0                                                        
      DO 10 I=1,IP                                                      
   10 W(I,J) = 0.D0                                                     
      IF(ID .EQ. 0) GO TO 997                                           
      DO 50 J=1,ID                                                      
      I = IDP1 - J - 1                                                  
      DO 50 K=1,J                                                       
      I = I + 1                                                         
   50 WW(I) = WW(I) - WW(I+1)                                           
  997 CONTINUE                                                          
      DO 500 J=1,IDAR                                                   
      W(1,J) = WW(J)                                                    
      JLX = MIN0(IAR,IDP1 - J)                                          
      IF(IAR .EQ. 0) GO TO 500                                          
      DO 400 JY=1,JLX                                                   
  400 W(1,J) = W(1,J) - AR(JY)*WW(J+JY)                                 
  500 CONTINUE                                                          
  998 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE  INIT(W,LENGTH,DOP,ISTEP)                              
      IMPLICIT REAL*8 (A-H,O-Z)                                         
cc      DIMENSION  W(100), DDOP(100), DOP(100)                            
      DIMENSION  W(LENGTH), DDOP(LENGTH), DOP((LENGTH-1)*ISTEP+1)
      J=1                                                               
      DO 1 I=1,LENGTH                                                   
      DDOP(I)=DOP(J)                                                    
    1 J=J+ISTEP                                                         
      DO 20 J=1,LENGTH                                                  
      SUM = 0.D0                                                        
      ITEM = 0                                                          
      DO 10 K=J,LENGTH                                                  
      ITEM=ITEM+1                                                       
   10 SUM=SUM-W(K)*DDOP(ITEM)                                           
   20 W(J)=SUM                                                          
      RETURN                                                            
      END                                                               
C                                                                       
C                                                                       
C                                                                       
      SUBROUTINE  EXHSLD(H1,N1,H2,N2,H3,N3,H4,M1,IPOS)                  
      IMPLICIT REAL*8(A-H,O-Z)                                          
      DIMENSION  H1(N1,1),H2(N2,1),H3(1),H4(1)                          
      DATA EPS/1.D-30/                                                  
      IF(IPOS .LE. M1) GO TO 30                                         
      M1 = IPOS                                                         
      DO 10 J=1,N1                                                      
   10 H1(J,M1) = 0.D0                                                   
      DO 20 J=1,N2                                                      
   20 H2(J,M1+N2) = 0.D0                                                
   30 CONTINUE                                                          
      IF(N3 .LT. 0) RETURN                                              
      M0 = IPOS - N3                                                    
      DO 100 J=1,N3                                                     
      MM = M0 + J                                                       
      IF(DABS(H3(J)).LT.EPS)GO TO 100                                   
      D = H1(1,MM)**2 + H3(J)**2                                        
      SQD = DSQRT(D)                                                    
      IF(H1(1,MM).GT.0.D0) SQD=-SQD                                     
      C = D - SQD*H1(1,MM)                                              
      D = H1(1,MM) - SQD                                                
      H1(1,MM) = SQD                                                    
      JP1=J+1                                                           
      IF(JP1.GT.N3)GO TO 60                                             
      M = 1                                                             
      DO 50 K=JP1,N3                                                    
      M = M + 1                                                         
      IF(M .GT. N1) GO TO 60                                            
      F = D*H1(M,MM) + H3(J)*H3(K)                                      
      F = F/C                                                           
      H1(M,MM) = H1(M,MM) - D*F                                         
   50 H3(K) = H3(K) - H3(J)*F                                           
   60 CONTINUE                                                          
      DO 70 K=1,N2                                                      
      F = D*H2(K,MM) + H3(J)*H4(K)                                      
      F = F/C                                                           
      H2(K,MM) = H2(K,MM) - D*F                                         
   70 H4(K) = H4(K) - H3(J)*F                                           
  100 CONTINUE                                                          
      DO 200 J=1,N2                                                     
      MM = M1 + J                                                       
      IF(DABS(H4(J)).LT.EPS)GO TO 200                                   
      D = H2(J,MM)**2 + H4(J)**2                                        
      SQD = DSQRT(D)                                                    
      IF(H2(J,MM).GT.0.D0) SQD=-SQD                                     
      C = D - SQD*H2(J,MM)                                              
      D = H2(J,MM) - SQD                                                
      H2(J,MM) = SQD                                                    
      IF(J .GE. N2) GO TO 200                                           
      JP1 = J + 1                                                       
      DO 150 K=JP1,N2                                                   
      F = D*H2(K,MM) + H4(J)*H4(K)                                      
      F = F/C                                                           
      H2(K,MM) = H2(K,MM) - D*F                                         
  150 H4(K) = H4(K) - H4(J)*F                                           
  200 CONTINUE                                                          
      RETURN                                                            
      END                                                               
C                                                                       
C                                                                       
C                                                                       
cc      SUBROUTINE SETX(H1,N1,H2,N2,M1,ICOUNT,FTRN,N,H,NH,WT,Y,NDATA,RLIM,
cc     *                 WEEK,IY,ALNDTD,F,DD,IART,ARFT)                   
      SUBROUTINE SETX(H1,N1,H2,N2,M1,ICOUNT,FTRN,N,H,NH,WT,Y,
     *                   NDATA,RLIM,WEEK,IY,ALNDTD,F,DD,IART,ARFT,
     *                   ALPHA,BETA,GAMMA,WTRD,DELTA,IP,ID,IS,YEAR)
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      INTEGER*4 YEAR                                                    
cc      DIMENSION H1(N1,1),H2(N2,1),Y(1),WEEK(IY,1),H3(200),H4(50)
cc      DIMENSION FTRN(1),H(NH,1),TI(10),T0(10),F(1),ARFT(1)              
      DIMENSION H1(N1,1),H2(N2,1),Y(1),WEEK(IY,1),H3(N1),H4(N2)
      DIMENSION FTRN(1),H(NH,1),TI(ID+IART),T0(2*(ID+IART+1)),
     *           F(1),ARFT(1)
cc      COMMON/IDATA/IP,ID,IS,YEAR                                        
cc      COMMON/ RDATA/ALPHA,BETA,GAMMA,DUMMY(4),WTRD,DELTA                
      DO 600 I=1,N2                                                     
  600 H2(I,N2)=0.D0                                                     
      IF(YEAR .EQ. 0) GO TO 200                                         
      TEMO=WTRD                                                         
      TEM=-TEMO/7.D0                                                    
      DO 300 J=1,7                                                      
      DO 400 I=1,7                                                      
  400 H2(I,J)=TEM                                                       
      H2(8,J)=TEM*DELTA                                                 
  300 H2(J,J)=H2(J,J)+TEMO                                              
cc      CALL HUSHLD(H2,8,8,8,0)                                           
      CALL BHUSHLD(H2,8,8,8,0)
      DO 500 I=1,N2                                                     
      DO 500 J=1,I                                                      
      TEM=H2(I,J)                                                       
      H2(I,J)=H2(J,I)                                                   
  500 H2(J,I)=TEM                                                       
      N2M1=N2-1                                                         
      DO 550 J=1,N2M1                                                   
      TEM=DABS(H2(J,J))                                                 
  550 ALNDTD=ALNDTD+DLOG(TEM)                                           
  200 CONTINUE                                                          
      ITEM=2                                                            
      IF(IP.EQ.1)ITEM=1                                                 
      CALL SETD(T0,ITEM,ID,WT,IART,ARFT)                                
      IDAR = ID + IART                                                  
      DO 10 I=1,IDAR                                                    
   10 TI(I) = FTRN(I)*ALPHA                                             
      CALL INIT(TI,IDAR,T0,ITEM)                                        
      ICOUNT=0                                                          
      M1=0                                                              
      ID0=IDAR*ITEM+1                                                   
      IPOS=0                                                            
      K=0                                                               
      DO 100 I=1,N                                                      
      IPOS = IPOS + 1                                                   
      N3 = IPOS                                                         
      IF(N3 .GT. ID0) N3 = ID0                                          
      JTEM = ID0 - N3                                                   
      DO 110 J=1,N3                                                     
      JTEM = JTEM + 1                                                   
  110 H3(J) = T0(JTEM)                                                  
      DO 120 J=1,N2                                                     
  120 H4(J) = 0.D0                                                      
      IF(I .GT. IDAR) GO TO 125                                         
      H4(N2)=TI(I)                                                      
      DO 121 J=1,N3                                                     
  121 H3(J)=H3(J)*ALPHA                                                 
  125 CONTINUE                                                          
      CALL EXHSLD(H1,N1,H2,N2,H3,N3,H4,M1,IPOS)                         
      N3 = -1                                                           
      IF(IP .GT. 1) IPOS=IPOS+1                                         
      IF(I .GT. NDATA) GO TO 90                                         
      IF(Y(I) .GT. RLIM .AND. RLIM .GT. 0.D0) GO TO 90                  
      ICOUNT=ICOUNT+1                                                   
      N3 = ITEM                                                         
      H3(1) = 1.D0                                                      
      IF(IP.NE.1)H3(2) = 1.D0                                           
      DO 30 J=1,N2                                                      
      IF(J .LT. N2) H4(J) = WEEK(J,I)                                   
      IF(J .EQ. N2) H4(J) = Y(I)                                        
   30 CONTINUE                                                          
   90 CALL EXHSLD(H1,N1,H2,N2,H3,N3,H4,M1,IPOS)                         
      IF(IP .EQ. 1) GO TO 100                                           
      IF(IPOS .LE. N1) GO TO 100                                        
      K=K+1                                                             
      J0=0                                                              
      DO 91 J=1,N1                                                      
      IF(MOD(J,2) .EQ. 0) GO TO 911                                     
      J0=J0+1                                                           
      H3(J)=H(J0,K)*DD                                                  
      GO TO 91                                                          
  911 H3(J)=0.D0                                                        
   91 CONTINUE                                                          
      DO 92 J=1,N2                                                      
   92 H4(J)=0.D0                                                        
      H4(N2)=F(K)*DD                                                    
      CALL EXHSLD(H1,N1,H2,N2,H3,N1,H4,M1,IPOS)                         
  100 CONTINUE                                                          
      IF(IP .EQ. 1) RETURN                                              
      NMK=N1                                                            
      KP1=K+1                                                           
      DO 95 K=KP1,N                                                     
      NMK=NMK-2                                                         
      J0=0                                                              
      DO 96 J=1,NMK                                                     
      IF(MOD(J,2).EQ.0)GO TO 93                                         
      J0=J0+1                                                           
      H3(J)=H(J0,K)*DD                                                  
      GO TO 96                                                          
   93 H3(J)=0.D0                                                        
   96 CONTINUE                                                          
      DO 97 L=1,N2                                                      
   97 H4(L)=0.D0                                                        
      H4(N2)=F(K)*DD                                                    
   95 CALL EXHSLD(H1,N1,H2,N2,H3,NMK,H4,M1,M1)                          
      RETURN                                                            
      END                                                               
C                                                                       
C                                                                       
C                                                                       
cc      SUBROUTINE SETDC(H1,N1,H2,M1,FSEAS,N,WS,WZ,IARS,ARFS,IARN,ARFN)   
cc     *                 ARFN,ALPHA,BETA,GAMMA,WTRD,DELTA,IP,ID,IS,YEAR)
      SUBROUTINE SETDC(H1,N1,H2,M1,FSEAS,N,WS,WZ,IARS,ARFS,IARN,
     *                  ARFN,ALPHA,BETA,GAMMA,WTRD,DELTA,IP,IS,YEAR)
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      INTEGER*4 YEAR                                                    
cc      DIMENSION  H1(N1,1),H2(1),FSEAS(1),H3(100),S0(100),ARFS(1),ARFN(1)
cc      DIMENSION  Z0(50),SI(100),ZI(50),H4(50)                           
      DIMENSION  H1(N1,1),H2(1),FSEAS(1),H3((IS+IARS)*IP+IARN+1),
     *            S0((IS+IARS+1)*IP+IARN+1),ARFS(1),ARFN(1)
      DIMENSION  Z0(IP),SI((IS+IARS)*IP+IARN),ZI(IP),H4(1)
cc      COMMON /RDATA/ ALPHA,BETA,GAMMA,DUMMY(4),WTRD,DELTA               
cc      COMMON /IDATA/ IP,ID,IS,YEAR                                      
      IPIS = (IS+IARS)*IP + IARN                                        
      IPM1 =IP - 1                                                      
      ITEM = IP*(IS-1) + 1                                              
      DO 30 I=1,IPM1                                                    
      ITEM = ITEM + 1                                                   
   30 ZI(I) = FSEAS(ITEM)*WZ*GAMMA                                      
      SUM = 0.D0                                                        
      ITEM = IP                                                         
      DO 40 I=1,IPM1                                                    
      ITEM = ITEM - 1                                                   
      SUM = SUM - ZI(ITEM)                                              
   40 ZI(ITEM) = SUM                                                    
      CALL SETD(S0,IP,IS,WS,IARS,ARFS)                                  
      IF(IARN .EQ. 0) GO TO 49                                          
      LENGTH = IPIS + 1                                                 
      LTEM = LENGTH - IARN                                              
      DO 41 I=1,LTEM                                                    
      S0(LENGTH)=S0(LENGTH - IARN)                                      
   41 LENGTH = LENGTH - 1                                               
      DO 42 I=1,IARN                                                    
   42 S0(I)=0.D0                                                        
      LENGTH = IPIS + 1                                                 
      DO 400 I=1,LENGTH                                                 
      DT=S0(I)                                                          
      DO 410 J=1,IARN                                                   
  410 IF(I+J .LE. LENGTH) DT=DT-ARFN(J)*S0(I+J)                         
  400 S0(I)=DT                                                          
   49 CONTINUE                                                          
      IF(IPIS .EQ. 0) GO TO 998                                         
      DO 20 I=1,IPIS                                                    
   20 SI(I) = FSEAS(I)*BETA                                             
      CALL INIT(SI,IPIS,S0,1)                                           
  998 DO 50 I=1,IP                                                      
      Z0(I)=WZ                                                          
   50 CONTINUE                                                          
   55 CONTINUE                                                          
      IS0 = IPIS + 1                                                    
      IZ0 = IP                                                          
      M1 = 0                                                            
      IPOS=0                                                            
      H2(1) = 0.D0                                                      
      DO 100 I=1,N                                                      
      IPOS = I                                                          
      N3 = IPOS                                                         
      IF(N3 .GT. IS0) N3 = IS0                                          
      ITEM = IS0 - N3                                                   
      DO 130 J=1,N3                                                     
      ITEM = ITEM + 1                                                   
  130 H3(J) = S0(ITEM)                                                  
      H4(1) = 0.D0                                                      
      IF(I .GT. IPIS) GO TO 145                                         
      H4(1)=SI(I)                                                       
      DO 141 J=1,N3                                                     
  141 H3(J)=H3(J)*BETA                                                  
  145 CONTINUE                                                          
      CALL EXHSLD(H1,N1,H2,1,H3,N3,H4,M1,IPOS)                          
      N3 = IPOS                                                         
      IF(N3 .GT. IZ0) N3 =IZ0                                           
      ITEM = IZ0 - N3                                                   
      DO 150 J=1,N3                                                     
      ITEM = ITEM + 1                                                   
  150 H3(J) = Z0(ITEM)                                                  
      H4(1) = 0.D0                                                      
      IF(I .GE. IP) GO TO 165                                           
      H4(1)=ZI(I)                                                       
      DO 161 J=1,N3                                                     
  161 H3(J)=H3(J)*GAMMA                                                 
  165 CONTINUE                                                          
      CALL EXHSLD(H1,N1,H2,1,H3,N3,H4,M1,IPOS)                          
  100 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE PARTAR(R,A,M)                                          
C **** PARCOR R TO AR ********                                          
      IMPLICIT REAL*8(A-H,O-Z)                                          
      DIMENSION R(M),A(M,M)                                             
      DO 10 I=1,M                                                       
      DO 10 J=1,I                                                       
   10 A(I,J)=0.D0                                                       
      A(1,1)=R(1)                                                       
      IF(M.LE.1) RETURN                                                 
      DO 40 I=2,M                                                       
      A(I,I)=R(I)                                                       
      IM1=I-1                                                           
      DO 30 J=1,IM1                                                     
   30 A(I,J)=A(I-1,J)-R(I)*A(I-1,I-J)                                   
   40 CONTINUE                                                          
      RETURN                                                            
      END                                                               
cc      SUBROUTINE SOLVE(H1,N1,H2,N2,A,M1,SQE,NANS,ERR)                   
      SUBROUTINE BSOLVE(H1,N1,H2,N2,A,M1,SQE,NANS,ERR)
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      DIMENSION  H1(N1,1), H2(N2,1), A(1), ERR(1)                       
C                                                                       
      DO 30 I=1,NANS                                                    
   30 ERR(I)=0.D0                                                       
      DO 100 LER = 1, NANS                                              
         K = M1 + N2                                                    
         KA = NANS                                                      
         JJ = NANS - 1                                                  
         IF( LER.EQ.NANS )  GO TO 44                                    
            SQE = 0.D0                                                  
            KAM1 = KA - 1                                               
            DO 40 KKA = 1, KAM1                                         
   40       A(KKA) = 0.D0                                               
            A(LER) = 1.D0                                               
            GO TO 48                                                    
   44    CONTINUE                                                       
            SQE = H2(N2,K)**2                                           
            KKA = KA                                                    
            KK  = K                                                     
            DO 46 J = 1, JJ                                             
               KK = KK - 1                                              
               KKA = KKA - 1                                            
               A(KKA) = H2(N2,KK)                                       
   46       CONTINUE                                                    
   48    CONTINUE                                                       
C                                                                       
         KCOPY = K                                                      
         DO 50 J = 1, JJ                                                
            KA = KA - 1                                                 
            IF( A(KA).EQ.0.D0 )  GO TO 50                               
            K = KCOPY - J                                               
            IF( J.GE.N2 )  GO TO 20                                     
C                                                                       
            A(KA) = A(KA)/H2(N2-J,K)                                    
            IF( LER.LT.NANS )  ERR(KA) = ERR(KA) + A(KA)**2             
            KM1 = KA - 1                                                
            IF( KM1.LE.0 )  GO TO 50                                    
            LTEM = K                                                    
            KATEM = KA                                                  
            AKA = A(KA)                                                 
            N2MJ = N2-J                                                 
            DO 10 L = 1, KM1                                            
               KATEM = KATEM - 1                                        
               LTEM = LTEM - 1                                          
               A(KATEM) = A(KATEM) - AKA*H2(N2MJ,LTEM)                  
   10       CONTINUE                                                    
            GO TO 50                                                    
C                                                                       
   20      A(KA) = A(KA)/H1(1,K)                                        
            IF( LER.LT.NANS )  ERR(KA) = ERR(KA) + A(KA)**2             
            LTEM = K                                                    
            L = KA                                                      
            IF( N1.LT.2 )  GO TO 50                                     
            DO 25 I = 2, N1                                             
               L = L - 1                                                
               LTEM = LTEM - 1                                          
               IF( L.LE.0 )  GO TO 50                                   
               A(L) = A(L) - A(KA)*H1(I,LTEM)                           
   25       CONTINUE                                                    
   50    CONTINUE                                                       
  100 CONTINUE                                                          
C                                                                       
      RETURN                                                            
      END                                                               
