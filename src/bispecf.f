      SUBROUTINE BISPECF(N,MH,CC0,C0,P1,P2,Q,A,BR,BI,RAT)
C
      INCLUDE 'timsac.h'
C
cc	PROGRAM BISPEC
C     PROGRAM 74.6.2.
C-----------------------------------------------------------------------
C     ** DESIGNED BY H. AKAIKE, THE INSTITUTE OF STATISTICAL MATHEMATICS
C     ** PROGRAMMED BY E. ARAHATA, THE INSTITUTE OF STATISTICAL MATHEMAT
C         TOKYO
C     ** DATE OF THE LATEST REVISION: MARCH 25, 1977
C     ** THIS PROGRAM WAS ORIGINALLY PUBLISHED IN
C         "TIMSAC-74 A TIME SERIES ANALYSIS AND CONTROL PROGRAM PACKAGE(2
C         BY H. AKAIKE, E. ARAHATA AND T. OZAKI, COMPUTER SCIENCE MONOGRA
C         NO.6 MARCH 1976, THE INSTITUTE OF STATISTICAL MATHEMATICS
C-----------------------------------------------------------------------
C     THIS PROGRAM COMPUTES BISPECTRUM USING THE DIRECT FOURIER TRANSFOR
C     OF SAMPLE THIRD ORDER MOMENTS.
C     THIS PROGRAM REQUIRES THE FOLLOWING INPUTS;
C     OUTPUTS OF THE PROGRAM THIRMO:
C      N; DATA LENGTH,
C      MH; MAXIMUM LAG,
C      CC(I); AUTOCOVARIANCES,
C      C(I,J); THIRD ORDER MOMENTS.
C
cc      !DEC$ ATTRIBUTES DLLEXPORT :: BISPECF
C
cxx      IMPLICIT REAL*8(A-H,O-Z)
cc      DIMENSION C(51,51),C1(102),S1(102)
cc      DIMENSION CL(51,51),SL(51,51)
cc      DIMENSION CA1(51,51),CA2(51,51)
cc      DIMENSION P(51)
cc      DIMENSION CC(51)
cxx      DIMENSION C(MH+1,MH+1),C1((MH+1)*2),S1((MH+1)*2)
cxx      DIMENSION CL(MH+1,MH+1),SL(MH+1,MH+1)
cxx      DIMENSION CA1(MH+1,MH+1),CA2(MH+1,MH+1)
cxx      DIMENSION P1(MH+1),P2(MH+1),Q(MH+1)
cxx      DIMENSION CC(MH+1)
cxx      DIMENSION CC0(MH+1),C0(MH+1,MH+1)
cxx      DIMENSION A(MH+1,MH+1),BR(MH+1,MH+1),BI(MH+1,MH+1)
      INTEGER N, MH
      DOUBLE PRECISION CC0(MH+1), C0(MH+1,MH+1), P1(MH+1), P2(MH+1),
     1                 Q(MH+1), A(MH+1,MH+1), BR(MH+1,MH+1),
     2                 BI(MH+1,MH+1), RAT
c local
      DOUBLE PRECISION C(MH+1,MH+1), C1((MH+1)*2), S1((MH+1)*2),
     1                 CL(MH+1,MH+1), SL(MH+1,MH+1), CA1(MH+1,MH+1),
     2                 CA2(MH+1,MH+1), CC(MH+1), CST0, CST1, CST2,
     3                 CST6, CST2I, CST6I, H, PI, PP, AI, T, C00, SL0,
     4                 TC1, TC2, TC3, TS1, TS3, CN0, SN0
C
C     INPUT / OUTPUT DATA FILE OPEN
cc	CALL SETWND
cc	CALL FLOPN2(NFL)
cc	IF (NFL.EQ.0) GO TO 999
C
      MH1 = MH+1
      A(1:MH1,1:MH1) = 0.0D-00
      BR(1:MH1,1:MH1) = 0.0D-00
      BI(1:MH1,1:MH1) = 0.0D-00
      CL(1:MH1,1:MH1) = 0.0D-00
      SL(1:MH1,1:MH1) = 0.0D-00
c
cc      MJ1=51
      CST0 =0.0D-00
      CST1 =1.0D-00
      CST2 =2.0D-00
      CST6=6.0D-00
      CST2I=1.0D-00/CST2
      CST6I=1.0D-00/CST6
C     INITIAL LOADING
C     AUTOCOVARIANCE AND THIRD ORDER MOMENT INPUT
cc	READ(5,1) N,MH
      L1=MH+1
cc	READ(5,2) (CC(I),I=1,L1)
      H=MH
cxx      DO 10 I=1,L1
      DO 11 I=1,L1
cc   10 READ(5,2) (C(I,J),J=1,I)
         CC(I)=CC0(I)
         DO 10 J=1,I
            C(I,J)=C0(I,J)
   10 CONTINUE
   11 CONTINUE
cc      WRITE(6,130)
cc      WRITE(6,129) N,MH
C     AUTOCOVARIANCES PRINT OUT
cc      WRITE(6,163)
cc      DO 164 I=1,L1
cc      IM1=I-1
cc  164 WRITE(6,165) IM1,CC(I)
C     POWER SPECTRUM COMPUTATION
cc      CALL SAUSP1(CC,P,N,L1,L1)
      CALL SAUSP1(CC,P1,P2,Q,N,L1,L1)
C     THIRD ORDER MOMENTS PRINT OUT
cc      WRITE(6,180)
cc      CALL MOMTPR (C,MJ1,L1)
C     END CORRECTIONS
      C(1,1)=C(1,1)*CST6I
      DO 112 I=2,L1
      C(I,1) =C(I,1)*CST2I
      C(I,I)=C(I,I)*CST2I
cxx  112 C(L1,I)=C(L1,I)*CST2I
      C(L1,I)=C(L1,I)*CST2I
  112 CONTINUE
      C(L1,1) =C(L1,1)*CST2I
C     FOURIER TRANSFORMATION
      PI=3.141592653
      PP=PI/H
      I2H=MH+MH
      DO 115 I=1,I2H
      AI=I
      T=PP*AI
      C1(I)=DCOS(T)
cxx  115 S1(I)=DSIN(T)
      S1(I)=DSIN(T)
  115 CONTINUE
C     CL(0,0) COMPUTATION
      C00 =CST6*C(1,1)
      SL0 =CST0
      IS=0
      ID=0
      IR=0
      TC1 =CST0
      DO 600 MMS1=1,L1
      TC2=CST0
      DO 610 NS1=1,MMS1
cxx  610 TC2=TC2+C(MMS1,NS1)
      TC2=TC2+C(MMS1,NS1)
  610 CONTINUE
      TC1=TC1+TC2
  600 CONTINUE
      CL(1,1)=CST6*TC1
C     CL(IS,0) COMPUTATION
      DO 120 IS=1,MH
      ID=IS
      TC3=CST0
      MS=0
      DO 122 MMS=1,MH
      MS=MS+IS
      IF(MS.LE.I2H) GO TO 461
      MS=MS-I2H
  461 MC1=MS
      MC2=MS
      MC3=0
      TC2 =C(MMS+1,1)*(C1(MC1)+C1(MC2)+CST1)
      TC1 =CST0
      DO 123 NS=1,MMS
      MC2=MC2-IS
      MC3=MC3+IS
      IF(MC2.GT.0) GO TO 561
      MC2=MC2+I2H
  561 IF(MC3.LE.I2H) GO TO 562
      MC3=MC3-I2H
  562 TC1=TC1+C(MMS+1,NS+1)*(C1(MC1)+C1(MC2)+C1(MC3))
  123 CONTINUE
      TC3=TC3+TC2+TC1
  122 CONTINUE
      IS1=IS+1
      CL(IS1,1)=C00+TC3+TC3
      CL(1,IS1)=CL(IS1,1)
      SL(IS1,1) =CST0
      SL(1,IS1)=SL(IS1,1)
  120 CONTINUE
C     CL(ID,IR),SL(ID,IR) COMPUTATION
      DO 20 IS=2,MH
      IRL=IS/2
      DO 21 IR=1,IRL
      TC3=CST0
      TS3 =CST0
      ID=IS-IR
      MD=0
      MS=0
      MR=0
      DO 22 MMS=1,MH
      MD=MD+ID
      MS=MS-IS
      MR=MR+IR
      IF(MD.LE.I2H) GO TO 401
      MD=MD-I2H
  401 IF(MS.GT.0) GO TO 412
      MS=MS+I2H
  412 IF(MR.LE.I2H) GO TO 413
      MR=MR-I2H
  413 MDR=MD
      MDS=MDR
      MSR=MS
      MSD=MSR
      MRD=MR
      MRS=MRD
      CN0=C(MMS+1,1)*(C1(MDR)+C1(MDS)+C1(MSR)+C1(MSD)+C1(MRD)+C1(MRS))
      SN0=C(MMS+1,1)*(S1(MDR)+S1(MDS)+S1(MSR)+S1(MSD)+S1(MRD)+S1(MRS))
      TC1 =CST0
      TS1=CST0
      DO 23 NS=1,MMS
      MSD=MSD+ID
      MRD=MRD+ID
      MDR=MDR+IR
      MSR=MSR+IR
      MDS=MDS-IS
      MRS=MRS-IS
      IF(MSD.LE.I2H) GO TO 501
      MSD=MSD-I2H
  501 IF(MRD.LE.I2H) GO TO 511
      MRD=MRD-I2H
  511 IF(MDR.LE.I2H) GO TO 521
      MDR=MDR-I2H
  521 IF(MSR.LE.I2H) GO TO 531
      MSR=MSR-I2H
  531 IF(MDS.GT.0) GO TO 542
      MDS=MDS+I2H
  542 IF(MRS.GT.0) GO TO 552
      MRS=MRS+I2H
  552 TC1=TC1+C(MMS+1,NS+1)*(C1(MSD)+C1(MRD)+C1(MDR)+C1(MSR)+C1(MDS)+C1(
     AMRS))
      TS1=TS1+C(MMS+1,NS+1)*(S1(MSD)+S1(MRD)+S1(MDR)+S1(MSR)+S1(MDS)+S1(
     AMRS))
   23 CONTINUE
      TC3=TC3+CN0+TC1
      TS3=TS3+SN0+TS1
   22 CONTINUE
      IR1=IR+1
      ID1=ID+1
      CL(ID1,IR1)=C00+TC3
      SL(ID1,IR1)=TS3
      CL(IR1,ID1)=CL(ID1,IR1)
      SL(IR1,ID1)=SL(ID1,IR1)
   21 CONTINUE
   20 CONTINUE
C     SMOOTHING OF THE RAW ESTIMATES
C     ROWWISE SMOOTHING
      CALL SUBCA(CL,CA1,MH,0)
      CALL SUBCA(SL,CA2,MH,1)
C     COLUMNWISE SMOOTHING
      CALL SUBCB(CA1,CL,MH)
      CALL SUBCB(CA2,SL,MH)
C     OBLIQUE SMOOTHING
cc      WRITE(6,160)
cc      CALL SUBCD(CL,CA1,MH)
cc      WRITE(6,161)
cc      CALL SUBCD(SL,CA2,MH)
C     Q1(I,J) COMPUTATION
cc      CALL SUBQ1(CA1,CA2,P,N,MH)
      CALL SUBCD(CL,CA1,MH,BR)
      CALL SUBCD(SL,CA2,MH,BI)
      CALL SUBQ1(CA1,CA2,P1,N,MH,A,RAT)
cc      CALL FLCLS2(NFL)
cc  999 CONTINUE
      RETURN
cxx    1 FORMAT(16I5)
cxx    2 FORMAT(4D20.10)
cxx  130 FORMAT(1H ,'PROGRAM 74.6.2. BISPEC / BISPECTRUM COMPUTATION.')
cxx  180 FORMAT(//1H ,4X,'C(I,J): THIRD ORDER MOMENTS'/)
cxx  129 FORMAT(1H ,'N=',I5,5X,'MH=',I5)
cxx  160 FORMAT(/1H ,'REAL PART OF BISPECTRUM B(I,J)'/)
cxx  161 FORMAT(/1H ,'IMAGINARY PART OF BISPECTRUM B(I,J)'/)
cxx  163 FORMAT(1H ,'CC(I):',1X,'AUTOCOVARIANCES')
cxx  165 FORMAT(1H ,I5,D16.5)
      END
C
      SUBROUTINE SUBCA(CL,CA,MH,ISW)
C     ROWWISE SMOOTHING
C     CA(I,J)=
C              (CL(I-1,J)+2.0*CL(I,J)+CL(I+1,J))/4.0
C     ISW=0: COS
C     ISW=1: SIN
cxx      IMPLICIT REAL*8(A-H,O-Z)
cc      DIMENSION CL(51,51),CA(51,51)
cxx      DIMENSION CL(MH+1,MH+1),CA(MH+1,MH+1)
      INTEGER MH, ISW
      DOUBLE PRECISION CL(MH+1,MH+1), CA(MH+1,MH+1)
c local
      DOUBLE PRECISION CST2, CST4, CST4I
      CST2=2.0D-00
      CST4=4.0D-00
      CST4I=1.0D-00/CST4
      L1=MH+1
C     ON AND ABOVE X-AXIS COMPUTATION
      JJL=MH/2+1
      DO 10 JJ=1,JJL
      J=JJ-1
      IF(J.LE.1) GO TO 11
      IJ=J
      GO TO 12
   11 IJ=2
   12 IL=MH-J
      DO 20 II=IJ,IL
      CA(II,JJ)=
     A        (CL(II-1,JJ)+CST2*CL(II,JJ)+CL(II+1,JJ))*CST4I
      I=II-1
   20 CONTINUE
   10 CONTINUE
C     BELOW X-AXIS
      DO 40 J1=1,2
      JJ=MH/2+J1+1
      J2=J1+2
      DO 41 II=J2,MH
      I3=II-J1
      IF(ISW.EQ.1) GO TO 42
      CA(II,JJ)=CA(I3,J1+1)
      I=II-1
      J=-J1
      GO TO 41
   42 CA(II,JJ)=-CA(I3,J1+1)
      I=II-1
      J=-J1
   41 CONTINUE
   40 CONTINUE
      RETURN
      END
C
      SUBROUTINE SUBCB(CA,CB,MH)
C     COLUMNWISE SMOOTHING
C     CB(I,J)=
C             (CA(I,J-1)+2.0*CA(I,J)+CA(I,J+1))/4.0
cxx      IMPLICIT REAL*8(A-H,O-Z)
cc      DIMENSION CA(51,51),CB(51,51)
cxx      DIMENSION CA(MH+1,MH+1),CB(MH+1,MH+1)
      INTEGER MH
      DOUBLE PRECISION CA(MH+1,MH+1), CB(MH+1,MH+1)
c local
      DOUBLE PRECISION CST2, CST4, CST4I
C     ON AND ABOVE 1-AXIS
      CST2=2.0D-00
      CST4=4.0D-00
      CST4I=1.0D-00/CST4
      L1=MH+1
      MHJ=MH/2
      DO 10 J=2,MHJ
      MHI=MH-J
      DO 11 I=J,MHI
      CB(I,J)=
     A      (CA(I,J-1)+CST2*CA(I,J)+CA(I,J+1))*CST4I
      IM1=I-1
      JM1=J-1
   11 CONTINUE
   10 CONTINUE
C     ON 0-AXIS
      MH1=MH-1
      JJ=MH/2+2
      DO 12 I=3,MH1
      CB(I,1)=
     A      (CA(I,JJ)+CST2*CA(I,1)+CA(I,2))*CST4I
      IM1=I-1
      JC=0
   12 CONTINUE
C     ON (-1)-AXIS
      JJ1=MHJ+1
      DO 14 I=4,MH
      CB(I,JJ1)=
     A        (CA(I,JJ+1)+CST2*CA(I,JJ)+CA(I,1))*CST4I
      IM1=I-1
      JC=-1
   14 CONTINUE
      RETURN
      END
C
cc      SUBROUTINE SUBCD(CB,CD,MH)
      SUBROUTINE SUBCD(CB,CD,MH,B)
C     OBLIQUE SMOOTHING
C     CD(I,J)=
C             (CB(I-1,J-1)+2.0*CB(I,J)+CB(I+1,J+1))/4.0
cxx      IMPLICIT REAL*8(A-H,O-Z)
cc      INTEGER KANA1 / ' B  ' /
cc      DIMENSION CB(51,51),CD(51,51)
cxx      DIMENSION CB(MH+1,MH+1),CD(MH+1,MH+1)
cxx      DIMENSION B(0:MH,0:MH)
      INTEGER MH
      DOUBLE PRECISION CB(MH+1,MH+1), CD(MH+1,MH+1), B(0:MH,0:MH)
c local
      DOUBLE PRECISION CST2, CST4, CST4I, CZ
      CST2=2.0D-00
      CST4=4.0D-00
      CST4I=1.0D-00/CST4
      CZ=0.0D0
cc      CALL BISPPR (KANA1,0,CZ,0)
      L1=MH+1
C     ON AND ABOVE 2-AXIS
      JL=MH/2-1
      DO 10 J=3,JL
      IL=MH-2-J
      DO 11 I=J,IL
      CD(I,J)=
     A      (CB(I-1,J-1)+CST2*CB(I,J)+CB(I+1,J+1))*CST4I
      IM1=I-1
      JM1=J-1
      CZ=CD(I,J)
cc      CALL BISPPR (IM1,JM1,CZ,1)
      B(IM1,JM1)=CZ
   11 CONTINUE
   10 CONTINUE
C     ON 1-AXIS
      MH4=MH-4
      DO 12 I=4,MH4
      CD(I,2)=
     A      (CB(I-1,1)+CST2*CB(I,2)+CB(I+1,3))*CST4I
      IM1=I-1
      JC=1
      CZ=CD(I,2)
cc      CALL BISPPR(IM1,JC,CZ,1)
      B(IM1,JC)=CZ
   12 CONTINUE
C     ON 0-AXIS
      MH3=MH-3
      JJ1=MH/2+1
      DO 14 I=5,MH3
      CD(I,1)=
     A      (CB(I-1,JJ1)+CST2*CB(I,1)+CB(I+1,2))*CST4I
      IM1=I-1
      JC=0
      CZ=CD(I,1)
cc      CALL BISPPR (IM1,JC,CZ,1)
      B(IM1,JC)=CZ
   14 CONTINUE
      RETURN
      END
C
cc      SUBROUTINE SUBQ1(CL,SL,P,N,MH)
      SUBROUTINE SUBQ1(CL,SL,P,N,MH,A,RAT)
C     A MEASURE OF COHERENCE BETWEEN X(I)X(J) AND X(I+J) IS GIVEN BY Q1(
C     WHICH IS DEFINED BY
C     Q1(I,J)=(CL(I,J)**2+SL(I,J)**2)/(P(I)P(J)P(I+J)MH).
cxx      IMPLICIT REAL*8(A-H,O-Z)
cc      INTEGER KANA2 / 'Q1  ' /
cc      DIMENSION CL(51,51),SL(51,51),P(51)
cxx      DIMENSION CL(MH+1,MH+1),SL(MH+1,MH+1),P(MH+1)
cxx      DIMENSION A(0:MH,0:MH)
      INTEGER N, MH
      DOUBLE PRECISION CL(MH+1,MH+1), SL(MH+1,MH+1), P(MH+1),
     1                 A(0:MH,0:MH), RAT
c local
      DOUBLE PRECISION H, CZ, AN, AS3, CST075
c
      L1=MH+1
      H = MH
      CZ =0.0D-00
cc      WRITE(6,185)
cc      CALL BISPPR (KANA2,0,CZ,0)
C     ON AND ABOVE 2-AXIS
      JL=MH/2-1
cxx      DO 10 J=3,JL
      DO 11 J=3,JL
      IL=MH-2-J
      DO 10 I=J,IL
      I4=I+J-1
      CL(I,J)=(CL(I,J)*CL(I,J)+SL(I,J)*SL(I,J))/P(I)/P(J)/P(I4)/H
      IM1=I-1
      JM1=J-1
      CZ=CL(I,J)
cc      CALL BISPPR (IM1,JM1,CZ,1)
      A(IM1,JM1)=CZ
   10 CONTINUE
   11 CONTINUE
C     ON 1-AXIS
      MH4=MH-4
      DO 12 I=4,MH4
      I4=I+1
      CL(I,2)=(CL(I,2)*CL(I,2)+SL(I,2)*SL(I,2))/P(I)/P(2)/P(I4)/H
      IM1=I-1
      JC=1
      CZ=CL(I,2)
cc      CALL BISPPR (IM1,JC,CZ,1)
      A(IM1,JC)=CZ
   12 CONTINUE
C     ON 0-AXIS
      MH3=MH-3
      JJ1=MH/2+1
      DO 14 I=5,MH3
      IM1=I-1
      CL(I,1)=(CL(I,1)*CL(I,1)+SL(I,1)*SL(I,1))/P(I)/P(I)/P(1)/H
      JC=0
      CZ=CL(I,1)
cc      CALL BISPPR (IM1,JC,CZ,1)
      A(IM1,JC)=CZ
   14 CONTINUE
      AN=N
      AS3 =3.0D-00
      CST075=0.75D-00
      RAT=(H/AN)*CST075*CST075/DSQRT(AS3)
cc      WRITE(6,186) RAT
      RETURN
cxx  181 FORMAT(1H ,2I5,D17.5)
cxx  185 FORMAT(/1H ,'Q1(I,J)=(CL(I,J)**2+SL(I,J)**2)/(P(I)P(J)P(I+J)MH);',
cxx     1'A MESURE OF COHERENCE BETWEEN X(I)*X(J) AND X(I+J)'/)
cxx  186 FORMAT(/1H ,'APPROXIMATE EXPECTED VALUE OF Q1(I,J) UNDER ',
cxx     A'GAUSSIAN ASSUMPTION;',D17.5,'.')
      END
C
cc      SUBROUTINE SAUSP1(CXX,P1,N,LAGH3,LAGH1)
      SUBROUTINE SAUSP1(CXX,P1,P2,Q,N,LAGH3,LAGH1)
C     THIS SUBROUTINE COMPUTES POWER SPECTRUM.
cxx      IMPLICIT REAL*8(A-H,O-Z)
cc      DIMENSION CXX(51),FC(51),P1(51),P2(51),Q(51),P3(51)
cc      DIMENSION A1(10),A2(10)
cxx      DIMENSION CXX(LAGH1),FC(LAGH1),P1(LAGH1),P2(LAGH1),Q(LAGH1)
cxx      DIMENSION A1(2),A2(3)
      INTEGER N, LAGH3, LAGH1
      DOUBLE PRECISION CXX(LAGH1), P1(LAGH1), P2(LAGH1), Q(LAGH1)
c local
      DOUBLE PRECISION FC(LAGH1), A1(2), A2(3)
C     WINDOW W1 DEFINITION
      MLA1=2
      A1(1)=0.5D-00
      A1(2) =0.25D-00
C     WINDOW W2 DEFINITION
      MLA2=3
      A2(1)=0.625D-00
      A2(2)=0.25D-00
      A2(3)=-0.0625D-00
      LAGH=LAGH1-1
      LAGH0=LAGH3-1
      DO 10 I=2,LAGH
cxx   10 CXX(I)=CXX(I)+CXX(I)
      CXX(I)=CXX(I)+CXX(I)
   10 CONTINUE
C     F-COS TRANSFORMATION
C     COMMON SUBROUTINE CALL
      CALL FGERCO(CXX,LAGH1,FC,LAGH1)
C     SPECTRUM SMOOTHING BY WINDOW W1
C     COMMON SUBROUTINE CALL
      CALL AUSP(FC,P1,LAGH1,A1,MLA1)
C     SPECTRUM SMOOTHING BY WINDOW W2
C     COMMON SUBROUTINE CALL
      CALL AUSP(FC,P2,LAGH1,A2,MLA2)
C     TEST STATISTICS COMPUTATION
C     COMMON SUBROUTINE CALL
      CALL SIGNIF(P1,P2,Q,LAGH1,N)
C     AUTO SPECTRUM AND TEST STATISTICS PRINT OUT
cc      WRITE(6,66) N,LAGH
cc      WRITE(6,63) A1(1),A1(2)
cc      WRITE(6,64) A2(1),A2(2),A2(3)
cc      WRITE(6,67)
C     COMMON SUBROUTINE CALL
cc      CALL PRCOL3(P1,P2,Q,1,LAGH1,1)
      RETURN
cxx   63 FORMAT(1H ,'WINDOW W1',5X,'A(0)=',F10.5,5X,'A(1)=',F10.5)
cxx   64 FORMAT(1H ,'WINDOW W2',5X,'A(0)=',F10.5,5X,'A(1)=',F10.5,5X,'A(2)=
cxx     1',F10.5)
cxx   66 FORMAT(//1H ,14HPOWER SPECTRUM,5X,2HN=,I5,5X,5HLAGH=,I5)
cxx   67 FORMAT(/1H ,4X,1HI,8X,8HPOWER W1,6X,8HPOWER W2,2X,12HSIGNIFICANCE)
      END
