C-----------------------------------------------------------------------
C     This  USDFLD is used for considering the combined hardening behavior 
C     of low-carbon steel!
C                 written by Wang Youde                   2024-05-19
C-----------------------------------------------------------------------
      SUBROUTINE USDFLD(FIELD,STATEV,PNEWDT,DIRECT,T,CELENT,TIME,DTIME,
     1 CMNAME,ORNAME,NFIELD,NSTATV,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,
     2 NDI,NSHR,COORD,JMAC,JMATYP,MATLAYO,LACCFLA)
C
      INCLUDE 'ABA_PARAM.INC'
C PLATEAU REGION STRAIN (EPT IS THE YIELD PLATEAU LENGTH, ST IS THE C IS THE MATERIAL PARAMETER, SET AS(0~0.5))
      PARAMETER(EPT=0.00834d0)
      PARAMETER(ST=0.5d0)
      PARAMETER(Qs=-129.0d0)
      PARAMETER(bs=600.d0)      
      PARAMETER(Qh=128.1d0)
      PARAMETER(bh=100.d0)
      PARAMETER(Q0=10.d0)
      PARAMETER(u=26.81d0)
C
C STRENGTH PARAMETERS
      PARAMETER(E=2.7183,TOL=0.00001,CLAMDA=0.178,ETAMONO=2.679)
      CHARACTER*80 CMNAME,ORNAME
      CHARACTER*3  FLGRAY(15)
      DIMENSION FIELD(NFIELD),STATEV(NSTATV),DIRECT(3,3),T(3,3),TIME(2)
      DIMENSION ARRAY(15),JARRAY(15),JMAC(*),JMATYP(*),COORD(*)
C
C****************************Constitutive Model**********************************************           
      CALL GETVRM('SDV',ARRAY,JARRAY,FLGRAY,JRCD,
     1    JMAC,JMATYP,MATLAYO,LACCFLA)
           Q11=ARRAY(1)
           Q22=ARRAY(2)
           Q33=ARRAY(3)
           Q12=ARRAY(4)
           Q23=ARRAY(5)
           Q31=ARRAY(6)
           P011=ARRAY(7)
           P022=ARRAY(8)
           P033=ARRAY(9)
           P012=ARRAY(10)
           P023=ARRAY(11)
           P031=ARRAY(12)
           PEEQ0=ARRAY(13)
           Tq=ARRAY(14)
           R=ARRAY(15)
           V=ARRAY(16)
           PE_TR=ARRAY(17)
           R_TR=ARRAY(18)
           V0=V
C           
      CALL GETVRM('PE',ARRAY,JARRAY,FLGRAY,JRCD,
     1    JMAC,JMATYP,MATLAYO,LACCFLA)
           P11=ARRAY(1)
           P22=ARRAY(2)
           P33=ARRAY(3)
           P12=ARRAY(4)
           P23=ARRAY(5)
           P31=ARRAY(6)
           PEEQ=ARRAY(7)
      IF ((PEEQ.EQ.0).OR.(PEEQ.EQ.PEEQ0)) THEN
           GOTO 10
      ELSE
C          
           DP11=(P11-P011)*(2.d0/3.d0)**0.5d0/(PEEQ-PEEQ0)
           DP22=(P22-P022)*(2.d0/3.d0)**0.5d0/(PEEQ-PEEQ0)
           DP33=(P33-P033)*(2.d0/3.d0)**0.5d0/(PEEQ-PEEQ0)
           DP12=(P12-P012)*(2.d0/3.d0)**0.5d0/(PEEQ-PEEQ0)
           DP23=(P23-P023)*(2.d0/3.d0)**0.5d0/(PEEQ-PEEQ0)
           DP31=(P31-P031)*(2.d0/3.d0)**0.5d0/(PEEQ-PEEQ0) 
C      
           PQ11=P11-Q11
           PQ22=P22-Q22
           PQ33=P33-Q33
           PQ12=P12-Q12
           PQ23=P23-Q23
           PQ31=P31-Q31
      Tq=SQRT(2.d0/3.d0*(PQ11*PQ11 + PQ22*PQ22 +PQ33*PQ33 
     1    + 2.d0*(PQ12*PQ12+PQ23*PQ23+PQ31*PQ31)))
           DQ11=PQ11*(2.d0/3.d0)**0.5d0/Tq
           DQ22=PQ22*(2.d0/3.d0)**0.5d0/Tq
           DQ33=PQ33*(2.d0/3.d0)**0.5d0/Tq
           DQ12=PQ12*(2.d0/3.d0)**0.5d0/Tq
           DQ23=PQ23*(2.d0/3.d0)**0.5d0/Tq
           DQ31=PQ31*(2.d0/3.d0)**0.5d0/Tq
C      
           g=Tq-R
      IF (g.LT.0) THEN
           Hg=0
      ELSE
           Hg=1
      ENDIF
           Tmn=DP11*DQ11 + DP22*DQ22 + DP33*DQ33
     1        + 2.d0*(DP12*DQ12+DP23*DQ23+DP31*DQ31)
      IF (Tmn.LT.0) THEN
           Tmn=0
      ENDIF      
C      
           R=R+0.5d0*Hg*Tmn*(PEEQ-PEEQ0)   
           Q11=Q11+(1.d0-0.5d0)*Hg*Tmn*(PEEQ-PEEQ0)*DQ11
           Q22=Q22+(1.d0-0.5d0)*Hg*Tmn*(PEEQ-PEEQ0)*DQ22
           Q33=Q33+(1.d0-0.5d0)*Hg*Tmn*(PEEQ-PEEQ0)*DQ33
           Q12=Q12+(1.d0-0.5d0)*Hg*Tmn*(PEEQ-PEEQ0)*DQ12
           Q23=Q23+(1.d0-0.5d0)*Hg*Tmn*(PEEQ-PEEQ0)*DQ23
           Q31=Q31+(1.d0-0.5d0)*Hg*Tmn*(PEEQ-PEEQ0)*DQ31     
C    
C UPDATE STATEV
      IF((R .GT. ST*EPT) .AND. (PEEQ .GT. EPT)) THEN
           V=1
      ENDIF
      IF((V0 .EQ. 0) .AND. (V .EQ. 1)) THEN
           PE_TR=PEEQ
           R_TR=R
      ENDIF
10       FIELD(1)=V    
      IF(V .EQ. 1)  THEN
           RR_TR=R-R_TR
           FIELD(2)=Qh+(Q0-Qh)*EXP(-u*RR_TR)
           FIELD(3)=bh*(PEEQ-PE_TR)/PEEQ
      ELSE
           RR_TR=0
           FIELD(2)=Qs
           FIELD(3)=bs
      ENDIF
C 
           STATEV(1)=Q11
           STATEV(2)=Q22
           STATEV(3)=Q33
           STATEV(4)=Q12
           STATEV(5)=Q23
           STATEV(6)=Q31
           STATEV(7)=P11
           STATEV(8)=P22
           STATEV(9)=P33
           STATEV(10)=P12
           STATEV(11)=P23
           STATEV(12)=P31
           STATEV(13)=PEEQ
           STATEV(14)=Tq      
           STATEV(15)=R
           STATEV(16)=V
           STATEV(17)=PE_TR
           STATEV(18)=R_TR   
C****************************CVGM**********************************************      
C GET TRIAXIALITY     
      CALL GETVRM('SINV',ARRAY,JARRAY,FLGRAY,JRCD,
     1     JMAC,JMATYP,MATLAYO,LACCFLA)
      IF (ARRAY(1).LE.TOL) THEN
            TRIAX=0
      ELSE
            TRIAX=-ARRAY(3)/ARRAY(1)
      ENDIF            
C GET STATEV OF LAST INC     
      CALL GETVRM('SDV',ARRAY,JARRAY,FLGRAY,JRCD,
     1     JMAC,JMATYP,MATLAYO,LACCFLA)
           PE_C=ARRAY(19)
           CVGI=ARRAY(20)
           RAT_MAX=ARRAY(23)
C CALCULATE COMPRESSIVE EQUIVALENT PLASTIC STRAIN            
      IF (TRIAX .LE. 0) THEN 
          PE_C=PE_C+PEEQ-PEEQ0!compress strain 
      ENDIF
C CALCULATE DEMAND PARAMETERS FOR CVGM 
      IF ((TRIAX .GT. 0) .AND. (TRIAX .LT. 10)) THEN
           CVGI=CVGI+(PEEQ-PEEQ0)*E**(ABS(1.5*TRIAX))
      ELSEIF ((TRIAX .LT. 0) .AND. (TRIAX .GT. -10))THEN
           CVGI=CVGI-0.85*(PEEQ-PEEQ0)*E**(ABS(1.5*TRIAX))
      ENDIF
      IF (CVGI.LT.0) THEN
            CVGI=0
      ENDIF
C CALCULATE CAPACITY PARAMETERS FOR CVGM 
      Yita_C=ETAMONO*E**(-CLAMDA*PE_C)
      IF (Yita_C .GT. TOL) THEN
           RAT=CVGI/Yita_C
      ENDIF
C      
      IF (RAT .GT. RAT_MAX)THEN
         RAT_MAX=RAT
      ENDIF
C ELEMENT DELETION     
      IF (RAT .GE. 1) THEN
           CRIT = 0
      ELSE 
           CRIT = 1
      ENDIF       
C UPDATE STATEV
      STATEV(19)=PE_C
      STATEV(20)=CVGI
      STATEV(21)=Yita_C
      STATEV(22)=RAT
      STATEV(23)=RAT_MAX
      STATEV(24)=CRIT     
C
      ENDIF   
      RETURN
      END