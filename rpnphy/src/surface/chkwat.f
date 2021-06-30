      SUBROUTINE CHKWAT(ISFC,PCPR,EVAP,RUNOFF,WLOST,RAICAN,SNOCAN,
     1                  RAICNI,SNOCNI,ZPOND,ZPONDI,THLIQ,THICE,
     2                  THLIQI,THICEI,ZSNOW,RHOSNO,XSNOW,SNOWI,
     3                  WSNOW,WSNOWI,FCS,FGS,FI,BAL,THPOR,THLMIN,
     4                  DELZW,ISAND,IG,ILG,IL1,IL2,JL,N,QDIS,QIN,Isat,
     5                  WTnew,DMLIQ,Ibed,GWTP,Isats,EXCW,
     6                  igwscheme)
C
C     * APR 28/10 - B.DUGAS.    INTRODUCE SEPARATE ACCURACY LIMITS
C     *                         FOR BAL AND FOR THE OTHER CHECKS.
C     * JUN 06/06 - D.VERSEGHY. MODIFY CHECK ON RUNOFF.
C     * APR 03/06 - D.VERSEGHY. ALLOW FOR PRESENCE OF WATER IN SNOW.
C     * SEP 26/05 - D.VERSEGHY. REMOVE HARD CODING OF IG=3 IN 300 LOOP.
C     * SEP 23/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.
C     * AUG 04/04 - D.VERSEGHY. RELAX ACCURACY LIMIT FOR WATER
C     *                         BALANCE CHECK, CONSISTENT WITH
C     *                         ROUNDOFF ERROR CONSTRAINTS.
C     * JUN 25/02 - D.VERSEGHY. RENAME VARIABLES FOR CLARITY; UPDATES
C     *                         CAUSED BY CONVERSION OF PONDING DEPTH
C     *                         TO A PROGNOSTIC VARIABLE; SHORTENED
C     *                         CLASS4 COMMON BLOCK.
C     * MAY 23/02 - D.VERSEGHY. MOVE CALCULATION OF "XSNOW" INTO
C     *                         THIS ROUTINE.
C     * JUN 20/97 - D.VERSEGHY. CLASS - VERSION 2.7.
C     *                         MODIFICATIONS TO ALLOW FOR VARIABLE
C     *                         SOIL PERMEABLE DEPTH.
C     * AUG 24/95 - D.VERSEGHY. CLASS - VERSION 2.4.
C     *                         RATIONALIZE USE OF WLOST.
C     *                         ALSO INTRODUCE NEW VALUE OF ACCLMT
C     *                         CORRESPONDING TO 3 MM/YR AS USED
C     *                         BY THE PILPS COMMUNITY.
C     * AUG 18/95 - D.VERSEGHY. REVISIONS TO ALLOW FOR INHOMOGENEITY
C     *                         BETWEEN SOIL LAYERS.
C     * JAN 31/94 - D.VERSEGHY. LOCAL VERSION FOR CLIMATE RESEARCH
C     *                         NETWORK: CHECK RAICAN AND SNOCAN
C     *                         AGAINST -ACCLMT INSTEAD OF AGAINST 0.
C     * AUG 16/93 - D.VERSEGHY. CLASS - VERSION 2.2.
C     *                         RETURN WATER BALANCE CHECK TO ROUTINE
C     *                         USE (COMMENTED OUT IN PREVIOUS VERSION)
C     *                         AND RENAME SUBROUTINE FROM "CHKVAL"
C     *                         TO "CHKWAT".
C     * MAY 15/92 - M.LAZARE.   CLASS - VERSION 2.1.
C     *                         MOISTURE BALANCE CHECKS EXTRACTED FROM
C     *                         "CLASSW" AND VECTORIZED.
C     * APR 11/89 - D.VERSEGHY. THE FOLLOWING MOISTURE BALANCE CHECKS 
C     *                         ARE CARRIED OUT: INTERCEPTED MOISTURE
C     *                         STORES AND LOCAL RUNOFF MUST BE .GE.0;
C     *                         LIQUID SOIL LAYER MOISTURE STORES MUST
C     *                         BE LESS THAN THE PORE VOLUME AND GREATER
C     *                         THAN THE LIMITING VALUE "THLMIN"; FROZEN
C     *                         SOIL LAYER MOISTURE STORES MUST BE LESS
C     *                         THAN THE MAXIMUM AVAILABLE VOLUME (THE
C     *                         PORE VOLUME - THLMIN) AND GE.0; AND THE
C     *                         MOISTURE BALANCE OF THE TOTAL CANOPY/
C     *                         SNOW/SOIL COLUMN MUST BE WITHIN A 
C     *                         SPECIFIED TOLERANCE.  THE TOLERANCE
C     *                         LEVEL ADOPTED IS DESIGNATED BY "ACCLMT".
C
      IMPLICIT NONE
C
C     * INTEGER CONSTANTS.
C
      INTEGER ISFC,IG,ILG,IL1,IL2,JL,I,J,K
C
      INTEGER IPTBAD,JPTBAD,KPTBAD,IPTBDI,JPTBDI,KPTBDI,LPTBDI,
     1        IPTBDJ,JPTBDJ,KPTBDJ,LPTBDJ,N,Ibed(ILG)
C
C     * INPUT FIELDS.

      REAL PCPR  (ILG),    EVAP  (ILG),    RUNOFF(ILG),    WLOST (ILG),    
     1     RAICAN(ILG),    SNOCAN(ILG),    RAICNI(ILG),    SNOCNI(ILG),    
     2     ZPOND (ILG),    ZPONDI(ILG),    ZSNOW (ILG),    RHOSNO(ILG),    
     3     XSNOW (ILG),    SNOWI (ILG),    WSNOW (ILG),    WSNOWI(ILG),
     4     FCS   (ILG),    FGS   (ILG),    FI    (ILG),    QDIS(ILG),
     5     QIN(ILG),  DMLIQ(ILG),EXCW(ILG)
C
      REAL THLIQ (ILG,IG), THICE (ILG,IG), 
     1     THLIQI(ILG,IG), THICEI(ILG,IG), GWTP(ILG) 
C
C     * WORK ARRAYS.
C
      REAL BAL   (ILG), WTnew(ILG)
C
C     * TEMPORARY VARIABLES.
C
      REAL ACCLMT,BALLMT,CANFAC,SNOFAC
C
C     * SOIL INFORMATION ARRAYS.
C
      REAL THPOR (ILG,IG), THLMIN(ILG,IG), DELZW (ILG,IG)
C
      INTEGER              ISAND (ILG,IG)   , Isat(ILG),
     1                     IGDR(ILG)        , Isats(ILG)
      integer igwscheme
C
C     * COMMON BLOCK PARAMETERS.
C
      REAL DELT,TFREZ,HCPW,HCPICE,HCPSOL,HCPOM,HCPSND,HCPCLY,SPHW,
     1     SPHICE,SPHVEG,SPHAIR,RHOW,RHOICE,TCGLAC,CLHMLT,CLHVAP
C
      COMMON /CLASS1/ DELT,TFREZ
      COMMON /CLASS4/ HCPW,HCPICE,HCPSOL,HCPOM,HCPSND,HCPCLY,
     1                SPHW,SPHICE,SPHVEG,SPHAIR,RHOW,RHOICE,
     2                TCGLAC,CLHMLT,CLHVAP

      if (igwscheme.ne.1) then
      DO 501 I=IL1,IL2
      DO 101 J=1,IG
         IF(ISAND(I,J).GT.-3) THEN
            IGDR(I)=J
         ENDIF
101   continue
501   continue
      endif


C
C      ACCLMT=3.0*DELT/3.1536E7
      ACCLMT=1.0E-2
      BALLMT=1.0E-1
C-----------------------------------------------------------------------
      IF(ISFC.EQ.1 .OR. ISFC.EQ.3)                                  THEN      
          IPTBAD=0
          JPTBAD=0
      ENDIF
      KPTBAD=0
      DO 100 I=IL1,IL2
          IF(FI(I).GT.0. .AND. ISAND(I,1).GT.-4)                   THEN
              IF(ISFC.EQ.1 .OR. ISFC.EQ.3)                   THEN      
                  IF(RAICAN(I).LT.(-1.0*ACCLMT)) IPTBAD=I
                  IF(SNOCAN(I).LT.(-1.0*ACCLMT)) JPTBAD=I
              ENDIF
              IF(RUNOFF(I).LT.(-1.0*ACCLMT)) KPTBAD=I
          ENDIF
  100 CONTINUE
C
      IF(ISFC.EQ.1 .OR. ISFC.EQ.3)                                  THEN      
          IF(IPTBAD.NE.0)                                    THEN
             WRITE(6,6100) IPTBAD,JL,ISFC,RAICAN(IPTBAD)
 6100        FORMAT('0AT (I,JL)=(',I3,',',I3,'), ISFC=',I2,' RAICAN = ',
     1               E13.5)
             print *, 'CHKWAT -1' ! CALL XIT('CHKWAT',-1)
          ENDIF
          IF(JPTBAD.NE.0)                                           THEN
             WRITE(6,6150) JPTBAD,JL,ISFC,SNOCAN(JPTBAD)
 6150        FORMAT('0AT (I,JL)=(',I3,',',I3,'), ISFC=',I2,' SNOCAN = ',
     1               E13.5)
             print *, 'CHKWAT -2' ! CALL XIT('CHKWAT',-2)
          ENDIF
      ENDIF
      IF(KPTBAD.NE.0)                                           THEN
         WRITE(6,6200) KPTBAD,JL,ISFC,RUNOFF(KPTBAD)
 6200    FORMAT('0AT (I,JL)=(',I3,',',I3,'), ISFC=',I2,' RUNOFF = ',
     1           E13.5)
         print *, 'CHKWAT -3' ! CALL XIT('CHKWAT',-3)
      ENDIF
C
      IPTBDI=0
      JPTBDI=0
      KPTBDI=0
      LPTBDI=0
      DO 150 J=1,min(IG,Isat(I))
      DO 150 I=IL1,IL2
          IF(FI(I).GT.0. .AND. ISAND(I,1).GT.-4)                   THEN
              IF((THLIQ(I,J)-THPOR(I,J)).GT.ACCLMT)           THEN
              WRITE(6,6009) I,J,N,THLIQ(I,J),THPOR(I,J),DELZW(I,J),
     1                      THICE(I,J),THLIQI(I,J),THICEI(I,J),
     1                      ZPOND(I),ZPONDI(I),ISAND(I,J)
6009          FORMAT(2X,3I6,8F16.8,I6)
              DO 145 K=1,IG
              WRITE(6,6008) K,THLIQ(I,K),THLIQI(I,K),THICE(I,K),
     1                      THICEI(I,K),DELZW(I,K),THPOR(I,K),ISAND(I,K)
6008          FORMAT(2X,I6,6F14.8,I6)
145           CONTINUE
                  IPTBDI=I
                  IPTBDJ=J
              ENDIF
              IF(THLIQ(I,J).LT.(THLMIN(I,J)-ACCLMT) .AND. 
     1                          ISAND(I,J).NE.-3)             THEN
                  JPTBDI=I
                  JPTBDJ=J
              ENDIF
              IF((THICE(I,J)*RHOICE/RHOW-THPOR(I,J)+THLMIN(I,J))
     1                      .GT.ACCLMT.AND.ISAND(I,J).NE.-3)  THEN 
                  KPTBDI=I
                  KPTBDJ=J
              ENDIF
              IF(THICE(I,J).LT.-1.*ACCLMT)                    THEN
                  LPTBDI=I
                  LPTBDJ=J
              ENDIF
          ENDIF
  150 CONTINUE
C
      IF(IPTBDI.NE.0)                                               THEN
          WRITE(6,6250) IPTBDI,JL,ISFC,N,THLIQ(IPTBDI,IPTBDJ),
     1                  THPOR(IPTBDI,IPTBDJ),IPTBDJ
 6250     FORMAT('0AT (I,JL)=(',I6,',',I6,'), ISFC=',I2,' STEP =',I8,
     1            ' THLIQ = ',E13.5,' THPOR = ',E13.5,' FOR J=',I2)
          print *, 'CHKWAT -4' ! CALL XIT('CHKWAT',-4)
      ENDIF
      IF(JPTBDI.NE.0)                                               THEN
          WRITE(6,6300) JPTBDI,JL,ISFC,THLIQ(JPTBDI,JPTBDJ),JPTBDJ
 6300     FORMAT('0AT (I,JL)=(',I3,',',I3,'), ISFC=',I2,' THLIQ = ',
     1            E13.5,' FOR J=',I2)
          print *, 'CHKWAT -5' ! CALL XIT('CHKWAT',-5)
      ENDIF
      IF(KPTBDI.NE.0)                                               THEN
          WRITE(6,6350) KPTBDI,JL,ISFC,THICE(KPTBDI,KPTBDJ),
     1                  THPOR(KPTBDI,KPTBDJ),KPTBDJ
 6350     FORMAT('0AT (I,JL)=(',I3,',',I3,'), ISFC=',I2,' THICE = ',
     1            E13.5,' THPOR = ',E13.5,' FOR J=',I2)
          print *, 'CHKWAT -6' ! CALL XIT('CHKWAT',-6)
      ENDIF
      IF(LPTBDI.NE.0)                                               THEN
          WRITE(6,6400) LPTBDI,JL,ISFC,THICE(LPTBDI,LPTBDJ),LPTBDJ
 6400     FORMAT('0AT (I,JL)=(',I3,',',I3,'), ISFC=',I2,' THICE = ',
     1            E13.5,' FOR J=',I2)
          print *, 'CHKWAT -7' ! CALL XIT('CHKWAT',-7)
      ENDIF
C
      IPTBAD=0
      IF(ISFC.EQ.1 .OR. ISFC.EQ.3)                                  THEN
          CANFAC=1.0
      ELSE
          CANFAC=0.0
      ENDIF
 
      DO 300 I=IL1,IL2
          IF(FI(I).GT.0. .AND. ZSNOW(I).GT.0.) XSNOW(I)=1.0
          IF(FI(I).GT.0. .AND. ISAND(I,1).GT.-4)                   THEN
             IF(ISFC.EQ.1 .OR. ISFC.EQ.2)         THEN
                SNOFAC=1.0/(FCS(I)+FGS(I))
             ELSE
                SNOFAC=0.0
             ENDIF
             IF(igwscheme.eq.1) then
               BAL(I)=PCPR(I)*DELT-                                               
     1                 EVAP(I)*RHOW*DELT-RUNOFF(I)*RHOW+WLOST(I)-
     2                 CANFAC*(RAICAN(I)-RAICNI(I)+SNOCAN(I)-
     3                 SNOCNI(I))-(ZPOND(I)-ZPONDI(I))*RHOW-                   
     4                 ZSNOW(I)*RHOSNO(I)+SNOFAC*SNOWI(I)-
     5                 WSNOW(I)+WSNOWI(I)
             ELSEIF  (Isats(I).EQ.Ibed(I)) THEN
                 BAL(I)=PCPR(I)*DELT-                                               
     1                 EVAP(I)*RHOW*DELT-RUNOFF(I)*RHOW+WLOST(I)-
     2                 CANFAC*(RAICAN(I)-RAICNI(I)+SNOCAN(I)-
     3                 SNOCNI(I))-(ZPOND(I)-ZPONDI(I))*RHOW-                   
     4                 ZSNOW(I)*RHOSNO(I)+SNOFAC*SNOWI(I)-
     5                 WSNOW(I)+WSNOWI(I)
     6                 -QIN(I)*RHOW*DELT
     8                 +QDIS(I)*RHOW*DELT     ! QDIS is added because of runoff that can be from FDT at the last layer that is baseflow
     7                 +DMLIQ(I)*RHOW+EXCW(I)*RHOW

              ELSE
                 BAL(I)=PCPR(I)*DELT-                                               
     1                 EVAP(I)*RHOW*DELT-RUNOFF(I)*RHOW+WLOST(I)-
     2                 CANFAC*(RAICAN(I)-RAICNI(I)+SNOCAN(I)-
     3                 SNOCNI(I))-(ZPOND(I)-ZPONDI(I))*RHOW-                   
     4                 ZSNOW(I)*RHOSNO(I)+SNOFAC*SNOWI(I)-
     5                 WSNOW(I)+WSNOWI(I)+
     6                 GWTP(I)*DELT+EXCW(I)*RHOW
C     6                 +QIN(I)*RHOW*DELT
C     7                 +QDIS(I)*RHOW*DELT
              ENDIF

!     7                 +MAX(WTnew-5.0,0.0)*RHOW
!              ELSE
!                 BAL(I)=PCPR(I)*DELT-                                               
!     1                 EVAP(I)*RHOW*DELT-RUNOFF(I)*RHOW+WLOST(I)-
!     2                 CANFAC*(RAICAN(I)-RAICNI(I)+SNOCAN(I)-
!     3                 SNOCNI(I))-(ZPOND(I)-ZPONDI(I))*RHOW-                   
!     4                 ZSNOW(I)*RHOSNO(I)+SNOFAC*SNOWI(I)-
!     5                 WSNOW(I)+WSNOWI(I)
c     6                 -QIN(I)*RHOW*DELT
c     7                 +QDIS(I)*RHOW*DELT        
!              ENDIF

                               ! qdis must be used as a positive values it has increased basedflow earlier
             DO 275 J=1,MIN(Isat(I),IG)
                 BAL(I)=BAL(I)-
     1                 (THLIQ(I,J)-THLIQI(I,J))*RHOW*DELZW(I,J)-                       
     2                 (THICE(I,J)-THICEI(I,J))*RHOICE*DELZW(I,J)                     
  275        CONTINUE
             IF(ABS(BAL(I)).GT.BALLMT)                           THEN
                 IPTBAD=I
             ENDIF  
          ENDIF
  300 CONTINUE

      IF(IPTBAD.NE.0)                                            THEN
          WRITE(6,6450) IPTBAD,JL,N,ISFC,BAL(IPTBAD)
          WRITE(6,6460) PCPR(IPTBAD)*DELT,EVAP(IPTBAD)*RHOW*DELT,
     1        RUNOFF(IPTBAD)*RHOW,WLOST(IPTBAD),
     2        RAICNI(IPTBAD)-RAICAN(IPTBAD),SNOCNI(IPTBAD)-
     3        SNOCAN(IPTBAD),(ZPOND(IPTBAD)-ZPONDI(IPTBAD))*RHOW
          WRITE(6,6460) RAICAN(IPTBAD),RAICNI(IPTBAD),
     1         SNOCAN(IPTBAD),SNOCNI(IPTBAD),ZPOND(IPTBAD),
     2         ZPONDI(IPTBAD)
          WRITE(6,6460) ZSNOW(IPTBAD)*RHOSNO(IPTBAD),
     1        SNOFAC*SNOWI(IPTBAD),WSNOW(IPTBAD),WSNOWI(IPTBAD),
     2        SNOFAC*SNOWI(IPTBAD)-ZSNOW(IPTBAD)*RHOSNO(IPTBAD),
     3        WSNOWI(IPTBAD)-WSNOW(IPTBAD)
          WRITE(6,6460) ZSNOW(IPTBAD),RHOSNO(IPTBAD),SNOFAC,
     1        SNOWI(IPTBAD)
          if (igwscheme.ne.1) then
             write(6,*) isat(iptbad),isats(iptbad),ibed(iptbad)
             WRITE(6,*) QIN(IPTBAD)*RHOW*DELT,DMLIQ(IPTBAD)*RHOW,
     1       QDIS(IPTBAD)*RHOW*DELT,GWTP(IPTBAD)*DELT,EXCW(IPTBAD)*RHOW
          endif
          DO 350 J=1,IG
              WRITE(6,6460) 
     1        THLIQ(IPTBAD,J),THLIQI(IPTBAD,J),
     2        THICE(IPTBAD,J),THICEI(IPTBAD,J),
     3        DELZW(IPTBAD,J),THPOR(IPTBAD,J),
     4        (THLIQ(IPTBAD,J)-THLIQI(IPTBAD,J))*RHOW*DELZW(IPTBAD,J),
     5        (THICE(IPTBAD,J)-THICEI(IPTBAD,J))*RHOICE*DELZW(IPTBAD,J)
  350     CONTINUE
          WRITE(6,6470) FCS(IPTBAD),FGS(IPTBAD)
6450      FORMAT('0AT (I,JL)=(',I8,',',I8,'),  TIME=',I8,' ISFC=',I2,
     1        ' BAL = ',E13.5)
6460      FORMAT(2X,8F15.8)
6470      FORMAT(2X,4E20.6)
          print *, 'CHKWAT -8' ! CALL XIT('CHKWAT',-8)
      ENDIF

      RETURN
      END  
