      SUBROUTINE TSOLVE(ISNOW,FI,
     1                  QSWNET,QLWOUT,QTRANS,QSENS,QEVAP,EVAP,
     2                  TZERO,QZERO,GZERO,QMELT,CDH,CDM,RIB,CFLUX,
     3                  FTEMP,FVAP,ILMO,UE,H,
     4                  QSWINV,QSWINI,QLWIN,TPOTA,QA,VA,PADRY,RHOAIR,
     5                  ALVISG,ALNIRG,CRIB,CPHCH,CEVAP,TVIRTA,
     6                  ZOSCLH,ZOSCLM,ZRSLFH,ZRSLFM,ZOH,ZOM,FCOR,
     7                  GCONST,GCOEFF,TSTART,TRSNOW,PCPR,
     8                  IWATER,IEVAP,ITERCT,ISAND,
     9                  ISLFD,ITG,ILG,IG,IL1,IL2,JL,
     A                  TSTEP,TVIRTS,EVBETA,Q0SAT,RESID,
     B                  DCFLXM,CFLUXM,WZERO,TRTOP,A,B,
     C                  LZZ0,LZZ0T,FM,FH,ITER,NITER,JEVAP,KF)
C
C     Purpose: Solution of surface energy balance for non-vegetated 
C     subareas.
C
C     * OCT 14/11 - D.VERSEGHY. FOR POST-ITERATION CLEANUP WITH N-R SCHEME,
C     *                         REMOVE CONDITION INVOLVING LAST ITERATION
C     *                         TEMPERATURE.
C     * DEC 07/09 - D.VERSEGHY. RESTORE EVAPORATION WHEN PRECIPITATION
C     *                         IS OCCURRING.
C     * MAR 13/09 - D.VERSEGHY. REPLACE SURFCON COMMON BLOCK WITH CLASSD2;
C     *                         REVISED CALL TO FLXSURFZ.
C     * JAN 06/09 - D.VERSEGHY/M.LAZARE. SPLIT IF CONDITIONS FRAMING
C     *                         300 LOOP.
C     * FEB 25/08 - D.VERSEGHY. STREAMLINE SOME CALCULATIONS; REMOVE
C     *                         "ILW" SWITCH; SUPPRESS WATER VAPOUR FLUX
C     *                         IF PRECIPITATION IS OCCURRING.
C     * MAY 17/06 - D.VERSEGHY. SUPPRESS EVAPORATION WHEN PONDED WATER
C     *                         IS FREEZING; ADD IL1 AND IL2 TO CALL TO
C     *                         FLXSURFZ; REMOVE JL FROM CALL TO DRCOEF.
C     * APR 13/05 - R.BROWN. ADD WINDLESS TRANFER COEFFICIENT TO QSENS
C     *                         CALCULATION FOR SNOW PACKS.
C     * DEC 17/04 - Y.DELAGE/D.VERSEGHY. ADD SWITCH TO USE EITHER SECANT/
C     *                         BISECTION OR NEWTON-RAPHSON ITERATION
C     *                         SCHEME (WITH NUMBER OF ITERATIONS LIMITED
C     *                         TO FIVE AND CORRECTION FOR REMAINING
C     *                         RESIDUAL).
C     * NOV 04/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.
C     * AUG 06/04 - Y.DELAGE/D.VERSEGHY. PROTECT SENSITIVE CALCULATIONS
C     *                         FROM ROUNDOFF ERRORS.
C     * NOV 07/02 - Y.DELAGE/D.VERSEGHY. NEW CALL TO FLXSURFZ.
C     * JUL 26/02 - D.VERSEGHY. SHORTENED CLASS4 COMMON BLOCK.
C     * MAR 28/02 - D.VERSEGHY. STREAMLINED SUBROUTINE CALL.
C     *                         BYPASS EVAPORATION EFFICIENCY PARAMETER 
C     *                         IN CASES OF CONDENSATION.
C     * JAN 18/02 - P.BARTLETT/D.VERSEGHY. NEW "BETA" FORMULATION FOR 
C     *                         BARE SOIL EVAPORATION BASED ON LEE AND
C     *                         PIELKE.
C     * APR 11/01 - M.LAZARE.   SHORTENED "CLASS2" COMMON BLOCK.
C     * OCT 06/00 - D.VERSEGHY. CONDITIONAL "IF" IN ITERATION SEQUENCE
C     *                         TO AVOID DIVIDE BY ZERO.
C     * DEC 07/99 - A.WU/D.VERSEGHY. NEW SOIL EVAPORATION FORMULATION.
C     * JUL 24/97 - D.VERSEGHY. CLASS - VERSION 2.7.
C     *                         REPLACE BISECTION METHOD IN SURFACE 
C     *                         TEMPERATURE ITERATION SCHEME WITH 
C     *                         SECANT METHOD FOR FIRST TEN ITERATIONS.
C     *                         PASS QZERO,QA,ZOMS,ZOHS TO REVISED
C     *                         DRCOEF (ZOMS AND ZOHS ALSO NEW WORK ARRAYS
C     *                         PASSED TO THIS ROUTINE).
C     * JUN 20/97 - D.VERSEGHY. PASS IN NEW "CLASS4" COMMON BLOCK.
C     * JAN 02/96 - D.VERSEGHY. CLASS - VERSION 2.5.
C     *                         COMPLETION OF ENERGY BALANCE 
C     *                         DIAGNOSTICS.  ALSO, PASS SWITCH "ILW"
C     *                         THROUGH SUBROUTINE CALL, SPECIFYING 
C     *                         WHETHER QLWIN REPRESENTS INCOMING
C     *                         (ILW=1) OR NET (ILW=2) LONGWAVE
C     *                         RADIATION ABOVE THE GROUND.
C     * NOV 30/94 - M.LAZARE.   CLASS - VERSION 2.3.
C     *                         NEW DRAG COEFFICIENT AND RELATED FIELDS,
C     *                         NOW DETERMINED IN ROUTINE "DRCOEF"
C     *                         "CFLUX" NOW WORK FIELD INSTEAD OF "CLIMIT".
C     * OCT 04/94 - D.VERSEGHY. CHANGE "CALL ABORT" TO "CALL XIT" TO
C     *                         ENABLE RUNNING ON PCS.
C     * JAN 24/94 - M.LAZARE.   UNFORMATTED I/O COMMENTED OUT IN LOOP 200.
C     * JUL 29/93 - D.VERSEGHY. CLASS - VERSION 2.2.
C     *                         REMOVE RE-DEFINITION OF QMELT NEAR END
C     *                         (SINCE DONE ELSEWHERE ALREADY) AND
C     *                         REDEFINE QSWNET FOR DIAGNOSTIC PURPOSES
C     *                         TO INCLUDE TRANSMISSION THROUGH 
C     *                         SNOWPACK.
C     * OCT 15/92 - D.VERSEGHY/M.LAZARE. CLASS - VERSION 2.1.
C     *                                  REVISED AND VECTORIZED CODE 
C     *                                  FOR MODEL VERSION GCM7.                  
C     * AUG 12/91 - D.VERSEGHY. CODE FOR MODEL VERSION GCM7U -
C     *                         CLASS VERSION 2.0 (WITH CANOPY).
C     * APR 11/89 - D.VERSEGHY. ITERATIVE SURFACE TEMPERATURE 
C     *                         CALCULATIONS FOR SNOW/SOIL.
C
      IMPLICIT NONE

C     * INTEGER CONSTANTS.
C
      INTEGER ISNOW !Flag indicating presence or absence of snow
      INTEGER ISLFD,ITG,ILG,IG,IL1,IL2,JL,I
C
      INTEGER NUMIT,NIT,IBAD,ITERMX
C
C     * OUTPUT ARRAYS.
C
      REAL QSWNET(ILG)  !Net shortwave radiation at surface [W m-2]    
      REAL QLWOUT(ILG)  !Upwelling longwave radiation at surface [W m-2]
                        ! (L)  
      REAL QTRANS(ILG)  !Shortwave radiation transmitted into surface 
                        ![W m-2]
      REAL QSENS (ILG)  !Sensible heat flux from surface [W m-2] (QH )  
      REAL QEVAP (ILG)  !Latent heat flux from surface [W m-2] (QE)  
      REAL EVAP  (ILG)  !Evaporation rate at surface [kg m-2 s-1] (E(0))  
      REAL TZERO (ILG)  !Temperature at surface [K] (T(0))  
      REAL QZERO (ILG)  !Specific humidity at surface [kg kg-1] (q(0)) 
      REAL GZERO (ILG)  !Heat flux into surface [W m-2] (G(0))  
      REAL QMELT (ILG)  !Heat available for melting snow or freezing 
                        !water at the surface [W m-2]
      REAL CDH   (ILG)  !Surface drag coefficient for heat [ ] (CDH)  
      REAL CDM   (ILG)  !Surface drag coefficient for momentum [ ]  
      REAL RIB   (ILG)  !Bulk Richardson number at surface [ ]  
      REAL CFLUX (ILG)  !Product of surface drag coefficient and wind 
                        !speed [m s-1]
      REAL FTEMP (ILG)  !Product of surface-air temperature gradient, 
                        !drag coefficient and wind speed [K m s-1]
      REAL FVAP  (ILG)  !Product of surface-air humidity gradient, drag 
                        !coefficient and wind speed [kg kg-1 m s-1]
      REAL ILMO  (ILG)  !Inverse of Monin-Obukhov roughness length (m-1]  
      REAL UE    (ILG)  !Friction velocity of air [m s-1]  
      REAL H     (ILG)  !Height of the atmospheric boundary layer [m]
C
C     * INPUT ARRAYS.
C
      REAL FI    (ILG)  !Fractional coverage of subarea in question on 
                        !modelled area [ ]
      REAL QSWINV(ILG)  !Visible radiation incident on horizontal 
                        !surface [W m-2]
      REAL QSWINI(ILG)  !Near-infrared radiation incident on horizontal 
                        !surface [W m -2]
      REAL QLWIN (ILG)  !Downwelling longwave radiation at bottom of 
                        !atmosphere [W m-2]
      REAL TPOTA (ILG)  !Potential temperature of air at reference 
                        !height [K] (Ta,pot)
      REAL QA    (ILG)  !Specific humidity at reference height [kg kg-1] 
                        !(qa)
      REAL VA    (ILG)  !Wind speed at reference height [m s-1] (va)  
      REAL PADRY (ILG)  !Partial pressure of dry air [Pa] (pdry)  
      REAL RHOAIR(ILG)  !Density of air [kg m-3] (rho_a)  
      REAL ALVISG(ILG)  !Visible albedo of ground surface [ ]  
      REAL ALNIRG(ILG)  !Near-IR albedo of ground surface [ ]  
      REAL CRIB  (ILG)  !Richardson number coefficient [K-1]  
      REAL CPHCH (ILG)  !Latent heat of vaporization at surface [J kg-1]  
      REAL CEVAP (ILG)  !Soil evaporation efficiency coefficient [ ] 
                        !(beta)
      REAL TVIRTA(ILG)  !Virtual potential temperature of air at 
                        !reference height [K]
      REAL ZOSCLH(ILG)  !Ratio of roughness length for heat to reference 
                        !height for temperature and humidity [ ]
      REAL ZOSCLM(ILG)  !Ratio of roughness length for momentum to 
                        !reference height for wind speed [ ]
      REAL ZRSLFH(ILG)  !Difference between reference height for 
                        !temperature and humidity and height at which 
                        !extrapolated wind speed goes to zero [m]
      REAL ZRSLFM(ILG)  !Difference between reference height for wind 
                        !speed and height at which extrapolated wind 
                        !speed goes to zero [m]
      REAL ZOH   (ILG)  !Surface roughness length for heat [m]  
      REAL ZOM   (ILG)  !Surface roughness length for momentum [m]  
      REAL GCONST(ILG)  !Intercept used in equation relating surface 
                        !heat flux to surface temperature [W m-2 ]
      REAL GCOEFF(ILG)  !Multiplier used in equation relating surface 
                        !heat flux to surface temperature [W m-2 K-1]
      REAL TSTART(ILG)  !Starting point for surface temperature 
                        !iteration [K]  
      REAL TRSNOW(ILG)  !Short-wave transmissivity of snow pack [ ]  
      REAL FCOR  (ILG)  !Coriolis parameter [s-1]  
      REAL PCPR  (ILG)  !Surface precipitation rate [kg m-2 s-1]

C
      INTEGER          IWATER(ILG)  !Flag indicating condition of 
                                    !surface (dry, water-covered or 
                                    !snow-covered)
      INTEGER          IEVAP (ILG)  !Flag indicating whether surface 
                                    !evaporation is occurring or not
      INTEGER          ITERCT(ILG,6,50) !Counter of number of iterations 
                                        !required to solve energy 
                                        !balance for four subareas
      INTEGER          ISAND(ILG,IG)    !Sand content flag
C
C     * INTERNAL WORK ARRAYS.
C
      REAL TSTEP (ILG),    TVIRTS(ILG),    EVBETA(ILG),    Q0SAT (ILG),
     1     RESID (ILG),    DCFLXM(ILG),    CFLUXM(ILG),    TRTOP (ILG),    
     2     A     (ILG),    B     (ILG),
     3     LZZ0  (ILG),    LZZ0T (ILG),    FM    (ILG),    FH    (ILG),
     4     WZERO (ILG)
C
      INTEGER              ITER  (ILG),    NITER (ILG),    JEVAP (ILG),
     1                     KF    (ILG)
C
C     * TEMPORARY VARIABLES.
C
      REAL QSWNV,QSWNI,DCFLUX,DRDT0,TZEROT,QEVAPT,BOWEN,EZERO
C
C     * COMMON BLOCK PARAMETERS.
C
      REAL DELT     !Time step [s]
      REAL TFREZ    !Freezing point of water [K]
      REAL RGAS     !Gas Constant [J kg-1 K-1]
      REAL RGASV    !Gas constant for water vapour [J kg-1 K-1]
      REAL GRAV     !Acceleration due to gravity [m s-1]
      REAL SBC      !Stefan-Boltzmann constant [W m-2 K-4]
      REAL VKC      !Von Karman constant (0.40)
      REAL CT       !Drag coefficient for water (1.15*10^-3)
      REAL VMIN     !Minimum wind speed (0.1) [m s-1]
      REAL HCPW     !Volumetric heat capacity of water (4.187*10^6) 
                    ![J m-3 K-1]
      REAL HCPICE   !Volumetric heat capacity of ice (1.9257*10^6) 
                    ![J m-3 K-1]
      REAL HCPSOL   !Volumetric heat capacity of mineral matter 
                    !(2.25*10^6) [J m-3 K-1]
      REAL HCPOM    !Volumetric heat capacity of organic matter 
                    !(2.50*10^6) [J m-3 K-1]
      REAL HCPSND   !Volumetric heat capacity of sand particles 
                    !(2.13*10^6) [J m-3 K-1]
      REAL HCPCLY   !Volumetric heat capacity of fine mineral particles 
                    !(2.38*10^6) [J m-3 K-1]
      REAL SPHW     !Specific heat of water (4.186*10^3) [J kg-1 K-1]
      REAL SPHICE   !Specific heat of ice (2.10*10^3) [J kg-1 K-1]
      REAL SPHVEG   !Specific heat of vegetation matter (2.70*10^3) 
                    ![J kg-1 K-1]
      REAL SPHAIR   !Specific heat of air [J kg-1 K-1]
      REAL RHOW     !Density of water (1.0*10^3) [kg m-3]
      REAL RHOICE   !Density of ice (0.917*10^3) [kg m-3]
      REAL TCGLAC   !Thermal conductivity of ice sheets (2.24) 
                    ![W m-1 K-1]
      REAL CLHMLT   !Latent heat of freezing of water (0.334*10^6) 
                    ![J kg-1]
      REAL CLHVAP   !Latent heat of vaporization of water (2.501*10^6) 
                    ![J kg-1]
      REAL DELTA,CGRAV,CKARM,CPD,AS,ASX,CI,BS,BETA,FACTN,HMIN,ANGMAX
C
      COMMON /CLASS1/ DELT,TFREZ                                                  
      COMMON /CLASS2/ RGAS,RGASV,GRAV,SBC,VKC,CT,VMIN
      COMMON /CLASS4/ HCPW,HCPICE,HCPSOL,HCPOM,HCPSND,HCPCLY,
     1                SPHW,SPHICE,SPHVEG,SPHAIR,RHOW,RHOICE,
     2                TCGLAC,CLHMLT,CLHVAP
      COMMON /PHYCON/ DELTA,CGRAV,CKARM,CPD
      COMMON /CLASSD2/ AS,ASX,CI,BS,BETA,FACTN,HMIN,ANGMAX
C-----------------------------------------------------------------------
      !
      !For the surface temperature iteration, two alternative schemes 
      !are offered: the bisection method (selected if the flag ITG = 1) 
      !and the Newton-Raphson method (selected if ITG = 2). In the first 
      !case, the maximum number of iterations ITERMX is set to 12, and 
      !in the second case it is set to 5. An optional windless transfer 
      !coefficient EZERO is available, which can be used, following the 
      !recommendations of Brown et al. (2006), to prevent the sensible 
      !heat flux over snow packs from becoming vanishingly small under 
      !highly stable conditions. If the snow cover flag ISNOW is zero 
      !(indicating bare ground), EZERO is set to zero; if ISNOW=1, EZERO 
      !is set to 2.0 W m -2 K-1.
      !
C     * INITIALIZATION AND PRE-ITERATION SEQUENCE.
C
      IF(ITG.LT.2) THEN
          ITERMX=12
      ELSE
          ITERMX=5
      ENDIF
C
C      IF(ISNOW.EQ.0) THEN
C          EZERO=0.0
C      ELSE
C          EZERO=2.0
C      ENDIF
       EZERO=0.0
C
      !
      !In the 50 loop, some preliminary calculations are done. The 
      !shortwave transmissivity at the surface, TRTOP, is set to zero in 
      !the absence of a snow pack, and to the transmissivity of snow 
      !TRSNOW otherwise. The net shortwave radiation at the surface, 
      !QSWNET, is calculated as the sum of the incoming visible and 
      !near-infrared shortwave radiation, weighted according to one 
      !minus their respective albedos. This average value is corrected 
      !for the amount of radiation transmitted into the surface, 
      !obtained using TRTOP. The initial value of the surface 
      !temperature TZERO is set to TSTART, which contains the value of 
      !TZERO from the previous time step, and the first step in the 
      !iteration sequence, TSTEP, is set to 1.0 K. The flag ITER is set 
      !to 1 for each element of the set of modelled areas, indicating 
      !that its surface temperature has not yet been found. The 
      !iteration counter NITER is initialized to 1 for each element. 
      !Initial values are assigned to several other variables.
      !
      DO 50 I=IL1,IL2
          IF(FI(I).GT.0.)                                          THEN
              IF(ISNOW.EQ.0)                      THEN
                  TRTOP(I)=0.
              ELSE
                  TRTOP(I)=TRSNOW(I)
              ENDIF
              QSWNV=QSWINV(I)*(1.0-ALVISG(I))   
              QSWNI=QSWINI(I)*(1.0-ALNIRG(I))  
              QSWNET(I)=QSWNV+QSWNI           
              QTRANS(I)=QSWNET(I)*TRTOP(I)   
              QSWNET(I)=QSWNET(I)-QTRANS(I) 
              TZERO(I)=TSTART(I)            
              TSTEP(I)=1.0
              ITER(I)=1
              NITER(I)=1
C
              QMELT(I)=0.0                   
              RESID(I)=999999.
              DCFLXM(I)=0.0
              CFLUX(I)=0.0
              IF(ISNOW.EQ.1)                      THEN
                  KF(I)=3
              ELSE
                  KF(I)=6
              ENDIF
          ENDIF
   50 CONTINUE
      !
      !The 100 continuation line marks the beginning of the surface 
      !temperature iteration sequence. First the flags NIT (indicating 
      !that there are still locations at the beginning of the current 
      !iteration step for which the surface temperature has not yet been 
      !found) and NUMIT (indicating that there are still locations at 
      !the end of the current iteration step for which the surface 
      !temperature has not yet been found) are set to zero. Loop 150 is 
      !then performed over the set of modelled areas. If ITER=1, NIT is 
      !incremented by one, and the initial value of the surface transfer 
      !coefficient CFLUXM for this iteration pass is set to its value 
      !from the previous pass. The virtual temperature at the surface, 
      !T(0)v, is obtained using the standard expression (see 
      !documentation for subroutine CLASST):
      !
      !T(0)v = TZERO*[1 + 0.61*QZERO]
      !
      !where TZERO is the surface temperature and QZERO is the specific 
      !humidity at the surface. The surface humidity can be obtained 
      !from the saturated specific humidity q(0)sat by making use of the 
      !definition of the surface evaporation efficiency CEVAP:
      !
      !CEVAP = [QZERO – QA]/[q(0)sat - QA]
      !
      !where QA is the specific humidity of the air at the reference 
      !height. This expression is inverted to obtain an expression for 
      !QZERO. The saturated specific humidity q(0)sat is determined from 
      !the mixing ratio at saturation, w(0)sat:
      !
      !q(0)sat = w(0)sat/[1 + w(0) sat]
      !
      !The saturation mixing ratio is a function of the saturation 
      !vapour pressure e(0)sat at the surface:
      !
      !w(0)sat = 0.622*e(0)sat/(PADRY)
      !
      !where PADRY is the partial pressure of dry air. A standard 
      !empirical equation for the saturation vapour pressure dependence 
      !on the temperature T is used:
      !
      !esat = 611.0*exp[17.269*(T – Tf)/(T – 35.86)]    T >= Tf
      !esat = 611.0*exp[21.874*(T – Tf)/(T – 7.66)]     T < Tf
      !
      !where Tf is the freezing point. If there is a snow cover or 
      !ponded water present on the surface (IWATER > 0), the surface 
      !evaporation efficiency EVBETA is set to 1 and QZERO is set to 
      !q(0)sat. Otherwise EVBETA is set to CEVAP, the value obtained in 
      !subroutine TPREP on the basis of ambient conditions, and QZERO is 
      !calculated as above. If QZERO > QA and the evaporation flag IEVAP 
      !has been set to zero, EVBETA is reset to zero and QZERO is reset 
      !to QA. Finally, T(0)v is determined using the equation above.
C
C     * ITERATION SECTION.
C     * LOOP IS REPEATED UNTIL SOLUTIONS HAVE BEEN FOUND FOR ALL POINTS 
C     * ON THE CURRENT LATITUDE CIRCLE(S). 
C  
  100 CONTINUE
C
      NUMIT=0
      NIT=0
      DO 150 I=IL1,IL2
          IF(FI(I).GT.0. .AND. ITER(I).EQ.1)                       THEN    
              NIT=NIT+1
              CFLUXM(I)=CFLUX(I)
              IF(TZERO(I).GE.TFREZ)                        THEN
                  A(I)=17.269       
                  B(I)=35.86       
              ELSE                
                  A(I)=21.874    
                  B(I)=7.66     
              ENDIF                       
              WZERO(I)=0.622*611.0*EXP(A(I)*(TZERO(I)-TFREZ)/
     1              (TZERO(I)-B(I)))/PADRY(I)           
              Q0SAT(I)=WZERO(I)/(1.0+WZERO(I))    
              IF(IWATER(I).GT.0)                              THEN
                  EVBETA(I)=1.0
                  QZERO(I)=Q0SAT(I)
              ELSE
                  EVBETA(I)=CEVAP(I)
                  QZERO(I)=EVBETA(I)*Q0SAT(I)+(1.0-EVBETA(I))*QA(I)
                  IF(QZERO(I).GT.QA(I) .AND. IEVAP(I).EQ.0) THEN
                      EVBETA(I)=0.0
                      QZERO(I)=QA(I)
                  ENDIF
              ENDIF
              TVIRTS(I)=TZERO(I)*(1.0+0.61*QZERO(I))
          ENDIF
  150 CONTINUE      
C
      !If NIT > 0, a subroutine is called to evaluate the stability-
      !corrected surface drag coefficients for heat and momentum. The 
      !subroutine selection is made on the basis of the flag ISLFD. If 
      !ISLFD=1, indicating that the calculations are to be consistent 
      !with CCCma conventions, subroutine DRCOEF is called; if ISLFD=2, 
      !indicating that the calculations are to be consistent with RPN 
      !conventions, subroutine FLXSURFZ is called.
      !
      IF(NIT.GT.0)                                                  THEN
C
C     * CALCULATE SURFACE DRAG COEFFICIENTS (STABILITY-DEPENDENT) AND
C     * OTHER RELATED QUANTITIES.
C
        IF(ISLFD.LT.2) THEN
            CALL DRCOEF (CDM,CDH,RIB,CFLUX,QZERO,QA,ZOSCLM,ZOSCLH,
     1                   CRIB,TVIRTS,TVIRTA,VA,FI,ITER,
     2                   ILG,IL1,IL2)
        ELSE
            CALL FLXSURFZ(CDM,CDH,CFLUX,RIB,FTEMP,FVAP,ILMO,
     1                    UE,FCOR,TPOTA,QA,ZRSLFM,ZRSLFH,VA,
     2                    TZERO,QZERO,H,ZOM,ZOH,
     3                    LZZ0,LZZ0T,FM,FH,ILG,IL1,IL2,FI,ITER,JL )
        ENDIF
C
C     * REMAINING CALCULATIONS.
C
        !In loop 175, the terms of the surface energy balance are 
        !evaluated. The energy balance equation is written as:
        !
        !K* + L* - QSENS – QEVAP – GZERO = 0
        !
        !where K* is the net shortwave radiation, L* is the net longwave 
        !radiation, QSENS is the sensible heat flux, QEVAP is the latent 
        !heat flux, and GZERO is the conduction into the surface. K* was 
        !evaluated earlier in loop 50. L* is obtained as the difference 
        !between the downwelling radiation QLWIN and the upwelling 
        !radiation QLWOUT, which in turn is determined using the Stefan-
        !Boltzmann equation:
        !
        !QLWOUT = SBC*TZERO^4
        !
        !where SBC is the Stefan-Boltzmann constant. (It is assumed that 
        !natural surfaces, because of their radiative complexity, act as 
        !effective black bodies, so that their emissivity can be taken 
        !to be 1.) The sensible heat flux is given by
        !
        !QH = [RHOAIR*SPHAIR*CDH*VA + epsilon_0]*[TZERO – TPOTA]
        !
        !where RHOAIR is the density of the air, SPHAIR is its specific 
        !heat, CDH is the surface drag coefficient for heat, and VA and 
        !TPOTA are the wind speed and potential temperature respectively 
        !at the reference height. (Note in the code that the variable 
        !CFLUX, evaluated in subroutine DRCOEF or FLXSURFZ, represents 
        !the product of CDH and va.) The windless transfer coefficient 
        !epsilon_0, evaluated at the beginning of the subroutine, is 
        !used only under stable conditions, i.e. if TZERO < TPOTA. The 
        !evaporation rate at the surface, EVAP, is calculated as
        !
        !EVAP = RHOAIR*CDH*VA [QZERO – QA]
        !
        !QE is obtained by multiplying EVAP by the latent heat of 
        !vaporization at the surface. The ground heat flux GZERO is 
        !determined as a linear function of TZERO (see documentation for 
        !subroutines TNPREP and TSPREP). It can be seen that each of the 
        !terms of the surface energy balance is a function of a single 
        !unknown, TZERO. The residual RESID of the energy balance is now 
        !evaluated on the basis of the current estimation for TZERO. If 
        !the absolute value of RESID is less than 5.0 W m-2, or if the 
        !absolute value of the iteration step TSTEP most recently used 
        !is less than 0.01 K, the surface temperature is deemed to have 
        !been found and ITER is set to 0. If the iteration counter NITER 
        !is equal to the maximum number and ITER is still 1, ITER is set 
        !to -1.
        !
        DO 175 I=IL1,IL2
          IF(FI(I).GT.0. .AND. ITER(I).EQ.1)                       THEN    
              QLWOUT(I)=SBC*TZERO(I)*TZERO(I)*TZERO(I)*TZERO(I)
              IF(TZERO(I).LT.TPOTA(I))                        THEN
                  QSENS(I)=(RHOAIR(I)*SPHAIR*CFLUX(I)+EZERO)*(TZERO(I)-
     1                TPOTA(I))
              ELSE
                  QSENS(I)=RHOAIR(I)*SPHAIR*CFLUX(I)*(TZERO(I)-
     1                TPOTA(I))
              ENDIF
              EVAP(I)=RHOAIR(I)*CFLUX(I)*(QZERO(I)-QA(I)) 
              QEVAP(I)=CPHCH(I)*EVAP(I)      
              GZERO(I)=GCOEFF(I)*TZERO(I)+GCONST(I)
              RESID(I)=QSWNET(I)+QLWIN(I)-QLWOUT(I)-QSENS(I)-QEVAP(I)-
     1                 GZERO(I)                             
              IF(ABS(RESID(I)).LT.5.0)                       ITER(I)=0
              IF(ABS(TSTEP(I)).LT. 1.0E-2)                   ITER(I)=0
              IF(NITER(I).EQ.ITERMX .AND. ITER(I).EQ.1)      ITER(I)=-1
          ENDIF
175     CONTINUE
      ENDIF
C
      !
      !If ITG = 2, then if NIT > 0 and if ITER for the array element in 
      !question is 1, the calculations for the Newton-Raphson method of 
      !iteration are performed. In this approach, the value xn+1 used at 
      !each iteration step is obtained from the value xn at the previous 
      !step as follows:
      !
      !xn+1 = xn – f(xn)/f′(xn)
      !
      !Identifying xn with TZERO and f(xn) with the surface energy 
      !balance equation, it can be seen that the second term on the 
      !right-hand side corresponds to TSTEP; the numerator is equal to 
      !RESID and the denominator to the first derivative of the energy 
      !balance equation evaluated at TZERO, which in turn is equal to 
      !the sum of the derivatives of the individual terms:
      !
      !d(QLWOUT)/dT = -4*SBC*TZERO^3
      !
      !d(QSENS)/dT = RHOAIR*SPHAIR*{CDH*VA + 
      !                             [TZERO – TPOTA]*d(CDH*VA)/dT}
      !
      !d(QEVAP)/dT = Lv*RHOAIR{CDH*VA*d(QZERO)/dT+ 
      !                        [QZERO – QA]*d(CDH*VA)/dT}
      !
      !and dG(0)/dT is equal to the coefficient multiplying TZERO in the 
      !equation for G(0). (Lv is the latent heat of vaporization at the 
      !surface.) At the end of the calculations the iteration counter 
      !NITER and the flag NUMIT are each incremented by one, and upon 
      !exiting the loop, if NUMIT > 0, the iteration cycle is repeated 
      !from line 100 on.
      !
      IF(ITG.LT.2) THEN
C
C     * OPTION #1: BISECTION ITERATION METHOD.
C
      IF(NIT.GT.0)                                                  THEN
        DO 180 I=IL1,IL2      
          IF(FI(I).GT.0. .AND. ITER(I).EQ.1)                       THEN    
              IF(NITER(I).EQ.1) THEN
                  IF(RESID(I).GT.0.0) THEN
                      TZERO(I)=TZERO(I)+1.0
                  ELSE
                      TZERO(I)=TZERO(I)-1.0
                  ENDIF
              ELSE
                  IF((RESID(I).GT.0. .AND. TSTEP(I).LT.0.) .OR.
     1                (RESID(I).LT.0. .AND. TSTEP(I).GT.0.))    THEN 
                      TSTEP(I)=-TSTEP(I)/2.0     
                  ENDIF
                  TZERO(I)=TZERO(I)+TSTEP(I)
              ENDIF
          ENDIF
C
          IF(FI(I).GT.0. .AND. ITER(I).EQ.1)                       THEN
              NITER(I)=NITER(I)+1
              NUMIT=NUMIT+1
          ENDIF
  180   CONTINUE
      ENDIF
C
c     DO 185 I=IL1,IL2
C         IF(FI(I).GT.0. .AND. ITER(I).EQ.-1)                      THEN 
C             WRITE(6,6250) I,JL,RESID(I),TZERO(I),RIB(I)
C6250         FORMAT('0GROUND ITERATION LIMIT',3X,2I3,3(F8.2,E12.4))            
C         ENDIF                              
c 185 CONTINUE
C
      IF(NUMIT.GT.0)                                    GO TO 100
C
      ELSE
C
C     * OPTION #2: NEWTON-RAPHSON ITERATION METHOD.
C
      IF(NIT.GT.0)                                                  THEN
        DO 190 I=IL1,IL2
          IF(FI(I).GT.0. .AND. ITER(I).EQ.1)                      THEN
              IF(NITER(I).GT.1)                                 THEN
                  DCFLUX=(CFLUX(I)-CFLUXM(I))/
     1                SIGN(MAX(.001,ABS(TSTEP(I))),TSTEP(I))
                  IF(ABS(TVIRTA(I)-TVIRTS(I)).LT.0.4)
     1                DCFLUX=MAX(DCFLUX,0.8*DCFLXM(I))
                  DCFLXM(I)=DCFLUX
              ELSE
                  DCFLUX=0.
              ENDIF
              DRDT0= -4.0*SBC*TZERO(I)**3
     1           -RHOAIR(I)*SPHAIR*(CFLUX(I)+MAX(0.,TZERO(I)-TPOTA(I))
     2           *DCFLUX) -GCOEFF(I)
     3           +CPHCH(I)*RHOAIR(I)*(CFLUX(I)*WZERO(I)*A(I)
     4           *EVBETA(I)*(B(I)-TFREZ)/((TZERO(I)-B(I))*
     5           (1.0+WZERO(I)))**2-(QZERO(I)-QA(I))*DCFLUX)
              TSTEP(I)=-RESID(I)/DRDT0
              TSTEP(I)=MAX(-10.,MIN(5.,TSTEP(I)))
              TZERO(I)=TZERO(I)+TSTEP(I)
              NITER(I)=NITER(I)+1
              NUMIT=NUMIT+1
          ENDIF
  190   CONTINUE
      ENDIF
C
      IF(NUMIT.GT.0)                                    GO TO 100
      !
      !After the iteration has been completed, NUMIT is reset to zero 
      !and a check is carried out to ascertain whether convergence has 
      !not been reached (i.e. whether ITER = -1) for any location. In 
      !such cases it is assumed that conditions of near-neutral 
      !stability at the surface are the cause of the difficulty in 
      !finding a solution. A trial value of TZERO is calculated using 
      !the virtual potential temperature of the air. If 
      !RESID > 50 W m-2, TZERO is set to this trial value. The values of 
      !QZERO and the components of the surface energy balance are 
      !recalculated as above, except that QSENS and QEVAP are assumed as 
      !a first approximation to be zero. RESID is determined on this 
      !basis. If RESID is positive, it is assigned to QEVAP; otherwise 
      !RESID is divided equally between QSENS and QEVAP, except in the 
      !case of an absolutely dry surface, in which case QEVAP is set to 
      !zero and QSENS to RESID. RESID is reset to zero, EVAP is obtained 
      !from QEVAP, and T(0)v is recalculated. The flag JEVAP for the 
      !location is set to 1, and NUMIT is incremented by 1. If NUMIT > 0 
      !at the end of this loop, DRCOEF or FLXSURF are called again, and 
      !their calculations are performed for any location where JEVAP is 
      !1, to ensure consistency with the new surface temperature and 
      !humidity.
C
C     * IF CONVERGENCE HAS NOT BEEN REACHED, CALCULATE TEMPERATURE AND
C     * FLUXES ASSUMING NEUTRAL STABILITY.
C
      DO 195 I=IL1,IL2
          NUMIT=0
          JEVAP(I)=0
          IF(FI(I).GT.0. .AND.ITER(I).EQ.-1)                       THEN
              TZEROT=TVIRTA(I)/(1.0+0.61*QZERO(I))
              IF(ABS(RESID(I)).GT.50.) THEN
                  TZERO(I)=TZEROT
                  IF(TZERO(I).GE.TFREZ)                        THEN
                      A(I)=17.269
                      B(I)=35.86
                  ELSE
                      A(I)=21.874
                      B(I)=7.66
                  ENDIF
                  WZERO(I)=0.622*611.0*EXP(A(I)*(TZERO(I)-TFREZ)/
     1                (TZERO(I)-B(I)))/PADRY(I)
                  Q0SAT(I)=WZERO(I)/(1.0+WZERO(I))
                  QZERO(I)=EVBETA(I)*Q0SAT(I)+(1.0-EVBETA(I))*QA(I)
                  QLWOUT(I)=SBC*TZERO(I)*TZERO(I)*TZERO(I)*TZERO(I)
                  GZERO(I)=GCOEFF(I)*TZERO(I)+GCONST(I)
                  RESID(I)=QSWNET(I)+QLWIN(I)-QLWOUT(I)-GZERO(I)
                  IF(RESID(I).GT.0.)                 THEN
                      QEVAP(I)=RESID(I)
                  ELSE
                      QEVAP(I)=RESID(I)*0.5
                  ENDIF
                  QSENS(I)=RESID(I)-QEVAP(I)
                  RESID(I)=0.
                  EVAP(I)=QEVAP(I)/CPHCH(I)
                  TVIRTS(I)=TZERO(I)*(1.0+0.61*QZERO(I))
                  JEVAP(I)=1
                  NUMIT=NUMIT+1
              ENDIF
          ENDIF
  195 CONTINUE
C
      IF(NUMIT.GT.0)                   THEN
        IF(ISLFD.LT.2) THEN
            CALL DRCOEF (CDM,CDH,RIB,CFLUX,QZERO,QA,ZOSCLM,ZOSCLH,
     1                   CRIB,TVIRTS,TVIRTA,VA,FI,JEVAP,
     2                   ILG,IL1,IL2)
        ELSE
            CALL FLXSURFZ(CDM,CDH,CFLUX,RIB,FTEMP,FVAP,ILMO,
     1                    UE,FCOR,TPOTA,QA,ZRSLFM,ZRSLFH,VA,
     2                    TZERO,QZERO,H,ZOM,ZOH,
     3                    LZZ0,LZZ0T,FM,FH,ILG,IL1,IL2,FI,JEVAP,JL )
        ENDIF
      ENDIF
C
      ENDIF
      !
      !At this point a check is performed for unphysical values of the 
      !surface temperature, i.e. for values greater than 100 C or less 
      !than -100 C. If such values are encountered, an error message is 
      !printed and a call to abort is carried out.
C
C     * CHECK FOR BAD ITERATION TEMPERATURES.
C
      IBAD=0
      DO 200 I=IL1,IL2
          IF(FI(I).GT.0. .AND. (TZERO(I).LT.173.16 .OR. 
     1                           TZERO(I).GT.373.16))               THEN 
              IBAD=I
          ENDIF  
 200  CONTINUE
C
      IF(IBAD.NE.0)                                                 THEN
          WRITE(6,6275) IBAD,JL,TZERO(IBAD),NITER(IBAD),ISNOW
 6275     FORMAT('0BAD ITERATION TEMPERATURE',3X,2I3,F16.2,2I4)
          WRITE(6,6280) QSWNET(IBAD),QLWIN(IBAD),QSENS(IBAD),
     1        QEVAP(IBAD),GZERO(IBAD),CFLUX(IBAD),RIB(IBAD)
 6280     FORMAT(2X,7F12.4)
          CALL XIT('TSOLVE',-1)
      ENDIF 
      !
      !Finally, a check is performed to ensure that TZERO is not less 
      !than 0 C if ponded water is present on the surface (IWATER = 1) 
      !or greater than 0 C if snow is present on the surface 
      !(IWATER = 2), or greater than zero if the surface is an ice sheet 
      !(ISAND = -4). If any of these cases is true, TZERO is reset to 
      !the freezing point, and QZERO and T(0)v are recalculated. DRCOEF 
      !or FLXSURFZ are called, and their calculations are performed for 
      !all locations meeting these criteria. The components of the 
      !surface energy balance are recalculated; the residual amount is 
      !assigned to the energy associated with phase change of water at 
      !the surface, QMELT, and RESID is set to zero. If QMELT < 0, i.e. 
      !if there is freezing taking place, the evaporation flux QEVAP is 
      !added to it and then reset to zero (since if ponded water is 
      !freezing, it is unavailable for evaporation).
C
C     * POST-ITERATION CLEAN-UP. 
C
      NIT=0
      DO 300 I=IL1,IL2
          IF(FI(I).GT.0.)                                          THEN
              IF(((IWATER(I).EQ.1 .AND. TZERO(I).LT.TFREZ) .OR. 
     1            (IWATER(I).EQ.2 .AND. TZERO(I).GT.TFREZ)) .OR.
     2            (ISAND(I,1).EQ.-4 .AND. TZERO(I).GT.TFREZ))   THEN 
                  TZERO(I)=TFREZ        
                  WZERO(I)=0.622*611.0/PADRY(I)
                  QZERO(I)=WZERO(I)/(1.0+WZERO(I))    
                  TVIRTS(I)=TZERO(I)*(1.0+0.61*QZERO(I))
                  ITER(I)=1
                  NIT=NIT+1 
              ELSE
                  ITER(I)=0
              ENDIF  
          ENDIF
  300 CONTINUE
C
      IF(NIT.GT.0)                                                  THEN 
C
C       * CALCULATE SURFACE DRAG COEFFICIENTS (STABILITY-DEPENDENT) AND
C       * OTHER RELATED QUANTITIES.
C
        IF(ISLFD.LT.2) THEN
            CALL DRCOEF (CDM,CDH,RIB,CFLUX,QZERO,QA,ZOSCLM,ZOSCLH,
     1                   CRIB,TVIRTS,TVIRTA,VA,FI,ITER,
     2                   ILG,IL1,IL2)
        ELSE
            CALL FLXSURFZ(CDM,CDH,CFLUX,RIB,FTEMP,FVAP,ILMO,
     1                    UE,FCOR,TPOTA,QA,ZRSLFM,ZRSLFH,VA,
     2                    TZERO,QZERO,H,ZOM,ZOH,
     3                    LZZ0,LZZ0T,FM,FH,ILG,IL1,IL2,FI,ITER,JL )
        ENDIF
      ENDIF
      !
      !In the last half of the loop, some final adjustments are made to 
      !a few variables. If the evaporation flux is vanishingly small, it 
      !is added to RESID and reset to zero. If an anomalous case has 
      !arisen in which QMELT < 0 over a snow-covered surface or 
      !QMELT > 0 over a snow-free surface, QMELT is added to the heat 
      !flux into the surface and then reset to zero. Any remaining 
      !residual flux is added to QH. The shortwave radiation transmitted 
      !into the surface is added back to the net shortwave radiation for 
      !diagnostic purposes. The surface vapour flux is converted into 
      !units of m s-1. Lastly, the iteration counter ITERCT is updated 
      !for the level corresponding to the subarea type and the value of 
      !NITER.
C
C     * REMAINING CALCULATIONS.
C
      DO 350 I=IL1,IL2 
          IF(FI(I).GT.0. .AND. ITER(I).EQ.1)                       THEN
              QLWOUT(I)=SBC*TZERO(I)*TZERO(I)*TZERO(I)*TZERO(I)
              IF(TZERO(I).LT.TPOTA(I))                        THEN
                  QSENS(I)=(RHOAIR(I)*SPHAIR*CFLUX(I)+EZERO)*(TZERO(I)-
     1                TPOTA(I))
              ELSE
                  QSENS(I)=RHOAIR(I)*SPHAIR*CFLUX(I)*(TZERO(I)-
     1                TPOTA(I))
              ENDIF
              EVAP(I)=RHOAIR(I)*CFLUX(I)*(QZERO(I)-QA(I)) 
              QEVAP(I)=CPHCH(I)*EVAP(I)       
              GZERO(I)=GCOEFF(I)*TZERO(I)+GCONST(I)
              QMELT(I)=QSWNET(I)+QLWIN(I)-QLWOUT(I)-QSENS(I)-QEVAP(I)-
     1                 GZERO(I)                             
              RESID(I)=0.0
              IF(QMELT(I).LT.0.0) THEN
                  QMELT(I)=QMELT(I)+QEVAP(I)
                  QEVAP(I)=0.0
                  EVAP(I) =0.0
              ENDIF
          ENDIF                              
C
          IF(FI(I).GT.0.)                                 THEN
              IF(ABS(EVAP(I)).LT.1.0E-8) THEN
                  RESID(I)=RESID(I)+QEVAP(I)
                  EVAP(I)=0.0
                  QEVAP(I)=0.0
              ENDIF
              IF((ISNOW.EQ.1 .AND. QMELT(I).LT.0.0) .OR.
     1            (ISNOW.EQ.0 .AND. QMELT(I).GT.0.0))     THEN
                  GZERO(I)=GZERO(I)+QMELT(I)
                  QMELT(I)=0.0
              ENDIF
C              QSENS(I)=QSENS(I)+0.5*RESID(I)
C              GZERO(I)=GZERO(I)+0.5*RESID(I)
              QSENS(I)=QSENS(I)+RESID(I)
              QSWNET(I)=QSWNET(I)+QTRANS(I)
              EVAP(I)=EVAP(I)/RHOW
              ITERCT(I,KF(I),NITER(I))=ITERCT(I,KF(I),NITER(I))+1
          ENDIF
  350 CONTINUE
C
      RETURN                                                                      
      END  
