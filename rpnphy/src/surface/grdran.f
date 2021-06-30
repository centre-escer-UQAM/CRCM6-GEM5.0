      SUBROUTINE GRDRAN(IVEG,THLIQ,THICE,TBARW,FDT,TFDT,BASFLW,TBASFL,
     1                  RUNOFF,TRUNOF,QFG,WLOST,FI,EVAP,R,ZPOND,DT,
     2                  WEXCES,THLMAX,THTEST,THPOR,THLRET,THLMIN,
     3                  BI,PSISAT,GRKSAT,THFC,DELZW,XDRAIN,ISAND,LZF,
     4                  IGRN,IGRD,IGDR,IG,IGP1,IGP2,ILG,IL1,IL2,JL,N,
     5                  WT,ZBOTW,QDIS,QIN,Isat,WTnew,DMLIQ,
     6                  igrinfl,KQIN,Ibed,DMLIQJ,Isats,xsi,ZWT,EXCW,
     7                  ARE,LEG,SLP,SANI,igwscheme)
C
C     Purpose: Quantify movement of liquid water between soil layers 
C     under non-infiltrating conditions, in response to gravity and 
C     tension forces.
C
c     * JUL 06/16 - A.GANJI.    ADD REVISED FROZEN SOIL SCHEME
c                               (IMPLEMENTED BY L.DUARTE)
C     * OCT 18/11 - M.LAZARE.   PASS IN "IGDR" AS AN INPUT FIELD 
C     *                         (ORIGINATING IN CLASSB) RATHER
C     *                         THAN REPEATING THE CALCULATION HERE
C     *                         AS AN INTERNAL WORK FIELD.
C     * SEP 11/11 - D.VERSEGHY. CHANGE IF CONDITION ON "J.LT.IG"
C     *                         TO "J.LT.IGDR(I)" IN LOOPS 400
C     *                         AND 500, TO BE CONSISTENT WITH
C     *                         OTHERS.
C     * DEC 10/10 - D.VERSEGHY. ALLOW DRAINAGE AT BEDROCK SURFACE
C     *                         ANYWHERE IN SOIL PROFILE.
C     * DEC 23/09 - V.FORTIN.   NEW CALCULATION OF BASEFLOW.
C     * MAR 31/09 - D.VERSEGHY. PASS IN LZF, AND ZERO OUT FLOWS AT
C     *                         TOP AND BOTTOM OF SOIL LAYERS DOWN
C     *                         TO LAYER CONTAINING WETTING FRONT 
C     *                         IN CASES WHERE INFILTRATION IS
C     *                         OCCURRING.
C     * JAN 06/09 - D.VERSEGHY. MODIFIED CALCULATION OF GRKSATF;
C     *                         ADJUSTMENTS TO WATER FLUX 
C     *                         CORRECTIONS IN 500 LOOP.
C     * MAR 27/08 - D.VERSEGHY. MOVE VISCOSITY ADJUSTMENT TO WPREP.
C     * OCT 31/06 - R.SOULIS.   ADJUST GRKSAT FOR VISCOSITY OF
C     *                         WATER AND PRESENCE OF ICE; ADJUST
C     *                         THPOR FOR PRESENCE OF ICE.
C     * JUN 06/06 - F.SEGLENIEKS. CHANGE CALCULATION OF GRSBND
C     *                           TO USE HARMONIC MEAN INSTEAD OF
C     *                           GEOMETRIC MEAN.
C     * MAY 17/06 - D.VERSEGHY. MODIFY CALCULATION OF THLMAX TO
C     *                         ALLOW FOR OVERSATURATED CONDITIONS.
C     * MAR 21/06 - D.VERSEGHY. PROTECT CALCULATIONS OF TBASFL AND
C     *                         TRUNOF AGAINST DIVISION BY ZERO.
C     * MAR 23/05 - R.SOULIS/D.VERSEGHY. CALCULATE GRKSAT AND
C     *                         PSISAT AT LAYER BOUNDARIES USING 
C     *                         GEOMETRIC MEAN; SET BASEFLOW TO
C     *                         GRKSAT IF THLIQ > THFC; ADD
C     *                         CALCULATION OF RUNOFF TEMPERATURE.
C     * MAR 16/05 - D.VERSEGHY. TREAT FROZEN SOIL WATER AS ICE
C     *                         VOLUME RATHER THAN AS EQUIVALENT
C     *                         LIQUID WATER VOLUME.
C     * SEP 24/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.
C     * JUL 27/04 - D.VERSEGHY,Y.DELAGE. PROTECT SENSITIVE 
C     *                         CALCULATIONS AGAINST ROUNDOFF ERRORS.
C     * DEC 03/03 - D.VERSEGHY. IMPROVE HANDLING OF CAPILLARY RISE
C     *                         VS. GRAVITY DRAINAGE (ESPECIALLY
C     *                         FOR ORGANIC SOILS).
C     * JUL 31/03 - D.VERSEGHY. ALWAYS CALCULATE THLMAX IN 100 LOOP.
C     * OCT 23/02 - D.VERSEGHY. REFINEMENT OF TEST IN 400 LOOP.
C     * JUN 21/02 - D.VERSEGHY. BUGFIX IN CALCULATION OF FDT'S IN
C     *                         400 LOOP; UPDATE SUBROUTINE CALL;
C     *                         SHORTENED CLASS4 COMMON BLOCK.
C     * MAY 21/02 - D.VERSEGHY. STREAMLINE CALCULATIONS FOR ORGANIC
C     *                         SOILS AND MODIFY CHECK ON EVAPORATION
C     *                         RATE.
C     * DEC 12/01 - D.VERSEGHY. ADD SEPARATE CALCULATION OF BASEFLOW
C     *                         AT BOTTOM OF SOIL COLUMN.
C     * SEP 28/00 - P.BARTLETT/D.VERSEGHY. BUG FIX IN CALCULATION
C     *                                    OF PSI IN LOOP 200.
C     * FEB 08/00 - D.VERSEGHY/L.SPACEK. MINOR BUG FIX IN LOOP 600
C     *                                  RE. ADDRESSING OF THLIQ.
C     * JUN 20/97 - D.VERSEGHY. CLASS - VERSION 2.7.
C     *                         MODIFICATIONS TO ALLOW FOR VARIABLE
C     *                         SOIL PERMEABLE DEPTH.
C     * DEC 30/96 - D.VERSEGHY. CLASS - VERSION 2.6.
C     *                         BUGFIX IN CALCULATION OF QFG.
C     * AUG 30/95 - D.VERSEGHY. CLASS - VERSION 2.4.
C     *                         ADDITIONAL DIAGNOSTIC CALCULATIONS.
C     * AUG 18/95 - D.VERSEGHY. REVISIONS TO ALLOW FOR INHOMOGENEITY
C     *                         BETWEEN SOIL LAYERS.
C     * APR 24/92 - D.VERSEGHY/M.LAZARE. CLASS - VERSION 2.1.
C     *                                  REVISED AND VECTORIZED CODE
C     *                                  FOR MODEL VERSION GCM7.
C     * AUG 12/91 - D.VERSEGHY. CODE FOR MODEL VERSION GCM7U -
C     *                         CLASS VERSION 2.0 (WITH CANOPY).
C     * APR 11/89 - D.VERSEGHY. UPDATE SOIL LAYER TEMPERATURES AND
C     *                         LIQUID MOISTURE CONTENTS FOR
C     *                         NON-INFILTRATING CONDITIONS (I.E.
C     *                         NO PONDED WATER AND NO RAINFALL 
C     *                         OCCURRING WITHIN CURRENT TIMESTEP).
C
      IMPLICIT NONE
C
C     * INTEGER CONSTANTS.
C
      integer ifrsoil
      INTEGER IVEG  !Subarea type flag  
      INTEGER IG,IGP1,IGP2,ILG,IL1,IL2,JL,I,J,K,IPTBAD,N

      integer IWT(ILG),KQIN(ILG),Ibed(ILG), Isat(ILG), Isats(ILG)
      INTEGER igrinfl
      real ZBOTW (ILG,IG),DMLIQ(ILG),DMLIQJ(ILG), xsi(ILG,IG), EXCW(ILG)
      real WT(ILG), QDIS(ILG), QIN(ILG), WTnew(ILG), ZWT(ILG)
      real QDIS2, QIN2
      real SANI(ILG), ARE(ILG), SLP(ILG),    LEG(ILG)
      integer igwscheme
      real l_k(IG), limit
C
C     * INPUT/OUTPUT FIELDS.
C
      REAL THLIQ (ILG,IG)   !Volumetric liquid water content of soil 
                            !layer [m3 m-3] (theta_l)
      REAL THICE (ILG,IG)   !Volumetric frozen water content of soil 
                            !layer [m3 m-3] (theta_i)
      REAL TBARW (ILG,IG)   !Temperature of water in soil layer [C]
      REAL FDT  (ILG,IGP1)  !Water flow at soil layer interfaces during 
                            !current time step [m]
      REAL TFDT  (ILG,IGP1) !Temperature of water flowing between soil 
                            !layers [C]
C                        
      REAL BASFLW(ILG)  !Base flow from bottom of soil column [m]    
      REAL TBASFL (ILG) !Temperature of base flow from bottom of soil 
                        !column [K]
      REAL RUNOFF(ILG)  !Total runoff from soil column [m]
      REAL TRUNOF (ILG) !Temperature of total runoff from soil column 
                        ![K]
      REAL QFG    (ILG) !Evaporation from soil surface (diagnostic) 
                        ![kg m-2 s-1]
      REAL WLOST (ILG)  !Residual amount of water that cannot be 
                        !supplied by surface stores [kg m-2]
C
C     * INPUT FIELDS.
C
      REAL FI    (ILG)  !Fractional coverage of subarea in question on 
                        !modelled area [ ] (Xi)
      REAL EVAP  (ILG)  !Evaporation rate from ground surface [m s-1]  
      REAL R     (ILG)  !Rainfall rate at ground surface [m s-1]  
      REAL ZPOND (ILG)  !Depth of ponded water on soil surface [m]  
      REAL DT    (ILG)  !Time period over which water movement takes 
                        !place [s]
C
      INTEGER              IGRN  (ILG)  !Flag to indicate whether 
                                        !calculations in subroutine 
                                        !GRINFL are done
      INTEGER              LZF   (ILG)  !Index of soil layer in which 
                                        !wetting front is located
      INTEGER              IGRD  (ILG)  !Flag to indicate whether 
                                        !calculations in this subroutine 
                                        !are to be done
      INTEGER              IGDR  (ILG)  !Index of soil layer in which 
                                        !bedrock is encountered
C
C     * WORK FIELDS.
C
      REAL WEXCES(ILG),    THLMAX(ILG,IG), THTEST(ILG,IG),
     1     GRKSATF(ILG,IG),THPORF(ILG,IG) 
C
C     * SOIL INFORMATION ARRAYS.
C
      REAL THPOR (ILG,IG)   !Pore volume in soil layer [m3 m-3] 
                            !(theta_p) 
      REAL THLRET(ILG,IG)   !Liquid water retention capacity for organic 
                            !soil [m3 m-3 ] (theta_l,ret)
      REAL THLMIN(ILG,IG)   !Residual soil liquid water content 
                            !remaining after freezing or evaporation 
                            ![m3 m-3] (theta_l,min)
      REAL BI    (ILG,IG)   !Clapp and Hornberger empirical “b” 
                            !parameter [ ] (b)
      REAL PSISAT(ILG,IG)   !Soil moisture suction at saturation [m] 
                            !(psi_sat)
      REAL GRKSAT(ILG,IG)   !Hydraulic conductivity of soil at 
                            !saturation [m s-1] (Ksat)
      REAL THFC  (ILG,IG)   !Field capacity [m3 m-3]
      REAL DELZW (ILG,IG)   !Permeable depth of soil layer [m] 
                            !(delta_zg,w)
      REAL XDRAIN(ILG)      !Drainage index for water flow at bottom of 
                            !soil profile [ ]
C  
      INTEGER              ISAND (ILG,IG)   !Sand content flag
C
C     * TEMPORARY VARIABLES.
C
      real alpha, balpha
      REAL THPBND,THLBND,DTHLDZ,DTHPDZ,BBND,GRSBND,PSSBND,GRK,PSI,
     1     WLIMIT,THSUBL,THLTHR,CCH,ASAT,ASATC,SATB
C
C     * COMMON BLOCK PARAMETERS.
C
      REAL DELT     !Time step [s]
      REAL TFREZ    !Freezing point of water [K]
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
C
      COMMON /CLASS1/ DELT,TFREZ
      COMMON /CLASS4/ HCPW,HCPICE,HCPSOL,HCPOM,HCPSND,HCPCLY,
     1                SPHW,SPHICE,SPHVEG,SPHAIR,RHOW,RHOICE,
     2                TCGLAC,CLHMLT,CLHVAP
      common /switches/ ifrsoil
C-----------------------------------------------------------------------
C     * DETERMINE POINTS WHICH SATISFY CONDITIONS FOR THESE CALCULATIONS
C     * AND STORE THEM AS HAVING NON-ZERO VALUES FOR WORK ARRAY "IGRD".
C     * NOTE THAT POINTS WHICH GO THROUGH THE ROUTINE "GRINFL" SHOULD
C     * NOT GO THROUGH THIS ROUTINE WHEN IT IS CALLED FROM CLASSW.
C     * THE INPUT ARRAY "IGRN" HANDLES THIS CONDITION (PASSED AS
C     * "IZERO" ARRAY WHEN CALLED FROM "WEND" OR THE END OF "GRINFL"). 
C
      !In loop 50, the flag IGRD is first set to 1 for all grid cells 
      !where the calculations in this subroutine are to be performed. 
      !The necessary conditions are: that the surface being modelled is 
      !not a glacier or ice sheet (ISAND > -4); that the time period DT 
      !is greater than zero; that the infiltration calculations in 
      !subroutine GRINFL are not simultaneously being performed 
      !(IGRN = 0); and that the rainfall rate and the depth of ponded 
      !water are both vanishingly small. If any of these conditions is 
      !not met, IGRD is set to zero.
      !

      QIN2=0.0
      if (igwscheme.eq.1) then
        do i=il1,il2
          EXCW(I)=0.0
        enddo
      endif
      DO 50 I=IL1,IL2
          IF(FI (I).GT.0. .AND. 
     1       ISAND(I,1).GT.-4 .AND.DT(I).GT.0. .AND.IGRN(I).EQ.0 .AND.
     2       (R(I).LT.1.0E-12 .AND. ZPOND(I).LT.1.0E-12))     THEN
              IGRD(I)=1
          ELSE
              IGRD(I)=0
          ENDIF
   50 CONTINUE
C
C     * CALCULATE MAXIMUM LIQUID WATER CONTENT OF EACH SOIL LAYER;
C     * ADJUST GRKSAT FOR VISCOSITY OF WATER AND PRESENCE OF ICE;
C     * ADJUST THPOR FOR PRESENCE OF ICE.
C
      !In loop 100, if the surface being modelled is not an ice sheet 
      !(ISAND = -4) and if the soil layer in question is not completely 
      !rock (ISAND = -3), the maximum possible water content of each 
      !soil layer, THLMAX, is calculated as the maximum of the available 
      !pore volume THPOR – THICE (where THPOR is the total pore volume 
      !and THICE is the ice content of the layer), the actual liquid 
      !water content THLIQ, and the minimum residual liquid water 
      !content THLMIN. The last two conditions are required because in 
      !the case of a saturated soil undergoing freezing, since water 
      !expands when frozen, the sum of the liquid and frozen volumetric 
      !water contents may be greater than the pore volume, and thus 
      !THLIQ or THLMIN may be greater than THPOR – THICE. An effective 
      !saturated hydraulic conductivity GRKSATF of the soil layer is 
      !also defined, applying an empirical correction for the presence 
      !of ice. This ice content factor, fice, is calculated from Zhao 
      !and Gray (1997) as:
      !
      !fice = [1.0 – min((THPOR - THLMIN)/THPOR, THICE/THPOR )]^2
      !
      !The pore volume of the soil layer is corrected for the presence 
      !of ice by defining the effective porevolume THPORF as equivalent 
      !to THLMAX.
      !
      if (ifrsoil.eq.2) then
        alpha=3.0
        balpha=8.0
      else
        alpha=0.0
        balpha=0.0
      endif
      DO 100 J=1,IG
      DO 100 I=IL1,IL2
        IF(IGRD(I).GT.0)                             THEN
          IF(ISAND(I,J).GT.-3)             THEN
              THLMAX(I,J)=MAX((THPOR(I,J)-THICE(I,J)-0.00001),
     1            THLIQ(I,J),THLMIN(I,J))
              if (ifrsoil.eq.2) then
                  GRKSATF(I,J)=GRKSAT(I,J)/(1.0+balpha*THICE(I,J))**4.0
     1                     *MIN(1.0,(1.0-MAX(0.0,
     2                     (EXP(-alpha*(MAX(0.0,1.0-THICE(I,J)/
     3                     THPOR(I,J))))-EXP(-alpha))))) 
              else
                  GRKSATF(I,J)=GRKSAT(I,J)*(1.0-MAX(0.0,MIN((THPOR(I,J)-
     1               THLMIN(I,J))/THPOR(I,J),THICE(I,J)/THPOR(I,J))))**2
              endif
              THPORF(I,J)=THLMAX(I,J)
          ELSE
              THLMAX(I,J)=0.0
              GRKSATF(I,J)=0.0
              THPORF(I,J)=0.0
          ENDIF
        ENDIF
  100 CONTINUE
C
C     * CALCULATE THEORETICAL FLOW RATES AT BOTTOM OF PERMEABLE SOIL
C     * DEPTH AND BETWEEN SOIL LAYERS.
C
      !In loops 150 and 200, the theoretical amounts of water FDT 
      !flowing across the soil layer boundaries are calculated. FDT is 
      !evaluated as the product of the flow rate F(z) at the given depth 
      !z, and the period of time (DT) over which the flow is occurring. 
      !At z=0, the water flux is simply equal to the soil surface 
      !evaporation rate. Within the soil, the flow rate F at a given 
      !depth z is obtained using equation 21 from Verseghy (1991):
      !
      !F(z) = K(z)*[-b*psi(z)/THLIQ(z) · d(THLIQ)/dz + 1]
      !
      !where K(z) is the hydraulic conductivity and psi(z) is the soil 
      !moisture suction at depth z, and b is an empirical parameter 
      !developed by Clapp and Hornberger (1978). K(z) and ψ(z) are 
      !calculated following Clapp and Hornberger as
      !
      !K(z) = GRKSAT*(THLIQ/THPOR)^(2b + 3)
      !psi(z) = PSISAT*(THLIQ/THPOR)^(-b )
      !
      !where GRKSAT and PSISAT are the values of K and psi respectively 
      !at saturation.
      !
      !At the bottom of the permeable soil depth, in layer IGDR, if the 
      !liquid water content of the soil is greater than the field 
      !capacity, the vertical flow out of the bottom of the soil profile 
      !is calculated using a relation derived from Soulis et al. (2010):
      !
      !F(zb) = GRKSAT · min{1, (THLIQ/THPOR)/[1 – 1/(2b + 3)]}^(2b + 3)
      !
      !This flow rate is multiplied by a drainage parameter XDRAIN, 
      !which is set to 0 if the soil is underlain by an impermeable 
      !layer (as in a bog), and to 1 otherwise.
      !
      DO 150 I=IL1,IL2
          IF(IGRD(I).GT.0)                                          THEN
             FDT(I,1)=-EVAP(I)*DT(I)                                                           
             IF(DELZW(I,IGDR(I)).GT.0.0001)                      THEN
               if (igwscheme.ne.1) then
                 DMLIQ(I)=0.0
                 EXCW(I)=0.0

                 CALL  GW(THLIQ,THICE,DT,THPOR,
     1           BI,PSISAT,GRKSAT,DELZW,
     2           IGRD,IG,IGP1,IGP2,ILG,IL1,IL2,JL,N,IGDR,
     3           ZBOTW,THPORF,I,GRKSATF,THLMIN,THLMAX,IWT,QDIS2,
     4           QIN2,WT,WTnew,DMLIQ,igrinfl,KQIN,Ibed,DMLIQJ,
     5           Isats,xsi,ZWT,TBARW,EXCW,ARE,LEG,SLP,SANI,igwscheme)
                 EXCW(I) = EXCW(I)*DT(I)/1000.0
                 IF (IWT(I).LE.IGDR(I)+1) THEN
                   IGDR(I)=Isat(I)    
                 ENDIF               
                 FDT(I,IGDR(I)+1)=QDIS2*DT(I)/1000.0
                 QDIS(I)=QDIS2/1000.0        ! (m/s)
                 QIN2=QIN2/1000.0
                 QIN(I)=QIN2
               else
                 IF(THLIQ(I,IGDR(I)).GT.THFC(I,IGDR(I)))      THEN
                     CCH=2.0*BI(I,IGDR(I))+3.0
                     ASATC=1.0-(1.0/CCH)
                     ASAT=THLIQ(I,IGDR(I))/THPORF(I,IGDR(I))
                     SATB=MIN(1.0,ASAT/ASATC)
                     FDT(I,IGDR(I)+1)=GRKSATF(I,IGDR(I))*DT(I)*
     1                   XDRAIN(I)*SATB**CCH
                 ELSE                                                                        
                     FDT(I,IGDR(I)+1)=0.0                                                           
                 ENDIF
               endif
             ELSE
               FDT(I,IGDR(I)+1)=0.0
             ENDIF
          ELSE
             if (igwscheme.ne.1) then
               IWT(I)=0
               QDIS2=0.0
             endif
          ENDIF
          if (igwscheme.ne.1 .and. QIN(I).NE.0.0.AND.
     &        IGRD(I).EQ.0.AND.DT(I).NE.0.0) QIN2=QIN(I)
  150 CONTINUE
C
      !
      !Between soil layers, values of K(z) and psi(z) must be 
      !determined. This requires the estimation of soil properties at 
      !the layer interfaces. For THLIQ, THPOR and d(THLIQ)/dz, slightly 
      !different approaches are followed if the permeable soil layer 
      !thickness DELZW increases with depth (the normal case), or if it 
      !decreases between one soil layer and the next (indicating the 
      !presence of an impermeable barrier at some depth in the second 
      !layer). In the first case, the pore volume THPBND, liquid water 
      !content THLBND, and liquid water gradient DTHLDZ are simply 
      !calculated as arithmetic averages of the values above and below 
      !the interface. This has the effect of weighting the interface 
      !values towards the upper layer values, roughly emulating a 
      !surfaceward exponential decay curve of liquid water content. In 
      !the second case, in order to avoid a spurious weighting towards 
      !the lower layer, the gradients of liquid water content and pore 
      !volume between the upper and lower layers are calculated as 
      !linear relations, and then solved for the values at the 
      !interface.
      !
      !The Clapp and Hornberger b parameter at the interface is 
      !calculated as a simple average of the values in the upper and 
      !lower soil layers. The saturated hydraulic conductivity GRSBND is 
      !calculated as a harmonic mean, and the saturated soil moisture 
      !suction PSSBND as a geometric mean, of the top and bottom 
      !layer values:
      !
      !GRSBND = GRKSAT,t*GRKSAT,b*(DELZW,t + DELZW,b)/
      !                           (GRKSAT,t DELZW,b+ GRKSAT,b DELZW,t )
      !
      !PSSBND = PSISAT,t^[DELZW,t/(DELZW,t + DELZW,b)]
      !         *PSISAT,b^[DELZW,b/(DELZW,t + DELZW,b)]
      !
      !Finally, K(z), psi(z) and F(z) at the interface are calculated 
      !from the equations given above. At the end of the loop, if 
      !infiltration is occurring, the flow rates between soil layers are 
      !set to zero behind the soil layer LZF containing the wetting 
      !front, and the flow rate at the bottom of LZF is constrained to 
      !be >= 0.
      !
      DO 200 J=1,IG-1                                                             
      DO 200 I=IL1,IL2
        if (igwscheme.eq.1 .or. j.le.Isat(i)) then
          IF(IGRD(I).GT.0)                                     THEN
            IF(J.LT.IGDR(I))                    THEN
              IF(THPOR(I,J).GT.0.0.AND.THPOR(I,J+1).GT.0.0.AND.
     1                  ISAND(I,J+1).GT.-3)            THEN
                  IF(DELZW(I,J+1).GT.DELZW(I,J)) THEN
                      THPBND=(THPORF(I,J)+THPORF(I,J+1))/2.0
                      THLBND=(THLIQ(I,J)+THLIQ(I,J+1))/2.0                                        
                      DTHLDZ=(THLIQ(I,J+1)-THLBND)/DELZW(I,J+1)+
     1                       (THLBND-THLIQ(I,J))/DELZW(I,J)
                  ELSE
                      DTHLDZ=2.0*(THLIQ(I,J+1)-THLIQ(I,J))/
     1                       (DELZW(I,J+1)+DELZW(I,J))
                      THLBND=THLIQ(I,J)+0.5*DTHLDZ*DELZW(I,J)
                      DTHPDZ=2.0*(THPORF(I,J+1)-THPORF(I,J))/
     1                       (DELZW(I,J+1)+DELZW(I,J))
                      THPBND=THPORF(I,J)+0.5*DTHPDZ*DELZW(I,J)
                  ENDIF
                  BBND=(BI(I,J)+BI(I,J+1))/2.0
C                  GRSBND=GRKSAT(I,J)**(DELZW(I,J)/(DELZW(I,J)+
C     1                DELZW(I,J+1)))*GRKSAT(I,J+1)**(DELZW(I,J+1)/
C     2                (DELZW(I,J)+DELZW(I,J+1)))
                  GRSBND=GRKSATF(I,J)*GRKSATF(I,J+1)*(DELZW(I,J)+
     1                DELZW(I,J+1))/(GRKSATF(I,J)*DELZW(I,J+1)+
     2                GRKSATF(I,J+1)*DELZW(I,J))
                  PSSBND=PSISAT(I,J)**(DELZW(I,J)/(DELZW(I,J)+
     1                DELZW(I,J+1)))*PSISAT(I,J+1)**(DELZW(I,J+1)/
     2                (DELZW(I,J)+DELZW(I,J+1)))
                  GRK=MIN(GRSBND*(THLBND/THPBND)**(2.*BBND+3.),
     1                   GRSBND)                     
                  PSI=MAX(PSSBND*(THLBND/THPBND)**(-BBND),PSSBND)
     1             *(1.0+balpha*THICE(I,J))**2.0
                  FDT(I,J+1)=GRK*DT(I)*((-BBND*PSI*DTHLDZ/THLBND)+1.)
              ELSE
                  FDT(I,J+1)=0.0
              ENDIF
              IF(ABS(THLIQ(I,J)-THLIQ(I,J+1)).LT.0.05 .AND. 
     1            FDT(I,J).LT.0.0) FDT(I,J+1)=0.0
              IF(LZF(I).GT.0 .AND. J.LT.LZF(I)) FDT(I,J+1)=0.0
              IF(LZF(I).GT.0 .AND. J.EQ.LZF(I) .AND. FDT(I,J+1)
     1                 .LT.0.0) FDT(I,J+1)=0.0
            ENDIF
          ENDIF
        endif
  200 CONTINUE 
C                               
C     * CHECK FOR SUSTAINABLE EVAPORATION RATE FROM TOP SOIL LAYER; IF
C     * LIQUID WATER SUPPLY IS INSUFFICIENT, TRY TO REMOVE WATER FROM 
C     * FROZEN SOIL MOISTURE.
C
      !In the next several loops, checks are carried out to ascertain 
      !whether the theoretically determined flow rates at the soil layer 
      !interfaces can be supported by the water contents in the layers. 
      !At the beginning of loop 250, the soil surface evaporation rate 
      !is addressed. A trial liquid water content THTEST is calculated 
      !for the first soil layer, reflecting the removal of the 
      !evaporated water. If THTEST is less than the minimum soil water 
      !content, F(0) is set to zero and THLIQ is set to THLMIN. The 
      !excess surface flux that was not met by the removed liquid water 
      !is converted to a frozen water content, and an attempt is made to 
      !remove it from the frozen water in the layer. (The energy 
      !involved in converting the required frozen water mass to liquid 
      !water is obtained by decreasing the water temperature of the 
      !layer.) The frozen water content of the layer is adjusted to 
      !reflect the removal of the required amount. If the demand exceeds 
      !the supply of frozen water available, the frozen water content is 
      !set to zero, the diagnostic evaporative flux QFG at the soil 
      !surface is adjusted, and the remaining water flux not able to be 
      !met by the first soil layer is assigned to the variable WLOST.
      !
      !In the remainder of the 250 loop, checks are carried out for soil 
      !layers with liquid water contents already effectively equal to 
      !THLMIN. If the flux at the top of the layer is upward and that at 
      !the bottom of the layer is downward, both are set to zero. If 
      !both are downward, only the flux at the bottom is set to zero; if 
      !both are upward, only the flux at the top is set to zero. If the 
      !soil layer is an organic one and the liquid water content is less 
      !than the layer retention capacity THLRET and the flux at the 
      !bottom of the layer is downward, it is set to zero.
      !
      IPTBAD=0                                        
      DO 250 J=1,IG                                                               
      DO 250 I=IL1,IL2
        if (igwscheme.eq.1 .or. j.le.Isat(i)) then
          IF(IGRD(I).GT.0 .AND. J.EQ.1 .AND.
     1                          DELZW(I,J).GT.0.0)    THEN
              !Huziy bring inside the if on FDT
            IF(FDT(I,J) .LT. 0) THEN
              
              THTEST(I,J)=THLIQ(I,J)+FDT(I,J)/DELZW(I,J)
              IF(THTEST(I,J).LT.THLMIN(I,J))             THEN
                  FDT(I,J)=FDT(I,J)+(THLIQ(I,J)-THLMIN(I,J))*
     1                DELZW(I,J)
                  THLIQ(I,J)=THLMIN(I,J)
                  WEXCES(I)=-FDT(I,J)                                                  
                  FDT(I,J)=0.0                                                      
                  THSUBL=WEXCES(I)*RHOW/(RHOICE*DELZW(I,J))
                  IF(THLIQ(I,J).GT.0.0)                 THEN
                      TBARW(I,J)=TBARW(I,J)-(CLHMLT*RHOICE*THSUBL)/ 
     1                           (HCPW*THLIQ(I,J)) 
                  ENDIF
                  IF(THSUBL.LE.THICE(I,J))              THEN
                      THICE(I,J)=THICE(I,J)-THSUBL
                  ELSE
                      THSUBL=THSUBL-THICE(I,J)
                      THICE(I,J)=0.0
                      QFG(I)=QFG(I)-FI(I)*THSUBL*RHOICE*DELZW(I,J)/
     1                       DELT
                      WLOST(I)=WLOST(I)+THSUBL*RHOICE*DELZW(I,J)
                  ENDIF
              ENDIF
              IF(THICE(I,J).LT.0.) IPTBAD=I
            ENDIF
          ENDIF
C                                                    
C     * ENSURE THAT CALCULATED WATER FLOWS BETWEEN SOIL LAYERS DO NOT
C     * CAUSE LIQUID MOISTURE CONTENT OF ANY LAYER TO FALL BELOW THE
C     * RESIDUAL VALUE OR TO EXCEED THE CALCULATED MAXIMUM.
C
          IF(IGRD(I).GT.0)                                        THEN
              if (igwscheme.ne.1) then
                limit=-0.0001
              else
                limit=0.001
              endif
c              IF(THLIQ(I,J).LE.(THLMIN(I,J)+0.001)
              IF(THLIQ(I,J).LE.(THLMIN(I,J)+limit)
     1            .AND. J.LE.IGDR(I))                            THEN
               if (igwscheme.ne.1 .and. J.EQ.IGDR(I)) then
                IF(KQIN(I).EQ.1)                                    THEN
                  IF(FDT(I,J).LE.0. .AND.(QIN(I)*DT(I)-DMLIQ(I)).GE.0.)
     1                                                           THEN
                    FDT(I,J)=0.0
c$$$                  ELSEIF(FDT(I,J).GE.0. .AND.
c$$$     1                  (QIN(I)*DT(I)-DMLIQ(I)).GT.0.) THEN
                  ELSEIF(FDT(I,J).LT.0. .AND.
     1                  (QIN(I)*DT(I)-DMLIQ(I)).LE.0.) THEN
                    FDT(I,J)=0.0
                  ENDIF
                ELSE
                  IF(FDT(I,J).LE.0. .AND. FDT(I,J+1).GE.0.)      THEN                          
                    FDT(I,J)=0.0
                  ELSE IF(FDT(I,J).GE.0. .AND. FDT(I,J+1).GT.0.) THEN                      
                  ELSE IF(FDT(I,J).LT.0. .AND. FDT(I,J+1).LE.0.) THEN                      
                    FDT(I,J)=0.0
                  ENDIF
                ENDIF
               else
                IF(FDT(I,J).LE.0. .AND. FDT(I,J+1).GE.0.)      THEN                          
                    FDT(I,J)=0.0    
                    FDT(I,J+1)=0.0   
                ELSE IF(FDT(I,J).GE.0. .AND. FDT(I,J+1).GT.0.) THEN                      
                    FDT(I,J+1)=0.0 
                ELSE IF(FDT(I,J).LT.0. .AND. FDT(I,J+1).LE.0.) THEN                      
                    FDT(I,J)=0.0    
                ENDIF
               endif
              ENDIF
          ENDIF   
          IF(IGRD(I).GT.0 .AND. ISAND(I,J).EQ.-2 .AND. 
     1                           THLIQ(I,J).LE.THLRET(I,J))    THEN
              IF(FDT(I,J+1).GT.0.0) FDT(I,J+1)=0.0
          ENDIF
        endif
  250 CONTINUE    
C
      IF(IPTBAD.NE.0)                                           THEN
          WRITE(6,6500) IPTBAD,JL,IVEG,THICE(IPTBAD,1)
 6500     FORMAT('0AT (I,J)=(',I3,',',I3,'), IVEG=',I2,' THICE(1)= ',
     1            E13.5)
          CALL XIT('GRDRAN',-1)
      ENDIF
C
      !In the 300 loop, checks are carried out for soil layers with 
      !liquid water contents already effectively equal to THLMAX. If the 
      !flux at the top of the layer is downward and that at the bottom 
      !of the layer is upward, both are set to zero. If both are 
      !downward, then if the top flux is greater than the bottom flux, 
      !it is set equal to the bottom flux. If both are upward, then if 
      !magnitude of the bottom flux is greater than that of the top 
      !flux, it is set equal to the top flux.
      !
      DO 300 J=IG,1,-1                                                            
      DO 300 I=IL1,IL2
c          IF(IGRD(I).GT.0)                                          THEN
          IF(J.LE.Isat(I).AND.IGRD(I).GT.0)  THEN
            if (igwscheme.ne.1) then
              limit=0.00001
            else
              limit=0.001
            endif
            IF(igwscheme.ne.1 .and. THLIQ(I,J).GE.(THLMAX(I,J)-limit) 
     1                      .AND. J.EQ.IGDR(I))                     THEN  
                  IF(KQIN(I).EQ.1. AND. FDT(I,J).GE.0. .AND. (QIN(I)
     1                   *DT(I)-DMLIQ(I)).LE.0.)                    THEN                     
                         FDT(I,J)=0.0                                                      
                  ELSE IF(KQIN(I).EQ.1. AND. FDT(I,J).GT.0. .AND.
     1                   (QIN(I)*DT(I)-DMLIQ(I)).GE.0.)             THEN                     
                         IF(FDT(I,J).GT.(QIN(I)*DT(I)-DMLIQ(I)))
     1                      FDT(I,J)=QIN(I)*DT(I)-DMLIQ(I)                                               
                  ELSEIF(KQIN(I).EQ.0. AND. FDT(I,J).GT.0. .AND.
     1                   FDT(I,J+1).GE.0.)                         THEN
                         IF(FDT(I,J).GT.FDT(I,J+1)) FDT(I,J)=FDT(I,J+1)
                  ENDIF
            ELSEIF(THLIQ(I,J).GE.(THLMAX(I,J)-limit)
     1                        .AND. J.LE.IGDR(I))                THEN
              IF(FDT(I,J).GE.0. .AND. FDT(I,J+1).LE.0.)      THEN                          
                  FDT(I,J)=0.0                                                      
                  FDT(I,J+1)=0.0                                                    
              ELSE IF(FDT(I,J).GT.0. .AND. FDT(I,J+1).GE.0.) THEN                      
                  IF(FDT(I,J).GT.FDT(I,J+1)) FDT(I,J)=FDT(I,J+1)                          
              ELSE IF(FDT(I,J).LE.0. .AND. FDT(I,J+1).LT.0.) THEN                      
                  IF(FDT(I,J+1).LT.FDT(I,J)) FDT(I,J+1)=FDT(I,J)                          
              ENDIF                                                               
            ENDIF                                                               
          ENDIF                                                                   
  300 CONTINUE
C
      !In the 400 loop, for each soil layer THTEST is recalculated as 
      !the liquid water content resulting from the updated boundary 
      !fluxes, and is compared with a residual value THLTHR. For organic 
      !soil layers deeper than the first layer, THLTHR is set to the 
      !minimum of THLRET and THLIQ, since only evapotranspiration can 
      !cause the soil moisture to fall below the retention capacity. 
      !(The first layer is excepted because surface evaporation can 
      !drive its moisture content below THLRET.) For mineral soils, 
      !THLTHR is set to THLMIN. If THTEST < THLTHR, then if the flow at 
      !the bottom of the soil layer is downward, it is recalculated as 
      !the sum of the flow at the top of the layer plus the amount of 
      !water that must be removed to make the liquid water content of 
      !the layer equal to THLTHR. If the flow at the bottom of the layer 
      !is upward, the flow at the top of the layer is recalculated as 
      !the sum of the flow at the bottom of the layer minus the amount 
      !of water that must be removed to make the liquid water content of 
      !the layer equal to THLTHR. THTEST of the current layer is then 
      !reset to THLTHR, and THTEST of the overlying and underlying 
      !layers (if any) are recalculated using the new interface fluxes.
      !
      DO 400 J=1,IG                                                               
      DO 400 I=IL1,IL2
       if (igwscheme.eq.1 .or. J.LE.Isat(I)) then
        IF(IGRD(I).GT.0)                                           THEN
          IF(J.LE.IGDR(I) .AND. ISAND(I,J).NE.-3)              THEN
            if (igwscheme.ne.1) then
              IF (J.EQ.IGDR(I).AND.KQIN(I).EQ.1) THEN 
                 THTEST(I,J)=THLIQ(I,J)+FDT(I,J)/DELZW(I,J) 
              ELSEIF (J.EQ.IGDR(I).AND.KQIN(I).EQ.0) THEN
                 THTEST(I,J)=THLIQ(I,J)+FDT(I,J)/DELZW(I,J)  
              ENDIF
              IF (J.LT.IGDR(I)) THEN
                 THTEST(I,J)=THLIQ(I,J)+(FDT(I,J)-FDT(I,J+1))/DELZW(I,J)
              ENDIF
            else
              THTEST(I,J)=THLIQ(I,J)+(FDT(I,J)-FDT(I,J+1))/DELZW(I,J)
            endif
              IF(ISAND(I,J) .EQ. -2 .AND. J .NE. 1)      THEN
                  THLTHR=MIN(THLRET(I,J),THLIQ(I,J))
              ELSE
                  THLTHR=THLMIN(I,J)
              ENDIF
              IF(THTEST(I,J).LT.THLTHR)                  THEN
                  IF(FDT(I,J+1).GT.0. .and. (igwscheme.eq.1 .or.
     &            J.LT.IGDR(I))) THEN                      
                      FDT(I,J+1)=FDT(I,J)+(THLIQ(I,J)-THLTHR)*
     1                    DELZW(I,J)  
                  ELSE
                    IF (igwscheme.eq.1 .or. KQIN(I).EQ.1) THEN
                      FDT(I,J)=FDT(I,J+1)-(THLIQ(I,J)-THLTHR)*
     1                    DELZW(I,J)
                    ELSEIF(igwscheme.ne.1.and.
     &                     KQIN(I).EQ.0.AND.J.LT.IWT(I)) THEN
                      FDT(I,J)=FDT(I,J+1)-(THLIQ(I,J)-THLTHR)*DELZW(I,J)
                    ELSEIF(igwscheme.ne.1.and.
     &                     KQIN(I).EQ.0.AND.J.GE.IWT(I)) THEN
                      FDT(I,J)=FDT(I,J)-(THLIQ(I,J)-THLTHR)*DELZW(I,J)
                    ENDIF
                  ENDIF
                  THTEST(I,J)=THLTHR
                  IF(J.LT.IGDR(I)) THEN
                      IF(DELZW(I,J+1).GT.0.0) THTEST(I,J+1)=THLIQ(I,J+1)
     1                    +(FDT(I,J+1)-FDT(I,J+2))/DELZW(I,J+1)
                  ENDIF
                  IF(J.GT.1) THEN
                      IF(DELZW(I,J-1).GT.0.0) THTEST(I,J-1)=THLIQ(I,J-1)
     1                    +(FDT(I,J-1)-FDT(I,J))/DELZW(I,J-1)
                  ENDIF
              ENDIF                                                                   
          ELSE
              THTEST(I,J)=0.0
          ENDIF
        ENDIF
       endif
  400 CONTINUE               
C
      !In the 500 loop, for each soil layer THTEST is compared with 
      !THLMAX. If THTEST > THLMAX, two temporary variables are defined: 
      !WLIMIT as the amount of water required to raise the liquid water 
      !content of the soil layer to THLMAX, and WEXCES as the amount of 
      !water representing the excess of THTEST over THLMAX. If the flux 
      !at the top of the layer is downward and the flux at the bottom of 
      !the layer is upward, then if the negative of the flux at the 
      !bottom of the layer is greater than WLIMIT, it is set equal to 
      !-WLIMIT and the flux at the top is set to zero; otherwise WEXCES 
      !is subtracted from the flux at the top. If both the flux at the 
      !top and the flux at the bottom are downward, then WEXCES is 
      !subtracted from the flux at the top. If both are upward, then 
      !WEXCES is added to the flux at the bottom, and a correction is 
      !applied to the flux at the bottom of the underlying layer. 
      !Finally, THTEST is recalculated for all the soil layers using the 
      !new interface fluxes.
      !
      DO 500 J=IG,1,-1
      DO 500 I=IL1,IL2
        IF(igwscheme.eq.1.or.J.LE.Isat(I))  THEN
          IF(IGRD(I).GT.0)  THEN
            IF(THTEST(I,J).GT.THLMAX(I,J) .AND. J.LE.IGDR(I))    THEN
              WLIMIT=MAX((THLMAX(I,J)-THLIQ(I,J)),0.0)*DELZW(I,J)                      
              WEXCES(I)=(THTEST(I,J)-THLMAX(I,J))*DELZW(I,J)                                
              IF(igwscheme.ne.1.and.J.EQ.IGDR(I)) THEN 
               IF(KQIN(I).EQ.1)   THEN
                IF(FDT(I,J).GT.0. .AND. 
     1            (QIN(I)*DT(I)-DMLIQ(I)).LE.0.) THEN
                    FDT(I,J)=FDT(I,J)-WEXCES(I)
                ELSE IF(FDT(I,J).GT.0. .AND. 
     1                 (QIN(I)*DT(I)-DMLIQ(I)).GE.0.) THEN                      
                    FDT(I,J)=FDT(I,J)-WEXCES(I)                                            
                ELSE IF(FDT(I,J).LE.0. .AND. 
     1                 (QIN(I)*DT(I)-DMLIQ(I)).LT.0.) THEN                      
                ENDIF
               ELSE
                IF(FDT(I,J).GT.0. .AND. FDT(I,J+1).EQ.0.)        THEN                          
                   FDT(I,J)=FDT(I,J)-WEXCES(I)
                ELSE IF(FDT(I,J).GT.0. .AND. FDT(I,J+1).GE.0.)   THEN                      
                  FDT(I,J)=FDT(I,J)-WEXCES(I)                                            
                ELSE IF(FDT(I,J).LE.0. .AND. FDT(I,J+1).LT.0.)   THEN                      
                ENDIF
               ENDIF
              ELSE
               IF(FDT(I,J).GT.0. .AND. FDT(I,J+1).LE.0.)        THEN                          
                IF(-FDT(I,J+1).GT.WLIMIT)          THEN
                    FDT(I,J+1)=-WLIMIT
                    FDT(I,J)=0.0
                ELSE
                    FDT(I,J)=FDT(I,J)-WEXCES(I)
                ENDIF
C                IF(FDT(I,J).GE.WLIMIT)             THEN                                       
C                   FDT(I,J)=WLIMIT                                               
C                   FDT(I,J+1)=0.0                                                
C                ELSE                                                            
C                   FDT(I,J+1)=FDT(I,J)-WLIMIT                                      
C                ENDIF                                                           
              ELSE IF(FDT(I,J).GT.0. .AND. FDT(I,J+1).GE.0.)   THEN                      
                FDT(I,J)=FDT(I,J)-WEXCES(I)                                            
              ELSE IF(FDT(I,J).LE.0. .AND. FDT(I,J+1).LT.0.)   THEN                      
                FDT(I,J+1)=FDT(I,J+1)+WEXCES(I)                                        
c                IF(J.LT.IGDR(I))                       THEN
                IF((igwscheme.eq.1.and.J.LT.IGDR(I)) .or.
     &            (igwscheme.ne.1.and.J.LT.IG.AND.(J+2).NE.IGDR(I)))THEN
                    IF(FDT(I,J+2).LT.0.) FDT(I,J+2)=0.0
                ENDIF
               ENDIF
              ENDIF
              DO 450 K=1,IG
                  IF(DELZW(I,K).GT.0.0)                            THEN
                   IF(igwscheme.ne.1.and.K.EQ.IGDR(I))  THEN
                      THTEST(I,IGDR(I))=THLIQ(I,IGDR(I))+(
     1                          FDT(I,IGDR(I)))/DELZW(I,IGDR(I))
                   ELSE
                      THTEST(I,K)=THLIQ(I,K)+(FDT(I,K)-FDT(I,K+1))/
     1                            DELZW(I,K)
                   ENDIF
                  ENDIF
  450         CONTINUE
            ENDIF                                                                   
          ENDIF                                                                   
        endif
  500 CONTINUE
C
      !In the first part of the 600 loop, a correction is performed in 
      !the unlikely event that as a result of the above adjustments, the 
      !flux at the bottom of the last permeable soil layer has become 
      !upward. If this occurs, all of the fluxes above are corrected by 
      !the inverse of this same amount. However, this will result in a 
      !small spurious upward water flux at the surface. Since the liquid 
      !water flows are now accounted for, an attempt is made, as in loop 
      !250, to remove the required water from the frozen water store in 
      !the first layer. The excess that cannot be supplied from this 
      !source is assigned to WLOST.
      !
      !In the second part of the 600 loop, the temperatures TFDT of the 
      !water fluxes at the top and bottom of the layers in the permeable 
      !soil profile are set equal to the temperatures of the water in 
      !the first and last layers respectively. The flow at the bottom of 
      !the profile and the temperature of the flow are used to update 
      !BASFLW and TBASFL, the baseflow from the current subarea and its 
      !temperature, and RUNOFF and TRUNOF, the total runoff from the 
      !grid cell in question and its temperature.
      !
      IPTBAD=0
      DO 600 I=IL1,IL2
          IF(IGRD(I).GT.0)                                        THEN
              IF(FDT(I,IGDR(I)+1).LT.0.)                     THEN
                  WEXCES(I)=-FDT(I,IGDR(I)+1)
                  DO 550 J=1,IGDR(I)+1
                    IF(J.LE.Isat(I)+1)  THEN
                      FDT(I,J)=FDT(I,J)+WEXCES(I)
                    ENDIF
  550             CONTINUE
                  THSUBL=WEXCES(I)*RHOW/(RHOICE*DELZW(I,1))                                     
                  IF(THLIQ(I,1).GT.0.0)               THEN
                      TBARW(I,1)=TBARW(I,1)-(CLHMLT*RHOICE*THSUBL)/
     1                           (HCPW*THLIQ(I,1))                
                  ENDIF
                  IF(THSUBL.LE.THICE(I,1))            THEN
                      THICE(I,1)=THICE(I,1)-THSUBL                                        
                  ELSE
                      THSUBL=THSUBL-THICE(I,1)
                      THICE(I,1)=0.0
                      QFG(I)=QFG(I)-FI(I)*THSUBL*RHOICE*DELZW(I,1)/
     1                       DELT
                      WLOST(I)=WLOST(I)+THSUBL*RHOICE*DELZW(I,1)
                  ENDIF
                  IF(THICE(I,1).LT.0.0) IPTBAD=I
              ENDIF                                                                       
C
C     * CALCULATE DRAINAGE FROM BOTTOM OF SOIL COLUMN AND RE-EVALUATE
C     * SOIL LAYER TEMPERATURES AND LIQUID MOISTURE CONTENTS AFTER
C     * WATER MOVEMENT.
C
              TFDT(I,1)=TBARW(I,1) 
              TFDT(I,IGDR(I)+1)=TBARW(I,IGDR(I))
c              IF(FDT(I,IGDR(I)+1).GT.0.0) THEN
              IF(FDT(I,IGDR(I)+1).GT.0.0 .OR. EXCW(I).GT.0.0) THEN
                  IF((BASFLW(I)+FI(I)*FDT(I,IGDR(I)+1)).GT.0.0) 
     1              TBASFL(I)=(TBASFL(I)*BASFLW(I)+(TFDT(I,IGDR(I)+1)+
     2                TFREZ)*FI(I)*FDT(I,IGDR(I)+1))/(BASFLW(I)+FI(I)*
     3                FDT(I,IGDR(I)+1))
                  BASFLW(I)=BASFLW(I)+FI(I)*FDT(I,IGDR(I)+1)
                  IF((RUNOFF(I)+FDT(I,IGDR(I)+1)).GT.0.0) 
c     1              TRUNOF(I)=(TRUNOF(I)*RUNOFF(I)+(TFDT(I,IGDR(I)+1)+
     1              TRUNOF(I)=(TRUNOF(I)*(RUNOFF(I)+EXCW(I))+
     2                        (TFDT(I,IGDR(I)+1)+
     2                TFREZ)*FDT(I,IGDR(I)+1))/(RUNOFF(I)+
     3                FDT(I,IGDR(I)+1)+EXCW(I))
                  RUNOFF(I)=RUNOFF(I)+FDT(I,IGDR(I)+1)+EXCW(I)
              ENDIF
          ENDIF
  600 CONTINUE
C
      IF(IPTBAD.NE.0)                                           THEN
          WRITE(6,6500) IPTBAD,JL,IVEG,THICE(IPTBAD,1)
          CALL XIT('GRDRAN',-3)
      ENDIF
C     
      !
      !In the 700 loop, for each successive soil layer the temperature 
      !of the water flux at the bottom of the layer is first assigned. 
      !If the flux is downward, TFDT is set to the temperature of the 
      !water in the current layer; if it is upward, TFDT is set to the 
      !temperature of the water in the layer below. The temperature of 
      !the water in the current layer is then updated using the FDT and 
      !TFDT values at the top and bottom of the layer, and the liquid 
      !water content of the layer is set to THTEST.
      !                                                 
      DO 700 J=1,IG
      DO 700 I=IL1,IL2
        IF(J.LE.Isat(I))  THEN
          IF(IGRD(I).GT.0)                                      THEN
            IF(J.LE.IGDR(I))                               THEN
              IF(J.LT.IGDR(I))                       THEN
                IF(FDT(I,J+1).GT.0.)          THEN                                                
                   TFDT(I,J+1)=TBARW(I,J)                                                  
                ELSE                                                                    
                   TFDT(I,J+1)=TBARW(I,J+1)                                                
                ENDIF                                                                   
              ENDIF
              IF(THTEST(I,J).GT.0.0 .AND. DELZW(I,J).GT.0.0)    THEN
                IF(igwscheme.ne.1.and.J.EQ.IGDR(I))  THEN
                 TBARW(I,J)=(THLIQ(I,J)*TBARW(I,J)+(FDT(I,J)*TFDT(I,J))
     1                      /DELZW(I,J))/THTEST(I,J)
                ELSE
                 TBARW(I,J)=(THLIQ(I,J)*TBARW(I,J)+(FDT(I,J)*TFDT(I,J)-
     1                      FDT(I,J+1)*TFDT(I,J+1))/DELZW(I,J))/
     2                      THTEST(I,J) 
                ENDIF                                      
              ENDIF
              THLIQ(I,J)=THTEST(I,J)
            ENDIF                                                      
          ENDIF
        endif
  700 CONTINUE                                                                    
C                                                                                  
      RETURN                                                                      
      END        
