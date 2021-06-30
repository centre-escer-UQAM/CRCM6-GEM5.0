      SUBROUTINE GRINFL(IVEG,THLIQ,THICE,TBARW,BASFLW,TBASFL,
     1                  RUNOFF,TRUNOF,ZFAV,LZFAV,THLINV,QFG,
     2                  WLOST,FI,EVAP,R,TR,TPOND,ZPOND,DT,
     3                  ZMAT,WMOVE,TMOVE,THLIQX,THICEX,TBARWX,
     4                  DELZX,ZBOTX,FDT,TFDT,PSIF,THLINF,GRKINF,
     5                  THLMAX,THTEST,ZRMDR,FDUMMY,TDUMMY,THLDUM,
     6                  THIDUM,TDUMW,TRMDR,ZF,FMAX,TUSED,RDUMMY,
     7                  ZERO,WEXCES,FDTBND,WADD,TADD,WADJ,TIMPND,
     8                  DZF,DTFLOW,THLNLZ,THLQLZ,DZDISP,WDISP,WABS,
     9                  THPOR,THLRET,THLMIN,BI,PSISAT,GRKSAT,
     A                  THLRAT,THFC,DELZW,ZBOTW,XDRAIN,DELZ,ISAND,
     B                  IGRN,IGRD,IFILL,IZERO,LZF,NINF,IFIND,ITER,
     C                  NEND,ISIMP,IGDR,
     D                  IG,IGP1,IGP2,ILG,IL1,IL2,JL,N,WT,
     d                  QDIS,QIN,Isat,WTnew,DMLIQ,KQIN,Ibed,DMLIQJ,
     e                  Isats,xsi,ZWT,EXCW,ARE,LEG,SLP,SANI,igwscheme)
C
C     Purpose: Quantify movement of liquid water between soil layers 
C     under conditions of infiltration.
C
c     * JUL 06/16 - A.GANJI.    ADD REVISED FROZEN SOIL SCHEME
c                               (IMPLEMENTED BY L.DUARTE)
C     * OCT 18/11 - M.LAZARE.   PASS IN "IGDR" AS AN INPUT FIELD 
C     *                         (ORIGINATING IN CLASSB) TO
C     *                         GRDRAN AND WEND.
C     * APR 04/11 - D.VERSEGHY. MODIFY TEST IN 150 LOOP TO USE DELZ
C     *                         INSTEAD OF XDRAIN.
C     * JAN 06/09 - D.VERSEGHY. MODIFY DELZX AND ZBOTX OF BOTTOM LAYER;
C     *                         ADDITIONAL THLIQ CHECK IN 350 LOOP;
C     *                         PASS ADDITIONAL VARIABLES TO WEND.
C     * MAR 27/08 - D.VERSEGHY. MOVE VISCOSITY ADJUSTMENT TO WPREP.
C     * OCT 31/06 - R.SOULIS.   ADJUST GRKSAT FOR VISCOSITY OF WATER
C     *                         AND PRESENCE OF ICE; ADJUST THPOR FOR
C     *                         PRESENCE OF ICE.
C     * MAR 22/06 - D.VERSEGHY. UNCONDITIONALLY DEFINE VARIABLES FOR
C     *                         ALL "IF" STATEMENTS.
C     * SEP 28/05 - D.VERSEGHY. REMOVE HARD CODING OF IG=3 IN 400 LOOP.
C     * MAR 23/05 - D.VERSEGHY.R.SOULIS. PASS ADDITIONAL VARIABLES 
C     *                         TO GRDRAN AND WEND; PASS OUT ZFAV, 
C     *                         LZFAV, THLINV; CALCULATE GRKTLD 
C     *                         INTERNALLY; REVISE CALCULATION OF 
C     *                         THLINF.
C     * MAR 16/04 - D.VERSEGHY. TREAT FROZEN SOIL WATER AS ICE
C     *                         VOLUME RATHER THAN AS EQUIVALENT
C     *                         LIQUID WATER VOLUME.
C     * SEP 24/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.
C     * JUL 27/04 - Y.DELAGE/D.VERSEGHY. PROTECT SENSITIVE
C     *                         CALCULATIONS AGAINST ROUNDOFF ERRORS.
C     * JUL 26/02 - D.VERSEGHY. SHORTENED CLASS4 COMMON BLOCK.
C     * DEC 12/01 - D.VERSEGHY. PASS NEW VARIABLE IN FOR CALCULATION
C     *                         OF BASEFLOW.
C     * JUN 20/97 - D.VERSEGHY. CLASS - VERSION 2.7.
C     *                         MODIFICATIONS TO ALLOW FOR VARIABLE
C     *                         SOIL PERMEABLE DEPTH.
C     * APR 17/96 - D.VERSEGHY. CLASS - VERSION 2.5.
C     *                         BUG FIX: INITIALIZE FDT AND TFDT
C     *                         TO ZERO.
C     * AUG 18/95 - D.VERSEGHY. CLASS - VERSION 2.4.
C     *                         REVISIONS TO ALLOW FOR INHOMOGENEITY
C     *                         BETWEEN SOIL LAYERS.
C     * APR 24/92 - D.VERSEGHY/M.LAZARE. CLASS - VERSION 2.1.
C     *                                  REVISED AND VECTORIZED CODE 
C     *                                  FOR MODEL VERSION GCM7.
C     * AUG 12/91 - D.VERSEGHY. CODE FOR MODEL VERSION GCM7U -
C     *                         CLASS VERSION 2.0 (WITH CANOPY).            
C     * APR 11/89 - D.VERSEGHY. UPDATE SOIL LAYER TEMPERATURES AND 
C     *                         LIQUID MOISTURE CONTENTS FOR 
C     *                         INFILTRATING CONDITIONS (I.E.
C     *                         PONDED WATER OR RAINFALL OCCURRING
C     *                         WITHIN CURRENT TIMESTEP).
C
      IMPLICIT NONE
C
C     * INTEGER CONSTANTS.
C
      integer ifrsoil
      INTEGER IVEG  !Subarea type flag
      INTEGER IG,IGP1,IGP2,ILG,IL1,IL2,JL,I,J,N
      integer KQIN(ILG),Ibed(ILG), Isat(ILG), Isats(ILG)
      real WT(ILG),WTnew(ILG),DMLIQ(ILG),QINold(ILG),DMLIQold(ILG),
     &     QDISold(ILG), DMLIQJ(ILG),ZWT(ILG)
      real QDIS(ILG), QIN(ILG),xsi(ILG,IG),EXCW(ILG),EXCWold(ILG),
     &     SANI(ILG), ARE(ILG), SLP(ILG),    LEG(ILG)
      integer igwscheme
C  
C     * INPUT/OUTPUT FIELDS.
C
      REAL THLIQ (ILG,IG)   !Volumetric liquid water content of soil 
                            !layer [m3 m-3] (theta_l)
      REAL THICE (ILG,IG)   !Volumetric frozen water content of soil 
                            !layer [m3 m-3] (theta_i)
      REAL TBARW (ILG,IG)   !Temperature of water in soil layer [C]
C                        
      REAL BASFLW(ILG)  !Base flow from bottom of soil column [m] 
      REAL TBASFL(ILG)  !Temperature of base flow from bottom of 
                        !soil column [K]  
      REAL RUNOFF(ILG)  !Total runoff from soil column [m]  
      REAL TRUNOF(ILG)  !Temperature of total runoff from soil column 
                        ![K]
      REAL QFG   (ILG)  !Evaporation from soil surface (diagnostic) 
                        ![kg m-2 s-1]
      REAL WLOST (ILG)  !Residual amount of water that cannot be 
                        !supplied by surface stores [kg m-2]
      REAL ZFAV  (ILG)  !Average depth of wetting front over current 
                        !time step [m]
      REAL THLINV(ILG)  !Liquid water content behind the wetting front 
                        ![m3 m-3]
C
      INTEGER              LZFAV (ILG)  !Soil layer index in which the 
                                        !average position of the wetting 
                                        !front occurs
C
C     * INPUT FIELDS.
C
      REAL FI    (ILG)  !Fractional coverage of subarea in question on 
                        !modelled area [ ]
      REAL EVAP  (ILG)  !Evaporation rate from ground surface [m s-1]  
      REAL R     (ILG)  !Rainfall rate at ground surface [m s-1]  
      REAL TR    (ILG)  !Temperature of rainfall [C]
      REAL TPOND (ILG)  !Temperature of ponded water [C]  
      REAL ZPOND (ILG)  !Depth of ponded water on soil surface [m]  
      REAL DT    (ILG)  !Time period over which water movement takes 
                        !place [s]
C
C     * WORK FIELDS (FOR ALL CALLED ROUTINES AS WELL).
C
      REAL ZMAT  (ILG,IGP2,IGP1)
C
      REAL WMOVE (ILG,IGP2),   TMOVE (ILG,IGP2),
     1     GRKSATF(ILG,IG),    THPORF(ILG,IG)
C
      REAL THLIQX(ILG,IGP1),   THICEX(ILG,IGP1),   TBARWX(ILG,IGP1),
     1     DELZX (ILG,IGP1),   ZBOTX (ILG,IGP1),   FDT   (ILG,IGP1),
     2     TFDT  (ILG,IGP1),   PSIF  (ILG,IGP1),   THLINF(ILG,IGP1),   
     3     GRKINF(ILG,IGP1),   THLMAX(ILG,IG),     THTEST(ILG,IG),     
     4     ZRMDR (ILG,IGP1),   FDUMMY(ILG,IGP1),   TDUMMY(ILG,IGP1),
     5     THLDUM(ILG,IG),     THIDUM(ILG,IG),     TDUMW (ILG,IG)
C
      REAL TRMDR (ILG),    ZF    (ILG),    FMAX  (ILG),    TUSED (ILG),
     1     RDUMMY(ILG),    ZERO  (ILG),    WEXCES(ILG),    FDTBND(ILG),    
     2     WADD  (ILG),    TADD  (ILG),    WADJ  (ILG),    TIMPND(ILG),    
     3     DZF   (ILG),    DTFLOW(ILG),    THLNLZ(ILG),    THLQLZ(ILG),    
     4     DZDISP(ILG),    WDISP (ILG),    WABS  (ILG)  
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
      REAL THLRAT(ILG,IG)   !Fractional saturation of soil behind the 
                            !wetting front [ ] (finf)
      REAL THFC  (ILG,IG)   !Field capacity [m3 m-3]
      REAL DELZW (ILG,IG)   !Permeable depth of soil layer [m] 
                            !(delta_zg,w)
      REAL ZBOTW (ILG,IG)   !Depth to permeable bottom of soil layer [m]
      REAL XDRAIN(ILG)      !Drainage index for water flow at bottom of 
                            !soil profile [ ]
      REAL DELZ(IG)         !Overall thickness of soil layer [m]
C
C     * TEMPORARY VARIABLES.
C
      real alpha, balpha
      REAL PSIINF,GRK,PSI
C
C     * VARIOUS INTEGER ARRAYS.
C
      INTEGER              ISAND (ILG,IG)   !Sand content flag 
      INTEGER              IGRN  (ILG)  !Flag to indicate whether 
                                        !calculations in this subroutine 
                                        !are to be done
      INTEGER              IGRD  (ILG)  !Flag to indicate whether 
                                        !calculations in subroutine 
                                        !GRDRAN are to be done
      INTEGER              IFILL (ILG),    IZERO (ILG),    LZF   (ILG),
     1                     NINF  (ILG),    IFIND (ILG),    ITER  (ILG),
     2                     NEND  (ILG),    ISIMP (ILG)    
      INTEGER              IGDR  (ILG)  !Index of soil layer in which 
                                        !bedrock is encountered
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
C     * AND STORE THEM AS HAVING NON-ZERO VALUES FOR WORK ARRAY "IGRN".
C
      !
      !In loop 50, the flag IGRN is first set to 1 for all grid cells 
      !where the calculations in this subroutine are to be performed. 
      !The necessary conditions are: that the surface being modelled is 
      !not a glacier or ice sheet (ISAND > -4); that the time period DT 
      !is greater than zero; and that either the rainfall rate is 
      !greater than zero or that ponded water exists on the surface. If 
      !any of these conditions is not met, IGRN is set to zero.
      !
      DO 50 I=IL1,IL2
          IF(FI(I).GT.0. .AND. 
     1       ISAND(I,1).GT.-4 .AND. DT(I).GT.0. .AND.
     2       (R(I).GT.0. .OR. ZPOND(I).GT.0.))                     THEN
              IGRN(I)=1
              RDUMMY(I)=0.
          ELSE
              IGRN(I)=0
          ENDIF
   50 CONTINUE
C
C     * ADJUST GRKSAT FOR VISCOSITY OF WATER AND PRESENCE OF ICE;
C     * ADJUST THPOR FOR PRESENCE OF ICE.
C     * INITIALIZATION; DETERMINATION OF SOIL HYDRAULIC CONDUCTIVITIES
C     * AND SOIL MOISTURE SUCTION ACROSS WETTING FRONT.
C
      !
      !In loop 100, a series of local arrays THLIQX, THICEX, TBARWX, 
      !DELZX and ZBOTX are defined for use in the subsequent 
      !infiltration calculations, with a “y” dimension of IG+1, where IG 
      !is the number of soil layers. The entries from 1 to IG are set 
      !equal to the corresponding values in THLIQ, THICE, TBARW, DELZW 
      !and ZBOTW. An effective saturated hydraulic conductivity GRKSATF 
      !is calculated for each soil layer, applying an empirical 
      !correction for the presence of ice. This ice content factor, 
      !fice, is calculated from Zhao and Gray (1997) as:
      !
      !fice = [1.0 – min(1.0, THICE/THPOR )]^2
      !
      !To further account for the presence of ice, a modified pore 
      !volume THPORF for the soil layer is calculated as the maximum of 
      !the available pore volume THPOR – THICE (where THPOR is the total 
      !pore volume and THICE is the ice content of the layer), the 
      !actual liquid water content THLIQ, and the minimum residual 
      !liquid water content THLMIN. (The last two conditions are 
      !required because in the case of a saturated soil undergoing 
      !freezing, since water expands when frozen, the sum of the liquid 
      !and frozen volumetric water contents may be greater than the pore 
      !volume, and thus THLIQ or THLMIN may be greater than 
      !THPOR – THICE.) Finally, the water content THLINF and the 
      !hydraulic conductivity GRKINF behind the wetting front are 
      !evaluated, following the analysis of Green and Ampt (1911) and 
      !treating the change in soil moisture due to infiltration as a 
      !downward-propagating square wave. THLINF is calculated as the 
      !maximum of finf(THPOR – THICE), THLIQ, and THLMIN, where finf 
      !represents the fractional saturation of the soil behind the 
      !wetting front, corresponding to a hydraulic conductivity of half 
      !the saturated value (this correction is applied in order to 
      !account for the fact that as infiltration occurs, a small amount 
      !of air is usually trapped in the soil). GRKINF is calculated from 
      !GRKSATF, THLINF and THPORF, using the classic Clapp and 
      !Hornberger (1978) equation
      !
      !K(z) = GRKSAT*(THLIQ/THPOR)^(2b + 3)
      !
      !where K(z) is the hydraulic conductivity at a depth z, GRKSAT is 
      !the hydraulic conductivity at saturation, THPOR is the pore 
      !volume and b is an empirical coefficient.
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
          IF(IGRN(I).GT.0)                                         THEN
              THLIQX(I,J)=THLIQ(I,J)                                                      
              THICEX(I,J)=THICE(I,J)                                                      
              TBARWX(I,J)=TBARW(I,J)                                                      
              DELZX(I,J)=DELZW(I,J)                                                        
              ZBOTX(I,J)=ZBOTW(I,J)                                                        
              FDT (I,J)=0.0
              TFDT(I,J)=0.0
              IF(ISAND(I,J).GT.-3)                             THEN
                if (ifrsoil.eq.2) then
                  GRKSATF(I,J)=GRKSAT(I,J)*(1.0-MAX(0.0,MIN(1.0,
     1             THICE(I,J)/THPOR(I,J))))**2/(1.0+balpha*THICE(I,J))
     2                **4.0*MIN(1.0,(1.0-MAX(0.0,
     3                (EXP(-alpha*(MAX(0.0,1.0-THICE(I,J)/
     4                THPOR(I,J))))-EXP(-alpha)))))
                else
                  GRKSATF(I,J)=GRKSAT(I,J)*(1.0-MAX(0.0,MIN(1.0,
     1                THICE(I,J)/THPOR(I,J))))**2
                endif
                  THPORF(I,J)=MAX((THPOR(I,J)-THICE(I,J)-0.00001),
     1                THLIQ(I,J),THLMIN(I,J))                
                  THLINF(I,J)=MAX(THLIQ(I,J),THLMIN(I,J),
     1                        THLRAT(I,J)*(THPOR(I,J)-
     2                        THICE(I,J)-0.00001))
                  GRKINF(I,J)=GRKSATF(I,J)*(THLINF(I,J)/THPORF(I,J))
     1                        **(2.*BI(I,J)+3.)
              ELSE
                  GRKSATF(I,J)=0.0
                  THPORF(I,J)=0.0
                  THLINF(I,J)=0.0
                  GRKINF(I,J)=0.0
              ENDIF
          ENDIF
  100 CONTINUE
C
      !
      !In loop 150, values for the IG+1 entries in each of the above 
      !arrays are assigned. If the permeable thickness DELZW of the IG 
      !layer is less than its overall thickness DELZ (indicating the 
      !presence of bedrock) the IG+1 entries are set to 0; otherwise 
      !they are set to the values for the IG layer in the case of 
      !THLIQX, THICEX, TBARWX and THLINF, and for DELZX and ZBOTX they 
      !are set to large positive numbers. GRKINF is set to the value for 
      !the IG layer multiplied by the drainage parameter XDRAIN, which 
      !takes a value of 0 if the soil is underlain by an impermeable 
      !layer (as in a bog), and a value of 1 otherwise.
      !
      DO 150 I=IL1,IL2
          IF(IGRN(I).GT.0)                                          THEN
              IF(DELZW(I,IG).LT.DELZ(IG))                 THEN
                  THLIQX(I,IG+1)=0.0
                  THICEX(I,IG+1)=0.0
                  TBARWX(I,IG+1)=0.0
                  DELZX(I,IG+1)=0.0
                  THLINF(I,IG+1)=0.0
                  GRKINF(I,IG+1)=0.0
              ELSE
                  THLIQX(I,IG+1)=THLIQX(I,IG)                                                     
                  THICEX(I,IG+1)=THICEX(I,IG)                                                     
                  TBARWX(I,IG+1)=TBARWX(I,IG)                                                     
                  DELZX(I,IG+1)=999999.                                                        
                  THLINF(I,IG+1)=THLINF(I,IG)                                                     
                  GRKINF(I,IG+1)=GRKINF(I,IG)*XDRAIN(I)
              ENDIF
              ZBOTX (I,IG+1)=ZBOTX(I,IG)+DELZX(I,IG+1)
              FDT   (I,IG+1)=0.0
              TFDT  (I,IG+1)=0.0
          ENDIF
  150 CONTINUE
C
      !
      !In loop 200, the soil water suction across the wetting front 
      !psi_f is calculated for each soil layer, using equation 25 in 
      !Verseghy (1991):
      !
      !psi_f = b[psi_inf*Kinf - psi(z)*K(z)]/[Kinf*(b+3)]

      !
      !where psi_inf and Kinf are the soil moisture suction and the 
      !hydraulic conductivity behind the wetting front, and psi(z) and 
      !K(z) are the soil moisture suction and the hydraulic conductivity 
      !ahead of the wetting front. The soil moisture suction values are 
      !obtained using the classic Clapp and Hornberger (1978) equation:
      !
      !psi(z) = PSISAT*(THLIQ/THPPOR)^(-b)
      !
      !where PSISAT is the soil moisture suction at saturation.
      !                                                 
      DO 200 J=1,IG
      DO 200 I=IL1,IL2
          IF(IGRN(I).GT.0)                                          THEN
             IF(THPOR(I,J).GT.0.0001)                     THEN
                 PSIINF=MAX(PSISAT(I,J)*(THLINF(I,J)/THPORF(I,J))**
     1                        (-BI(I,J)),PSISAT(I,J))
     2                        *(1.0+balpha*THICE(I,J))**2.0 
                 GRK=MIN(GRKSATF(I,J)*(THLIQ(I,J)/THPORF(I,J))**
     1                     (2.*BI(I,J)+3.),GRKSATF(I,J))                   
                 PSI=MAX(PSISAT(I,J)*(THLIQ(I,J)/THPORF(I,J))**
     1                     (-BI(I,J)),PSISAT(I,J))
     2                      *(1.0+balpha*THICE(I,J))**2.0
             ELSE
                 PSIINF=PSISAT(I,J)*(1.0+balpha*THICE(I,J))**2.0
                 GRK=GRKSATF(I,J)
                 PSI=PSISAT(I,J)*(1.0+balpha*THICE(I,J))**2.0
             ENDIF
             IF(THLINF(I,J).GT.THLIQ(I,J))                  THEN 
                PSIF(I,J)=MAX(BI(I,J)*(GRKINF(I,J)*PSIINF-GRK*PSI)/
     1                    (GRKINF(I,J)*(BI(I,J)+3.)), 0.0) 
             ELSE                                                                    
                PSIF(I,J)=0.0                                                         
             ENDIF                                                                   
          ENDIF
  200 CONTINUE
C
      DO 250 I=IL1,IL2
          IF(IGRN(I).GT.0)                                          THEN
             PSIF(I,IG+1)=PSIF(I,IG)      
             TRMDR(I)=DELT
          ELSE
             TRMDR(I)=0. 
          ENDIF
  250 CONTINUE
C    
      DO 300 J=1,IGP2
      DO 300 I=IL1,IL2
          IF(IGRN(I).GT.0)                                          THEN 
             WMOVE(I,J)=0.0                 
             TMOVE(I,J)=0.0
          ENDIF
  300 CONTINUE
C
C     * DETERMINE STARTING POSITION OF WETTING FRONT; INITIALIZATION
C     * FOR SATURATED INFILTRATION.
C
      !
      !In loop 400, a test is carried out to determine whether saturated 
      !conditions are already present in the soil. Generally it is 
      !assumed that a period of unsaturated flow occurs initially, so 
      !the flag IFILL is set to 1, the depth of the wetting front ZF is 
      !set to zero, and the flag LZF indicating the index of the soil 
      !layer in which the wetting front occurs is set to 1. If a pond is 
      !present on the soil surface or if the hydraulic conductivity if 
      !the first soil layer is very small, it is deemed that saturated 
      !flow begins immediately, so IFILL is set to zero. Then each soil 
      !layer is checked in turn, to ascertain whether the liquid water 
      !content is greater than THLINF. If so, it is concluded that 
      !saturated flow is occurring in this layer; ZF is set to ZBOTW, 
      !the depth of the bottom of the soil layer, LZF is set to the 
      !index of the next layer, and IFILL is set to zero. For each layer 
      !found to be undergoing saturated flow, its water content is 
      !stored in the J+1 level of the water movement matrix WMOVE, and 
      !its water temperature is stored in the J+1 level of the matrix 
      !TMOVE. The counter NINF is set to the number of soil layers 
      !behind and including the one with the wetting front, plus 2.
      !
      DO 400 I=IL1,IL2
          IF(IGRN(I).GT.0)                                          THEN 
              IFILL(I)=1
              ZF(I)=0.0
              LZF(I)=1
              IF(ZPOND(I).GT.0. .OR. GRKINF(I,1).LT.1.0E-12)   THEN                  
                  NINF(I)=2
                  TMOVE(I,2)=TBARWX(I,1)
                  IFILL(I)=0
              ENDIF
              DO 350 J=1,MIN(IG,Isat(I))
                  IF(THLIQ(I,J).GE.(THLINF(I,J)-1.0E-6) .AND.
     1                    THLIQ(I,J).GT.0.0001 .AND. LZF(I).EQ.J)  THEN
                      ZF(I)=ZBOTW(I,J)
                      LZF(I)=J+1
                      NINF(I)=J+2
                      WMOVE(I,J+1)=THLIQ(I,J)*DELZW(I,J)
                      TMOVE(I,J+1)=TBARWX(I,J)
                      TMOVE(I,J+2)=TBARWX(I,J+1)
                      IFILL(I)=0
                  ENDIF
  350         CONTINUE
          ELSE
              IFILL(I)=0
              LZF(I)=0
              NINF(I)=0
          ENDIF

  400 CONTINUE
      !
      !A series of subroutines is called to complete the infiltration 
      !calculations. Subroutine WFILL treats the period of unsaturated 
      !flow up the point when saturated flow begins. Subroutine WFLOW 
      !treats the period of saturated flow. Subroutine WEND reassigns 
      !liquid water contents and temperatures at the conclusion of the 
      !infiltration period. THLIQ, THICE and TBARW are then set to the 
      !corresponding updated values of THLIQX, THICEX and TBARWX. The 
      !average depth of the wetting front in the soil layer in which it 
      !occurs, the index of the soil layer, and the moisture content 
      !behind the wetting front are stored in output variables ZFAV, 
      !LZFAV and THLINV respectively. Finally, if there is any time 
      !remaining in the current time step following the infiltration 
      !period, subroutine GRDRAN is called to treat the movement of soil 
      !water over the remainder of the time step.
C
C     * IF SATURATED INFILTRATION CONDITIONS ARE NOT PRESENT AT ONCE
C     * (IFILL=1), CALL "WFILL" TO DO PROCESSING FOR PERIOD OF
C     * UNSATURATED INFILTRATION.
C                                                            
      CALL WFILL(WMOVE,TMOVE,LZF,NINF,ZF,TRMDR,R,TR,
     1           PSIF,GRKINF,THLINF,THLIQX,TBARWX,
     2           DELZX,ZBOTX,DZF,TIMPND,WADJ,WADD,
     3           IFILL,IFIND,IGP1,IGP2,ILG,IL1,IL2,Isat )
C
C     * CALL "WFLOW" TO DO PROCESSING FOR PERIOD OF SATURATED
C     * INFILTRATION.
C
      CALL WFLOW(WMOVE,TMOVE,LZF,NINF,TRMDR,TPOND,ZPOND,
     1           R,TR,EVAP,PSIF,GRKINF,THLINF,THLIQX,TBARWX,
     2           DELZX,ZBOTX,FMAX,ZF,DZF,DTFLOW,THLNLZ,
     3           THLQLZ,DZDISP,WDISP,WABS,ITER,NEND,ISIMP,
     4           IGRN,IG,IGP1,IGP2,ILG,IL1,IL2,JL,N,Isat )
C
C     * RECALCULATE TEMPERATURES AND LIQUID MOISTURE CONTENTS OF
C     * SOIL LAYERS FOLLOWING INFILTRATION.
C
      CALL WEND(THLIQX,THICEX,TBARWX,ZPOND,TPOND,
     1          BASFLW,TBASFL,RUNOFF,TRUNOF,FI,
     2          WMOVE,TMOVE,LZF,NINF,TRMDR,THLINF,DELZX,
     3          ZMAT,ZRMDR,FDTBND,WADD,TADD,FDT,TFDT,
     4          THLMAX,THTEST,THLDUM,THIDUM,TDUMW,
     5          TUSED,RDUMMY,ZERO,WEXCES,XDRAIN,
     6          THPOR,THLRET,THLMIN,BI,PSISAT,GRKSAT,
     7          THFC,DELZW,ISAND,IGRN,IGRD,IGDR,IZERO,
     8          IVEG,IG,IGP1,IGP2,ILG,IL1,IL2,JL,N,WT,ZBOTW,
     9          QDIS,QIN,Isat,WTnew,DMLIQ,KQIN,Ibed,DMLIQJ,
     a          Isats,xsi,ZWT,EXCW,ARE,LEG,SLP,SANI,igwscheme)
C
      DO 800 J=1,IG
      DO 800 I=IL1,IL2
          IF(IGRN(I).GT.0)                                          THEN
            IF(J.LT.Isat(I)+1) THEN
              THLIQ(I,J)=THLIQX(I,J)                                                      
              THICE(I,J)=THICEX(I,J)                                                      
              TBARW(I,J)=TBARWX(I,J)      
            ELSE
              THLIQ(I,J)=THPOR(I,J)-THICEX(I,J)                                                      
              THICE(I,J)=THICEX(I,J)                                                      
              TBARW(I,J)=TBARWX(I,J)
            ENDIF     
          ENDIF                                                
  800 CONTINUE
C
      DO 850 I=IL1,IL2
c          IF(IGRN(I).GT.0 .AND. LZF(I).LT.IG+1)                     THEN
          IF(IGRN(I).GT.0 .AND. LZF(I).LT.MIN(IG+1,Isat(I)+1))   THEN
              ZFAV(I)=(ZF(I)+MAX(ZBOTW(I,LZF(I))-DELZW(I,LZF(I)),0.0))/
     1                2.0
              LZFAV(I)=LZF(I)
              THLINV(I)=THLINF(I,LZF(I))
          ELSE
              ZFAV(I)=0.0
              LZFAV(I)=0
              THLINV(I)=0.0
          ENDIF                                                
  850 CONTINUE
C
C     * IF TIME REMAINS IN THE CURRENT MODEL STEP AFTER INFILTRATION
C     * HAS CEASED (TRMDR>0), CALL "GRDRAN" TO CALCULATE WATER FLOWS
C     * BETWEEN LAYERS FOR THE REMAINDER OF THE TIME STEP.
C
      if (igwscheme.ne.1) then
      DO 854 I=IL1,IL2
        IF(TUSED(I).GT.0.0.AND.TRMDR(I).GT.0.0.AND.IGRD(I).EQ.1) THEN
           QINold(I)=QIN(I)
           DMLIQold(I)=DMLIQ(I)
           QDISold(I)=QDIS(I)
           EXCWold(I)=EXCW(I)
        ENDIF
854   CONTINUE
      endif

      CALL GRDRAN(IVEG,THLIQ,THICE,TBARW,FDUMMY,TDUMMY,BASFLW,
     1            TBASFL,RUNOFF,TRUNOF,QFG,WLOST,FI,EVAP,ZERO,ZERO,
     2            TRMDR,WEXCES,THLMAX,THTEST,THPOR,THLRET,THLMIN,
     3            BI,PSISAT,GRKSAT,THFC,DELZW,XDRAIN,ISAND,IZERO,
     4            IZERO,IGRD,IGDR,
     4            IG,IGP1,IGP2,ILG,IL1,IL2,JL,N,WT,ZBOTW,
     5            QDIS,QIN,Isat,WTnew,DMLIQ,1,KQIN,Ibed,DMLIQJ,Isats,
     6            xsi,ZWT,EXCW,ARE,LEG,SLP,SANI,igwscheme)     
C
      if (igwscheme.ne.1) then
       DO 851 J=1,IG
              DO 852 I=IL1,IL2
                 IF(TRMDR(I).GT.0.AND.IGRD(I).EQ.0.AND.
     1           TUSED(I).GT.0.) THEN
                    IF (J.EQ.MIN(IG,Isat(I))) THEN
                       THLIQ(I,J)=THLIQ(I,J)-(QIN(I)*TRMDR(I)*KQIN(I)
     2                 +DMLIQ(I))/DELZW(I,J)
                    ENDIF
                 ENDIF 
852           CONTINUE
851    CONTINUE

       DO 853 I=IL1,IL2
         IF (TRMDR(I).GT.0.AND.TUSED(I).EQ.0.) THEN
             QIN(I)=QIN(I)*TRMDR(I)/1800.
             QDIS(I)=QDIS(I)*TRMDR(I)/1800.
         ENDIF
853    CONTINUE

       DO 855 I=IL1,IL2
          IF(TUSED(I).GT.0.0.AND.TRMDR(I).GT.0.0.AND.IGRD(I).EQ.1) THEN   
             QIN(I)=QIN(I)*TRMDR(I)/1800.+
     1              QINold(I)*(1800.-TRMDR(I))/1800.            
             DMLIQ(I)=DMLIQ(I)+DMLIQold(I)
             QDIS(I)=QDIS(I)*TRMDR(I)/1800.+
     1               QDISold(I)*(1800.-TRMDR(I))/1800.              
             EXCW(I)=EXCW(I)+EXCWold(I)
          ENDIF
855    CONTINUE
      endif
        
      RETURN                                                                      
      END       
