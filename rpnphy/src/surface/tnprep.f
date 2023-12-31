      SUBROUTINE TNPREP(A1,A2,B1,B2,C2,GDENOM,GCOEFF,
     1                  GCONST,CPHCHG,IWATER, 
     2                  TBAR,TCTOP,TCBOT,
     3                  FI,ZPOND,TBAR1P,DELZ,TCSNOW,ZSNOW,
     4                  ISAND,ILG,IL1,IL2,IG                     )
C
C     Purpose: Calculate coefficients for solution of heat conduction 
C     into soil.
C
C     * SEP 05/12 - J.MELTON.   REMOVED UNUSED VARS, JL AND J
C     * MAR 03/08 - D.VERSEGHY. ASSIGN TCTOP3 AND TCBOT3 ON THE BASIS
C     *                         OF SUBAREA VALUES FROM TPREP; REPLACE
C     *                         THREE-LEVEL TEMPERATURE AND THERMAL
C     *                         CONDUCTIVITY VECTORS WITH STANDARD
C     *                         VALUES.
C     * AUG 16/06 - D.VERSEGHY. REMOVE TSTART.
C     * MAY 24/05 - D.VERSEGHY. LIMIT DELZ3 TO <= 4.1 M.
C     * OCT 04/05 - D.VERSEGHY. USE THREE-LAYER TBAR,TCTOP,TCBOT.
C     * NOV 04/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.
C     * AUG 06/02 - D.VERSEGHY. SHORTENED CLASS3 COMMON BLOCK,
C     * JUN 17/02 - D.VERSEGHY. USE NEW LUMPED SOIL AND PONDED
C     *                         WATER TEMPERATURE FOR FIRST LAYER;
C     *                         SHORTENED COMMON BLOCK.
C     * MAR 28/02 - D.VERSEGHY. CHANGE POND THRESHOLD VALUE FOR
C     *                         CALCULATION OF "IWATER".
C     * JUN 20/97 - D.VERSEGHY. CLASS - VERSION 2.7.
C     *                         INCORPORATE EXPLICITLY CALCULATED
C     *                         THERMAL CONDUCTIVITIES AT TOPS AND
C     *                         BOTTOMS OF SOIL LAYERS.
C     * SEP 27/96 - D.VERSEGHY. CLASS - VERSION 2.5.
C     *                         SURFACE TREATED AS WATER ONLY IF
C     *                         ZPOND > 1 MM.
C     * AUG 18/95 - D.VERSEGHY. CLASS - VERSION 2.4.
C     *                         REVISIONS TO ALLOW FOR INHOMOGENEITY
C     *                         BETWEEN SOIL LAYERS AND FRACTIONAL
C     *                         ORGANIC MATTER CONTENT.
C     * NOV 28/94 - M. LAZARE.  CLASS - VERSION 2.3.
C     *                         TCSATW,TCSATI DECLARED REAL(16).
C     * APR 10/92 - M. LAZARE.  CLASS - VERSION 2.1.
C     *                         DIVIDE PREVIOUS SUBROUTINE "T3LAYR"
C     *                         INTO "TNPREP" AND "TNPOST" AND
C     *                         VECTORIZE.
C     * APR 11/89 - D.VERSEGHY. CALCULATE COEFFICIENTS FOR GROUND HEAT
C     *                         FLUX, EXPRESSED AS A LINEAR FUNCTION
C     *                         OF SURFACE TEMPERATURE.  COEFFICIENTS
C     *                         ARE CALCULATED FROM LAYER TEMPERATURES,
C     *                         THICKNESSES AND THERMAL CONDUCTIVITIES,
C     *                         ASSUMING A QUADRATIC VARIATION OF
C     *                         TEMPERATURE WITH DEPTH WITHIN EACH
C     *                         SOIL LAYER. SET THE SURFACE LATENT 
C     *                         HEAT OF VAPORIZATION OF WATER AND 
C     *                         THE STARTING TEMPERATURE FOR THE 
C     *                         ITERATION IN "TSOLVC"/"TSOLVE".
C                                                                                 
      IMPLICIT NONE
C
C     * INTEGER CONSTANTS.
C
      INTEGER ILG,IL1,IL2,IG,I
C
C     * OUTPUT ARRAYS.
C
      !Work arrays used in calculation of GCONST and GCOEFF
      REAL A1    (ILG),    A2    (ILG),    B1    (ILG),
     1     B2    (ILG),    C2    (ILG)  
C  
      REAL GDENOM(ILG)  !Work array used in calculation of GCONST and 
                        !GCOEFF    
      REAL GCOEFF(ILG)  !Multiplier used in equation relating ground 
                        !surface heat flux to surface temperature 
                        ![W m-2 K-1]
      REAL GCONST(ILG)  !Intercept used in equation relating ground 
                        !surface heat flux to surface temperature 
                        ![W m-2 ]
      REAL CPHCHG(ILG)  !Latent heat of sublimation [J kg-1]
C
      INTEGER              IWATER(ILG)  !Flag indicating condition of 
                                        !surface (dry, water-covered or 
                                        !snow-covered)
C
C     * INPUT ARRAYS.
C
      REAL TBAR  (ILG,IG)   !Temperatures of soil layers, averaged over 
                            !modelled area [K] 
      REAL TCTOP (ILG,IG)   !Thermal conductivity of soil at top of 
                            !layer [W m-1 K-1] (lambda_t)
      REAL TCBOT (ILG,IG)   !Thermal conductivity of soil at bottom of 
                            !layer [W m-1 K-1] (lambda_b)

C
      REAL FI    (ILG)  !Fractional coverage of subarea in question on 
                        !modelled area [ ]
      REAL ZPOND (ILG)  !Depth of ponded water on surface [m]  
      REAL TBAR1P(ILG)  !Lumped temperature of ponded water and first 
                        !soil layer [K]
      REAL TCSNOW(ILG)  !Thermal conductivity of snow [W m-1 K-1]  
      REAL ZSNOW (ILG)  !Depth of snow pack [m]
C
      INTEGER              ISAND (ILG,IG)   !Sand content flag
C
      REAL DELZ  (IG)   !Overall thickness of soil layer [m] (delta_z)
C
C     * TEMPORARY VARIABLES.
C
      REAL DELZ1,A3,B3,C3,TCZERO
C
C     * COMMON BLOCK PARAMETERS.
C
      REAL TCW      !Thermal conductivity of water (0.57) [W m-1 K-1]
      REAL TCICE    !Thermal conductivity of ice (2.24) [W m-1 K-1]
      REAL TCSAND   !Thermal conductivity of sand particles (2.5) 
                    ![W m-1 K-1]
      REAL TCCLAY   !Thermal conductivity of fine mineral particles 
                    !(2.5) [W m-1 K-1]
      REAL TCOM     !Thermal conductivity of organic matter (0.25) 
                    ![W m-1 K-1]
      REAL TCDRYS   !Thermal conductivity of dry mineral soil (0.275) 
                    ![W m-1 K-1]
      REAL RHOSOL   !Density of soil mineral matter (2.65*10^3) [kg m-3]
      REAL RHOOM    !Density of soil organic matter (1.30*10^3) [kg m-3]
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
      COMMON /CLASS3/ TCW,TCICE,TCSAND,TCCLAY,TCOM,TCDRYS,
     1                RHOSOL,RHOOM
      COMMON /CLASS4/ HCPW,HCPICE,HCPSOL,HCPOM,HCPSND,HCPCLY,
     1                SPHW,SPHICE,SPHVEG,SPHAIR,RHOW,RHOICE,
     2                TCGLAC,CLHMLT,CLHVAP
C-----------------------------------------------------------------------
C     * INITIALIZATION OF ARRAYS.
C
      !In this subroutine, coefficients are derived for an equation 
      !relating the heat flux at the ground surface to the ground 
      !surface temperature, using the average temperatures and the 
      !thermal conductivities of the underlying first three soil layers. 
      !It is assumed that the variation of temperature T with depth z 
      !within each soil layer can be modelled by using a quadratic 
      !equation:
      !
      !T(z) = (1⁄2)*a*z^2 + b*z + c
      !
      !By substituting 0 for z in the above equation and in the 
      !expressions for its first and second derivatives, it can be shown 
      !that a = T′′(0), b = T′(0), and c = T(0). The term T′′(0) can be 
      !evaluated from the expression for the first derivative evaluated 
      !at the bottom of the soil layer, T(DELZ):
      !
      !T′′(0) = [T′(DELZ) - T′(0)]/DELZ
      !
      !The temperature gradient T′(0) at the top of each layer is 
      !related to the heat flux G(0) through the thermal conductivity 
      !TCTOP; and the temperature gradient and heat flux at the bottom 
      !of the layer, G(DELZ) and T(DELZ), are similarly related through 
      !the bottom thermal conductivity TCBOT: 
      !
      !G(0) = - TCTOP*T′(0)
      !G(DELZ) = - TCBOT*T′(DELZ)
      !
      !The average soil layer temperature, Tav(DELZ), can be obtained by 
      !integrating the resulting equation for T(z) between 0 and DELZ. 
      !Making use of all of the above expressions, recognizing that the 
      !heat fluxes and temperatures at the bottoms of layers 1 and 2 
      !must equal the heat fluxes and temperatures at the tops of layers 
      !2 and 3 respectively, and neglecting as a first approximation the 
      !heat flux at the bottom of the third layer, a linear equation can 
      !be derived relating G(0) to T(0) at the soil surface, where the 
      !slope and intercept of the equation are functions only of the 
      !average temperatures, thicknesses, and top and bottom thermal 
      !conductivities of the three soil layers.
      !
      !In the subroutine loop, first the depth corresponding to TBAR1P 
      !(the lumped temperature of the first soil layer and the ponded 
      !water) is calculated, as the sum of the first soil layer 
      !thickness and the ponded water depth. If the ponded water depth 
      !is not vanishingly small, the surface water flag IWATER is set to 
      !1; otherwise it is set to 0 for soils and 2 for ice sheets 
      !(indicated by ISAND = -4). If IWATER = 2, indicating a frozen 
      !water surface, the latent heat of vaporization, CPHCHG, is set to 
      !the value for sublimation (by adding the latent heat of melting 
      !to the latent heat of vaporization). If there is a snow pack 
      !present, the thermal conductivity at the top of the ground 
      !surface is calculated as the harmonic mean of the thermal 
      !conductivity at the top of the first soil layer and that of the 
      !snow pack. Finally, a series of work arrays is evaluated and is 
      !used to calculate the slope and intercept, GCOEFF and GCONST, of 
      !the equation relating G(0) to T(0) at the ground surface.
      !
      DO 100 I=IL1,IL2
          IF(FI(I).GT.0.)                                          THEN
              DELZ1=DELZ(1)+ZPOND(I)                                                         
              IF(ZPOND(I).GT.0.5E-3)                          THEN
                  IWATER(I)=1
              ELSE                                                                        
                  IF(ISAND(I,1).GT.-4)                THEN
                      IWATER(I)=0
                  ELSE
                      IWATER(I)=2
                  ENDIF
              ENDIF    
C
              IF(IWATER(I).EQ.2)                                    THEN
                  CPHCHG(I)=CLHVAP+CLHMLT
              ELSE                                                                        
                  CPHCHG(I)=CLHVAP
              ENDIF                                                                   
C
              IF(ZSNOW(I).GT.0.0) THEN
                  TCZERO=1.0/(0.5/TCSNOW(I)+0.5/TCTOP(I,1))
              ELSE
                  TCZERO=TCTOP(I,1)
              ENDIF
              A1(I)=DELZ1/(3.0*TCZERO)
              A2(I)=DELZ1/(2.0*TCZERO)
              A3=A2(I)           
              B1(I)=DELZ1/(6.0*TCBOT(I,1)) 
              B2(I)=DELZ1/(2.0*TCBOT(I,1))+DELZ(2)/(3.0*TCTOP(I,2))
              B3=DELZ1/(2.0*TCBOT(I,1))+DELZ(2)/(2.0*TCTOP(I,2))
              C2(I)=DELZ(2)/(6.0*TCBOT(I,2))
              C3=DELZ(2)/(2.0*TCBOT(I,2))+DELZ(3)/(3.0*TCTOP(I,3))
              GDENOM(I)=A1(I)*(B2(I)*C3-B3*C2(I))-B1(I)*(A2(I)*C3-
     1                  A3*C2(I))                                    
              GCOEFF(I)=(B2(I)*C3-B3*C2(I)-B1(I)*(C3-C2(I)))/GDENOM(I) 
              GCONST(I)=(-TBAR1P(I)*(B2(I)*C3-B3*C2(I))+
     1                    TBAR(I,2)*B1(I)*C3-
     2                    TBAR(I,3)*B1(I)*C2(I))/GDENOM(I)           
          ENDIF                                                            
  100 CONTINUE
C                                                                                  
      RETURN                                                                      
      END       
