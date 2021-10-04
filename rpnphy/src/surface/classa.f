      SUBROUTINE CLASSA(FC,     FG,     FCS,    FGS,    ALVSCN, ALIRCN,
     1                  ALVSG,  ALIRG,  ALVSCS, ALIRCS, ALVSSN, ALIRSN,           
     2                  ALVSGC, ALIRGC, ALVSSC, ALIRSC, TRVSCN, TRIRCN, 
     3                  TRVSCS, TRIRCS, FSVF,   FSVFS,  
     4                  RAICAN, RAICNS, SNOCAN, SNOCNS, FRAINC, FSNOWC, 
     5                  FRAICS, FSNOCS, DISP,   DISPS,  ZOMLNC, ZOMLCS, 
     6                  ZOELNC, ZOELCS, ZOMLNG, ZOMLNS, ZOELNG, ZOELNS, 
     7                  CHCAP,  CHCAPS, CMASSC, CMASCS, CWLCAP, CWFCAP,
     8                  CWLCPS, CWFCPS, RC,     RCS,    RBCOEF, FROOT,  
     9                  ZPLIMC, ZPLIMG, ZPLMCS, ZPLMGS, TRSNOW, ZSNOW,  
     A                  WSNOW,  ALVS,   ALIR,   HTCC,   HTCS,   HTC,    
     B                  WTRC,   WTRS,   WTRG,   CMAI,   FSNOW,
     C                  FCANMX, ZOLN,   ALVSC,  ALIRC,  PAIMAX, PAIMIN, 
     D                  CWGTMX, ZRTMAX, RSMIN,  QA50,   VPDA,   VPDB,
     E                  PSIGA,  PSIGB,  PAIDAT, HGTDAT, ACVDAT, ACIDAT, 
     F                  ASVDAT, ASIDAT, AGVDAT, AGIDAT, ALGWET, ALGDRY, 
     G                  THLIQ,  THICE,  TBAR,   RCAN,   SNCAN,  TCAN,   
     H                  GROWTH, SNO,    TSNOW,  RHOSNO, ALBSNO, ZBLEND,
     I                  Z0ORO,  Z0VEG,  SNOLIM, ZPLMG0, ZPLMS0, 
     J                  FCLOUD, TA,     VPD,    RHOAIR, COSZS,  
     K                  QSWINV, RADJ,   DLON,   RHOSNI, DELZ,   DELZW,  
     L                  ZBOTW,  THPOR,  THLMIN, PSISAT, BI,     PSIWLT, 
     M                  HCPS,   ISAND,  
     N                  FCANCMX,ICTEM,  ICTEMMOD, RMATC, ZOLNC,CMASVEGC,
     O                  AILC,   PAIC,   L2MAX,  NOL2PFTS, SLAIC,
     P                  AILCG,  AILCGS, FCANC,  FCANCS,
     Q                  IDAY,   ILG,    IL1,    IL2,    
     R                  JL,N,   IC,     ICP1,   IG,     IDISP,  IZREF,
     S                  IWF,    IPAI,   IHGT,   IALC,   IALS,   IALG,
     T                  ALVSCTM, ALIRCTM )

C
C     Purpose: Organize calculation of radiation-related and other 
C     surface parameters.
C
C===================== CTEM =====================================/
C
C     * OCT 04/21 - K.WINGER.   ADD Z0VEG
C     * SEP 05/12 - J.MELTON.   REMOVED UNUSED VAR, CWCPAV
C     * NOV 14/11 - M.LAZARE.   IMPLEMENT CTEM SUPPORT, PRIMARILY
C     *                         INVOLVING ADDITIONAL FIELDS TO PASS
C     *                         IN/OUT OF NEW APREP ROUTINE. THIS 
C     *                         INCLUDES NEW INPUT ARRAY "PAIC".
C     * NOV 30/06 - D.VERSEGHY. CONVERT RADJ TO REGULAR PRECISION.
C     * APR 13/06 - D.VERSEGHY. SEPARATE GROUND AND SNOW ALBEDOS FOR 
C     *                         OPEN AND CANOPY-COVERED AREAS; KEEP
C     *                         FSNOW AS OUTPUT ARRAY.
C     * APR 06/06 - D.VERSEGHY. INTRODUCE MODELLING OF WSNOW.
C     * MAR 14/05 - D.VERSEGHY. RENAME SCAN TO SNCAN (RESERVED NAME
C     *                         IN F90); CHANGE SNOLIM FROM CONSTANT
C     *                         TO VARIABLE.
C     * NOV 03/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.
C     * DEC 05/02 - D.VERSEGHY. NEW PARAMETERS FOR APREP.
C     * JUL 31/02 - D.VERSEGHY. MODIFICATIONS ASSOCIATED WITH NEW
C     *                         CALCULATION OF STOMATAL RESISTANCE.
C     *                         SHORTENED CLASS3 COMMON BLOCK.
C     * JUL 23/02 - D.VERSEGHY. MODIFICATIONS TO MOVE ADDITION OF AIR
C     *                         TO CANOPY MASS INTO APREP; SHORTENED
C     *                         CLASS4 COMMON BLOCK.
C     * MAR 18/02 - D.VERSEGHY. NEW CALLS TO ALL SUBROUTINES TO ENABLE
C     *                         ASSIGNMENT OF USER-SPECIFIED VALUES TO
C     *                         ALBEDOS AND VEGETATION PROPERTIES; NEW
C     *                         "CLASS8" COMMON BLOCK; MOVE CALCULATION 
C     *                         OF "FCLOUD" INTO CLASS DRIVER.
C     * SEP 19/00 - D.VERSEGHY. PASS ADDITIONAL ARRAYS TO APREP IN COMMON 
C     *                         BLOCK CLASS7, FOR CALCULATION OF NEW 
C     *                         STOMATAL RESISTANCE COEFFICIENTS USED
C     *                         IN TPREP.
C     * APR 12/00 - D.VERSEGHY. RCMIN NOW VARIES WITH VEGETATION TYPE:
C     *                         PASS IN BACKGROUND ARRAY "RCMINX".
C     * DEC 16/99 - D.VERSEGHY. ADD "XLEAF" ARRAY TO CLASS7 COMMON BLOCK
C     *                         AND CALCULATION OF LEAF DIMENSION PARAMETER
C     *                         "DLEAF" IN APREP.
C     * NOV 16/98 - M.LAZARE.   "DLON" NOW PASSED IN AND USED DIRECTLY
C     *                         (INSTEAD OF INFERRING FROM "LONSL" AND 
C     *                         "ILSL" WHICH USED TO BE PASSED) TO PASS
C     *                         TO APREP TO CALCULATE GROWTH INDEX. THIS
C     *                         IS DONE TO MAKE THE PHYSICS PLUG COMPATIBLE
C     *                         FOR USE WITH THE RCM WHICH DOES NOT HAVE
C     *                         EQUALLY-SPACED LONGITUDES.
C     * JUN 20/97 - D.VERSEGHY. CLASS - VERSION 2.7.
C     *                         MODIFICATIONS TO ALLOW FOR VARIABLE
C     *                         SOIL PERMEABLE DEPTH.
C     * SEP 27/96 - D.VERSEGHY. CLASS - VERSION 2.6.
C     *                         FIX BUG TO CALCULATE GROUND ALBEDO
C     *                         UNDER CANOPIES AS WELL AS OVER BARE
C     *                         SOIL.
C     * JAN 02/96 - D.VERSEGHY. CLASS - VERSION 2.5.
C     *                         COMPLETION OF ENERGY BALANCE
C     *                         DIAGNOSTICS.
C     *                         ALSO, PASS IDISP TO SUBROUTINE APREP.
C     * AUG 30/95 - D.VERSEGHY. CLASS - VERSION 2.4.
C     *                         VARIABLE SURFACE DETENTION CAPACITY
C     *                         IMPLEMENTED.
C     * AUG 16/95 - D.VERSEGHY. THREE NEW ARRAYS TO COMPLETE WATER
C     *                         BALANCE DIAGNOSTICS.
C     * OCT 14/94 - D.VERSEGHY. CLASS - VERSION 2.3.
C     *                         REVISE CALCULATION OF FCLOUD TO
C     *                         HANDLE CASES WHERE INCOMING SOLAR
C     *                         RADIATION IS ZERO AT LOW SUN ANGLES.
C     * NOV 24/92 - M.LAZARE.   CLASS - VERSION 2.1.
C     *                         MODIFIED FOR MULTIPLE LATITUDES.
C     * OCT 13/92 - D.VERSEGHY/M.LAZARE. REVISED AND VECTORIZED CODE
C     *                                  FOR MODEL VERSION GCM7.
C     * AUG 12/91 - D.VERSEGHY. CODE FOR MODEL VERSION GCM7U -
C     *                         CLASS VERSION 2.0 (WITH CANOPY).
C     * APR 11/89 - D.VERSEGHY. VISIBLE AND NEAR-IR ALBEDOS AND 
C     *                         TRANSMISSIVITIES FOR COMPONENTS OF
C     *                         LAND SURFACE.
C
      IMPLICIT NONE
C
C     * INTEGER CONSTANTS.
C
      INTEGER IDAY,ILG,IL1,IL2,JL,IC,ICP1,IG,IDISP,IZREF,IWF,
     1        IPAI,IHGT,IALC,IALS,IALG,I,J,N
C
C     * OUTPUT ARRAYS.
C
      REAL FC    (ILG)  !Subarea fractional coverage of modelled area 
                        ![ ]
      REAL FG    (ILG)  !Subarea fractional coverage of modelled area 
                        ![ ]
      REAL FCS   (ILG)  !Subarea fractional coverage of modelled area 
                        ![ ]
      REAL FGS   (ILG)  !Subarea fractional coverage of modelled area 
                        ![ ]
      REAL ALVSCN(ILG)  !Visible albedo of vegetation over bare ground 
                        ![ ]
      REAL ALIRCN(ILG)  !Near-IR albedo of vegetation over bare ground 
                        ![ ]
      REAL ALVSG (ILG)  !Visible albedo of open bare ground [ ]
      REAL ALIRG (ILG)  !Near-IR albedo of open bare ground [ ]   
      REAL ALVSCS(ILG)  !Visible albedo of vegetation over snow [ ]
      REAL ALIRCS(ILG)  !Near-IR albedo of vegetation over snow [ ]
      REAL ALVSSN(ILG)  !Visible albedo of open snow cover [ ]
      REAL ALIRSN(ILG)  !Near-IR albedo of open snow cover [ ]
      REAL ALVSGC(ILG)  !Visible albedo of bare ground under vegetation 
                        ![ ]
      REAL ALIRGC(ILG)  !Near-IR albedo of bare ground under vegetation 
                        ![ ]
      REAL ALVSSC(ILG)  !Visible albedo of snow under vegetation [ ] 
      REAL ALIRSC(ILG)  !Near-IR albedo of snow under vegetation [ ]
      REAL TRVSCN(ILG)  !Visible transmissivity of vegetation over bare 
                        !ground [ ]
      REAL TRIRCN(ILG)  !Near-IR transmissivity of vegetation over bare 
                        !ground [ ] 
      REAL TRVSCS(ILG)  !Visible transmissivity of vegetation over snow 
                        ![ ]
      REAL TRIRCS(ILG)  !Near-IR transmissivity of vegetation over snow 
                        ![ ]
      REAL FSVF  (ILG)  !Sky view factor for bare ground under canopy 
                        ![ ]
      REAL FSVFS (ILG)  !Sky view factor for snow under canopy [ ] 
      REAL RAICAN(ILG)  !Intercepted liquid water stored on canopy over 
                        !bare ground [kg m-2]
      REAL RAICNS(ILG)  !Intercepted liquid water stored on canopy over 
                        !snow [kg m-2]
      REAL SNOCAN(ILG)  !Intercepted frozen water stored on canopy over 
                        !bare soil [kg m-2]
      REAL SNOCNS(ILG)  !Intercepted frozen water stored on canopy over 
                        !snow [kg m-2]
      REAL FRAINC(ILG)  !Fractional coverage of canopy by liquid water 
                        !over snow-free subarea [ ]
      REAL FSNOWC(ILG)  !Fractional coverage of canopy by frozen water 
                        !over snow-free subarea [ ]
      REAL FRAICS(ILG)  !Fractional coverage of canopy by liquid water 
                        !over snow-covered subarea [ ]
      REAL FSNOCS(ILG)  !Fractional coverage of canopy by frozen water 
                        !over snow-covered subarea [ ]
      REAL DISP  (ILG)  !Displacement height of vegetation over bare 
                        !ground [m]
      REAL DISPS (ILG)  !Displacement height of vegetation over snow [m] 
      REAL ZOMLNC(ILG)  !Logarithm of roughness length for momentum of 
                        !vegetation over bare ground [ ]
      REAL ZOMLCS(ILG)  !Logarithm of roughness length for momentum of 
                        !vegetation over snow [ ]
      REAL ZOELNC(ILG)  !Logarithm of roughness length for heat of 
                        !vegetation over bare ground [ ]
      REAL ZOELCS(ILG)  !Logarithm of roughness length for heat of 
                        !vegetation over snow [ ]
      REAL ZOMLNG(ILG)  !Logarithm of roughness length for momentum of 
                        !bare ground [ ]
      REAL ZOMLNS(ILG)  !Logarithm of roughness length for momentum of 
                        !snow [ ]
      REAL ZOELNG(ILG)  !Logarithm of roughness length for heat of bare 
                        !ground [ ]
      REAL ZOELNS(ILG)  !Logarithm of roughness length for heat of snow 
                        ![ ]
      REAL CHCAP (ILG)  !Heat capacity of canopy over bare ground 
                        ![J m-2 K-1]
      REAL CHCAPS(ILG)  !Heat capacity of canopy over snow [J m-2 K-1] 
      REAL CMASSC(ILG)  !Mass of canopy over bare ground [kg m-2] 
      REAL CMASCS(ILG)  !Mass of canopy over snow [kg m-2] 
      REAL CWLCAP(ILG)  !Storage capacity of canopy over bare ground for 
                        !liquid water [kg m-2]
      REAL CWFCAP(ILG)  !Storage capacity of canopy over bare ground for 
                        !frozen water [kg m-2]
      REAL CWLCPS(ILG)  !Storage capacity of canopy over snow for liquid 
                        !water [kg m-2]
      REAL CWFCPS(ILG)  !Storage capacity of canopy over snow for frozen 
                        !water [kg m-2]
      REAL RC    (ILG)  !Stomatal resistance of vegetation over bare 
                        !ground [s m-1]
      REAL RCS   (ILG)  !Stomatal resistance of vegetation over snow 
                        ![s m-1]
      REAL ZPLIMC(ILG)  !Maximum water ponding depth for ground under 
                        !canopy [m]
      REAL ZPLIMG(ILG)  !Maximum water ponding depth for bare ground [m] 
      REAL ZPLMCS(ILG)  !Maximum water ponding depth for ground under 
                        !snow under canopy [m]
      REAL ZPLMGS(ILG)  !Maximum water ponding depth for ground under 
                        !snow [m]
      REAL RBCOEF(ILG)  !Parameter for calculation of leaf boundary 
                        !resistance
      REAL TRSNOW(ILG)  !Short-wave transmissivity of snow pack [ ] 
      REAL ZSNOW (ILG)  !Depth of snow pack [m] (zs) 
      REAL WSNOW (ILG)  !Liquid water content of snow pack [kg m-2]
      REAL ALVS  (ILG)  !Diagnosed total visible albedo of land surface 
                        ![ ]
      REAL ALIR  (ILG)  !Diagnosed total near-infrared albedo of land 
                        !surface [ ]
      REAL HTCC  (ILG)  !Diagnosed internal energy change of vegetation 
                        !canopy due to conduction and/or change in mass 
                        ![W m -2]
      REAL HTCS  (ILG)  !Diagnosed internal energy change of snow pack 
                        !due to conduction and/or change in mass [W m-2]
      REAL WTRC  (ILG)  !Diagnosed residual water transferred off the 
                        !vegetation canopy [kg m-2 s-1]
      REAL WTRS  (ILG)  !Diagnosed residual water transferred into or 
                        !out of the snow pack [kg m-2 s-1]
      REAL WTRG  (ILG)  !Diagnosed residual water transferred into or 
                        !out of the soil [kg m-2 s-1]
      REAL CMAI  (ILG)  !Aggregated mass of vegetation canopy [kg m-2]
      REAL FSNOW (ILG)  !Diagnosed fractional snow coverage [ ]
C
      REAL FROOT (ILG,IG)   !Fraction of total transpiration contributed 
                            !by soil layer [ ]  
      REAL HTC   (ILG,IG)   !Diagnosed internal energy change of soil 
                            !layer due to conduction and/or change in 
                            !mass [W m-2]
C
C     * INPUT ARRAYS DEPENDENT ON LONGITUDE.
C  
      REAL FCANMX(ILG,ICP1) !Maximum fractional coverage of modelled 
                            !area by vegetation category [ ] 
      REAL ZOLN  (ILG,ICP1) !Natural logarithm of maximum roughness 
                            !length of vegetation category [ ]
      REAL ALVSC (ILG,ICP1) !Background average visible albedo of 
                            !vegetation category [ ]
      REAL ALIRC (ILG,ICP1) !Background average near-infrared albedo of 
                            !vegetation category [ ]
      REAL PAIMAX(ILG,IC)   !Maximum plant area index of vegetation 
                            !category [ ]
      REAL PAIMIN(ILG,IC)   !Minimum plant area index of vegetation 
                            !category [ ]
      REAL CWGTMX(ILG,IC)   !Maximum canopy mass for vegetation category 
                            ![kg m-2]
      REAL ZRTMAX(ILG,IC)   !Maximum rooting depth of vegetation 
                            !category [m]
      REAL RSMIN (ILG,IC)   !Minimum stomatal resistance of vegetation 
                            !category [s m-1]
      REAL QA50  (ILG,IC)   !Reference value of incoming shortwave 
                            !radiation for vegetation category (used in 
                            !stomatal resistance calculation) [W m-2]
      REAL VPDA  (ILG,IC)   !Vapour pressure deficit coefficient for 
                            !vegetation category (used in stomatal 
                            !resistance calculation) [ ]
      REAL VPDB  (ILG,IC)   !Vapour pressure deficit coefficient for 
                            !vegetation category (used in stomatal 
                            !resistance calculation) [ ]
      REAL PSIGA (ILG,IC)   !Soil moisture suction coefficient for 
                            !vegetation category (used in stomatal
                            !resistance calculation) [ ]
      REAL PSIGB (ILG,IC)   !Soil moisture suction coefficient for 
                            !vegetation category (used in stomatal
                            !resistance calculation) [ ]
      REAL PAIDAT(ILG,IC)   !Optional user-specified value of plant area 
                            !indices of vegetation categories to 
                            !override CLASS-calculated values [ ]
      REAL HGTDAT(ILG,IC)   !Optional user-specified values of height of 
                            !vegetation categories to override CLASS-
                            !calculated values [m]
      REAL ACVDAT(ILG,IC)   !Optional user-specified value of canopy 
                            !visible albedo to override CLASS-calculated 
                            !value [ ]
      REAL ACIDAT(ILG,IC)   !Optional user-specified value of canopy 
                            !near-infrared albedo to override CLASS-
                            !calculated value [ ]
      REAL THLIQ (ILG,IG)   !Volumetric liquid water content of soil 
                            !layers [m3 m-3]
      REAL THICE (ILG,IG)   !Volumetric frozen water content of soil 
                            !layers [m3 m-3]
      REAL TBAR  (ILG,IG)   !Temperature of soil layers [K]
C
      REAL ASVDAT(ILG)  !Optional user-specified value of snow visible 
                        !albedo to override CLASS-calculated value [ ]
      REAL ASIDAT(ILG)  !Optional user-specified value of snow near-
                        !infrared albedo to override CLASS-calculated 
                        !value [ ]
      REAL AGVDAT(ILG)  !Optional user-specified value of ground visible 
                        !albedo to override CLASS-calculated value [ ]
      REAL AGIDAT(ILG)  !Optional user-specified value of ground near-
                        !infrared albedo to override CLASS-calculated 
                        !value [ ]
      REAL ALGWET(ILG)  !Reference albedo for saturated soil [ ] 
      REAL ALGDRY(ILG)  !Reference albedo for dry soil [ ] 
      REAL RHOSNI(ILG)  !Density of fresh snow [kg m-3] 
      REAL Z0ORO (ILG)  !Orographic roughness length [m]
      REAL Z0VEG (ILG)  !Vegetation roughness length [m]
      REAL RCAN  (ILG)  !Intercepted liquid water stored on canopy 
                        ![kg m-2]
      REAL SNCAN (ILG)  !Intercepted frozen water stored on canopy 
                        ![kg m-2]
      REAL TCAN  (ILG)  !Vegetation canopy temperature [K] 
      REAL GROWTH(ILG)  !Vegetation growth index [ ]    
      REAL SNO   (ILG)  !Mass of snow pack [kg m-2] (Ws) 
      REAL TSNOW (ILG)  !Snowpack temperature [K] 
      REAL RHOSNO(ILG)  !Density of snow [kg m-3] (rho_s)
      REAL ALBSNO(ILG)  !Snow albedo [ ]
      REAL FCLOUD(ILG)  !Fractional cloud cover [ ] 
      REAL TA    (ILG)  !Air temperature at reference height [K] 
      REAL VPD   (ILG)  !Vapour pressure deficit of air [mb] 
      REAL RHOAIR(ILG)  !Density of air [kg m-3]
      REAL COSZS (ILG)  !Cosine of solar zenith angle [ ] 
      REAL QSWINV(ILG)  !Visible radiation incident on horizontal 
                        !surface [W m-2]
      REAL DLON  (ILG)  !Longitude of grid cell (east of Greenwich) 
                        ![degrees]
      REAL ZBLEND(ILG)  !Atmospheric blending height for surface 
                        !roughness length averaging [m]
      REAL SNOLIM(ILG)  !Limiting snow depth below which coverage is 
                        !< 100% [m] (zs,lim)
      REAL ZPLMG0(ILG)  !Maximum water ponding depth for snow-free 
                        !subareas (user-specified when running MESH 
                        !code) [m]
      REAL ZPLMS0(ILG)  !Maximum water ponding depth for snow-covered 
                        !subareas (user-specified when running MESH 
                        !code) [m]
      REAL RADJ  (ILG)  !Latitude of grid cell (positive north of 
                        !equator) [rad]
C
C    * SOIL PROPERTY ARRAYS.
C
      REAL DELZW (ILG,IG)   !Permeable thickness of soil layer [m]  
      REAL ZBOTW (ILG,IG)   !Depth to permeable bottom of soil layer [m]
      REAL THPOR (ILG,IG)   !Pore volume in soil layer [m3 m-3]
      REAL THLMIN(ILG,IG)   !Residual soil liquid water content 
                            !remaining after freezing or evaporation
                            ![m3 m-3]
      REAL PSISAT(ILG,IG)   !Soil moisture suction at saturation [m]
      REAL BI    (ILG,IG)   !Clapp and Hornberger empirical “b” 
                            !parameter [ ]
      REAL PSIWLT(ILG,IG)   !Soil moisture suction at wilting point [m]
      REAL HCPS  (ILG,IG)   !Volumetric heat capacity of soil particles 
                            ![J m-3]
C
      INTEGER   ISAND (ILG,IG)  !Sand content flag
C
C     * OTHER DATA ARRAYS WITH NON-VARYING VALUES.
C
      REAL GROWYR(18,4,2),  DELZ  (IG)  !Soil layer thickness [m]     
      REAL ZORAT (4),       CANEXT(4),       XLEAF (4)
C
C     * CTEM-RELATED FIELDS.

C     * AILCG  - GREEN LAI FOR USE IN PHOTOSYNTHESIS
C     * AILCGS - GREEN LAI FOR CANOPY OVER SNOW SUB-AREA
C     * AILCMIN- MIN. LAI FOR CTEM PFTs
C     * AILCMAX- MAX. LAI FOR CTEM PFTs
C     * L2MAX  - MAXIMUM NUMBER OF LEVEL 2 CTEM PFTs
C     * NOL2PFTS - NUMBER OF LEVEL 2 CTEM PFTs
C     * FCANC  - FRACTION OF CANOPY OVER GROUND FOR CTEM's 9 PFTs
C     * FCANCS - FRACTION OF CANOPY OVER SNOW FOR CTEM's 9 PFTs
C     * SEE BIO2STR SUBROUTINE FOR DEFINITION OF OTHER CTEM VARIABLES
C
      REAL FCANCMX(ILG,ICTEM),   RMATC(ILG,IC,IG),
     1      AILC  (ILG,IC),      PAIC (ILG,IC),
     2      AILCG (ILG,ICTEM),   AILCGS(ILG,ICTEM),
     3      FCANC(ILG,ICTEM),    FCANCS(ILG,ICTEM),
     4      ZOLNC(ILG,IC),       CMASVEGC(ILG,IC),
     5      SLAIC(ILG,IC),       ALVSCTM(ILG,IC),
     6      ALIRCTM(ILG,IC)      

      INTEGER ICTEM, ICTEMMOD, L2MAX, NOL2PFTS(IC)
C                                                                                 
C     * INTERNAL WORK ARRAYS FOR THIS AND ASSOCIATED SUBROUTINES.
C
      REAL RMAT (ILG,IC,IG),H     (ILG,IC),  HS    (ILG,IC),
     1     PAI   (ILG,IC),  PAIS  (ILG,IC),  FCAN  (ILG,IC),  
     2     FCANS (ILG,IC),  CXTEFF(ILG,IC),  AIL   (ILG,IC),
     3     RCACC (ILG,IC),  RCG   (ILG,IC),  RCV   (ILG,IC)
C
      REAL PSIGND(ILG),     
     1     GROWA (ILG),     GROWN (ILG),     GROWB (ILG),     
     2     RRESID(ILG),     SRESID(ILG),     FRTOT (ILG),
     3     TRVS  (ILG),     TRIR  (ILG),     RCT   (ILG),     
     4     GC    (ILG),     PAICAN(ILG),     PAICNS(ILG) 
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
      REAL PI       !Pi
      REAL ZOLNG    !Natural log of roughness length of soil (-4.605)
      REAL ZOLNS    !Natural log of roughness length of snow (-6.908)
      REAL ZOLNI    !Natural log of roughness length of ice (-6.215)
      REAL ZORATG   !Ratio of soil roughness for momentum to roughness
                    !length for heat (3.0)
      REAL ALVSI    !Visible albedo of ice (0.95)
      REAL ALIRI    !Near-infrared albedo of ice (0.73)
      REAL ALVSO    !Visible albedo of organic matter (0.05)
      REAL ALIRO    !Near-infrared albedo of organic matter (0.30)
      REAL ALBRCK   !Albedo of rock (0.27)

      COMMON /CLASS1/ DELT,TFREZ                                               
      COMMON /CLASS2/ RGAS,RGASV,GRAV,SBC,VKC,CT,VMIN
      COMMON /CLASS3/ TCW,TCICE,TCSAND,TCCLAY,TCOM,TCDRYS,
     1                RHOSOL,RHOOM
      COMMON /CLASS4/ HCPW,HCPICE,HCPSOL,HCPOM,HCPSND,HCPCLY,
     1                SPHW,SPHICE,SPHVEG,SPHAIR,RHOW,RHOICE,
     2                TCGLAC,CLHMLT,CLHVAP
      COMMON /CLASS6/ PI,GROWYR,ZOLNG,ZOLNS,ZOLNI,ZORAT,ZORATG                   
      COMMON /CLASS7/ CANEXT,XLEAF
      COMMON /CLASS8/ ALVSI,ALIRI,ALVSO,ALIRO,ALBRCK
                                                                                  
C------------------------------------------------------------------
      !
      !In the first loop, the depth of snow ZSNOW is calculated from the 
      !snow mass SNO and density RHOSNO as:
      !
      !ZSNOW = SNO/RHOSNO.
      !
      !If the calculated value of ZSNOW is less than the limiting snow 
      !depth SNOLIM, the snow cover is deemed to be discontinuous. The 
      !fractional snow coverage FSNOW of the modelled area is evaluated 
      !as
      !
      !FSNOW = ZSNOW/SNOLIM
      !
      !and the snow depth is reset to SNOLIM. The water content of the 
      !snow pack is corrected according to the new snow fractional area.
      !
      !The subarea albedo and transmissivity arrays (for canopy, bare 
      !ground, canopy over snow and snow over bare ground) are next 
      !initialized to zero, and the four CLASSA subsidiary subroutines 
      !are called in turn: APREP to evaluate various model parameters 
      !for the four subareas, SNOALBA to calculate the snow albedo and 
      !transmissivity, GRALB to calculate the ground surface albedo, and 
      !CANALB to calculate the canopy albedo, transmissivity and 
      !stomatal resistance. Finally, the overall visible and near-
      !infrared albedos for the modelled area are determined as weighted 
      !averages over the four subareas.
      !
C     * CALCULATION OF SNOW DEPTH ZSNOW AND FRACTIONAL SNOW COVER
C     * FSNOW; INITIALIZATION OF COMPUTATIONAL ARRAYS. 
C                                                                                  
      DO 100 I=IL1,IL2                                                            
          IF(SNO(I).GT.0.0) THEN
              ZSNOW(I)=SNO(I)/RHOSNO(I)                                       
              IF(ZSNOW(I).GE.(SNOLIM(I)-0.00001)) THEN                                     
                  FSNOW(I)=1.0                                                   
              ELSE                                                            
                  FSNOW(I)=ZSNOW(I)/SNOLIM(I)
                  ZSNOW(I)=SNOLIM(I)
                  WSNOW(I)=WSNOW(I)/FSNOW(I)
              ENDIF                                                           
          ELSE                                                                
              ZSNOW(I)=0.0                                                    
              FSNOW(I)=0.0                                                       
          ENDIF
C
          ALVSCN(I)=0.0                                                   
          ALIRCN(I)=0.0                                                   
          ALVSCS(I)=0.0  
          ALIRCS(I)=0.0    
          TRVSCN(I)=0.0                                                   
          TRIRCN(I)=0.0                                                   
          TRVSCS(I)=0.0                                                   
          TRIRCS(I)=0.0
          ALVSSN(I)=0.0                                                   
          ALIRSN(I)=0.0                                                   
          ALVSG (I)=0.0
          ALIRG (I)=0.0
          ALVSGC(I)=0.0
          ALIRGC(I)=0.0
          ALVSSC(I)=0.0
          ALIRSC(I)=0.0
          TRSNOW(I)=0.0                                                       
  100 CONTINUE
C
C ===================== CTEM =====================================\
C     IF USING DYNAMIC VEGETATION COMPONENT OF CTEM, REPLACE ALBEDOS
C     THAT ARE BASED ON CTEM.

      IF(ICTEMMOD.EQ.1)THEN
        DO J = 1, IC
          DO I = IL1, IL2
            ALVSC(I,J)=ALVSCTM(I,J)
            ALIRC(I,J)=ALIRCTM(I,J)
          ENDDO
        ENDDO
      ENDIF
C===================== CTEM =====================================/
C
C     * PREPARATION.
C
      CALL APREP (FC,FG,FCS,FGS,PAICAN,PAICNS,FSVF,FSVFS, 
     1            FRAINC,FSNOWC,FRAICS,FSNOCS,RAICAN,RAICNS,SNOCAN,
     2            SNOCNS,DISP,DISPS,ZOMLNC,ZOMLCS,ZOELNC,ZOELCS,
     3            ZOMLNG,ZOMLNS, ZOELNG,ZOELNS,CHCAP,CHCAPS,CMASSC,
     4            CMASCS,CWLCAP,CWFCAP,CWLCPS,CWFCPS,RBCOEF,FROOT,
     5            ZPLIMC,ZPLIMG,ZPLMCS,ZPLMGS,HTCC,HTCS,HTC,WTRC,
     6            WTRS,WTRG,CMAI,PAI,PAIS,AIL,FCAN,FCANS,PSIGND,
     7            FCANMX,ZOLN,PAIMAX,PAIMIN,CWGTMX,ZRTMAX,
     8            PAIDAT,HGTDAT,THLIQ,THICE,TBAR,RCAN,SNCAN,
     9            TCAN,GROWTH,ZSNOW,TSNOW,FSNOW,RHOSNO,SNO,Z0ORO,Z0VEG,
     A            ZBLEND,ZPLMG0,ZPLMS0,
     B            TA,RHOAIR,RADJ,DLON,RHOSNI,DELZ,DELZW,ZBOTW,
     C            THPOR,THLMIN,PSISAT,BI,PSIWLT,HCPS,ISAND,
     D            ILG,IL1,IL2,JL,IC,ICP1,IG,IDAY,IDISP,IZREF,IWF,
     E            IPAI,IHGT,RMAT,H,HS,GROWA,GROWN,GROWB,
     F            RRESID,SRESID,FRTOT,
     G            FCANCMX,ICTEM,ICTEMMOD,RMATC,
     H            AILC,PAIC,AILCG,L2MAX,NOL2PFTS,
     I            AILCGS,FCANCS,FCANC,ZOLNC,CMASVEGC,SLAIC)

C     * SNOW ALBEDOS AND TRANSMISSIVITY.
C 
      CALL SNOALBA(ALVSSN,ALIRSN,ALVSSC,ALIRSC,ALBSNO,TRSNOW,
     1             ZSNOW,FSNOW,ASVDAT,ASIDAT,
     2             ILG,IG,IL1,IL2,JL,IALS)
C
C     * BARE SOIL ALBEDOS.
C
      CALL GRALB(ALVSG,ALIRG,ALVSGC,ALIRGC,
     1            ALGWET,ALGDRY,THLIQ,ALVSC(1,5),ALIRC(1,5),
     2            FCANMX(1,5),AGVDAT,AGIDAT,ISAND,
     3            ILG,IG,IL1,IL2,JL,IALG) 
C
C     * CANOPY ALBEDOS AND TRANSMISSIVITIES, AND VEGETATION
C     * STOMATAL RESISTANCE.
C
      CALL CANALB(ALVSCN,ALIRCN,ALVSCS,ALIRCS,TRVSCN,TRIRCN,
     1            TRVSCS,TRIRCS,RC,RCS,
     2            ALVSC,ALIRC,RSMIN,QA50,VPDA,VPDB,PSIGA,PSIGB,
     3            FC,FCS,FSNOW,FSNOWC,FSNOCS,FCAN,FCANS,PAI,PAIS,
     4            AIL,PSIGND,FROOT,FCLOUD,COSZS,QSWINV,VPD,TA,
     5            ACVDAT,ACIDAT,ALVSGC,ALIRGC,ALVSSC,ALIRSC,
     6            ILG,IL1,IL2,JL,IC,ICP1,IG,IALC,
     7            CXTEFF,TRVS,TRIR,RCACC,RCG,RCV,RCT,GC) 
C
C     * EFFECTIVE WHOLE-SURFACE VISIBLE AND NEAR-IR ALBEDOS.
C
      DO 500 I=IL1,IL2
          ALVS(I)=FC(I)*ALVSCN(I)+FG(I)*ALVSG(I)+FCS(I)*ALVSCS(I)+
     1            FGS(I)*ALVSSN(I)                                                
          ALIR(I)=FC(I)*ALIRCN(I)+FG(I)*ALIRG(I)+FCS(I)*ALIRCS(I)+            
     1            FGS(I)*ALIRSN(I)                                                
  500 CONTINUE
C
      RETURN                                                                      
      END
