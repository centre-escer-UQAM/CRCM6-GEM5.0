      SUBROUTINE CLASSW(THLIQ,  THICE,  TBAR,   TCAN,   RCAN,   SNCAN,
     1                  RUNOFF, TRUNOF, SNO,    TSNOW,  RHOSNO, ALBSNO,
     2                  WSNOW,  ZPOND,  TPOND,  GROWTH, TBASE,  GFLUX,
     3                  PCFC,   PCLC,   PCPN,   PCPG,   QFCF,   QFCL,
     4                  QFN,    QFG,    QFC,    HMFC,   HMFG,   HMFN,
     5                  HTCC,   HTCS,   HTC,    ROFC,   ROFN,   ROVG,
     6                  WTRS,   WTRG,   OVRFLW, SUBFLW, BASFLW,
     7                  TOVRFL, TSUBFL, TBASFL, EVAP,
     8                  TBARC,  TBARG,  TBARCS, TBARGS, THLIQC, THLIQG,
     9                  THICEC, THICEG, HCPC,   HCPG,   RPCP,   TRPCP,
     A                  SPCP,   TSPCP,  PCPR,   TA,     RHOSNI, GGEO,
     B                  FC,     FG,     FCS,    FGS,    TPONDC, TPONDG,
     C                  TPNDCS, TPNDGS, EVAPC,  EVAPCG, EVAPG,  EVAPCS,
     D                  EVPCSG, EVAPGS, QFREZC, QFREZG, QMELTC, QMELTG,
     E                  RAICAN, SNOCAN, RAICNS, SNOCNS, FROOT,  FSVF,
     F                  FSVFS,  CWLCAP, CWFCAP, CWLCPS, CWFCPS, TCANO,
     G                  TCANS,  CHCAP,  CHCAPS, CMASSC, CMASCS, ZSNOW,
     H                  GZEROC, GZEROG, GZROCS, GZROGS, G12C,   G12G,
     I                  G12CS,  G12GS,  G23C,   G23G,   G23CS,  G23GS,
     J                  TSNOCS, TSNOGS, WSNOCS, WSNOGS, RHOSCS, RHOSGS,
     K                  ZPLIMC, ZPLIMG, ZPLMCS, ZPLMGS, TSFSAV,
     L                  TCTOPC, TCBOTC, TCTOPG, TCBOTG,
     M                  THPOR,  THLRET, THLMIN, BI,     PSISAT, GRKSAT,
     N                  THLRAT, THFC,   XDRAIN, HCPS,   DELZ,
     O                  DELZW,  ZBOTW,  XSLOPE, GRKFAC, WFSURF, WFCINT,
     P                  ISAND,  IGDR,
     Q                  IWF,    ILG,    IL1,    IL2,    N,
     R                  JL,     IC,     IG,     IGP1,   IGP2,
c     S                  NLANDCS,NLANDGS,NLANDC, NLANDG, NLANDI,
     S                  NLANDCS,NLANDGS,NLANDC, NLANDG, NLANDI, WT,
     S                  QDIST, QINT, QINC,QING,QINGS,QINCS,QDISC, Isat,
     T                  WTnew,DMLIQT,Ibed,Isats,ZWT,EXCWT,SLP,ARE,LEG,
     +                  XSANI,igwscheme)
!!!     T                  DD, INTFL,SANI, INTF_FLAG) !Huziy added DD, INTFL - drainage density for interflow model

C                                                                        
C     Purpose: Call subroutines to perform surface water budget 
C     calculations,
C
C     * OCT 18/11 - M.LAZARE.   PASS IN IGDR THROUGH CALLS TO
C     *                         GRDRAN/GRINFL (ORIGINATES NOW
C     *                         IN CLASSB - ONE CONSISTENT
C     *                         CALCULATION).

C     * APR 04/11 - D.VERSEGHY. ADD DELZ TO GRINFL CALL.
C     * DEC 07/09 - D.VERSEGHY. ADD RADD AND SADD TO WPREP CALL.
C     * JAN 06/09 - D.VERSEGHY. INCREASE LIMITING SNOW AMOUNT.
C     * FEB 25/08 - D.VERSEGHY. MODIFICATIONS REFLECTING CHANGES
C     *                         ELSEWHERE IN CODE.
C     * MAR 23/06 - D.VERSEGHY. CHANGES TO ADD MODELLING OF WSNOW;
C     *                         PASS IN GEOTHERMAL HEAT FLUX.
C     * MAR 21/06 - P.BARTLETT. PASS ADDITIONAL VARIABLES TO WPREP.
C     * DEC 07/05 - D.VERSEGHY. REVISIONS TO CALCULATION OF TBASE.
C     * OCT 05/05 - D.VERSEGHY. MODIFICATIONS TO ALLOW OPTION OF SUB-
C     *                         DIVIDING THIRD SOIL LAYER.
C     * MAR 23/05 - D.VERSEGHY. ADD VARIABLES TO SUBROUTINE CALLS.
C     * MAR 14/05 - D.VERSEGHY. RENAME SCAN TO SNCAN (RESERVED NAME
C     *                         IN F90).
C     * NOV 04/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.
C     * JUL 08/04 - D.VERSEGHY. NEW LOWER LIMITS FOR RCAN, SCAN, ZPOND
C     *                         AND SNOW.
C     * DEC 09/02 - D.VERSEGHY. SWITCH CALLING ORDER OF TFREEZ AND
C     *                         SNOVAP FOR CONSISTENCY WITH DIAGNOSTICS.
C     * SEP 26.02 - D.VERSEGHY. CHANGED CALL TO SUBCAN.
C     * AUG 01/02 - D.VERSEGHY. ADD CALL TO WATROF, NEW SUBROUTINE
C     *                         CONTAINING WATERLOO OVERLAND FLOW
C     *                         AND INTERFLOW CALCULATIONS.
C     *                         SHORTENED CLASS3 COMMON BLOCK.
C     * JUL 03/02 - D.VERSEGHY. STREAMLINE SUBROUTINE CALLS; MOVE
C     *                         CALCULATION OF BACKGROUND SOIL
C     *                         PROPERTIES INTO "CLASSB"; CHANGE
C     *                         RHOSNI FROM CONSTANT TO VARIABLE.
C     * OCT 04/01 - M.LAZARE.   REMOVE SEVERAL OLD DIAGNOSTIC FIELDS
C     *                         AND ADD NEW FIELD "ROVG".
C     * MAY 14/01 - M.LAZARE.   ADD CALLS TO SUBROUTINE "SNOVAP" FOR
C     *                         FC AND FG SUBAREAS OF GRID CELL.
C     * OCT 20/00 - D.VERSEGHY. ADD WORK ARRAY "RHOMAX" FOR SNOALBW.
C     * JUN 20/97 - D.VERSEGHY. CLASS - VERSION 2.7.
C     *                         CHANGES RELATED TO VARIABLE SOIL DEPTH
C     *                         (MOISTURE HOLDING CAPACITY) AND DEPTH-
C     *                         VARYING SOIL PROPERTIES.
C     * JAN 02/96 - D.VERSEGHY. CLASS - VERSION 2.5.
C     *                         COMPLETION OF ENERGY BALANCE
C     *                         DIAGNOSTICS; INTRODUCE CALCULATION OF
C     *                         OVERLAND FLOW.
C     * AUG 30/95 - D.VERSEGHY. CLASS - VERSION 2.4.
C     *                         VARIABLE SURFACE DETENTION CAPACITY
C     *                         IMPLEMENTED.
C     * AUG 24/95 - D.VERSEGHY. UPDATE ARRAY "EVAP" TO TAKE INTO
C     *                         ACCOUNT "WLOST"; RATIONALIZE
C     *                         CALCULATION OF THE LATTER.
C     *                         COMPLETION OF WATER BUDGET DIAGNOSTICS.
C     * AUG 18/95 - D.VERSEGHY. REVISIONS TO ALLOW FOR INHOMOGENEITY
C     *                         BETWEEN SOIL LAYERS AND FRACTIONAL
C     *                         ORGANIC MATTER CONTENT.
C     * DEC 22/94 - D.VERSEGHY. CLASS - VERSION 2.3.
C     *                         CHANGES TO SUBROUTINE CALLS ASSOCIATED
C     *                         WITH REVISIONS TO DIAGNOSTICS.
C     *                         ALLOW SPECIFICATION OF LIMITING POND
C     *                         DEPTH "PNDLIM" (PARALLEL CHANGES MADE
C     *                         SIMULTANEOUSLY IN TMCALC).
C     * DEC 16/94 - D.VERSEGHY. TWO NEW DIAGNOSTIC FIELDS.
C     * NOV 18/93 - D.VERSEGHY. LOCAL VERSION WITH INTERNAL WORK ARRAYS
C     *                         HARD-CODED FOR USE ON PCS.
C     * NOV 01/93 - D.VERSEGHY. CLASS - VERSION 2.2.
C     *                         REVISIONS ASSOCIATED WITH NEW VERSION
C     *                         OF TMCALC.
C     * JUL 30/93 - D.VERSEGHY/M.LAZARE. NUMEROUS NEW DIAGNOSTIC FIELDS.
C     * MAY 06/93 - D.VERSEGHY/M.LAZARE. CORRECT BUG IN CALL TO TMCALC
C     *                                  FOR CANOPY-SNOW CASE, WHERE
C     *                                  SHOULD BE PASSING "HCPCS"
C     *                                  INSTEAD OF "HCPGS".
C     * MAY 15/92 - D.VERSEGHY/M.LAZARE. REVISED AND VECTORIZED CODE
C     *                                  FOR MODEL VERSION GCM7.
C     * AUG 12/91 - D.VERSEGHY. CODE FOR MODEL VERSION GCM7U -
C                               CLASS VERSION 2.0 (WITH CANOPY).
C     * APR 11/89 - D.VERSEGHY. LAND SURFACE WATER BUDGET CALCULATIONS.
C
      IMPLICIT NONE

      integer igwscheme
C     TEMPORARY VARS
      REAL QOUT, DMLIQG(ILG), DMLIQGS(ILG), DMLIQC(ILG), DMLIQCS(ILG),
     &     DMLIQT(ILG), xsi(ILG,IG), ZWT(ILG),
     1     DMLIQJG(ILG), DMLIQJGS(ILG), DMLIQJC(ILG), DMLIQJCS(ILG),
     2     EXCWT(ILG),EXCWCS(ILG),EXCWGS(ILG),EXCWG(ILG),EXCWC(ILG)

C     * INTEGER CONSTANTS.
C
      INTEGER IWF   !Flag governing lateral soil water flow calculations
      INTEGER ILG,IL1,IL2,JL,IC,IG,IGP1,IGP2,I,J
      INTEGER NLANDCS   !Number of modelled areas that contain subareas 
                        !of canopy over snow
      INTEGER NLANDGS   !Number of modelled areas that contain subareas 
                        !of snow
      INTEGER NLANDC    !Number of modelled areas that contain subareas 
                        !of canopy over bare ground
      INTEGER NLANDG    !Number of modelled areas that contain subareas 
                        !of bare ground
      INTEGER NLANDI    !Number of modelled areas that are ice sheets 
                        ![ ]
      INTEGER IPTBAD,JPTBAD,KPTBAD,LPTBAD,N
      integer Ibed(ilg),KQIN(ILG)
C
C     * MAIN OUTPUT FIELDS.
C                                                                                  
      REAL THLIQ (ILG,IG)   !Volumetric liquid water content of soil 
                            !layers [m3 m-3]
      REAL THICE (ILG,IG)   !Volumetric frozen water content of soil 
                            !layers [m3 m-3]
      REAL TBAR  (ILG,IG)   !Temperature of soil layers [K]
      REAL GFLUX (ILG,IG)   !Heat flux at interfaces between soil layers 
                            ![W m-2]
C
      REAL TCAN  (ILG)  !Vegetation canopy temperature [K]   
      REAL RCAN  (ILG)  !Intercepted liquid water stored on canopy 
                        ![kg m-2]
      REAL SNCAN (ILG)  !Intercepted frozen water stored on canopy 
                        ![kg m-2]
      REAL RUNOFF(ILG)  !Total runoff from soil [m or kg m-2 s-1]
      REAL SNO   (ILG)  !Mass of snow pack [kg m-2]  
      REAL TSNOW (ILG)  !Snowpack temperature [K]  
      REAL RHOSNO(ILG)  !Density of snow [kg m-3]  
      REAL ALBSNO(ILG)  !Snow albedo [ ]
      REAL ZPOND (ILG)  !Depth of ponded water on surface [m]
      REAL TPOND (ILG)  !Temperature of ponded water [K]
      REAL GROWTH(ILG)  !Vegetation growth index [ ]  
      REAL TBASE (ILG)  !Temperature of bedrock in third soil layer [K]
      REAL TRUNOF(ILG)  !Temperature of total runoff [K]  
      REAL WSNOW (ILG)  !Liquid water content of snow pack [kg m-2]
C
      real QDIS(ILG),      QING(ILG),
     4     QINCS(ILG),     QINGS(ILG),     QINC(ILG),      QDISG(ILG),
     5     QDISCS(ILG),    QDISGS(ILG),    QDISC(ILG),   QDIST(ILG), 
     6     QINT(ILG),      GWTPC(ILG),     GWTPG(ILG),   GWTPCS(ILG),
     7     GWTPGS(ILG)
C
      REAL INTFL(ILG, IG) !HUZIY
C
C     * DIAGNOSTIC ARRAYS.
C
      REAL PCFC  (ILG)  !Frozen precipitation intercepted by vegetation 
                        ![kg m-2 s-1]
      REAL PCLC  (ILG)  !Liquid precipitation intercepted by vegetation 
                        ![kg m-2 s-1]
      REAL PCPN  (ILG)  !Precipitation incident on snow pack 
                        ![kg m-2 s-1]
      REAL PCPG  (ILG)  !Precipitation incident on ground [kg m-2 s-1]  
      REAL QFCF  (ILG)  !Sublimation from frozen water on vegetation 
                        ![kg m-2 s-1]
      REAL QFCL  (ILG)  !Evaporation from liquid water on vegetation 
                        ![kg m-2 s-1]
      REAL QFN   (ILG)  !Sublimation from snow pack [kg m-2 s-1]  
      REAL QFG   (ILG)  !Evaporation from ground [kg m-2 s-1]  
      REAL HMFC  (ILG)  !Diagnosed energy associated with phase change 
                        !of water on vegetation [W m-2]
      REAL HMFN  (ILG)  !Diagnosed energy associated with phase change 
                        !of water in snow pack [W m-2]
      REAL HTCC  (ILG)  !Diagnosed internal energy change of vegetation 
                        !canopy due to conduction and/or change in mass 
                        ![W m -2]
      REAL HTCS  (ILG)  !Diagnosed internal energy change of snow pack 
                        !due to conduction and/or change in mass [W m-2]
      REAL ROFC  (ILG)  !Liquid/frozen water runoff from vegetation 
                        ![kg m-2 s-1]
      REAL ROFN  (ILG)  !Liquid water runoff from snow pack [kg m-2 s-1]  
      REAL ROVG  (ILG)  !Liquid/frozen water runoff from vegetation to 
                        !ground surface [kg m-2 s-1]
      REAL WTRS  (ILG)  !Diagnosed residual water transferred into or 
                        !out of the snow pack [kg m-2 s-1]
      REAL WTRG  (ILG)  !Diagnosed residual water transferred into or 
                        !out of the soil [kg m-2 s-1]
      REAL OVRFLW(ILG)  !Overland flow from top of soil column 
                        ![m or kg m-2 s-1]  
      REAL SUBFLW(ILG)  !Interflow from sides of soil column 
                        ![kg m-2 s-1]  
      REAL BASFLW(ILG)  !Base flow from bottom of soil column 
                        ![m or kg m-2 s-1]
      REAL TOVRFL(ILG)  !Temperature of overland flow from top of soil 
                        !column [K]
      REAL TSUBFL(ILG)  !Temperature of interflow from sides of soil 
                        !column [K]
      REAL TBASFL(ILG)  !Temperature of base flow from bottom of soil 
                        !column [K]
      REAL EVAP  (ILG)  !Diagnosed total surface water vapour flux over 
                        !modelled area [kg m -2 s-1]
C
      REAL QFC   (ILG,IG)   !Water removed from soil layers by 
                            !transpiration [kg m-2 s-1]
      REAL HMFG  (ILG,IG)   !Diagnosed energy associated with phase 
                            !change of water in soil layers [W m-2]
      REAL HTC   (ILG,IG)   !Diagnosed internal energy change of soil 
                            !layer due to conduction and/or change in 
                            !mass [W m-2]
C
C     * I/O FIELDS PASSED THROUGH CLASS.
C
      !(In composite definitions, suffix C or CO = vegetation over 
      !ground; G or GO = bare ground; CS = vegetation over snow cover; 
      !GS = bare snow cover.)
C
      REAL RPCP  (ILG)  !Rainfall rate over modelled area [m s-1]
      REAL TRPCP (ILG)  !Rainfall temperature over modelled area [C] 
      REAL SPCP  (ILG)  !Snowfall rate over modelled area [m s-1] 
      REAL TSPCP (ILG)  !Snowfall temperature over modelled area [C]
      REAL PCPR  (ILG)  !Surface precipitation rate [kg m-2 s-1] 
      REAL TA    (ILG)  !Air temperature at reference height [K]
C
      REAL TBARC(ILG,IG), TBARG(ILG,IG), TBARCS(ILG,IG),TBARGS(ILG,IG)
           !Subarea temperatures of soil layers [C]
C                            
      REAL THLIQC(ILG,IG)   !Liquid water content of soil layers under 
                            !vegetation [m3 m-3]
      REAL THLIQG(ILG,IG)   !Liquid water content of soil layers in bare 
                            !areas [m3 m-3]
      REAL THICEC(ILG,IG)   !Frozen water content of soil layers under 
                            !vegetation [m3 m-3]
      REAL THICEG(ILG,IG)   !Frozen water content of soil layers in bare 
                            !areas [m3 m-3]                
      REAL HCPC  (ILG,IG)   !Heat capacity of soil layers under 
                            !vegetation [J m-3 K-1]
      REAL HCPG  (ILG,IG)   !Heat capacity of soil layers in bare areas 
                            ![J m-3 K-1]
      REAL TCTOPC(ILG,IG)   !Thermal conductivity of soil at top of 
                            !layer (vegetation over ground) [W m-1 K-1]
      REAL TCBOTC(ILG,IG)   !Thermal conductivity of soil at bottom of 
                            !layer (vegetation over ground) [W m-1 K-1]    
      REAL TCTOPG(ILG,IG)   !Thermal conductivity of soil at top of 
                            !layer (bare ground) [W m-1 K-1]
      REAL TCBOTG(ILG,IG)   !Thermal conductivity of soil at bottom of 
                            !layer (bare ground) [W m-1 K-1]
      REAL FROOT (ILG,IG)   !Fraction of total transpiration contributed 
                            !by soil layer [ ]
      REAL TSFSAV(ILG,4)    !Ground surface temperature over subarea [K]
C
      REAL FC    (ILG),   FG    (ILG),   FCS   (ILG),   FGS   (ILG)
           !Subarea fractional coverage of modelled area [ ]
C
      REAL TPONDC(ILG),   TPONDG(ILG),   TPNDCS(ILG),   TPNDGS(ILG)
           !Subarea temperature of surface ponded water [C]
C          
      REAL EVAPC (ILG)  !Evaporation from vegetation over ground [m s-1]   
      REAL EVAPCG(ILG)  !Evaporation from ground under vegetation 
                        ![m s-1]
      REAL EVAPG (ILG)  !Evaporation from bare ground [m s-1] 
      REAL EVAPCS(ILG)  !Evaporation from vegetation over snow [m s-1]             
      REAL EVPCSG(ILG)  !Evaporation from snow under vegetation [m s-1] 
      REAL EVAPGS(ILG)  !Evaporation from snow on bare ground [m s-1] 
      REAL QFREZC(ILG)  !Heat sink to be used for freezing water on 
                        !ground under canopy [W m-2]
      REAL QFREZG(ILG)  !Heat sink to be used for freezing water on bare 
                        !ground [W m-2]
      REAL QMELTC(ILG)  !Heat to be used for melting snow under canopy 
                        ![W m-2]
      REAL QMELTG(ILG)  !Heat to be used for melting snow on bare ground 
                        ![W m-2] 
      REAL RAICAN(ILG)  !Intercepted liquid water stored on canopy over 
                        !ground [kg m-2]
      REAL SNOCAN(ILG)  !Intercepted frozen water stored on canopy over 
                        !ground [kg m-2] 
      REAL RAICNS(ILG)  !Intercepted liquid water stored on canopy over 
                        !snow [kg m-2]
      REAL SNOCNS(ILG)  !Intercepted frozen water stored on canopy over 
                        !snow [kg m-2]
      REAL FSVF  (ILG)  !Sky view factor of ground under vegetation 
                        !canopy [ ]
      REAL FSVFS (ILG)  !Sky view factor of snow under vegetation canopy 
                        ![ ]
      REAL CWLCAP(ILG)  !Storage capacity of canopy over bare ground for 
                        !liquid water [kg m-2]
      REAL CWFCAP(ILG)  !Storage capacity of canopy over bare ground for 
                        !frozen water [kg m-2]
      REAL CWLCPS(ILG)  !Storage capacity of canopy over snow for liquid 
                        !water [kg m-2]
      REAL CWFCPS(ILG)  !Storage capacity of canopy over snow for frozen 
                        !water [kg m-2]
      REAL TCANO (ILG)  !Temperature of canopy over ground [K]   
      REAL TCANS (ILG)  !Temperature of canopy over snow [K] 
      REAL CHCAP (ILG)  !Heat capacity of canopy over bare ground 
                        ![J m-2 K-1]
      REAL CHCAPS(ILG)  !Heat capacity of canopy over snow [J m-2 K-1] 
      REAL CMASSC(ILG)  !Mass of canopy over bare ground [kg m-2]   
      REAL CMASCS(ILG)  !Mass of canopy over snow [kg m-2] 
      REAL ZSNOW (ILG)  !Depth of snow pack [m] 
      REAL RHOSNI(ILG)  !Density of fresh snow [kg m-3]          
      REAL GZEROC(ILG),   GZEROG(ILG),   GZROCS(ILG),   GZROGS(ILG)
           !Subarea heat flux at soil surface [W m-2]
C            
      REAL G12C  (ILG),   G12G  (ILG),   G12CS (ILG),   G12GS (ILG)
           !Subarea heat flux between first and second soil layers 
           ![W m-2]
C               
      REAL G23C  (ILG),   G23G  (ILG),   G23CS (ILG),   G23GS (ILG)
           !Subarea heat flux between second and third soil layers 
           ![W m-2]
C
      REAL TSNOCS(ILG)  !Temperature of snow pack under vegetation [K]   
      REAL TSNOGS(ILG)  !Temperature of snow pack in bare areas [K] 
      REAL WSNOCS(ILG)  !Liquid water content of snow pack under 
                        !vegetation [kg m-2]
      REAL WSNOGS(ILG)  !Liquid water content of snow pack in bare areas 
                        ![kg m-2]
      REAL RHOSCS(ILG)  !Density of snow under vegetation [kg m-3]   
      REAL RHOSGS(ILG)  !Density of snow in bare areas [kg m-3] 
      REAL ZPLIMC(ILG)  !Subarea maximum ponding depth [m] 
      REAL ZPLIMG(ILG)  !Subarea maximum ponding depth [m]
      REAL ZPLMCS(ILG)  !Subarea maximum ponding depth [m]
      REAL ZPLMGS(ILG)  !Subarea maximum ponding depth [m] 
      REAL GGEO  (ILG)  !Geothermal heat flux at bottom of soil profile 
                        ![W m-2]
C
C     * SOIL PROPERTY ARRAYS.
C
      REAL THPOR (ILG,IG)   !Pore volume in soil layer [m3 m-3]
      REAL THLRET(ILG,IG)   !Liquid water retention capacity for organic 
                            !soil [m3 m-3 ]
      REAL THLMIN(ILG,IG)   !Residual soil liquid water content 
                            !remaining after freezing or evaporation 
                            ![m3 m-3]
      REAL BI    (ILG,IG)   !Clapp and Hornberger empirical “b” 
                            !parameter [ ]
      REAL GRKSAT(ILG,IG)   !Saturated hydraulic conductivity of soil 
                            !layer [m s-1]
      REAL PSISAT(ILG,IG)   !Soil moisture suction at saturation [m]
      REAL THLRAT(ILG,IG)   !Fractional saturation of soil behind the 
                            !wetting front [ ]
      REAL THFC  (ILG,IG)   !Field capacity [m3 m-3]
      REAL HCPS  (ILG,IG)   !Heat capacity of soil material [J m-3 K-1]
      REAL DELZW (ILG,IG)   !Overall thickness of soil layer [m]
      REAL DELZZ (ILG,IG)   !Permeable thickness of soil layer [m]
      REAL ZBOTW (ILG,IG)   !Depth to permeable bottom of soil layer [m]
      REAL XDRAIN(ILG)      !Drainage index at bottom of soil profile 
                            ![ ]   
      REAL XSLOPE(ILG)      !Surface slope (used when running MESH code) 
                            ![degrees]
      REAL GRKFAC(ILG)      !WATROF parameter used when running MESH 
                            !code [ ]
      REAL WFSURF(ILG)      !WATROF parameter used when running MESH 
                            !code [ ]
      REAL WFCINT(ILG)      !WATROF parameter used when running MESH 
                            !code [ ]
      REAL DELZ  (IG)       !Overall thickness of soil layer [m]
      real XSANI(ILG),ARE(ILG),SLP(ILG),    LEG(ILG)
C
      INTEGER             ISAND(ILG,IG) !Sand content flag
      INTEGER             IGDR  (ILG)   !Index of soil layer in which 
                                        !bedrock is encountered
C
C     * INTERNAL WORK ARRAYS USED THROUGHOUT CLASSW.
C
      REAL TBARWC(ILG,IG),TBARWG(ILG,IG),TBRWCS(ILG,IG),TBRWGS(ILG,IG),
     1     THLQCO(ILG,IG),THLQGO(ILG,IG),THLQCS(ILG,IG),THLQGS(ILG,IG),
     2     THICCO(ILG,IG),THICGO(ILG,IG),THICCS(ILG,IG),THICGS(ILG,IG),
     3     HCPCO (ILG,IG),HCPGO (ILG,IG),HCPCS (ILG,IG),HCPGS (ILG,IG),
     4     GRKSC (ILG,IG),GRKSG (ILG,IG),GRKSCS(ILG,IG),GRKSGS(ILG,IG),
     5     GFLXC (ILG,IG),GFLXG (ILG,IG),GFLXCS(ILG,IG),GFLXGS(ILG,IG)
C
      REAL SPCC  (ILG),   SPCG  (ILG),   SPCCS (ILG),   SPCGS (ILG),
     1     TSPCC (ILG),   TSPCG (ILG),   TSPCCS(ILG),   TSPCGS(ILG),
     2     RPCC  (ILG),   RPCG  (ILG),   RPCCS (ILG),   RPCGS (ILG),
     3     TRPCC (ILG),   TRPCG (ILG),   TRPCCS(ILG),   TRPCGS(ILG),
     4     EVPIC (ILG),   EVPIG (ILG),   EVPICS(ILG),   EVPIGS(ILG),
     5     ZPONDC(ILG),   ZPONDG(ILG),   ZPNDCS(ILG),   ZPNDGS(ILG),
     6     XSNOWC(ILG),   XSNOWG(ILG),   XSNOCS(ILG),   XSNOGS(ILG),
     7     ZSNOWC(ILG),   ZSNOWG(ILG),   ZSNOCS(ILG),   ZSNOGS(ILG),
     8     ALBSC (ILG),   ALBSG (ILG),   ALBSCS(ILG),   ALBSGS(ILG),
     9     RHOSC (ILG),   RHOSG (ILG),
     A     HCPSC (ILG),   HCPSG (ILG),   HCPSCS(ILG),   HCPSGS(ILG),
     B     RUNFC (ILG),   RUNFG (ILG),   RUNFCS(ILG),   RUNFGS(ILG),
     C     TRUNFC(ILG),   TRUNFG(ILG),   TRNFCS(ILG),   TRNFGS(ILG),
     D     TBASC (ILG),   TBASG (ILG),   TBASCS(ILG),   TBASGS(ILG)
C

!Huziy (interflows for different subareas)
      REAL INTFL_C(ILG, IG), INTFL_G(ILG, IG),
     1     INTFL_CS(ILG, IG), INTFL_GS(ILG, IG)

      REAL DD(ILG)
      REAL SANI(ILG)
      LOGICAL INTF_FLAG
      parameter (INTF_FLAG=.false.)


      REAL SUBLC (ILG),   SUBLCS(ILG),   WLOSTC(ILG),   WLOSTG(ILG),
     1     WLSTCS(ILG),   WLSTGS(ILG),   RAC   (ILG),   RACS  (ILG),
     2     SNC   (ILG),   SNCS  (ILG),   TSNOWC(ILG),   TSNOWG(ILG),
     3     DT    (ILG),   ZERO  (ILG),   RALB  (ILG),   ZFAV  (ILG),
     4     THLINV(ILG)
C
      INTEGER             LZFAV (ILG)
C
C     * INTERNAL WORK ARRAYS FOR WPREP AND CANADD.
C
      REAL RADD  (ILG),    SADD  (ILG)
C
C     * INTERNAL WORK FIELDS FOR GRINFL/GRDRAN (AND THEIR CALLED
C     * ROUTINES (I.E. WFILL,WFLOW,WEND) AND ICEBAL.
C
      REAL ZMAT  (ILG,IGP2,IGP1)
C
      REAL WMOVE (ILG,IGP2),   TMOVE (ILG,IGP2)
C
      REAL THLIQX(ILG,IGP1),   THICEX(ILG,IGP1),   TBARWX(ILG,IGP1),
     1     DELZX (ILG,IGP1),   ZBOTX (ILG,IGP1),   FDT   (ILG,IGP1),
     2     TFDT  (ILG,IGP1),   PSIF  (ILG,IGP1),   THLINF(ILG,IGP1),
     3     GRKINF(ILG,IGP1),   FDUMMY(ILG,IGP1),   TDUMMY(ILG,IGP1),
     4     ZRMDR (ILG,IGP1)
C
      REAL THLMAX(ILG,IG),     THTEST(ILG,IG),     THLDUM(ILG,IG),
     1     THIDUM(ILG,IG),     TDUMW (ILG,IG)
C
      REAL TRMDR (ILG),    ZF    (ILG),    FMAX  (ILG),    TUSED (ILG),
     1     RDUMMY(ILG),    WEXCES(ILG),    FDTBND(ILG),    WADD  (ILG),
     2     TADD  (ILG),    WADJ  (ILG),    TIMPND(ILG),    DZF   (ILG),
     3     DTFLOW(ILG),    THLNLZ(ILG),    THLQLZ(ILG),    DZDISP(ILG),
     4     WDISP (ILG),    WABS  (ILG),    ZMOVE (ILG),    TBOT  (ILG)
C
      real WTnew(ILG),WTnewC(ILG),WTnewCS(ILG),WTnewG(ILG),WTnewGS(ILG)
C
      INTEGER              IGRN  (ILG),    IGRD  (ILG),    IZERO (ILG),
     1                     IFILL (ILG),    LZF   (ILG),    NINF  (ILG),
     2                     IFIND (ILG),    ITER  (ILG),    NEND  (ILG),
     3                     ISIMP (ILG),    ICONT (ILG)
      integer Isat(ILG),   Isats(ILG)
C
C     * INTERNAL WORK ARRAYS FOR CANVAP AND SNOALBW.
C
      REAL EVLOST(ILG),    RLOST (ILG),    RHOMAX(ILG)
C
      INTEGER              IROOT (ILG)
C
C     * INTERNAL WORK ARRAYS FOR WATROF.
C
      REAL THCRIT(ILG,IG), DODRN (ILG),     DOVER (ILG),
     1     DIDRN (ILG,IG), DIDRNMX(ILG,IG)
C
C     * INTERNAL WORK ARRAYS FOR CHKWAT.
C
      REAL BAL   (ILG)
      real WT(ILG),WTC(ILG),WTCS(ILG),WTG(ILG), WTGS(ILG)
C
C     * COMMON BLOCK PARAMETERS.
C
      REAL DELT     !Time step [s]
      REAL TFREZ    !Freezing point of water [K]
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
      COMMON /CLASS1/ DELT,TFREZ
      COMMON /CLASS3/ TCW,TCICE,TCSAND,TCCLAY,TCOM,TCDRYS,
     1                RHOSOL,RHOOM
      COMMON /CLASS4/ HCPW,HCPICE,HCPSOL,HCPOM,HCPSND,HCPCLY,
     1                SPHW,SPHICE,SPHVEG,SPHAIR,RHOW,RHOICE,
     2                TCGLAC,CLHMLT,CLHVAP
C
C-----------------------------------------------------------------------
      !
      !First, subroutine WPREP is called to initialize various arrays 
      !and produce parameters for the four subareas of canopy over snow 
      !(CS), snow on ground (GS), canopy over ground (C) and bare ground 
      !(G). Then, for each of the four subareas, if the number of 
      !modelled areas containing that subarea is greater than zero, a 
      !series of subroutines is called. The subroutines associated with 
      !each subarea are listed in the table below.
      !
      !----------------------------------------------------------------|
      ! CANVAP  | Evaporation/sublimation of water from    |   CS,C    |
      !         | vegetation canopy                        |           |
      !----------------------------------------------------------------|
      ! CANADD  | Addition of rainfall/snowfall to canopy; |   CS,C    |
      !         | throughfall and drip                     |           |
      !----------------------------------------------------------------|
      ! CWCALC  | Freezing/thawing of liquid/frozen water  |   CS,C    |
      !         | on canopy                                |           |
      !----------------------------------------------------------------|
      ! SUBCAN  | Precipitaiton and condensation under     |   CS,C    |
      !         | canopy                                   |           |
      !----------------------------------------------------------------|
      ! TWCALC  | Freezing/thawing of liquid/frozen water  | CS,GS,C,G |
      !         | in soil                                  |           |
      !----------------------------------------------------------------|
      ! SNOVAP  | Sublimaiton from snow pack               | CS,GS,C,G |
      !         |                                          |           |
      !----------------------------------------------------------------|
      ! TFREEZ  | Freezing of ponded water on soil         | CS,GS,C,G |
      !         |                                          |           |
      !----------------------------------------------------------------|
      ! TMELT   | Melting of snow pack                     |   CS,GS   |
      !         |                                          |           |
      !----------------------------------------------------------------|
      ! SNOADD  | Accumulation of snow on ground           | CS,GS,C,G |
      !         |                                          |           |
      !----------------------------------------------------------------|
      ! SNINFL  | Infiltration of rain into snow pack      |   CS,GS   |
      !         |                                          |           |
      !----------------------------------------------------------------|
      ! ICEBAL  | Energy and water budget of ice sheets    |   GS,G    |
      !         |                                          |           |
      !----------------------------------------------------------------|
      ! GRINFL  | Infiltraiton of water into soil          | CS,GS,C,G |
      !         |                                          |           |
      !----------------------------------------------------------------|
      ! GRDRAN  | Soil water movement in response to       | CS,GS,C,G |
      !         | gravity and suction forces               |           |
      !----------------------------------------------------------------|
      ! TMCALC  | Step ahead soil layer temperatures,      | CS,GS,C,G |
      !         | check for freezing/thawing               |           |
      !----------------------------------------------------------------|
      ! CHKWAT  | Check subarea moisture balances for      | CS,GS,C,G |
      !         | closure                                  |           |
      !----------------------------------------------------------------|
      ! SNOALBW | Temporal variation of snow albedo and    |   CS,GS   |
      !         | density                                  |           |
      !----------------------------------------------------------------| 
      !
C     * PREPARATION.
C
      CALL WPREP(THLQCO, THLQGO, THLQCS, THLQGS, THICCO, THICGO,
     1           THICCS, THICGS, HCPCO,  HCPGO,  HCPCS,  HCPGS,
     2           GRKSC,  GRKSG,  GRKSCS, GRKSGS,
     3           SPCC,   SPCG,   SPCCS,  SPCGS,  TSPCC,  TSPCG,
     4           TSPCCS, TSPCGS, RPCC,   RPCG,   RPCCS,  RPCGS,
     5           TRPCC,  TRPCG,  TRPCCS, TRPCGS, EVPIC,  EVPIG,
     6           EVPICS, EVPIGS, ZPONDC, ZPONDG, ZPNDCS, ZPNDGS,
     7           XSNOWC, XSNOWG, XSNOCS, XSNOGS, ZSNOWC, ZSNOWG,
     8           ZSNOCS, ZSNOGS, ALBSC,  ALBSG,  ALBSCS, ALBSGS,
     9           RHOSC,  RHOSG,  HCPSC,  HCPSG,  HCPSCS, HCPSGS,
     A           RUNFC,  RUNFG,  RUNFCS, RUNFGS,
     B           TRUNFC, TRUNFG, TRNFCS, TRNFGS, TBASC,  TBASG,
     C           TBASCS, TBASGS, GFLXC,  GFLXG,  GFLXCS, GFLXGS,
     D           SUBLC,  SUBLCS, WLOSTC, WLOSTG, WLSTCS, WLSTGS,
     E           RAC,    RACS,   SNC,    SNCS,   TSNOWC, TSNOWG,
     F           OVRFLW, SUBFLW, BASFLW, TOVRFL, TSUBFL, TBASFL,
     G           PCFC,   PCLC,   PCPN,   PCPG,   QFCF,   QFCL,
     H           QFN,    QFG,    QFC,    HMFG,
     I           ROVG,   ROFC,   ROFN,   TRUNOF,
     J           THLIQX, THICEX, THLDUM, THIDUM,
     K           DT,     RDUMMY, ZERO,   IZERO,  DELZZ,
     L           FC,     FG,     FCS,    FGS,
     M           THLIQC, THLIQG, THICEC, THICEG, HCPC,   HCPG,
     N           TBARC,  TBARG,  TBARCS, TBARGS, TBASE,  TSFSAV,
     O           FSVF,   FSVFS,  RAICAN, SNOCAN, RAICNS, SNOCNS,
     P           EVAPC,  EVAPCG, EVAPG,  EVAPCS, EVPCSG, EVAPGS,
     Q           RPCP,   TRPCP,  SPCP,   TSPCP,  RHOSNI, ZPOND,
     R           ZSNOW,  ALBSNO, WSNOCS, WSNOGS, RHOSCS, RHOSGS,
     S           THPOR,  HCPS,   GRKSAT, ISAND,  DELZW,  DELZ,
     T           ILG,    IL1,    IL2,    JL,     IG,     IGP1,
     U           NLANDCS,NLANDGS,NLANDC, NLANDG, RADD,   SADD,
     V           WT,WTGS,WTG,WTC,WTCS,WTnew,WTnewGS,WTnewG,
     W           WTnewC,WTnewCS,DMLIQC,DMLIQCS,DMLIQG,
     X           DMLIQGS,KQIN,GWTPC,GWTPG,GWTPCS,GWTPGS,
     Y           EXCWCS,EXCWGS,EXCWC,EXCWG,QINCS,QINGS,QINC,QING,
     Z           qdist,qdiscs,qdisgs,qdisc,qdisg,
     *           DMLIQJG,DMLIQJGS, DMLIQJC, DMLIQJCS,xsi,
     +           igwscheme)



      !Huziy initialize interflow
      IF (INTF_FLAG) THEN
        INTFL_C = 0
        INTFL_CS = 0
        INTFL_G = 0
        INTFL_GS = 0
      ENDIF

C
C
C     * CALCULATIONS FOR CANOPY OVER SNOW.
C
      IF(NLANDCS.GT.0)                                              THEN
          CALL CANVAP(EVAPCS,SUBLCS,RAICNS,SNOCNS,TCANS,THLQCS,
     1                TBARCS,ZSNOCS,WLSTCS,CHCAPS,QFCF,QFCL,QFN,QFC,
     2                HTCC,HTCS,HTC,FCS,CMASCS,TSNOCS,HCPSCS,RHOSCS,
     3                FROOT,THPOR,THLMIN,DELZW,EVLOST,RLOST,IROOT,
     4                IG,ILG,IL1,IL2,JL,N  )
          CALL CANADD(2,RPCCS,TRPCCS,SPCCS,TSPCCS,RAICNS,SNOCNS,
     1                TCANS,CHCAPS,HTCC,ROFC,ROVG,PCPN,PCPG,
     2                FCS,FSVFS,CWLCPS,CWFCPS,CMASCS,RHOSNI,
     3                TSFSAV(1,1),RADD,SADD,ILG,IL1,IL2)
          CALL CWCALC(TCANS,RAICNS,SNOCNS,RDUMMY,RDUMMY,CHCAPS,
     1                HMFC,HTCC,FCS,CMASCS,ILG,IL1,IL2,JL)
          CALL SUBCAN(2,RPCCS,TRPCCS,SPCCS,TSPCCS,RHOSNI,EVPCSG,
     1                QFN,QFG,PCPN,PCPG,FCS,ILG,IL1,IL2,JL)
          CALL TWCALC(TBARCS,THLQCS,THICCS,HCPCS,TBRWCS,HMFG,HTC,
     1                FCS,ZERO,THPOR,THLMIN,HCPS,DELZW,DELZZ,ISAND,
     2                IG,ILG,IL1,IL2,JL)
          CALL SNOVAP(RHOSCS,ZSNOCS,HCPSCS,TSNOCS,EVPCSG,QFN,QFG,
     1                HTCS,WLSTCS,TRNFCS,RUNFCS,TOVRFL,OVRFLW,
     2                FCS,RPCCS,SPCCS,RHOSNI,WSNOCS,ILG,IL1,IL2,JL)
          CALL TFREEZ(ZPNDCS,TPNDCS,ZSNOCS,TSNOCS,ALBSCS,
     1                RHOSCS,HCPSCS,GZROCS,HMFG,HTCS,HTC,
     2                WTRS,WTRG,FCS,ZERO,WSNOCS,TA,TBARCS,
     3                ISAND,IG,ILG,IL1,IL2,JL)
          CALL TMELT(ZSNOCS,TSNOCS,QMELTC,RPCCS,TRPCCS,
     1               GZROCS,RALB,HMFN,HTCS,HTC,FCS,HCPSCS,
     2               RHOSCS,WSNOCS,ISAND,IG,ILG,IL1,IL2,JL)
          CALL SNOADD(ALBSCS,TSNOCS,RHOSCS,ZSNOCS,
     1                HCPSCS,HTCS,FCS,SPCCS,TSPCCS,RHOSNI,WSNOCS,
     2                ILG,IL1,IL2,JL)
          CALL SNINFL(RPCCS,TRPCCS,ZSNOCS,TSNOCS,RHOSCS,HCPSCS,
     1                WSNOCS,HTCS,HMFN,PCPG,ROFN,FCS,ILG,IL1,IL2,JL)
          CALL GRINFL(1,THLQCS,THICCS,TBRWCS,BASFLW,TBASFL,RUNFCS,
     1                TRNFCS,ZFAV,LZFAV,THLINV,QFG,WLSTCS,
     2                FCS,EVPCSG,RPCCS,TRPCCS,TPNDCS,ZPNDCS,
     3                DT,ZMAT,WMOVE,TMOVE,THLIQX,THICEX,TBARWX,
     4                DELZX,ZBOTX,FDT,TFDT,PSIF,THLINF,GRKINF,
     5                THLMAX,THTEST,ZRMDR,FDUMMY,TDUMMY,THLDUM,
     6                THIDUM,TDUMW,TRMDR,ZF,FMAX,TUSED,RDUMMY,
     7                ZERO,WEXCES,FDTBND,WADD,TADD,WADJ,TIMPND,
     8                DZF,DTFLOW,THLNLZ,THLQLZ,DZDISP,WDISP,WABS,
     9                THPOR,THLRET,THLMIN,BI,PSISAT,GRKSCS,
     A                THLRAT,THFC,DELZW,ZBOTW,XDRAIN,DELZ,ISAND,
     B                IGRN,IGRD,IFILL,IZERO,LZF,NINF,IFIND,ITER,
     C                NEND,ISIMP,IGDR,
     D                IG,IGP1,IGP2,ILG,IL1,IL2,JL,N,WTCS,QDISCS,QINCS,
     D                Isat,WTnewCS,DMLIQCS,KQIN,Ibed,DMLIQJCS,
     E                Isats,xsi,ZWT,EXCWCS,ARE,LEG,SLP,XSANI,igwscheme)
          CALL GRDRAN(1,THLQCS,THICCS,TBRWCS,FDUMMY,TDUMMY,
     1                BASFLW,TBASFL,RUNFCS,TRNFCS,
     2                QFG,WLSTCS,FCS,EVPCSG,RPCCS,ZPNDCS,
     3                DT,WEXCES,THLMAX,THTEST,THPOR,THLRET,THLMIN,
     4                BI,PSISAT,GRKSCS,THFC,DELZW,XDRAIN,ISAND,
     5                IZERO,IGRN,IGRD,IGDR,
     6                IG,IGP1,IGP2,ILG,IL1,IL2,JL,N,
     7                WTCS, ZBOTW,QDISCS,QINCS,Isat,WTnewCS,
     8                DMLIQCS,2,KQIN,Ibed,DMLIQJCS,Isats,xsi,ZWT,EXCWCS,
     9                ARE,LEG,SLP,XSANI,igwscheme)


!Huziy

c$$$        IF (INTF_FLAG) THEN
c$$$         CALL INTERFLOW(INTFL_CS,
c$$$     &         BI, THPOR, THLMIN, THLQCS, THICCS, RUNFCS, TRNFCS,
c$$$     &         TBRWCS,
c$$$     &         DD,  GRKSCS(:,1),! horizontal conductivity = vertical cond.
c$$$     &         DELZW, IG,ILG ,IL1,IL2, THFC,!<- not sure where to take the bulk_fc from?
c$$$     &         ISAND, XSLOPE, DELT)
c$$$        ENDIF

          CALL TMCALC(TBARCS,THLQCS,THICCS,HCPCS,TPNDCS,ZPNDCS,
     1                TSNOCS,ZSNOCS,ALBSCS,RHOSCS,HCPSCS,TBASCS,
     2                OVRFLW,TOVRFL,RUNFCS,TRNFCS,HMFG,HTC,HTCS,
     3                WTRS,WTRG,FCS,TBRWCS,GZROCS,G12CS,
     4                G23CS,GGEO,TA,WSNOCS,TCTOPC,TCBOTC,GFLXCS,
     5                ZPLMCS,THPOR,THLMIN,HCPS,DELZW,DELZZ,DELZ,
     6                ISAND,IWF,IG,ILG,IL1,IL2,JL,N,Isat)
          CALL CHKWAT(1,PCPR,EVPICS,RUNFCS,WLSTCS,RAICNS,SNOCNS,
     1                RACS,SNCS,ZPNDCS,ZPOND,THLQCS,THICCS,
     2                THLIQC,THICEC,ZSNOCS,RHOSCS,XSNOCS,SNO,
     3                WSNOCS,WSNOW,FCS,FGS,FCS,BAL,THPOR,THLMIN,
     4                DELZW,ISAND,IG,ILG,IL1,IL2,JL,N,QDISCS,QINCS,Isat,
     5                WTnewCS,DMLIQCS,Ibed,GWTPCS,Isats,EXCWCS,
     6                igwscheme)
          CALL SNOALBW(ALBSCS,RHOSCS,ZSNOCS,HCPSCS,
     1                 TSNOCS,FCS,SPCCS,RALB,WSNOCS,RHOMAX,
     2                 ISAND,ILG,IG,IL1,IL2,JL)
      ENDIF
C
C     * CALCULATIONS FOR SNOW-COVERED GROUND.
C
      IF(NLANDGS.GT.0)                                              THEN
          CALL TWCALC(TBARGS,THLQGS,THICGS,HCPGS,TBRWGS,HMFG,HTC,
     1                FGS,ZERO,THPOR,THLMIN,HCPS,DELZW,DELZZ,ISAND,
     2                IG,ILG,IL1,IL2,JL)
          CALL SNOVAP(RHOSGS,ZSNOGS,HCPSGS,TSNOGS,EVAPGS,QFN,QFG,
     1                HTCS,WLSTGS,TRNFGS,RUNFGS,TOVRFL,OVRFLW,
     2                FGS,RPCGS,SPCGS,RHOSNI,WSNOGS,ILG,IL1,IL2,JL)
          CALL TFREEZ(ZPNDGS,TPNDGS,ZSNOGS,TSNOGS,ALBSGS,
     1                RHOSGS,HCPSGS,GZROGS,HMFG,HTCS,HTC,
     2                WTRS,WTRG,FGS,ZERO,WSNOGS,TA,TBARGS,
     3                ISAND,IG,ILG,IL1,IL2,JL)
          CALL TMELT(ZSNOGS,TSNOGS,QMELTG,RPCGS,TRPCGS,
     1               GZROGS,RALB,HMFN,HTCS,HTC,FGS,HCPSGS,
     2               RHOSGS,WSNOGS,ISAND,IG,ILG,IL1,IL2,JL)
          CALL SNOADD(ALBSGS,TSNOGS,RHOSGS,ZSNOGS,
     1                HCPSGS,HTCS,FGS,SPCGS,TSPCGS,RHOSNI,WSNOGS,
     2                ILG,IL1,IL2,JL)
          CALL SNINFL(RPCGS,TRPCGS,ZSNOGS,TSNOGS,RHOSGS,HCPSGS,
     1                WSNOGS,HTCS,HMFN,PCPG,ROFN,FGS,ILG,IL1,IL2,JL)
          IF(NLANDI.NE.0)                                       THEN
              CALL ICEBAL(TBARGS,TPNDGS,ZPNDGS,TSNOGS,RHOSGS,ZSNOGS,
     1                    HCPSGS,ALBSGS,HMFG,HTCS,HTC,WTRS,WTRG,GFLXGS,
     2                    RUNFGS,TRNFGS,OVRFLW,TOVRFL,ZPLMGS,GGEO,
     3                    FGS,EVAPGS,RPCGS,TRPCGS,GZROGS,G12GS,G23GS,
     4                    HCPGS,QMELTG,WSNOGS,ZMAT,TMOVE,WMOVE,ZRMDR,
     5                    TADD,ZMOVE,TBOT,DELZ,ISAND,ICONT,
     6                    IWF,IG,IGP1,IGP2,ILG,IL1,IL2,JL,N )
          ENDIF
          CALL GRINFL(2,THLQGS,THICGS,TBRWGS,BASFLW,TBASFL,RUNFGS,
     1                TRNFGS,ZFAV,LZFAV,THLINV,QFG,WLSTGS,
     2                FGS,EVAPGS,RPCGS,TRPCGS,TPNDGS,ZPNDGS,
     3                DT,ZMAT,WMOVE,TMOVE,THLIQX,THICEX,TBARWX,
     4                DELZX,ZBOTX,FDT,TFDT,PSIF,THLINF,GRKINF,
     5                THLMAX,THTEST,ZRMDR,FDUMMY,TDUMMY,THLDUM,
     6                THIDUM,TDUMW,TRMDR,ZF,FMAX,TUSED,RDUMMY,
     7                ZERO,WEXCES,FDTBND,WADD,TADD,WADJ,TIMPND,
     8                DZF,DTFLOW,THLNLZ,THLQLZ,DZDISP,WDISP,WABS,
     9                THPOR,THLRET,THLMIN,BI,PSISAT,GRKSGS,
     A                THLRAT,THFC,DELZW,ZBOTW,XDRAIN,DELZ,ISAND,
     B                IGRN,IGRD,IFILL,IZERO,LZF,NINF,IFIND,ITER,
     C                NEND,ISIMP,IGDR,
     D                IG,IGP1,IGP2,ILG,IL1,IL2,JL,N,WTGS,QDISGS,QINGS,
     D                Isat,WTnewGS,DMLIQGS,KQIN,Ibed,DMLIQJGS,
     E                Isats,xsi,ZWT,EXCWGS,ARE,LEG,SLP,XSANI,igwscheme)
          CALL GRDRAN(2,THLQGS,THICGS,TBRWGS,FDUMMY,TDUMMY,
     1                BASFLW,TBASFL,RUNFGS,TRNFGS,
     2                QFG,WLSTGS,FGS,EVAPGS,RPCGS,ZPNDGS,
     3                DT,WEXCES,THLMAX,THTEST,THPOR,THLRET,THLMIN,
     4                BI,PSISAT,GRKSGS,THFC,DELZW,XDRAIN,ISAND,
     5                IZERO,IGRN,IGRD,IGDR,
     6                IG,IGP1,IGP2,ILG,IL1,IL2,JL,N,
     6                WTGS,ZBOTW,QDISGS,QINGS,Isat,WTnewGS,
     7                DMLIQGS,2,KQIN,Ibed,DMLIQJGS,Isats,xsi,ZWT,EXCWGS,
     8                ARE,LEG,SLP,XSANI,igwscheme)

!Huziy call interflow

c$$$        IF (INTF_FLAG) THEN
c$$$          CALL INTERFLOW(INTFL_GS,
c$$$     &         BI, THPOR, THLMIN, THLQGS, THICGS, RUNFGS, TRNFGS,
c$$$     &         TBRWGS,
c$$$     &         DD,  GRKSGS(1,:),! horizontal conductivity = vertical cond
c$$$     &         DELZW, IG,ILG,IL1,IL2, THFC,!<- not sure where to take the bulk_fc from?
c$$$     &         ISAND, XSLOPE, DELT)
c$$$        ENDIF

          CALL TMCALC(TBARGS,THLQGS,THICGS,HCPGS,TPNDGS,ZPNDGS,
     1                TSNOGS,ZSNOGS,ALBSGS,RHOSGS,HCPSGS,TBASGS,
     2                OVRFLW,TOVRFL,RUNFGS,TRNFGS,HMFG,HTC,HTCS,
     3                WTRS,WTRG,FGS,TBRWGS,GZROGS,G12GS,
     4                G23GS,GGEO,TA,WSNOGS,TCTOPG,TCBOTG,GFLXGS,
     5                ZPLMGS,THPOR,THLMIN,HCPS,DELZW,DELZZ,DELZ,
     6                ISAND,IWF,IG,ILG,IL1,IL2,JL,N,Isat)
          CALL CHKWAT(2,PCPR,EVPIGS,RUNFGS,WLSTGS,RAICNS,SNOCNS,
     1                RACS,SNCS,ZPNDGS,ZPOND,THLQGS,THICGS,
     2                THLIQG,THICEG,ZSNOGS,RHOSGS,XSNOGS,SNO,
     3                WSNOGS,WSNOW,FCS,FGS,FGS,BAL,THPOR,THLMIN,
     4                DELZW,ISAND,IG,ILG,IL1,IL2,JL,N,QDISGS,QINGS,Isat,
     5                WTnewGS,DMLIQGS,Ibed,GWTPGS,Isats,EXCWGS,
     6                igwscheme) 
          CALL SNOALBW(ALBSGS,RHOSGS,ZSNOGS,HCPSGS,
     1                 TSNOGS,FGS,SPCGS,RALB,WSNOGS,RHOMAX,
     2                 ISAND,ILG,IG,IL1,IL2,JL)
      ENDIF
C
C     * CALCULATIONS FOR CANOPY OVER BARE GROUND.
C
      IF(NLANDC.GT.0)                                               THEN
          CALL CANVAP(EVAPC,SUBLC,RAICAN,SNOCAN,TCANO,THLQCO,
     1                TBARC,ZSNOWC,WLOSTC,CHCAP,QFCF,QFCL,QFN,QFC,
     2                HTCC,HTCS,HTC,FC,CMASSC,TSNOWC,HCPSC,RHOSC,
     3                FROOT,THPOR,THLMIN,DELZW,EVLOST,RLOST,IROOT,
     4                IG,ILG,IL1,IL2,JL,N  )
          CALL CANADD(1,RPCC,TRPCC,SPCC,TSPCC,RAICAN,SNOCAN,
     1                TCANO,CHCAP,HTCC,ROFC,ROVG,PCPN,PCPG,
     2                FC,FSVF,CWLCAP,CWFCAP,CMASSC,RHOSNI,
     3                TSFSAV(1,3),RADD,SADD,ILG,IL1,IL2)    
          CALL CWCALC(TCANO,RAICAN,SNOCAN,RDUMMY,RDUMMY,CHCAP,
     1                HMFC,HTCC,FC,CMASSC,ILG,IL1,IL2,JL)
          CALL SUBCAN(1,RPCC,TRPCC,SPCC,TSPCC,RHOSNI,EVAPCG,
     1                QFN,QFG,PCPN,PCPG,FC,ILG,IL1,IL2,JL)
          CALL TWCALC(TBARC,THLQCO,THICCO,HCPCO,TBARWC,HMFG,HTC,
     1                FC,EVAPCG,THPOR,THLMIN,HCPS,DELZW,DELZZ,
     2                ISAND,IG,ILG,IL1,IL2,JL)
          CALL SNOVAP(RHOSC,ZSNOWC,HCPSC,TSNOWC,EVAPCG,QFN,QFG,
     1                HTCS,WLOSTC,TRUNFC,RUNFC,TOVRFL,OVRFLW,
     2                FC,RPCC,SPCC,RHOSNI,ZERO,ILG,IL1,IL2,JL)
          CALL TFREEZ(ZPONDC,TPONDC,ZSNOWC,TSNOWC,ALBSC,
     1                RHOSC,HCPSC,GZEROC,HMFG,HTCS,HTC,
     2                WTRS,WTRG,FC,QFREZC,ZERO,TA,TBARC,
     3                ISAND,IG,ILG,IL1,IL2,JL)
          CALL SNOADD(ALBSC,TSNOWC,RHOSC,ZSNOWC,
     1                HCPSC,HTCS,FC,SPCC,TSPCC,RHOSNI,ZERO,
     2                ILG,IL1,IL2,JL)
          CALL GRINFL(3,THLQCO,THICCO,TBARWC,BASFLW,TBASFL,RUNFC,
     1                TRUNFC,ZFAV,LZFAV,THLINV,QFG,WLOSTC,
     2                FC,EVAPCG,RPCC,TRPCC,TPONDC,ZPONDC,
     3                DT,ZMAT,WMOVE,TMOVE,THLIQX,THICEX,TBARWX,
     4                DELZX,ZBOTX,FDT,TFDT,PSIF,THLINF,GRKINF,
     5                THLMAX,THTEST,ZRMDR,FDUMMY,TDUMMY,THLDUM,
     6                THIDUM,TDUMW,TRMDR,ZF,FMAX,TUSED,RDUMMY,
     7                ZERO,WEXCES,FDTBND,WADD,TADD,WADJ,TIMPND,
     8                DZF,DTFLOW,THLNLZ,THLQLZ,DZDISP,WDISP,WABS,
     9                THPOR,THLRET,THLMIN,BI,PSISAT,GRKSC,
     A                THLRAT,THFC,DELZW,ZBOTW,XDRAIN,DELZ,ISAND,
     B                IGRN,IGRD,IFILL,IZERO,LZF,NINF,IFIND,ITER,
     C                NEND,ISIMP,IGDR,
     D                IG,IGP1,IGP2,ILG,IL1,IL2,JL,N,WTC,
     D                QDISC,QINC,Isat,WTnewC,DMLIQC,KQIN,Ibed,DMLIQJC,
     E                Isats,xsi,ZWT,EXCWC,ARE,LEG,SLP,XSANI,igwscheme)
          CALL GRDRAN(3,THLQCO,THICCO,TBARWC,FDUMMY,TDUMMY,
     1                BASFLW,TBASFL,RUNFC,TRUNFC,
     2                QFG,WLOSTC,FC,EVAPCG,RPCC,ZPONDC,
     3                DT,WEXCES,THLMAX,THTEST,THPOR,THLRET,THLMIN,
     4                BI,PSISAT,GRKSC,THFC,DELZW,XDRAIN,ISAND,
     5                IZERO,IGRN,IGRD,IGDR,
     6                IG,IGP1,IGP2,ILG,IL1,IL2,JL,N,
     6                WTC,ZBOTW,QDISC,QINC,Isat,WTnewC,
     7                DMLIQC,2,KQIN,Ibed,DMLIQJC,Isats,xsi,ZWT,EXCWC,
     8                ARE,LEG,SLP,XSANI,igwscheme)

        !Huziy call interflow

c$$$        IF (INTF_FLAG) THEN
c$$$          CALL INTERFLOW(INTFL_C,
c$$$     &          BI, THPOR, THLMIN, THLIQC, THICEC, RUNFC, TRUNFC,
c$$$     &          TBARWC,
c$$$     &          DD,  GRKSC(:,1),!horizontal conductivity = vertical cond
c$$$     &          DELZW, IG,ILG,IL1,IL2, THFC,!<- not sure where to take the bulk_fc from?
c$$$     &          ISAND, XSLOPE, DELT)
c$$$        ENDIF

          CALL TMCALC(TBARC,THLQCO,THICCO,HCPCO,TPONDC,ZPONDC,
     1                TSNOWC,ZSNOWC,ALBSC,RHOSC,HCPSC,TBASC,
     2                OVRFLW,TOVRFL,RUNFC,TRUNFC,HMFG,HTC,HTCS,
     3                WTRS,WTRG,FC,TBARWC,GZEROC,G12C,
     4                G23C,GGEO,TA,ZERO,TCTOPC,TCBOTC,GFLXC,
     5                ZPLIMC,THPOR,THLMIN,HCPS,DELZW,DELZZ,DELZ,
     6                ISAND,IWF,IG,ILG,IL1,IL2,JL,N,Isat)
          CALL CHKWAT(3,PCPR,EVPIC,RUNFC,WLOSTC,RAICAN,SNOCAN,
     1                RAC,SNC,ZPONDC,ZPOND,THLQCO,THICCO,
     2                THLIQC,THICEC,ZSNOWC,RHOSC,XSNOWC,SNO,
     3                ZERO,ZERO,FCS,FGS,FC,BAL,THPOR,THLMIN,
     4                DELZW,ISAND,IG,ILG,IL1,IL2,JL,N,QDISC,QINC,Isat,
     5                WTnewC,DMLIQC,Ibed,GWTPC,Isats,EXCWC,
     6                igwscheme) 
C
      ENDIF
C
C     * CALCULATIONS FOR BARE GROUND.
C
      IF(NLANDG.GT.0)                                               THEN
          CALL TWCALC(TBARG,THLQGO,THICGO,HCPGO,TBARWG,HMFG,HTC,
     1                FG,EVAPG,THPOR,THLMIN,HCPS,DELZW,DELZZ,
     2                ISAND,IG,ILG,IL1,IL2,JL)
          CALL SNOVAP(RHOSG,ZSNOWG,HCPSG,TSNOWG,EVAPG,QFN,QFG,
     1                HTCS,WLOSTG,TRUNFG,RUNFG,TOVRFL,OVRFLW,
     2                FG,RPCG,SPCG,RHOSNI,ZERO,ILG,IL1,IL2,JL)
          CALL TFREEZ(ZPONDG,TPONDG,ZSNOWG,TSNOWG,ALBSG,
     1                RHOSG,HCPSG,GZEROG,HMFG,HTCS,HTC,
     2                WTRS,WTRG,FG,QFREZG,ZERO,TA,TBARG,
     3                ISAND,IG,ILG,IL1,IL2,JL)
          CALL SNOADD(ALBSG,TSNOWG,RHOSG,ZSNOWG,
     1                HCPSG,HTCS,FG,SPCG,TSPCG,RHOSNI,ZERO,
     2                ILG,IL1,IL2,JL)
          IF(NLANDI.NE.0)                                       THEN
              CALL ICEBAL(TBARG,TPONDG,ZPONDG,TSNOWG,RHOSG,ZSNOWG,
     1                    HCPSG,ALBSG,HMFG,HTCS,HTC,WTRS,WTRG,GFLXG,
     2                    RUNFG,TRUNFG,OVRFLW,TOVRFL,ZPLIMG,GGEO,
     3                    FG,EVAPG,RPCG,TRPCG,GZEROG,G12G,G23G,
     4                    HCPGO,QFREZG,ZERO,ZMAT,TMOVE,WMOVE,ZRMDR,
     5                    TADD,ZMOVE,TBOT,DELZ,ISAND,ICONT,
     6                    IWF,IG,IGP1,IGP2,ILG,IL1,IL2,JL,N )
          ENDIF
          CALL GRINFL(4,THLQGO,THICGO,TBARWG,BASFLW,TBASFL,RUNFG,
     1                TRUNFG,ZFAV,LZFAV,THLINV,QFG,WLOSTG,
     2                FG,EVAPG,RPCG,TRPCG,TPONDG,ZPONDG,
     3                DT,ZMAT,WMOVE,TMOVE,THLIQX,THICEX,TBARWX,
     4                DELZX,ZBOTX,FDT,TFDT,PSIF,THLINF,GRKINF,
     5                THLMAX,THTEST,ZRMDR,FDUMMY,TDUMMY,THLDUM,
     6                THIDUM,TDUMW,TRMDR,ZF,FMAX,TUSED,RDUMMY,
     7                ZERO,WEXCES,FDTBND,WADD,TADD,WADJ,TIMPND,
     8                DZF,DTFLOW,THLNLZ,THLQLZ,DZDISP,WDISP,WABS,
     9                THPOR,THLRET,THLMIN,BI,PSISAT,GRKSG,
     A                THLRAT,THFC,DELZW,ZBOTW,XDRAIN,DELZ,ISAND,
     B                IGRN,IGRD,IFILL,IZERO,LZF,NINF,IFIND,ITER,
     C                NEND,ISIMP,IGDR,
     D                IG,IGP1,IGP2,ILG,IL1,IL2,JL,N,WTG,
     D                QDISG,QING,Isat,WTnewG,DMLIQG,KQIN,Ibed,DMLIQJG,
     E                Isats,xsi,ZWT,EXCWG,ARE,LEG,SLP,XSANI,igwscheme)
          CALL GRDRAN(4,THLQGO,THICGO,TBARWG,FDUMMY,TDUMMY,
     1                BASFLW,TBASFL,RUNFG,TRUNFG,
     2                QFG,WLOSTG,FG,EVAPG,RPCG,ZPONDG,
     3                DT,WEXCES,THLMAX,THTEST,THPOR,THLRET,THLMIN,
     4                BI,PSISAT,GRKSG,THFC,DELZW,XDRAIN,ISAND,
     5                IZERO,IGRN,IGRD,IGDR,
     6                IG,IGP1,IGP2,ILG,IL1,IL2,JL,N,
     6                WTG,ZBOTW,QDISG,QING,Isat,WTnewG,
     7                DMLIQG,2,KQIN,Ibed,DMLIQJG,Isats,xsi,ZWT,EXCWG,
     8                ARE,LEG,SLP,XSANI,igwscheme)


        !Huziy call interflow
c$$$       IF (INTF_FLAG) THEN
c$$$        CALL INTERFLOW(INTFL_G,
c$$$     &         BI, THPOR, THLMIN, THLIQG, THICEG, RUNFG, TRUNFG,
c$$$     &         TBARWG,
c$$$     &         DD, GRKSG(:,1) ,! horizontal conductivity = vertical conductivity
c$$$     &         DELZW, IG,ILG,IL1,IL2, THFC,!<- not sure where to take the bulk_fc from?
c$$$     &         ISAND, XSLOPE, DELT)
c$$$       ENDIF

          CALL TMCALC(TBARG,THLQGO,THICGO,HCPGO,TPONDG,ZPONDG,
     1                TSNOWG,ZSNOWG,ALBSG,RHOSG,HCPSG,TBASG,
     2                OVRFLW,TOVRFL,RUNFG,TRUNFG,HMFG,HTC,HTCS,
     3                WTRS,WTRG,FG,TBARWG,GZEROG,G12G,
     4                G23G,GGEO,TA,ZERO,TCTOPG,TCBOTG,GFLXG,
     5                ZPLIMG,THPOR,THLMIN,HCPS,DELZW,DELZZ,DELZ,
     6                ISAND,IWF,IG,ILG,IL1,IL2,JL,N,Isat)
          CALL CHKWAT(4,PCPR,EVPIG,RUNFG,WLOSTG,RAICAN,SNOCAN,
     1                RAC,SNC,ZPONDG,ZPOND,THLQGO,THICGO,
     2                THLIQG,THICEG,ZSNOWG,RHOSG,XSNOWG,SNO,
     3                ZERO,ZERO,FCS,FGS,FG,BAL,THPOR,THLMIN,
     4                DELZW,ISAND,IG,ILG,IL1,IL2,JL,N,QDISG,QING,Isat,
     5                WTnewG,DMLIQG,Ibed,GWTPG,Isats,EXCWG,
     6                igwscheme)
C
      ENDIF
      !
      !After these calls have been done, average values of the main 
      !prognostic variables over the modelled area are determined by 
      !performing weighted averages over the four subareas, and checks 
      !are carried out to identify and remove vanishingly small values. 
      !First the bedrock temperature in the third soil layer, the total 
      !runoff and the runoff temperature are calculated. Then the total 
      !runoff and the overland flow, interflow and baseflow are 
      !converted from units of m to kg m-2 s-1. The total surface water 
      !vapour flux over the modelled area is updated to account for the 
      !residual amounts of evaporative demand over the four subareas 
      !that could not be supplied by surface stores (WLSTCS, WLSTGS, 
      !WLOSTC and WLOSTG, variables that are defined internally in this 
      !subroutine).
      !
      !The temperature of the vegetation canopy TCAN and the amount of 
      !intercepted liquid water RCAN are calculated as weighted averages 
      !over the two canopy subareas. A flag is set to trigger a call to 
      !abort if TCAN is less than -100 C or greater than 100 C. If RCAN 
      !is vanishingly small, it is added to the overland flow and to the 
      !total runoff, and their respective temperatures are recalculated. 
      !The diagnostic arrays ROFC, ROVG, PCPG and HTCC are updated, and 
      !RCAN is set to zero. The amount of intercepted snow SNCAN is 
      !likewise calculated as a weighted average over the two canopy 
      !subareas. If SNCAN is vanishingly small, it is added to the 
      !overland flow and to the total runoff, and their respective 
      !temperatures are recalculated. The diagnostic arrays ROFC, ROVG, 
      !PCPG and HTCC are updated, and SNCAN is set to zero. If there is 
      !no canopy present, TCAN is set to zero.
      !
      !At the end of the 600 loop, the depth of ponded water ZPOND and 
      !its temperature TPOND over the modelled area are calculated as 
      !weighted averages over the four subareas. If ZPOND is vanishingly 
      !small, then as in the case of intercepted water, it is added to 
      !the overland flow and to the total runoff, and their respective 
      !temperatures are recalculated. The diagnostic array HTC is 
      !updated, and ZPOND and TPOND are set to zero.
      !
C
C     * AVERAGE RUNOFF AND PROGNOSTIC VARIABLES OVER FOUR GRID CELL
C     * SUBAREAS.
C
      JPTBAD=0
      KPTBAD=0
      LPTBAD=0
      DO 600 I=IL1,IL2
          TBASE (I)=FCS(I)*(TBASCS(I)+TFREZ) +
     1              FGS(I)*(TBASGS(I)+TFREZ) +
     2              FC (I)*(TBASC (I)+TFREZ) +
     3              FG (I)*(TBASG (I)+TFREZ)
          RUNOFF(I)=FCS(I)*RUNFCS(I) + FGS(I)*RUNFGS(I) +
     1              FC (I)*RUNFC (I) + FG (I)*RUNFG (I)

      if (igwscheme.ne.1) then
          QINT(I)=FCS(I)*QINCS(I) + FGS(I)*QINGS(I) +
     1            FC (I)*QINC (I) + FG (I)*QING (I) 
          QDIST(I)=FCS(I)*QDISCS(I) + FGS(I)*QDISGS(I) +
     1             FC (I)*QDISC (I) + FG (I)*QDISG (I)  

          DMLIQT(I)=FCS(I)*DMLIQCS(I) + FGS(I)*DMLIQGS(I) +
     1              FC (I)*DMLIQC (I) + FG (I)*DMLIQG (I)

          WT(I)=FCS(I)*WTCS(I) + FGS(I)*WTGS(I) +
     1          FC (I)*WTC (I) + FG (I)*WTG (I)

          WTnew(I)=FCS(I)*WTnewCS(I) + FGS(I)*WTnewGS(I) +
     1             FC (I)*WTnewC (I) + FG (I)*WTnewG (I)
          
          EXCWT(I)=FCS(I)*EXCWCS(I) + FGS(I)*EXCWGS(I) +
     1             FC (I)*EXCWC (I) + FG (I)*EXCWG (I)
      endif
 
          !Huziy calculate average interflow
c$$$       IF (INTF_FLAG) THEN
c$$$         INTFL(I,:) = FCS(I) * INTFL_CS(I,:) + FGS(I) * INTFL_GS(I,:) +
c$$$     1                 FC(I) * INTFL_C(I,:) + FG(I) * INTFL_G(I,:)
c$$$       ENDIF

          IF(RUNOFF(I).GT.0.0)
     1        TRUNOF(I)=(FCS(I)*RUNFCS(I)*TRNFCS(I) +
     2                   FGS(I)*RUNFGS(I)*TRNFGS(I) +
     3                   FC (I)*RUNFC (I)*TRUNFC(I) +
     4                   FG (I)*RUNFG (I)*TRUNFG(I))/RUNOFF(I)
          RUNOFF(I)=RUNOFF(I)*RHOW/DELT
          OVRFLW(I)=OVRFLW(I)*RHOW/DELT
          SUBFLW(I)=SUBFLW(I)*RHOW/DELT
          BASFLW(I)=BASFLW(I)*RHOW/DELT
          EVAP  (I)=EVAP(I)-(FCS(I)*WLSTCS(I)+FGS(I)*WLSTGS(I)+
     1              FC(I)*WLOSTC(I)+FG(I)*WLOSTG(I))/DELT
          IF((FC(I)+FCS(I)).GT.0.)                                  THEN
              TCAN(I)=(FCS(I)*TCANS(I)*CHCAPS(I)+FC(I)*TCANO(I)*
     1                CHCAP(I))/(FCS(I)*CHCAPS(I)+FC(I)*CHCAP(I))
              RCAN(I)= FCS(I)*RAICNS(I) + FC (I)*RAICAN(I)
              IF(TCAN(I).LT.173.16 .OR. TCAN(I).GT.373.16) JPTBAD=I
              IF(RCAN(I).LT.0.0) RCAN(I)=0.0
              IF(RCAN(I).LT.1.0E-5 .AND. RCAN(I).GT.0.0) THEN
                  TOVRFL(I)=(TOVRFL(I)*OVRFLW(I)+TCAN(I)*RCAN(I)/
     1                DELT)/(OVRFLW(I)+RCAN(I)/DELT)
                  OVRFLW(I)=OVRFLW(I)+RCAN(I)/DELT
                  TRUNOF(I)=(TRUNOF(I)*RUNOFF(I)+TCAN(I)*RCAN(I)/
     1                DELT)/(RUNOFF(I)+RCAN(I)/DELT)
                  RUNOFF(I)=RUNOFF(I)+RCAN(I)/DELT
                  ROFC(I)=ROFC(I)+RCAN(I)/DELT
                  ROVG(I)=ROVG(I)+RCAN(I)/DELT
                  PCPG(I)=PCPG(I)+RCAN(I)/DELT
                  HTCC(I)=HTCC(I)-TCAN(I)*SPHW*RCAN(I)/DELT
                  RCAN(I)=0.0
              ENDIF
              SNCAN  (I)=FCS(I)*SNOCNS(I) + FC (I)*SNOCAN(I)
              IF(SNCAN(I).LT.0.0) SNCAN(I)=0.0
              IF(SNCAN(I).LT.1.0E-5 .AND. SNCAN(I).GT.0.0) THEN
                  TOVRFL(I)=(TOVRFL(I)*OVRFLW(I)+TCAN(I)*SNCAN(I)/
     1                DELT)/(OVRFLW(I)+SNCAN(I)/DELT)
                  OVRFLW(I)=OVRFLW(I)+SNCAN(I)/DELT
                  TRUNOF(I)=(TRUNOF(I)*RUNOFF(I)+TCAN(I)*SNCAN(I)/
     1                DELT)/(RUNOFF(I)+SNCAN(I)/DELT)
                  RUNOFF(I)=RUNOFF(I)+SNCAN(I)/DELT
                  ROFC(I)=ROFC(I)+SNCAN(I)/DELT
                  ROVG(I)=ROVG(I)+SNCAN(I)/DELT
                  PCPG(I)=PCPG(I)+SNCAN(I)/DELT
                  HTCC(I)=HTCC(I)-TCAN(I)*SPHICE*SNCAN(I)/DELT
                  SNCAN(I)=0.0
              ENDIF
          ELSE
              TCAN(I)=0.0
          ENDIF
          IF(ZPNDCS(I).GT.0. .OR. ZPNDGS(I).GT.0. .OR.
     1                ZPONDC(I).GT.0. .OR. ZPONDG(I).GT.0.)    THEN
              ZPOND(I)=(FCS(I)*ZPNDCS(I)+FGS(I)*ZPNDGS(I)+
     1                  FC (I)*ZPONDC(I)+FG (I)*ZPONDG(I))
              TPOND(I)=(FCS(I)*(TPNDCS(I)+TFREZ)*ZPNDCS(I)+
     1                  FGS(I)*(TPNDGS(I)+TFREZ)*ZPNDGS(I)+
     2                  FC (I)*(TPONDC(I)+TFREZ)*ZPONDC(I)+
     3                  FG (I)*(TPONDG(I)+TFREZ)*ZPONDG(I))/
     4                  ZPOND(I)
              IF(ZPOND(I).LT.0.0) ZPOND(I)=0.0
              IF(ZPOND(I).LT.1.0E-8 .AND. ZPOND(I).GT.0.0) THEN
                  TOVRFL(I)=(TOVRFL(I)*OVRFLW(I)+TPOND(I)*RHOW*
     1                ZPOND(I)/DELT)/(OVRFLW(I)+RHOW*ZPOND(I)/DELT)
                  OVRFLW(I)=OVRFLW(I)+RHOW*ZPOND(I)/DELT
                  TRUNOF(I)=(TRUNOF(I)*RUNOFF(I)+TPOND(I)*RHOW*
     1                ZPOND(I)/DELT)/(RUNOFF(I)+RHOW*ZPOND(I)/DELT)
                  RUNOFF(I)=RUNOFF(I)+RHOW*ZPOND(I)/DELT
                  HTC(I,1)=HTC(I,1)-TPOND(I)*HCPW*ZPOND(I)/DELT
                  TPOND(I)=0.0
                  ZPOND(I)=0.0
              ENDIF
         ELSE
              ZPOND(I)=0.0
              TPOND(I)=0.0
         ENDIF
  600 CONTINUE
C
      !
      !In the 650 loop, values of the snow prognostic variables are 
      !calculated as weighted averages over the four subareas. The 
      !weightings for the subareas include the four internally-defined 
      !CLASSW variables XSNOCS, XSNOGS, XSNOWC and XSNOWG, which are set 
      !in subroutine CHKWAT to 1 if the subarea snow depth is greater 
      !than zero, and to zero otherwise. If the snow depth over the CS 
      !and GS subareas is greater than zero (meaning that there was a 
      !pre-existing snow cover at the beginning of the time step), the 
      !average snow albedo ALBSNO is preferentially set to the average 
      !over these two subareas. Otherwise ALBSNO is set to the average 
      !over the C and G subareas (where snow has just been added in the 
      !current time step). The snow temperature TSNOW and density RHOSNO 
      !are set to weighted averages over the four subareas, using the 
      !internally-defined subarea volumetric heat capacities 
      !HCPSCS/GS/C/G and RHOSCS/GS/C/G. Finally the snow depth ZSNOW is 
      !calculated from the subarea depths; the liquid water content of 
      !the snow pack WSNOW is obtained as a weighted average over the CS 
      !and GS subareas (assuming that freshly fallen snow does not yet 
      !contain liquid water); and the snow mass is determined from ZSNOW 
      !and RHOSNO. As in the case of intercepted and ponded water, if 
      !the snow mass is vanishingly small it and its liquid water 
      !content are added to the overland flow and to the total runoff, 
      !and their respective temperatures are recalculated. The 
      !diagnostic arrays ROFN, PCPG and HTCS are updated, and TSNOW, 
      !RHOSNO, SNO and WSNOW are set to zero. Flags are set to trigger 
      !calls to abort if TSNOW is less than 0 K or greater than 0.001 C. 
      !Finally, the three abort flags set thus far are checked, and 
      !calls to abort are performed if they are greater than zero.
      !
      DO 650 I=IL1,IL2     
          IF(ZSNOCS(I).GT.0. .OR. ZSNOGS(I).GT.0. .OR.
     1       ZSNOWC(I).GT.0. .OR. ZSNOWG(I).GT.0.)              THEN
              IF(ZSNOCS(I).GT.0. .OR. ZSNOGS(I).GT.0.)    THEN
                  ALBSNO(I)=(FCS(I)*ALBSCS(I)*XSNOCS(I)+
     1                       FGS(I)*ALBSGS(I)*XSNOGS(I))/
     2                      (FCS(I)*XSNOCS(I)+FGS(I)*XSNOGS(I))
              ELSE
                  ALBSNO(I)=(FC (I)*ALBSC(I)*XSNOWC(I) +
     1                       FG (I)*ALBSG(I)*XSNOWG(I))/
     2                      (FC (I)*XSNOWC(I)+FG (I)*XSNOWG(I))
              ENDIF
              TSNOW(I)=(FCS(I)*(TSNOCS(I)+TFREZ)*HCPSCS(I)*
     1                  ZSNOCS(I)*XSNOCS(I) +
     2                  FGS(I)*(TSNOGS(I)+TFREZ)*HCPSGS(I)*
     3                  ZSNOGS(I)*XSNOGS(I) +
     4                  FC (I)*(TSNOWC(I)+TFREZ)*HCPSC(I)*
     5                  ZSNOWC(I)*XSNOWC(I) +
     6                  FG (I)*(TSNOWG(I)+TFREZ)*HCPSG(I)*
     7                  ZSNOWG(I)*XSNOWG(I))/
     8                 (FCS(I)*HCPSCS(I)*ZSNOCS(I)*XSNOCS(I) +
     9                  FGS(I)*HCPSGS(I)*ZSNOGS(I)*XSNOGS(I) +
     A                  FC (I)*HCPSC(I)*ZSNOWC(I)*XSNOWC(I) +
     B                  FG (I)*HCPSG(I)*ZSNOWG(I)*XSNOWG(I))
              RHOSNO(I)=(FCS(I)*RHOSCS(I)*ZSNOCS(I)*XSNOCS(I) +
     1                   FGS(I)*RHOSGS(I)*ZSNOGS(I)*XSNOGS(I) +
     2                   FC (I)*RHOSC(I)*ZSNOWC(I)*XSNOWC(I) +
     3                   FG (I)*RHOSG(I)*ZSNOWG(I)*XSNOWG(I))/
     4                  (FCS(I)*ZSNOCS(I)*XSNOCS(I) +
     5                   FGS(I)*ZSNOGS(I)*XSNOGS(I) +
     6                   FC (I)*ZSNOWC(I)*XSNOWC(I) +
     7                   FG (I)*ZSNOWG(I)*XSNOWG(I))
              ZSNOW(I)=FCS(I)*ZSNOCS(I) + FGS(I)*ZSNOGS(I) +
     1                 FC (I)*ZSNOWC(I) + FG (I)*ZSNOWG(I)
              WSNOW(I)=FCS(I)*WSNOCS(I) + FGS(I)*WSNOGS(I)
              SNO(I)=ZSNOW(I)*RHOSNO(I)
              IF(SNO(I).LT.0.0) SNO(I)=0.0
              IF(SNO(I).LT.1.0E-2 .AND. SNO(I).GT.0.0) THEN
                  TOVRFL(I)=(TOVRFL(I)*OVRFLW(I)+TSNOW(I)*(SNO(I)+
     1                WSNOW(I))/DELT)/(OVRFLW(I)+(SNO(I)+WSNOW(I))/
     2                DELT)
                  OVRFLW(I)=OVRFLW(I)+(SNO(I)+WSNOW(I))/DELT
                  TRUNOF(I)=(TRUNOF(I)*RUNOFF(I)+TSNOW(I)*(SNO(I)+
     1                WSNOW(I))/DELT)/(RUNOFF(I)+(SNO(I)+WSNOW(I))/
     2                DELT)
                  RUNOFF(I)=RUNOFF(I)+(SNO(I)+WSNOW(I))/DELT
                  ROFN(I)=ROFN(I)+(SNO(I)+WSNOW(I))/DELT
                  PCPG(I)=PCPG(I)+(SNO(I)+WSNOW(I))/DELT
                  HTCS(I)=HTCS(I)-TSNOW(I)*(SPHICE*SNO(I)+SPHW*
     1                WSNOW(I))/DELT
                  TSNOW(I)=0.0
                  RHOSNO(I)=0.0
                  SNO(I)=0.0
                  WSNOW(I)=0.0
              ENDIF
          ELSE
              TSNOW(I)=0.0
              RHOSNO(I)=0.0
              SNO(I)=0.0
              WSNOW(I)=0.0
          ENDIF
C
          IF(TSNOW(I).LT.0.0) KPTBAD=I
          IF((TSNOW(I)-TFREZ).GT.1.0E-3) LPTBAD=I
  650 CONTINUE
C
      IF(JPTBAD.NE.0)                                               THEN
          WRITE(6,6625) JPTBAD,JL,TCAN(JPTBAD)
 6625     FORMAT('0AT (I,J)= (',I3,',',I3,'), TCAN = ',F10.5)
          CALL XIT('CLASSW2',-2)
      ENDIF
C
      IF(KPTBAD.NE.0)                                               THEN
          WRITE(6,6626) KPTBAD,JL,TSNOW(KPTBAD)
 6626     FORMAT('0AT (I,J)= (',I3,',',I3,'), TSNOW = ',F10.5)
          CALL XIT('CLASSW2',-3)
      ENDIF
C
      IF(LPTBAD.NE.0)                                               THEN
          WRITE(6,6626) LPTBAD,JL,TSNOW(LPTBAD)
          CALL XIT('CLASSW2',-4)
      ENDIF
C
      !
      !In the 700 loop, the temperature of each soil layer is calculated 
      !as a weighted average over the four subareas. In the case of the 
      !third soil layer., if the standard three-layer configuration is 
      !being modelled (with a very thick third soil layer of 3.75 m), 
      !the subarea layer temperatures TBARCS/GS/C/G and the layer heat 
      !capacities HCPCS/GS/C/G apply to the permeable depth DELZW of the 
      !layer, and the bedrock temperature TBASE and the rock heat 
      !capacity HCPSND to the remainder, DELZ-DELZW. The averaging is 
      !carried out accordingly. In all other soil layers, the layer 
      !temperature applies to the whole thickness, whose heat capacity 
      !is a weighted average of HCPCS/GS/C/G over DELZW and HCPSND over 
      !DELZ-DELZW. The volumetric liquid water content THLIQ, the 
      !volumetric frozen water content THICE, and the heat flux at the 
      !soil layer interfaces GFLUX are calculated as simple weighted 
      !averages over the subareas. A flag is set to trigger a call to 
      !abort if the soil layer temperature is less than -100 C or 
      !greater than 100 C, and after the end of the loop, a call to 
      !abort is performed if the flag is greater than zero.
      !
      IPTBAD=0
      DO 700 J=1,IG
      DO 700 I=IL1,IL2
          IF(IG.EQ. 3 .AND. J.EQ.IG .AND. ISAND(I,1).GT.-4)    THEN
              TBAR(I,J)=((FCS(I)*(TBARCS(I,J)+TFREZ)*HCPCS(I,J) +
     1                   FGS(I)*(TBARGS(I,J)+TFREZ)*HCPGS(I,J) +
     2                   FC (I)*(TBARC (I,J)+TFREZ)*HCPCO(I,J) +
     3                   FG (I)*(TBARG (I,J)+TFREZ)*HCPGO(I,J))*
     4                   DELZW(I,J)+TBASE(I)*HCPSND*
     5                   (DELZ(J)-DELZW(I,J)))/
     4                  ((FCS(I)*HCPCS(I,J) + FGS(I)*HCPGS(I,J) +
     5                   FC (I)*HCPCO(I,J) + FG (I)*HCPGO(I,J))*
     8                   DELZW(I,J)+HCPSND*(DELZ(J)-DELZW(I,J)))
          ELSE
              TBAR(I,J)=(FCS(I)*(TBARCS(I,J)+TFREZ)*(DELZW(I,J)*
     1                   HCPCS(I,J)+(DELZ(J)-DELZW(I,J))*HCPSND)+
     2                   FGS(I)*(TBARGS(I,J)+TFREZ)*(DELZW(I,J)*
     3                   HCPGS(I,J)+(DELZ(J)-DELZW(I,J))*HCPSND)+
     4                   FC (I)*(TBARC (I,J)+TFREZ)*(DELZW(I,J)*
     5                   HCPCO(I,J)+(DELZ(J)-DELZW(I,J))*HCPSND)+
     6                   FG (I)*(TBARG (I,J)+TFREZ)*(DELZW(I,J)*
     7                   HCPGO(I,J)+(DELZ(J)-DELZW(I,J))*HCPSND))/
     8                  (FCS(I)*(DELZW(I,J)*HCPCS(I,J)+
     9                   (DELZ(J)-DELZW(I,J))*HCPSND) +
     A                   FGS(I)*(DELZW(I,J)*HCPGS(I,J)+
     B                   (DELZ(J)-DELZW(I,J))*HCPSND) +
     C                   FC (I)*(DELZW(I,J)*HCPCO(I,J)+
     D                   (DELZ(J)-DELZW(I,J))*HCPSND) +
     E                   FG (I)*(DELZW(I,J)*HCPGO(I,J)+
     F                   (DELZ(J)-DELZW(I,J))*HCPSND))
          ENDIF
          THLIQ(I,J)=FCS(I)*THLQCS(I,J)+FGS(I)*THLQGS(I,J)+
     1               FC (I)*THLQCO(I,J)+FG (I)*THLQGO(I,J)
          THICE(I,J)=FCS(I)*THICCS(I,J)+FGS(I)*THICGS(I,J)+
     1               FC (I)*THICCO(I,J)+FG (I)*THICGO(I,J)
          GFLUX(I,J)=FCS(I)*GFLXCS(I,J)+FGS(I)*GFLXGS(I,J)+
     1               FC (I)*GFLXC (I,J)+FG (I)*GFLXG (I,J)
C     ipy test
C          IF(THLIQ(I,J).GT.THFC(I,J))                               THEN
C              BASFLW(I)=BASFLW(I)+(THLIQ(I,J)-THFC(I,J))*DELZW(I,J)*
C     1            RHOW/DELT
C              RUNOFF(I)=RUNOFF(I)+(THLIQ(I,J)-THFC(I,J))*DELZW(I,J)*
C     1            RHOW/DELT
C              HTC(I,J)=HTC(I,J)-TBAR(I,J)*(THLIQ(I,J)-THFC(I,J))*
C     1            HCPW*DELZW(I,J)/DELT
C              THLIQ(I,J)=THFC(I,J)
C          ENDIF
          IF(TBAR(I,1).LT.173.16 .OR. TBAR(I,1).GT.373.16) IPTBAD=I
  700 CONTINUE
C
      IF(IPTBAD.NE.0)                                               THEN
          WRITE(6,6600) IPTBAD,JL,TBAR(IPTBAD,1)
 6600     FORMAT('0AT (I,J)= (',I3,',',I3,'), TBAR(1) = ',F10.5)
          CALL XIT('CLASSW2',-1)
      ENDIF
C
      !Finally, subroutine CGROW is called to update the vegetation 
      !growth index.
      !
      CALL CGROW(GROWTH,TBAR,TA,FC,FCS,ILG,IG,IL1,IL2,JL)
C
      RETURN
      END
