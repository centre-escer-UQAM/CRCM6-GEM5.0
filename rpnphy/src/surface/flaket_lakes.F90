! File %M% from Library %Q%
! Version %I% from %G% extracted: %H%
!------------------------------------------------------------------------------

SUBROUTINE flaket_lakes ( dMsnowdt_inR, I_atm_inR, Q_atm_lw_inR, height_u_inR, height_tq_inR,     &
               U_a_inR, T_a_inR, q_a_inR, P_a_inR,etaR,                                     &
               depth_wR, fetchR, depth_bsR, T_bsR, par_CoriolisR, del_timeR,           &
               T_snow_inR,  T_ice_inR,  T_mnw_inR,  T_wML_inR,  T_bot_inR,  T_B1_inR,  &
               C_T_inR,  h_snow_inR,  h_ice_inR,  h_ML_inR,  H_B1_inR, T_sfc_pR,aR, albedoR,        &
               T_snow_outR, T_ice_outR, T_mnw_outR, T_wML_outR, T_bot_outR, T_B1_outR,  & 
               h_snow_outR, h_ice_outR, h_ML_outR, H_B1_outR,Q_sensibleR,    &
               pR,Q_latentR, T_sfc_nR, C_T_outR,trnch,n,m,i,frac )

!------------------------------------------------------------------------------
!
! Description:
!
!  The FLake interface is
!  a communication routine between "flake_driver"
!  and a prediction system that uses FLake.
!  It assigns the FLake variables at the previous time step 
!  to their input values given by the driving model,
!  calls a number of routines to compute the heat and radiation fluxes,
!  calls "flake_driver",
!  and returns the updated FLake variables to the driving model.
!  The "flake_interface" does not contain any Flake physics. 
!  It only serves as a convenient means to organize calls of "flake_driver"
!  and of external routines that compute heat and radiation fluxes.
!  The interface may (should) be changed so that to provide 
!  the most convenient use of FLake.
!  Within a 3D atmospheric prediction system,
!  "flake_driver" may be called in a DO loop within "flake_interface" 
!  for each grid-point where a lake is present.
!  In this way, the driving atmospheric model should call "flake_interface"
!  only once, passing the FLake variables to "flake_interface" as 2D fields. 
!
!  Lines embraced with "!_tmp" contain temporary parts of the code.
!  These should be removed prior to using FLake in applications.
!  Lines embraced/marked with "!_dev" may be replaced
!  as improved parameterizations are developed and tested.
!  Lines embraced/marked with "!_dm" are DM's comments
!  that may be helpful to a user.
!
!
! Current Code Owner: DWD, Dmitrii Mironov
!  Phone:  +49-69-8062 2705
!  Fax:    +49-69-8062 3721
!  E-mail: dmitrii.mironov@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! N.NN       YYYY/MM/DD Dmitrii Mironov
!  Initial release
! !VERSION!  !DATE!     <Your name>
!  <Modification comments>
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================
!
! Declarations:
!
! Modules used:

USE flake_data_parameters , ONLY : &
    ireals,                  & ! KIND-type parameter for real variables
    iintegers                  ! KIND-type parameter for "normal" integer variables

USE flake_derivedtypes         ! Definitions of several derived TYPEs

USE flake_parameters , ONLY :   &
  tpl_T_f                     , & ! Fresh water freezing point [K]
  tpl_rho_w_r                 , & ! Maximum density of fresh water [kg m^{-3}]
  h_Snow_min_flk              , & ! Minimum snow thickness [m]
  h_Ice_min_flk                   ! Minimum ice thickness [m]

USE flake_paramoptic_ref       ! Reference values of the optical characteristics
                               ! of the lake water, lake ice and snow 

USE flakem_albedo_ref           ! Reference values the albedo for the lake water, lake ice and snow

USE flaket           , ONLY :    &
  flake_driver                , & ! Subroutine, FLake driver
  flake_radflux               , & ! Subroutine, computes radiation fluxes at various depths
                                  !
  T_snow_p_flk, T_snow_n_flk  , & ! Temperature at the air-snow interface [K]
  T_ice_p_flk, T_ice_n_flk    , & ! Temperature at the snow-ice or air-ice interface [K]
  T_mnw_p_flk, T_mnw_n_flk    , & ! Mean temperature of the water column [K]
  T_wML_p_flk, T_wML_n_flk    , & ! Mixed-layer temperature [K]
  T_bot_p_flk, T_bot_n_flk    , & ! Temperature at the water-bottom sediment interface [K]
  T_B1_p_flk, T_B1_n_flk      , & ! Temperature at the bottom of the upper layer of the sediments [K]
  C_T_p_flk, C_T_n_flk        , & ! Shape factor (thermocline)
  h_snow_p_flk, h_snow_n_flk  , & ! Snow thickness [m]
  h_ice_p_flk, h_ice_n_flk    , & ! Ice thickness [m]
  h_ML_p_flk, h_ML_n_flk      , & ! Thickness of the mixed-layer [m]
  H_B1_p_flk, H_B1_n_flk      , & ! Thickness of the upper layer of bottom sediments [m]
                                  !
  Q_snow_flk                  , & ! Heat flux through the air-snow interface [W m^{-2}]
  Q_ice_flk                   , & ! Heat flux through the snow-ice or air-ice interface [W m^{-2}]
  Q_w_flk                     , & ! Heat flux through the ice-water or air-water interface [W m^{-2}]
  Q_bot_flk                   , & ! Heat flux through the water-bottom sediment interface [W m^{-2}]
  I_atm_flk                   , & ! Radiation flux at the lower boundary of the atmosphere [W m^{-2}],
                                  ! i.e. the incident radiation flux with no regard for the surface albedo
  I_snow_flk                  , & ! Radiation flux through the air-snow interface [W m^{-2}]
  I_ice_flk                   , & ! Radiation flux through the snow-ice or air-ice interface [W m^{-2}]
  I_w_flk                     , & ! Radiation flux through the ice-water or air-water interface [W m^{-2}]
  I_h_flk                     , & ! Radiation flux through the mixed-layer-thermocline interface [W m^{-2}]
  I_bot_flk                   , & ! Radiation flux through the water-bottom sediment interface [W m^{-2}]
  I_intm_0_h_flk              , & ! Mean radiation flux over the mixed layer [W m^{-1}]
  I_intm_h_D_flk              , & ! Mean radiation flux over the thermocline [W m^{-1}]
  Q_star_flk                  , & ! A generalized heat flux scale [W m^{-2}]
  u_star_w_flk                , & ! Friction velocity in the surface layer of lake water [m s^{-1}]
  w_star_sfc_flk              , & ! Convective velocity scale, using a generalized heat flux scale [m s^{-1}]
  dMsnowdt_flk                    ! The rate of snow accumulation [kg m^{-2} s^{-1}]


USE flakem_SfcFlx          , ONLY :    &
  flake_SfcFlx_lwradwsfc            , & ! Function, returns the surface long-wave radiation flux
  flake_SfcFlx_momsenlat                ! Subroutine, computes fluxes of momentum and of sensible and latent heat

!==============================================================================

IMPLICIT NONE

!==============================================================================
!
! Declarations

!  Input (procedure arguments)

REAL, INTENT(INOUT) ::   &
  dMsnowdt_inR                       , & ! The rate of snow accumulation [kg m^{-2} s^{-1}]
  I_atm_inR                          , & ! Solar radiation flux at the surface [W m^{-2}]
  Q_atm_lw_inR                       , & ! Long-wave radiation flux from the atmosphere [W m^{-2}]
  height_u_inR                       , & ! Height above the lake surface where the wind speed is measured [m]
  height_tq_inR                      , & ! Height where temperature and humidity are measured [m]
  U_a_inR                            , & ! Wind speed at z=height_u_in [m s^{-1}]
  T_a_inR                            , & ! Air temperature at z=height_tq_in [K]
  q_a_inR                            , & ! Air specific humidity at z=height_tq_in
  P_a_inR                            , & ! Surface air pressure [N m^{-2} = kg m^{-1} s^{-2}]
  etaR                                   ! The SW extinction coefficient in water [m-1]
REAL, INTENT(INOUT) ::   &
  depth_wR                           , & ! The lake depth [m]
  fetchR                             , & ! Typical wind fetch [m]
  depth_bsR                          , & ! Depth of the thermally active layer of the bottom sediments [m]
  T_bsR                              , & ! Temperature at the outer edge of 
                                        ! the thermally active layer of the bottom sediments [K]
  par_CoriolisR                      , & ! The Coriolis parameter [s^{-1}]
  del_timeR                              ! The model time step [s]

REAL, INTENT(INOUT)  :: &
  T_snow_inR                        , & ! Temperature at the air-snow interface [K] 
  T_ice_inR                         , & ! Temperature at the snow-ice or air-ice interface [K]
  T_mnw_inR                         , & ! Mean temperature of the water column [K]
  T_wML_inR                         , & ! Mixed-layer temperature [K]
  T_bot_inR                         , & ! Temperature at the water-bottom sediment interface [K]
  T_B1_inR                          , & ! Temperature at the bottom of the upper layer of the sediments [K]
  C_T_inR                           , & ! Shape factor (thermocline)
  h_snow_inR                        , & ! Snow thickness [m]
  h_ice_inR                         , & ! Ice thickness [m]
  h_ML_inR                          , & ! Thickness of the mixed-layer [m]
  H_B1_inR                          , & ! Thickness of the upper layer of bottom sediments [m]
  T_sfc_pR                              ! Surface temperature at the previous time step [K]  

!  Input/Output (procedure arguments)

REAL (KIND = ireals)                 :: &
  albedo_water                        , & ! Water surface albedo with respect to the solar radiation
  albedo_ice                          , & ! Ice surface albedo with respect to the solar radiation
  albedo_snow                             ! Snow surface albedo with respect to the solar radiation

TYPE (opticpar_medium)               :: & 
  opticpar_water                       , & ! Optical characteristics of water
  opticpar_ice                         , & ! Optical characteristics of ice
  opticpar_snow                            ! Optical characteristics of snow 

!  Output (procedure arguments)

REAL,  INTENT(INOUT)  :: &
  T_snow_outR                        , & ! Temperature at the air-snow interface [K] 
  T_ice_outR                         , & ! Temperature at the snow-ice or air-ice interface [K]
  T_mnw_outR                         , & ! Mean temperature of the water column [K]
  T_wML_outR                         , & ! Mixed-layer temperature [K]
  T_bot_outR                         , & ! Temperature at the water-bottom sediment interface [K]
  T_B1_outR                          , & ! Temperature at the bottom of the upper layer of the sediments [K]
  C_T_outR                           , & ! Shape factor (thermocline)
  h_snow_outR                        , & ! Snow thickness [m]
  h_ice_outR                         , & ! Ice thickness [m]
  h_ML_outR                          , & ! Thickness of the mixed-layer [m]
  H_B1_outR                          , & ! Thickness of the upper layer of bottom sediments [m]
  T_sfc_nR                           , & ! Updated surface temperature [K]  
  albedoR                            ,&    ! Effective albedo
  Q_sensibleR             , & ! Sensible heat flux [W m^{-2}]
  Q_latentR               ,& ! Latent heat flux [W m^{-2}]
  pR,aR ,frac                        ! empty parameter


INTEGER,  INTENT(IN)  :: trnch,n,m,i

!  Local variables of type REAL


REAL (KIND = ireals) ::    &
  Q_momentum             , & ! Momentum flux [N m^{-2}]
  Q_sensible             , & ! Sensible heat flux [W m^{-2}]
  Q_latent               , & ! Latent heat flux [W m^{-2}]
  Q_watvap               , &    ! Flux of water vapour [kg m^{-2} s^{-1}]
  albedo


REAL (KIND = ireals) ::   &
  dMsnowdt_in                       , & ! The rate of snow accumulation [kg m^{-2} s^{-1}]
  I_atm_in                          , & ! Solar radiation flux at the surface [W m^{-2}]
  Q_atm_lw_in                       , & ! Long-wave radiation flux from the atmosphere [W m^{-2}]
  height_u_in                       , & ! Height above the lake surface where the wind speed is measured [m]
  height_tq_in                      , & ! Height where temperature and humidity are measured [m]
  U_a_in                            , & ! Wind speed at z=height_u_in [m s^{-1}]
  T_a_in                            , & ! Air temperature at z=height_tq_in [K]
  q_a_in                            , & ! Air specific humidity at z=height_tq_in
  P_a_in                            , & ! Surface air pressure [N m^{-2} = kg m^{-1} s^{-2}]
  eta                                   ! The SW extinction coefficient in water [m-1]

REAL (KIND = ireals) ::   &
  depth_w                           , & ! The lake depth [m]
  fetch                             , & ! Typical wind fetch [m]
  depth_bs                          , & ! Depth of the thermally active layer of the bottom sediments [m]
  T_bs                              , & ! Temperature at the outer edge of 
                                        ! the thermally active layer of the bottom sediments [K]
  par_Coriolis                      , & ! The Coriolis parameter [s^{-1}]
  del_time                              ! The model time step [s]

REAL (KIND = ireals) :: &
  T_snow_in                        , & ! Temperature at the air-snow interface [K] 
  T_ice_in                         , & ! Temperature at the snow-ice or air-ice interface [K]
  T_mnw_in                         , & ! Mean temperature of the water column [K]
  T_wML_in                         , & ! Mixed-layer temperature [K]
  T_bot_in                         , & ! Temperature at the water-bottom sediment interface [K]
  T_B1_in                          , & ! Temperature at the bottom of the upper layer of the sediments [K]
  C_T_in                           , & ! Shape factor (thermocline)
  h_snow_in                        , & ! Snow thickness [m]
  h_ice_in                         , & ! Ice thickness [m]
  h_ML_in                          , & ! Thickness of the mixed-layer [m]
  H_B1_in                          , & ! Thickness of the upper layer of bottom sediments [m]
  T_sfc_p                              ! Surface temperature at the previous time step [K]  


!  Output (procedure arguments)

REAL (KIND = ireals)  :: &
  T_snow_out                        , & ! Temperature at the air-snow interface [K] 
  T_ice_out                         , & ! Temperature at the snow-ice or air-ice interface [K]
  T_mnw_out                         , & ! Mean temperature of the water column [K]
  T_wML_out                         , & ! Mixed-layer temperature [K]
  T_bot_out                         , & ! Temperature at the water-bottom sediment interface [K]
  T_B1_out                          , & ! Temperature at the bottom of the upper layer of the sediments [K]
  C_T_out                           , & ! Shape factor (thermocline)
  h_snow_out                        , & ! Snow thickness [m]
  h_ice_out                         , & ! Ice thickness [m]
  h_ML_out                          , & ! Thickness of the mixed-layer [m]
  H_B1_out                          , & ! Thickness of the upper layer of bottom sediments [m]
  T_sfc_n                               ! Updated surface temperature [K]  




 dMsnowdt_in =  (dMsnowdt_inR)
 I_atm_in =  (I_atm_inR)
 Q_atm_lw_in =  (Q_atm_lw_inR)
 height_u_in =  (height_u_inR)
 height_tq_in =  (height_tq_inR)
 U_a_in =  (U_a_inR)
 T_a_in =  (T_a_inR)
 q_a_in =  (q_a_inR)
 P_a_in =  (P_a_inR)
 eta    =  (etaR)
 depth_w =  (depth_wR)
 fetch =  (fetchR)
 depth_bs =  (depth_bsR)
 T_bs =  (T_bsR)
 par_Coriolis =  (par_CoriolisR)
 del_time =  (del_timeR)
 T_snow_in =  (T_snow_inR)
 T_ice_in =  (T_ice_inR)
 T_mnw_in =  (T_mnw_inR)
 T_wML_in =  (T_wML_inR)
 T_bot_in =  (T_bot_inR)
 T_B1_in =  (T_B1_inR)
 C_T_in =  (C_T_inR)
 h_snow_in =  (h_snow_inR)
 h_ice_in =  (h_ice_inR)
 h_ML_in =  (h_ML_inR)
 H_B1_in =  (H_B1_inR)
 T_sfc_p =  (T_sfc_pR)
 !albedo =  (albedoR   ) 


!==============================================================================
!  Start calculations
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!  Set albedos of the lake water, lake ice and snow
!  Provided by CLASS
!------------------------------------------------------------------------------

! Use default value 
albedo_water = albedo_water_ref !AM albedo from CLASS
! Use empirical formulation proposed by Mironov and Ritter (2004) for GME 
!_nu albedo_ice   = albedo_whiteice_ref 
albedo_ice   = EXP(-c_albice_MR*(tpl_T_f-T_sfc_p)/tpl_T_f)
albedo_ice   = albedo_whiteice_ref*(1._ireals-albedo_ice) + albedo_blueice_ref*albedo_ice
! Snow is not considered
albedo_snow  = albedo_ice  

	if(trnch.eq.-1.and.n.eq.27.and.m.eq.27.and.i.eq.18) then
write(*,4507)albedo_water,albedo_ice,c_albice_MR,tpl_T_f,T_sfc_p,albedo_whiteice_ref, &
albedo_blueice_ref
endif

!------------------------------------------------------------------------------
!  Set optical characteristics of the lake water, lake ice and snow
!------------------------------------------------------------------------------

! Use default values
opticpar_water = opticpar_water_ref
opticpar_ice   = opticpar_ice_opaque   ! Opaque ice
opticpar_snow  = opticpar_snow_opaque  ! Opaque snow

! Externally-specified water transparency
opticpar_water%extincoef_optic(1) = eta


!------------------------------------------------------------------------------
!  Set initial values
!------------------------------------------------------------------------------

T_snow_p_flk = T_snow_in     
T_ice_p_flk  = T_ice_in         
T_mnw_p_flk  = T_mnw_in        
T_wML_p_flk  = T_wML_in       
T_bot_p_flk  = T_bot_in      
T_B1_p_flk   = T_B1_in       
C_T_p_flk    = C_T_in        
h_snow_p_flk = h_snow_in      
h_ice_p_flk  = h_ice_in       
h_ML_p_flk   = h_ML_in        
H_B1_p_flk   = H_B1_in       

!------------------------------------------------------------------------------
!  Set the rate of snow accumulation
!------------------------------------------------------------------------------

dMsnowdt_flk = dMsnowdt_in  

!------------------------------------------------------------------------------
!  Compute solar radiation fluxes (positive downward)
!------------------------------------------------------------------------------

I_atm_flk = I_atm_in
CALL flake_radflux ( depth_w, albedo_water, albedo_ice, albedo_snow, &
                     opticpar_water, opticpar_ice, opticpar_snow )

!------------------------------------------------------------------------------
!  Compute long-wave radiation fluxes (positive downward)
!------------------------------------------------------------------------------

Q_w_flk = Q_atm_lw_in                          ! Radiation of the atmosphere 
Q_w_flk = Q_w_flk - flake_SfcFlx_lwradwsfc(T_sfc_p)  ! Radiation of the surface (notice the sign)

!------------------------------------------------------------------------------
!  Compute the surface friction velocity and fluxes of sensible and latent heat 
!------------------------------------------------------------------------------


CALL flake_SfcFlx_momsenlat ( height_u_in, height_tq_in, fetch,                      &
                        U_a_in, T_a_in, q_a_in, T_sfc_p, P_a_in, h_ice_p_flk,  &
                        Q_momentum, Q_sensible, Q_latent, Q_watvap, trnch,i,frac )
u_star_w_flk = SQRT(-Q_momentum/tpl_rho_w_r)

!------------------------------------------------------------------------------
!  Compute heat fluxes Q_snow_flk, Q_ice_flk, Q_w_flk
!------------------------------------------------------------------------------

Q_w_flk = Q_w_flk - Q_sensible - Q_latent  ! Add sensible and latent heat fluxes (notice the signs)


!if(trnch.eq.1.and.n.eq.34.and.m.eq.34.and.i.eq.7) then
!write(6,3507)frac,I_atm_in,Q_atm_lw_in ,flake_SfcFlx_lwradwsfc(T_sfc_p),  &
!     -Q_sensible,-Q_latent,Q_w_flk
!endif


IF(h_ice_p_flk.GE.h_Ice_min_flk) THEN            ! Ice exists
  IF(h_snow_p_flk.GE.h_Snow_min_flk) THEN        ! There is snow above the ice
    Q_snow_flk = Q_w_flk
    Q_ice_flk  = 0._ireals
    Q_w_flk    = 0._ireals
  ELSE                                           ! No snow above the ice
    Q_snow_flk = 0._ireals
    Q_ice_flk  = Q_w_flk
    Q_w_flk    = 0._ireals
  END IF
ELSE                                             ! No ice-snow cover
    Q_snow_flk = 0._ireals
    Q_ice_flk  = 0._ireals
END IF

!------------------------------------------------------------------------------
!  Advance FLake variables
!------------------------------------------------------------------------------

CALL flake_driver ( depth_w, depth_bs, T_bs, par_Coriolis,         &
                    opticpar_water%extincoef_optic(1),             &
                    del_time, T_sfc_p, T_sfc_n )

!------------------------------------------------------------------------------
!  Set output values
!------------------------------------------------------------------------------

T_snow_out = T_snow_n_flk  
T_ice_out  = T_ice_n_flk      
T_mnw_out  = T_mnw_n_flk     
T_wML_out  = T_wML_n_flk    
T_bot_out  = T_bot_n_flk   
T_B1_out   = T_B1_n_flk    
C_T_out    = C_T_n_flk     
h_snow_out = h_snow_n_flk   
h_ice_out  = h_ice_n_flk    
h_ML_out   = h_ML_n_flk     
H_B1_out   = H_B1_n_flk    


 T_snow_outR = (T_snow_out)
 T_ice_outR = (T_ice_out)
 T_mnw_outR = (T_mnw_out)
 T_wML_outR =  (T_wML_out)
 T_bot_outR =  (T_bot_out)
 T_B1_outR =  (T_B1_out)
 C_T_outR =  (C_T_out)
! write(6,*)"CTOUT1 ",C_T_out
! write(6,*)"CTOUT2 ",C_T_outR
 h_snow_outR =  (h_snow_out)
 h_ice_outR =  (h_ice_out)
 h_ML_outR =  (h_ML_out)
 H_B1_outR =  (H_B1_out)
 T_sfc_nR  =  (T_sfc_n)
 Q_sensibleR =  (Q_sensible)
 Q_latentR =  (Q_latent)

 IF(h_ice_out.GT.0.0) THEN
 albedo=albedo_ice
 ELSE
 albedo=albedo_water
 ENDIF
 albedoR= (albedo)


2604     FORMAT('IFLAKEAFTER 11 ',8(e14.8,' '))
3604     FORMAT('IFLAKEAFTER 12 ',8(e14.8,' '))

2605     FORMAT('IFLAKEAFTER 21 ',8(e14.8,' '))
3507     FORMAT('FLAKEFLUX ',8(e14.8,' '))
3605     FORMAT('IFLAKEAFTER 22 ',8(e14.8,' '))
3606     FORMAT('IFLAKEAFTER 23 ',8(e14.8,' '))
4507     FORMAT('AFLAKEAFTER 01 ',8(e14.8,' '))

 


!------------------------------------------------------------------------------
!  End calculations
!==============================================================================

END SUBROUTINE flaket_lakes

