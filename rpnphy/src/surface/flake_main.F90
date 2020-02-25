!-------------------------------------- LICENCE BEGIN ------------------------------------
!Environment Canada - Atmospheric Science and Technology License/Disclaimer,
!                     version 3; Last Modified: May 7, 2008.
!This is free but copyrighted software; you can use/redistribute/modify it under the terms
!of the Environment Canada - Atmospheric Science and Technology License/Disclaimer
!version 3 or (at your option) any later version that should be found at:
!http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html
!
!This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
!without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!See the above mentioned License/Disclaimer for more details.
!You should have received a copy of the License/Disclaimer along with this software;
!if not, you can write to: EC-RPN COMM Group, 2121 TransCanada, suite 500, Dorval (Quebec),
!CANADA, H9P 1J3; or send e-mail to service.rpn@ec.gc.ca
!-------------------------------------- LICENCE END --------------------------------------
   !**
   !
SUBROUTINE flake_main ( bus, bussiz,          &
                        ptsurf, ptsurfsiz,dt, &
                        trnch, kount,         &
                        n, m, nk )
   use sfclayer_mod, only: sl_prelim,sl_sfclayer,SL_OK
   use sfc_options
   use sfcbus_mod
   implicit none
#include <arch_specific.hf>
   !
   !
   integer bussiz, kount, trnch
   real,target :: bus(bussiz)
   integer ptsurfsiz
   integer ptsurf(ptsurfsiz)
   integer n, m, nk
   real dt
   !
   !Author
   !          A. Martynov (UQAM)
   !
   !Revisions
   ! 001      A. Martynov (Summer 2011) - Initial version
   ! 002      K. Winger   (Spring 2017) - Adapt for physics 5.8.3
   !
   !Arguments
   !             - Input/Output -
   ! BUS         Bus for the lake surface scheme
   !
   !             - Input -
   ! BUSSIZ      dimension of bus
   ! PTSURF      surface pointers
   ! PTSURFSIZ   dimension of ptsurf
   ! TRNCH       row number
   ! KOUNT       timestep number
   ! DT          length of timestep
   ! N           horizontal dimension (row length)
   ! M           horiatzontal dimensions of fields
   !             (not used for the moment)
   ! NK          vertical dimension
   !
   !
   !Notes
   !          Z0 = BETA*USTAR**2/GRAV (in metres) with minimum value
   !          Z0MIN and a maximum value Z0MAX
   !
   !IMPLICITES
   !
   !#include "consphy.cdk"
   !
   !#include "clefcon.cdk"
   !#include "options.cdk"
   !#include "lakes.cdk"
   include "thermoconsts.inc"
   include "dintern.inc"
   include "fintern.inc"
   !
   integer SURFLEN
#define x(fptr,fj,fk) ptsurf(vd%fptr%i)+(fk-1)*surflen+fj-1

   integer, parameter :: INDX_SFC = INDX_LAKE
   !
   logical, parameter :: LAKE_TDIAGLIM=.false.
   !
   INTEGER I, J, K
   !
   !
   !***************************************************
   !     AUTOMATIC ARRAYS
   !***************************************************
   !
   real,dimension(n) :: EMIST, VMOD, VDIR
   real,dimension(n) :: SCR1, SCR2, SCR3, SCR4, SCR5
   !
   !***************************************************
   !
   REAL RHO
   !

   !AMTV I/O VARIABLES

   real,pointer,dimension(:) ::  ALVIS_AVG, FC_AVG
   real,pointer,dimension(:) ::  FV_AVG, HU, ZSNODP
   real,pointer,dimension(:) ::  CMU_AVG,CTU_AVG,Z0H_AVG,Z0M_AVG
   real,pointer,dimension(:) ::  HST_AVG, ILMO_AVG,QLWIN,FSOL
   real,pointer,dimension(:) ::  PS, QS, TH, TT, UU, VV
   real,pointer,dimension(:) ::  ZALFAQ, ZALFAT, ZDLAT, ZFCOR,ZDLON
   real,pointer,dimension(:) ::  ZFTEMP, ZFVAP, ZQDIAG, ZTDIAG
   real,pointer,dimension(:) ::  ZTSURF, ZTSRAD, ZUDIAG, ZVDIAG
   real,pointer,dimension(:) ::  ZFRV_AVG, ZZUSL, ZZTSL
   real,pointer,dimension(:) ::  ZRAINRATE,ZSNOWRATE, zrunofftot
   real,pointer,dimension(:) ::  ztt2m

   !AMTV OPEN WATER VARIABLES

   real,pointer,dimension(:) :: ZFRV_WAT
   real, dimension(n) ::  ALVIS_WAT, CMU_WAT, CTU_WAT, FC_WAT
   real, dimension(n) ::  FV_WAT
   real, dimension(n) ::  HST_WAT, ILMO_WAT
   real, dimension(n) ::  Z0H_WAT, Z0M_WAT
   real, dimension(n) ::  ZALFAQ_WAT, ZALFAT_WAT,  ZFCOR_WAT
   real, dimension(n) ::  ZFTEMP_WAT, ZFVAP_WAT , ZQDIAG_WAT, ZTDIAG_WAT
   real, dimension(n) ::  ZTSURF_WAT, ZTSRAD_WAT, ZUDIAG_WAT, ZVDIAG_WAT
   real, dimension(n) ::  WA_WAT
   real, dimension(n) ::  QS_WAT,TS_WAT
   real WINDCRIT, ZFRVMIN
!Huziy uncommented
   real, dimension(n) ::   ZRUNOFF_WAT

   !AMTV ICE-COVER VARIABLES

   real,pointer,dimension(:) :: ZFRV_ICE
   real, dimension(n) ::  ALVIS_ICE, CMU_ICE, CTU_ICE, FC_ICE
   real, dimension(n) ::  FV_ICE
   real, dimension(n) ::  HST_ICE, ILMO_ICE
   real, dimension(n) ::  Z0H_ICE, Z0M_ICE
   real, dimension(n) ::  ZALFAQ_ICE, ZALFAT_ICE,  ZFCOR_ICE
   real, dimension(n) ::  ZFTEMP_ICE, ZFVAP_ICE,  ZQDIAG_ICE, ZTDIAG_ICE
   real, dimension(n) ::  ZTSURF_ICE, ZTSRAD_ICE, ZUDIAG_ICE, ZVDIAG_ICE
   real, dimension(n) ::  WA_ICE
   real, dimension(n) ::  QS_ICE, TS_ICE, DQS_ICE
!Huziy uncommented
   real, dimension(n) ::   ZRUNOFF_ICE

   !AMTV LAKE MODEL VARIABLES

   real,pointer,dimension(:)   :: LDEPTH,LFICE,LHICE,LHML,LFR
   real,pointer,dimension(:)   :: LTBOT,LTICE,LTMNW,LTRANSP,LTWML,LCT
   real FDEPTH

   real WIND, T_sfc_in, T_sfc_out, SWD, LHMIN, IDAY1, DAY1
   real julien1, HSNOW, T_ice_out, yz, XPI
   real, dimension(n) :: ZQDIAG_TMP, WA_AVG, PCPR, coszrs

   real JULIAND,HZ1,DAYLNT,albedo
   EXTERNAL JULIAND

   !FLAKE VARIABLES
   real CORIO,omega_earth,t_ice_in,t_snow_in,t_mnw_in
   real t_wml_in,t_bot_in,c_t_in,h_ml_in,t_b1_in
   real h_b1_in,t_air_in,tpl_t_r,t_bot_out
   real t_snow_out,t_mnw_out,t_wml_out,t_b1_out
   real c_t_out,h_snow_out,tpl_t_f
   real h_ice_out,h_ml_out,h_b1_out


   real CON1,CON2,CON3,CON4,CON5,CON6
   real CON7,CON8,CON9,CON10,CON11,CON12

   real FI0,CONDFI,TFRZW,TMELI,TMELS

   real ALBOW,ALBDI,ALBMI,ALBDS,ALBMS,EMISI,EMISNO,EMISW

   real COEFCOND,COEFHCAP,COEFEXT

   real ROICE,ROSNOW(2),ROSWTR

   real Z0ICE,Z0W, HICE_ALB(N)

   real BASEHF

   real HCAPI,HCAPW,VHFICE,VHFBAS,VHFSNO

   real HSMIN,DMIX


   DATA   CON1     , CON2   , CON3  , CON4  &
         /2.845E-6 , 2.7E-4 , 233.0 , 0.2   /
   DATA   CON5     , CON6   , CON7  , CON8  &
         /92.88    , 7.364  , 3.2   , 14.24 /
   DATA   CON9     , CON10 , CON11 , CON12  &
         /19.39    , 0.1   , 0.44  , 0.075  /
   DATA  TFRZW , TMELI , TMELS /  271.2 , 273.05 , 273.15  /
   DATA   ALBOW   ,  ALBDI  ,  ALBMI ,  ALBDS  ,  ALBMS &
         /0.08    ,  0.57   ,  0.50  ,  0.83   ,  0.77  /
   DATA   FI0   ,  CONDFI , COEFCOND , COEFHCAP , COEFEXT &
         /0.17  ,  2.034  , 0.1172   , 1.715E+7 , 1.5     /
   DATA  EMISI , EMISNO , EMISW / 0.99 , 0.99 , 0.97 /
   DATA  ROICE ,  ROSWTR  / 913.0 ,  1025.0  /
   DATA  ROSNOW / 330.0 , 450.0 /
   DATA  Z0ICE ,  Z0W  / 1.6E-4 , 3.2E-5 /
   DATA  BASEHF / 2.0 /
   DATA   HCAPI    , HCAPW    , VHFICE   , VHFBAS   , VHFSNO   &
         /2.062E+3 , 4.088E+3 , 3.014E+8 , 2.679E+8 , 1.097E+8 /
   DATA  HSMIN , DMIX / 0.010 , 30.0 /
   DATA  XPI /3.1415927/
   DATA  DAYLNT /86400.0/
   !

!Huziy
   real, dimension(n) :: LHICE_prev, HSNOW_prev !ice thickness and snow depth at the previous time step (needed for runoff calculation)
   real, dimension(n) :: HSNOW_current, LFICE_prev
   real delta_thickness

   !parameters taken from lakepar.cdk
   real, parameter :: rhowat  = 1000.             ! density of water
   real, parameter :: rhosnow =  330.             ! density of snow
   real, parameter :: rhoice  =  917.             ! density of ice
   !
   !
   !
   real Z0MAX, QSENS,QLAT
   !      SAVE Z0MAX

   real, dimension(n) :: dummy1,dummy2,dummy3,dummy4
   real dummy5,dummy6
   !
   !**
   !
   !**  WARNING ------ the value for Z0MAX needs to be increased when coupled with WAM
   !
   !**   DATA Z0MAX / 5.E-2 /
   DATA Z0MAX / 5.E-3 /
   !
   !** ------------------------------------------------------------------
   !

   !# In offline mode the t-step 0 is (correctly) not performed
!   if (FLUVERT.eq.'SURFACE'.and.KOUNT.eq.0) return

!print *,'bussiz, kount, trnch:',bussiz, kount, trnch
!print *,'ptsurfsiz, ptsurf(ptsurfsiz):',ptsurfsiz, ptsurf(ptsurfsiz)
!print *,'n, m, nk, dt:',n, m, nk, dt
!print *,'bus(bussiz):',bus(1:bussiz)

   !
   SURFLEN = M
   !
   ! Input fields
   hu         (1:n) => bus( x(humoins,   1,nk)       : ) ! Input
   ps         (1:n) => bus( x(pmoins,    1,1)        : ) ! Input
   th         (1:n) => bus( x(thetaa,    1,1)        : ) ! Input
   tt         (1:n) => bus( x(tmoins,    1,nk)       : ) ! Input
   uu         (1:n) => bus( x(umoins,    1,nk)       : ) ! Input
   vv         (1:n) => bus( x(vmoins,    1,nk)       : ) ! Input
   zdlat      (1:n) => bus( x(dlat,      1,1)        : ) ! Input
   zdlon      (1:n) => bus( x(dlon,      1,1)        : ) ! Not used
   zfcor      (1:n) => bus( x(fcor,      1,1)        : ) ! Input
   zrainrate  (1:n) => bus( x(rainrate,  1,1)        : ) ! Input
   zsnowrate  (1:n) => bus( x(snowrate,  1,1)        : ) ! Input
   zzusl      (1:n) => bus( x(zusl,      1,1)        : ) ! Input
   zztsl      (1:n) => bus( x(ztsl,      1,1)        : ) ! Input
   if (radslope)  then
      fsol  (1:n)   => bus( x(fluslop,   1,1)        : ) ! Input
   else
      fsol  (1:n)   => bus( x(flusolis,  1,1)        : ) ! Input
   endif
   qlwin      (1:n) => bus( x(fdsi,      1,1)        : ) ! Input
   ldepth     (1:n) => bus( x(lakedepth, 1,1)        : ) ! Input
   lfr        (1:n) => bus( x(lakefr,    1,1)        : ) ! Input
   ltransp    (1:n) => bus( x(laketransp,1,1)        : ) ! Input

   ! Input/output
   lct        (1:n) => bus( x(lakect,    1,1)        : ) ! Input/Ouput
   lfice      (1:n) => bus( x(lakefice,  1,1)        : ) ! Input/Ouput
   lhice      (1:n) => bus( x(lakehice,  1,1)        : ) ! Input/Ouput
   lhml       (1:n) => bus( x(lakehml,   1,1)        : ) ! Input/Ouput
   ltbot      (1:n) => bus( x(laketbot,  1,1)        : ) ! Input/Ouput
   ltice      (1:n) => bus( x(laketice,  1,1)        : ) ! Input/Ouput
   ltmnw      (1:n) => bus( x(laketmnw,  1,1)        : ) ! Input/Ouput
   ltwml      (1:n) => bus( x(laketwml,  1,1)        : ) ! Input/Ouput
   zsnodp     (1:n) => bus( x(snodp,     1,indx_sfc) : ) ! Input/Ouput

   ! Output
   alvis_avg  (1:n) => bus( x(alvis,     1,indx_sfc) : ) ! Output
   cmu_avg    (1:n) => bus( x(bm,        1,1)        : ) ! Output
   ctu_avg    (1:n) => bus( x(bt,        1,indx_sfc) : ) ! Output
   fc_avg     (1:n) => bus( x(fc,        1,indx_sfc) : ) ! Output
   fv_avg     (1:n) => bus( x(fv,        1,indx_sfc) : ) ! Output
   hst_avg    (1:n) => bus( x(hst,       1,indx_sfc) : ) ! Output
   ilmo_avg   (1:n) => bus( x(ilmo,      1,indx_sfc) : ) ! Output
   z0h_avg    (1:n) => bus( x(z0t,       1,indx_sfc) : ) ! Output
   z0m_avg    (1:n) => bus( x(z0,        1,indx_sfc) : ) ! Output
   zalfaq     (1:n) => bus( x(alfaq,     1,1)        : ) ! Output
   zalfat     (1:n) => bus( x(alfat,     1,1)        : ) ! Output
   zftemp     (1:n) => bus( x(ftemp,     1,indx_sfc) : ) ! Output
   zfvap      (1:n) => bus( x(fvap,      1,indx_sfc) : ) ! Output
   zqdiag     (1:n) => bus( x(qdiag,     1,1)        : ) ! Output
   zrunofftot (1:n) => bus( x(runofftot, 1,indx_sfc) : ) ! Output
   ztsurf     (1:n) => bus( x(tsurf,     1,indx_sfc) : ) ! Output
   qs         (1:n) => bus( x(qsurf,     1,indx_sfc) : ) ! Output
   ztsrad     (1:n) => bus( x(tsrad,     1,1)        : ) ! Output
   ztdiag     (1:n) => bus( x(tdiag,     1,1)        : ) ! Output
   zudiag     (1:n) => bus( x(udiag,     1,1)        : ) ! Output
   zvdiag     (1:n) => bus( x(vdiag,     1,1)        : ) ! Output
   zfrv_avg   (1:n) => bus( x(frv,       1,indx_sfc) : ) ! Output
   zfrv_ice   (1:n) => bus( x(frv_li,    1,1)        : ) ! Output
   zfrv_wat   (1:n) => bus( x(frv_lw,    1,1)        : ) ! Output

!   ztt2m      (1:n) => bus( x(tt2m,      1,indx_sfc) : ) ! Output


!* 1.     Preliminaries
!  --------------------
!  Calculate windspeed (vmod)

   i = sl_prelim(tt,hu,uu,vv,ps,zzusl,spd_air=vmod,dir_air=vdir, &
                 min_wind_speed=sqrt(6.25))
!                 min_wind_speed=sqrt(vamin))
   if (i /= SL_OK) then
      print*, 'Aborting in flake_main() because of error returned by sl_prelim()'
      stop
   endif

!print *,'vmod,vdir:',vmod,vdir

   !
   !
   !********************************************************************************
   !       A)     SURFACE FLUXES AND ALBEDO FOR OPEN WATER AND ICE COVERED FRACTIONS
   !********************************************************************************

!Huziy remember ice thickness and snow height on the ice for runoff calculation
   LHICE_prev = LHICE
   HSNOW_prev = ZSNODP
   LFICE_prev = LFICE


   ! A.1    Saturated specific humidity at the water and ice surface
   ! ---------------------------------------------------------------
   DO I=1,N

      IF(LFICE(I).LE.0.01) LTICE(I)=TMELS

      TS_WAT(I) = LTWML(I) ! Mixed layer temperature
      TS_ICE(I) = LTICE(I) ! Ice surface temperature

      ! Saturated specific humidity at surface
      QS_WAT(I) = FOQST (TS_WAT(I), PS(I))
      QS_ICE(I) = FOQST (TS_ICE(I), PS(I))

   END DO


   ! A.2    Calculate roughness lengths based on generalized Charnock's relation (open water)
   ! ------------------------------------------------------------------------------------------

   beta     =  0.018
   windcrit = 12.5
   ZFRVMIN  =  0.01


   ! min of zfrv to have no impact on z0m is given by sqrt(g*z0min/beta)
   ! will also solve the division by zero peroblem for z0h

   ! Beginning ---New formulation of Z0M from Daniel Deacu and for winds
   !              stronger than 12.5 m/s new formulation of Z0M from Moon.

   DO I=1,N
      ZFRV_WAT(I)=MAX(ZFRVMIN,  ZFRV_WAT(I))
      if (VMOD(I) .le. windcrit) then
         Z0M_WAT(I) =MIN( BETA*ZFRV_WAT(I)**2/GRAV + 1.65e-06/ ZFRV_WAT(I) ,Z0MAX )
      else
         Z0M_WAT(I)=MAX((0.085*(-0.56*ZFRV_WAT(I)**2+20.255*ZFRV_WAT(I)+2.458)-0.58)/1000.,Z0MIN )
      endif
   END DO

   ! End ---New formulation of Z0M from Deacu and Moon.

   !     Note:  For |lat| >= Z0TLAT(2)  Charnock's relation is used
   !            For |lat| <= Z0TLAT(1)  Z0HCON is used.
   !            For Z0TLAT(1) < |lat| < Z0TLAT(2)
   !            we do a linear interpolation between Charnock and Z0HCON.

   ! Beginning --New formulation of Z0H from Daniel Deacu
   DO I=1,N
      Z0H_WAT(I) = MIN(2.e-05/ZFRV_WAT(I), 1.e-04)
   END DO
   ! End --New formulation of Z0H from Daniel Deacu

   DO I=1,N
      IF (ABS(ZDLAT(I)) .GE. Z0TLAT(2)) THEN
         Z0H_WAT(I) = Z0H_WAT(I)
      ELSE IF (ABS(ZDLAT(I)) .LE. Z0TLAT(1)) THEN
         Z0H_WAT(I) = Z0HCON
      ELSE
      Z0H_WAT(I)=( ((ABS(ZDLAT(I))-Z0TLAT(1))/(Z0TLAT(2)-Z0TLAT(1))) &
                    *(Z0H_WAT(I)-Z0HCON) ) + Z0HCON
      ENDIF
   END DO

   ! A.3    Calculate the surface transfer coefficient and fluxes (open water)
   !              and estimate diagnostic-level quantities (open water)
   ! -------------------------------------------------------------------------

   i = sl_sfclayer(th,hu,vmod,vdir,zzusl,zztsl,ts_wat,qs_wat,z0m_wat,z0h_wat,zdlat,zfcor, &
                   coefm=cmu_wat, coeft=ctu_wat, flux_t=zftemp_wat, flux_q=zfvap_wat, &
                   ilmo=ilmo_wat, ue=zfrv_wat, h=hst_wat, &
                   hghtm_diag=zu, hghtt_diag=zt, tdiaglim=LAKE_TDIAGLIM, &
                   t_diag=ztdiag_wat, q_diag=zqdiag_wat, u_diag=zudiag_wat, v_diag=zvdiag_wat)

!print *,'cmu_wat:',cmu_wat
!print *,'ctu_wat:',ctu_wat
!print *,'zftemp_wat:',zftemp_wat
!print *,'zfvap_wat:',zfvap_wat
!print *,'ilmo_wat:',ilmo_wat
!print *,'zfrv_wat:',zfrv_wat
!print *,'hst_wat:',hst_wat
!print *,'ztdiag_wat:',ztdiag_wat
!print *,'zqdiag_wat:',zqdiag_wat
!print *,'zudiag_wat:',zudiag_wat
!print *,'zvdiag_wat:',zvdiag_wat

   if (i /= SL_OK) then
      print*, 'Aborting in flake_main() because of error returned by sl_sfclayer()'
      stop
   endif


!   CALL FLXSURF3( CMU_WAT, CTU_WAT, SCR1, ZFTEMP_WAT, ZFVAP_WAT, &
!                  ILMO_WAT, ZFRV_WAT, ZFCOR, TH, HU,             &
!                  ZZUSL, ZZTSL, VMOD, TS_WAT, QS_WAT, HST_WAT,   &
!                  Z0M_WAT, Z0H_WAT,SCR2, SCR3, SCR4, SCR5, N )

!   CALL DIASURF2(ZUDIAG_WAT, ZVDIAG_WAT, ZTDIAG_WAT, ZQDIAG_WAT, &
!                 N, UU, VV, TS_WAT, QS_WAT,                      &
!                 Z0M_WAT, Z0H_WAT, ILMO_WAT, ZZUSL,              &
!                 HST_WAT, ZFRV_WAT, ZFTEMP_WAT, ZFVAP_WAT,       &
!                 ZUN, ZTN, ZDLAT)

   ! 4.     Finalize the fluxes (open water)
   ! ---------------------------------------

   !VDIR NODEP
   DO I=1,N

      ZTSURF_WAT   (I) = TS_WAT (I)
      ZTSRAD_WAT   (I) = TS_WAT (I)

      ZALFAT_WAT   (I) = - CTU_WAT(I) * ( TS_WAT(I)-TH(I) )
      ZALFAQ_WAT   (I) = - CTU_WAT(I) * ( QS_WAT(I)-HU(I) )
      IF (.NOT.IMPFLX) CTU_WAT(I) = 0.0
      RHO = PS(I)/(RGASD * ZTDIAG_WAT(I)*(1.+DELTA*ZQDIAG_WAT(I)))
      FC_WAT(I) = -CPD *RHO*ZALFAT_WAT(I)
      FV_WAT(I) = -CHLC*RHO*ZALFAQ_WAT(I)

      IF (IMPFLX) THEN
         ZALFAT_WAT   (I) = - CTU_WAT(I) *  TS_WAT(I)
         ZALFAQ_WAT   (I) = - CTU_WAT(I) *  QS_WAT(I)
      ENDIF

   END DO

   !*******************************************************************************
   !       B)      ICE-COVERED LAKE FRACTION
   !*******************************************************************************

   !       B.1    Preliminaries (ice)
   !       --------------------------

   DO I=1,N

      IF( ICEMELT ) THEN
         IF( LHICE(I) .GE. HIMIN ) THEN
            Z0M_ICE(I) = Z0ICE ! Roughness lengths for the surface (ice/water)
         ELSE
            Z0M_ICE(I) = Z0W
   !        ZSNODP(I) = 0.0 ! Remove snow if ice is too thin
         ENDIF
      ELSE
         HICE_ALB(I) = MAX ( LHICE(I) , HIMIN ) ! minimum ice thickness
   !     ZSNODP(I) = 0.0
         Z0M_ICE(I) = Z0ICE
      ENDIF

      Z0H_ICE(I) = Z0M_ICE(I)

   END DO

   !       B.2    Calculate the drag and heat coefficients (ice)
   !              and estimate diagnostic-level quantities (ice)
   !       -----------------------------------------------------

   i = sl_sfclayer(th,hu,vmod,vdir,zzusl,zztsl,ts_ice,qs_ice,z0m_ice,z0h_ice,zdlat,zfcor, &
                   coefm=cmu_ice, coeft=ctu_ice, flux_t=zftemp_ice, flux_q=zfvap_ice, &
                   ilmo=ilmo_ice, ue=zfrv_ice, h=hst_ice, &
                   hghtm_diag=zu, hghtt_diag=zt, tdiaglim=LAKE_TDIAGLIM, &
                   t_diag=ztdiag_ice, q_diag=zqdiag_ice, u_diag=zudiag_ice, v_diag=zvdiag_ice)
   if (i /= SL_OK) then
      print*, 'Aborting in flake_main() because of error returned by sl_sfclayer()'
      stop
   endif


!   CALL FLXSURF3( CMU_ICE, CTU_ICE, SCR1, ZFTEMP_ICE, ZFVAP_ICE, ILMO_ICE, &
!                  ZFRV_ICE, ZFCOR, TH, HU, ZZUSL, ZZTSL, VMOD, TS_ICE,     &
!                  QS_ICE, HST_ICE, Z0M_ICE, Z0H_ICE,                       &
!                  SCR2, SCR3, SCR4, SCR5, N)

!   CALL DIASURF2(ZUDIAG_ICE, ZVDIAG_ICE, ZTDIAG_ICE, ZQDIAG_ICE,           &
!                 N, UU, VV, TS_ICE, QS_ICE,                                &
!                 Z0M_ICE, Z0H_ICE, ILMO_ICE, ZZUSL,                        &
!                 HST_ICE, ZFRV_ICE, ZFTEMP_ICE, ZFVAP_ICE,                 &
!                 ZUN, ZTN, ZDLAT)

!VDIR NODEP
   DO I=1,N

   !     IF( HICE_ALB(I).LT.HIMIN ) ZSNODP(I) = 0.0 ! remove snow if ice is too thin
      ZTSURF_ICE   (I) = TS_ICE(I)
      ZTSRAD_ICE   (I) = TS_ICE(I)

      ZALFAT_ICE   (I) = - CTU_ICE(I) * ( TS_ICE (I) - TH(I) )
      ZALFAQ_ICE   (I) = - CTU_ICE(I) * ( QS_ICE (I) - HU(I) )
      IF (.NOT.IMPFLX) CTU_ICE (I) = 0.
      RHO = PS(I)/(RGASD * ZTDIAG_ICE(I)*(1.+DELTA*ZQDIAG_ICE(I)))
      FC_ICE(I) = -CPD *RHO*ZALFAT_ICE(I)
      FV_ICE(I) = -(CHLC+CHLF)*RHO*ZALFAQ_ICE(I)

      IF (IMPFLX) THEN
         ZALFAT_ICE   (I) = - CTU_ICE(I) *  TS_ICE(I)
         ZALFAQ_ICE   (I) = - CTU_ICE(I) *  QS_ICE(I)
      ENDIF

   END DO

   !	     SURFACE-AVERAGED VALUES

   DO I=1,N
      WA_AVG(I)=(((ZUDIAG_WAT(I)**2+(ZVDIAG_WAT(I))**2)*(1.-LFICE(I)))  &
               + ((ZUDIAG_ICE(I)**2+(ZVDIAG_ICE(I))**2)*    LFICE(I)))**0.5

      ZQDIAG_TMP(I)=ZQDIAG_WAT(I)*(1.-LFICE(I))+ZQDIAG_ICE(I)*LFICE(I)
   END DO


   !     FLAKE MODEL
   !     ---------------------------

   DO I=1,N                  !* INPUT DATA FORMATION

      SWD = MAX(0.0,FSOL(I)) !* SW RADIATION DOWNWARD, AS IN FSS
      WIND=WA_AVG(I)         !* WIND FORCE @10m
   !  HICE=LHICE(I)          !* ICE THICKNESS
      HSNOW=0.0              !* NO SNOW ON ICE

      omega_earth = 7.29E-05
      CORIO=2.*omega_earth*SIN(zdlat(I))

   !    SURFACE TEMPERATURE: Ice temp. if ice present, else GT

      IF(LHICE(I).GT.0.0) THEN
         T_sfc_in  = LTICE(I)
      ELSE
         T_sfc_in  = LTWML(I)
      ENDIF

   !    FLAKE PARAMETERS READ FROM ARRAYS

      T_ice_in = LTICE(I)
      T_snow_in= T_ice_in
      T_mnw_in = LTMNW(I)
      T_wML_in = LTWML(I)
      T_bot_in = LTBOT(I)
      C_T_in   = LCT(I)
      h_ML_in  = LHML(I)

      IF(LHICE(I).GT.0.0) THEN
         LFICE(I) =1.0
         FC_AVG(I)=FC_ICE(I)
         FV_AVG(I)=FV_ICE(I)
         T_air_in =ZTDIAG_ICE(I)
         ZQDIAG_TMP(I)=ZQDIAG_ICE(I)
      ELSE
         LFICE(I) =0.0
         FC_AVG(I)=FC_WAT(I)
         FV_AVG(I)=FV_WAT(I)
         T_air_in =ZTDIAG_WAT(I)
         ZQDIAG_TMP(I)=ZQDIAG_WAT(I)
      ENDIF

   !    sediment layer(not used)
      T_B1_in = TMELS
      H_B1_in = 10.0

      ! Always keep snow depth at zero in FLake
      HSNOW    = 0.0


   !  Running FLake

      QSENS=0.0
      QLAT=0.0
      albedo=0.0
      dummy5=0.0
      dummy6=484.0
      FDEPTH=LDEPTH(I)

!print*,'FLake kount,trnch,i:',kount,trnch,i
!print*,'T_wML_in:',T_wML_in
!print*,'T_bot_in:',T_bot_in
!print*,'T_mnw_in:',T_mnw_in
!print*,'h_ML_in:',h_ML_in
!print*,'C_T_in:',C_T_in
!print*,'HSNOW:',HSNOW
!print*,'LHICE(I):',LHICE(I)
!print*,'T_ice_in:',T_ice_in
!print*,'T_snow_in:',T_snow_in
!print*,'T_sfc_in:',T_sfc_in
!print*,'LTRANSP(i):',LTRANSP(i)
!print*,'FDEPTH:',FDEPTH

!      CALL flaket_lakes ( HSNOW, SWD, QLWIN(I),zu, zt, WIND,        &
!      CALL flaket_lakes ( HSNOW, SWD, QLWIN(I),zu, 2., WIND,        &
!                          T_air_in,ZQDIAG_TMP(I),PS(I),LTRANSP(I),  &
      CALL flaket_lakes ( HSNOW, SWD, QLWIN(I),zzusl,zztsl, vmod(i),&
                          tt(i),hu(I),PS(I),LTRANSP(I),  &
                          FDEPTH,1.0E3,10.0,TMELS+4.0,CORIO,DELT,   &
                          T_snow_in, T_ice_in, T_mnw_in, T_wML_in,  &
                          T_bot_in,T_B1_in,C_T_in,0.0,LHICE(I),     &
                          h_ML_in,H_B1_in,T_sfc_in,dummy6,albedo,   &
                          T_snow_out,T_ice_out,T_mnw_out,T_wML_out, &
                          T_bot_out,T_B1_out,h_snow_out, h_ice_out, &
                          h_ML_out, H_B1_out, QSENS, dummy5,QLAT,   &
                          T_sfc_out, C_T_out,trnch,n,m,i,LFR(I) )

!  if (trnch == 131 .and. i == 32) then ! debug !
!     print *,"Apres appel a flaket_lakes pour trnch=131 et i=32..."
!     print *,"HSNOW, SWD, QLWIN(I),WIND",HSNOW, SWD, QLWIN(I),WIND
!     print *,"T_air_in,ZQDIAG_TMP(I),PS(I),LTRANSP(I)",T_air_in,ZQDIAG_TMP(I),PS(I),LTRANSP(I)
!     print *,"FDEPTH,TMELS+4.0,CORIO,DELT",FDEPTH,TMELS+4.0,CORIO,DELT
!     print *,"T_snow_in, T_ice_in, T_mnw_in, T_wML_in",T_snow_in, T_ice_in, T_mnw_in, T_wML_in
!     print *,"T_bot_in,T_B1_in,C_T_in,LHICE(I)",T_bot_in,T_B1_in,C_T_in,LHICE(I)
!     print *,"h_ML_in,H_B1_in,T_sfc_in,dummy6,albedo",h_ML_in,H_B1_in,T_sfc_in,dummy6,albedo
!     print *,"T_snow_out,T_ice_out,T_mnw_out,T_wML_out",T_snow_out,T_ice_out,T_mnw_out,T_wML_out
!     print *,"T_bot_out,T_B1_out,h_snow_out, h_ice_out",T_bot_out,T_B1_out,h_snow_out, h_ice_out
!     print *,"h_ML_out, H_B1_out, QSENS, dummy5,QLAT",h_ML_out, H_B1_out, QSENS, dummy5,QLAT
!     print *,"T_sfc_out, C_T_out,trnch,n,m,i,LFR(I)",T_sfc_out, C_T_out,trnch,n,m,i,LFR(I)
!  endif

!       Output parameters

      FC_AVG(I)=QSENS
      FV_AVG(I)=QLAT

      LTICE(I)  = T_ice_out
      LTMNW(I)  = T_mnw_out
      LTWML(I)  = T_wML_out
      LTBOT(I)  = T_bot_out
      LCT(I)    = C_T_out
      LHML(I)   = h_ML_out

      HSNOW=0.0 !* No snow in FLake AMTV

      LHICE(I)=h_ice_out
      ALVIS_AVG(I)=albedo

      IF(LHICE(I).GT.0.0) THEN
         LFICE(I)=1.0
         ALVIS_ICE(I)=albedo
         ALVIS_WAT(I)=0.199
         FC_ICE(I)=FC_AVG(I)
         FV_ICE(I)=FV_AVG(I)
      ELSE
         LFICE(I)=0.0
         ALVIS_WAT(I)=albedo
         ALVIS_ICE(I)=0.199
         FC_WAT(I)=FC_AVG(I)
         FV_WAT(I)=FV_AVG(I)
      ENDIF

      !Huziy
      HSNOW_current(I) = HSNOW

   ENDDO

   ! FLAKE MODEL FINISHED
   ! --------------------


   ! SURFACE TEMPERATURES
   ! --------------------
!Sasha: set temperatures from NEMO, where required
!Sasha all  the temperatures here should be in Kelvins
   DO I=1,N
      TS_WAT(I) = LTWML(I)
      TS_ICE(I) = LTICE(I)
   ENDDO

   !*************************************************************
   !       D)     AFTER THE LAKE MODEL: SURFACE FLUXES AND ALBEDO
   !              FOR OPEN WATER AND ICE COVERED FRACTIONS
   !*************************************************************


   !       D.1    Saturated specific humidity at the water surface (open water)
   !       --------------------------------------------------------------------

   DO I=1,N

      IF(LFICE(I).EQ.0.0) LTICE(I)=TMELS

      QS_WAT(I)= FOQST(TS_WAT(I),PS(I))
   !  QS_WAT(I)= QS(I)

   END DO

   !       D.2    Calculate the surface transfer coefficient and fluxes (open water)
   !              and estimate diagnostic-level quantities (open water)
   !       -------------------------------------------------------------------------


   i = sl_sfclayer(th,hu,vmod,vdir,zzusl,zztsl,ts_wat,qs_wat,z0m_wat,z0h_wat,zdlat,zfcor, &
                   coefm=cmu_wat, coeft=ctu_wat, flux_t=zftemp_wat, flux_q=zfvap_wat, &
                   ilmo=ilmo_wat, ue=zfrv_wat, h=hst_wat, &
                   hghtm_diag=zu, hghtt_diag=zt, tdiaglim=LAKE_TDIAGLIM, &
                   t_diag=ztdiag_wat, q_diag=zqdiag_wat, u_diag=zudiag_wat, v_diag=zvdiag_wat)
   if (i /= SL_OK) then
      print*, 'Aborting in flake_main() because of error returned by sl_sfclayer()'
      stop
   endif


!   CALL FLXSURF3( CMU_WAT, CTU_WAT, SCR1, ZFTEMP_WAT, ZFVAP_WAT, &
!                  ILMO_WAT, ZFRV_WAT, ZFCOR, TH, HU,             &
!                  ZZUSL, ZZTSL, VMOD, TS_WAT, QS_WAT, HST_WAT,   &
!                  Z0M_WAT, Z0H_WAT,SCR2, SCR3, SCR4, SCR5, N )


!   CALL DIASURF2(ZUDIAG_WAT, ZVDIAG_WAT, ZTDIAG_WAT, ZQDIAG_WAT, &
!                 N, UU, VV, TS_WAT, QS_WAT,                      &
!                 Z0M_WAT, Z0H_WAT, ILMO_WAT, ZZUSL,              &
!                 HST_WAT, ZFRV_WAT, ZFTEMP_WAT, ZFVAP_WAT,       &
!                 ZUN, ZTN, ZDLAT)

   !       D.3    Finalize the fluxes (open water)
   !       ---------------------------------------

   !VDIR NODEP
   DO I=1,N

      ZTSURF_WAT   (I) = TS_WAT (I)
      ZTSRAD_WAT   (I) = TS_WAT (I)

   !  ZALFAT_WAT   (I) = - CTU_WAT(I) * ( TS_WAT(I)-TH(I) )
   !  ZALFAQ_WAT   (I) = - CTU_WAT(I) * ( QS_WAT(I)-HU(I) )
      IF (.NOT.IMPFLX) CTU_WAT(I) = 0.
      RHO = PS(I)/(RGASD * ZTDIAG_WAT(I)*(1.+DELTA*ZQDIAG_WAT(I)))
   !  FC_WAT(I) = -CPD *RHO*ZALFAT_WAT(I)
   !  FV_WAT(I) = -CHLC*RHO*ZALFAQ_WAT(I)

      ZALFAT_WAT(I)=-FC_WAT(I)/(CPD *RHO)
      ZALFAQ_WAT(I)=-FV_WAT(I)/(CHLC*RHO)
!Huziy uncommented RUNOFF_WAT (convert to kg/(m**2 * s))
      ZRUNOFF_WAT (I) = (rhowat * (ZRAINRATE(I) + ZSNOWRATE(I)) &
                      + RHO*ZALFAQ_WAT(I))

     IF (IMPFLX) THEN
         ZALFAT_WAT   (I) = - CTU_WAT(I) *  TS_WAT(I)
         ZALFAQ_WAT   (I) = - CTU_WAT(I) *  QS_WAT(I)
   !     ZRUNOFF_WAT (I) = (1000.*(ZRAINRATE(I) +ZSNOWRATE(I)) &
   !     + RHO*ZALFAQ_WAT(I))
      ENDIF

   END DO

   !    D.4    Preliminaries (ICE-COVERED LAKE FRACTION)
   !    ------------------------------------------------

   DO I=1,N

      IF( ICEMELT ) THEN
         IF( LHICE(I) .GE. HIMIN ) THEN !* Roughness lengths for the surface (ice/water)
            Z0M_ICE(I) = Z0ICE
         ELSE
            Z0M_ICE(I) = Z0W
   !        ZSNODP(I) = 0.0 !* Remove snow if ice is too thin
         ENDIF
      ELSE
         HICE_ALB(I) = MAX ( LHICE(I) , HIMIN ) !* Minimum ice thickness
         Z0M_ICE(I) = Z0ICE
      ENDIF

      Z0H_ICE(I) = Z0M_ICE(I)
   !  ZSNODP(I) = 0.0

   END DO

   !    D.5    Calculate the drag and heat coefficients (ice)
   !           and estimate diagnostic-level quantities (ice)
   !    -----------------------------------------------------

   DO I=1,N
   !                               Saturated specific humidity at surface
      QS_ICE(I)  = FOQST ( TS_ICE  (I), PS(I) )
   END DO

   i = sl_sfclayer(th,hu,vmod,vdir,zzusl,zztsl,ts_ice,qs_ice,z0m_ice,z0h_ice,zdlat,zfcor, &
                   coefm=cmu_ice, coeft=ctu_ice, flux_t=zftemp_ice, flux_q=zfvap_ice, &
                   ilmo=ilmo_ice, ue=zfrv_ice, h=hst_ice, &
                   hghtm_diag=zu, hghtt_diag=zt, tdiaglim=LAKE_TDIAGLIM, &
                   t_diag=ztdiag_ice, q_diag=zqdiag_ice, u_diag=zudiag_ice, v_diag=zvdiag_ice)
   if (i /= SL_OK) then
      print*, 'Aborting in flake_main() because of error returned by sl_sfclayer()'
      stop
   endif


!   CALL FLXSURF3( CMU_ICE, CTU_ICE, SCR1, ZFTEMP_ICE, ZFVAP_ICE, ILMO_ICE, &
!                  ZFRV_ICE, ZFCOR, TH, HU, ZZUSL, ZZTSL, VMOD, TS_ICE,     &
!                  QS_ICE, HST_ICE, Z0M_ICE, Z0H_ICE,                       &
!                  SCR2, SCR3, SCR4, SCR5, N)

!   CALL DIASURF2(ZUDIAG_ICE, ZVDIAG_ICE, ZTDIAG_ICE, ZQDIAG_ICE,           &
!                 N, UU, VV, TS_ICE, QS_ICE,                                &
!                 Z0M_ICE, Z0H_ICE, ILMO_ICE, ZZUSL,                        &
!                 HST_ICE, ZFRV_ICE, ZFTEMP_ICE, ZFVAP_ICE,                 &
!                 ZUN, ZTN, ZDLAT)

!VDIR NODEP
   DO I=1,N

      ZTSURF_ICE   (I) = TS_ICE(I)
      ZTSRAD_ICE   (I) = TS_ICE(I)

      IF (.NOT.IMPFLX) CTU_ICE (I) = 0.
      RHO = PS(I)/(RGASD * ZTDIAG_ICE(I)*(1.+DELTA*ZQDIAG_ICE(I)))

!      FC_ICE(I) = -CPD *RHO*ZALFAT_ICE(I)
!      FV_ICE(I) = -(CHLC+CHLF)*RHO*ZALFAQ_ICE(I)

      ZALFAT_ICE(I)=-FC_ICE(I)/(CPD *RHO)
      ZALFAQ_ICE(I)=-FV_ICE(I)/((CHLC+CHLF)*RHO) !Huziy fix: was (CHLC+CHLF*RHO)

!      ZRUNOFF_ICE (I) = (1000.*(ZRAINRATE(I) +ZSNOWRATE(I)) &
!                         + RHO*ZALFAQ_ICE(I))

!Huziy
      ZRUNOFF_ICE (I) = (1000.*(ZRAINRATE(I)) + RHO*ZALFAQ_ICE(I))
!        ZRUNOFF_ICE (I) = 0.0

      IF (IMPFLX) THEN
         ZALFAT_ICE   (I) = - CTU_ICE(I) *  TS_ICE(I)
         ZALFAQ_ICE   (I) = - CTU_ICE(I) *  QS_ICE(I)
!         ZRUNOFF_ICE (I) = (1000.*(ZRAINRATE(I) +ZSNOWRATE(I)) + RHO*ZALFAQ_ICE(I))
      ENDIF


!Huziy account for the ice and snow melting, using ROSNOW(2), since runoff happens only when it melts
!and the melting snow is of bigger density
!print *,'HSNOW_current,HSNOW_prev:',HSNOW_current(I),HSNOW_prev(I)
!print *,'ROSNOW,rhowat:',ROSNOW(2),rhowat
     delta_thickness = (HSNOW_current(I) - HSNOW_prev(I)) * (ROSNOW(2) / rhowat)
!print *,'delta_thickness:',delta_thickness
     if (delta_thickness < 0) then
       !print*, "delta_thickness (snow) = ", delta_thickness
       ZRUNOFF_ICE (I) =  ZRUNOFF_ICE (I) - 1000.0 * delta_thickness / dt
     endif

     delta_thickness = (LHICE(I) - LHICE_prev(I)) * (ROICE / rhowat)
     !print*, "delta_thickness (ice) = ", delta_thickness
     !*1000.0 - converts m/s to mm/s
     ZRUNOFF_ICE (I) =  ZRUNOFF_ICE (I) - 1000.0 * delta_thickness / dt

   END DO

   ! OPEN WATER + ICE: AGGREGATION

   DO I=1,N

   ! LINEAR AGGREGATION (OPEN WATER + ICE):

      ZALFAQ  (I)=ZALFAQ_WAT(I)*(1.-LFICE(I))+ZALFAQ_ICE(I)*LFICE(I) !*ALFAQ
      ZALFAT  (I)=ZALFAT_WAT(I)*(1.-LFICE(I))+ZALFAT_ICE(I)*LFICE(I) !*ALFAT
      CTU_AVG (I)=CTU_WAT   (I)*(1.-LFICE(I))+CTU_ICE   (I)*LFICE(I) !*BT
      FC_AVG  (I)=FC_WAT    (I)*(1.-LFICE(I))+FC_ICE    (I)*LFICE(I) !*FC
      ZFRV_AVG(I)=ZFRV_WAT  (I)*(1.-LFICE(I))+ZFRV_ICE  (I)*LFICE(I) !*FRV
      ZFTEMP  (I)=ZFTEMP_WAT(I)*(1.-LFICE(I))+ZFTEMP_ICE(I)*LFICE(I) !*FTEMP
      FV_AVG  (I)=FV_WAT    (I)*(1.-LFICE(I))+FV_ICE    (I)*LFICE(I) !*FV
      ZFVAP   (I)=ZFVAP_WAT (I)*(1.-LFICE(I))+ZFVAP_ICE (I)*LFICE(I) !*FVAP
      HST_AVG (I)=HST_WAT   (I)*(1.-LFICE(I))+HST_ICE   (I)*LFICE(I) !*HST
      ILMO_AVG(I)=ILMO_WAT  (I)*(1.-LFICE(I))+ILMO_ICE  (I)*LFICE(I) !*ILMO
      QS      (I)=QS_WAT    (I)*(1.-LFICE(I))+QS_ICE    (I)*LFICE(I) !*QSURF
      ZSNODP  (I)=ZSNODP    (I) !*SNODP: NO SNOW OVER ICE (YET) AND OPEN WATER
      ZTSURF  (I)=ZTSURF_WAT(I)*(1.-LFICE(I))+ZTSURF_ICE(I)*LFICE(I) !*TSURF
      ZQDIAG  (I)=ZQDIAG_WAT(I)*(1.-LFICE(I))+ZQDIAG_ICE(I)*LFICE(I) !*QDIAG
      CMU_AVG (I)=CMU_WAT   (I)*(1.-LFICE(I))+CMU_ICE   (I)*LFICE(I) !*BM
      ZTDIAG  (I)=ZTDIAG_WAT(I)*(1.-LFICE(I))+ZTDIAG_ICE(I)*LFICE(I) !*TDIAG
      ZVDIAG  (I)=ZVDIAG_WAT(I)*(1.-LFICE(I))+ZVDIAG_ICE(I)*LFICE(I) !*VDIAG
      ZUDIAG  (I)=ZUDIAG_WAT(I)*(1.-LFICE(I))+ZUDIAG_ICE(I)*LFICE(I) !*UDIAG
!Huziy changed, was =0.0, maybe uncomment later

!print *,'ZRUNOFF_WAT,ZRUNOFF_ICE,LFICE:',ZRUNOFF_WAT(I),ZRUNOFF_ICE(I),LFICE(I)

      ZRUNOFFTOT(I)=ZRUNOFF_WAT(I)*(1.-LFICE(I))+ZRUNOFF_ICE(I)*LFICE(I)

!    NON-LINEAR AGGREGATION:

      ZTSRAD(I)=((ZTSRAD_WAT(I))**4.0*(1.-LFICE(I))+(ZTSRAD_ICE(I))**4.0*LFICE(I))**0.25 !*TSRAD

      Z0M_AVG(I) = ALOG(MAX(1.e-10,Z0M_WAT(I)))*(1.-LFICE(I)) &
                 + ALOG(MAX(1.e-10,Z0M_ICE(I)))*LFICE(I)      !*Z0M
      Z0M_AVG(I) = EXP( Z0M_AVG(I))

      Z0H_AVG(I) = ALOG(MAX(1.e-10,Z0H_WAT(I)))*(1.-LFICE(I)) &
                 + ALOG(MAX(1.e-10,Z0H_ICE(I)))*LFICE(I)      !*Z0T
      Z0H_AVG(I) = EXP( Z0H_AVG(I))

   END DO


!  Copy surface averaged fields into fields for each surface fraction
!   zTT2m(:)  = ztdiag(:)                ! instantaneos

!    FILL THE ARRAYS TO BE AGGREGATED LATER IN S/R AGREGE
   CALL FILLAGG( BUS,BUSSIZ, PTSURF,PTSURFSIZ, INDX_LAKE, SURFLEN )

   RETURN

1001     FORMAT('AMTVA BEFORE01 ',5(i8,' '),8(e14.8,' '))
1002     FORMAT('AMTVA BEFORE02 ',5(i8,' '),8(e14.8,' '))
1003     FORMAT('AMTVA BEFORE03 ',5(i8,' '),8(e14.8,' '))
1004     FORMAT('AMTVA BEFORE04 ',5(i8,' '),8(e14.8,' '))
1005     FORMAT('AMTVA BEFORE05 ',5(i8,' '),8(e14.8,' '))
1006     FORMAT('AMTVA BEFORE06 ',5(i8,' '),8(e14.8,' '))
1007     FORMAT('AMTVA BEFORE07 ',5(i8,' '),8(e14.8,' '))
1008     FORMAT('AMTVA BEFORE08 ',5(i8,' '),8(e14.8,' '))
1009     FORMAT('AMTVA BEFORE09 ',5(i8,' '),8(e14.8,' '))
1010     FORMAT('AMTVA BEFORE10 ',5(i8,' '),8(e14.8,' '))
1011     FORMAT('AMTVA BEFORE11 ',5(i8,' '),8(e14.8,' '))
1012     FORMAT('AMTVA BEFORE12 ',6(i8,' '),15(e14.8,' '))
1013     FORMAT('AMTVA BEFORE13 ',5(i8,' '),8(e14.8,' '))



1101     FORMAT('AMTVI AFTER 01 ',5(i8,' '),8(e14.8,' '))
1102     FORMAT('AMTVA AFTER 02 ',5(i8,' '),8(e14.8,' '))
1103     FORMAT('AMTVA AFTER 03 ',5(i8,' '),8(e14.8,' '))
1104     FORMAT('AMTVA AFTER 04 ',5(i8,' '),8(e14.8,' '))
1105     FORMAT('AMTVA AFTER 05 ',5(i8,' '),8(e14.8,' '))
1106     FORMAT('AMTVA AFTER 06 ',5(i8,' '),8(e14.8,' '))
1107     FORMAT('AMTVA AFTER 07 ',5(i8,' '),8(e14.8,' '))
1108     FORMAT('AMTVA AFTER 08 ',5(i8,' '),8(e14.8,' '))
1109     FORMAT('AMTVA AFTER 09 ',5(i8,' '),8(e14.8,' '))
1110     FORMAT('AMTVA AFTER 10 ',5(i8,' '),8(e14.8,' '))
1111     FORMAT('AMTVA AFTER 11 ',5(i8,' '),8(e14.8,' '))
1112     FORMAT('AMTVA AFTER 12 ',6(i8,' '),15(e14.8,' '))
1113     FORMAT('AMTVA AFTER 13 ',5(i8,' '),8(e14.8,' '))

1100     FORMAT('AMTV WATER1 ',7(e14.8,' '),4(i8,' '))
1200     FORMAT('AMTV WATER2 ',7(e14.8,' '),4(i8,' '))
865      FORMAT('AMTV OUT1 ',5(i8,' '),7(e14.8,' '))


1501     FORMAT('FLAKE BEFORE01 ',5(i8,' '),8(e14.8,' '))
1502     FORMAT('FLAKE BEFORE02 ',5(i8,' '),8(e14.8,' '))
1503     FORMAT('FLAKE BEFORE03 ',5(i8,' '),8(e14.8,' '))
1504     FORMAT('FLAKE BEFORE04 ',5(i8,' '),8(e14.8,' '))
1505     FORMAT('FLAKE BEFORE05 ',5(i8,' '),8(e14.8,' '))
1506     FORMAT('FLAKE BEFORE06 ',5(i8,' '),8(e14.8,' '))

1601     FORMAT('FLAKE AFTER_01 ',5(i8,' '),8(e14.8,' '))
1602     FORMAT('FLAKE AFTER_02 ',5(i8,' '),8(e14.8,' '))
1603     FORMAT('FLAKE AFTER_03 ',5(i8,' '),8(e14.8,' '))
1604     FORMAT('FLAKE AFTER_04 ',5(i8,' '),8(e14.8,' '))
1605     FORMAT('FLAKE AFTER_05 ',5(i8,' '),8(e14.8,' '))
1606     FORMAT('FLAKE AFTER_06 ',5(i8,' '),8(e14.8,' '))

1609     FORMAT('FLAKEAGG1 ',1(i8,' '),8(e14.8,' '))
1610     FORMAT('FLAKEAGG2 ',1(i8,' '),8(e14.8,' '))

1621     FORMAT('FLAKE AFTER_01 ',1(i8,' '),18(e14.8,' '))

4604     FORMAT('FLAKEIAFTER 05 ',8(e14.8,' '))
4605     FORMAT('FLAKEIAFTER 06 ',10(e14.8,' '))

END
