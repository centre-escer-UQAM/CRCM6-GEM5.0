!---------------------------------- LICENCE BEGIN -------------------------------
! GEM - Library of kernel routines for the GEM numerical atmospheric model
! Copyright (C) 1990-2010 - Division de Recherche en Prevision Numerique
!                       Environnement Canada
! This library is free software; you can redistribute it and/or modify it
! under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, version 2.1 of the License. This library is
! distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
! without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
! PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
! You should have received a copy of the GNU Lesser General Public License
! along with this library; if not, write to the Free Software Foundation, Inc.,
! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
!---------------------------------- LICENCE END ---------------------------------

!  s/r set_zeta    - Generates A and B of the hybrid coordinate
!                    Also sets Z and other related vertical parameters.
!
      subroutine set_zeta2( F_hybuser, Nk )
      use vGrid_Descriptors, only: vgrid_descriptor,vgd_new,vgd_get,vgd_put,&
                                   vgd_levels,VGD_OK,VGD_ERROR,vgd_print,vgd_free
      use vgrid_wb, only: vgrid_wb_put
      use gmm_pw
      use grid_options
      use gem_options
      use glb_ld
      use cstv
      use lun
      use dimout
      use out_mod
      use levels
      use ver
      use wb_itf_mod
      implicit none
#include <arch_specific.hf>

      integer Nk
      real, dimension(Nk) :: F_hybuser        !user-specified hybrid coordinate values
!
! authors
!      A. Plante & C. Girard - CMC - janvier 2008
!
! revision
!
! v4_00 - Plante & Girard   - Log-hydro-pressure coord on Charney-Phillips grid
! v4_4  - Plante - add standard pressure profils for physics.
!
! object
!    To return A, B parameters for momentum and thermodynamic levels
!    These levels are used for DYNAMICAL CALCULATIONS in the models !
!
!    Also to return other parameters related to the vertical discretization
!
!         Z, dZ, 1/dZ, dBdZ, etc
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! n staggered MOMENTUM & THERMO **LEVELS** !!!!!!!!!!!!!!!!!!!!!!!!!
!
!               ttttttttttttttttt Cstv_ztop_8=Ver_z_8%m(0)=Ver_z_8%t(0)=Ver_z_8%x(0)
!
!               - - - - - - - - - Ver_z_8%m(1)
!
!               ================= Ver_z_8%t(1)=Ver_z_8%x(1) = ( Ver_z_8%m(2) + Ver_z_8%m(1) ) / 2
!
!               - - - - - - - - - Ver_z_8%m(2)
!
!                      ...
!
!               - - - - - - - - - Ver_z_8%m(k)
!
!               ================= Ver_x_8%t(k)=Ver_z_8%x(k) = ( Ver_z_8%m(k+1) + Ver_z_8%m(k) )/2
!
!               - - - - - - - - - Ver_z_8%m(k+1)
!
!                      ...
!
!               - - - - - - - - - Ver_z_8%m(n)
!
!               ================= Ver_z_8%t(n)
!
!               sssssssssssssssss Cstv_zsrf_8=Ver_z_8%m(n+1)=Ver_z_8%x(n)
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! n staggered MOMENTUM & THERMO **LAYERS** !!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!               Cstv_ztop_8 ttttttttttttttttttttttttttttttttttttttttttt Ver_z_8%x(0)
!                                                     \
!               Ver_z_8%m(1)- - - - - - - - - - - - - - Ver_dz_8%m(1)=Ver_z_8%x(1)-Cstv_ztop_8
!                                          /          /
!   Ver_z_8%m(2)-Ver_z_8%m(1)=Ver_dz_8%t(1)    ======================== Ver_z_8%x(1)
!                                          \          \
!               Ver_z_8%m(2)- - - - - - - - - - - - - - Ver_dz_8%m(2)=Ver_z_8%x(2)-Ver_z_8%x(1)
!                                                     /
!                           =========================================== Ver_z_8%x(2)
!
!                                                     ...
!
!                           =========================================== Ver_z_8%x(k-1)
!                                                     \
!               Ver_z_8%m(k)- - - - - - - - - - - - - - Ver_dz_8%m(k)=Ver_z_8%x(k)-Ver_z_8%x(k-1)
!                                          /          /
! Ver_z_8%m(k+1)-Ver_z_8%m(k)=Ver_dz_8%t(k)    ======================== Ver_z_8%x(k)
!                                          \
!             Ver_z_8%m(k+1)- - - - - - - - - - - - - - - - - - - - - -
!
!                                                     ...
!
!                           =========================================== Ver_z_8%x(n-1)
!                                                   \
!                                                    \
!                                                     \
!                           - - - - - - - - - - - - - - Ver_dz_8%m(n)=Ver_z_8%x(n)-Ver_z_8%x(n-1)
!                                          /          /
!    Cstv_zsrf_8-Ver_z_8%m(n)=Ver_dz_8%t(n)          /
!                                          \        /
!               Cstv_zsrf_8 sssssssssssssssssssssssssssssssssssssssssss Ver_z_8%x(n)
!
! arguments
! none
!

      character(len=32), parameter  :: VGRID_M_S  = 'ref-m'
      character(len=32), parameter  :: VGRID_T_S  = 'ref-t'

      type(vgrid_descriptor) :: vcoord
      integer k,istat,options_readwrite,options_readonly
      integer, dimension(:), pointer :: wkpti
      real, dimension(:), pointer :: std_p_prof=>null(),wkpt
      real*8  wk_8
      real*8, parameter :: zero=0.d0, one=1.d0, half=0.5d0
      real*8, dimension(:), pointer :: wkpt8
      character(len=32) :: REFP0_S, REFP0_LS_S
!     __________________________________________________________________
!
      allocate(   Ver_hyb%m(G_nk+1),       Ver_hyb%t(G_nk+1), &
           Ver_std_p_prof%m(G_nk+1),Ver_std_p_prof%t(G_nk+1), &
                  Ver_ip1%m(G_nk+1),       Ver_ip1%t(G_nk+1), &
                  Ver_a_8%m(G_nk+1),       Ver_a_8%t(G_nk+1), &
                  Ver_b_8%m(G_nk+1),       Ver_b_8%t(G_nk+1), &
                  Ver_c_8%m(G_nk+1),       Ver_c_8%t(G_nk+1), &
                Ver_z_8%m(0:G_nk+1),     Ver_z_8%t(0:G_nk  ), &
                                         Ver_z_8%x(0:G_nk  ), &
                 Ver_dz_8%m(G_nk  ),      Ver_dz_8%t(G_nk  ), &
                Ver_idz_8%m(G_nk  ),     Ver_idz_8%t(G_nk  ), &
               Ver_dbdz_8%m(G_nk  ),    Ver_dbdz_8%t(G_nk  ), &
               Ver_dcdz_8%m(G_nk  ),    Ver_dcdz_8%t(G_nk  ), &
                 Ver_wp_8%m(G_nk  ),      Ver_wp_8%t(G_nk  ), &
                 Ver_wm_8%m(G_nk  ),      Ver_wm_8%t(G_nk  ), &
              Ver_Tstar_8%m(G_nk+1),   Ver_Tstar_8%t(G_nk  ), &
                                          Ver_gama_8(G_nk  ), &
                 Ver_epsi_8(G_nk  ),     Ver_FIstr_8(G_nk+1), &
                  Ver_bzz_8(G_nk  ),     Ver_onezero(G_nk+1), &
               Ver_wpstar_8(G_nk  ),    Ver_wmstar_8(G_nk  ), &
                  Ver_wpA_8(G_nk  ),       Ver_wmA_8(G_nk  ), &
                  Ver_wpM_8(G_nk  ),       Ver_wmM_8(G_nk  ), &
                  Ver_wpC_8(G_nk  ),       Ver_wmC_8(G_nk  ), &
                 Ver_czz_8(G_nk  ))

      Cstv_pref_8 = 100000.d0
      Ver_code    = 6

      ! Construct vertical coordinate
      schm_sleve_L = .true.
      if( Hyb_rcoef(3) == -1 .or. Hyb_rcoef(4) == -1 )then
         if( Hyb_rcoef(3) /= -1 .or. Hyb_rcoef(4) /= -1 )then
            call handle_error_l(.true.,'set_zeta','Incorrect rcoef 3 and/or 4')
         endif
          Schm_sleve_L = .false.
      endif
      if(schm_sleve_L)then
         istat = vgd_new ( vcoord, kind=5, version=100, hyb=F_hybuser    , &
              rcoef1=Hyb_rcoef(1), rcoef2=Hyb_rcoef(2)                , &
              rcoef3=Hyb_rcoef(3), rcoef4=Hyb_rcoef(4)                , &
              ptop_out_8=wk_8, pref_8=Cstv_pref_8,  &
              dhm=0., dht=0. , avg_L=Schm_bcavg_L)
      else
         istat = vgd_new ( vcoord, kind=5, version=5, hyb=F_hybuser    , &
              rcoef1=Hyb_rcoef(1), rcoef2=Hyb_rcoef(2)                , &
              ptop_out_8=wk_8, pref_8=Cstv_pref_8,  &
              dhm=0., dht=0. , avg_L=Schm_bcavg_L)
      endif

      Cstv_ptop_8=wk_8

      if (Lun_debug_L) istat = vgd_print(vcoord)

      Cstv_Zsrf_8 = log(Cstv_pSref_8)
      Cstv_Ztop_8 = log(Cstv_ptop_8)
      Cstv_Sstar_8 = log(Cstv_pref_8/Cstv_pSref_8)
      Ver_zmin_8 = Cstv_Ztop_8
      Ver_zmax_8 = Cstv_Zsrf_8

      call handle_error_l(istat==VGD_OK,'set_zeta','coordinate construction failed')

      ! Retrieve information required to fill model arrays
      nullify(wkpt,wkpti,wkpt8)
      if (vgd_get(vcoord,'CA_M - vertical A coefficient (m)',wkpt8) /= VGD_OK) istat = VGD_ERROR
      Ver_a_8%m = wkpt8(1:size(Ver_a_8%m)); deallocate(wkpt8); nullify(wkpt8)
      if (vgd_get(vcoord,'CB_M - vertical B coefficient (m)',wkpt8) /= VGD_OK) istat = VGD_ERROR
      Ver_b_8%m = wkpt8(1:size(Ver_b_8%m)); deallocate(wkpt8); nullify(wkpt8)
      if (vgd_get(vcoord,'CA_T - vertical A coefficient (t)',wkpt8) /= VGD_OK) istat = VGD_ERROR
      Ver_a_8%t = wkpt8(1:size(Ver_a_8%t)); deallocate(wkpt8); nullify(wkpt8)
      if (vgd_get(vcoord,'CB_T - vertical B coefficient (t)',wkpt8) /= VGD_OK) istat = VGD_ERROR
      Ver_b_8%t = wkpt8(1:size(Ver_b_8%t)); deallocate(wkpt8); nullify(wkpt8)
      if(Schm_sleve_L)then
         if (vgd_get(vcoord,'CC_M - vertical C coefficient (m)',wkpt8) /= VGD_OK) istat = VGD_ERROR
         Ver_c_8%m = wkpt8(1:size(Ver_c_8%m)); deallocate(wkpt8); nullify(wkpt8)
         if (vgd_get(vcoord,'CC_T - vertical C coefficient (t)',wkpt8) /= VGD_OK) istat = VGD_ERROR
         Ver_c_8%t = wkpt8(1:size(Ver_c_8%t)); deallocate(wkpt8); nullify(wkpt8)
      else
         Ver_c_8%m=zero; Ver_c_8%t=zero
      endif
      if (vgd_get(vcoord,'VCDM - vertical coordinate (m)'   ,wkpt) /= VGD_OK) istat = VGD_ERROR
      Ver_hyb%m = wkpt(1:size(Ver_hyb%m)); deallocate(wkpt); nullify(wkpt)
      if (vgd_get(vcoord,'VCDT - vertical coordinate (t)'   ,wkpt) /= VGD_OK) istat = VGD_ERROR
      Ver_hyb%t = wkpt(1:size(Ver_hyb%t)); deallocate(wkpt); nullify(wkpt)
      if (vgd_get(vcoord,'VIPM - level ip1 list (m)'        ,wkpti) /= VGD_OK) istat = VGD_ERROR
      Ver_ip1%m = wkpti(1:size(Ver_ip1%m)); deallocate(wkpti); nullify(wkpti)
      if (vgd_get(vcoord,'VIPT - level ip1 list (t)'        ,wkpti) /= VGD_OK) istat = VGD_ERROR
      Ver_ip1%t = wkpti(1:size(Ver_ip1%t)); deallocate(wkpti); nullify(wkpti)

      call handle_error_l(istat==VGD_OK,'set_zeta','retrieving coordinate info')
      if(Schm_sleve_l)then
         istat = vgd_levels(vcoord,Ver_ip1%m,std_p_prof,sfc_field=100000.,sfc_field_ls=100000.,in_log=.false.)
         call handle_error_l(istat==VGD_OK,'set_zeta','problem getting standard pressure profile for m levels')
         Ver_std_p_prof%m=std_p_prof
         deallocate(std_p_prof)
         istat = vgd_levels(vcoord,Ver_ip1%t,std_p_prof,sfc_field=100000.,sfc_field_ls=100000.,in_log=.false.)
         call handle_error_l(istat==VGD_OK,'set_zeta','problem getting standard pressure profile for t levels')
         Ver_std_p_prof%t=std_p_prof
      else
         istat = vgd_levels(vcoord,Ver_ip1%m,std_p_prof,100000.,in_log=.false.)
         call handle_error_l(istat==VGD_OK,'set_zeta','problem getting standard pressure profile for m levels')
         Ver_std_p_prof%m=std_p_prof
         deallocate(std_p_prof)
         istat = vgd_levels(vcoord,Ver_ip1%t,std_p_prof,100000.,in_log=.false.)
         call handle_error_l(istat==VGD_OK,'set_zeta','problem getting standard pressure profile for t levels')
         Ver_std_p_prof%t=std_p_prof
      endif

!     -------------------------------
!     Define zeta(m/t/x) from A(m/t):
!     -------------------------------

!     Ver_a_8%m(1:G_nk+1):       G_nk   momentum levels + surface
!     Ver_a_8%t(1:G_nk+1):       G_nk   thermo   levels + surface
!     Ver_z_8%m(0:G_nk+1): top + G_nk   momentum levels + surface
!     Ver_z_8%t(0:G_nk  ): top + G_nk   thermo   levels
!     Ver_z_8%x(0:G_nk  ): top + G_nk-1 thermo   levels + surface

      Ver_z_8%m(0) = Cstv_Ztop_8
      do k = 1, G_nk+1
         Ver_z_8%m(k) = Ver_a_8%m(k)-Ver_b_8%m(k)*Cstv_Sstar_8
      enddo

     !Define the positions of true thermo levels
      Ver_z_8%t(0) = Cstv_Ztop_8
      do k = 1, G_nk
         Ver_z_8%t(k) = Ver_a_8%t(k)-Ver_b_8%t(k)*Cstv_Sstar_8
      enddo
      if( .not.Schm_trapeze_L .or. Schm_autobar_L ) Ver_z_8%t(G_nk)=Cstv_Zsrf_8

     !Define the positions of zeta_dot
      Ver_z_8%x(0) = Cstv_Ztop_8
      do k = 1, G_nk-1
         Ver_z_8%x(k) = Ver_a_8%t(k)-Ver_b_8%t(k)*Cstv_Sstar_8
      enddo
      Ver_z_8%x(G_nk)=Cstv_Zsrf_8

!     ----------------------
!     Compute dZ, 1/dZ, dBdZ
!     ----------------------

      do k=1,G_nk
           Ver_dz_8%m(k) = Ver_z_8%x(k) - Ver_z_8%x(k-1)
      enddo

      do k=1,G_nk
          Ver_idz_8%m(k) = one/Ver_dz_8%m(k)
           Ver_dz_8%t(k) = Ver_z_8%m(k+1) - Ver_z_8%m(k)
          Ver_idz_8%t(k) = one/Ver_dz_8%t(k)
         Ver_dbdz_8%t(k) = (Ver_b_8%m(k+1)-Ver_b_8%m(k))*Ver_idz_8%t(k)
         Ver_dcdz_8%t(k) = (Ver_c_8%m(k+1)-Ver_c_8%m(k))*Ver_idz_8%t(k)
         if(k == 1) then
            Ver_dbdz_8%m(k) = (Ver_b_8%t(k)-zero)*Ver_idz_8%m(k)
            Ver_dcdz_8%m(k) = (Ver_c_8%t(k)-zero)*Ver_idz_8%m(k)
         elseif(k == G_nk) then
            Ver_dbdz_8%m(k) = (one -Ver_b_8%t(k-1))*Ver_idz_8%m(k)
            Ver_dcdz_8%m(k) = (zero-Ver_c_8%t(k-1))*Ver_idz_8%m(k)
         else
            Ver_dbdz_8%m(k) = (Ver_b_8%t(k)-Ver_b_8%t(k-1))*Ver_idz_8%m(k)
            Ver_dcdz_8%m(k) = (Ver_c_8%t(k)-Ver_c_8%t(k-1))*Ver_idz_8%m(k)
         endif
      enddo

      if(Schm_autobar_L) then
         do k=1,G_nk
            Ver_dbdz_8%m(k) = one/(Cstv_Zsrf_8-Ver_z_8%m(1))
            Ver_dcdz_8%m(k) = one/(Cstv_Zsrf_8-Ver_z_8%m(1))
         enddo
      endif

!     -------------------------------------------------------
!     Compute AVERGING WEIGHTS FROM THERMO TO MOMENTUM LEVELS
!     -------------------------------------------------------

      do k=1,G_nk
         Ver_wmM_8(k) = (Ver_z_8%t(k)-Ver_z_8%m(k))/(Ver_z_8%t(k)-Ver_z_8%t(k-1))
         Ver_wpM_8(k) = one - Ver_wmM_8(k)

         Ver_wmC_8(k) = (Ver_z_8%x(k)-Ver_z_8%m(k))/(Ver_z_8%x(k)-Ver_z_8%x(k-1))
         Ver_wpC_8(k) = one - Ver_wmC_8(k)

         Ver_wmA_8(k) = (Ver_z_8%m(k)-Ver_z_8%x(k-1))/(Ver_z_8%x(k)-Ver_z_8%x(k-1))
         Ver_wpA_8(k) = one - Ver_wmA_8(k)

         Ver_wm_8%m(k)= Ver_wmA_8(k)
         Ver_wp_8%m(k)= Ver_wpA_8(k)
      enddo

!     -------------------------------------------------------
!     Compute AVERGING WEIGHTS FROM MOMENTUM TO THERMO LEVELS
!     -------------------------------------------------------

      do k=1,G_nk-1
         Ver_wp_8%t(k) = half
         Ver_wm_8%t(k) = half
      enddo
      Ver_wp_8%t(G_nk) = one
      Ver_wm_8%t(G_nk) = zero
!
!     SPECIAL WEIGHTS due to last thermo level
!
      Ver_wmstar_8 = zero
      Ver_wpstar_8 = one

      if(Schm_trapeze_L .and. .not.Schm_autobar_L ) then
         Ver_wmstar_8(G_nk)=half*Ver_dz_8%t(G_nk)/Ver_dz_8%m(G_nk)
         Ver_wpstar_8(G_nk)=one-Ver_wmstar_8(G_nk)
         Ver_wp_8%m(G_nk) = Ver_wpstar_8(G_nk) * Ver_wpA_8(G_nk)
         Ver_wm_8%m(G_nk) = one - Ver_wp_8%m(G_nk)
      endif
!     -------------------------------------------------------
!     Compute VERTICAL AVERGING OF B from thermo to momentum
!     -------------------------------------------------------

      do k=1,G_nk
         if(k == 1) then
            Ver_bzz_8(k) = Ver_wp_8%m(k)*Ver_b_8%t(k) &
                         + Ver_wm_8%m(k)*zero
            Ver_czz_8(k) = Ver_wp_8%m(k)*Ver_c_8%t(k) &
                         + Ver_wm_8%m(k)*zero
         elseif(k == G_nk) then
            Ver_bzz_8(k) = Ver_wp_8%m(k)*one &
                         + Ver_wm_8%m(k)*Ver_b_8%t(k-1)
            Ver_czz_8(k) = Ver_wp_8%m(k)*zero &
                         + Ver_wm_8%m(k)*Ver_c_8%t(k-1)
         else
            Ver_bzz_8(k) = Ver_wp_8%m(k)*Ver_b_8%t(k) &
                         + Ver_wm_8%m(k)*Ver_b_8%t(k-1)
            Ver_czz_8(k) = Ver_wp_8%m(k)*Ver_c_8%t(k) &
                         + Ver_wm_8%m(k)*Ver_c_8%t(k-1)
         endif
      enddo

!     -------------------------------------------------------
!     Initialize Ver_onezero
!     -------------------------------------------------------

      Ver_onezero=1.
      Ver_onezero(1)=0.

!     ----------------------------------------------------------
!     Save vcoord and ip1m/t for output
!     ----------------------------------------------------------
      REFP0_S = 'PW_P0:P'  !# gmmk_pw_p0_plus_s !NOTE: could gmmk_* be defined as parameters in a .cdk, this way it could be used here and would be more consistent
      REFP0_LS_S = ' '
      if (Schm_sleve_L) REFP0_LS_S = 'PW_P0_LS'  !# gmmk_pw_p0_ls_s
      istat = vgrid_wb_put(VGRID_M_S, vcoord, Ver_ip1%m,  &
           REFP0_S, REFP0_LS_S, F_overwrite_L=.true.)
      istat = vgrid_wb_put(VGRID_T_S, vcoord, Ver_ip1%t,  &
           REFP0_S, REFP0_LS_S, F_overwrite_L=.true.)
      istat = vgd_free(vcoord)

      options_readwrite = WB_IS_LOCAL
      options_readonly = options_readwrite + WB_REWRITE_NONE

      istat= wb_put('model/Vgrid/size-hybm',size(Ver_hyb%m),options_readonly)
      istat= wb_put('model/Vgrid/size-hybt',size(Ver_hyb%t),options_readonly)
      istat= wb_put('model/Vgrid/hybm'     ,Ver_hyb%m      ,options_readwrite)
      istat= wb_put('model/Vgrid/hybt'     ,Ver_hyb%t      ,options_readwrite)
      istat= wb_put('model/Vgrid/am'       ,Ver_a_8%m      ,options_readonly)
      istat= wb_put('model/Vgrid/bm'       ,Ver_b_8%m      ,options_readonly)
      istat= wb_put('model/Vgrid/at'       ,Ver_a_8%t      ,options_readonly)
      istat= wb_put('model/Vgrid/bt'       ,Ver_b_8%t      ,options_readonly)
      istat= wb_put('model/Vgrid/rcoef'    ,Hyb_rcoef      ,options_readonly)
      istat= wb_put('model/Vgrid/ptop'     ,Cstv_ptop_8    ,options_readonly)
      istat= wb_put('model/Vgrid/pref'     ,Cstv_pref_8    ,options_readonly)
      istat= wb_put('model/Vgrid/vcode'    ,Ver_code       ,options_readonly)

      if (Lun_debug_L) call prgenab()
!     __________________________________________________________________
!
      return
      end
