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

!**s/r canonical_cases - Various actions related to canonical cases (Williamson/DCMIP)

      subroutine canonical_cases (F_action_S)

      use canonical
      use dcmip_options
      use wil_options
      use step_options
      use gmm_vt0
      use gmm_vt1
      use gem_options
      use tdpack, only : rgasd_8, cpd_8

      use glb_ld
      use cstv
      use lun
      use tr3d
      use ver
      use gmm_itf_mod
      use var_gmm
      implicit none

      !arguments
      character(len=*), intent(in) :: F_action_S

      !object
      !     F_action_S ='SET_GEOM': Print dcmip_HEIGHTS
      !     F_action_S ='SET_VT'  : Initialize gmm variables
      !     F_action_S ='BAC'     : Back subtitution (WINDS)
      !     F_action_S ='PHY'     : Physics
      !     F_action_S ='RAYLEIGH': Rayleigh Friction damping
      !     F_action_S ='ERR'     : Evaluate Errors/Diagnostics
      !     F_action_S ='OUT'     : Output dependent variables on standard file
      !     F_action_S ='XCHNG'   : Exchange staggered U-V
      !================================================================================

#include "gmm_gem_flags.hf"
#include <msg.h>

#define SET_GMMUSR_FLAG(MYMETA,MYFLAG) gmm_metadata(MYMETA%l,gmm_attributes(MYMETA%a%key,ior(MYMETA%a%uuid1,MYFLAG),MYMETA%a%uuid2,MYMETA%a%initmode,MYMETA%a%flags))

      !-------------------------------------------------------------------------------

      integer k,istat,pnip1,flag_r_n,i,j,n
      real, pointer, dimension(:,:,:) :: hu,cl,cl2,tr,tr_r,tr_e
      integer*8 :: flag_m_t,flag_m_u,flag_m_v,flag_s_f
      character*8 dumc
      real dcmip_height,dcmip_heightp1
      real bidon(l_minx:l_maxx,l_miny:l_maxy,G_nk),th(l_minx:l_maxx,l_miny:l_maxy,G_nk)
      type(gmm_metadata) :: mymeta3d_nk_u, mymeta3d_nk_v, mymeta3d_nk_t, mymeta2d_s
      real*8 :: pr_8
      logical Terminator_L
      real*4, parameter :: CLY_REF = 4.*10.**(-6)

      !-------------------------------------------------------------------------------

      if (.not. Schm_testcases_L) return

      Terminator_L = Dcmip_Terminator_L.or.Williamson_Terminator_L

      !-------------------
      !Print dcmip_HEIGHTS
      !-------------------
      if (F_action_S=="SET_GEOM") then

         if (Dcmip_case==0) return

         if (Lun_out > 0) then
            write (Lun_out,1005) G_nk,Hyb_rcoef

            do k=1,G_nk
               dcmip_height  =-8780.2*alog(Ver_hyb%m(k))
               if (k < G_nk)&
               dcmip_heightp1 =-8780.2*alog(Ver_hyb%m(k+1))
               if (k == G_nk) dcmip_heightp1 = 0.
               call convip(pnip1,Ver_hyb%m(k),5,1,dumc,.false.)
               write (Lun_out,1006) k,Ver_hyb%m(k),dcmip_height, &
                                    dcmip_height-dcmip_heightp1,pnip1
            end do

         endif

      !------------------------
      !Initialize gmm variables
      !------------------------
      elseif (F_action_S=="SET_VT") then

         flag_r_n = GMM_FLAG_RSTR+GMM_FLAG_IZER

         flag_m_t = FLAG_LVL_T                  !thermo
         flag_m_u = FLAG_LVL_M+GMM_FLAG_STAG_X  !momentum u
         flag_m_v = FLAG_LVL_M+GMM_FLAG_STAG_Y  !momentum v
         flag_s_f = 0                           !2d-surf

         mymeta3d_nk_u  = SET_GMMUSR_FLAG(meta3d_nk  ,flag_m_u)
         mymeta3d_nk_v  = SET_GMMUSR_FLAG(meta3d_nk  ,flag_m_v)
         mymeta3d_nk_t  = SET_GMMUSR_FLAG(meta3d_nk  ,flag_m_t)
         mymeta2d_s     = SET_GMMUSR_FLAG(meta2d     ,flag_s_f)

         gmmk_pth_s   = 'PTH'
         gmmk_thbase_s= 'THBA'
         gmmk_dtv_s   = 'DTV'

         gmmk_cly_s  = 'CLY'
         gmmk_acl_s  = 'ACL'
         gmmk_acl2_s = 'ACL2'
         gmmk_acly_s = 'ACLY'

         gmmk_irt_s  = 'IRT'
         gmmk_art_s  = 'ART'
         gmmk_wrt_s  = 'WRT'

         gmmk_uref_s = 'UREF'
         gmmk_vref_s = 'VREF'
         gmmk_wref_s = 'WREF'
         gmmk_zdref_s= 'ZDRF'
         gmmk_qvref_s= 'QVRF'
         gmmk_qcref_s= 'QCRF'
         gmmk_qrref_s= 'RWRF'
         gmmk_thref_s= 'THRF'
         gmmk_fcu_s  = 'FUC'
         gmmk_fcv_s  = 'FVC'
         gmmk_fcw_s  = 'FWC'

         gmmk_q1ref_s= 'QR1'
         gmmk_q2ref_s= 'QR2'
         gmmk_q3ref_s= 'QR3'
         gmmk_q4ref_s= 'QR4'

         gmmk_q1err_s= 'QE1'
         gmmk_q2err_s= 'QE2'
         gmmk_q3err_s= 'QE3'
         gmmk_q4err_s= 'QE4'

         gmmk_clyref_s= 'CLYR'
         gmmk_clyerr_s= 'CLYE'

         istat = GMM_OK

         istat = min(gmm_create(gmmk_pth_s,   pth,   mymeta3d_nk_t, flag_r_n),istat)
         istat = min(gmm_create(gmmk_thbase_s,thbase,mymeta3d_nk_t, flag_r_n),istat)
         istat = min(gmm_create(gmmk_dtv_s,   dtv,   mymeta3d_nk_t, flag_r_n),istat)

         istat = min(gmm_create(gmmk_cly_s,   cly,   mymeta3d_nk_t, flag_r_n),istat)
         istat = min(gmm_create(gmmk_acl_s,   acl,   mymeta2d_s   , flag_r_n),istat)
         istat = min(gmm_create(gmmk_acl2_s,  acl2,  mymeta2d_s   , flag_r_n),istat)
         istat = min(gmm_create(gmmk_acly_s,  acly,  mymeta2d_s   , flag_r_n),istat)

         istat = min(gmm_create(gmmk_irt_s,   irt,   mymeta2d_s   , flag_r_n),istat)
         istat = min(gmm_create(gmmk_art_s,   art,   mymeta2d_s   , flag_r_n),istat)
         istat = min(gmm_create(gmmk_wrt_s,   wrt,   mymeta2d_s   , flag_r_n),istat)

         istat = min(gmm_create(gmmk_uref_s,  uref,  mymeta3d_nk_u, flag_r_n),istat)
         istat = min(gmm_create(gmmk_vref_s,  vref,  mymeta3d_nk_v, flag_r_n),istat)
         istat = min(gmm_create(gmmk_wref_s,  wref,  mymeta3d_nk_t, flag_r_n),istat)
         istat = min(gmm_create(gmmk_zdref_s, zdref, mymeta3d_nk_t, flag_r_n),istat)
         istat = min(gmm_create(gmmk_qvref_s, qvref, mymeta3d_nk_t, flag_r_n),istat)
         istat = min(gmm_create(gmmk_qcref_s, qcref, mymeta3d_nk_t, flag_r_n),istat)
         istat = min(gmm_create(gmmk_qrref_s, qrref, mymeta3d_nk_t, flag_r_n),istat)
         istat = min(gmm_create(gmmk_thref_s, thref, mymeta3d_nk_t, flag_r_n),istat)
         istat = min(gmm_create(gmmk_fcu_s,   fcu,   mymeta3d_nk_u, flag_r_n),istat)
         istat = min(gmm_create(gmmk_fcv_s,   fcv,   mymeta3d_nk_v, flag_r_n),istat)
         istat = min(gmm_create(gmmk_fcw_s,   fcw,   mymeta3d_nk_t, flag_r_n),istat)

         istat = min(gmm_create(gmmk_q1ref_s, q1ref, mymeta3d_nk_t, flag_r_n),istat)
         istat = min(gmm_create(gmmk_q2ref_s, q2ref, mymeta3d_nk_t, flag_r_n),istat)
         istat = min(gmm_create(gmmk_q3ref_s, q3ref, mymeta3d_nk_t, flag_r_n),istat)
         istat = min(gmm_create(gmmk_q4ref_s, q4ref, mymeta3d_nk_t, flag_r_n),istat)

         istat = min(gmm_create(gmmk_q1err_s, q1err, mymeta3d_nk_t, flag_r_n),istat)
         istat = min(gmm_create(gmmk_q2err_s, q2err, mymeta3d_nk_t, flag_r_n),istat)
         istat = min(gmm_create(gmmk_q3err_s, q3err, mymeta3d_nk_t, flag_r_n),istat)
         istat = min(gmm_create(gmmk_q4err_s, q4err, mymeta3d_nk_t, flag_r_n),istat)

         istat = min(gmm_create(gmmk_clyref_s,clyref,mymeta3d_nk_t, flag_r_n),istat)
         istat = min(gmm_create(gmmk_clyerr_s,clyerr,mymeta3d_nk_t, flag_r_n),istat)

         if (GMM_IS_ERROR(istat)) &
             call msg(MSG_ERROR,'set_vt ERROR at gmm_create(V_TEST*)')

      !------------------------
      !Back subtitution (WINDS)
      !------------------------
      elseif (F_action_S=="BAC") then

         istat = gmm_get(gmmk_ut0_s, ut0)
         istat = gmm_get(gmmk_vt0_s, vt0)
         istat = gmm_get(gmmk_zdt0_s,zdt0)

         if (Williamson_case==1) then

            call wil_uvcase1 (ut0,vt0,l_minx,l_maxx,l_miny,l_maxy,G_nk,.true.,Lctl_step)

            return

         endif

         if (Dcmip_case>=11.and.Dcmip_case<=13) then

             if (Dcmip_case==11) call dcmip_tracers11_transport (ut0,vt0,zdt0,bidon,bidon,bidon,bidon,bidon,bidon,bidon,bidon, &
                                                                 l_minx,l_maxx,l_miny,l_maxy,G_nk)

             if (Dcmip_case==12) call dcmip_tracers12_transport (ut0,vt0,zdt0,bidon,bidon,bidon,bidon,bidon, &
                                                                 l_minx,l_maxx,l_miny,l_maxy,G_nk)

             if (Dcmip_case==13) call dcmip_tracers13_transport (ut0,vt0,zdt0,bidon,bidon,bidon,bidon,bidon,bidon,bidon,bidon, &
                                                                 l_minx,l_maxx,l_miny,l_maxy,G_nk)

             return

         endif

      !-------
      !Physics
      !-------
      elseif (F_action_S=="PHY") then

         if (Dcmip_prec_type/=-1.or.Dcmip_pbl_type/=-1) call dcmip_2016_physics()

         if (Terminator_L) call canonical_Terminator()

         if (Dcmip_vrd_L) call dcmip_vrd_main()

      !------------------------------------------------------------------------
      !Rayleigh Friction damping for Mountain waves over a Schaer-type mountain
      !------------------------------------------------------------------------
      elseif (F_action_S=="RAYLEIGH") then

         if (Dcmip_case>0.and.Dcmip_Rayleigh_friction_L) call dcmip_Rayleigh_friction

      !---------------------------
      !Evaluate Errors/Diagnostics
      !---------------------------
      elseif (F_action_S=="ERR") then

         if (Williamson_case==1) call wil_diagnostics (Lctl_step)

         if (Dcmip_case>0) call dcmip_diagnostics (Lctl_step)

      !-------------------------------------------
      !Output dependent variables on standard file
      !-------------------------------------------
      elseif (F_action_S=="OUT") then

         if (Williamson_case==1) then

            do n=1,Tr3d_ntr

               if (.NOT.((Tr3d_name_S(n)(1:2)=='Q1').or. &
                         (Tr3d_name_S(n)(1:2)=='Q2').or. &
                         (Tr3d_name_S(n)(1:2)=='Q3').or. &
                         (Tr3d_name_S(n)(1:2)=='Q4'))) cycle

               if (Lctl_step/=Step_total    .and.(Williamson_NAIR==1.or.Williamson_NAIR==2)) cycle

               if (Tr3d_name_S(n)(1:2)/='Q1'.and.(Williamson_NAIR==0.or.Williamson_NAIR==3)) cycle

               istat = gmm_get('TR/'//trim(Tr3d_name_S(n))//':P',tr)

               if (Tr3d_name_S(n)(1:2)=='Q1') istat = gmm_get(gmmk_q1ref_s,tr_r)
               if (Tr3d_name_S(n)(1:2)=='Q2') istat = gmm_get(gmmk_q2ref_s,tr_r)
               if (Tr3d_name_S(n)(1:2)=='Q3') istat = gmm_get(gmmk_q3ref_s,tr_r)
               if (Tr3d_name_S(n)(1:2)=='Q4') istat = gmm_get(gmmk_q4ref_s,tr_r)

               if (Tr3d_name_S(n)(1:2)=='Q1') istat = gmm_get(gmmk_q1err_s,tr_e)
               if (Tr3d_name_S(n)(1:2)=='Q2') istat = gmm_get(gmmk_q2err_s,tr_e)
               if (Tr3d_name_S(n)(1:2)=='Q3') istat = gmm_get(gmmk_q3err_s,tr_e)
               if (Tr3d_name_S(n)(1:2)=='Q4') istat = gmm_get(gmmk_q4err_s,tr_e)

               !Initialize REFERENCE at TIME>0
               !------------------------------
               if (Williamson_Nair==0) call wil_case1(tr_r,l_minx,l_maxx,l_miny,l_maxy,G_nk,0,Lctl_step)
               if (Williamson_Nair==3) call wil_case1(tr_r,l_minx,l_maxx,l_miny,l_maxy,G_nk,5,Lctl_step)

               !Initialize ERROR
               !----------------
               tr_e(1:l_ni,1:l_nj,1:G_nk) = tr(1:l_ni,1:l_nj,1:G_nk) - tr_r(1:l_ni,1:l_nj,1:G_nk)

            enddo

         endif

         if (Terminator_L) then

            istat = gmm_get ('TR/CL:P' , cl )
            istat = gmm_get ('TR/CL2:P', cl2)

            !Initialize CLY
            !--------------
            istat = gmm_get (gmmk_cly_s, cly)

            cly(1:l_ni,1:l_nj,1:G_nk) = cl(1:l_ni,1:l_nj,1:G_nk) + 2.0d0 * cl2(1:l_ni,1:l_nj,1:G_nk)

            !Initialize CLY REFERENCE
            !------------------------
            istat = gmm_get (gmmk_clyref_s, clyref)

            clyref(1:l_ni,1:l_nj,1:G_nk) =  CLY_REF

            !Initialize CLY ERROR
            !--------------------
            istat = gmm_get (gmmk_clyerr_s, clyerr)

            clyerr(1:l_ni,1:l_nj,1:G_nk) = cly(1:l_ni,1:l_nj,1:G_nk) - clyref(1:l_ni,1:l_nj,1:G_nk)

            !Initialize Average Column Integrated of CL/CL2/CLY
            !--------------------------------------------------
            istat = gmm_get (gmmk_acl_s,  acl )
            istat = gmm_get (gmmk_acl2_s, acl2)
            istat = gmm_get (gmmk_acly_s, acly)

            call dcmip_avg_column_integrated (acl ,cl, l_minx,l_maxx,l_miny,l_maxy,G_nk)
            call dcmip_avg_column_integrated (acl2,cl2,l_minx,l_maxx,l_miny,l_maxy,G_nk)
            call dcmip_avg_column_integrated (acly,cly,l_minx,l_maxx,l_miny,l_maxy,G_nk)

         endif

         if (Dcmip_case>=11.and.Dcmip_case<=13) then

            do n=1,Tr3d_ntr

               if (.NOT.((Tr3d_name_S(n)(1:2)=='Q1').or. &
                         (Tr3d_name_S(n)(1:2)=='Q2').or. &
                         (Tr3d_name_S(n)(1:2)=='Q3').or. &
                         (Tr3d_name_S(n)(1:2)=='Q4'))) cycle

               if (Tr3d_name_S(n)(1:2)/='Q1'.and.Dcmip_case==12) cycle

               if (Lctl_step/=Step_total) cycle

               istat = gmm_get('TR/'//trim(Tr3d_name_S(n))//':P',tr)

               if (Tr3d_name_S(n)(1:2)=='Q1') istat = gmm_get(gmmk_q1ref_s,tr_r)
               if (Tr3d_name_S(n)(1:2)=='Q2') istat = gmm_get(gmmk_q2ref_s,tr_r)
               if (Tr3d_name_S(n)(1:2)=='Q3') istat = gmm_get(gmmk_q3ref_s,tr_r)
               if (Tr3d_name_S(n)(1:2)=='Q4') istat = gmm_get(gmmk_q4ref_s,tr_r)

               if (Tr3d_name_S(n)(1:2)=='Q1') istat = gmm_get(gmmk_q1err_s,tr_e)
               if (Tr3d_name_S(n)(1:2)=='Q2') istat = gmm_get(gmmk_q2err_s,tr_e)
               if (Tr3d_name_S(n)(1:2)=='Q3') istat = gmm_get(gmmk_q3err_s,tr_e)
               if (Tr3d_name_S(n)(1:2)=='Q4') istat = gmm_get(gmmk_q4err_s,tr_e)

               !Initialize ERROR
               !----------------
               tr_e(1:l_ni,1:l_nj,1:G_nk) = tr(1:l_ni,1:l_nj,1:G_nk) - tr_r(1:l_ni,1:l_nj,1:G_nk)

            enddo

         endif

         !Prepare Perturbation of real potential Temperature
         !--------------------------------------------------
         if (Dcmip_case==31.or.Dcmip_case==163) then

            istat = gmm_get(gmmk_pth_s,   pth)
            istat = gmm_get(gmmk_thbase_s,thbase)
            istat = gmm_get(gmmk_tt1_s,   tt1)
            istat = gmm_get(gmmk_st1_s,   st1)
            istat = gmm_get('TR/'//'HU'//':P',hu)

            do k=1,G_nk
               do j=1,l_nj
                  do i=1,l_ni

                     !Real potential temperature
                     !--------------------------
                     pr_8       = exp(Ver_a_8%t(k) + Ver_b_8%t(k)*st1(i,j))
                     th (i,j,k) = (tt1(i,j,k) / (1.d0 + 0.608d0 * hu(i,j,k))) * (Cstv_pref_8/pr_8) ** (rgasd_8/cpd_8)

                     pth(i,j,k) = th(i,j,k) - thbase(i,j,k)

                  enddo
               enddo
            enddo

         endif

         !Prepare Perturbation of Virtual Temperature
         !-------------------------------------------
         istat = gmm_get(gmmk_dtv_s, dtv)
         istat = gmm_get(gmmk_tt1_s, tt1)

         dtv(1:l_ni,1:l_nj,1:G_nk) = tt1(1:l_ni,1:l_nj,1:G_nk) - Cstv_Tstr_8

      !----------------------
      !Exchange staggered U-V
      !----------------------
      elseif (F_action_S=="XCHNG") then

         istat = gmm_get(gmmk_ut1_s , ut1)
         istat = gmm_get(gmmk_vt1_s , vt1)

         call yyg_nestuv (ut1,vt1,l_minx,l_maxx,l_miny,l_maxy,G_nk)

         call pw_update_UV

      else

         call handle_error(-1,'CANONICAL_CASES','F_action_S unknown')

      endif

      !-------------------------------------------------------------------------------

      return

 1005 format (/'STAGGERED VERTICAL LAYERING ON',I4,' MOMENTUM HYBRID LEVELS WITH ', &
               'Hyb_rcoef= ',2f7.2,':'/ &
               2x,'level',10x,'HYB',2x,'~dcmip_HEIGHTS',2x,'~dcmip_DELTA_Z',2x,'IP1')
 1006 format (1x,i4,3x,es16.4,2(6x,f6.0),4x,i10)

      end
