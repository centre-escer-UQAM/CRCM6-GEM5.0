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

   subroutine itf_phy_geom4 (F_istat)
   use dcst
   use iso_c_binding
   use nest_blending, only: nest_blend
   use gmm_geof
   use gem_options
   use geomh
   use tdpack
      use glb_ld
      use cstv
      use gmm_itf_mod
   implicit none
#include <arch_specific.hf>

   integer F_istat

#include <rmnlib_basics.hf>
#include <msg.h>

   logical :: nest_it
   integer :: i,j,istat,flag_r_n
   real, pointer :: wrk1(:,:)
   real*8 :: deg2rad_8
   real :: w1(l_minx:l_maxx,l_miny:l_maxy,2),&
           w2(l_minx:l_maxx,l_miny:l_maxy,2)
   type(gmm_metadata) :: mymeta
!
!-------------------------------------------------------------------
!
   F_istat = RMN_OK

   mymeta = GMM_NULL_METADATA
   mymeta%l(1) = gmm_layout(1,l_ni,0,0,l_ni)
   mymeta%l(2) = gmm_layout(1,l_nj,0,0,l_nj)

   flag_r_n = GMM_FLAG_RSTR+GMM_FLAG_IZER

   deg2rad_8 = acos(-1.D0)/180.D0

   nullify(wrk1)
   istat = gmm_create('DLAT',wrk1,mymeta)
   if (RMN_IS_OK(istat)) then
      wrk1 = deg2rad_8*geomh_latrx
   else
      F_istat = RMN_ERR
      call msg(MSG_ERROR,'(itf_phy_geom) Problem creating DLAT')
   endif

   nullify(wrk1)
   istat = gmm_create('DLON',wrk1,mymeta)
   if (RMN_IS_OK(istat)) then
      where(geomh_lonrx >= 0)
         wrk1 = deg2rad_8*geomh_lonrx
      elsewhere
         wrk1 = deg2rad_8*(geomh_lonrx+360.)
      endwhere
   else
      F_istat = RMN_ERR
      call msg(MSG_ERROR,'(itf_phy_geom) Problem creating DLON')
   endif

   nullify(wrk1)
   istat = gmm_create('DXDY',wrk1,mymeta)
   if (RMN_IS_OK(istat)) then
      do j = 1,l_nj
         do i = 1,l_ni
             wrk1(i,j) = geomh_hx_8 * geomh_hy_8 * Dcst_rayt_8**2 *geomh_cy_8(j)
         end do
      end do
   else
      F_istat = RMN_ERR
      call msg(MSG_ERROR,'(itf_phy_geom) Problem creating DXDY')
   endif

   nullify(wrk1)
   istat = gmm_create('TDMASK',wrk1,mymeta,flag_r_n)
   if (RMN_IS_OK(istat)) then
      w1 = 1.
      nest_it = ( Lam_0ptend_L .and. &
                ((Lam_blend_Hx > 0).or.(Lam_blend_Hy > 0)) )
      if ( nest_it ) then
         w2 = 0.
         call nest_blend (w1(:,:,1),w2(:,:,1),l_minx,l_maxx,l_miny,l_maxy,'M',level=G_nk+1)
      endif
      wrk1(1:l_ni,1:l_nj) = w1(1:l_ni,1:l_nj,1)
   else
      F_istat= RMN_ERR
      call msg(MSG_ERROR,'(itf_phy_geom) Problem creating TDMASK')
   endif
!
!-------------------------------------------------------------------
!
   return
   end subroutine itf_phy_geom4
