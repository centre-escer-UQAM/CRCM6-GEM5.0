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
!
      subroutine adv_tracers (F_before_psadj_L)
      use adv_pos
      use grid_options
      use gem_options
      use glb_ld
      use tr3d
      use gmm_itf_mod
      use adv
      use tracers
      use gmm_tracers
      implicit none
#include <arch_specific.hf>

      logical, intent(in) :: F_before_psadj_L

   !@author  Michel Desgagne
   !@revisions
   ! v4_70 - Desgagne          - Initial Version
   ! v4_80 - Tanguay M.        - GEM4 Mass-Conservation
   !@objective Perform advection of all tracers


      logical qw_L,tr_not_before_psadj_L,tr_not_after_psadj_L
      integer  n,count,jext,erra,i,j,k,err
      integer nind, nind_s , num
      integer i0,j0,in,jn,k0        ! scope of advection operations
      integer i0_s,j0_s,in_s,jn_s   ! scope of advection operations on CORE
      integer, dimension(:),  allocatable :: ii , ii_s      ! pre-computed index used in the tricubic lagrangian interp loop
      type(gmm_metadata) :: mymeta
      real, pointer, dimension (:,:,:) :: fld_in, fld_out
!
!---------------------------------------------------------------------
!
      i0 = 1
      in = l_ni
      j0 = 1
      jn = l_nj

      jext=1
      if (Grd_yinyang_L) jext=2
      if (l_west)  i0 =        pil_w - jext  ;  i0_s = pil_w + 1
      if (l_south) j0 =        pil_s - jext  ;  j0_s = pil_s + 1
      if (l_east)  in = l_ni - pil_e + jext  ;  in_s = l_ni - pil_e
      if (l_north) jn = l_nj - pil_n + jext  ;  jn_s = l_nj - pil_n

      k0=Lam_gbpil_t+1

      num    =  l_ni*l_nj*l_nk
      nind   = (in-i0+1)*(jn-j0+1)*(l_nk-k0+1)
      nind_s = (in_s-i0_s+1)*(jn_s-j0_s+1)*(l_nk-k0+1)

      allocate (ii(4*nind), ii_s(4*nind_s))

      !Pre-compute indices ii used in: adv_tricub_lag3d_loop
      !-----------------------------------------------------
      call adv_get_indices(ii, pxt, pyt, pzt, num ,nind, &
                           i0, in, j0, jn, k0 , l_nk, 't')

      count=0

      Adv_component_S = 'INTP_TR'
      call timing_start2 (27, 'ADV_INTP_TR', 10)

      if (.NOT.F_before_psadj_L) Adv_done_precompute_L = .FALSE.

      do n=1,Tr3d_ntr

         qw_L= Tr3d_wload(n) .or. Tr3d_name_S(n)(1:2) == 'HU'

         if (qw_L .and. Schm_psadj==2 .and. (Tr3d_mono(n)>1 .or. Tr3d_mass(n)>0)) &
            call handle_error (-1,'ADV_TRACERS','CAUTION: Conservation TRACER (Water) + PSADJ DRY')

         if (F_before_psadj_L) then

            tr_not_before_psadj_L= .not.qw_L .or. (qw_L .and. Schm_psadj/=2)

            if (tr_not_before_psadj_L) cycle

         else

            tr_not_after_psadj_L= qw_L .and. Schm_psadj==2

            if (tr_not_after_psadj_L) cycle

         endif

         err= gmm_get('TR/'//trim(Tr3d_name_S(n))//':P' ,fld_in ,mymeta)
         err= gmm_get('TR/'//trim(Tr3d_name_S(n))//':M' ,fld_out,mymeta)

         if (Tr_extension_L) then

             !Pre-compute indices ii_s used in: adv_tricub_lag3d_loop when Bermejo-Conde or SLICE when LAM
             !--------------------------------------------------------------------------------------------
             if (count==0) call adv_get_indices(ii_s, pxt, pyt, pzt, num, nind_s, &
                                                i0_s, in_s, j0_s, jn_s, k0, l_nk, 't')

             call adv_cubic ('TR/'//trim(Tr3d_name_S(n))//':M', fld_out , fld_in, pxt, pyt, pzt, &
                             pxmu_s, pymu_s, pzmu_s, pxmv_s, pymv_s, pzmv_s, &
                             l_ni, l_nj, l_nk, l_minx, l_maxx, l_miny, l_maxy, &
                             nind_s, ii_s, i0_s, in_s, j0_s, jn_s, k0,'t', Tr3d_mono(n), Tr3d_mass(n) )

             count=1

         else

             call adv_cubic ('TR/'//trim(Tr3d_name_S(n))//':M', fld_out, fld_in, pxt, pyt, pzt, &
                              pxmu_s, pymu_s, pzmu_s, pxmv_s, pymv_s, pzmv_s, &
                              l_ni, l_nj, l_nk, l_minx, l_maxx, l_miny, l_maxy, &
                              nind, ii, i0, in, j0, jn, k0,'t', Tr3d_mono(n), Tr3d_mass(n) )
         endif

    end do

    if (.NOT.F_before_psadj_L) then 

       err = gmm_get (gmmk_pxto_s, pxto)
       err = gmm_get (gmmk_pyto_s, pyto)
       err = gmm_get (gmmk_pzto_s, pzto)

       do k=1,l_nk
       do j=1,l_nj
       do i=1,l_ni
          n = (k-1)*l_ni*l_nj + (j-1)*l_ni + i
          pxto(n) = pxt(i,j,k)
          pyto(n) = pyt(i,j,k)
          pzto(n) = pzt(i,j,k)
       end do
       end do
       end do

    end if

    call timing_stop (27)

    deallocate(ii,ii_s)
!
!---------------------------------------------------------------------
!
      return
      end subroutine adv_tracers
