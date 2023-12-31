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
!**s/r sol_fgmres - FGMRES based iterative elliptic solver
!
      subroutine sol_fgmres ( F_w2_8, F_w1_8, nil, njl,&
                              Minx, Maxx, Miny, Maxy, F_gnk,&
                              F_nk, print_conv_L, conv, its )
      use gem_options
      use glb_ld
      use lun
      use prec
      use sol
      implicit none
#include <arch_specific.hf>

      logical print_conv_L
      integer nil, njl, Minx, Maxx, Miny, Maxy, F_gnk, F_nk, its
      real*8 F_w1_8 (Minx:Maxx,Miny:Maxy,F_gnk), &
             F_w2_8 (Minx:Maxx,Miny:Maxy,F_gnk), conv

!author
!       Abdessamad Qaddouri -  2013
!
!revision
! v4_70 - Qaddouri A.       - initial version
!

      integer niloc,njloc,nloc
      integer halox, haloy, minx1, maxx1, minx2, maxx2
      real*8, dimension(:  ), allocatable :: wk11,wk22,rhs1,sol1
      real*8, dimension(:,:), allocatable :: vv_8,ww_8
      real*8, dimension(:), pointer :: sol2
      integer icode,ischmi
      integer :: lastdt = -1
      save sol2, lastdt
!
!     ---------------------------------------------------------------
!
      niloc = (Nil-pil_e)-(1+pil_w)+1
      njloc = (Njl-pil_n)-(1+pil_s)+1
      nloc  = niloc*njloc*F_nk

      if (lastdt == -1) then  ! allocate and initialize sol2
          allocate (sol2(nloc))
          sol2 = 0.
          lastdt = 0
      endif

      allocate (vv_8(nloc,sol_im+1),ww_8(nloc,sol_im),wk11(nloc), &
                wk22(nloc),rhs1(nloc),sol1(nloc))

! really necessary??? will kill the performance...
      vv_8=0.d0
      ww_8=0.d0
      wk11=0.d0
      wk22=0.d0
      rhs1=0.d0

      halox = 1
      haloy = halox

      minx1 = 1-halox
      maxx1 = nil+halox
      minx2 = 1-haloy
      maxx2 = njl+haloy

      sol_i0 = 1   + pil_w
      sol_in = nil - pil_e
      sol_j0 = 1   + pil_s
      sol_jn = njl - pil_n

      call tab_vec ( F_w1_8, Minx,Maxx,Miny,Maxy, F_nk, &
                     rhs1  , sol_i0,sol_in,sol_j0,sol_jn, +1 )

!     Initialize sol1 to last solution
      sol1  = sol2
      icode = 0
      conv  = 1.d0
      its   = 0
!
!-----------------------------------------------------------------------
!     F G M R E S   L O O P
!-----------------------------------------------------------------------
 1    continue
!-----------------------------------------------------------------------

      call fgmres2( nloc,sol_im,rhs1,sol1,ischmi,vv_8,ww_8,wk11,wk22, &
                           sol_fgm_eps,sol_fgm_maxits,its,conv,icode )

      sol2= sol1  ! Save solution sol1 into sol2

      if (icode == 1) then

         if (sol2D_precond_S == 'JACOBI')   then
               call pre_jacobi2D ( wk22,wk11,Prec_xevec_8,niloc,njloc,&
                                   F_nk,Prec_ai_8,Prec_bi_8,Prec_ci_8 )
         else
            call dcopy (nloc, wk11, 1, wk22, 1)
         endif

         goto 1

      else

         if (icode >= 2) then

            if (Lun_debug_L.and.print_conv_L) write(lun_out, 199) conv,its

            call sol_matvec ( wk22, wk11, Minx, Maxx, Miny, Maxy, &
                           nil,njl, F_nk, minx1,maxx1,minx2,maxx2 )

            goto 1

         endif

      endif

      call tab_vec ( F_w2_8 , Minx,Maxx,Miny,Maxy, F_nk, &
                     sol1   , sol_i0,sol_in,sol_j0,sol_jn, -1 )

 199  format (3x,'Iterative FGMRES solver convergence criteria: ',1pe14.7,' at iteration', i3)

      deallocate (wk11,wk22,rhs1,sol1,vv_8,ww_8)
!
!     ---------------------------------------------------------------
!
      return
      end
