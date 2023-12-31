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

      logical function decomp3 ( F_npts, F_min, F_max, F_lni, F_npartiel, F_halo, F_start, &
                                 F_alongx_L, F_fill_L, F_npe, F_lowestsize, F_checkparti_L,&
                                 F_relax)
      implicit none
#include <arch_specific.hf>

      logical F_alongx_L, F_fill_L, F_checkparti_L
      integer F_npts, F_min, F_max, F_lni, F_npartiel, F_halo, F_start, F_npe
      integer F_lowestsize, F_relax

      logical  check_parti
      integer  RPN_COMM_limit_2,rpn_comm_topo_2
      external RPN_COMM_limit_2,rpn_comm_topo_2,check_parti

      integer istat,my_id
      integer count(F_npe),depl(F_npe)
!
!-------------------------------------------------------------------
!
      decomp3 = .false.

      if (F_checkparti_L ) then

         istat= RPN_COMM_limit_2 (0, F_npe, 1, F_npts, F_min, F_max, &
                                  count, depl, F_relax)
         if (istat >= 0) &
         decomp3 = ( minval(count) >= F_lowestsize )

      else

         istat= rpn_comm_topo_2( F_npts, F_min, F_max, F_lni, F_npartiel, &
                  F_halo, F_start, F_alongx_L, F_fill_L, F_relax, .false. )

         decomp3 = .not.(istat < 0)
         if (F_lowestsize > 0) decomp3= (decomp3) .and. (F_lni >= F_lowestsize)

      endif

      if (.not.decomp3) write (6,1001) F_npts,F_npe
 1001 format(/' DECOMP: illegal partitionning ====> ',i7,' / ',i7)
!
!-------------------------------------------------------------------
!
      return
      end

