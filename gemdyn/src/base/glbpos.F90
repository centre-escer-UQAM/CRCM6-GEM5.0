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

!**s/r glbpos - sets up global indices for each PE
!

!
      subroutine glbpos
!
      use glb_ld
      use ptopo
      implicit none
#include <arch_specific.hf>
!
!author
!     M. Desgagne - V. Lee
!
!revision
! v2_00 - Desgagne/Lee       - initial MPI version
! v2_10 - Desgagne M.        - remove partitioning check
! v2_21 - Desgagne M.        - rpn_comm stooge for MPI
! v2_31 - Desgagne M.        - remove stkmemw
!
!object
!
!arguments
!     None
!

!
      integer dim, err, gindx(6,Ptopo_numproc)
!
!----------------------------------------------------------------------
!
      gindx = 0
!
      gindx(1,Ptopo_myproc+1) = l_i0
      gindx(2,Ptopo_myproc+1) = l_i0 + l_ni - 1
      gindx(3,Ptopo_myproc+1) = l_j0
      gindx(4,Ptopo_myproc+1) = l_j0 + l_nj - 1
      gindx(5,Ptopo_myproc+1) = 1
      gindx(6,Ptopo_myproc+1) = G_nk
!
      allocate (Ptopo_gindx(6,Ptopo_numproc))
      dim = 6*Ptopo_numproc
      call rpn_comm_ALLREDUCE (gindx,Ptopo_gindx,dim,"MPI_INTEGER", &
                                             "MPI_BOR","grid",err)
!
!----------------------------------------------------------------------
      return
      end

