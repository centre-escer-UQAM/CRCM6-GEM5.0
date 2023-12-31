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

!**s/r ptopo_nml - Read namelist ptopo
!
      integer function ptopo_nml (F_namelistf_S)
      use lun
      use ptopo
      implicit none
#include <arch_specific.hf>
!
      character(len=*) F_namelistf_S
!
!author
!     Michel Desgagne - Summer 2006
!
!revision
! v3_30 - Desgagne M.       - initial version
! v3_31 - Lee V.            - binding is restricted to when SMT is equal
!                             or twice Ptopo_npeOpenMP
!
!object
!  Default configuration and reading namelist ptopo
!

      integer,external :: fnom

      integer unf
!
!-------------------------------------------------------------------
!
      ptopo_nml = -1

      if ((F_namelistf_S == 'print').or.(F_namelistf_S == 'PRINT')) then
         ptopo_nml = 0
         if (Lun_out > 0) write (6  ,nml=resources)
         return
      endif

! Defaults values for ptopo namelist variables

      Ptopo_npex   =  1
      Ptopo_npey   =  1
      Ptopo_nthreads_dyn = 0
      Ptopo_nthreads_phy = 0
      Ptopo_bind_L = .false.

      unf=0
      if (fnom (unf, F_namelistf_S, 'SEQ+OLD' , 0) == 0) then
         rewind(unf)
         read (unf, nml=resources, end=7110, err=9110)
 7110    call fclos (unf)
         ptopo_nml = 1
         goto 7777
      else
         goto 9220
      endif

 9110 write (6, 8150) trim( F_namelistf_S )
      call fclos (unf)
      goto 7777

 9220 write (6, 8155) trim( F_namelistf_S )

 8150 format (/,' NAMELIST resources    INVALID IN FILE: ',a)
 8155 format (/,' FILE: ',a,' NOT AVAILABLE')
!
!-------------------------------------------------------------------
!
 7777 return
      end
