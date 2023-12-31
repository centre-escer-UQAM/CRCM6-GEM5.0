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

!**s/r dynkernel_nml - Read namelist dyn_kernel

      integer function dynkernel_nml (F_namelistf_S)
      use dynkernel_options
      use lun
      use clib_itf_mod
      implicit none
#include <arch_specific.hf>

      character(len=*), intent(in):: F_namelistf_S

#include <rmnlib_basics.hf>

      integer unf,err
!
!-------------------------------------------------------------------
!
      dynkernel_nml = -1

      if ((F_namelistf_S == 'print').or.(F_namelistf_S == 'PRINT')) then
         dynkernel_nml = 0
         if (Lun_out > 0) then
            write (Lun_out  ,nml=dyn_kernel)
         endif
         return
      endif

      if (F_namelistf_S /= '') then
         unf = 0
         if (fnom (unf,F_namelistf_S, 'SEQ+OLD', 0) /= 0) then
            if (Lun_out >= 0) write (Lun_out, 7050) trim( F_namelistf_S )
            goto 9999
         endif
         rewind(unf)
         read (unf, nml=dyn_kernel, iostat=err)
         err = -1 * err
         if (err < 0) goto 9130
         goto 9000
      endif

 9130 if (Lun_out >= 0) write (Lun_out, 7070) trim( F_namelistf_S )
      goto 9999

 9000 err = clib_toupper ( Dynamics_Kernel_S )

      dynkernel_nml = 1

 9999 err = fclos (unf)

 7050 format (/,' FILE: ',A,' NOT AVAILABLE'/)
 7070 format (/,' NAMELIST &dynKernel_S IS INVALID IN FILE: ',a/)
!
!-------------------------------------------------------------------
!
      return
      end function dynkernel_nml
