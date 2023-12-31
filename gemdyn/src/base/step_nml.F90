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

!**s/r step_nml - Read namelist time

      integer function step_nml (F_namelistf_S)
      use timestr_mod
      use step_options
      use grid_options
      use gem_options
      use lun
      use rstr
      use clib_itf_mod
      implicit none
#include <arch_specific.hf>

      character(len=*) F_namelistf_S

!authors    Michel Desgagne - Spring 2011
!
!revision
! v4_40 - Desgagne M.       - initial MPI version
! v5_00 - Winger K. (ESCER/UQAM) - write Fcst_rstrt_S as i8 
!
!object
!  Default configuration and reading namelist 'step'

#include <rmnlib_basics.hf>

      integer unf,err
      real*8 nesdt,nsteps
!
!-------------------------------------------------------------------
!
      step_nml = -1

      if ((F_namelistf_S == 'print').or.(F_namelistf_S == 'PRINT')) then
         step_nml = 0
         if (Lun_out > 0) then
            write (Lun_out  ,nml=step)
            write (Lun_out,8000) Step_alarm
         endif
         return
      endif

! Defaults values for ptopo namelist variables

      if (F_namelistf_S /= '') then
         unf = 0
         if (fnom (unf,F_namelistf_S, 'SEQ+OLD', 0) /= 0) then
            if (Lun_out >= 0) write (Lun_out, 7050) trim( F_namelistf_S )
            goto 9999
         endif
         rewind(unf)
         read (unf, nml=step, end = 9120, err=9130)
         goto 9000
      endif

 9120 if (Lun_out >= 0) write (Lun_out, 7060) trim( F_namelistf_S )
      goto 9999
 9130 if (Lun_out >= 0) write (Lun_out, 7070) trim( F_namelistf_S )
      goto 9999

 9000 if (Step_dt < 0.) then
         if (Lun_out > 0) write(Lun_out,*)  &
                    ' Step_dt must be specified in namelist &step'
         goto 9999
      endif

      err= 0

      if ( Fcst_start_S  == '' ) Fcst_start_S = '0H'
      if ( Fcst_end_S    == '' ) Fcst_end_S   = Fcst_start_S

      err= min( timestr2step (Step_initial, Fcst_start_S, Step_dt), err)
      err= min( timestr2step (Step_total  , Fcst_end_S  , Step_dt), err)
! transforming the Step_total into actual number of timesteps
      Step_total= Step_total - Step_initial

! Fcst_nesdt_S is transformed into a number of secondes (into Step_nesdt)
      nesdt= 1.d0
      err= min( timestr2step (nsteps, Fcst_nesdt_S, nesdt), err)
      Step_nesdt= dble(nsteps)

      if ( Fcst_rstrt_S  == '' ) then
         write(Fcst_rstrt_S,'(a,i8)') 'step,',Step_total+1
      else
         err= timestr_check ( Fcst_rstrt_S )
      endif

      Step_bkup_additional= Step_total+1
      err = clib_toupper ( Fcst_bkup_additional_S )
      if ( Fcst_bkup_additional_S /= 'NIL' ) then
         if (Fcst_bkup_additional_S == 'END' ) then
            Step_bkup_additional= Step_total
         else
            err= min( timestr2step (Step_bkup_additional, &
                      Fcst_bkup_additional_S, Step_dt), err)
         endif
      endif

      err = clib_toupper ( Fcst_bkup_S )
      if ( Fcst_bkup_S == 'END' ) Fcst_bkup_S= Fcst_end_S
      if ( Fcst_bkup_S /= 'NIL' ) then
         err = timestr_check (Fcst_bkup_S)
      endif

      if ( Fcst_gstat_S  == '' ) then
         Step_gstat= Step_total-Step_initial+1
      else
         err= min( timestr2step (Step_gstat, Fcst_gstat_S, Step_dt), err)
      endif
      if ( Fcst_spinphy_S  == '' ) then
         Step_spinphy= Step_total-Step_initial+1
      else
         err= min( timestr2step (Step_spinphy, Fcst_spinphy_S, Step_dt), err)
      endif

      if (err < 0) goto 9999

      Step_delay= Step_initial

      if (.not.Rstri_rstn_L) Lctl_step= Step_initial

      step_nml = 1

 7050 format (/,' FILE: ',A,' NOT AVAILABLE'/)
 7060 format (/,' Namelist &step NOT AVAILABLE in FILE: ',a/)
 7070 format (/,' NAMELIST &step IS INVALID IN FILE: ',a/)
 8000 format (/,' MODEL ALARM SET TO: ',i8,' secondes'/)

 9999 err = fclos (unf)
!
!-------------------------------------------------------------------
!
      return
      end function step_nml
