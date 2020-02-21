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

      subroutine out_outdir
      use iso_c_binding
      use timestr_mod, only: timestr_prognum,timestr_unitfact,timestr_date_S
      use step_options
      use grid_options
      use gem_options
      use cstv
      use lun
      use out_mod
      use out3
      use path
      use clib_itf_mod
      use ptopo
      use out_listes        ! needed to only create directories when needed (KW)
      implicit none
#include <arch_specific.hf>

!AUTHOR   Michel Desgagne  - Summer 2015
!
!REVISION
! v4_80  - Desgagne M.      - Initial version
! v4_812 - Winger K. (ESCER/UQAM) - Create laststep directories as laststep_YYYYMMDD.hhmmss
!                                 - Make sure 00Z goes to the previous day

#include <rmnlib_basics.hf>
      include "rpn_comm.inc"

      character(len=1024),save :: dirstep_S=' ', diryy_S=' ', dirbloc_S=' ', &
                                  FMT=' ', last_S=' '
      character*15 postjob_S
      character*7  blocxy_S
      character*15 out_date_S                ! YYYYMMDD.hhmmss
      character*15, save :: last_out_date_S = "YYYYMMDD.hhmmss"
      integer      err, stepno, prognum, prognum1
      integer      date, year, month, time, stamp, stamp_last, t
      integer      pp_int_step, next_out_step
      integer, save :: last_out_step = -1
      real :: interval
      real*8       nhours_8
!
!----------------------------------------------------------------------
!

!print *,'In out_outdir'
!print *,'out_outdir Step_kount    = ',Step_kount
!print *,'out_outdir Step_total    = ',Step_total
!print *,'out_outdir Out_dateo     = ',Out_dateo
!print *,'out_outdir Out_endstepno = ',Out_endstepno

!print *,'Step_initial = ',Step_initial
!print *,'Step_delay   = ',Step_delay
!      upperlimit = Step_total + Step_initial

      Out_post_L = .false.

      call out_steps

      ! if there is no output for this time step => return (KW)
      if ( outd_sorties(0,Lctl_step) .eq. 0 .and.  &
           outp_sorties(0,Lctl_step) .eq. 0 .and.  &
           outc_sorties(0,Lctl_step) .eq. 0 ) return



      if ( Init_mode_L .and. (Step_kount.ge.Init_halfspan) ) return

      write (blocxy_S,'(I3.3,"-",I3.3)') Ptopo_mycol, Ptopo_myrow

!      if ( Out3_unit_S(1:3) .eq. 'MON' ) Out3_close_interval = 1.
      Out3_close_interval = 1   ! KW
      Out3_postproc_fact  = max(1, Out3_postproc_fact)

!print *,'out_outdir Out_endstepno:',Out_endstepno
      interval = Out3_close_interval * Out3_postproc_fact
      stepno = max(Step_kount,1)
      err = timestr_prognum(prognum ,Out3_unit_S,interval,Out_dateo,&
                            float(Out_deet),stepno  ,Out_endstepno)
      err = timestr_prognum(prognum1,Out3_unit_S,interval,Out_dateo,&
                            float(Out_deet),stepno+1,Out_endstepno)
      Out_post_L = (prognum1 > prognum .or. stepno == Out_endstepno)



      ! Always create monthly laststep directories as laststep_YYYYMMDD.hhmmss (KW)
      ! And make sure 00Z goes to the previous day (KW)

      ! Determine output directory name/date
      Out_post_L = .false.
      pp_int_step = 0
      ! Use Step_kount-1 to make sure 00Z goes to the previous day/month
      out_date_S = timestr_date_S(Out_dateo,sngl(Cstv_dt_8),Step_kount-1)

      ! Determine post processing interval in steps
      select case (Out3_unit_S(1:3))
         case ('MON') ! Monthly files
                      out_date_S( 7:15) = '01.000000'    ! First of the month at 00Z

         case ('DAY') ! Daily files
                      out_date_S(10:15) =    '000000'    ! 00Z of the day
                      pp_int_step =       int(interval * 86400. / Cstv_dt_8)

         case ('HOU') ! Hourly files
                      out_date_S(12:15) =      '0000'    ! 00m of the hour
                      pp_int_step = max(1,int(interval *  3600. / Cstv_dt_8))

         case ('MIN') ! Minutely files
                      out_date_S(14:15) =        '00'    ! 00s of the minute
                      pp_int_step = max(1,int(interval *    60. / Cstv_dt_8))

         case ('SEC') ! Secondly files
                      pp_int_step = max(1,int(interval          / Cstv_dt_8))

         case ('STE') ! Step files
                      pp_int_step = max(1,int(interval                     ))
      end select

!print *,'out_outdir Step_job, pp_int_step:',Step_job, pp_int_step
!print *,'out_outdir Step_kount, Step_job:',Step_kount, Step_job

!print *,'out_outdir Step_kount,mod(Step_kount,pp_int_step):',Step_kount,mod(Step_kount,pp_int_step)

      ! Determine if post processing should get started
      if ( Out3_unit_S(1:3) .ne. 'MON' ) then
        if ( mod(Step_job,pp_int_step) == 0 ) Out_post_L = .true.
      end if

      ! Always start post processing at the end of a job
      if ( stepno == Out_endstepno ) Out_post_L = .true.

!print *,'out_outdir next_out_step,Out_post_L:', next_out_step,Out_post_L

      postjob_S = out_date_S

      ! Always create new output directory at time steps 0 and 1
      if ( Step_kount .eq. 0 .or. Step_kount .eq. 1 ) last_out_date_S = "YYYYMMDD.hhmmss"

      ! If output directory has not been created yet => create it
      if ( Step_kount .eq. 0) next_out_step = 0

print *,'out_outdir t,last_out_date_S,out_date_S:',Step_kount,last_out_date_S,out_date_S

      if ( last_out_date_S .ne. out_date_S ) then

        last_out_date_S = out_date_S

        ! Give special name for time step 0
        if ( Step_kount.eq.0 ) then
           out_date_S = timestr_date_S(Out_dateo,sngl(Cstv_dt_8),Step_kount)
           postjob_S = '00000000.' // out_date_S(1:6)
           last_out_date_S = "YYYYMMDD.hhmmss"
        endif

        Out_laststep_S = 'laststep_'//postjob_S
        Out_dirname_S  = trim(Path_output_S)//'/'//trim(Out_laststep_S) !  Add 'trim' (KW)
!print *,'outdir Step_kount,Out_laststep_S: ',Step_kount,Out_laststep_S
!print *,'outdir Out_dirname_S :',Out_dirname_S
!print *,'outdir Out_post_L:',Out_post_L
        ! PE0 is responsible for creating shared subdir structure
        if (dirstep_S /= Out_dirname_S) then
           dirstep_S = Out_dirname_S        
           if (Ptopo_myproc == 0 .and. Ptopo_couleur == 0) then
              err = clib_mkdir(trim(Out_dirname_S))
              if (Lun_out>0) write(Lun_out,1001) trim(Out_laststep_S),Step_kount
           endif
        endif
      
        ! Wait for Grid PE0 to be finished subdir creation
        call rpn_comm_barrier (RPN_COMM_ALLGRIDS, err)

        ! Each io pe now creates an independent subdir for outputs
        Out_dirname_S = trim(Out_dirname_S)//'/'//blocxy_S
        err = CLIB_OK
        if (Out3_iome .ge. 0 .and. dirbloc_S /= Out_dirname_S &
                             .and. Ptopo_couleur == 0) then
           dirbloc_S = Out_dirname_S
           err = clib_mkdir ( trim(Out_dirname_S) )
           err = clib_isdir ( trim(Out_dirname_S) )
        endif

        call gem_error (err,'out_outdir','unable to create output directory structure')

      endif ! Create output directory
      
 1001 format (' OUT_OUTDIR: DIRECTORY output/',a,' was created at timestep: ',i9)

!
!----------------------------------------------------------------------
!
      return
      end
