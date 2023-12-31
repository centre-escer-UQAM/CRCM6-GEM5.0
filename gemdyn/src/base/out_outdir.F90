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
      character*15 curr_date_S                ! YYYYMMDD.hhmmss
      integer      err, stepno, prognum, prognum1
      integer      date, year, month, time, stamp, stamp_last, t
      integer      pp_int_step, next_out_step
      integer, save :: last_out_step = -1
      integer :: interval
      real*8       nhours_8

      logical, save :: new_dir_L = .true.
      integer, save :: out_int   = -1
!
!----------------------------------------------------------------------
!

!print *,'In out_outdir'
!print *,'out_outdir Step_kount    = ',Step_kount
!print *,'out_outdir startofroutine: out_int = ',out_int

!print *,'out_outdir Step_total    = ',Step_total
!print *,'out_outdir Out_dateo     = ',Out_dateo
!print *,'out_outdir Out_endstepno = ',Out_endstepno

!print *,'Step_initial = ',Step_initial
!print *,'Step_delay   = ',Step_delay
!      upperlimit = Step_total + Step_initial

      Out_post_L = .false.
      if (Step_kount == 0 .or. Step_kount == 1) out_int = -1

      call out_steps

      ! if there is no output for this time step => return (KW)
      if ( outd_sorties(0,Lctl_step) .eq. 0 .and.  &
           outp_sorties(0,Lctl_step) .eq. 0 .and.  &
           outc_sorties(0,Lctl_step) .eq. 0 ) return

      if ( Init_mode_L .and. (Step_kount.ge.Init_halfspan) ) return

!      if ( Out3_unit_S(1:3) .eq. 'MON' ) Out3_close_interval = 1.
      Out3_close_interval = 1   ! KW
      Out3_postproc_fact  = max(1, Out3_postproc_fact)

!print *,'out_outdir Out_endstepno:',Out_endstepno
      interval = int(Out3_close_interval * Out3_postproc_fact)

      ! Decide if new directory needs to get created
      if ( out_int    == interval .or. &
           out_int    == -1 ) then
        new_dir_L = .true.

        ! Reset counter
        out_int = 0

      else
        new_dir_L = .false.
      endif

!print *,'out_outdir interval,out_int,new_dir_L:',interval,out_int,new_dir_L




      ! Make sure 00Z goes to the previous day (KW)

      ! Current date
      curr_date_S = timestr_date_S(Out_dateo,sngl(Cstv_dt_8),Step_kount)
!print *,'out_outdir curr_date_S:',curr_date_S

      ! Determine saving interval
      select case (Out3_unit_S(1:3))
         case ('MON') ! Monthly files : First of the month at 00Z
                      if (curr_date_S( 7:15) == '01.000000') out_int = out_int + 1

         case ('DAY') ! Daily files   : 00Z of the day
                      if (curr_date_S(10:15) ==    '000000') out_int = out_int + 1

         case ('HOU') ! Hourly files  : 00m of the hour
                      if (curr_date_S(12:15) ==      '0000') out_int = out_int + 1

         case ('MIN') ! Minutely files: 00s of the minute
                      if (curr_date_S(14:15) ==        '00') out_int = out_int + 1

         case ('SEC') ! Secondly files
                      out_int = out_int + 1

         case ('STE') ! Step files
                      out_int = out_int + 1
      end select



      ! Determine if post processing should get started
      !   Create new directory after Out3_postproc_fact
      !   Always create new output directory at time steps 0 and 1
      !   Always create new output directory at beginning of a month
      !   Always start post processing at the end of a job
      if ( out_int    == interval            .or. &
           Step_kount == 0                   .or. &
           Step_kount == 1                   .or. &
           curr_date_S( 7:15) == '01.000000' .or. &
           stepno == Out_endstepno) then
        Out_post_L = .true.
      else
        Out_post_L = .false.
      endif


      ! Always create a new directory at the beginning of a month
      ! Set flag for next timestep
      if (curr_date_S( 7:15) == '01.000000') out_int = -1


!print *,'out_outdir interval,out_int,Out_post_L:',interval,out_int,Out_post_L


!print *,'out_outdir Step_job, pp_int_step:',Step_job, pp_int_step
!print *,'out_outdir Step_kount, Step_job:',Step_kount, Step_job

!print *,'out_outdir Step_kount,mod(Step_kount,pp_int_step):',Step_kount,mod(Step_kount,pp_int_step)

!print *,'out_outdir next_out_step,Out_post_L:', next_out_step,Out_post_L


      ! Create new output directory
      if ( new_dir_L ) then

        ! Give special name for time step 0
        if ( Step_kount.eq.0 ) then
           postjob_S = '00000000.' // curr_date_S(1:6)
        else
           postjob_S = curr_date_S
        endif

        Out_laststep_S = 'firststep_'//postjob_S
        Out_dirname_S  = trim(Path_output_S)//'/'//trim(Out_laststep_S) !  Add 'trim' (KW)
!print *,'out_outdir Step_kount,Out_laststep_S: ',Step_kount,Out_laststep_S
!print *,'out_outdir Out_dirname_S :',Out_dirname_S
!print *,'out_outdir Out_post_L:',Out_post_L
        ! PE0 is responsible for creating shared subdir structure
        if (dirstep_S /= Out_dirname_S) then
           dirstep_S = Out_dirname_S        
           if (Ptopo_myproc == 0 .and. Ptopo_couleur == 0) then
              err = clib_mkdir(trim(Out_dirname_S))
print *,'out_outdir: Create ',trim(Out_dirname_S)
!print *,'out_outdir: err: ',err
              if (Lun_out>0) write(Lun_out,1001) trim(Out_laststep_S),Step_kount
           endif
        endif
      
        ! Wait for Grid PE0 to be finished subdir creation
        call rpn_comm_barrier (RPN_COMM_ALLGRIDS, err)

        ! Each io pe now creates an independent subdir for outputs
        write (blocxy_S,'(I3.3,"-",I3.3)') Ptopo_mycol, Ptopo_myrow
        Out_dirname_S = trim(Out_dirname_S)//'/'//blocxy_S
        err = CLIB_OK
        if (Out3_iome .ge. 0 .and. dirbloc_S /= Out_dirname_S &
                             .and. Ptopo_couleur == 0) then
           dirbloc_S = Out_dirname_S
           err = clib_mkdir ( trim(Out_dirname_S) )
           err = clib_isdir ( trim(Out_dirname_S) )
        endif

        call gem_error (err,'out_outdir','unable to create output directory structure')

        ! Reset counter
!        out_int = 0

      endif ! Create output directory

!print *,'out_outdir end of routine: out_int = ',out_int
      
 1001 format (' OUT_OUTDIR: DIRECTORY output/',a,' was created at timestep: ',i9)

!
!----------------------------------------------------------------------
!
      return
      end
