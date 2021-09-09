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

!**s/r out_open_file - open an output FST file

subroutine out_open_file ( F_prefix_S )
      use timestr_mod, only: timestr_prognum, timestr_date_S
      use cstv
      use gem_options
      use out3
      use out_mod
      use ptopo
      use step_options

      implicit none
#include <arch_specific.hf>

      character*(*), intent(in ) :: F_prefix_S

!author 
!     Michel Desgagne  -  summer 2015
!revision
! v4_8   - Desgagne M.       - initial version
! v4_812 - Winger K. (ESCER/UQAM) - Write pilot files always as DAY-files
!                                 - Make sure 00Z oes to previous day
!                                 - Write time step 0 always as *_000000p
!                                 - Write monthly files as *_YYYYMMn

#include <rmnlib_basics.hf>

      character*4    unit_ext
      character*8    my_block
      character*16   datev,fdate
      character*1024 filen,myformat_S,my_hour
      integer prognum,err,i,indx,len0,len1
      real*8, parameter :: OV_day = 1.0d0/86400.0d0
      real*8  dayfrac
      character*4 curr_unit_S
      character*15 date_S, prognum_S
      integer date, time, stamp, stamp_prev
      integer interval, pp_int_step, next_out_step
!
!------------------------------------------------------------------
!
      if ( (Ptopo_couleur.gt.0) .or. (Out3_iome .lt. 0) )  return

      call datf2p (fdate, Out3_date)

      if ( Out_unf .gt. 0 ) return

      write(my_block,'(a,i3.3,a,i3.3)') '-',Ptopo_mycol,'-',Ptopo_myrow

      if ( trim(F_prefix_S) == 'casc' ) then

         dayfrac = dble(lctl_step-Step_delay) * Cstv_dt_8 * OV_day
         call incdatsd (datev,Step_runstrt_S,dayfrac)
         filen = trim(F_prefix_S)//'_'//trim(datev)//my_block

      else

! KW
         curr_unit_S = Out3_unit_S(1:3)
         interval    = max( int(Out3_close_interval), 1 )

         ! Write pilot files always as DAY or subdaily files (KW)
         if ( trim(F_prefix_S) == 'nm' ) then
            interval = 1
            if ( curr_unit_S(1:3) == 'MON' ) curr_unit_S = 'DAY'
         endif


         ! Unit extension
         unit_ext = ' '
         if (curr_unit_S(1:3) == 'SEC') unit_ext = 's'
         if (curr_unit_S(1:3) == 'MIN') unit_ext = 'm'
         if (curr_unit_S(1:3) == 'HOU') unit_ext = 'h'
         if (curr_unit_S(1:3) == 'DAY') unit_ext = 'd'
         if (curr_unit_S(1:3) == 'STE') unit_ext = 's'  ! originally 'p'  (KW)
         if (curr_unit_S(1:3) == 'MON') unit_ext = 'n'
         if (trim(F_prefix_S) == 'nm' ) unit_ext = ''   ! no extension for pilot files (KW)

! KW
         ! Monthly files
         if ( curr_unit_S(1:3) == 'MON' ) then
           ! Use Step_kount-1 to make sure 00Z goes to the previous day/month
           date_S = timestr_date_S(Out_dateo,sngl(Cstv_dt_8),Step_kount-1)

         ! Daily and sub-daily files
         else
           ! Determine file interval in steps
           select case (curr_unit_S(1:3))
              case ('DAY') ; pp_int_step =       int(interval * 86400. / Cstv_dt_8)
              case ('HOU') ; pp_int_step = max(1,int(interval *  3600. / Cstv_dt_8))
              case ('MIN') ; pp_int_step = max(1,int(interval *    60. / Cstv_dt_8))
              case ('SEC') ; pp_int_step = max(1,int(interval          / Cstv_dt_8))
              case ('STE') ; pp_int_step = max(1,int(interval                     ))
           end select
!print *,'out_outdir pp_int_step:',pp_int_step

           next_out_step = Step_kount
           if ( mod(Step_kount,pp_int_step) .ne. 0 ) then
             next_out_step = min(Step_kount + pp_int_step - mod(Step_kount,pp_int_step), Out_endstepno)
           else
             next_out_step = Step_kount
           endif

           ! For 'DAY' make sure 00Z goes to previous day
           if ( curr_unit_S(1:3) == 'DAY' ) next_out_step = next_out_step - 1

           ! Output date
           date_S = timestr_date_S(Out_dateo,sngl(Cstv_dt_8),next_out_step)

         endif

!print*,'out_open_file new date:',date_S


!print *,'out_open_file'
!print *,'out_open_file date_S:',date_S

         ! For hourly, minutely and secondly files set 00Z00m00s to 24Z00m00s of the previous day
         if ( curr_unit_S(1:3) == 'HOU' .or. &
              curr_unit_S(1:3) == 'MIN' .or. &
              curr_unit_S(1:3) == 'SEC' ) then
            read (date_S (10:15),'(i6.6)') time
            if ( time .eq. 0 ) then
              read (date_S ( 1: 8),'(i8.8)') date
              err = newdate (stamp, date, 0, 3)
              call incdatr (stamp_prev, stamp, -1.0_8)
!print *,'out_open_file stamp, stamp_prev:',stamp, stamp_prev
              err = newdate (stamp_prev, date, time, -3)
!print *,'out_open_file date, time:',date, time
              write (date_S,'(i8.8,a,i6.6)') date, '.', 240000
!print *,'out_open_file date_S:',date_S
            endif
         endif

         prognum_S = ''
         if ( lctl_step == 0 ) then
            prognum_S = '00000000'
            unit_ext  = 'p'
         elseif ( curr_unit_S(1:3) == 'MON' ) then
            ! Write monthly  files as *_YYYYMMn (KW)
            prognum_S = date_S(1:6)
         elseif ( curr_unit_S(1:3) == 'DAY' ) then
            ! Write daily    files as *_YYYYMMDDd (KW)
            prognum_S = date_S(1:8)
         elseif ( curr_unit_S(1:3) == 'HOU' ) then
            ! Write hourly   files as *_YYYYMMDD.hhh   and pilot files as *_YYYYMMDD.hh0000 (KW)
            prognum_S = date_S(1:11)
         elseif ( curr_unit_S(1:3) == 'MIN' ) then
            ! Write minutely files as *_YYYYMMDD.hhmmm and pilot files as *_YYYYMMDD.hhmm00 (KW)
            prognum_S = date_S(1:13)
         else
            ! Write secondly files and files per steps as *_YYYYMMDD.hhmmsss(KW)
            prognum_S = date_S(1:15)
         endif

         filen= trim(F_prefix_S)//fdate(1:8)//fdate(10:11)// &
                     my_block//'_'//trim(prognum_S)//trim(unit_ext)
!print *,'out_open_file filen:',trim(filen)
!if (lctl_step.eq.50) stop
      endif

      filen= trim(Out_dirname_S)//'/'//trim(filen)

      err= fnom  ( Out_unf, trim(filen), 'STD+RND', 0 )
      err= fstouv( Out_unf, 'RND' )
!
!------------------------------------------------------------------
!
      return
      end

