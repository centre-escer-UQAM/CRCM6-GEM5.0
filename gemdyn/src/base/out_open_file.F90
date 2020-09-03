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
      character*15 date_S, dateo_S, prognum_S
      integer date, time, stamp, stamp_prev, hour, minute, second
      integer interval, pp_int_step, pp_add_step, next_out_step
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
!         interval    = max( int(Out3_close_interval), 1 )
         interval    = 1

         ! Unit extension
         unit_ext = ''

         ! Write pilot files always as DAY or subdaily files (KW)
         if ( trim(F_prefix_S) == 'nm' ) then
            interval = 1
            curr_unit_S = Out3_pilot_unit_S(1:3)
         endif


         ! Use "Step_kount-1" to make sure 00Z goes to the previous month/day/hour/minute
         select case (curr_unit_S(1:3)) 
            case ('MON') ! Monthly files
                         unit_ext = 'n'
                         date_S = timestr_date_S(Out_dateo,sngl(Cstv_dt_8),Step_kount-1)
                         prognum_S = date_S(1: 6) ! Write as *_YYYYMMn

            case ('DAY') ! Daily files
                         unit_ext = 'd'
                         date_S = timestr_date_S(Out_dateo,sngl(Cstv_dt_8),Step_kount-1)
                         prognum_S = date_S(1: 8) ! Write as *_YYYYMMDDd

            case ('HOU') ! Hourly files
                         unit_ext = 'h'
                         date_S = timestr_date_S(Out_dateo,sngl(Cstv_dt_8),Step_kount-1+int(3600/sngl(Cstv_dt_8)))
                         prognum_S = date_S(1:11) ! Write as *_YYYYMMDD.hhh

            case ('MIN') ! Minutely files
                         unit_ext = 'm'
                         if ( Cstv_dt_8 >= 60. ) then
                            date_S = timestr_date_S(Out_dateo,sngl(Cstv_dt_8),Step_kount)
                         else
                            date_S = timestr_date_S(Out_dateo,sngl(Cstv_dt_8),Step_kount-1+int(60/sngl(Cstv_dt_8)))
                         endif
                         prognum_S = date_S(1:13) ! Write as *_YYYYMMDD.hhmmm

            case default ! Secondly & step files
                         unit_ext = 's'
                         date_S = timestr_date_S(Out_dateo,sngl(Cstv_dt_8),Step_kount)
                         prognum_S = date_S(1:15) ! Write as *_YYYYMMDD.hhmmsss
         end select


!print*,'out_open_file new date:',date_S
!print *,'out_open_file next_out_step:', next_out_step
!print *,'out_open_file date_S:',date_S

         ! No unit extension for monthly or daily pilot files (KW)
         if ( trim(F_prefix_S) == 'nm' ) then
            if ( curr_unit_S(1:3) == "MON" .or. curr_unit_S(1:3) == "DAY" ) then
               unit_ext = ''
            end if
         endif


         if ( lctl_step == 0 ) then
            prognum_S = '00000000'
            unit_ext  = 'p'
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

