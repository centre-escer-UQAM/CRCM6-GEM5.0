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

!s/r wrrstrt - Write the restart file
!
      subroutine wrrstrt ()
use iso_c_binding
      use phy_itf, only: phy_restart
      use step_options
      use gem_options
      use lun
      use psadjust
      use gmm_itf_mod
      use wb_itf_mod
      use timestr_mod, only: timestr_date_S  ! Needed to write last date (KW)
      use out_mod
      use cstv

      implicit none

!revision
! v5_00 - Winger K. (ESCER/UQAM) - Write restart file twice
!                                - Write last date to be able to remove eventually
!                                  existing output files from previous attempts

#include <arch_specific.hf>
#include <clib_interface_mu.hf>  ! Needed to write double restart file (KW)


      include "rpn_comm.inc"

      integer, external :: fnom,fclos
      integer ier,gmmstat,me,howmany,newcomm,i
      character* 2 :: double_restart_S  ! Needed to write double restart file (KW)
      character*15 :: date_S            ! Needed to write last date - 19790101.120000 (KW)
      integer oun
!
!     ---------------------------------------------------------------
!
      if (Lun_out > 0) write(Lun_out,2000) Lctl_step

      ! Write file with last date (KW)
      if (Lun_out > 0) then
         date_S  = timestr_date_S(Out_dateo,sngl(Cstv_dt_8),Step_kount)
         oun = 0
         ier = fnom (oun,'../last_date','SEQ',0)
         write(oun,'(a15)') date_S
         ier = fclos(oun)
print *,'In wrrstrt last_date: ', date_S
      end if

      call split_on_hostid (RPN_COMM_comm('GRID'),me,howmany,newcomm)

      call timing_start2 ( 33, 'RESTART', 34 )
      do i=0,howmany-1
         if (i == me) then

            Lun_rstrt = 0
            ier = fnom (Lun_rstrt,'gem_restart','SEQ+UNF',0)

            write(Lun_rstrt) Lctl_step,Step_kount,Init_mode_L
            write(Lun_rstrt) PSADJ_g_avg_ps_initial_8,PSADJ_scale_8,PSADJ_fact_8

            ier = fclos(Lun_rstrt)

            !        Write Gmm-files

            gmmstat = gmm_checkpoint_all(GMM_WRIT_CKPT)
            ier = wb_checkpoint()

            ier = phy_restart ('W', .false.)

         endif

         call mpi_barrier (newcomm,ier)

      end do
      call timing_stop (33)

 2000 format(/,'WRITING A RESTART FILE AT TIMESTEP #',I8, &
             /,'=========================================')
!
!     ---------------------------------------------------------------
!
      return
      end
