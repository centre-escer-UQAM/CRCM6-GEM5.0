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

!**s/r gem_ctrl - initiate the forward integration of the model
!
      subroutine gem_ctrl
      use gem_options
      use glb_ld
      use lun
      use rstr
      implicit none
#include <arch_specific.hf>

!author
!     Michel Roch - rpn - june 1993
!revision
! v2_00 - Desgagne M.       - initial MPI version (from gefcn2 v1_03)
! v2_10 - Tanguay M.        - cancel indata and p_main when 4D-Var
! v2_20 - Desgagne M.       - correction for Step_total=0
! v2_21 - Dugas B.          - adapt to climate mode
! v2_30 - Desgagne M.       - entry vertical interpolator in gemdm
! v2_31 - Tanguay M.        - initialization done_once and
!                             remove call to adv_set
! v3_21 - Tanguay M.        - do filtering at t=0 in v4d_predat
! v3_21 - valcke, S.        - add oasis coupling initialization
! v3_30 - Desgagne M.       - moved calls set_adw,set_opr to set_dync
!object
!     Beginning of the integration. This subroutine
!     reads the data and performs initialization if required.
!     It then initiates the forward intergration of the model.

      logical :: rstrt_L= .false.
      integer err
!
!     ---------------------------------------------------------------
!
      call gemtim4 ( Lun_out, 'GEM_CTRL: START', .false. )

      if ( .not. Rstri_rstn_L ) then

         call indata

      else

         call set_dync ( .true., err )

      endif

      call spn_init ()
      call set_smago ()

      call gemtim4 ( Lun_out, 'GEM_CTRL: INIT COMPLETED', .false. )
      call timing_stop ( 2 )

      if (  Init_mode_L ) call initial (rstrt_L)

      if ( .not.rstrt_L ) call gem_run (rstrt_L)

      if (Lun_out > 0) write(Lun_out,3000) Lctl_step

 3000 format(/,'GEM_CTRL: END OF CURRENT TIME SLICE AT TIMESTEP',I8, &
             /,'===================================================')
!
!     ---------------------------------------------------------------
!
      return
      end
