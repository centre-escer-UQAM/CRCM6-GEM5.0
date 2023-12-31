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

!**s/r init_component

      subroutine init_component()
      use clib_itf_mod
      use component_mod
      use dcst
      use iso_c_binding
      use glb_ld
      use grid_options
      use path
      use ptopo
      use step_options
      use tdpack
      use version
      implicit none
#include <arch_specific.hf>

!author   Michel Desgagne -- Spring 2013 --
!
!revision
! v4_50 - M. Desgagne       - Initial version
! v4_80 - M. Desgagne       - remove *bloc*, introduce RPN_COMM_*_io_*
! v5_00 - K. Winger (ESCER/UQAM) - set Step_alarm= 60 resp. 600
!
      include "rpn_comm.inc"

      external init_ndoms, pe_zero_topo
      integer, external :: model_timeout_alarm

      character(len=256) :: my_dir
      integer :: ierr,mydomain,ngrids
!
!--------------------------------------------------------------------
!
      call gemtim4 ( 6, '', .false. )
      Step_alarm= 600
      ierr= model_timeout_alarm (Step_alarm)
      call rpn_comm_mydomain (init_ndoms, mydomain)

      write(my_dir,'(a,i4.4)') 'cfg_',mydomain

      ierr = clib_getenv ('TASK_BASEDIR',Path_basedir_S)
      ierr = clib_getenv ('TASK_WORK'   ,Path_work_S   )
      ierr = clib_getenv ('TASK_INPUT'  ,Path_input_S  )
      ierr = clib_getenv ('TASK_OUTPUT' ,Path_output_S )

      Path_input_S  = trim(Path_input_S ) // '/' // trim(my_dir)
      Path_work_S   = trim(Path_work_S  ) // '/' // trim(my_dir)
      Path_output_S = trim(Path_output_S) // '/' // trim(my_dir)

      ierr = clib_chdir (trim(Path_work_S))

      COMPONENT= 'MOD'

      ngrids=1
      if (Grd_yinyang_L) ngrids=2

      Ptopo_couleur= RPN_COMM_init_multi_level (             &
                   pe_zero_topo, Ptopo_myproc,Ptopo_numproc, &
                   Ptopo_npex,Ptopo_npey,Grd_ndomains,ngrids )

      Step_alarm= 600
      ierr= model_timeout_alarm (Step_alarm)

      call timing_init2 ( Ptopo_myproc, COMPONENT )
      call timing_start2 ( 1, 'GEMDM', 0)

      if (Grd_yinyang_L) then
         Ptopo_intracomm = RPN_COMM_comm ('GRID')
         Ptopo_intercomm = RPN_COMM_comm ('GRIDPEERS')
         call RPN_COMM_size ('MULTIGRID',Ptopo_world_numproc,ierr)
         call RPN_COMM_rank ('MULTIGRID',Ptopo_world_myproc ,ierr)

         if (Ptopo_couleur == 0) Grd_yinyang_S = 'YIN'
         if (Ptopo_couleur == 1) Grd_yinyang_S = 'YAN'
         ierr= clib_chdir(trim(Grd_yinyang_S))
         Ptopo_ncolors = 2
      else
         Ptopo_couleur = 0
         Ptopo_ncolors = 1
      endif

!      call RPN_COMM_size ( RPN_COMM_GRIDPEERS, Ptopo_nodes, ierr )

      call msg_set_can_write (Ptopo_myproc == 0)

      Dcst_rayt_8      = rayt_8
      Dcst_inv_rayt_8  = 1.d0 / rayt_8
      Dcst_omega_8     = omega_8

      call pe_all_topo

      ierr= 0

      G_periodx= .false.  ;  G_periody= .false.

      l_west  = (0 == Ptopo_mycol)
      l_east  = (Ptopo_npex-1 == Ptopo_mycol)
      l_south = (0 == Ptopo_myrow)
      l_north = (Ptopo_npey-1 == Ptopo_myrow)

      north= 0 ; south= 0 ; east= 0 ; west= 0
      if (l_north) north = 1
      if (l_south) south = 1
      if (l_east ) east  = 1
      if (l_west ) west  = 1
!
!--------------------------------------------------------------------
!
      return
      end

