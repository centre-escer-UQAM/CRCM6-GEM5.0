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

!**s/r pe_all_topo - First initialization steps

      subroutine pe_all_topo
use iso_c_binding
      use step_options
      use gem_options
      use lun
      use rstr
      use path
      use clib_itf_mod
      use ptopo
      implicit none
#include <arch_specific.hf>

!author
!     Michel Desgagne - Summer 2006

      include "rpn_comm.inc"

      integer,parameter :: BUFSIZE=10000
      integer, external :: fnom,wkoffit,OMP_get_max_threads
      integer, external :: generate_unique_name

      character(len=3)   :: mycol_S,myrow_S
      character(len=33)  :: buffer
      character(len=1024):: fn, pwd_S, scratch_dir
      integer :: err,unf,nc
      integer, dimension(4) :: bcast_ptopo
      integer :: bufnml(BUFSIZE),bufoutcfg  (BUFSIZE), &
                 bufinphycfg(BUFSIZE)
!
!-------------------------------------------------------------------
!
      lun_out    = -1
      Lun_debug_L=.false.
      if (Ptopo_myproc == 0) lun_out= 6

      call gemtim4 ( Lun_out, 'STARTING GEMDM', .false. )
      call timing_start2 ( 2, 'INIT_GEM', 1)
!
! Broadcasts processor topology
!
      bcast_ptopo(1) = Ptopo_npex
      bcast_ptopo(2) = Ptopo_npey
      bcast_ptopo(3) = Ptopo_nthreads_dyn
      bcast_ptopo(4) = Ptopo_nthreads_phy

      call RPN_COMM_bcast (bcast_ptopo, 4, "MPI_INTEGER",0,"grid",err)

      Ptopo_npex         = bcast_ptopo(1)
      Ptopo_npey         = bcast_ptopo(2)
      Ptopo_nthreads_dyn = bcast_ptopo(3)
      Ptopo_nthreads_phy = bcast_ptopo(4)

      call RPN_COMM_bcast (Ptopo_bind_L, 1, "MPI_LOGICAL",0,"grid",err)

      Ptopo_npeOpenMP = OMP_get_max_threads()

      if (Ptopo_nthreads_dyn < 1) Ptopo_nthreads_dyn=Ptopo_npeOpenMP
      if (Ptopo_nthreads_phy < 1) Ptopo_nthreads_phy=Ptopo_npeOpenMP

      if (Lun_out > 0) then
         write (Lun_out, 8255) Ptopo_npex, Ptopo_npey, &
             Ptopo_npeOpenMP,Ptopo_nthreads_dyn, Ptopo_nthreads_phy
         write (Lun_out, 8256) trim(Path_work_S)
      endif
!
! Initializes Path_nml_S and Path_outcfg_S
!
      Path_nml_S    = trim(Path_work_S)//'/model_settings.nml'
      Path_outcfg_S = trim(Path_work_S)//'/output_settings'
      Path_phyincfg_S = trim(Path_input_S)//'/physics_input_table'
!
      err= rpn_comm_mype (Ptopo_myproc, Ptopo_mycol, Ptopo_myrow)
!
! Initializes OpenMP
!
      call set_num_threads ( Ptopo_nthreads_dyn, 0 )
!
! Reading namelist file Path_nml_S (blind read)
!
      if (Ptopo_myproc == 0) then
         call array_from_file(bufnml,size(bufnml),Path_nml_S)
         call array_from_file(bufoutcfg,size(bufoutcfg),Path_outcfg_S)
         call array_from_file(bufinphycfg,size(bufinphycfg),Path_phyincfg_S)
      endif

! Changing directory to local Ptopo_mycol_Ptopo_myrow

      write(mycol_S,10) Ptopo_mycol
      write(myrow_S,10) Ptopo_myrow
10    format(i3.3)

      if (Ptopo_myproc == 0) call mkdir_gem ('./',Ptopo_npex,Ptopo_npey)
      call rpn_comm_barrier("grid", err)

      err= clib_chdir(mycol_S//'-'//myrow_S)
!
! Writing local namelist file Path_nml_S (blind write)
!
      if (clib_getenv ('GEM_scratch_dir',scratch_dir) < 0) scratch_dir='/tmp'
      nc = generate_unique_name(buffer,len(buffer))
      Path_nml_S      = trim(scratch_dir)//'/model_settings_'//buffer(1:nc)
      Path_outcfg_S   = trim(scratch_dir)//'/output_settings_'//buffer(1:nc)
      Path_phyincfg_S = trim(scratch_dir)//'/physics_input_table_'//buffer(1:nc)

      call RPN_COMM_bcast(bufnml,size(bufnml),"MPI_INTEGER",0, &
                                                 "grid",err )
      call RPN_COMM_bcast(bufoutcfg,size(bufoutcfg),"MPI_INTEGER",0, &
                                                 "grid",err )
      call RPN_COMM_bcast(bufinphycfg,size(bufinphycfg),"MPI_INTEGER",0, &
                                                 "grid",err )

      call array_to_file (bufnml,size(bufnml),trim(Path_nml_S))
      call array_to_file (bufoutcfg,size(bufoutcfg),trim(Path_outcfg_S))
      call array_to_file (bufinphycfg,size(bufinphycfg),trim(Path_phyincfg_S))

      Lun_rstrt = 0 ; Rstri_rstn_L= .false. ; Step_kount= 0
      err = wkoffit('gem_restart')
      if (err >= -1) then
         err= fnom( Lun_rstrt,'gem_restart','SEQ+UNF+OLD',0 )
         if (err >= 0) then
           Rstri_rstn_L = .true.
           if (lun_out > 0) write (lun_out,1001)
           call rdrstrt
         endif
      endif
!
! Determine theoretical mode with presence of file ${TASK_WORK}/theoc
!
      unf=0
      Schm_theoc_L = .false.
      fn=trim(Path_work_S)//'/theoc'
      if (wkoffit(fn) > -3) then
         if (Ptopo_myproc == 0) write (Lun_out,*) &
                                'Assume Theoretical case'
         Schm_theoc_L = .true.
      else
         call fclos (unf)
      endif

 1001 format (/' RESTART DETECTED'/)
 8255 format (/," MPI CONFIG (npex x npey): ",i4,' x ',i3,/, &
          " OMP CONFIG (npeOpenMP x nthreads_dyn x nthreads_phy): ",&
            i4,' x ',i3,' x ',i3)
 8256 format (/," WORKING DIRECTORY:"/a/)
!
!-------------------------------------------------------------------
!
      return
      end

      subroutine mkdir_gem (F_path_S,F_npex,F_npey)
      use lun
      use path
      use clib_itf_mod
      implicit none

      character(len=*), intent(IN) :: F_path_S
      integer      , intent(IN) :: F_npex,F_npey

      character(len=2048) ici
      integer :: i,j,status,err
      integer :: mk_gem_dir

      err= clib_getcwd(ici)
      err= clib_chdir(trim(F_path_S))

      do j=0,F_npey-1
      do i=0,F_npex-1
         status = mk_gem_dir(i,j)
      enddo
      enddo

      err= clib_chdir(ici)

return
end

integer function mk_gem_dir(x,y)
use ISO_C_BINDING
      use lun
      use path
implicit none
integer, intent(IN) :: x, y
character (len=7) :: name

character (len=1), dimension(8), target :: temp

interface
integer(C_INT) function mkdir(path,mode) bind(C,name='mkdir')
use ISO_C_BINDING
type(C_PTR), value :: path
integer(C_INT), value :: mode
end function mkdir
end interface
write(name,100)x,'-',y
100 format(I3.3,A1,I3.3)

temp = transfer( trim(name)//achar(0) , temp )
mk_gem_dir = mkdir(C_LOC(temp),511)

return
end function mk_gem_dir
