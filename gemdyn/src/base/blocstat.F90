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

!**s/r blocstat  - Performs 2D-3D statistics on model fields

      subroutine blocstat (F_forceit_L)
      use step_options
      use gem_options
      use glb_ld
      use gmm_itf_mod
      use ptopo
      implicit none
#include <arch_specific.hf>

      logical F_forceit_L


      logical,save :: done = .false.
      logical :: flag
      type(gmm_metadata) :: tmp_meta,mymeta
      integer i0,in,j0,jn,inn,jnn,n,indx,istat
      real, pointer, dimension(:,:  ) :: wk2d
      real, pointer, dimension(:,:,:) :: wk3d
      character(len=GMM_MAXNAMELENGTH) :: stat_varname
!
!     ---------------------------------------------------------------
!
      call timing_start2 ( 99, 'BLOCSTAT', 1 )

      if (.not.done) call set_statliste

      flag = .false.
      if (Step_gstat > 0) flag = (mod(Lctl_step,Step_gstat) == 0)
      flag = flag .or. F_forceit_L

      if (flag) then

         if (Ptopo_myproc == 0) write(6,1000) Lctl_step

         i0 = 1 ; in = G_ni
         j0 = 1 ; jn = G_nj

         do n=1,stat_nombre
            istat = gmm_getmeta(stat_liste(n),tmp_meta)
            if (istat == 0) then
               indx = index(stat_liste(n),"TR/")*3 + 1
               stat_varname = stat_liste(n)(indx:)
               inn= 0
               jnn= 0
               if ( stat_liste(n)(1:3)=='URT' ) inn=1
               if ( stat_liste(n)(1:3)=='VRT' ) jnn=1
               nullify (wk2d,wk3d)
               if (tmp_meta%l(3)%high == tmp_meta%l(3)%low) then
                  istat = gmm_get(stat_liste(n), wk2d, mymeta)
                  call glbstat2 (wk2d,stat_varname,'', &
                                 tmp_meta%l(1)%low,tmp_meta%l(1)%high, &
                                 tmp_meta%l(2)%low,tmp_meta%l(2)%high, &
                                 1,1, &
                                 i0,in-inn,j0,jn-jnn, 1,1)
               else
                  istat = gmm_get(stat_liste(n), wk3d, mymeta)
                  call glbstat2 (wk3d,stat_varname,'', &
                                tmp_meta%l(1)%low,tmp_meta%l(1)%high,&
                                tmp_meta%l(2)%low,tmp_meta%l(2)%high,&
                                tmp_meta%l(3)%low,tmp_meta%l(3)%high,&
                                i0,in-inn,j0,jn-jnn                 ,&
                                tmp_meta%l(3)%low,tmp_meta%l(3)%high)
               endif
            endif
         end do

!         if ((Ptopo_numproc == 1).and.(Schm_theoc_L)) then
!            call lipschitz(u, v, zd, LDIST_DIM, G_Nk, i0,in,j0,jn)
!         endif

!     Print max courrant numbers
         if (Ptopo_myproc==0 .and. Lctl_step>0) call adv_cfl_print()

         if (Ptopo_myproc == 0) write(6,1001)

      endif

      done = .true.
      call timing_stop ( 99 )
      call flush(6)

 1000 format (/ 19('#'),' BLOC STAT ',i6,1X,19('#'))
 1001 format (  19('#'),' BLOC STAT ...done')
!     ---------------------------------------------------------------
!
      return
      end

!**s/r set_statliste
!
      subroutine set_statliste
      use gem_options
      use tr3d
      use gmm_itf_mod
      implicit none
#include <arch_specific.hf>
!
!author
!     M. Desgagne  	- fall 2009
!
!revision
! v4_10 - Desgagne M.       - initial version
!

      character(len=GMM_MAXNAMELENGTH) :: tmp_liste (1000), dumc_S
      integer k,n,cnt
!     ---------------------------------------------------------------
!
      if (stat_liste(1) == '') then
         stat_liste(1) = 'ALL_DYN_T1'
         stat_liste(2) = 'TR/HU:P'
      endif

      cnt = 0
      do k=1,stat_maxn
         if (stat_liste(k) == '') exit
         call low2up  (stat_liste(k),dumc_S)
         stat_liste(k) = dumc_S

         if ((stat_liste(k) == 'ALL_DYN_T1') .or. (stat_liste(k) == 'ALL')) then
            cnt = cnt + 1
            tmp_liste(cnt  ) = 'URT1'
            tmp_liste(cnt+1) = 'VRT1'
            tmp_liste(cnt+2) = 'ZDT1'
            tmp_liste(cnt+3) = 'WT1'
            tmp_liste(cnt+4) = 'TT1'
            tmp_liste(cnt+5) = 'ST1'
            cnt = cnt + 5
            if(.not.Schm_hydro_L)then
               tmp_liste(cnt+1) = 'QT1'
               cnt = cnt + 1
            endif
            if (stat_liste(k) /= 'ALL') cycle
         endif
         if (stat_liste(k) == 'ALL_DYN_T0') then
            cnt = cnt + 1
            tmp_liste(cnt  ) = 'URT0'
            tmp_liste(cnt+1) = 'VRT0'
            tmp_liste(cnt+2) = 'ZDT0'
            tmp_liste(cnt+3) = 'WT0'
            tmp_liste(cnt+4) = 'TT0'
            tmp_liste(cnt+5) = 'ST0'
            cnt = cnt + 5
            if(.not.Schm_hydro_L)then
               tmp_liste(cnt+1) = 'QT0'
               cnt = cnt + 1
            endif
            cycle
         endif
         if ((stat_liste(k) == 'ALL_TR_T1') .or. (stat_liste(k) == 'ALL')) then
            do n=1,tr3d_ntr
               cnt = cnt + 1
               tmp_liste(cnt) = 'TR/'//trim(Tr3d_name_S(n))//':P'
            end do
            cycle
         endif
         if (stat_liste(k) == 'ALL_TR_T0') then
            do n=1,tr3d_ntr
               cnt = cnt + 1
               tmp_liste(cnt) = 'TR/'//trim(Tr3d_name_S(n))//':M'
            end do
            cycle
         endif
         cnt = cnt + 1
         tmp_liste(cnt) = stat_liste(k)
      end do

      stat_liste (1:cnt) = tmp_liste (1:cnt)
      stat_nombre        = cnt
!
!     ---------------------------------------------------------------
!
      return
      end
