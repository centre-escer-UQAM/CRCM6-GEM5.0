!----------------------------------LICENCE BEGIN -------------------------------
!     GEM - Library of kernel routines for the GEM numerical atmospheric model
!     Copyright (C) 1990-2010 - Division de Recherche en Prevision Numerique
!     Environnement Canada
!     This library is free software; you can redistribute it and/or modify it
!     under the terms of the GNU Lesser General Public License as published by
!     the Free Software Foundation, version 2.1 of the License. This library is
!     distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
!     without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!     PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
!     You should have received a copy of the GNU Lesser General Public License
!     along with this library; if not, write to the Free Software Foundation, Inc.,
!     59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
!----------------------------------LICENCE END ---------------------------------

!**s/r hzd_exp5p  - 5 points explicit horizontal diffusion operator
!                   which includes geometric terms (hzd_geom_*)

      subroutine hzd_exp5p ( F_champ, F_temp, Minx,Maxx,Miny,Maxy, NK,&
                                        F_coef_8, F_arakawa_S,mm,dpwr )
      use gem_options
      use hzd_mod
      use glb_ld
      implicit none
#include <arch_specific.hf>

      character*1 F_arakawa_S
      integer Minx,Maxx,Miny,Maxy,NK,mm,dpwr
      real F_champ(Minx:Maxx,Miny:Maxy,NK)
      real F_temp(Minx:Maxx,Miny:Maxy,NK)
      real*8 F_coef_8(NK)

!author
!    Abdessamad Qaddouri - summer 2015
!
!revision
! v4_80 - Qaddouri A.      - initial version
! v4_80 - Desgagne & Lee   - optimization
!


      integer i,j,k,i0,in,j0,jn
      real*8 wk_8 (Minx:Maxx,Miny:Maxy,Nk)
      real*8, dimension(:,:), pointer :: stencils => null()
!
!---------------------------------------------------------------------
!
      i0  = 2        - G_halox*(1-west )
      j0  = 2        - G_haloy*(1-south)
      in  = l_ni - 1 + G_halox*(1-east )
      jn  = l_nj - 1 + G_haloy*(1-north)

      select case (F_arakawa_S)
      case ('M')
         stencils => Hzd_geom_q
      case ('U')
         stencils => Hzd_geom_u
         in  = l_niu - 1 + G_halox*(1-east )
      case ('V')
         stencils => Hzd_geom_v
         jn  = l_njv - 1 + G_haloy*(1-north)
      end select

      if (mm == 1) then
          F_champ(:,:,1:NK) = F_temp(:,:,1:NK)
      else if (mm == 2) then
          F_temp(:,:,1:NK) = F_champ(:,:,1:NK) - F_temp(:,:,1:NK)
      else
          F_temp(:,:,1:NK) = F_champ(:,:,1:NK) - F_temp(:,:,1:NK)
      endif

!     del2 explicit
!$omp parallel private(i,j,k) shared (i0,in,j0,jn,wk_8,stencils)
!$omp do
       do k=1,NK
          do j= j0, jn
          do i= i0, in
             wk_8(i,j,k)= stencils(j,1)*F_temp(i  ,j  ,k) + &
                          stencils(j,2)*F_temp(i-1,j  ,k) + &
                          stencils(j,3)*F_temp(i+1,j  ,k) + &
                          stencils(j,4)*F_temp(i  ,j-1,k) + &
                          stencils(j,5)*F_temp(i  ,j+1,k)
          enddo
          enddo
       enddo
!$omp end do
!$omp end parallel

      if (mm == dpwr) then

!     Apply diffusion operator
!$omp parallel private(i,j,k) shared(i0,in,j0,jn,wk_8)
!$omp do
          do k=1,NK
             do j= j0, jn
             do i= i0, in
                F_champ(i,j,k)= F_champ(i,j,k) + F_coef_8(k)*wk_8(i,j,k)
             enddo
             enddo
          enddo
!$omp end do
!$omp end parallel
          F_temp(:,:,1:NK) = F_champ(:,:,1:NK)

      else
!$omp parallel private(i,j,k) shared(i0,in,j0,jn,wk_8)
!$omp do
          do k=1,NK
             do j= j0, jn
             do i= i0, in
                F_temp(i,j,k)= F_champ(i,j,k) + F_coef_8(k)*wk_8(i,j,k)
             enddo
             enddo
          enddo
!$omp end do
!$omp end parallel

      endif
!
!----------------------------------------------------------------------
!
      return
      end

