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

!**s/r  pre_jacobio3D -  BlocJacobi_3D additive-Schwarz preconditioner with
!                        local solver = "separable" approximate elliptic problem
!
      subroutine pre_jacobi3D ( Sol,Rhs,evec_local,Ni,Nj,Nk,ai,bi,ci )
      use glb_ld
      use opr
      use ptopo
      implicit none
#include <arch_specific.hf>

      integer Ni,Nj,Nk
      real*8 Rhs(Ni,Nj,Nk),Sol(Ni,Nj,Nk)
      real*8  ai(Ni,Nj,Nk), bi(Ni,Nj,Nk), ci(Ni,Nj,Nk)
      real*8 evec_local(Ni,Ni)

!author
!       Abdessamad Qaddouri -  2013
!
!revision
! v4_70 - Qaddouri A.       - initial version


      integer i,j,k,jr,offi,offj
      real*8 fdg(Ni,Nj,Nk), w2_8(Ni,Nj,Nk)
!
!     ---------------------------------------------------------------

      offi = Ptopo_gindx(1,Ptopo_myproc+1)-1
      offj = Ptopo_gindx(3,Ptopo_myproc+1)-1




      do j=1,Nj
         call dgemm ( 'N','N', Ni, nk, nk, 1.0D0,     &
                      rhs(1,j,1), Ni*Nj, Opr_lzevec_8,&
                      G_nk,0.0d0, w2_8 (1,j,1), Ni*Nj )
         do k=1,Nk
            do i = 1, Ni
               w2_8(i,j,k)= Opr_opsxp0_8(G_ni+offi+i+pil_w) * &
                    Opr_opsyp0_8(G_nj+offj+j+pil_s) * w2_8(i,j,k)
            enddo
         end do
     end do



      do k=1,Nk
         call dgemm ( 'T','N',Ni,Nj,Ni,1.0d0,evec_local,Ni,&
                      w2_8(1,1,k),Ni,0.0d0,fdg(1,1,k),Ni)
         do j =2, Nj
            jr =  j - 1
            do i=1,Ni
               fdg(i,j,k) = fdg(i,j,k) - ai(i,j,k)*fdg(i,jr,k)
            enddo
         enddo
         j = Nj
         do i=1,Ni
            fdg(i,j,k) = fdg(i,j,k)/bi(i,j,k)
         enddo
         do j = Nj-1, 1, -1
            jr =  j + 1
            do i=1 , Ni
               fdg(i,j,k)=(fdg(i,j,k)-ci(i,j,k)*fdg(i,jr,k))/bi(i,j,k)
            enddo
         enddo

         call dgemm ( 'N','N',Ni,Nj,Ni,1.0d0,evec_local,Ni,&
                      fdg(1,1,k),Ni,0.d0,w2_8(1,1,k),Ni )
      enddo



      do j=1,Nj
         call dgemm ( 'N','T', Ni, nk, nk, 1.0D0,     &
                      w2_8(1,j,1), Ni*Nj, Opr_zevec_8,&
                      G_nk,0.0d0, Sol (1,j,1), Ni*Nj  )
      enddo



!
!     ---------------------------------------------------------------
!
      return
      end

