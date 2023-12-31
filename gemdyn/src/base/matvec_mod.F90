module matvec
   ! Matrix-vector product subroutines
   !
   ! Author
   !     Stephane Gaudreault / Abdessamad Qaddouri -- June 2014
   !
   ! Revision
   !     v4_70 - Gaudreault/Qaddouri      - initial version
   !

   use cstv
   use geomh
   use glb_ld
   use ldnh
   use opr
   use sol
   implicit none
   private

#include <arch_specific.hf>

   integer, parameter :: IDX_POINT=1, IDX_WEST=2, IDX_EAST=3, IDX_NORTH=4, IDX_SOUTH=5, IDX_TOP=6, IDX_BOTTOM=7
   real*8, dimension(:,:,:,:), allocatable :: matrix

   public :: matvec_init, matvec_3d

contains

   subroutine matvec_init()
      implicit none

      real*8  :: di_8
      real*8  :: xxx, yyy
      integer :: i, j, k, jj, ii

      allocate (matrix(7,ldnh_minx:ldnh_maxx, ldnh_miny:ldnh_maxy, l_nk))

      xxx = - Cstv_hco2_8
      yyy = - Cstv_hco1_8



      do k=1, l_nk
         do j=1+sol_pil_s, l_nj-sol_pil_n
            jj=j+l_j0-1
            di_8 = Opr_opsyp0_8(G_nj+jj) * geomh_invcy2_8(j)
            do i=1+sol_pil_w, l_ni-sol_pil_e
               ii=i+l_i0-1

               matrix(IDX_POINT,i,j,k) = Cstv_hco0_8 * (Opr_opszp2_8(G_nk+k) + Opr_opszpl_8(G_nk+k) &
                           + xxx * Opr_opszpm_8(G_nk+k) + yyy * Opr_opszp0_8(G_nk+k)) &
                           + Opr_opszp0_8(G_nk+k) * (Opr_opsxp2_8(G_ni+ii) * di_8     &
                           + Opr_opsxp0_8(G_ni+ii) * Opr_opsyp2_8(G_nj+jj))           &
                           / (Opr_opsxp0_8(G_ni+ii) * Opr_opsyp0_8(G_nj+jj))

               matrix(IDX_WEST,i,j,k) = Opr_opsxp2_8(ii) * Opr_opszp0_8(G_nk+k) * Opr_opsyp0_8(G_nj+jj) &
                           / (cos( G_yg_8 (jj) )**2) / (Opr_opsxp0_8(G_ni+ii)*Opr_opsyp0_8(G_nj+jj))

               matrix(IDX_EAST,i,j,k) = Opr_opsxp2_8(2*G_ni+ii) * Opr_opszp0_8(G_nk+k) * Opr_opsyp0_8(G_nj+jj) &
                          / (cos( G_yg_8 (jj) )**2) / (Opr_opsxp0_8(G_ni+ii) * Opr_opsyp0_8(G_nj+jj))

               matrix(IDX_SOUTH,i,j,k) = Opr_opsxp0_8(G_ni+ii) * Opr_opsyp2_8(jj) * Opr_opszp0_8(G_nk+k) &
                          / (Opr_opsxp0_8(G_ni+ii) * Opr_opsyp0_8(G_nj+jj))

               matrix(IDX_NORTH,i,j,k) = Opr_opsxp0_8(G_ni+ii) * Opr_opsyp2_8(2*G_nj+jj) * Opr_opszp0_8(G_nk+k) &
                           / (Opr_opsxp0_8(G_ni+ii) * Opr_opsyp0_8(G_nj+jj))

               matrix(IDX_TOP,i,j,k) = Cstv_hco0_8 * (Opr_opszp2_8(k) + Opr_opszpl_8(k) + xxx * Opr_opszpm_8(k))

               matrix(IDX_BOTTOM,i,j,k) = Cstv_hco0_8 * (Opr_opszp2_8(2*G_nk+k) + Opr_opszpl_8(2*G_nk+k) + xxx * Opr_opszpm_8(2*G_nk+k))

            enddo
         enddo
      enddo


   end subroutine matvec_init


   subroutine matvec_3d(vec, prod)
      implicit none
      real*8, dimension(ldnh_minx:ldnh_maxx, ldnh_miny:ldnh_maxy, l_nk), intent(in) :: vec
      real*8, dimension(ldnh_minx:ldnh_maxx, ldnh_miny:ldnh_maxy, l_nk), intent(out) :: prod

      integer :: i, j, k
      real*8 :: vector(0:l_ni+1, 0:l_nj+1, 0:l_nk+1)

      vector = 0.0d0
      prod = 0.0d0



      do k = 1, l_nk
         vector(:,:,k) = 0.0d0
         do j=1+sol_pil_s, l_nj-sol_pil_n
            do i=1+sol_pil_w, l_ni-sol_pil_e
               vector(i,j,k) = vec(i,j,k)
            enddo
         enddo
      enddo



      call rpn_comm_xch_halon (vector, 0, l_ni+1, 0, l_nj+1, l_ni, l_nj, l_nk+2, &
                               1, 1, G_periodx, G_periody, l_ni, 0, 2)



      do k=1,l_nk
         do j=1+sol_pil_s, l_nj-sol_pil_n
            do i=1+sol_pil_w, l_ni-sol_pil_e

               prod(i,j,k) = matrix(IDX_POINT,i,j,k)  * vector(i  ,j  ,k  ) + &
                             matrix(IDX_WEST,i,j,k)   * vector(i-1,j  ,k  ) + &
                             matrix(IDX_EAST,i,j,k)   * vector(i+1,j  ,k  ) + &
                             matrix(IDX_NORTH,i,j,k)  * vector(i  ,j+1,k  ) + &
                             matrix(IDX_SOUTH,i,j,k)  * vector(i  ,j-1,k  ) + &
                             matrix(IDX_TOP,i,j,k)    * vector(i  ,j  ,k-1) + &
                             matrix(IDX_BOTTOM,i,j,k) * vector(i  ,j  ,k+1)
            enddo
         enddo
      enddo


   end subroutine matvec_3d

end module matvec
