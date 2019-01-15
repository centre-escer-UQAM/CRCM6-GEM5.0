!-------------------------------------- LICENCE BEGIN ------------------------------------
!Environment Canada - Atmospheric Science and Technology License/Disclaimer,
!                     version 3; Last Modified: May 7, 2008.
!This is free but copyrighted software; you can use/redistribute/modify it under the terms
!of the Environment Canada - Atmospheric Science and Technology License/Disclaimer
!version 3 or (at your option) any later version that should be found at:
!http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html
!
!This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
!without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!See the above mentioned License/Disclaimer for more details.
!You should have received a copy of the License/Disclaimer along with this software;
!if not, you can write to: EC-RPN COMM Group, 2121 TransCanada, suite 500, Dorval (Quebec),
!CANADA, H9P 1J3; or send e-mail to service.rpn@ec.gc.ca
!-------------------------------------- LICENCE END --------------------------------------

      subroutine lightning (d,f,v,dsiz,fsiz,vsiz,ni,nk)
      use phy_options
      use phybus
      implicit none
#include <arch_specific.hf>

      integer,      intent(in)     :: dsiz,fsiz,vsiz,ni,nk
      real, target, intent(in)     :: d(dsiz)
      real, target, intent(inout)  :: f(fsiz),v(vsiz)

!Author
!          Anna Glazer   (Sept 2014)
!
!Revisions
!
!          Anna Glazer   (Avril 2015)
!          Pointers for version 4.8-a1
!          Anna Glazer   (Mai 2015)
!          Adapted to 4.8.a5
!
!Object
!          To calculate the lightning threat expressed as
!          flash-rate density (number of flashes per sec and m2)
!
!          f3 = r1*f1 + r2*f2   ! f3 = f(foudre)
!
!          ref: McCaul et al., Wea. Forecasting 2009, vol. 40, pp. 709-729
!
!Arguments
!
!          - Input/Output -
! D        dynamic bus
! F        permanent bus
!
!          - Input -
! V        volatile (output) bus
!
!          - Input -
! DSIZ     dimension of d
! FSIZ     dimension of f
! VSIZ     dimension of v
! NI       horizontal running length
! NK       vertical dimension
!
#include "thermoconsts.inc"

      integer i, k, ik
      integer index_t15(ni), array_t15(ni,nk-1), k_t15(ni)
      real t15, r1, r2, k1, k2, dsg, dpsg, invconv, tempo, i_grav
      real iiwc(ni,nk), iiwp(ni), f1(ni)
      real :: tmp1,tmp2,tmp3

      real, dimension(ni,nk) :: q_grpl, q_totalice

      real, parameter :: conv = 3.e+8    !to convert f3 unit from 5min*km2 to sec*m2
      real, parameter :: fmin = 1.e-12   !s-1 m-2, flash-rate density threshold for output
      real, parameter :: seuil = 1.e-6   !kg kg-1, mixing ratio threshold for graupel flux calculation


      real, pointer, dimension(:)   :: zfoudre, zpmoins
      real, pointer, dimension(:,:) :: zqiplus, zqgplus, zqnplus, zsigm, ztplus, zwplus, &
                                       a_qi_1, a_qi_2, a_qi_3, a_qi_4, a_qi_5, a_qi_6

      zfoudre(1:ni)      => v( foudre:)
      zpmoins(1:ni)      => f( pmoins:)
      ztplus (1:ni,1:nk) => d( tplus:)
      zsigm  (1:ni,1:nk) => d( sigm:)
      zwplus (1:ni,1:nk) => d( wplus:)

      t15 = tcdk - 15.

      r1 = 0.95
      r2 = 0.05
      k1 = 0.042
      k2 = 0.20
      invconv = 1.0/conv
      i_grav = 1.*i_grav

      if (stcond(1:5)=='MP_P3') then

        !-- diagnostic ice types (from P3)
         a_qi_1(1:ni,1:nk) => v(qi_1:)  !small ice crystals
         a_qi_2(1:ni,1:nk) => v(qi_2:)  !unrimed snow crystals
         a_qi_3(1:ni,1:nk) => v(qi_3:)  !lightly rimed snow
         a_qi_4(1:ni,1:nk) => v(qi_4:)  !graupel
         a_qi_5(1:ni,1:nk) => v(qi_5:)  !hail
         a_qi_6(1:ni,1:nk) => v(qi_6:)  !ice pellets
         q_grpl(:,:) = a_qi_4(:,:) + a_qi_5(:,:)
         iiwc(:,:)   = a_qi_1(:,:) + a_qi_1(:,:) + a_qi_1(:,:) + a_qi_4(:,:) + a_qi_5(:,:) + a_qi_6(:,:)

      elseif (stcond(1:6)=='MP_MY2') then

        !-- explicit ice types (from MY2)
         zqiplus(1:ni,1:nk) => d( qiplus:)
         zqgplus(1:ni,1:nk) => d( qgplus:)
         zqnplus(1:ni,1:nk) => d( qnplus:)
         q_grpl(:,:) = zqgplus(:,:)
         iiwc(:,:)   = zqiplus(:,:) + zqnplus(:,:) + zqgplus(:,:)

      endif

!  Graupel flux (at -15 deg C level) contribution
!  f1 = k1 * (wplus*QG); vertical velocity and graupel mixing ratio , both at -15 deg C level

!  Find level (k=k_t15(i))whose temperature most closely matches -15 deg C
!  if there is none then index_t15(i)=0 and k_t15(i)=0
!
      do i=1,ni
          index_t15(i) = 0
          k_t15(i)     = 0
          f1(i)        = 0.
          iiwp(i)      = 0.
      enddo
!
      do 15 k=1,nk-1
         do i=1,ni
         if ((ztplus(i,k) .ge. t15 .and. &
              ztplus(i,k+1) .lt. t15) .or. &
             (ztplus(i,k) .le. t15 .and. &
              ztplus(i,k+1) .gt. t15)) then

               index_t15(i)=index_t15(i)+1
               if(abs(ztplus(i,k)-t15) .lt. &
                  abs(ztplus(i,k+1)-t15)) then
                 array_t15(i,index_t15(i))=k
               else
                 array_t15(i,index_t15(i))=k+1
               endif
         endif
         enddo
 15   continue

      do i=1,ni
       if (index_t15(i) .ge. 1) then
          k_t15(i) = array_t15(i,index_t15(i))
       endif
      enddo

      do i=1,ni
       if (index_t15(i) .ge. 1 .and. &
           zwplus(i,k_t15(i)) .gt. 0. .and. &
           q_grpl(i,k_t15(i)) .gt. seuil) then
         f1(i) =1.e3*k1*zwplus(i,k_t15(i)) * q_grpl(i,k_t15(i))

       endif
      enddo

!  Contribution from solid species (mixing ratios):    iiwc = QI + QG + QS
!  f2 = k2 * iiwp (iiwp = vertical integral of total solid condensate, here iiwc)

      do i=1,ni
          dsg= 0.5 * ( zsigm(i,2) - zsigm(i,1) )
          dpsg= zpmoins(i)*dsg*i_grav
          iiwp(i) = iiwp(i) + max(iiwc(i,1) , 0. ) * dpsg
      enddo
!
     do k=2,nk-1
         do i=1,ni
            dsg= 0.5 * ( zsigm(i,k+1) - zsigm(i,k-1) )
            dpsg= zpmoins(i)*dsg*i_grav
            iiwp(i) = iiwp(i) + max(iiwc(i,k) , 0. ) * dpsg
         enddo
      enddo

      do i=1,ni
          dsg= 0.5 * ( zsigm(i,nk) - zsigm(i,nk-1) )
          dpsg= zpmoins(i)*dsg*i_grav
          iiwp(i) = iiwp(i) + max(iiwc(i,nk) , 0. ) * dpsg
      enddo

!  Combine all together and find lighting threat f3 (variable foudre)

      do i=1,ni
         tempo = r1*f1(i) + r2*k2*iiwp(i)
         if (invconv * tempo .gt. fmin) then
            zfoudre(i) = invconv * tempo
         else
            zfoudre(i) = 0.
         endif
      enddo

      return
      end

