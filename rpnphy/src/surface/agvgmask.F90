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
!**S/P AGGVEG
!
      SUBROUTINE AGVGMASK( FRACT, VGTYPE, AGGF, NI, NC, NCLASS, \
                           NCLASSURB, NMOS, K)
!
!
      implicit none
#include <arch_specific.hf>
!
!
!      INTEGER NI, NC, NCLASS,  VGTYPE(NCLASS)
      INTEGER NI, NC, NCLASS,NCLASSURB,  VGTYPE(NCLASS), NMOS,K
!      REAL AGGF(NI,NC), FRACT(NI,NCLASS)
      REAL AGGF(NI,NC,NMOS+1), FRACT(NI,NCLASS+NCLASSURB)
!
!
!Author
!        Yves Delage
!Revision
!
!Object
!        Aggregation of parameters for CLASS that depend on the vegetation
!        fraction masks.
!
!
!Arguments
!
!            - Input -
! FRACT      Fraction of vegetation (masks)
!
!            - Output -
! AGGF       Aggretated geophysical field representative of an entire
!            grid area
!
!            - Input -
! NI         Horizontal dimension
! NC         Number of vegetation classes in CLASS
! NCLASS     Number of landuse classes
!
!
!
      INTEGER I,J,M
      INTEGER vgtype1,vgtype2
!
!
!
      do j=1,nc
         DO i=1,ni
!           aggf(i,j) = 0.0
           aggf(i,j,k) = 0.0
         END DO
!
!
         DO i=1,ni
              DO m=4,nclass
!                if(vgtype(m).eq.j) aggf(i,j) = aggf(i,j) &
                if(vgtype(m).eq.j) aggf(i,j,k) = aggf(i,j,k) &
!                    + fract(i,m)
!  workaround for values of fract that are negative on input (LD)
                    + max(0.,fract(i,m))
!
                if (vgtype(m).gt.10 .and. vgtype(m).le.99) then
                  vgtype1 = vgtype(m) / 10
                  vgtype2 = mod(vgtype(m) ,10)
                  if(vgtype1.eq.j .or. vgtype2.eq.j) &
!                       aggf(i,j) = aggf(i,j) + fract(i,m)/2.
                       aggf(i,j,k) = aggf(i,j,k) + max(0.,fract(i,m)/2.)
                end if
!
              END DO
         END DO
      end do
!
!
!
      RETURN
      END
