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
!**S/P AGVGMASK2

      SUBROUTINE AGVGMASKC( FRACT, VGTYPE, AGGF, NI, NC, NCC, \
                            NCLASS,NCLASSURB, NMOS, K)

      implicit none
#include <arch_specific.hf>

      INTEGER NI, NC, NCC, NCLASS,NCLASSURB,  VGTYPE(NCLASS), NMOS,K
      REAL AGGF(NI,NCC,NMOS+1), FRACT(NI,NCLASS+NCLASSURB)

!Object
!        Aggregation of parameters for CTEM that depend on the vegetation
!        fraction masks.

!Arguments

!            - Input -
! FRACT      Fraction of vegetation (masks)

!            - Output -
! AGGF       Aggregated geophysical field representative of an entire
!            grid area

!            - Input -
! NI         Horizontal dimension
! NC         Number of vegetation classes in CLASS
! NCC        Number of vegetation classes in CTEM
! NCLASS     Number of landuse classes
! NCLASSURB  Number of urban landuse classes
! NMOS       Number of mosaic levels -1
! K          Mosaic level to aggregate

#include "surface.cdk"

      INTEGER I,J,M,N,IND,vgtmp

      INTEGER NOL2PFT(4)
      DATA NOL2PFT /2,3,2,2/
      INTEGER FIRSTPFT(4)
      DATA FIRSTPFT /1,3,6,8/

      REAL VGCTEM(26,3)
      DATA VGCTEM/ &
                    0.0    , 0.0    , 0.0    , 1.0    , 1.0    , &
                    0.0    , 0.0    , 0.5    , 0.0    , 0.5    , &
                    0.5    , 0.5    , 0.5    , 0.5    , 0.5    , &
                    1.0    , 0.0    , 0.0    , 0.0    , 0.5    , &
                    0.0    , 1.0    , 0.5    , 0.0    , 0.34   , &
                    0.5    , &
!  
                    0.0    , 0.0    , 0.0    , 0.0    , 0.0    , &
                    1.0    , 1.0    , 0.0    , 0.0    , 0.5    , &
                    0.5    , 0.5    , 0.5    , 0.5    , 0.5    , &
                    0.0    , 1.0    , 1.0    , 1.0    , 0.5    , &
                    0.0    , 0.0    , 0.5    , 0.0    , 0.33   , &
                    0.5    , &
!  
                    0.0    , 0.0    , 0.0    , 0.0    , 0.0    , &
                    0.0    , 0.0    , 0.5    , 1.0    , 0.0    , &
                    0.0    , 0.0    , 0.0    , 0.0    , 0.0    , &
                    0.0    , 0.0    , 0.0    , 0.0    , 0.0    , &
                    0.0    , 0.0    , 0.0    , 0.0    , 0.33   , &
                    0.0    / 
!  

      do j=1,ncc
         DO i=1,ni
           aggf(i,j,k) = 0.0
         END DO
      end do

      DO i=1,ni
         DO m=4,nclass
            if(vgtype(m).gt.0.and.vgtype(m).lt.5) then
               DO n=1,nol2pft(vgtype(m))
                  IND=firstpft(vgtype(m))+n-1
!                  aggf(i,ind,k) = aggf(i,ind,k) + fract(i,m)*VGCTEM(m,n)
!  workaround for values of fract that are negative on input (LD)
                  aggf(i,ind,k) = aggf(i,ind,k) + max(0.,fract(i,m)*VGCTEM(m,n))
               END DO
            end if 

!            if(vgtype(m).gt.10.and.vgtype(m).lt.99) then
!               vgtmp = vgtype(m)/10
!               DO n=1,nol2pft(vgtmp)
!                  IND=firstpft(vgtmp)+n-1
!                  aggf(i,ind,k) = aggf(i,ind,k) + fract(i,m)*VGCTEM(m,n)/2.0
!               END DO
!               vgtmp = mod(vgtype(m) ,10)
!               DO n=1,nol2pft(vgtmp)
!                  IND=firstpft(vgtmp)+n-1
!                  aggf(i,ind,k) = aggf(i,ind,k) + fract(i,m)*VGCTEM(m,n)/2.0
!               END DO
!            end if 
            if(vgtype(m).gt.10.and.vgtype(m).lt.99) then
               vgtmp = vgtype(m)/10
               DO n=1,nol2pft(vgtmp)
                  IND=firstpft(vgtmp)+n-1
!                  aggf(i,ind,k) = aggf(i,ind,k) + fract(i,m)*VGCTEM(m,n)/2.0
                  aggf(i,ind,k) = aggf(i,ind,k) + max(0.,fract(i,m)/(2.0*nol2pft(vgtmp)))
               END DO
               vgtmp = mod(vgtype(m) ,10)
               DO n=1,nol2pft(vgtmp)
                  IND=firstpft(vgtmp)+n-1
!                  aggf(i,ind,k) = aggf(i,ind,k) + fract(i,m)*VGCTEM(m,n)/2.0
                  aggf(i,ind,k) = aggf(i,ind,k) + max(0.,fract(i,m)/(2.0*nol2pft(vgtmp)))
               END DO
            end if 
         END DO
      END DO

      RETURN
      END
