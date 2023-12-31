!-------------------------------------- LICENCE BEGIN ------------------------
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
!-------------------------------------- LICENCE END --------------------------

!*@/
SUBROUTINE OZOREF2(O3F,LREF,DLAT,NP,NMAX,LBL,NLAT,ALAT, &
           PREF,F)
   use tdpack, only: PI
   implicit none
#include <arch_specific.hf>
      INTEGER LREF,NP,NMAX,IB,I,J,K,L,NLAT
      REAL SLOPE,XLATI,PREF(LREF)
      REAL DL,B,F1,F2,A1,A2,OZON,FAC,CONV
      REAL O3F(NMAX,LREF)
      INTEGER LBL(NP)
      REAL DLAT(NP)
      REAL F(NLAT,LREF),ALAT(NLAT)
!
!Author
!          L.Garand (1997)
!          rewritten from original CCRN code
!
!Revisions
! 001      B. Bilodeau (Jan 2000) - Exit if latitudes out of bounds
!
!Object
!          to calculate the ozone mixing ratio (kg/kg) at ozone
!          climatological levels for desired array of latitudes
!
!Arguments
!
!          - Output -
! O3F      ozone (kg O3/kg air) at  each reference
!          climatological level for DLAT latitudes
!
!          - Input -
! PREF     reference ozone pressure (N/m2) levels of field F
! LREF     number of climatological ozone level
! DLAT     latitude in radians of model points
! NP       number of points to process
! NMAX     number of maximum points permitted
! LBL      work field
! NLAT     number of latitude climatological bands
! ALAT     climatological ozone latitudes in degrees
! F        climatological ozone field in PPMV
!
!*@/

      FAC = 180./PI
!
      DO 145 J=1,NP
        IB=0
        XLATI= FLOAT( NINT(DLAT(J)*FAC) )
!
        DO 140 I=1,NLAT
          IF( XLATI.LT.ALAT(I) .AND. IB.EQ.0 ) IB=I
  140   CONTINUE
!
        IF ( XLATI.EQ. 90.0 ) IB=NLAT
!
        IF(IB.LE.1) THEN
          WRITE(6,6030) XLATI
          WRITE(6,*) (ALAT(I),i=1,NLAT)
          call physeterror('ozoref', 'O3 Inpol out bounds in latitude')
          RETURN
 6030     FORMAT(1X,' O3 INPOL OUT BOUNDS IN LATITUDE:',E12.4)
        ENDIF
        LBL(J)=IB-1
  145 CONTINUE
!
!  interpolate to desired latitudes
!  and transform into kg/kg using CONV:
!  1.E-6 converts PPMV to PPV , 48.0=M(O3), 28.964= M(dry air)
!
      CONV = 1.E-6*48./28.964

!
      DO 160 L=1,LREF
      DO 160 J=1,NP
      K=LBL(J)
      A1=ALAT(K)
      A2=ALAT(K+1)
      F1=F(K,L)
      F2=F(K+1,L)
      DL=DLAT(J)  *FAC
      SLOPE = (F2-F1)/(A2-A1)
      B = F1 - SLOPE*A1
      OZON = SLOPE*DL + B
      O3F(J,L)=OZON*CONV
  160 CONTINUE
!
!
      RETURN
      END
