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

subroutine inicover2(kount, ni, trnch)
   use mu_jdate_mod, only: jdate_day_of_year
   use sfc_options
   use sfcbus_mod
   use class_configs
   implicit none
#include <arch_specific.hf>
#include <rmnlib_basics.hf>
   include "sfcinput.cdk"

   integer ni, kount, trnch

   !@Author Bernard Bilodeau and Stephane Belair (May 2000)
   !@Revision
   ! 001      see version 5.5.0 for previous history
   ! 002      K. Winger (UQAM/ESCER) Jun 2020 - Section for CLASS added
   !@Object Initialize vegetation fields for the surface schemes
   !@Arguments
   !       - Input -
   ! kount    current timestep number
   ! ni       horizontal slice dimension
   !
   !@Notes    inisurf has been split in two subroutines:
   !          inisurf and inicover. the former calls the latter.
   !
   !     the geophysical fields determined from vegetation
   !     are done so using the following classification:
   !
   !     Class       Vegetation type
   !     =====       ===============
   !       1         (salt) water
   !       2         ice
   !       3         inland lake
   !       4         evergreen needleleaf trees
   !       5         evergreen broadleaf trees
   !       6         deciduous needleleaf trees
   !       7         deciduous broadleaf trees
   !       8         tropical broadleaf trees
   !       9         drought deciduous trees
   !       10        evergreen broadleaf shrub
   !       11        deciduous shrubs
   !       12        thorn shrubs
   !       13        short grass and forbs
   !       14        long grass
   !       15        crops
   !       16        rice
   !       17        sugar
   !       18        maize
   !       19        cotton
   !       20        irrigated crops
   !       21        urban
   !       22        tundra
   !       23        swamp
   !       24        desert
   !       25        mixed wood forests
   !       26        mixed shrubs

   !********************************************************************
   ! Tables for the veg characteristics for each veg type
   !********************************************************************

   real aldat(nclass), d2dat(nclass), rsminxdat(nclass)
   real laidat(nclass), vegdat(nclass)
   real cvdat(nclass), rgldat(nclass), gammadat(nclass)
   real alvsdat(nclass), alnidat(nclass),rsmindat(nclass)
   real qa50dat(nclass),vpdadat(nclass),vpdbdat(nclass)
   real psgadat(nclass),psgbdat(nclass),z0mdat(nclass),ln_z0mdat(nclass)
   real laimxdat(nclass),laimndat(nclass),vgmasdat(nclass)
   real rootdat(nclass),fveg(4)
   integer vgclass(nclass),vg000(nclass)

   integer :: nmos
   
   data d2dat/ &
        1.0    , 1.0    , 1.0    , 3.0    , 3.0    , &
        1.0    , 3.0    , 5.0    , 5.0    , 2.0    , &
        2.0    , 2.0    , 1.5    , 2.0    , 2.0    , &
        1.2    , 1.0    , 1.5    , 2.0    , 1.5    , &
        1.0    , 1.0    , 2.0    , 1.0    , 2.0    , &
        2.0    /

   data rsminxdat/ &
        500.   , 500.   , 500.   , 250.   , 250.   , &
        250.   , 250.   , 250.   , 250.   , 150.   , &
        150.   , 150.   ,  40.   ,  40.   ,  40.   , &
        40.   ,  40.   ,  40.   ,  40.   , 150.   , &
        150.   , 150.   , 150.   , 500.   , 250.   , &
        150.   /
   data laidat/ &
        0.00   , 0.00   , 0.00   , 5.00   , 6.00   , &
        -99.    , -99.   , 6.00   , 4.00   , 3.00   , &
        -99.    , 3.00   , 1.00   , -99.   , -99.   , &
        -99.    , -99.   , -99.   , -99.   , 1.00   , &
        1.00   , -99.   , 4.00   , 0.00   , -99.   , &
        -99.    /
   data vegdat/ &
        0.00   , 0.00   , 0.00   , 0.90   , 0.99   , &
        0.90   , 0.90   , 0.99   , 0.90   , 0.50   , &
        0.50   , 0.50   , 0.85   , 0.30   , -99.   , &
        -99.   , -99.   , -99.   , -99.   , 0.85   , &
        0.10   , 0.50   , 0.60   , 0.00   , 0.90   , &
        0.90   /
   data cvdat/ &
        2.0E-5 , 2.0E-5 , 2.0E-5 , 1.0E-5 , 1.0E-5 , &
        1.0E-5 , 1.0E-5 , 1.0E-5 , 1.0E-5 , 2.0E-5 , &
        2.0E-5 , 2.0E-5 , 2.0E-5 , 2.0E-5 , 2.0E-5 , &
        2.0E-5 , 2.0E-5 , 2.0E-5 , 2.0E-5 , 2.0E-5 , &
        2.0E-5 , 2.0E-5 , 2.0E-5 , 2.0E-5 , 2.0E-5 , &
        2.0E-5 /
   data rgldat/ &
        100.   , 100.   , 100.   , 30.    , 30.    , &
        30.    , 30.    , 30.    , 30.    , 100.   , &
        100.   , 100.   , 100.   , 100.   , 100.   , &
        100.   , 100.   , 100.   , 100.   , 100.   , &
        100.   , 100.   , 100.   , 100.   , 100.   , &
        100.   /
   data gammadat/ &
        0.    , 0.     , 0.     , 0.04   , 0.04   , &
        0.04  , 0.04   , 0.04   , 0.04   , 0.     , &
        0.    , 0.     , 0.     , 0.     , 0.     , &
        0.    , 0.     , 0.     , 0.     , 0.     , &
        0.    , 0.     , 0.     , 0.     , 0.     , &
        0.    /

   data rsmindat/ &
        0.0    , 0.0    , 0.0    , 250.   , 250.   , &
        263.   , 130.   , 130.   , 122.   , 323.   , &
        855.   , 500.   , 150.   , 150.   , 100.   , &
        120.   , 278.   , 90.    , 112.   , 86.    , &
        0.0    , 200.   , 200.   , 0.0    , 165.   , &
        855.   /

   data qa50dat/ &
        30.    , 30.    , 30.    , 30.    , 30.    , &
        30.    , 50.    , 30.    , 30.    , 30.    , &
        30.    , 30.    , 35.    , 30.    , 30.    , &
        30.    , 30.    , 30.    , 30.    , 30.    , &
        30.    , 50.    , 30.    , 30.    , 30.    , &
        30.    /

   data vpdadat/ &
        0.0    , 0.0    , 0.0    , 0.57   , 0.5    , &
        0.5    , 0.60   , 0.45   , 0.5    , 0.5    , &
        0.5    , 0.5    , 0.5    , 0.5    , 0.5    , &
        0.5    , 0.5    , 0.5    , 0.5    , 0.5    , &
        0.5    , 0.62   , 0.5    , 0.5    , 0.40   , &
        0.5    /

   data vpdbdat/ &
        1.0    , 1.0    , 1.0    , 1.0    , 1.0    , &
        1.0    , 0.5    , 0.0    , 1.0    , 1.0    , &
        1.0    , 1.0    , 1.0    , 1.0    , 1.0    , &
        1.0    , 1.0    , 1.0    , 1.0    , 1.0    , &
        1.0    , 0.4    , 1.0    , 1.0    , 0.6    , &
        1.0    /

   data psgadat/ &
        100    , 100    , 100    , 100    , 100    , &
        100    , 100    , 100    , 100    , 100    , &
        100    , 100    , 100    , 100    , 100    , &
        100    , 100    , 100    , 100    , 100    , &
        100    , 100    , 100    , 100    , 100    , &
        100    /

   data psgbdat/ &
        5.    , 5.     , 5.     , 5.     , 5.     , &
        5.    , 5.     , 5.     , 5.     , 5.     , &
        5.    , 5.     , 5.     , 5.     , 5.     , &
        5.    , 5.     , 5.     , 5.     , 5.     , &
        5.    , 5.     , 5.     , 5.     , 5.     , &
        5.    /

   data aldat/ &
        0.13   , 0.70   , 0.13   , 0.14   , 0.12   , &
        0.14   , 0.18   , 0.13   , 0.17   , 0.14   , &
        0.18   , 0.19   , 0.20   , 0.19   , 0.20   , &
        0.21   , 0.18   , 0.18   , 0.25   , 0.18   , &
        0.12   , 0.17   , 0.12   , 0.30   , 0.15   , &
        0.15   /

   data alvsdat/ &
        0.0    , 0.0    , 0.0    , 0.03   , 0.03   , &
        0.03   , 0.05   , 0.03   , 0.05   , 0.03   , &
        0.05   , 0.06   , 0.06   , 0.05   , 0.06   , &
        0.06   , 0.05   , 0.05   , 0.07   , 0.06   , &
        0.09   , 0.05   , 0.03   , 0.30   , 0.04   , &
        0.04   /

   data alnidat/ &
        0.0    , 0.0    , 0.0    , 0.19   , 0.23   , &
        0.19   , 0.29   , 0.23   , 0.29   , 0.19   , &
        0.29   , 0.32   , 0.34   , 0.31   , 0.34   , &
        0.36   , 0.31   , 0.31   , 0.43   , 0.36   , &
        0.15   , 0.29   , 0.25   , 0.30   , 0.26   , &
        0.26   /

   data z0mdat / &
        0.001  , 0.001  , 0.001  , 1.5    , 3.5    , &
        1.0    , 2.0    , 3.0    , 0.8    , 0.05   , &
        0.15   , 0.15   , 0.02   , 0.08   , 0.08   , &
        0.08   , 0.35   , 0.25   , 0.10   , 0.08   , &
        1.35   , 0.01   , 0.05   , 0.05   , 1.5    , &
        0.05   /

   DATA ln_z0mdat / &
        -6.91  , -6.91  , -6.91  ,  0.405 ,  1.25  , &
         0.0   ,  0.693 ,  1.10  , -0.223 , -3.0   , &
        -1.9   , -1.9   , -3.91  , -2.53  , -2.53  , &
        -2.53  , -1.05  , -1.39  , -2.30  , -2.53  , &
         0.3   , -4.61  , -3.0   , -3.0   ,  0.405 , &
        -3.0   /

   data laimxdat/ &
        0.0    , 0.0    , 0.0    , 2.0    , 10.    , &
        2.0    , 6.0    , 10.    , 4.0    , 2.0    , &
        4.0    , 3.0    , 3.0    , 4.0    , 4.0    , &
        6.5    , 5.0    , 4.0    , 5.0    , 4.0    , &
        0.0    , 1.5    , 1.5    , 0.0    , 5.5    , &
        3.0    /

   data laimndat/ &
        0.0    , 0.0    , 0.0    , 1.6    , 10.    , &
        0.5    , 0.5    , 10.    , 4.0    , 2.0    , &
        0.5    , 3.0    , 3.0    , 4.0    , 0.0    , &
        0.0    , 0.0    , 0.0    , 0.0    , 0.0    , &
        0.0    , 1.5    , 1.5    , 0.0    , 1.0    , &
        3.0    /

   data vgmasdat/ &
        0.0    , 0.0    , 0.0    , 25.    , 50.    , &
        15.    , 20.    , 40.    , 15.    , 2.     , &
        8.     , 8.     , 1.5    , 3.     , 2.     , &
        2.     , 5.     , 5.     , 2.     , 2.     , &
        0.     , 0.2    , 1.0    , 0.     , 20.    , &
        8.     /

   data rootdat/ &
        0.0    , 0.0    , 0.0    , 1.0    , 5.0    , &
        1.0    , 2.0    , 5.0    , 5.0    , 0.2    , &
        1.0    , 5.0    , 1.2    , 1.2    , 1.2    , &
        1.2    , 1.0    , 1.5    , 2.0    , 5.0    , &
        0.     , 0.1    , 5.0    , 0.     , 1.2    , &
        1.2    /

! KW:
! Move "mixed wood forests (25)" from "broadleaf (2)" to "needle- and broadleaf (12)"
! Move "deciduous shrubs (11)" from "grass (4) to "broadleaf (2)"
! Move "desert (24)" from "urbain (5)" to "bare soil (6)"
   data vgclass/ &
        0      , 0      , 0      , 1      , 2      , &
        1      , 2      , 2      , 2      , 4      , &
        2      , 4      , 4      , 4      , 3      , &
        3      , 3      , 3      , 3      , 3      , &
        5      , 4      , 4      , 6      , 12     , &
        4      /

   data vg000 / &
        1      , 1      , 1      , 1      , 1      , &
        1      , 1      , 1      , 1      , 1      , &
        1      , 1      , 1      , 1      , 1      , &
        1      , 1      , 1      , 1      , 1      , &
        1      , 1      , 1      , 1      , 1      , &
        1      /

! All values set to 1.00 so that no bare soil is added artificially (KW)
! Although, this parameter does not seem to get used anymore (KW)
!   data fveg/ 0.90,  0.90,   0.70,  0.60 /
   data fveg/ 1.00,  1.00,   1.00,  1.00 /

   !********************************************************************
   !                tables describing the annual evolution of veg fields
   !********************************************************************

   real, save :: vegcrops(13)

   data vegcrops/ &
        0.05   , 0.05   , 0.05   , 0.10   , 0.20   , &
        0.40   , 0.80   , 0.80   , 0.90   , 0.05   , &
        0.05   , 0.05   , 0.05                      /

   real, save :: lai6(13), lai7(13), lai11(13), lai14(13), lai15(13), &
        lai16(13), lai17(13), lai18(13), lai19(13), lai22(13), &
        lai25(13), lai26(13)

   data lai6 / &
        0.1   , 0.1   , 0.5   , 1.0   , 2.0   , &
        4.0   , 5.0   , 5.0   , 4.0   , 2.0   , &
        1.0   , 0.1   , 0.1                      /
   data lai7 / &
        0.1   , 0.1   , 0.5   , 1.0   , 2.0   , &
        4.0   , 5.0   , 5.0   , 4.0   , 2.0   , &
        1.0   , 0.1   , 0.1                      /
   data lai11/ &
        0.5   , 0.5   , 1.0   , 1.0   , 1.5   , &
        2.0   , 3.0   , 3.0   , 2.0   , 1.5   , &
        1.0   , 0.5   , 0.5                      /
   data lai14/ &
        0.5   , 0.5   , 0.5   , 0.5   , 0.5   , &
        0.5   , 1.0   , 2.0   , 2.0   , 1.5   , &
        1.0   , 1.0   , 0.5                      /
   data lai15/ &
        0.1   , 0.1   , 0.1   , 0.5   , 1.0   , &
        2.0   , 3.0   , 3.5   , 4.0   , 0.1   , &
        0.1   , 0.1   , 0.1                      /
   data lai16/ &
        0.1   , 0.1   , 0.1   , 0.5   , 1.0   , &
        2.5   , 4.0   , 5.0   , 6.0   , 0.1   , &
        0.1   , 0.1   , 0.1                      /
   data lai17/ &
        0.1   , 0.1   , 0.1   , 0.5   , 1.0   , &
        3.0   , 4.0   , 4.5   , 5.0   , 0.1   , &
        0.1   , 0.1   , 0.1                      /
   data lai18/ &
        0.1   , 0.1   , 0.1   , 0.5   , 1.0   , &
        2.0   , 3.0   , 3.5   , 4.0   , 0.1   , &
        0.1   , 0.1   , 0.1                      /
   data lai19/ &
        0.1   , 0.1   , 0.1   , 0.5   , 1.0   , &
        3.0   , 4.0   , 4.5   , 5.0   , 0.1   , &
        0.1   , 0.1   , 0.1                      /
   data lai22/ &
        1.0   , 1.0   , 0.5   , 0.1   , 0.1   , &
        0.1   , 0.1   , 1.0   , 2.0   , 1.5   , &
        1.5   , 1.0   , 1.0                      /
   data lai25/ &
        3.0   , 3.0   , 3.0   , 4.0   , 4.5   , &
        5.0   , 5.0   , 5.0   , 4.0   , 3.0   , &
        3.0   , 3.0   , 3.0                      /
   data lai26/ &
        3.0   , 3.0   , 3.0   , 4.0   , 4.5   , &
        5.0   , 5.0   , 5.0   , 4.0   , 3.0   , &
        3.0   , 3.0   , 3.0                      /

   !********************************************************************

   integer(IDOUBLE), parameter :: MU_JDATE_HALFDAY = 43200 !#TODO: use value from my_jdate_mod
   real, external :: interpveg

   integer :: i,j,k
   real :: julien, juliens
   real :: fcansum

   real, dimension(nclass) :: aldatd, cvdatd, d2datd, gammadatd, &
        laidatdn, laidatds, rgldatd, rsmindatd, &
        vegdatdn, vegdatds
   real, pointer, dimension (:,:) :: zfcanmx, zfcancmx, zvegf


   IF_ISBA: if (schmsol == 'ISBA') then

      ! Determine the current julian day
      julien = real(jdate_day_of_year(jdateo + kount*int(delt) + MU_JDATE_HALFDAY))

      ! Do the aggregation
      do i=1,nclass
         aldatd(i)    = aldat(i)
         d2datd(i)    = d2dat(i)
         rsmindatd(i) = rsminxdat(i)
         laidatdn(i)  = laidat(i)
         laidatds(i)  = laidat(i)
         vegdatdn(i)  = vegdat(i)
         vegdatds(i)  = vegdat(i)
         cvdatd(i)    = cvdat(i)
         rgldatd(i)   = rgldat(i)
         gammadatd(i) = gammadat(i)
      end do

      ! Fill the laidatd and vegdatd fields for
      ! land use classes varying with seasons
      ! (i.e., replace the -99 values in the table
      ! with temporal interpolations from the tables above)

      ! tables for northern hemisphere

      laidatdn( 6)  = interpveg(julien , lai6 )
      laidatdn( 7)  = interpveg(julien , lai7 )
      laidatdn(11)  = interpveg(julien , lai11)
      laidatdn(14)  = interpveg(julien , lai14)
      laidatdn(15)  = interpveg(julien , lai15)
      laidatdn(16)  = interpveg(julien , lai16)
      laidatdn(17)  = interpveg(julien , lai17)
      laidatdn(18)  = interpveg(julien , lai18)
      laidatdn(19)  = interpveg(julien , lai19)
      laidatdn(22)  = interpveg(julien , lai22)
      laidatdn(25)  = interpveg(julien , lai25)
      laidatdn(26)  = interpveg(julien , lai26)

      vegdatdn(15)  = interpveg(julien , vegcrops)
      vegdatdn(16)  = interpveg(julien , vegcrops)
      vegdatdn(17)  = interpveg(julien , vegcrops)
      vegdatdn(18)  = interpveg(julien , vegcrops)
      vegdatdn(19)  = interpveg(julien , vegcrops)

      !  tables for southern hermisphere
      juliens = julien  - 183
      if (juliens < 0.) juliens = juliens + 366.

      laidatds( 6)  = interpveg(juliens, lai6 )
      laidatds( 7)  = interpveg(juliens, lai7 )
      laidatds(11)  = interpveg(juliens, lai11)
      laidatds(14)  = interpveg(juliens, lai14)
      laidatds(15)  = interpveg(juliens, lai15)
      laidatds(16)  = interpveg(juliens, lai16)
      laidatds(17)  = interpveg(juliens, lai17)
      laidatds(18)  = interpveg(juliens, lai18)
      laidatds(19)  = interpveg(juliens, lai19)
      laidatds(22)  = interpveg(juliens, lai22)
      laidatds(25)  = interpveg(juliens, lai25)
      laidatds(26)  = interpveg(juliens, lai26)

      vegdatds(15)  = interpveg(juliens, vegcrops)
      vegdatds(16)  = interpveg(juliens, vegcrops)
      vegdatds(17)  = interpveg(juliens, vegcrops)
      vegdatds(18)  = interpveg(juliens, vegcrops)
      vegdatds(19)  = interpveg(juliens, vegcrops)

#define PTR1D(NAME2) busptr(vd%NAME2%i)%ptr(1,trnch)

      call aggcovernat(PTR1D(vegf), laidatdn, laidatds, PTR1D(lai), &
           PTR1D(dlat), ni, nclass)
      call aggcovernat(PTR1D(vegf), vegdatdn, vegdatds, PTR1D(vegfrac), &
           PTR1D(dlat), ni, nclass)
      call aggcovernat(PTR1D(vegf), aldatd, aldatd, PTR1D(alveg), &
           PTR1D(dlat), ni, nclass)
      call aggcovernat(PTR1D(vegf), d2datd, d2datd , PTR1D(rootdp), &
           PTR1D(dlat), ni, nclass)
      call aggcovernat(PTR1D(vegf), rsmindatd, rsmindatd, PTR1D(stomr), &
           PTR1D(dlat), ni, nclass)
      call aggcovernat(PTR1D(vegf), cvdatd, cvdatd, PTR1D(cveg), &
           PTR1D(dlat), ni, nclass)
      call aggcovernat(PTR1D(vegf), rgldatd, rgldatd, PTR1D(rgl), &
           PTR1D(dlat), ni, nclass)
      call aggcovernat(PTR1D(vegf), gammadatd , gammadatd, PTR1D(gamveg), &
           PTR1D(dlat), ni, nclass)

   endif IF_ISBA


   IF_CLASS: if (schmsol == 'CLASS') then
       nmos = 0
       call agvgclas(PTR1D(vegf),ALVSDAT, VGCLASS,PTR1D(ALVSC), NI,class_ic+1, &
                     NCLASS,0,nmos,1)
       call agvgclas(PTR1D(vegf),ALNIDAT, VGCLASS,PTR1D(ALIRC), NI,class_ic+1, &
                     NCLASS,0,nmos,1)
       call agvgclas(PTR1D(vegf),LAIMXDAT,VGCLASS,PTR1D(LAIMAX),NI,class_ic  , &
                     NCLASS,0,nmos,1)
       call agvgclas(PTR1D(vegf),LAIMNDAT,VGCLASS,PTR1D(LAIMIN),NI,class_ic  , &
                     NCLASS,0,nmos,1)
       call agvgclas(PTR1D(vegf),VGMASDAT,VGCLASS,PTR1D(VEGMA), NI,class_ic  , &
                     NCLASS,0,nmos,1)
       call agvgclas(PTR1D(vegf),ROOTDAT, VGCLASS,PTR1D(ROOTDP),NI,class_ic  , &
                     NCLASS,0,nmos,1)
       call agvgclas(PTR1D(vegf),LN_Z0MDAT,VGCLASS,PTR1D(ZOLN), NI,class_ic+1, &
                     NCLASS,0,nmos,1)
       call agvgclas(PTR1D(vegf),RSMINDAT,VGCLASS,PTR1D(STOMR), NI,class_ic  , &
                     NCLASS,0,nmos,1)
       call agvgclas(PTR1D(vegf),QA50DAT, VGCLASS,PTR1D(QA50),  NI,class_ic  , &
                     NCLASS,0,nmos,1)
       call agvgclas(PTR1D(vegf),VPDADAT, VGCLASS,PTR1D(VPDA),  NI,class_ic  , &
                     NCLASS,0,nmos,1)
       call agvgclas(PTR1D(vegf),VPDBDAT, VGCLASS,PTR1D(VPDB),  NI,class_ic  , &
                     NCLASS,0,nmos,1)
       call agvgclas(PTR1D(vegf),PSGADAT, VGCLASS,PTR1D(PSIGA), NI,class_ic  , &
                     NCLASS,0,nmos,1)
       call agvgclas(PTR1D(vegf),PSGBDAT, VGCLASS,PTR1D(PSIGB), NI,class_ic  , &
                     NCLASS,0,nmos,1)
       call agvgclas(PTR1D(vegf),ROOTDAT, VG000  ,PTR1D(SDEPTH),NI,1         , &
                     NCLASS,0,nmos,1)
       call agvgmask(PTR1D(vegf),         VGCLASS,PTR1D(FCANMX),NI,class_ic+1, &
                     NCLASS,0,nmos,1)
       IF (ctem_mode.gt.0) &
       call agvgmaskc(PTR1D(vegf), VGCLASS,PTR1D(FCANCMX), NI, class_ic, &
                      ctem_ICC, NCLASS,0,nmos,1)

#define MKPTR2D(NAME1,NAME2) nullify(NAME1); if (vd%NAME2%i > 0 .and.  associated(busptr(vd%NAME2%i)%ptr)) NAME1(1:ni,1:vd%NAME2%mul*vd%NAME2%niveaux) => busptr(vd%NAME2%i)%ptr(:,trnch)

   MKPTR2D(zfcancmx,fcancmx)
   MKPTR2D(zfcanmx,fcanmx)
   MKPTR2D(zvegf,vegf)

!
!       Normalize FCANMX
        DO i=1,ni
!          IF (vege_fields.ne.'CTEM') then
          IF (.not.any('fcancmx'==phyinread_list_s(1:phyinread_n))) then

!         Sum up CLASS vegetation fractions
          fcansum = 0.
!         Do not count too small fractions - their other fields were not initialized!
          DO J=1,class_ic+1
            if (ZFCANMX(i,j).lt.critmask) then
              IF (ctem_mode.gt.0.and.j.lt.class_ic+1) then
                DO k=1,nol2pft(j)
                  ZFCANCMX(i,(firstpft(j)+k-1))=0.
                enddo
              end if
              ZFCANMX(i,j)=0.
            end if
            fcansum = fcansum + ZFCANMX(i,j)
          ENDDO
!         Add bare soil (desert, VF(24))
!          fcansum = fcansum + zvegf(i,24)
          fcansum = fcansum + max(0.,zvegf(i,24))
!         Normalize
          if ( fcansum .ge. critmask ) then
            DO J=1,class_ic+1
              ZFCANMX(i,j) = ZFCANMX(i,j) / fcansum
            ENDDO
            IF (ctem_mode.gt.0) then
              DO J=1,ctem_icc
                ZFCANCMX(i,j) = ZFCANCMX(i,j) / fcansum
              ENDDO
            end if
          end if

          end if
        ENDDO
   endif IF_CLASS

   return
end subroutine inicover2
