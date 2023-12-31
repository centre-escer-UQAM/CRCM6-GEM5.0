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
      subroutine compute_re(REI, REW, &
                            lwc, iwc, aird, &
                            cloud, mg, ml, &
                            lmx, nk)
!
      implicit none
#include <arch_specific.hf>
!
#include "nbsnbl.cdk"
!
      integer lmx, m, nk, ioptrew
      real cloud(lmx,nk), mg(lmx), ml(lmx)
      real lwc(lmx,nk), iwc(lmx,nk)
      real aird(lmx,nk)
      real rew(lmx,nk), rei(lmx,nk)
      parameter(ioptrew=1)
!
!Author
!          Paul Vaillancourt / J. Li (Oct 2003)
!
!Revisions
! 001      Jason Cole (April 2006) - Modified to just compute the
!                    liquid and ice effective radii
!
!Object
!          A new scheme to calculate cloud optical properties,
!          cloud microphysical input based on cldoptx4
!
!Arguments
!          - Output -
! rei      effective radius for ice clouds
! rew      effective radius for liquid water clouds
!
!          - Input -
! LWC      in-cloud liquid water content in g/m^3
! IWC      in-cloud ice water content in  g/m^3
! CLOUD    layer cloud amount (0. or 1.) (LMX,NK)
! MG       ground cover (ocean=0.0,land <= 1.)  (LMX)
! ML       fraction of lakes (0.-1.) (LMX)
! LMX      number of profiles to compute
! NK       number of layers
!
!***********************************************************************

!...  AUTOMATIC ARRAYS
!

!
      integer i, j, k
      real rec_grav
      real cut, zrieff(lmx,nk), zrieff_vs(lmx,nk), third
      real epsilon,epsilon2,betan,betad
      real rec_cdd(lmx,nk), vs1(lmx,nk)
!
#include "cldop.cdk"
#include "tdpack_const.hf"
!
      data third/0.3333333/
      save third
      rec_grav=1./grav

!
! COMPUTE THE EFFECTIVE RADIUS FOR LIQUID AND ICE CLOUD PARTICLES
!

      DO K = 1, NK
          DO I =1, LMX
!
!...        EFFECTIVE RADIUS FOR WATER CLOUDS, SET NUMBER OF DROPS PER
!           CM^3, 100 FOR WATER AND 500 FOR LAND
!
            IF (MG(I) .LE. 0.5 .AND. ML(I) .LE. 0.5)                THEN
!              CDD=50.
              REC_CDD(I,K) =  0.01
           ELSE
!              CDD=250.
              REC_CDD(I,K) =  0.002
            ENDIF
!
!           calcul prealable au vspown de ZRIEFF
!           Units of IWC must be g/m3
!           for parameterization of REI
!
            ZRIEFF_VS(I,K) = IWC(I,K) !1000. * IWC(I,K) * AIRD(I,K)
          END DO
      END DO
!

      IF (IOPTREW .EQ. 1)                                           THEN
        DO K = 1, NK
          DO I = 1, LMX
             VS1(I,K) = (LWC(I,K)/1000.0) * REC_CDD(I,K) !* AIRD(I,K) * REC_CDD(I,K)
          ENDDO
        ENDDO
        CALL VSPOWN1(REW, VS1, THIRD, NK * LMX)
!
!...    THIS PARAMETERIZATION FROM H. BARKER, BASED ON AIRCRAFT DATA
!       RANGE 4-17 MICRON IS THAT SPECIFIED BY SLINGO FOR
!       PARAMETERIZATIONS
!
        DO K = 1, NK
          DO I = 1, LMX
            REW(I,K) = MIN ( MAX (4., 754.6 * REW(I,K)), 17.0)
          END DO
        END DO
      ELSEIF (IOPTREW .EQ. 2)                                       THEN
        DO K = 1, NK
          DO I = 1, LMX
             VS1(I,K) = (1.0 + LWC(I,K)/AIRD(I,K) * 1.E4) &
                         * (LWC(I,K)/1000.0) * REC_CDD(I,K) !* AIRD(I,K) * REC_CDD(I,K)
          ENDDO
        ENDDO
        CALL VSPOWN1(REW, VS1,THIRD, NK * LMX)
!
        DO K = 1, NK
          DO I = 1, LMX
            REW(I,K) =  3000. * REW(I,K)
            REW(I,K) =  MIN (MAX (2.5, REW(I,K)), 50.0)
          END DO
        END DO
      ELSEIF (IOPTREW .EQ. 3)                                       THEN
       DO K = 1, NK
          DO I = 1, LMX
            EPSILON =  1.0 - 0.7 * EXP(- 0.001 / REC_CDD(I,K))
            EPSILON2 =  EPSILON * EPSILON
            BETAD =  1.0 + EPSILON2
            BETAN =  BETAD + EPSILON2
!            REW(I,K) = 620.3504944*((BETAN*BETAN*LWC(i,k)*aird(i,k))
!     1                  / (BETAD / REC_CDD(I,K)) )**THIRD
            REW(I,K) = 620.3504944*((BETAN*BETAN*LWC(i,k)) &
                        / (BETAD / REC_CDD(I,K)) )**THIRD
            REW(I,K) =  MIN (MAX (2.5, REW(I,K)), 17.0)
!            REW(I,K) =  MIN (MAX (4.0, REW(I,K)), 17.0)
          END DO
        END DO
      ENDIF
!
!...        EFFECTIVE RADIUS FOR ICE CLOUDS
!

        CALL VSPOWN1(ZRIEFF, ZRIEFF_VS, 0.216, NK * LMX)
        DO K = 1, NK
          DO I = 1, LMX
            IF (IWC(I,K) .GE. 1.E-9)                              THEN
               ZRIEFF(I,K) = 83.8 * ZRIEFF(I,K)
            ELSE
               ZRIEFF(I,K) = 20.
            ENDIF
            REI(I,K) =  MAX (MIN (ZRIEFF(I,K), 50.0), 20.0)
         END DO
      END DO


      RETURN
      END
