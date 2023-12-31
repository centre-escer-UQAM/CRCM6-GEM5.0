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
!**S/P ISCCP_SIM_DRIVER - PART OF THE ISCCP CLOUD SIMULATOR PACKAGE
!
      SUBROUTINE ISCCP_SIM_DRIVER(CLDPT, CLDTP, CLDTAU, CLDEP, CLDFRAC,   &! OUTPUT
                                  SUNLIT, &
                                  CLW_SUB, CIC_SUB, PRESSG, SHJ, SHTJ,    &! INPUT
                                  IL1, IL2, ILG, LAY, LEV, &
                                  COSZ, GT, T, Q, MG, ML)


      implicit none
#include <arch_specific.hf>
!
#include "nbsnbl.cdk"
#include "tdpack_const.hf"
!
#include "mcica.cdk"

      EXTERNAL compute_re
!
!Author
!
!
!        Jason Cole Dec. 16, 2005
!
!Revisions
!
! 001
!
!Object
!
!---------------------------------------------------------------------------------
! This code generates is a slightly modified version of the ISCCP simulator.
! Rather than use the cloud generator provided with the simulator this version
! uses as input subcolumns generated by the stochastic cloud generator of Barker
! and Rasiasen.  I.e., the ISCCP simulator works on the same clouds as the radiation
! in the Monte Carlo Independent Column Approximation.
!---------------------------------------------------------------------------------

!
!
! INPUT DATA
!

      REAL, INTENT(IN) :: &
       CLW_SUB(ILG,LAY,NX_LOC), &
       CIC_SUB(ILG,LAY,NX_LOC), &
       PRESSG(ILG), &
       SHJ(ILG,LAY), &
       SHTJ(ILG,LEV), &
       COSZ(ILG), &
       GT(ILG), &
       T(ILG,LAY), &
       Q(ILG,LAY), &
       MG(ILG), &
       ML(ILG)

      INTEGER, INTENT(IN) :: &
       IL1, &
       IL2, &
       ILG, &
       LAY, &
       LEV

!
! OUTPUT DATA
!

      REAL, INTENT(OUT) :: &
       CLDPT(ILG,NTAU*NPTOP), &
       CLDTP(ILG), &
       CLDTAU(ILG), &
       CLDEP(ILG), &
       CLDFRAC(ILG), &
       SUNLIT(ILG)

!
! WORK DATA
!

      REAL :: &
       CLDMASK(ILG,LAY), &
       TAU_VIS(ILG,LAY), &
       EMS_CLD(ILG,LAY), &
       PFULL(ILG,LAY), &
       CIC(ILG,LAY), &
       CLW(ILG,LAY), &
       DZ(ILG,LAY), &
       PHALF(ILG,LEV), &
       EMS_SFC(ILG), &
       TAU_SUB(ILG), &
       LOGTAU_SUB(ILG), &
       PTOP_SUB(ILG), &
       RAN_NUM(ILG), &
       REI(ILG,LAY), &
       REL(ILG,LAY), &
       AIRD(ILG,LAY)

      INTEGER :: &
       RAN_INDEX(ILG), &
       SUNLIT_LOC(ILG)

!
! LOCAL DATA (SCALAR)
!

      INTEGER :: &
       ICOL, &
       IND1, &
       IL, &
       ILAY, &
       ILEV, &
       ILAYP1, &
       ITAU, &
       IPRES, &
       IT, &
       IP, &
       ITP

      REAL :: &
       INVREL, &
       INVREI, &
       ROG, &
       TAULIQVIS, &
       TAULIQIR, &
       TAUICEVIS, &
       TAUICEIR, &
       TAU_IR, &
       TEMP, &
       A, &
       B, &
       C, &
       DP1, &
       DP2, &
       DP3, &
       DP_LOC

!
! PARAMETERS
!

      REAL, PARAMETER :: COSZ_MIN = 0.2          ! CUTOFF FOR DAYLIGHT PER ISCCP
      REAL, PARAMETER :: THIRD    = 0.3333333333333
      REAL, PARAMETER :: RU       = 1.6487213

! LOWER BOUNDS CLOUD TAU AND CLOUD TOP PRESSURE BINS
      REAL, PARAMETER :: P_TOP_BNDS(NPTOP) = (/0.0  ,180.0,310.0,440.0, &
                                               560.0, 680.0, 800.0/)
      REAL, PARAMETER :: TAU_BNDS(NTAU) = (/0.0 ,0.3 ,1.3,3.6,9.4, &
                                            23.0,60.0/)

      REAL, PARAMETER :: TAUCHK = 0.000001 !-1.*log(0.9999999)
!
! ZERO OUT SOME ARRAYS
!

      DO ILAY = 1, LAY
         DO IL = IL1, IL2
            TAU_VIS(IL,ILAY) = 0.0
         END DO ! IL
      END DO ! ILAY

      DO IL = IL1, IL2
         CLDTP(IL)   = 0.0
         CLDTAU(IL)  = 0.0
         CLDEP(IL)   = 0.0
         CLDFRAC(IL) = 0.0
         SUNLIT(IL)  = 0.0
      END DO ! IL

      DO ITP = 1, NPTOP*NTAU
         DO IL = IL1, IL2
            CLDPT(IL,ITP) = 0.0
         END DO ! IL
      END DO ! ITP

!
!----------------------------------------------------------------------C
!     COMPUTE THE AIR DENSITY IN KG/M^3 AND THE LAYER THICKNESS IN M   C
!----------------------------------------------------------------------C
!

      DO ILAY = 1, LAY
         DO IL = IL1, IL2
            AIRD(IL,ILAY) =  SHJ(IL,ILAY) * PRESSG(IL) &
                          / ( T(IL,ILAY) * RGASD )
            DP1=0.5*(SHJ(IL,MIN(ILAY+1,LAY))-SHJ(IL,MAX(ILAY-1,1)))
            DP2=0.5*(SHJ(IL,1)+SHJ(IL,2))
            DP3=0.5*(1.-SHJ(IL,LAY))
            IF (ILAY .EQ. 1) THEN
               DP_LOC = DP2
            ELSE IF (ILAY .EQ. LAY) THEN
               DP_LOC = DP3
            ELSE
               DP_LOC = DP1
            ENDIF

            DP_LOC= MAX(DP_LOC*PRESSG(IL),0.)
            DZ(IL,ILAY)= DP_LOC/(AIRD(IL,ILAY)*GRAV)
         END DO ! IL
      END DO ! ILAY

! SET SURFACE EMISSIVITY TO 0.9999999
! SHOULD BE PASSED IN WHEN SURFACE EMISSIVITY IS COMPUTED BY
! GCM

      DO IL = IL1, IL2
         EMS_SFC(IL) = 0.9999999
      END DO ! IL

! SUNLIT ARRAY
      DO IL = IL1, IL2
         IF (COSZ(IL) .GE. COSZ_MIN) THEN
            SUNLIT_LOC(IL) = 1
         ELSE
            SUNLIT_LOC(IL) = 0
         END IF
      END DO

! COMPUTE PRESSURES

      DO IL = IL1, IL2
         PHALF(IL,LEV) = SHTJ(IL,LEV)*PRESSG(IL)
      END DO ! IL

      DO ILAY = 1, LAY
         DO IL = IL1, IL2
            PHALF(IL,ILAY) = SHTJ(IL,ILAY)*PRESSG(IL)
            PFULL(IL,ILAY) = SHJ(IL,ILAY)*PRESSG(IL)
         END DO ! IL
      END DO ! ILAY

! LOOP OVER THE SUBCOLUMNS TO SAMPLE FROM, RANDOMLY SELECT SUBCOLUMNS
! AND APPLY THE ISCCP SIMULATOR ALGORIGHTM TO THE SUBCOLUMNS.
! THE RESULTS FOR EACH SUBCOLUMN ARE THEN ACCUMULATED TO GIVE A
! GCM COLUMN MEAN.

      DO ICOL = 1, NSUBCOL

! RANDOMLY SELECT A SUBCOLUMN

         CALL RANDOM_NUMBER(RAN_NUM)
         DO Il = IL1, IL2
            IND1 = INT(RAN_NUM(IL)*REAL(NX_LOC))
            IF (IND1 .GT. NX_LOC) IND1 = NX_LOC
            IF (IND1 .LT. 1)      IND1 = 1
            RAN_INDEX(IL) = IND1
         END DO ! IL
         DO IL = IL1, IL2
            IND1 = RAN_INDEX(IL)
            DO ILAY = 1, LAY
               CLW(IL,ILAY) = CLW_SUB(IL,ILAY,IND1)
               CIC(IL,ILAY) = CIC_SUB(IL,ILAY,IND1)
               IF ((CLW(IL,ILAY)+CIC(IL,ILAY)) .GT. 1.0e-9) THEN
                  CLDMASK(IL,ILAY) = 1.0
               ELSE
                  CLDMASK(IL,ILAY) = 0.0
               END IF
            END DO ! ILAY
         END DO ! IL

! COMPUTE THE CLOUD OPTICAL THICKNESS AT ~0.6 AND ~10.5 MICRONS
! THESE VALUES ARE FROM J. LI (JAN. 22, 2006)

! COMPUTE THE CLOUD EFFECTIVE RADII

         CALL COMPUTE_RE(REI, REL, &
                         CLW, CIC, AIRD, &
                         CLDMASK, MG, ML, ILG, LAY)

!         WRITE(22,*) MAXVAL(REI),MAXVAL(REL),MINVAL(REI),MINVAL(REL)

! COMPUTE THE CLOUD OPTICAL THICKNESS AT ~0.6 AND ~10.5 MICRONS
! THESE CO-EFFICIENTS ARE FROM J. LI (JAN. 22, 2006)

         DO ILAY = 1, LAY
            DO IL = IL1, IL2

               IF (CLDMASK(IL,ILAY) .GT. 0.0) THEN

                  INVREL = 1.0/REL(IL,ILAY)
                  INVREI = 1.0/REI(IL,ILAY)

! COMPUTE OPTICAL THICKNESS WATER CLOUDS
                  IF (CLW(IL,ILAY) .GT. 0.0) THEN
                   TAULIQVIS = CLW(IL,ILAY)*(4.483e-04 + INVREL * (1.501 &
                             + INVREL*(7.441e-01 - INVREL * 9.620e-01)))
                   TAULIQIR  = CLW(IL,ILAY)*(0.14532e-01 - 0.47449e-03 &
                             * REL(IL,ILAY) + INVREL*(0.22898e+01 - INVREL &
                             * (0.92402e+01 - INVREL * 0.100999e+02)))
                  ELSE
                     TAULIQVIS = 0.0
                     TAULIQIR  = 0.0
                  END IF

! COMPUTE OPTICAL THICKNESS ICE CLOUDS
                  IF (CIC(IL,ILAY) .GT. 0.0) THEN
                     TAUICEVIS = CIC(IL,ILAY) &
                               * (-0.303108e-04 + 0.251805e+01 * INVREI)

                     TAUICEIR = CIC(IL,ILAY) * (-7.627102e-03 + 3.406420 &
                              * INVREI - 1.732583e+01 * (INVREI**2))
                  ELSE
                     TAUICEVIS = 0.0
                     TAUICEIR  = 0.0
                  END IF

                  TAU_VIS(IL,ILAY) = (TAULIQVIS + TAUICEVIS)*DZ(IL,ILAY)
                  TAU_IR           = (TAULIQIR  + TAUICEIR)*DZ(IL,ILAY)
                  EMS_CLD(IL,ILAY) = 1.0-EXP(-RU*TAU_IR)
               ELSE
                  TAU_VIS(IL,ILAY)  = 0.0
                  EMS_CLD(IL,ILAY) = 0.0
               END IF

            END DO              ! IL
         END DO                 ! ILAY

! CALL THE ISCCP SIMULATOR
         CALL ISCCP_SIM(TAU_SUB, PTOP_SUB,                                 &! OUTPUT
                        IL1, IL2,ILG, LAY, CLD_HGT,                        &! INPUT
                        PFULL, PHALF, Q, CLDMASK, TAU_VIS, EMS_CLD, T, &
                        GT, EMS_SFC, SUNLIT_LOC)

! CONVERT PTOP FROM Pa TO hPa

         DO IL = IL1, IL2
            PTOP_SUB(IL) = PTOP_SUB(IL)/100.0
         END DO

! COMPUTE THE VARIOUS ISCCP QUANTITIES

         DO IL = IL1, IL2
! ADD CONTRIBUTIONS TO THE HISTOGRAM OF CLOUD TOP PRESSURE AND CLOUD OPTICAL THICKNESS
! DO SO IF SUN IS UP
            IF (SUNLIT_LOC(IL) .GT. 0) THEN

               IPRES = 0

               IP = 1
               IF (PTOP_SUB(IL) .GT. P_TOP_BNDS(IP) .AND. &
                    PTOP_SUB(IL) .LT. P_TOP_BNDS(IP+1)) THEN
                  IPRES = IP
               END IF

               DO IP = 2, NPTOP-1
                  IF (PTOP_SUB(IL) .GE. P_TOP_BNDS(IP) .AND. &
                       PTOP_SUB(IL) .LT. P_TOP_BNDS(IP+1)) THEN
                     IPRES = IP
                  END IF
               END DO ! IP

               IF (PTOP_SUB(IL) .GE. P_TOP_BNDS(NPTOP)) IPRES=NPTOP

               ITAU = 0
               IT = 1
                  IF (TAU_SUB(IL) .GT. TAUCHK .AND. &
                       TAU_SUB(IL) .LT. TAU_BNDS(IT+1)) THEN
                     ITAU = IT
                  END IF

               DO IT = 2, NTAU-1
                  IF (TAU_SUB(IL) .GE. TAU_BNDS(IT) .AND. &
                       TAU_SUB(IL) .LT. TAU_BNDS(IT+1)) THEN
                     ITAU = IT
                  END IF
               END DO ! IT

               IF (TAU_SUB(IL) .GE. TAU_BNDS(NPTOP)) ITAU=NTAU

               IF (IPRES .GT. 0 .AND. ITAU .GT. 0) THEN
                  ITP = (IPRES-1)*NTAU+ITAU
                  CLDPT(IL,ITP) = CLDPT(IL,ITP) + 1.0
               END IF
            END IF
         END DO ! IL

! ACCUMULATE THE DOMAIN-MEAN PROPERTIES

         DO IL = IL1, IL2
            IF (SUNLIT_LOC(IL) .GT. 0   .AND. &
                 PTOP_SUB(IL)  .GT. 0.0 .AND. &
                 TAU_SUB(IL)   .GT. TAUCHK) THEN

! LINEAR TAU
               CLDTAU(IL) = CLDTAU(IL) + TAU_SUB(IL)

! CLOUD LOG MEAN VISIBLE OPTICAL THICKNESS (ACTUALLY ALBEDO (RADIATIVE) MEAN)
               CLDEP(IL) = CLDEP(IL) &
                      + REAL(INVTAU(MIN(NINT(100.0*TAU_SUB(IL)),45000)))

! CLOUD MEAN CLOUD TOP PRESSURE
               CLDTP(IL) = CLDTP(IL) + PTOP_SUB(IL)

! TOTAL CLOUD FRACTION
               CLDFRAC(IL) = CLDFRAC(IL) + 1.0
            END IF
         END DO ! IL
      END DO ! ICOL

! NOW THAT WE HAVE WORKED OVER THE ALL OF THE SUBCOLUMNS COMPUTE THE
! PROPER MEANS

      DO IL = IL1, IL2
! VARIABLES RELATED TO CLOUD OPTICAL THICKNESS AND VARIABLITY
         IF (SUNLIT_LOC(IL) .GT. 0) THEN
            IF (CLDFRAC(IL) .GT. 0.0) THEN
               CLDTAU(IL)  = CLDTAU(IL)/CLDFRAC(IL)
               A           = CLDEP(IL)/CLDFRAC(IL)
               B           = TAUTAB(MIN(255,MAX(1,NINT(A))))
               C           = MIN(B,CLDTAU(IL))
               CLDEP(IL)   = 1.0-C/CLDTAU(IL)
               CLDTAU(IL)  = B ! USE RADIATIVELY AVERAGED CLOUD OPTICAL THICKNESS
               CLDTP(IL)   = CLDTP(IL)/CLDFRAC(IL)
               CLDFRAC(IL) = CLDFRAC(IL)/REAL(NSUBCOL)
               DO ITP = 1, NPTOP*NTAU
                 CLDPT(IL,ITP) = CLDPT(IL,ITP)/REAL(NSUBCOL)
               END DO ! ITP
            ELSE
               CLDTAU(IL)  = -999.0
               CLDEP(IL)   = -999.0
               CLDTP(IL)   = -999.0
               CLDFRAC(IL) = 0.0
            END IF
         ELSE
            CLDTAU(IL)  = -999.0
            CLDEP(IL)   = -999.0
            CLDTP(IL)   = -999.0
            CLDFRAC(IL) = -999.0
            DO ITP = 1, NPTOP*NTAU
               CLDPT(IL,ITP) = -999.0
            END DO ! ITP
         END IF
         SUNLIT(IL) = REAL(SUNLIT_LOC(IL))
      END DO ! IL

      END

