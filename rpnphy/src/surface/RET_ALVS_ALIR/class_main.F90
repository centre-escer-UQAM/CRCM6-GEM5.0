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
!** S/P CLASSIC main_driver
!
!#include "phy_macros_f.h"
subroutine class_main (BUS, BUSSIZ, &
                       PTSURF, PTSURFSIZ, &
                       DT, KOUNT, TRNCH, &
                       N, NK, IG)
!
  use mu_jdate_mod, only: jdate_day_of_year, mu_js2ymdhms
  use sfc_options
  use sfcbus_mod
  use phy_options, only: MU_JDATE_HALFDAY, RAD_NUVBRANDS
!  use ctem_params, only : initpftpars, nlat, ilg, ican, ignd
  use class_configs, only : CLASS_IC, CTEM_ICC
  use classicParams, only : nlat, ilg, ignd

  implicit none

  ! Input parameters
  integer :: BUSSIZ  
  real, target :: BUS(BUSSIZ)
  integer :: PTSURFSIZ
  integer :: PTSURF(PTSURFSIZ)
  real    :: DT
  integer :: KOUNT, TRNCH, N, NK
  integer :: IG  ! number of soil layers

  integer :: LONCLEF, VSIZ
  integer :: SURFLEN

  integer, parameter :: NBS = RAD_NUVBRANDS

!
!
!Author
!          Y. Delage (November 2002)
!Revisions
!001       Y. Delage  (Jul 2004) Add calculations at kount=0
!002       Y. Delage  (Sep 2004) Replace ZA by ZUSL and ZTSL and
!                                UE2 by FRV
!003       V. Fortin  (Nov 2006) Use RAINRATE and SNOWRATE estimated
!                                by SURF_PRECIP (instead of TSS)
!                                to obtain total precipitation
!004       R. Larocque(Apr 2006) Mosaic driver
!005       V. Fortin  (Nov 2006) Adapt driver to CLASS 3.2
!006       V. Fortin  (Apr 2007) Adapt driver to CLASS 3.3
!007       J.P. Paquin(Aug 2008) -Level modification in LOCBUS to 0
!                                 for coupled use of CLASS in GEM
!                                -Full CLASS calculations at KOUNT=0
!                                 and replacement of integer :: value in calls
!                                -Add outputs for FTEMP, FVAP, RIB to GEM's
!                                 bus
!                                -Add calculations for roughness lenghts
!                                -ADD IMPFLX adjustment on ALFAT & ALFAT
!                                 (necessary for coupled execution)
!                                -ADD CALCULATIONS FOR BM OUTSIDE CLASS
!                                 ORIGINAL CODE.
!008       L. Duarte  (Oct 2008) Adapt driver to CLASS 3.4
!009       L. Duarte  (Dec 2008) Implement variable number of soil layers
!                                -IG now a subroutine argument
!                                -DELZ and ZBOT initialised according to
!                                 SCHMSOL_LEV
!                                -Add BDEPTH
!010       B. Dugas   (Jan 2009) Use VAMIN from options.cdk
!011       L. Duarte  (Feb 2009) -Removed BDEPTH
!                                -Check if ROOTDP less or equal than SDEPTH
!012       K. Winger  (Mar 2009) -Update snow depth in surface bus (ZSNOW -> SNODP)
!          Winger & Paquin       -Make sure RAINRATE, SNOWRATE, and CANG
!                                 do not get changed in surface bus
!          B. Dugas              -Use CANG instead of COSZ from surface bus
!          Dugas  & Winger       -Declare 2-D arrays as allocatable
!          K. Winger             -Add dimensions to pointer declarations
!                                -Always set ISAND,THLIQ,THICE according to bedrock
!013       L. Duarte  (Apr 2009) -Use LOCBUS_MOS instead of LOCBUS with
!                                 Z0H and Z0M
!014       L. Duarte  (May 2009) -Make calculations on certain variables
!                                 (TFLUX, QFLUX, BM, Z0M & Z0T) only
!                                 on points covered by the current mosaic level
!015       J.Toviessi (Aug 2009) -Adding the option of radiation
!                                 along the slopes (RADSLOPE)
!017       D. Deacu   (Mar 2010) -Set Z0ORO = f(Z0)
!018       K. Winger  (Nov 2001) -Change GGEO to ZGGEO and set it to value from namelist
!019       L. Duarte  (May 2010) -Add code to run CTEM in coupled mode
!020       K. Winger  (Nov 2012) -Pass FSNOW to output, renamed to ZFSNOW
!021       L. Duarte  (    2014) -Update for use with CLASS 3.6/CTEM 2.0
!022       K. Winger  (Aug 2016) -Initialize if soil just appeard due to melting glacier
!023       K. Winger  (Jun 2020) -Adapted for GEM5 and CLASSIC
!
!Object
!          Driver of the surface scheme CLASS
!
!Arguments
!
!               - Input/Output -
! BUS           bus of surface variables
!
!               - Input -
! BUSSIZ        size of the surface bus
! PTSURF        surface pointers
! PTSURFSIZ     dimension of ptsurf
! KOUNT         number of timestep
! TRNCH         row number
! DT            timestep
! N             running length
! NK            vertical dimension
! IG            number of soil layers
!
!
!*
!
!

  integer, parameter :: INDX_SFC = INDX_SOIL

  real   , parameter :: ggeo = 0.
  integer, parameter :: ctem_stdaln = 1
  integer, parameter :: iwf = 0
  logical, parameter :: INTERFLOW_L = .false.

  integer :: I
  integer :: IC,ICP1,IPAI,IHGT,IALC,IALS,IALG,IPCP, ICC,ICCP1
  integer :: IDISP,IZREF,ITC,ITCG,ITG,ISLFD,NMIM       !,ILW
  integer :: NLANDCS,NLANDGS,NLANDC, NLANDG, NLANDI
  logical :: DOTILE, DOPREC
  logical :: kount0  !to determine if we should run the initializations
  !inside KOUNT_EQ_0
!
!
!
  integer :: ptr, x
!
  integer :: k,j,ik,iday,ij,c,g
  real    :: juliand
!
!#include "locbus.cdk"
!  integer :: QUELNIVO(MAXVARSURF)
  integer :: QUELNIVO(NVARSURF)
!
! WARNING !!!! x in bus(x(varname,1,1)) is defined in the line below
! it is now case sensitive
#define x(fptr,fj,fk) ptsurf(vd%fptr%i)+(fk-1)*surflen+fj-1
!
!
#include "thermoconsts.inc"
#include "dintern.inc"
#include "fintern.inc"
!
!
!******************************************************
!     AUTOMATIC ARRAYS
!******************************************************
!
      real :: &
           ALVSCN(N),   ALIRCN(N),   ALVSG (N),   ALIRG (N), &
           ALVSCS(N),   ALIRCS(N),   ALVSSN(N),   ALIRSN(N), &
           TRVSCN(N),   TRIRCN(N),   TRVSCS(N),   TRIRCS(N), &
           FSVF  (N),   FSVFS (N), &
           RAICAN(N),   RAICNS(N),   SNOCAN(N),   SNOCNS(N), &
           FRAINC(N),   FSNOWC(N),   FRAICS(N),   FSNOCS(N), &
           DISP  (N),   DISPS (N),   ZOMLNC(N),   ZOMLCS(N), &
           ZOELNC(N),   ZOELCS(N),   ZOMLNG(N),   ZOMLNS(N), &
           ZOELNG(N),   ZOELNS(N),   CHCAP (N),   CHCAPS(N), &
           CMASSC(N),   CMASCS(N),   RC    (N),   RCS   (N), &
           ZPLIMC(N),   ZPLIMG(N),   ZPLMCS(N),   ZPLMGS(N), &
           QLWAVG(N),   QSOL  (N), &
           FCLOUD(N),   VPD   (N), &
           RHOAIR(N),   TADP  (N),   QSWINV(N),   QSWINI(N), &
           PADRY (N),   ZBLEND(N),   ZUN   (N),   ZTN   (N), &
           CWLCAP(N),   CWFCAP(N),   CWLCPS(N),   CWFCPS(N), &
           RBCOEF(N), &
           WSNOCS(N),   WSNOGS(N),   RHOSCS(N),   RHOSGS(N), &
           CDRAG (N), &
           ALVSGC(N),   ALIRGC(N),   ALVSSC(N),   ALIRSC(N), &
! roughness modif 1
!     N     Z0ORO (N),   SNOLIM(N),   ZPLMG0(N),   ZPLMS0(N),
           Z0M(N),       SNOLIM(N),   ZPLMG0(N),   ZPLMS0(N), &
           PCPR  (N),   ZGGEO (N),   VMOD  (N), &
           SRATE (N),   RRATE (N),   COSZS (N)

! New local fields added for CLASSIC (KW)
      integer :: isnoalb       ! 0: original two-band snow albedo algorithms are used. 
                               ! 1: the new four-band routines are used. 
                               ! At present, the four band algorithm should NOT be used offline.
      real :: BCSNO(N)         ! Black carbon mixing ratio [kg/m^3]
      real :: TRSNOWC(N)       ! Transmissivity of snow under vegetation to shortwave radiation
      real :: TRSNOWG(N,NBS)   ! Transmissivity of snow in bare areas to shortwave radiation
      real :: ALTG(N,NBS)      
      real :: ALSNO(N,NBS)     ! Albedo of snow in each modelled wavelength band []
      real :: ZTHLW(n,ig)      ! Soil water content at wilting point $[m^3/m^3]
      real :: TCSNOW(N)        ! Thermal conductivity of snow  $[W/m/K ]
      real :: GSNOW(N)         ! Heat conduction into surface of snow pack [W/m^2 ]
      real :: GTBS(N)          ! Surface temperature for CCCma black carbon scheme  [K]
      real :: SFCUBS(N)        ! Zonal surface wind velocity for CCCma black carbon scheme [m/s]
      real :: SFCVBS(N)        ! Meridional surface wind velocity for CCCma black carbon scheme [m/s]
      real :: USTARBS(N)       ! Friction velocity for CCCma black carbon scheme [m/s]
      ! Peatland variables, not used when ipeatland=0, just place holders
      real :: DAYL(N)          ! DAYLENGTH FOR THAT LOCATION
      real :: DAYL_MAX(N)      ! MAXIMUM DAYLENGTH FOR THAT LOCATION
      real :: ancsmoss(n), angsmoss(n), ancmoss(n), angmoss(n), rmlcsmoss(n), rmlgsmoss(n), &
              rmlcmoss(n), rmlgmoss(n), Cmossmas(n), dmoss(n), pdd(n)


!
!      real :: ZTHRC (N,IG,2), ZTHRG (N,IG,2),
!     1     ZTHRCS(N,IG,2), ZTHRGS(N,IG,2)
!
      real :: GZEROC(N),   GZEROG(N),   GZROCS(N),   GZROGS(N), &
           G12C  (N),   G12G  (N),   G12CS (N),   G12GS (N), &
           G23C  (N),   G23G  (N),   G23CS (N),   G23GS (N), &
           QFREZC(N),   QFREZG(N),   QMELTC(N),   QMELTG(N), &
           EVAPC (N),   EVAPCG(N),   EVAPG (N),   EVAPCS(N), &
           EVPCSG(N),   EVAPGS(N),   TCANO (N),   TCANS (N), &
           TPONDC(N),   TPONDG(N), &
           TPNDCS(N),   TPNDGS(N),   TSNOCS(N),   TSNOGS(N), &
                        WTABLE(N), &
           EVPPOT(N),   EVAPB (N)
!
      real :: ASVDAT(N), ASIDAT(N), AGVDAT(N), AGIDAT(N)
!
!     * DIAGNOSTIC ARRAYS USED FOR CHECKING ENERGY AND WATER
!     * BALANCES (CLASSZ).
      real :: CTVSTP(N),   CTSSTP(N),   CT1STP(N),   CT2STP(N), &
           CT3STP(N),   WTVSTP(N),   WTSSTP(N),   WTGSTP(N)
!
      real :: RHOSNI(N)
      real :: RPCP  (N),   TRPCP (N),   SPCP  (N),   TSPCP (N)
      integer ::           ILAND (N),   ITER  (N),   NITER (N)

! Peatland field for bog and fen peatlands
! Joe is till working on this so for now this field does not get read in but set to 0
      integer :: ipeatland(n)


! Thickness and depth of each soil layer in meters
      real :: DELZ(IG),ZBOT(IG)
!     DATA  DELZ    /0.10,0.25,3.75/
!     DATA  ZBOT    /0.10,0.35,4.10/
!
!******************************************************
!
!
  integer, DIMENSION(:,:,:), ALLOCATABLE :: ITERCT
  integer, DIMENSION(:,:),   ALLOCATABLE :: ISAND,IORG
  integer, DIMENSION(:),     ALLOCATABLE :: SOCI
  REAL, DIMENSION(:,:), ALLOCATABLE :: PAIDAT,HGTDAT,ACVDAT,ACIDAT
  REAL, DIMENSION(:,:), ALLOCATABLE :: TBARC ,TBARG ,TBARCS,TBARGS, &
                                       THLIQC,THLIQG,THICEC,THICEG, &
                                       HCPC  ,HCPG  ,FROOT , FROOTS, GFLUX, &
                                       TCTOPC,TCBOTC,TCTOPG,TCBOTG
!
  real :: XDIFFUS(N)
  real :: zlyglfmas(N,CTEM_ICC)

!
  real,pointer,dimension(:)   :: ZALVS, ZALIR
  real,pointer,dimension(:)   :: FSOLUACC, FIRUACC, SU, SV, ST, SQ, ALVIS_SOL
  real,pointer,dimension(:)   :: EVAPO, QSENS, QEVAP, HBL, ZILMO, ZFRV, PS, QS
  real,pointer,dimension(:)   :: Z0H, Z0ORO, ZTSURF, ZTSRAD, UA, VA, TA, QA
  real,pointer,dimension(:)   :: TFLUX, QFLUX, ZCANG, FLUSOL, ZHUAIRCAN
  real,pointer,dimension(:)   :: ZDLAT, ZDLON, ZSDEPTH, ZZTSL, ZZUSL, QLWIN
  real,pointer,dimension(:)   :: ZTAIRCAN, ZTSNOW, ZTBASE, ZTPOND, ZZPOND
  real,pointer,dimension(:)   :: ZRHOSNO, ZRUNOFF
  real,pointer,dimension(:)   :: ZSCAN, ZINISOIL, XSNO, ZALBSNO, ZGROWTH
  real,pointer,dimension(:)   :: ZXDRAIN, ZXSLOPE, ZGRKFAC, ZWFSURF, ZWFCINT
  real,pointer,dimension(:)   :: ZCMAI, ZFSGV, ZFSGS, ZFSGG, ZFSNOW, ZFLGV
  real,pointer,dimension(:)   :: ZFLGS, ZFLGG, ZHFSC, ZHFSS, ZHFSG, ZHEVC
  real,pointer,dimension(:)   :: ZHEVS, ZHEVG, ZHMFC, ZHTCC, ZHTCS, ZPCFC, ZPCLC
  real,pointer,dimension(:)   :: ZPCPG, ZQFCF, ZQFCL, ZQFG, ZQFN
  real,pointer,dimension(:)   :: ZTBASFL, ZTOVRFL, ZTRUNOFF, ZTSUBFL, ZWSNOW
  real,pointer,dimension(:)   :: ZWTRC, ZWTRS, ZWTRG, ZROFC, ZROFN, ZROVG
  real,pointer,dimension(:)   :: ZOVRFLW, ZSUBFLW, ZBASFLW, ZTCAN, ZRCAN
  real,pointer,dimension(:)   :: ZCFLUX, ZPCPN, ZHMFN, ZALGWN, ZALGWV, ZALGDN, ZALGDV
  real,pointer,dimension(:)   :: FFC, FCS, FG, FGS, ZFL, ZRAINRATE, ZSNOWRATE
  real,pointer,dimension(:)   :: ZSNOW, ZSOILCOL

  real,pointer,dimension(:)   :: ZFTEMP, ZFVAP, ZRIB, ZCDH, ZCDM, ZBM

  real,pointer,dimension(:)   :: ztt2m

  real,pointer,dimension(:,:) :: THLIQ, THICE, TS, ZFCANMX, ZSAND, ZCLAY
  real,pointer,dimension(:,:) :: ZORGM, ZDELZW, ZZBOTW, ZTHPOR, ZTHLMIN
  real,pointer,dimension(:,:) :: ZTHLRET, ZPSISAT, ZBI, ZPSIWLT, ZHCPS
  real,pointer,dimension(:,:) :: ZTCS, ZTSFS, ZTHFC
  real,pointer,dimension(:,:) :: ZGRKSAT, ZTHLRAT, ZHTC
  real,pointer,dimension(:,:) :: ZZOLN, ZALVSC, ZALIRC, ZPAIMAX, ZPAIMIN
  real,pointer,dimension(:,:) :: ZCWGTMX, ZZRTMAX, ZRSMIN, ZQA50, ZVPDA, ZVPDB
  real,pointer,dimension(:,:) :: ZPSIGA, ZPSIGB, ZQFC, ZHMFG

  real,pointer,dimension(:,:) :: ZFCANCMX
  real,pointer,dimension(:,:) :: ZAILCG, ZAILCB, ZAILC, ZZOLNC, ZSLAI, ZBMASVEG
  real,pointer,dimension(:,:) :: ZCMASVEGC, ZVEGHGHT, ZROOTDPTH, ZALVSCTM, ZALIRCTM 
  real,pointer,dimension(:,:) :: ZPAIC, ZSLAIC
  real,pointer,dimension(:,:,:) :: ZRMATC, ZRMATCTEM

! New fields for CLASSIC (KW)
  real,pointer,dimension(:)   :: REFSNO, groundHeatFlux
  real,pointer,dimension(:,:) :: ZFSSB     !  input Total   solar radiation in each modelled wavelength band [W/m^2]
  real,pointer,dimension(:,:) :: ZFSDB     !  input Direct  solar radiation in each modelled wavelength band [W/m^2]
  real,pointer,dimension(:,:) :: ZFSFB     !  input Diffuse solar radiation in each modelled wavelength band [W/m^2]


! New arrays for coupling CLASS and CTEM
      REAL,POINTER,DIMENSION(:,:) :: ZGLEAFMAS
      REAL,POINTER,DIMENSION(:,:) :: ZBLEAFMAS
      REAL,POINTER,DIMENSION(:,:) :: ZSTEMMASS
      REAL,POINTER,DIMENSION(:,:) :: ZROOTMASS
      REAL,POINTER,DIMENSION(:,:) :: ZLITRMASS
      REAL,POINTER,DIMENSION(:,:) :: ZSOILCMAS
      REAL,POINTER,DIMENSION(:,:) :: ZMLIGHTNG
      REAL,POINTER,DIMENSION(:) :: ZPRBFRHUC
      REAL,POINTER,DIMENSION(:) :: ZEXTNPROB
      REAL,POINTER,DIMENSION(:,:) :: ZLFSTATUR
      REAL,POINTER,DIMENSION(:,:) :: ZPANDAYR
!      POINTER (IZSTDALNR  , ZSTDALNR   (M) )

      REAL,POINTER,DIMENSION(:) :: ZCO2CONC
      REAL,POINTER,DIMENSION(:,:) :: ZCO2I1CS
      REAL,POINTER,DIMENSION(:,:) :: ZCO2I1CG
      REAL,POINTER,DIMENSION(:,:) :: ZCO2I2CS
      REAL,POINTER,DIMENSION(:,:) :: ZCO2I2CG

      REAL,POINTER,DIMENSION(:) :: ZUVACCGAT
      REAL,POINTER,DIMENSION(:) :: ZVVACCGAT

      REAL,POINTER,DIMENSION(:) :: ZFSNOWACC
      REAL,POINTER,DIMENSION(:) :: ZTCANOACCGAT
      REAL,POINTER,DIMENSION(:) :: ZTCANSACC
      REAL,POINTER,DIMENSION(:) :: ZTAACCGAT
      REAL,POINTER,DIMENSION(:,:) :: ZANCSVGAC
      REAL,POINTER,DIMENSION(:,:) :: ZANCGVGAC
      REAL,POINTER,DIMENSION(:,:) :: ZRMLCSVGA
      REAL,POINTER,DIMENSION(:,:) :: ZRMLCGVGA
      REAL,POINTER,DIMENSION(:,:) :: ZTODFRAC

      REAL,POINTER,DIMENSION(:,:) :: ZTBARACCGAT
      REAL,POINTER,DIMENSION(:,:) :: ZTBARCACC
      REAL,POINTER,DIMENSION(:,:) :: ZTBARCSACC
      REAL,POINTER,DIMENSION(:,:) :: ZTBARGACC
      REAL,POINTER,DIMENSION(:,:) :: ZTBARGSACC
      REAL,POINTER,DIMENSION(:,:) :: ZTHLIQCACC
      REAL,POINTER,DIMENSION(:,:) :: ZTHLIQGACC
      REAL,POINTER,DIMENSION(:,:) :: ZTHICECACC

      REAL,POINTER,DIMENSION(:,:) :: ZDVDFCAN

      REAL,POINTER,DIMENSION(:) :: ZVGBIOMAS
      REAL,POINTER,DIMENSION(:) :: ZGAVGLAI
      REAL,POINTER,DIMENSION(:) :: ZGAVGLTMS
      REAL,POINTER,DIMENSION(:) :: ZGAVGSCMS
      REAL,POINTER,DIMENSION(:,:) :: ZCOLDDAYR

      REAL,POINTER,DIMENSION(:,:) :: ZGRWTHEFF
      REAL,POINTER,DIMENSION(:,:) :: ZFLHRLOSS
      REAL,POINTER,DIMENSION(:,:) :: ZSTMHRLOS
      REAL,POINTER,DIMENSION(:,:) :: ZROTHRLOS
      REAL,POINTER,DIMENSION(:,:) :: ZLYSTMMAS
      REAL,POINTER,DIMENSION(:,:) :: ZLYROTMAS
      REAL,POINTER,DIMENSION(:,:) :: ZTYMAXLAI


      REAL,POINTER,DIMENSION(:) :: ZCFLUXCG
      REAL,POINTER,DIMENSION(:) :: ZCFLUXCS

      REAL,POINTER,DIMENSION(:,:) :: ZPFCANCMX
      REAL,POINTER,DIMENSION(:,:) :: ZNFCANCMX

      REAL,POINTER,DIMENSION(:) :: ZNPP
      REAL,POINTER,DIMENSION(:) :: ZNEP
      REAL,POINTER,DIMENSION(:) :: ZHETRORES
      REAL,POINTER,DIMENSION(:) :: ZAUTORES
      REAL,POINTER,DIMENSION(:) :: ZSOILRESP
      REAL,POINTER,DIMENSION(:) :: ZRM
      REAL,POINTER,DIMENSION(:) :: ZRG
      REAL,POINTER,DIMENSION(:) :: ZNBP
      REAL,POINTER,DIMENSION(:) :: ZLITRES
      REAL,POINTER,DIMENSION(:) :: ZSOCRES
      REAL,POINTER,DIMENSION(:) :: ZGPP
      REAL,POINTER,DIMENSION(:) :: ZDSTCEMLS
      REAL,POINTER,DIMENSION(:) :: ZLITRFALL
      REAL,POINTER,DIMENSION(:) :: ZHUMIFTRS
      REAL,POINTER,DIMENSION(:) :: ZRML
      REAL,POINTER,DIMENSION(:) :: ZRMS
      REAL,POINTER,DIMENSION(:) :: ZRMR
      REAL,POINTER,DIMENSION(:,:) :: ZTLTRLEAF
      REAL,POINTER,DIMENSION(:,:) :: ZTLTRSTEM
      REAL,POINTER,DIMENSION(:,:) :: ZTLTRROOT
      REAL,POINTER,DIMENSION(:,:) :: ZLEAFLITR
      REAL,POINTER,DIMENSION(:,:) :: ZROOTTEMP
      REAL,POINTER,DIMENSION(:,:) :: ZAFRLEAF
      REAL,POINTER,DIMENSION(:,:) :: ZAFRSTEM
      REAL,POINTER,DIMENSION(:,:) :: ZAFRROOT
      REAL,POINTER,DIMENSION(:,:) :: ZWTSTATUS
      REAL,POINTER,DIMENSION(:,:) :: ZLTSTATUS
      REAL,POINTER,DIMENSION(:) :: ZBURNAREA
      REAL,POINTER,DIMENSION(:) :: ZPROBFIRE
      REAL,POINTER,DIMENSION(:) :: ZLUCEMCOM
      REAL,POINTER,DIMENSION(:) :: ZLUCLTRIN
      REAL,POINTER,DIMENSION(:) :: ZLUCSOCIN
      REAL,POINTER,DIMENSION(:,:) :: ZNPPVEG
      REAL,POINTER,DIMENSION(:) :: ZGRCLAREA
      REAL,POINTER,DIMENSION(:) :: ZDSTCEMLS3
      REAL,POINTER,DIMENSION(:,:) :: ZRMLVEGACC
      REAL,POINTER,DIMENSION(:,:) :: ZRMSVEG
      REAL,POINTER,DIMENSION(:,:) :: ZRMRVEG
      REAL,POINTER,DIMENSION(:,:) :: ZRGVEG
      REAL,POINTER,DIMENSION(:,:) :: ZVGBIOMAS_VEG
      REAL,POINTER,DIMENSION(:,:) :: ZGPPVEG
      REAL,POINTER,DIMENSION(:,:) :: ZNEPVEG

      REAL,POINTER,DIMENSION(:,:,:) :: ZRMAT
      REAL,POINTER,DIMENSION(:,:) :: ZAIL
      REAL,POINTER,DIMENSION(:,:) :: ZPAI

      REAL,POINTER,DIMENSION(:) :: ZFSINACC
      REAL,POINTER,DIMENSION(:) :: ZFLINACC
      REAL,POINTER,DIMENSION(:) :: ZFLUTACC
      REAL,POINTER,DIMENSION(:) :: ZALSWACC
      REAL,POINTER,DIMENSION(:) :: ZALLWACC
      REAL,POINTER,DIMENSION(:) :: ZPREACC

      REAL,POINTER,DIMENSION(:) :: ZTCURM
      REAL,POINTER,DIMENSION(:) :: ZSRPCURYR
      REAL,POINTER,DIMENSION(:) :: ZDFTCURYR
      REAL,POINTER,DIMENSION(:,:) :: ZTMONTHB
      REAL,POINTER,DIMENSION(:) :: ZANPCPCUR
      REAL,POINTER,DIMENSION(:) :: ZANPECUR
      REAL,POINTER,DIMENSION(:) :: ZGDD5CUR
      REAL,POINTER,DIMENSION(:) :: ZSURMNR
      REAL,POINTER,DIMENSION(:) :: ZDEFMNR
      REAL,POINTER,DIMENSION(:) :: ZSRPLSCUR
      REAL,POINTER,DIMENSION(:) :: ZDEFCTCUR
!      REAL,POINTER,DIMENSION(:,:) :: ZLYGLFMAS
      REAL,POINTER,DIMENSION(:,:) :: ZGEREMORT
      REAL,POINTER,DIMENSION(:,:) :: ZINTRMORT
      REAL,POINTER,DIMENSION(:,:) :: ZLAMBDA
      REAL,POINTER,DIMENSION(:,:) :: ZPFTEXISTR
      REAL,POINTER,DIMENSION(:) :: ZTWARMM
      REAL,POINTER,DIMENSION(:) :: ZTCOLDM
      REAL,POINTER,DIMENSION(:) :: ZGDD5
      REAL,POINTER,DIMENSION(:) :: ZARIDITY
      REAL,POINTER,DIMENSION(:) :: ZSRPLSMON
      REAL,POINTER,DIMENSION(:) :: ZDEFCTMON
      REAL,POINTER,DIMENSION(:) :: ZANNDEFCT
      REAL,POINTER,DIMENSION(:) :: ZANNSRPLS
      REAL,POINTER,DIMENSION(:) :: ZANNPCP
      REAL,POINTER,DIMENSION(:) :: ZANPOTEVP
      REAL,POINTER,DIMENSION(:) :: ZDRYSLEN
      REAL,POINTER,DIMENSION(:,:) :: ZWDMINDEX
      REAL,POINTER,DIMENSION(:,:) :: ZBURNVEGF
      REAL,POINTER,DIMENSION(:,:) :: ZPSTEMMASS
      REAL,POINTER,DIMENSION(:,:) :: ZPGLEAFMASS
      REAL,POINTER,DIMENSION(:,:) :: ZCOLRATE
      REAL,POINTER,DIMENSION(:,:) :: ZMORTRATE
      REAL,POINTER,DIMENSION(:,:) :: ZAUTRESVEG
      REAL,POINTER,DIMENSION(:,:) :: ZHETRESVEG
      REAL,POINTER,DIMENSION(:,:) :: ZLITRESVEG
      REAL,POINTER,DIMENSION(:,:) :: ZNBPVEG
      REAL,POINTER,DIMENSION(:,:) :: ZSOCRESVEG


! Variables for coupling CLASS and CTEM
  real    :: DELTAT, ABSZERO
  logical :: CTEM_ON, compete, mosaic, inibioclim, dofire, LNDUSEON
  logical :: dowetlands, obswetf
  integer :: L2MAX, SPINFAST
  integer :: MODELPFT(12), NOL2PFTS(CLASS_IC)
  integer :: ISUMC, K1C, K2C, ICOUNT
  real    :: CSUM
  integer :: MONTH1, MONTH2, XDAY

  integer :: COLDDAYS(n,2),LFSTATUS(n,CTEM_ICC),PANDAYS(n,CTEM_ICC)
  real    :: LIGHTNG(n)
! for memory tests
  real    :: ZRH(n),ZANCSVEG(n,CTEM_ICC),ZANCGVEG(n,CTEM_ICC)
  real    :: ZRMLCSVEG(n,CTEM_ICC),ZRMLCGVEG(n,CTEM_ICC),ZCANRES(n)
  real    :: ZAILCGS(n,CTEM_ICC),ZFCANCS(n,CTEM_ICC),ZFCANC(n,CTEM_ICC)

  real    :: srh(n), popdin
  real    :: fieldsm(n,ig), wiltsm(n,ig),bterm(n,ig)

! Competition related variables
  real    :: altot_gat, fsstar_gat, flstar_gat, netrad_gat(n)
!
  real    :: tmonth(12,n)
  integer :: surmncur(n), defmncur(n)
  integer, dimension(n,12) :: wet_dry_mon_index
  logical :: pftexist(n,CTEM_ICC)

! Fire-related variables
! Not used for now (LD)
  real    :: emit_co2gat(n,CTEM_ICC), &
             emit_cogat(n,CTEM_ICC), &
             emit_ch4gat(n,CTEM_ICC), &
             emit_nmhcgat(n,CTEM_ICC), &
             emit_h2gat(n,CTEM_ICC), &
             emit_noxgat(n,CTEM_ICC), &
             emit_n2ogat(n,CTEM_ICC), &
             emit_pm25gat(n,CTEM_ICC), &
             emit_tpmgat(n,CTEM_ICC), &
             emit_tcgat(n,CTEM_ICC), &
             emit_ocgat(n,CTEM_ICC), &
             emit_bcgat(n,CTEM_ICC), &
!             burnvegfgat(n,CTEM_ICC), &
             btermgat(n), &
             ltermgat(n), &
             mtermgat(n)

! Methane(wetland) related variables
  real    :: WETFRACGRD(n), wetfrac_sgrd(n,8), &
             CH4WET1GAT(n), &
             CH4WET2GAT(n), &
             WETFDYNGAT(n), &
             CH4DYN1GAT(n), &
             CH4DYN2GAT(n)             !, wetfrac_mon(n,12)
!
  integer :: nml, igdr(n),    ilmos(n), jlmos(n)




! Huziy for debug
  logical :: ok_status

  integer :: Isat(n),Isats(n),Ibed(n)
  real    :: thpf(n,ig),RW(n,ig),DMLIQT(n),XWT(n)
  real    :: QDIST(n), QINT(n), QINC(n),QING(n),QINGS(n),QINCS(n),QDISC(n)
  integer :: jji
  real    :: dumw, sdlz
  real, parameter ::  maxd = 0.3

  logical :: prints_L

!===================== RLi =====================================\
! CONSTANTS FOR COUPLING CLASS AND CTEM
!
  DATA L2MAX/3/
!
! SEPARATION OF THESE PFTs INTO LEVEL 1 (FOR CLASS) AND LEVEL 2 (FOR CTEM) PFTs.
  DATA MODELPFT/1,     1,     0, &    ! CLASS PFT 1 NDL
!              EVG    DCD
                1,     1,     1, &    ! CLASS PFT 2 BDL
!              EVG  DCD-CLD DCD-DRY   ! NOTE 2 TYPES OF BDL DCD - COLD & DRY
                1,     1,     0, &    ! CLASS PFT 3 CROP
!              C3      C4
                1,     1,     0/    ! CLASS PFT 4 GRASS
!              C3      C4
!
! CTEM's TIME STEP IN DAYS
  DATA DELTAT/1.0/
!
  DATA ABSZERO/1E-06/

! ==============================================================================
! ==============================================================================

  prints_L = .false.
!if (kount==3) stop

  IC   = CLASS_IC
  ICP1 = CLASS_IC+1
  ICC  = CTEM_ICC
  ICCP1= CTEM_ICC+1
  ignd = ig
  ilg  = N

! Options and switches - likely to be set along with the other options in GEM
!
! CTEM1=(CTEM_MODE.gt.0)
! CTEM2=(CTEM_MODE.ge.2)
  ctem_on=(ctem_mode.gt.0)
  compete=ctem_compete
  dowetlands=.false.
  obswetf=.false.

  mosaic=.false.  ! mosaic option for competition purposes
                       ! always set it to false, even when mosaic is used! (LD)
  dofire=.false.
  LNDUSEON=.false.

!
  ALLOCATE( PAIDAT(N,IC),HGTDAT(N,IC),ACVDAT(N,IC),ACIDAT(N,IC), &
            TBARC (N,IG),TBARG (N,IG),TBARCS(N,IG),TBARGS(N,IG), &
            THLIQC(N,IG),THLIQG(N,IG),THICEC(N,IG),THICEG(N,IG), &
            HCPC  (N,IG),HCPG  (N,IG),FROOT (N,IG),FROOTS(N,IG), GFLUX (N,IG), &
            TCTOPC(N,IG),TCBOTC(N,IG),TCTOPG(N,IG),TCBOTG(N,IG), &
            ISAND (N,IG),IORG  (N,IG),ITERCT(N,6,50), SOCI(N) )
!
  SURFLEN = N
!  IDAY = JULIAND( DT , KOUNT, DATE )
  IDAY = real(jdate_day_of_year(jdateo + kount*int(dt) + MU_JDATE_HALFDAY))
! Initialize thickness and depth of soil layers
  DELZ(1) = SCHMSOL_LEV(1)
  ZBOT(1) = SCHMSOL_LEV(1)
  DO J=2,IG
     DELZ(J) = SCHMSOL_LEV(J)
     ZBOT(J) = SCHMSOL_LEV(J) + ZBOT(J-1)
  ENDDO
!
!     EQUIVALENCES
!
!cat sfc_businit.F90 ../../../rpnphy/src/base/phyvar.hf | grep -iw VN=VEGGRO
!grep -iw intent *.[fF]* | grep -iw PCFC
  PS        (1:N)      => bus( x(PMOINS ,1,1) : )  !  input surface pressure at t-dt (actually t+dt currently)
  ZCLAY     (1:N,1:IG) => bus( x(CLAY   ,1,1) : )  !  input percentage of clay in soil
  ZORGM     (1:N,1:IG) => bus( x(ORGM   ,1,1) : )  !  input percentage of organic matter in soil
  ZCANG     (1:N)      => bus( x(CANG   ,1,1) : )  !  input cosine of the solar angle
  ZSAND     (1:N,1:IG) => bus( x(SAND   ,1,1) : )  !  input percentage of sand in soil
  TA        (1:N)      => bus( x(TMOINS ,1,nk) : ) !  input temperature at t-dt
  QA        (1:N)      => bus( x(HUMOINS,1,nk) : ) !  input specific humidity at t-dt
  ZDLAT     (1:N)      => bus( x(DLAT   ,1,1) : )  !  input latitude
  ZDLON     (1:N)      => bus( x(DLON   ,1,1) : )  !  input longitude
  UA        (1:N)      => bus( x(UMOINS ,1,nk) : ) !  input wind speed along-x at t-dt
  VA        (1:N)      => bus( x(VMOINS ,1,nk) : ) !  input wind speed along-y at t-dt
  ZZTSL     (1:N)      => bus( x(ZTSL   ,1,1) : )  !  input temperature level at top of surface layer
  ZZUSL     (1:N)      => bus( x(ZUSL   ,1,1) : )  !  input wind level at top of surface layer
  IF (RADSLOPE) THEN
     FLUSOL (1:N)      => bus( x(FLUSLOP,1,1) : )  !  input SW flux twrds slp ground
  ELSE
     FLUSOL (1:N)      => bus( x(FLUSOLIS,1,1) : ) !  input SW flux towards ground
  ENDIF
  QLWIN     (1:N)      => bus( x(FDSI   ,1,1) : )  !  input accum. IR energy flux towards ground
  ZXDRAIN   (1:N)      => bus( x(XDRAIN ,1,1) : )  !  input drainage factor in CLASS
  ZGRKFAC   (1:N)      => bus( x(GRKFAC ,1,1) : )  !  input WATROF par. for MESH code - not used
  ZWFSURF   (1:N)      => bus( x(WFSURF ,1,1) : )  !  input WATROF par. for MESH code - not used
  ZWFCINT   (1:N)      => bus( x(WFCINT ,1,1) : )  !  input WATROF par. for MESH code - not used
  ZRAINRATE (1:N)      => bus( x(RAINRATE,1,1) : ) !  input liquid precip. rate
  ZSNOWRATE (1:N)      => bus( x(SNOWRATE,1,1) : ) !  input solid  precip. rate
  ZSNOW     (1:N)      => bus( x(SNODP  ,1,indx_sfc ) : )  !  inout snow depth

! New CLASSIC fields (KW)
  REFSNO    (1:N)      => bus( x(SNOWSIZE,1,1) : ) !  inout Snow grain size [m]
  groundHeatFlux(1:N)  => bus( x(grdhflx,1,1) : )  ! output Heat flux at soil surface [W/m^2]
  ZFSSB     (1:N,1:NBS) => bus( x(FATB,1,1) : )    !  input Total   solar radiation in each modelled wavelength band [W/m^2]
  ZFSDB     (1:N,1:NBS) => bus( x(FADB,1,1) : )    !  input Direct  solar radiation in each modelled wavelength band [W/m^2]
  ZFSFB     (1:N,1:NBS) => bus( x(FAFB,1,1) : )    !  input Diffuse solar radiation in each modelled wavelength band [W/m^2]


!print *,'class_main 000: ZSAND(1,:):', ZSAND(1,:)
!print *,'class_main 000: trnch,ZSNOW(:):', trnch,ZSNOW(:)

  igdr = 1
  DO J=1,IG
     DO I=1,N
        ISAND(I,J) = NINT(ZSAND(I,J))
        IF (ISAND(I,J).GT.-3) IGDR(I) = J  ! Index of soil layer in which bedrock is encountered
        bterm(i,j) = 0.159*zclay(i,j)+2.91
     ENDDO
  ENDDO
!print *,'class_main 999: ISAND(1,:):', ISAND(1,:)

! Peatland field for bog and fen peatlands
! Joe is till working on this so for now this field does not get read in but set to 0
  ipeatland = 0


  DO I=1,N
     ZBLEND(I) = ZZUSL(I)
     ILAND(I)  = I
     ZUN(I)    = ZU
     ZTN(I)    = ZT
     QSWINV(I) = 0.5*FLUSOL(I)
     QSWINI(I) = 0.5*FLUSOL(I)
     RRATE(I)  = ZRAINRATE(I)*1000.
     SRATE(I)  = ZSNOWRATE(I)*1000.
     PCPR(I)   = SRATE(I)+RRATE(I)
!-----------------------------------------------------------------
!  correctif pour le cas ou COSZS n'est pas disponible
     COSZS(I)  = ZCANG(I)
     if (FLUSOL(I).gt.20. .and. COSZS(I).le.0.) &
        COSZS(I) = FLUSOL(I)*.0006
!----------------------------------------------------------------
     IF (ABS(COSZS(I)).LT.0.10) COSZS(I) = 0.10
     IF (PCPR(I).GT.0.) THEN
        XDIFFUS(I) = 1.0
     ELSE
        XDIFFUS(I) = MIN(1.0-0.9*COSZS(I),1.)
     ENDIF
     FCLOUD(I) = XDIFFUS(I)
     QSOL(I)   = MAX(QSWINV(I)/COSZS(I),0.)
! roughness modif 4
!     Z0ORO(I)=0
     SNOLIM(I) = 0.10
     ZPLMG0(I) = 0.10
     ZPLMS0(I) = 0.10
     ZGGEO(I)  = GGEO
     vmod(i)   = sqrt(ua(i)*ua(i)+va(i)*va(i))
  END DO
!
!
!  DO_MOSAIC : do mos=dmos,nmos
    ZALVS     (1:N) => bus( x( ALVS   ,1,1 ) : )         ! output Diagnosed total visible albedo of land surface []
    ZALIR     (1:N) => bus( x( ALIR   ,1,1 ) : )         ! output Diagnosed total near-infrared albedo of land surface []
    Z0H       (1:N) => bus( x( Z0T    ,1,indx_sfc ) : )  !  input thermal  roughness length
    Z0ORO     (1:N) => bus( x( Z0     ,1,indx_sfc ) : )  !  input momentum roughness length
    QFLUX     (1:N) => bus( x( ALFAQ  ,1,1 ) : )         ! output Product of surface drag coefficient, wind speed and surface-air specific humidity difference [m/s] (inhomog. term for Q diff.)
    TFLUX     (1:N) => bus( x( ALFAT  ,1,1 ) : )         ! output Product of surface drag coefficient, wind speed and surface-air (inhomog. term for T diff.)
    ZCFLUX    (1:N) => bus( x( BT     ,1,1 ) : )         ! output Diagnosed product of drag coefficient and wind speed over modelled area [m/s] (homog. term for T,Q diffu.)
    ZBM       (1:N) => bus( x( BM     ,1,1 ) : )         ! output homog. term for U,V diffu.
    EVAPO     (1:N) => bus( x( WFLUX  ,1,1 ) : )         ! output Diagnosed total surface water vapour flux over modelled area [kg/m^2/s]
    QSENS     (1:N) => bus( x( FC     ,1,indx_sfc ) : )  ! output surface sensible heat flux
    QEVAP     (1:N) => bus( x( FV     ,1,indx_sfc ) : )  ! output surface latent   heat flux
    FFC       (1:N) => bus( x( FCOVC  ,1,1 ) : )         ! output fractional coverage for canopy
    FCS       (1:N) => bus( x( FCOVCS ,1,1 ) : )         ! output fractional coverage for canopy+snow
    FG        (1:N) => bus( x( FCOVG  ,1,1 ) : )         ! output fractional coverage for bare ground
    FGS       (1:N) => bus( x( FCOVGS ,1,1 ) : )         ! output fractional coverage for snow
    ZFRV      (1:N) => bus( x( FRV    ,1,indx_sfc ) : )  ! output Friction velocity of air [m/s]
    HBL       (1:N) => bus( x( HST    ,1,indx_sfc ) : )  ! output Height of the atmospheric boundary layer [m]
    ZILMO     (1:N) => bus( x( ILMO   ,1,indx_sfc ) : )  ! output Inverse of Monin-Obukhov roughness length [1/m]
    TS        (1:N,1:IG) => bus( x( TSOIL  ,1,1 ) : )    !  inout Temperature of soil layers [K]
    THLIQ     (1:N,1:IG) => bus( x( WSOIL  ,1,1 ) : )    !  inout Volumetric liquid water content of soil layers [m^3/m^3]
    THICE     (1:N,1:IG) => bus( x( ISOIL  ,1,1 ) : )    !  inout Volumetric frozen water content of soil layers [m^3/m^3]
    QS        (1:N) => bus( x( QSURF  ,1,indx_sfc ) : )  ! output skin specific humidity [kg/kg]
    ZTSURF    (1:N) => bus( x( TSURF  ,1,indx_sfc ) : )  ! output skin temperature [K]
    SQ        (1:N) => bus( x( QDIAG  ,1,1 ) : )         ! output screen level specific humidity
    ST        (1:N) => bus( x( TDIAG  ,1,1 ) : )         ! output screen level temperature
    SU        (1:N) => bus( x( UDIAG  ,1,1 ) : )         ! output screen level X-component of wind
    SV        (1:N) => bus( x( VDIAG  ,1,1 ) : )         ! output screen level X-component of wind
    XSNO      (1:N) => bus( x( SNOMA  ,1,1 ) : )         !  inout Mass of snow pack \f$[kg m^{-2}] (W_s)
    ZALBSNO   (1:N) => bus( x( SNOAL  ,1,1 ) : )         !  inout Snow albedo
    ZALGDN    (1:N) => bus( x( ALGDN  ,1,1 ) : )         !  inout Near-IR albedo of dry soil for modelled area
    ZALGDV    (1:N) => bus( x( ALGDV  ,1,1 ) : )         !  inout Visible albedo of dry soil for modelled area
    ZALGWN    (1:N) => bus( x( ALGWN  ,1,1 ) : )         !  inout Near-IR albedo of wet soil for modelled area
    ZALGWV    (1:N) => bus( x( ALGWV  ,1,1 ) : )         !  inout Visible albedo of wet soil for modelled area
    ZBASFLW   (1:N) => bus( x( DRAIN  ,1,1 ) : )         ! output Base flow from bottom of soil column [m]
    ZBI       (1:N,1:IG) => bus( x( BBI    ,1,1 ) : )    !  inout Clapp and Hornberger empirical "b" parameter
    ZCMAI     (1:N) => bus( x( CMAI   ,1,1 ) : )         ! output Aggregated mass of vegetation canopy [kg/m^2]
    ZDELZW    (1:N,1:IG) => bus( x( DELZW  ,1,1 ) : )    !  inout Permeable thickness of soil layer [m]
    ZFLGG     (1:N) => bus( x( FLGG   ,1,1 ) : )         ! output Diagnosed net longwave radiation at soil surface [W/m^2]
    ZFLGS     (1:N) => bus( x( FLGS   ,1,1 ) : )         ! output Diagnosed net longwave radiation at snow surface [W/m^2]
    ZFLGV     (1:N) => bus( x( FLGV   ,1,1 ) : )         ! output Diagnosed net longwave radiation on vegetation canopy [W/m^2]
    ZFSGG     (1:N) => bus( x( FSGG   ,1,1 ) : )         ! output Diagnosed net shortwave radiation at soil surface [W/m^2]
    ZFSGS     (1:N) => bus( x( FSGS   ,1,1 ) : )         ! output Diagnosed net shortwave radiation at snow surface [W/m^2]
    ZFSGV     (1:N) => bus( x( FSGV   ,1,1 ) : )         ! output Diagnosed net shortwave radiation on vegetation canopy [W/m^2]
    ZFSNOW    (1:N) => bus( x( FSNOW  ,1,1 ) : )         !  inout Diagnosed fractional snow coverage
    ZGRKSAT   (1:N,1:IG) => bus( x( GRKSAT ,1,1 ) : )    !  inout Hydraulic conductivity of soil at saturation [m/s]
    ZGROWTH   (1:N) => bus( x( VEGGRO ,1,1 ) : )         !  inout Vegetation growth index []
    ZHCPS     (1:N,1:IG) => bus( x( HCPS   ,1,1 ) : )    !  inout Volumetric heat capacity of soil matter $[J/m^3/K]
    ZHEVC     (1:N) => bus( x( HEVC   ,1,1 ) : )         ! output Diagnosed latent heat flux on vegetation canopy [W/m^2]
    ZHEVG     (1:N) => bus( x( HEVG   ,1,1 ) : )         ! output Diagnosed latent heat flux at soil surface [W/m^2]
    ZHEVS     (1:N) => bus( x( HEVS   ,1,1 ) : )         ! output Diagnosed latent heat flux at snow surface [W/m^2]
    ZHFSC     (1:N) => bus( x( HFSC   ,1,1 ) : )         ! output Diagnosed sensible heat flux on vegetation canopy [W/m^2]
    ZHFSG     (1:N) => bus( x( HFSG   ,1,1 ) : )         ! output Diagnosed sensible heat flux at soil surface [W/m^2]
    ZHFSS     (1:N) => bus( x( HFSS   ,1,1 ) : )         ! output Diagnosed sensible heat flux at snow surface [W/m^2]
    ZHMFC     (1:N) => bus( x( HMFC   ,1,1 ) : )         ! output Diagnosed energy associated with phase change
    ZHMFG     (1:N,1:IG) => bus( x( HMFG   ,1,1 ) : )    ! output Diagnosed energy associated with phase change of water in soil layers [W/m^2]
    ZHMFN     (1:N) => bus( x( HMFN   ,1,1 ) : )         ! output Diagnosed energy associated with freezing or thawing of water in the snow pack [W/m^2]
    ZHTC      (1:N,1:IG) => bus( x( HTC    ,1,1 ) : )    ! output Diagnosed internal energy change of soil layer due to conduction and/or change in mass [W/m^2]
    ZHTCC     (1:N) => bus( x( HTCC   ,1,1 ) : )         ! output Diagnosed internal energy change of canopy due to changes in temperature and/or mass [W/m^2]
    ZHTCS     (1:N) => bus( x( HTCS   ,1,1 ) : )         ! output Diagnosed internal energy change of snow pack due to conduction and/or change in mass [W/m^2]
    ZOVRFLW   (1:N) => bus( x( OVERFL ,1,1 ) : )         ! output Overland flow from top of soil column [m] 
    ZPCFC     (1:N) => bus( x( PCFC   ,1,1 ) : )         ! output Diagnosed frozen precipitation intercepted by vegetation [kg/m^2/s]
    ZPCLC     (1:N) => bus( x( PCLC   ,1,1 ) : )         ! output Diagnosed liquid precipitation intercepted by vegetation [kg/m^2/s]
    ZPCPG     (1:N) => bus( x( PCPG   ,1,1 ) : )         ! output Diagnosed precipitation incident on ground [kg/m^2/s]
    ZPCPN     (1:N) => bus( x( PCPN   ,1,1 ) : )         ! output Diagnosed precipitation incident on snow pack [kg/m^2/s]
    ZPSISAT   (1:N,1:IG) => bus( x( PSISAT ,1,1 ) : )    !  inout Soil moisture suction at saturation [m]
    ZPSIWLT   (1:N,1:IG) => bus( x( PSIWLT ,1,1 ) : )    !  inout Soil moisture suction at wilting point [m]
    ZQFC      (1:N,1:IG) => bus( x( QFC    ,1,1 ) : )    ! output Diagnosed water removed from soil layers by transpiration [kg/m^2/s]
    ZQFCF     (1:N) => bus( x( QFCF   ,1,1 ) : )         ! output Diagnosed sublimation from frozen water on vegetation [kg/m^2/s]
    ZQFCL     (1:N) => bus( x( QFCL   ,1,1 ) : )         ! output Diagnosed evaporation from liquid water on vegetation [kg/m^2/s]
    ZQFG      (1:N) => bus( x( QFG    ,1,1 ) : )         ! output Diagnosed evaporation from ground [kg/m^2/s]
    ZQFN      (1:N) => bus( x( QFN    ,1,1 ) : )         ! output Diagnosed sublimation from snow pack [kg/m^2/s]
    ZRCAN     (1:N) => bus( x( WVEG   ,1,1 ) : )         !  inout Intercepted liquid water stored on canopy [kg/m^2]
    ZRHOSNO   (1:N) => bus( x( SNODEN ,1,1 ) : )         !  inout Density of snow pack [kg/m^3]
    ZROFC     (1:N) => bus( x( ROFC   ,1,1 ) : )         ! output Liquid/frozen water runoff from vegetation [kg/m^2/s]
    ZROFN     (1:N) => bus( x( ROFN   ,1,1 ) : )         ! output Liquid water runoff from snow pack [kg/m^2/s]
    ZROVG     (1:N) => bus( x( ROVG   ,1,1 ) : )         ! output Liquid/frozen water runoff from vegetation to ground surface [kg/m^2/s]
    ZRUNOFF   (1:N) => bus( x( RUNOFF ,1,1 ) : )         ! output Total runoff from soil column [m]
    ZSCAN     (1:N) => bus( x( IVEG   ,1,1 ) : )         !  inout Intercepted frozen water stored on canopy [kg/m^2]
    ZSUBFLW   (1:N) => bus( x( SUBFLW ,1,1 ) : )         ! output Interflow from sides of soil column [m]
    ZTBASE    (1:N) => bus( x( TBASE  ,1,1 ) : )         !  inout Temperature of bedrock in third soil layer (if only three layers are being modelled) [K]
    ZTCAN     (1:N) => bus( x( TVEG   ,1,1 ) : )         !  inout Vegetation canopy temperature [K]
    ZTCS      (1:N,1:IG) => bus( x( TCS    ,1,1 ) : )    !  inout Thermal conductivity of soil [W/m/K]
    ZTHFC     (1:N,1:IG) => bus( x( THFC   ,1,1 ) : )    !  inout Field capacity [m^3/m^3]
    ZTHLMIN   (1:N,1:IG) => bus( x( THLMIN ,1,1 ) : )    !  inout Residual soil liquid water content remaining after freezing or evaporation [m^3/m^3]
    ZTHLRAT   (1:N,1:IG) => bus( x( THLRAT ,1,1 ) : )    !  inout Fractional saturation of soil behind the wetting front []
    ZTHLRET   (1:N,1:IG) => bus( x( THLRET ,1,1 ) : )    !  inout Liquid water retention capacity for organic soil [m^3/m^3]
    ZTHPOR    (1:N,1:IG) => bus( x( THPOR  ,1,1 ) : )    !  inout Pore volume in soil layer [m^3/m^3]
    ZTPOND    (1:N) => bus( x( TPOND  ,1,1 ) : )         !  inout Temperature of ponded water [C]
    ZTSNOW    (1:N) => bus( x( TSNO   ,1,1 ) : )         !  inout Temperature of the snow pack [C]
    ZTSRAD    (1:N) => bus( x( TSRAD  ,1,1 ) : )         !  inout Diagnosed effective surface black-body temperature [K]
    ZWSNOW    (1:N) => bus( x( WSNOW  ,1,1 ) : )         !  inout Liquid water content of snow pack [kg/m^2]
    ZWTRC     (1:N) => bus( x( WTRC   ,1,1 ) : )         ! output Diagnosed residual water transferred off the vegetation canopy [kg/m^2/s]
    ZWTRG     (1:N) => bus( x( WTRG   ,1,1 ) : )         ! output Diagnosed water transferred into or out of the ice [kg/m^2/s]
    ZWTRS     (1:N) => bus( x( WTRS   ,1,1 ) : )         ! output Diagnosed water transferred into or out of the snow pack [kg/m^2/s]
    ZXSLOPE   (1:N) => bus( x( XSLOPE ,1,1 ) : )         !  input Surface slope (used when running MESH code) [degrees]
    ZZBOTW    (1:N,1:IG) => bus( x( ZBOTW  ,1,1 ) : )    !  inout Depth to permeable bottom of soil layer [m]
    ZZPOND    (1:N) => bus( x( ZPOND  ,1,1 ) : )         !  inout Depth of ponded water on surface [m]
    FSOLUACC  (1:N) => bus( x( FSOLUPAF,1,1 ) : )        !  inout acc. of soil surf. upward solar flux
    FIRUACC   (1:N) => bus( x( FIRUPAF,1,1 ) : )         !  inout acc. of soil surf. upward infrared flux
    ALVIS_SOL (1:N) => bus( x( ALVIS  ,1,1 ) : )         ! output visible surface albedo 
    ZFL       (1:N) => bus( x( FL     ,1,1 ) : )         ! output heat flux at soil surface
       !
    ZSDEPTH   (1:N) => bus( x( SDEPTH ,1,1 ) : )         !  input Permeable depth of soil column (depth to bedrock) [m]
    ZSOILCOL  (1:N) => bus( x( SOILCOL ,1,1 ) : )        !  input soil color for albedo lookup table
    ZFCANMX   (1:N,1:ICP1) => bus( x( FCANMX ,1,1 ) : )  !  input fract. coverage of vegetation classes
    ZZOLN     (1:N,1:ICP1) => bus( x( ZOLN   ,1,1 ) : )  !  input ln of roughness length for each veg. class
    ZALVSC    (1:N,1:ICP1) => bus( x( ALVSC  ,1,1 ) : )  !  input canopy albedo (visible)
    ZALIRC    (1:N,1:ICP1) => bus( x( ALIRC  ,1,1 ) : )  !  input canopy albedo (near i.r.)
    ZPAIMAX   (1:N,1:IC) => bus( x( LAIMAX ,1,1 ) : )    !  input maximum leaf area index (LAI)
    ZPAIMIN   (1:N,1:IC) => bus( x( LAIMIN ,1,1 ) : )    !  input minimum leaf area index (LAI)
    ZCWGTMX   (1:N,1:IC) => bus( x( VEGMA  ,1,1 ) : )    !  input standing mass of canopy
    ZZRTMAX   (1:N,1:IC) => bus( x( ROOTDP ,1,1 ) : )    !  input rooting soil depth
    ZRSMIN    (1:N,1:IC) => bus( x( STOMR  ,1,1 ) : )    !  input minimum stomatal resistance
    ZQA50     (1:N,1:IC) => bus( x( QA50   ,1,1 ) : )    !  input Reference value of incoming shortwave radiation for vegetation category (used in stomatal resistance calculation) [W m-2]
    ZVPDA     (1:N,1:IC) => bus( x( VPDA   ,1,1 ) : )    !  input Vapour pressure deficit coefficient for vegetation category (used in stomatal resistance calculation) []
    ZVPDB     (1:N,1:IC) => bus( x( VPDB   ,1,1 ) : )    !  input Vapour pressure deficit coefficient for vegetation category (used in stomatal resistance calculation) []
    ZPSIGA    (1:N,1:IC) => bus( x( PSIGA  ,1,1 ) : )    !  input Soil moisture suction coefficient for vegetation category (used in stomatal resistance calculation) []
    ZPSIGB    (1:N,1:IC) => bus( x( PSIGB  ,1,1 ) : )    !  input Soil moisture suction coefficient for vegetation category (used in stomatal resistance calculation) []
    ZHUAIRCAN (1:N) => bus( x( HUAIRCAN,1,1 ) : )        !  inout Specific humidity of air within vegetation canopy space [kg/kg]
    ZTAIRCAN  (1:N) => bus( x( TAIRCAN,1,1 ) : )         !  inout Temperature of air within vegetation canopy [K]
    ZTBASFL   (1:N) => bus( x( TBASFL ,1,1 ) : )         ! output Temperature of base flow from bottom of soil column [K]
    ZTOVRFL   (1:N) => bus( x( TOVRFL ,1,1 ) : )         ! output Temperature of overland flow from top of soil column [K]
    ZTSUBFL   (1:N) => bus( x( TSUBFL ,1,1 ) : )         ! output Temperature of interflow from sides of soil column [K]
    ZTRUNOFF  (1:N) => bus( x( TRUNOFF,1,1 ) : )         ! output Temperature of total runoff from soil column [K]
    ZTSFS     (1:N,1:4) => bus( x( TSURFSA,1,1 ) : )     !  inout Ground surface temperature over subarea [K]
    ztt2m     (1:n) => bus( x(tdiagtyp,1,indx_sfc)  : )  ! ourput screen level temperature for each sfc type
!
    ZFTEMP    (1:N) => bus( x( FTEMP  ,1,indx_sfc ) : )  !  inout sfc. temperature flux (T-DT)
    ZFVAP     (1:N) => bus( x( FVAP   ,1,indx_sfc ) : )  !  inout surf. vapour flux (T-DT)
    ZRIB      (1:N) => bus( x( RIB    ,1,1 ) : )         ! output Bulk Richardson number [-10,5]
    ZCDH      (1:N) => bus( x( CDH    ,1,1 ) : )         ! output Surface drag coefficient for heat []
    ZCDM      (1:N) => bus( x( CDM    ,1,1 ) : )         ! output Surface drag coefficient for momentum []

! CTEM fields
    ZGLEAFMAS    (1:N,1:ICC)      => bus( x( GLEAFMAS,1,1 ) : ) !CTEM input
    ZBLEAFMAS    (1:N,1:ICC)      => bus( x( BLEAFMAS,1,1 ) : ) !CTEM input
    ZSTEMMASS    (1:N,1:ICC)      => bus( x( STEMMASS,1,1 ) : ) !CTEM input
    ZROOTMASS    (1:N,1:ICC)      => bus( x( ROOTMASS,1,1 ) : ) !CTEM input
    ZLITRMASS    (1:N,1:ICCP1)    => bus( x( LITRMASS,1,1 ) : ) !CTEM input
    ZSOILCMAS    (1:N,1:ICCP1)    => bus( x( SOILCMAS,1,1 ) : ) !CTEM input
    ZMLIGHTNG    (1:N,1:12)       => bus( x( MLIGHTNG,1,1 ) : ) !CTEM input
    ZPRBFRHUC    (1:N)            => bus( x( PRBFRHUC,1,1 ) : ) !CTEM input
    ZEXTNPROB    (1:N)            => bus( x( EXTNPROB,1,1 ) : ) !CTEM input
    ZLFSTATUR    (1:N,1:ICC)      => bus( x( LFSTATUR,1,1 ) : ) !CTEM input
    ZPANDAYR     (1:N,1:ICC)      => bus( x( PANDAYR ,1,1 ) : ) !CTEM input
    ZCO2CONC     (1:N)            => bus( x( CO2CONC,1,1 ) : ) !including this variable in the mosaic might be unnecessary
    ZCO2I1CS     (1:N,1:ICC)      => bus( x( CO2I1CS,1,1 ) : ) !CTEM
    ZCO2I1CG     (1:N,1:ICC)      => bus( x( CO2I1CG,1,1 ) : ) !CTEM
    ZCO2I2CS     (1:N,1:ICC)      => bus( x( CO2I2CS,1,1 ) : ) !CTEM
    ZCO2I2CG     (1:N,1:ICC)      => bus( x( CO2I2CG,1,1 ) : ) !CTEM
    ZUVACCGAT    (1:N)            => bus( x( UVACCGAT,1,1 ) : ) !CTEM
    ZVVACCGAT    (1:N)            => bus( x( VVACCGAT,1,1 ) : ) !CTEM
    ZFSNOWACC    (1:N)            => bus( x( FSNOWACC,1,1 ) : ) !CTEM
    ZTCANOACCGAT (1:N)            => bus( x( TCANOACCGAT,1,1 ) : ) !CTEM
    ZTCANSACC    (1:N)            => bus( x( TCANSACC,1,1 ) : ) !CTEM
    ZTAACCGAT    (1:N)            => bus( x( TAACCGAT,1,1 ) : ) !CTEM
    ZANCSVGAC    (1:N,1:ICC)      => bus( x( ANCSVGAC,1,1 ) : ) !CTEM
    ZANCGVGAC    (1:N,1:ICC)      => bus( x( ANCGVGAC,1,1 ) : ) !CTEM
    ZRMLCSVGA    (1:N,1:ICC)      => bus( x( RMLCSVGA,1,1 ) : ) !CTEM
    ZRMLCGVGA    (1:N,1:ICC)      => bus( x( RMLCGVGA,1,1 ) : ) !CTEM
    ZTODFRAC     (1:N,1:ICC)      => bus( x( TODFRAC ,1,1 ) : ) !CTEM
    ZTBARACCGAT  (1:N,1:IG)       => bus( x( TBARACCGAT,1,1 ) : ) !CTEM
    ZTBARCACC    (1:N,1:IG)       => bus( x( TBARCACC ,1,1 ) : ) !CTEM
    ZTBARCSACC   (1:N,1:IG)       => bus( x( TBARCSACC,1,1 ) : ) !CTEM
    ZTBARGACC    (1:N,1:IG)       => bus( x( TBARGACC ,1,1 ) : ) !CTEM
    ZTBARGSACC   (1:N,1:IG)       => bus( x( TBARGSACC,1,1 ) : ) !CTEM
    ZTHLIQCACC   (1:N,1:IG)       => bus( x( THLIQCACC,1,1 ) : ) !CTEM
    ZTHLIQGACC   (1:N,1:IG)       => bus( x( THLIQGACC,1,1 ) : ) !CTEM
    ZTHICECACC   (1:N,1:IG)       => bus( x( THICECACC,1,1 ) : ) !CTEM
    ZFCANCMX     (1:N,1:ICC)      => bus( x( FCANCMX,1,1 ) : ) !CTEM
    ZDVDFCAN     (1:N,1:ICC)      => bus( x( DVDFCAN,1,1 ) : ) !CTEM
    ZVGBIOMAS    (1:N)            => bus( x( VGBIOMAS ,1,1 ) : ) !CTEM
    ZGAVGLAI     (1:N)            => bus( x( GAVGLAI  ,1,1 ) : ) !CTEM
    ZGAVGLTMS    (1:N)            => bus( x( GAVGLTMS ,1,1 ) : ) !CTEM
    ZGAVGSCMS    (1:N)            => bus( x( GAVGSCMS ,1,1 ) : ) !CTEM
    ZCOLDDAYR    (1:N,1:2)        => bus( x( COLDDAYR ,1,1 ) : ) !CTEM
    ZGRWTHEFF    (1:N,1:ICC)      => bus( x( GRWTHEFF ,1,1 ) : ) !CTEM
    ZFLHRLOSS    (1:N,1:ICC)      => bus( x( FLHRLOSS ,1,1 ) : ) !CTEM
    ZSTMHRLOS    (1:N,1:ICC)      => bus( x( STMHRLOS ,1,1 ) : ) !CTEM
    ZROTHRLOS    (1:N,1:ICC)      => bus( x( ROTHRLOS ,1,1 ) : ) !CTEM
    ZLYSTMMAS    (1:N,1:ICC)      => bus( x( LYSTMMAS ,1,1 ) : ) !CTEM
    ZLYROTMAS    (1:N,1:ICC)      => bus( x( LYROTMAS ,1,1 ) : ) !CTEM
    ZTYMAXLAI    (1:N,1:ICC)      => bus( x( TYMAXLAI ,1,1 ) : ) !CTEM
    ZAILCG       (1:N,1:ICC)      => bus( x( AILCG    ,1,1 ) : ) !CTEM (BIO2STR output)
    ZAILCB       (1:N,1:ICC)      => bus( x( AILCB    ,1,1 ) : ) !CTEM
    ZAILC        (1:N,1:IC)       => bus( x( AILC     ,1,1 ) : ) !CTEM
    ZZOLNC       (1:N,1:IC)       => bus( x( ZOLNC    ,1,1 ) : ) !CTEM
    ZRMATC       (1:N,1:IC,1:IG)  => bus( x( RMATC    ,1,1 ) : ) !CTEM
    ZRMATCTEM    (1:N,1:ICC,1:IG) => bus( x( RMATCTEM ,1,1 ) : ) !CTEM
    ZSLAI        (1:N,1:ICC)      => bus( x( SLAI     ,1,1 ) : ) !CTEM
    ZBMASVEG     (1:N,1:ICC)      => bus( x( BMASVEG  ,1,1 ) : ) !CTEM
    ZCMASVEGC    (1:N,1:IC)       => bus( x( CMASVEGC ,1,1 ) : ) !CTEM
    ZVEGHGHT     (1:N,1:ICC)      => bus( x( VEGHGHT  ,1,1 ) : ) !CTEM
    ZROOTDPTH    (1:N,1:ICC)      => bus( x( ROOTDPTH ,1,1 ) : ) !CTEM
    ZALVSCTM     (1:N,1:IC)       => bus( x( ALVSCTM  ,1,1 ) : ) !CTEM
    ZALIRCTM     (1:N,1:IC)       => bus( x( ALIRCTM  ,1,1 ) : ) !CTEM
    ZPAIC        (1:N,1:IC)       => bus( x( PAIC     ,1,1 ) : ) !CTEM
    ZSLAIC       (1:N,1:IC)       => bus( x( SLAIC    ,1,1 ) : ) !CTEM
    ZCFLUXCG     (1:N)            => bus( x( CFLUXCG  ,1,1 ) : ) !CTEM (CLASST input)
    ZCFLUXCS     (1:N)            => bus( x( CFLUXCS  ,1,1 ) : )
    ZPFCANCMX    (1:N,1:ICC)      => bus( x( PFCANCMX ,1,1 ) : )
    ZNFCANCMX    (1:N,1:ICC)      => bus( x( NFCANCMX ,1,1 ) : )
    ZNPP         (1:N)            => bus( x( NPP      ,1,1 ) : ) !CTEM output
    ZNEP         (1:N)            => bus( x( NEP      ,1,1 ) : ) !CTEM output
    ZHETRORES    (1:N)            => bus( x( HETRORES ,1,1 ) : ) !CTEM output
    ZAUTORES     (1:N)            => bus( x( AUTORES  ,1,1 ) : ) !CTEM output
    ZSOILRESP    (1:N)            => bus( x( SOILRESP ,1,1 ) : ) !CTEM output
    ZRM          (1:N)            => bus( x( RM       ,1,1 ) : ) !CTEM output
    ZRG          (1:N)            => bus( x( RG       ,1,1 ) : ) !CTEM output
    ZNBP         (1:N)            => bus( x( NBP      ,1,1 ) : ) !CTEM output
    ZLITRES      (1:N)            => bus( x( LITRES   ,1,1 ) : ) !CTEM output
    ZSOCRES      (1:N)            => bus( x( SOCRES   ,1,1 ) : ) !CTEM output
    ZGPP         (1:N)            => bus( x( GPP      ,1,1 ) : ) !CTEM output
    ZDSTCEMLS    (1:N)            => bus( x( DSTCEMLS ,1,1 ) : ) !CTEM output
    ZLITRFALL    (1:N)            => bus( x( LITRFALL ,1,1 ) : ) !CTEM output
    ZHUMIFTRS    (1:N)            => bus( x( HUMIFTRS ,1,1 ) : ) !CTEM output
    ZRML         (1:N)            => bus( x( RML      ,1,1 ) : ) !CTEM output
    ZRMS         (1:N)            => bus( x( RMS      ,1,1 ) : ) !CTEM output
    ZRMR         (1:N)            => bus( x( RMR      ,1,1 ) : ) !CTEM output
    ZTLTRLEAF    (1:N,1:ICC)      => bus( x( TLTRLEAF ,1,1 ) : ) !CTEM output
    ZTLTRSTEM    (1:N,1:ICC)      => bus( x( TLTRSTEM ,1,1 ) : ) !CTEM output
    ZTLTRROOT    (1:N,1:ICC)      => bus( x( TLTRROOT ,1,1 ) : ) !CTEM output
    ZLEAFLITR    (1:N,1:ICC)      => bus( x( LEAFLITR ,1,1 ) : ) !CTEM output
    ZROOTTEMP    (1:N,1:ICC)      => bus( x( ROOTTEMP ,1,1 ) : ) !CTEM output
    ZAFRLEAF     (1:N,1:ICC)      => bus( x( AFRLEAF  ,1,1 ) : ) !CTEM output
    ZAFRSTEM     (1:N,1:ICC)      => bus( x( AFRSTEM  ,1,1 ) : ) !CTEM output
    ZAFRROOT     (1:N,1:ICC)      => bus( x( AFRROOT  ,1,1 ) : ) !CTEM output
    ZWTSTATUS    (1:N,1:ICC)      => bus( x( WTSTATUS ,1,1 ) : ) !CTEM output
    ZLTSTATUS    (1:N,1:ICC)      => bus( x( LTSTATUS ,1,1 ) : ) !CTEM output
    ZBURNAREA    (1:N)            => bus( x( BURNAREA ,1,1 ) : ) !CTEM output
    ZPROBFIRE    (1:N)            => bus( x( PROBFIRE ,1,1 ) : ) !CTEM output
    ZLUCEMCOM    (1:N)            => bus( x( LUCEMCOM ,1,1 ) : ) !CTEM output
    ZLUCLTRIN    (1:N)            => bus( x( LUCLTRIN ,1,1 ) : ) !CTEM output
    ZLUCSOCIN    (1:N)            => bus( x( LUCSOCIN ,1,1 ) : ) !CTEM output
    ZNPPVEG      (1:N,1:ICC)      => bus( x( NPPVEG   ,1,1 ) : ) !CTEM output
    ZGRCLAREA    (1:N)            => bus( x( GRCLAREA ,1,1 ) : ) !CTEM output
    ZDSTCEMLS3   (1:N)            => bus( x( DSTCEMLS3,1,1 ) : ) !CTEM output
    ZRMLVEGACC   (1:N,1:ICC)      => bus( x( RMLVEGACC,1,1 ) : ) !CTEM output
    ZRMSVEG      (1:N,1:ICC)      => bus( x( RMSVEG   ,1,1 ) : ) !CTEM output
    ZRMRVEG      (1:N,1:ICC)      => bus( x( RMRVEG   ,1,1 ) : ) !CTEM output
    ZRGVEG       (1:N,1:ICC)      => bus( x( RGVEG    ,1,1 ) : ) !CTEM output
    ZVGBIOMAS_VEG(1:N,1:ICC)      => bus( x( VGBIOMAS_VEG,1,1 ) : ) !CTEM output
    ZGPPVEG      (1:N,1:ICC)      => bus( x( GPPVEG   ,1,1 ) : ) !CTEM output
    ZNEPVEG      (1:N,1:ICCP1)    => bus( x( NEPVEG   ,1,1 ) : ) !CTEM output

    ZFSINACC   (1:N)         => bus( x( FSINACC  ,1,1 ) : )
    ZFLINACC   (1:N)         => bus( x( FLINACC  ,1,1 ) : )
    ZFLUTACC   (1:N)         => bus( x( FLUTACC  ,1,1 ) : )
    ZALSWACC   (1:N)         => bus( x( ALSWACC  ,1,1 ) : )
    ZALLWACC   (1:N)         => bus( x( ALLWACC  ,1,1 ) : )
    ZPREACC    (1:N)         => bus( x( PREACC   ,1,1 ) : )

    ZTCURM     (1:N)         => bus( x( TCURM    ,1,1 ) : )
    ZSRPCURYR  (1:N)         => bus( x( SRPCURYR ,1,1 ) : )
    ZDFTCURYR  (1:N)         => bus( x( DFTCURYR ,1,1 ) : )
    ZTMONTHB   (1:N,1:12)    => bus( x( TMONTHB  ,1,1 ) : )
    ZANPCPCUR  (1:N)         => bus( x( ANPCPCUR ,1,1 ) : )
    ZANPECUR   (1:N)         => bus( x( ANPECUR  ,1,1 ) : )
    ZGDD5CUR   (1:N)         => bus( x( GDD5CUR  ,1,1 ) : )
    ZSURMNR    (1:N)         => bus( x( SURMNR   ,1,1 ) : )
    ZDEFMNR    (1:N)         => bus( x( DEFMNR   ,1,1 ) : )
    ZSRPLSCUR  (1:N)         => bus( x( SRPLSCUR ,1,1 ) : )
    ZDEFCTCUR  (1:N)         => bus( x( DEFCTCUR ,1,1 ) : )
    ZGEREMORT  (1:N,1:ICC)   => bus( x( GEREMORT ,1,1 ) : )
    ZINTRMORT  (1:N,1:ICC)   => bus( x( INTRMORT ,1,1 ) : )
    ZLAMBDA    (1:N,1:ICC)   => bus( x( LAMBDA   ,1,1 ) : )
    ZPFTEXISTR (1:N,1:ICC)   => bus( x( PFTEXISTR,1,1 ) : )
    ZTWARMM    (1:N)         => bus( x( TWARMM   ,1,1 ) : )
    ZTCOLDM    (1:N)         => bus( x( TCOLDM   ,1,1 ) : )
    ZGDD5      (1:N)         => bus( x( GDD5     ,1,1 ) : )
    ZARIDITY   (1:N)         => bus( x( ARIDITY  ,1,1 ) : )
    ZSRPLSMON  (1:N)         => bus( x( SRPLSMON ,1,1 ) : )
    ZDEFCTMON  (1:N)         => bus( x( DEFCTMON ,1,1 ) : )
    ZANNDEFCT  (1:N)         => bus( x( ANNDEFCT ,1,1 ) : )
    ZANNSRPLS  (1:N)         => bus( x( ANNSRPLS ,1,1 ) : )
    ZANNPCP    (1:N)         => bus( x( ANNPCP   ,1,1 ) : )
    ZANPOTEVP  (1:N)         => bus( x( ANPOTEVP ,1,1 ) : )
    ZDRYSLEN   (1:N)         => bus( x( DRYSLEN  ,1,1 ) : )
    ZWDMINDEX  (1:N,1:12)    => bus( x( WDMINDEX ,1,1 ) : )
    ZBURNVEGF  (1:N,1:ICC)   => bus( x( BURNVEGF ,1,1 ) : )
    ZPSTEMMASS (1:N,1:ICC)   => bus( x( PSTEMMASS,1,1 ) : )
    ZPGLEAFMASS(1:N,1:ICC)   => bus( x( PGLEAFMASS,1,1 ) : )
    ZCOLRATE   (1:N,1:ICC)   => bus( x( COLRATE  ,1,1 ) : )
    ZMORTRATE  (1:N,1:ICC)   => bus( x( MORTRATE ,1,1 ) : )
    ZAUTRESVEG (1:N,1:ICC)   => bus( x( AUTRESVEG,1,1 ) : )
    ZHETRESVEG (1:N,1:ICCP1) => bus( x( HETRESVEG,1,1 ) : )
    ZLITRESVEG (1:N,1:ICCP1) => bus( x( LITRESVEG,1,1 ) : )
    ZNBPVEG    (1:N,1:ICCP1) => bus( x( NBPVEG   ,1,1 ) : )
    ZSOCRESVEG (1:N,1:ICCP1) => bus( x( SOCRESVEG,1,1 ) : )
! End CTEM fields

if (prints_L) then
print*,'class_main'
print*,'class_main Geophysical fields'
do j=1,IG
  print*,'class_main ZSAND:',j,minval(ZSAND(:,j)),maxval(ZSAND(:,j)),sum(ZSAND(:,j))/(N)
  print*,'class_main ZCLAY:',j,minval(ZCLAY(:,j)),maxval(ZCLAY(:,j)),sum(ZCLAY(:,j))/(N)
  print*,'class_main ZORGM:',j,minval(ZORGM(:,j)),maxval(ZORGM(:,j)),sum(ZORGM(:,j))/(N)
enddo
print*,'class_main ZXDRAIN        :',minval(ZXDRAIN),maxval(ZXDRAIN),sum(ZXDRAIN)/(N)

print*,'class_main'
print*,'class_main Physics input'
print*,'class_main QA             :',minval(QA),maxval(QA),sum(QA)/(N)
print*,'class_main PS             :',minval(PS),maxval(PS),sum(PS)/(N)
print*,'class_main TA             :',minval(TA),maxval(TA),sum(TA)/(N)
print*,'class_main UA             :',minval(UA),maxval(UA),sum(UA)/(N)
print*,'class_main VA             :',minval(VA),maxval(VA),sum(VA)/(N)
print*,'class_main ZCANG          :',minval(ZCANG),maxval(ZCANG),sum(ZCANG)/(N)
print*,'class_main ZDLAT          :',minval(ZDLAT),maxval(ZDLAT),sum(ZDLAT)/(N)
print*,'class_main ZDLON          :',minval(ZDLON),maxval(ZDLON),sum(ZDLON)/(N)
print*,'class_main ZZTSL          :',minval(ZZTSL),maxval(ZZTSL),sum(ZZTSL)/(N)
print*,'class_main ZZUSL          :',minval(ZZUSL),maxval(ZZUSL),sum(ZZUSL)/(N)
print*,'class_main FLUSOL         :',minval(FLUSOL),maxval(FLUSOL),sum(FLUSOL)/(N)
print*,'class_main QLWIN          :',minval(QLWIN),maxval(QLWIN),sum(QLWIN)/(N)
print*,'class_main ZRAINRATE      :',minval(ZRAINRATE),maxval(ZRAINRATE),sum(ZRAINRATE)/(N)
print*,'class_main ZSNOWRATE      :',minval(ZSNOWRATE),maxval(ZSNOWRATE),sum(ZSNOWRATE)/(N)
do j=1,NBS
  print*,'class_main ZFSSB:',j,minval(ZFSSB(:,j)),maxval(ZFSSB(:,j)),sum(ZFSSB(:,j))/(N)
  print*,'class_main ZFSDB:',j,minval(ZFSDB(:,j)),maxval(ZFSDB(:,j)),sum(ZFSDB(:,j))/(N)
  print*,'class_main ZFSFB:',j,minval(ZFSFB(:,j)),maxval(ZFSFB(:,j)),sum(ZFSFB(:,j))/(N)
enddo
!print*,'class_main ZGRKFAC        :',minval(ZGRKFAC),maxval(ZGRKFAC),sum(ZGRKFAC)/(N)
!print*,'class_main ZWFSURF        :',minval(ZWFSURF),maxval(ZWFSURF),sum(ZWFSURF)/(N)
!print*,'class_main ZWFCINT        :',minval(ZWFCINT),maxval(ZWFCINT),sum(ZWFCINT)/(N)
!print*,'class_main ZRAINRATE      :',minval(ZRAINRATE),maxval(ZRAINRATE),sum(ZRAINRATE)/(N)
!print*,'class_main ZSNOWRATE      :',minval(ZSNOWRATE),maxval(ZSNOWRATE),sum(ZSNOWRATE)/(N)
!print*,'class_main ZSNOW          :',minval(ZSNOW),maxval(ZSNOW),sum(ZSNOW)/(N)
!print*,'class_main REFSNO         :',minval(REFSNO),maxval(REFSNO),sum(REFSNO)/(N)
!print*,'class_main groundHeatFlux :',minval(groundHeatFlux),maxval(groundHeatFlux),sum(groundHeatFlux)/(N)
!print*,'class_main ZWFSURF        :',minval(ZWFSURF),maxval(ZWFSURF),sum(ZWFSURF)/(N)
!print*,'class_main Z0H            :',minval(Z0H),maxval(Z0H),sum(Z0H)/(N)
!print*,'class_main Z0ORO          :',minval(Z0ORO),maxval(Z0ORO),sum(Z0ORO)/(N)
!print*,'class_main QFLUX          :',minval(QFLUX),maxval(QFLUX),sum(QFLUX)/(N)
!print*,'class_main TFLUX          :',minval(TFLUX),maxval(TFLUX),sum(TFLUX)/(N)
!print*,'class_main ZCFLUX         :',minval(ZCFLUX),maxval(ZCFLUX),sum(ZCFLUX)/(N)
!print*,'class_main ZBM            :',minval(ZBM),maxval(ZBM),sum(ZBM)/(N)
!print*,'class_main EVAPO          :',minval(EVAPO),maxval(EVAPO),sum(EVAPO)/(N)
!print*,'class_main FFC            :',minval(FFC),maxval(FFC),sum(FFC)/(N)
!print*,'class_main QSENS          :',minval(QSENS),maxval(QSENS),sum(QSENS)/(N)
!print*,'class_main FCS            :',minval(FCS),maxval(FCS),sum(FCS)/(N)
!print*,'class_main FG             :',minval(FG),maxval(FG),sum(FG)/(N)
!print*,'class_main FGS            :',minval(FGS),maxval(FGS),sum(FGS)/(N)
!print*,'class_main ZFRV           :',minval(ZFRV),maxval(ZFRV),sum(ZFRV)/(N)
!print*,'class_main QEVAP          :',minval(QEVAP),maxval(QEVAP),sum(QEVAP)/(N)
!print*,'class_main HBL            :',minval(HBL),maxval(HBL),sum(HBL)/(N)
!print*,'class_main ZILMO          :',minval(ZILMO),maxval(ZILMO),sum(ZILMO)/(N)
!print*,'class_main THICE          :',minval(THICE),maxval(THICE),sum(THICE)/(N*IG)
!print*,'class_main SQ             :',minval(SQ),maxval(SQ),sum(SQ)/(N)
!print*,'class_main QS             :',minval(QS),maxval(QS),sum(QS)/(N)
!print*,'class_main ST             :',minval(ST),maxval(ST),sum(ST)/(N)
!print*,'class_main TS             :',minval(TS),maxval(TS),sum(TS)/(N*IG)
!print*,'class_main SU             :',minval(SU),maxval(SU),sum(SU)/(N)
!print*,'class_main SV             :',minval(SV),maxval(SV),sum(SV)/(N)
!print*,'class_main THLIQ          :',minval(THLIQ),maxval(THLIQ),sum(THLIQ)/(N*IG)
!print*,'class_main XSNO           :',minval(XSNO),maxval(XSNO),sum(XSNO)/(N)
!print*,'class_main ZALBSNO        :',minval(ZALBSNO),maxval(ZALBSNO),sum(ZALBSNO)/(N)
!print*,'class_main ZALGDN         :',minval(ZALGDN),maxval(ZALGDN),sum(ZALGDN)/(N)
!print*,'class_main ZALGDV         :',minval(ZALGDV),maxval(ZALGDV),sum(ZALGDV)/(N)
!print*,'class_main ZALGWN         :',minval(ZALGWN),maxval(ZALGWN),sum(ZALGWN)/(N)
!print*,'class_main ZALGWV         :',minval(ZALGWV),maxval(ZALGWV),sum(ZALGWV)/(N)
!print*,'class_main ZBASFLW        :',minval(ZBASFLW),maxval(ZBASFLW),sum(ZBASFLW)/(N)
!print*,'class_main ZBI            :',minval(ZBI),maxval(ZBI),sum(ZBI)/(N*IG)
!print*,'class_main ZCMAI          :',minval(ZCMAI),maxval(ZCMAI),sum(ZCMAI)/(N)
!print*,'class_main ZDELZW         :',minval(ZDELZW),maxval(ZDELZW),sum(ZDELZW)/(N*IG)
!print*,'class_main ZFLGG          :',minval(ZFLGG),maxval(ZFLGG),sum(ZFLGG)/(N)
!print*,'class_main ZFLGS          :',minval(ZFLGS),maxval(ZFLGS),sum(ZFLGS)/(N)
!print*,'class_main ZFLGV          :',minval(ZFLGV),maxval(ZFLGV),sum(ZFLGV)/(N)
!print*,'class_main ZFSGG          :',minval(ZFSGG),maxval(ZFSGG),sum(ZFSGG)/(N)
!print*,'class_main ZFSGS          :',minval(ZFSGS),maxval(ZFSGS),sum(ZFSGS)/(N)
!print*,'class_main ZFSGV          :',minval(ZFSGV),maxval(ZFSGV),sum(ZFSGV)/(N)
!print*,'class_main ZFSNOW         :',minval(ZFSNOW),maxval(ZFSNOW),sum(ZFSNOW)/(N)
!print*,'class_main ZGRKSAT        :',minval(ZGRKSAT),maxval(ZGRKSAT),sum(ZGRKSAT)/(N*IG)
!print*,'class_main ZGROWTH        :',minval(ZGROWTH),maxval(ZGROWTH),sum(ZGROWTH)/(N)
!print*,'class_main ZHCPS          :',minval(ZHCPS),maxval(ZHCPS),sum(ZHCPS)/(N*IG)
!print*,'class_main ZHEVC          :',minval(ZHEVC),maxval(ZHEVC),sum(ZHEVC)/(N)
!print*,'class_main ZHEVG          :',minval(ZHEVG),maxval(ZHEVG),sum(ZHEVG)/(N)
!print*,'class_main ZHEVS          :',minval(ZHEVS),maxval(ZHEVS),sum(ZHEVS)/(N)
!print*,'class_main ZHFSC          :',minval(ZHFSC),maxval(ZHFSC),sum(ZHFSC)/(N)
!print*,'class_main ZHFSG          :',minval(ZHFSG),maxval(ZHFSG),sum(ZHFSG)/(N)
!print*,'class_main ZHFSS          :',minval(ZHFSS),maxval(ZHFSS),sum(ZHFSS)/(N)
!print*,'class_main ZHMFC          :',minval(ZHMFC),maxval(ZHMFC),sum(ZHMFC)/(N)
!print*,'class_main ZHMFG          :',minval(ZHMFG),maxval(ZHMFG),sum(ZHMFG)/(N*IG)
!print*,'class_main ZHMFN          :',minval(ZHMFN),maxval(ZHMFN),sum(ZHMFN)/(N)
!print*,'class_main ZHTC           :',minval(ZHTC),maxval(ZHTC),sum(ZHTC)/(N*IG)
!print*,'class_main ZHTCC          :',minval(ZHTCC),maxval(ZHTCC),sum(ZHTCC)/(N)
!print*,'class_main ZHTCS          :',minval(ZHTCS),maxval(ZHTCS),sum(ZHTCS)/(N)
!print*,'class_main ZOVRFLW        :',minval(ZOVRFLW),maxval(ZOVRFLW),sum(ZOVRFLW)/(N)
!print*,'class_main ZPCFC          :',minval(ZPCFC),maxval(ZPCFC),sum(ZPCFC)/(N)
!print*,'class_main ZPCLC          :',minval(ZPCLC),maxval(ZPCLC),sum(ZPCLC)/(N)
!print*,'class_main ZPCPG          :',minval(ZPCPG),maxval(ZPCPG),sum(ZPCPG)/(N)
!print*,'class_main ZPCPN          :',minval(ZPCPN),maxval(ZPCPN),sum(ZPCPN)/(N)
!print*,'class_main ZPSISAT        :',minval(ZPSISAT),maxval(ZPSISAT),sum(ZPSISAT)/(N*IG)
!print*,'class_main ZPSIWLT        :',minval(ZPSIWLT),maxval(ZPSIWLT),sum(ZPSIWLT)/(N*IG)
!print*,'class_main ZQFC           :',minval(ZQFC),maxval(ZQFC),sum(ZQFC)/(N*IG)
!print*,'class_main ZQFCF          :',minval(ZQFCF),maxval(ZQFCF),sum(ZQFCF)/(N)
!print*,'class_main ZQFCL          :',minval(ZQFCL),maxval(ZQFCL),sum(ZQFCL)/(N)
!print*,'class_main ZQFG           :',minval(ZQFG),maxval(ZQFG),sum(ZQFG)/(N)
!print*,'class_main ZQFN           :',minval(ZQFN),maxval(ZQFN),sum(ZQFN)/(N)
!print*,'class_main ZRCAN          :',minval(ZRCAN),maxval(ZRCAN),sum(ZRCAN)/(N)
!print*,'class_main ZRHOSNO        :',minval(ZRHOSNO),maxval(ZRHOSNO),sum(ZRHOSNO)/(N)
!print*,'class_main ZROFC          :',minval(ZROFC),maxval(ZROFC),sum(ZROFC)/(N)
!print*,'class_main ZROFN          :',minval(ZROFN),maxval(ZROFN),sum(ZROFN)/(N)
!print*,'class_main ZROVG          :',minval(ZROVG),maxval(ZROVG),sum(ZROVG)/(N)
!print*,'class_main ZRUNOFF        :',minval(ZRUNOFF),maxval(ZRUNOFF),sum(ZRUNOFF)/(N)
!print*,'class_main ZSCAN          :',minval(ZSCAN),maxval(ZSCAN),sum(ZSCAN)/(N)
!print*,'class_main ZSUBFLW        :',minval(ZSUBFLW),maxval(ZSUBFLW),sum(ZSUBFLW)/(N)
!print*,'class_main ZTBASE         :',minval(ZTBASE),maxval(ZTBASE),sum(ZTBASE)/(N)
!print*,'class_main ZTCAN          :',minval(ZTCAN),maxval(ZTCAN),sum(ZTCAN)/(N)
!print*,'class_main ZTCS           :',minval(ZTCS),maxval(ZTCS),sum(ZTCS)/(N*IG)
!print*,'class_main ZTHFC          :',minval(ZTHFC),maxval(ZTHFC),sum(ZTHFC)/(N*IG)
!print*,'class_main ZTHLMIN        :',minval(ZTHLMIN),maxval(ZTHLMIN),sum(ZTHLMIN)/(N*IG)
!print*,'class_main ZTHLRAT        :',minval(ZTHLRAT),maxval(ZTHLRAT),sum(ZTHLRAT)/(N*IG)
!print*,'class_main ZTHLRET        :',minval(ZTHLRET),maxval(ZTHLRET),sum(ZTHLRET)/(N*IG)
!print*,'class_main ZTHPOR         :',minval(ZTHPOR),maxval(ZTHPOR),sum(ZTHPOR)/(N*IG)
!print*,'class_main ZTPOND         :',minval(ZTPOND),maxval(ZTPOND),sum(ZTPOND)/(N)
!print*,'class_main ZTSNOW         :',minval(ZTSNOW),maxval(ZTSNOW),sum(ZTSNOW)/(N)
!print*,'class_main ZTSRAD         :',minval(ZTSRAD),maxval(ZTSRAD),sum(ZTSRAD)/(N)
!print*,'class_main ZTSURF         :',minval(ZTSURF),maxval(ZTSURF),sum(ZTSURF)/(N)
!print*,'class_main ZWSNOW         :',minval(ZWSNOW),maxval(ZWSNOW),sum(ZWSNOW)/(N)
!print*,'class_main ZWTRC          :',minval(ZWTRC),maxval(ZWTRC),sum(ZWTRC)/(N)
!print*,'class_main ZWTRG          :',minval(ZWTRG),maxval(ZWTRG),sum(ZWTRG)/(N)
!print*,'class_main ZWTRS          :',minval(ZWTRS),maxval(ZWTRS),sum(ZWTRS)/(N)
!print*,'class_main ZXSLOPE        :',minval(ZXSLOPE),maxval(ZXSLOPE),sum(ZXSLOPE)/(N)
!print*,'class_main ZZBOTW         :',minval(ZZBOTW),maxval(ZZBOTW),sum(ZZBOTW)/(N*IG)
!print*,'class_main ZZPOND         :',minval(ZZPOND),maxval(ZZPOND),sum(ZZPOND)/(N)
!print*,'class_main FSOLUACC       :',minval(FSOLUACC),maxval(FSOLUACC),sum(FSOLUACC)/(N)
!print*,'class_main FIRUACC        :',minval(FIRUACC),maxval(FIRUACC),sum(FIRUACC)/(N)
!print*,'class_main ALVIS_SOL      :',minval(ALVIS_SOL),maxval(ALVIS_SOL),sum(ALVIS_SOL)/(N)
!print*,'class_main ZFL            :',minval(ZFL),maxval(ZFL),sum(ZFL)/(N)
!print*,'class_main ZSDEPTH        :',minval(ZSDEPTH),maxval(ZSDEPTH),sum(ZSDEPTH)/(N)
!print*,'class_main ZSOILCOL       :',minval(ZSOILCOL),maxval(ZSOILCOL),sum(ZSOILCOL)/(N)
!print*,'class_main ZFCANMX        :',minval(ZFCANMX),maxval(ZFCANMX),sum(ZFCANMX)/(N*ICP1)
!print*,'class_main ZZOLN          :',minval(ZZOLN),maxval(ZZOLN),sum(ZZOLN)/(N*ICP1)
!print*,'class_main ZALVSC         :',minval(ZALVSC),maxval(ZALVSC),sum(ZALVSC)/(N*ICP1)
!print*,'class_main ZALIRC         :',minval(ZALIRC),maxval(ZALIRC),sum(ZALIRC)/(N*ICP1)
!print*,'class_main ZPAIMAX        :',minval(ZPAIMAX),maxval(ZPAIMAX),sum(ZPAIMAX)/(N*IC)
!print*,'class_main ZPAIMIN        :',minval(ZPAIMIN),maxval(ZPAIMIN),sum(ZPAIMIN)/(N*IC)
!print*,'class_main ZCWGTMX        :',minval(ZCWGTMX),maxval(ZCWGTMX),sum(ZCWGTMX)/(N*IC)
!print*,'class_main ZZRTMAX        :',minval(ZZRTMAX),maxval(ZZRTMAX),sum(ZZRTMAX)/(N*IC)
!print*,'class_main ZRSMIN         :',minval(ZRSMIN),maxval(ZRSMIN),sum(ZRSMIN)/(N*IC)
!print*,'class_main ZQA50          :',minval(ZQA50),maxval(ZQA50),sum(ZQA50)/(N*IC)
!print*,'class_main ZVPDA          :',minval(ZVPDA),maxval(ZVPDA),sum(ZVPDA)/(N*IC)
!print*,'class_main ZVPDB          :',minval(ZVPDB),maxval(ZVPDB),sum(ZVPDB)/(N*IC)
!print*,'class_main ZPSIGA         :',minval(ZPSIGA),maxval(ZPSIGA),sum(ZPSIGA)/(N*IC)
!print*,'class_main ZPSIGB         :',minval(ZPSIGB),maxval(ZPSIGB),sum(ZPSIGB)/(N*IC)
!print*,'class_main ZHUAIRCAN      :',minval(ZHUAIRCAN),maxval(ZHUAIRCAN),sum(ZHUAIRCAN)/(N)
!print*,'class_main ZTAIRCAN       :',minval(ZTAIRCAN),maxval(ZTAIRCAN),sum(ZTAIRCAN)/(N)
!print*,'class_main ZTBASFL        :',minval(ZTBASFL),maxval(ZTBASFL),sum(ZTBASFL)/(N)
!print*,'class_main ZTOVRFL        :',minval(ZTOVRFL),maxval(ZTOVRFL),sum(ZTOVRFL)/(N)
!print*,'class_main ZTSUBFL        :',minval(ZTSUBFL),maxval(ZTSUBFL),sum(ZTSUBFL)/(N)
!print*,'class_main ZTRUNOFF       :',minval(ZTRUNOFF),maxval(ZTRUNOFF),sum(ZTRUNOFF)/(N)
!print*,'class_main ZTSFS          :',minval(ZTSFS),maxval(ZTSFS),sum(ZTSFS)/(N*4)
!print*,'class_main ztt2m          :',minval(ztt2m),maxval(ztt2m),sum(ztt2m)/(n)
!print*,'class_main ZFTEMP         :',minval(ZFTEMP),maxval(ZFTEMP),sum(ZFTEMP)/(N)
!print*,'class_main ZFVAP          :',minval(ZFVAP),maxval(ZFVAP),sum(ZFVAP)/(N)
!print*,'class_main ZRIB           :',minval(ZRIB),maxval(ZRIB),sum(ZRIB)/(N)
!print*,'class_main ZCDH           :',minval(ZCDH),maxval(ZCDH),sum(ZCDH)/(N)
!print*,'class_main ZCDM           :',minval(ZCDM),maxval(ZCDM),sum(ZCDM)/(N)

!stop
endif


!*************************************************************
!     CLASS'S EXECUTION OPTION AND COMMENTS
!     (FROM RICHARD'S SWITCHES.F90)
!
!    IDISP (CLASSA, APREP)              DEFAULT=0   [R=1]
! If idisp=0, vegetation displacement heights are ignored,
!    because the atmospheric model considers these to be part
!    of the "terrain".
! If idisp=1, vegetation displacement heights are calculated.
!
!    IZREF (CLASSA, APREP, CLASST)      DEFAULT=2   [R=1]
! If izref=1, the bottom of the atmospheric model is taken
!    to lie at the ground surface.
! If izref=2, the bottom of the atmospheric model is taken
!    to lie at the local roughness height.
!***
! N.B. COMBINATION EITHER IDISP=0 AND IZREF=2 OR
!                         IDISP=1 AND IZREF=1
!***
!    ISLFD (CLASST, TSOLC, TSOLVE)      DEFAULT=2   [R=2]
! If islfd=0, drcoef is called for surface stability corrections
!    and the original gcm set of screen-level diagnostic calculations
!    is done.
! If islfd=1, drcoef is called for surface stability corrections
!    and sldiag is called for screen-level diagnostic calculations.
! If islfd=2, flxsurfz is called for surface stability corrections
!    and diasurf is called for screen-level diagnostic calculations.
!
!    ITC   (CLASST, TSOLVC)             DEFAULT=1   [R=2]
!    ITCG  (CLASST, TSOLVC)             DEFAULT=1   [R=2]
!    ITG   (CLASST, TSOLVE)             DEFAULT=1   [R=2]
! itc, itcg and itg are switches to choose the iteration scheme to
! be used in calculating the canopy or ground surface temperature
! respectively.
! If the switch is set to 1, a combination of secant and bisection
! methods is used;
! If to 2, the newton-raphson method is used.
!
!    ILW   (CLASST, TSOLVC, TSOLVE)     DEFAULT=1   [R=?]
! If ilw=1, QLWIN is the incoming longwave
! If ilw=2, QLWIN is the      net longwave
! This option is no longer used (LD)
!
!    IWF   (CLASSA, CLASSW)             DEFAULT=0   [R=0]
! OPTION DEFINED IN option.cdk & SET IN phy_opt.ftn
! If iwf=0, only overland flow and baseflow are modelled, and
!    the ground surface slope is not modelled.
! If iwf=n (0<n<4), the watflood calculations of overland flow
!     and interflow are performed; interflow is drawn from the top
!     n soil layers.
!
!    IPAI  (CLASSA)                     DEFAULT=0   [R=0]
!    IHGT  (CLASSA)                     DEFAULT=0   [R=0]
!    IALC  (CLASSA)                     DEFAULT=0   [R=0]
!    IALS  (CLASSA)                     DEFAULT=0   [R=0]
!    IALG  (CLASSA)                     DEFAULT=0   [R=0]
! If =0 the values of leaf are index, vegetation height, canopy albedo,
!     snow albedo and soil albedo respectively calculated by class are used.
! If =1 the value of the corresponding parameter calculated by class is
!      overridden by a user-supplied input value.
!
!    IPCP  (CLASSI)                     DEFAULT=4   [R=3]
! If ipcp=1, the rainfall-snowfall cutoff is taken to lie at 0 c.
! If ipcp=2, a linear partitioning of precipitation betweeen
!       rainfall and snowfall is done between 0 c and 2 c.
! If ipcp=3, rainfall and snowfall are partitioned according to
!       a polynomial curve between 0 c and 6 c.
! If ipcp=4, ONLY IN GEM. Calculates PCPR
!
    IDISP = 0
    IZREF = 2
    ISLFD = 2
    ITC   = 1
    ITCG  = 1
    ITG   = 1
!   ILW   = 1
    NMIM  = 1

    IPAI  = 0
    IHGT  = 0
    IALC  = 0
    IALS  = 0
    IALG  = 0
    IPCP  = 4

    isnoalb = 0    ! New in CLASSIC (KW)
                   ! 0: original two-band snow albedo algorithms are used. 
                   ! 1: the new four-band routines are used.
                   ! At present, the four band algorithm should NOT be used offline.
    REFSNO = 0.001 ! Only for now!!! Later it needs to get read in!!! ToDo


! Convert soil color index to integer. New in CLASSIC (KW)
  DO I=1,N
     SOCI(I) = NINT(ZSOILCOL(I))
  ENDDO



!*************************************************************
!
!
! Initializing variables that are passed to the CTEM subroutine, but not used
! for the moment. These should go inside the KOUNT_EQ_0 loop once the features
! using them are implemented (LD)
         do i=1,n
           do j=1,8
            wetfrac_sgrd(i,j) = 0.   !setting this to 0 until wetlands are implemented (LD)
           enddo
         enddo
!
!
!    kount0 = (KOUNT.EQ.0 .and. FLUVERT.ne.'SURFACE') .or. &
!             (KOUNT.eq.1 .and. FLUVERT.eq.'SURFACE')
    kount0 = (KOUNT.EQ.0)

!    KOUNT_EQ_0 : IF(KOUNT.EQ.0 .or. maxval(zinisoil).gt.0.5) THEN
    KOUNT_EQ_0 : IF(kount0) THEN
!
!
! Check if ROOTDP less or equal than SDEPTH
!
        do i=1,N
           do j=1,ic
              ZZRTMAX(I,J) = min(ZZRTMAX(I,J), ZSDEPTH(I))
           enddo
        enddo
!
!       Initialize the soil characteristics
!       using the soil texture
!
!print *,'class_main AAA: ISAND(1,:):', ISAND(1,:)
        call soilProperties(ZTHPOR, ZTHLRET, ZTHLMIN, ZBI, ZPSISAT, ZGRKSAT, & ! Formerly CLASSB
                            ZTHLRAT, ZHCPS, ZTCS, ZTHFC, ZTHLW, ZPSIWLT,     &
                            ZDELZW, ZZBOTW, ZALGWV, ZALGWN, ZALGDV, ZALGDN,  &
                            ZSAND, ZCLAY, ZORGM, SOCI, DELZ, ZBOT, ZSDEPTH,  &
                            ISAND, igdr, N, NMIM, 1, N, NMIM, IG, ipeatland)

        ! Save negative ISAND indo ZSAND for following time steps.
        do k=1,IG
           do i=1,N
              if (ISAND(i,k) < 0) ZSAND(i,k) = ISAND(i,k)
           enddo
        enddo

!do i=1,N
!  print *,'class_main: i,ZTHPOR:',ZTHPOR(i,:)
!enddo

!       Set THLIQ, THICE according to ISAND
        do J=1,IG
           do I=1,N
              if (ISAND(i,j) == -3) then   ! Bedrock
                 THLIQ(i,j) = 0. ; THICE(i,j) = 0.
              endif
           enddo
        enddo

!

!print *,'class_main BBB: ISAND(1,:):', ISAND(1,:)

!CC     This statement should insure that the changes to
!CC     ISAND in CLASSB are saved. It does not seem to do
!CC     do so. The loop after this IF block corrects this
!CC     ZSAND = ISAND
           do I=1,N
              ztsurf(i)=ta(i)
              qs(i)    =qa(i)
! Initialize air temperature and specific humidity inside canopy
! from air temperature and specific humidity over the grid point
              ztaircan(i)=ta(i)
              zhuaircan(i)=qa(i)
! Initialize surface temperature for each subareas
! Surface temperature over snow-covered subareas set to zero
              do j=1,2
                 ztsfs(i,j)=TCDK
              enddo
! Surface temperature over snow-free areas
! set to temperature of first soil layer
              do j=3,4
                 ztsfs(i,j)=TS(i,1)
              enddo
           enddo

    ENDIF KOUNT_EQ_0



    do I=1,N
       Isat(i)=ig
       Isats(i)=ig
    enddo


!   Prescribe Black Carbon when not read in
!   New for CLASSIC (KW)
!    do I=1,N
!       BCSNO(I) = 0.25 - 0.2*ct * abs(ZDLAT(i))   ! Parameters need to get adjusted
!                                                  ! Height should get taken into account as well
!    enddo
    BCSNO = 0.

!
!
    doprec = .false.
!
!print *,'class_main CCC: ISAND(1,:):', ISAND(1,:)
    CLASS_MAIN_LOOP : do k=1,n+1
!
      dotile = (k.ne.n+1)
      if (dotile.and.doprec) then
         j = k
      elseif (dotile.and..not.doprec) then
         j = k
         i = k
      elseif (.not.dotile.and.doprec) then

!if (trnch==50) then
!  print *,'class_main: TBAR(28,:):',TS(28,:)
!endif
         CALL atmosphericVarsCalc(VPD, TADP, PADRY, RHOAIR, RHOSNI, RPCP, TRPCP, & ! Formerly CLASSI
                                  SPCP, TSPCP, TA, QA, PCPR, RRATE, SRATE,       &
                                  PS, IPCP, N, I, J)


      ! ToDo: Add call to energyWaterBalanceCheck (0,...)
      !> energyWaterBalanceCheck does the initial calculations for the energy and water balance checks
!      call energyWaterBalanceCheck(0, CTVSTP, CTSSTP, CT1STP, CT2STP, CT3STP, ...)  ! Formerly CLASSZ


!
!print *,'class_main: kount,trnch,TBAR(31,:):',kount,trnch,TS(31,:)
!print *,'class_main 000: trnch,ZSNOW(:):', trnch,ZSNOW(:)
!print *,'class_main 000: trnch,XSNO(:):', trnch,XSNO(:)
!print *,'class_main 000: trnch,ZRHOSNO(:):', trnch,ZRHOSNO(:)
         CALL  radiationDriver(FFC, FG, FCS, FGS, ALVSCN, ALIRCN, & ! Formerly CLASSA
                               ALVSG,  ALIRG,  ALVSCS, ALIRCS, ALVSSN, ALIRSN, &
                               ALVSGC, ALIRGC, ALVSSC, ALIRSC, TRVSCN, TRIRCN, &
                               TRVSCS, TRIRCS, FSVF,   FSVFS, &
                               RAICAN, RAICNS, SNOCAN, SNOCNS, FRAINC, FSNOWC, &
                               FRAICS, FSNOCS, DISP,   DISPS,  ZOMLNC, ZOMLCS, &
                               ZOELNC, ZOELCS, ZOMLNG, ZOMLNS, ZOELNG, ZOELNS, &
                               CHCAP,  CHCAPS, CMASSC, CMASCS, CWLCAP, CWFCAP, &
                               CWLCPS, CWFCPS, RC,     RCS,    RBCOEF, FROOT, &
                               FROOTS, ZPLIMC, ZPLIMG, ZPLMCS, ZPLMGS, ZSNOW, &
                               ZWSNOW, ZALVS,  ZALIR,  ZHTCC,  ZHTCS,  ZHTC, &
                               ALTG,   ALSNO,  TRSNOWC,TRSNOWG, &
                               ZWTRC,  ZWTRS,  ZWTRG,  ZCMAI,  ZFSNOW, &
                               ZFCANMX,ZZOLN,  ZALVSC, ZALIRC, ZPAIMAX,ZPAIMIN, &
                               ZCWGTMX,ZZRTMAX,ZRSMIN, ZQA50,  ZVPDA,  ZVPDB, &
                               ZPSIGA, ZPSIGB, PAIDAT, HGTDAT, ACVDAT, ACIDAT, &
                               ASVDAT, ASIDAT, AGVDAT, AGIDAT, &
                               ZALGWV, ZALGWN, ZALGDV, ZALGDN, &
                               THLIQ,  THICE,  TS,     ZRCAN,  ZSCAN,  ZTCAN, &
                               ZGROWTH,XSNO,   ZTSNOW, ZRHOSNO,ZALBSNO,ZBLEND, &
                               Z0ORO,  SNOLIM, ZPLMG0, ZPLMS0, &
                               FCLOUD, TA,     VPD,    RHOAIR, COSZS, &
                               ZFSDB, ZFSFB, REFSNO, BCSNO, &
                               QSWINV, ZDLAT,  ZDLON,  RHOSNI, DELZ,   ZDELZW, &
                               ZZBOTW, ZTHPOR, ZTHLMIN,ZPSISAT,ZBI,    ZPSIWLT, &
                               ZHCPS,  ISAND, &
                               ZFCANCMX,ICC,ctem_on,ZRMATC,ZZOLNC,ZCMASVEGC, &
                               ZAILC,  ZPAIC,  NOL2PFTS,ZSLAIC, &
                               ZAILCG, ZAILCGS,ZFCANC, ZFCANCS, &
                               IDAY,   N,      I,      J,      NBS, &
                               TRNCH,kount,  IC,     ICP1,   IG,     IDISP,  IZREF, &
                               IWF,    IPAI,   IHGT,   IALC,   IALS,   IALG, &
                               ISNOALB,ZALVSCTM, ZALIRCTM, ipeatland)
!
!
!print *,'class_main DDD: ISAND(1,:):', ISAND(1,:)
         ZTSURF = 0.  ! Mike Lazard had taken the 'TSURF' out of CLASS. 
                      ! K. Winger put the calculation back in 'energyBudgetDriver'
         call energyBudgetDriver(TBARC, TBARG, TBARCS, TBARGS, THLIQC, THLIQG, & ! Formerly CLASST
                                 THICEC, THICEG, HCPC,   HCPG,   TCTOPC, TCBOTC, TCTOPG, TCBOTG, &
                                 GZEROC, GZEROG, GZROCS, GZROGS, G12C,   G12G,   G12CS,  G12GS, &
                                 G23C,   G23G,   G23CS,  G23GS,  QFREZC, QFREZG, QMELTC, QMELTG, &
                                 EVAPC,  EVAPCG, EVAPG,  EVAPCS, EVPCSG, EVAPGS, TCANO,  TCANS, &
                                 RAICAN, SNOCAN, RAICNS, SNOCNS, CHCAP,  CHCAPS, TPONDC, TPONDG, &
                                 TPNDCS, TPNDGS, TSNOCS, TSNOGS, WSNOCS, WSNOGS, RHOSCS, RHOSGS, &
                                 ITERCT, ZCDH,   ZCDM,   QSENS,  TFLUX,  QEVAP,  EVAPO, &
                                 EVPPOT, ZCFLUX, EVAPB,  ZTSRAD, QS, ZTSURF, &
                                 ST,     SU,     SV,     SQ,     srh, &
                                 GTBS, SFCUBS, SFCVBS, USTARBS, &
                                 ZFSGV,  ZFSGS,  ZFSGG,  ZFLGV,  ZFLGS,  ZFLGG, &
                                 ZHFSC,  ZHFSS,  ZHFSG,  ZHEVC,  ZHEVS,  ZHEVG,  ZHMFC,  ZHMFN, &
                                 ZHTCC,  ZHTCS,  ZHTC,   ZQFCF,  ZQFCL,  CDRAG,  WTABLE, ZILMO, &
                                 ZFRV,   HBL, ZTAIRCAN, ZHUAIRCAN, ZZUSL, ZZTSL, ZUN,    ZTN, &
                                 VPD,    TADP,   RHOAIR, QSWINV, QSWINI, QLWIN,  UA,     VA, &
                                 TA,     QA,     PADRY,  FFC,    FG,     FCS,    FGS,    RBCOEF, &
                                 FSVF,   FSVFS,  PS,     vmod,   ALVSCN, ALIRCN, ALVSG,  ALIRG, &
                                 ALVSCS, ALIRCS, ALVSSN, ALIRSN, ALVSGC, ALIRGC, ALVSSC, ALIRSC, &
                                 TRVSCN, TRIRCN, TRVSCS, TRIRCS, RC,     RCS,    ZWTRG,  groundHeatFlux, qlwavg, &
                                 FRAINC, FSNOWC, FRAICS, FSNOCS, CMASSC, CMASCS, DISP,   DISPS, &
                                 ZOMLNC, ZOELNC, ZOMLNG, ZOELNG, ZOMLCS, ZOELCS, ZOMLNS, ZOELNS, &
                                 TS,     THLIQ,  THICE,  ZTPOND, ZZPOND, ZTBASE, ZTCAN,  ZTSNOW, &
                                 ZSNOW,  ZRHOSNO,ZWSNOW, ZTHPOR, ZTHLRET,ZTHLMIN,ZTHFC,  ZTHLW, &
                                 TRSNOWC,TRSNOWG,ALSNO,  ZFSSB,  FROOT,  FROOTS, &
                                 ZDLAT,  PCPR,   ZHCPS,  ZTCS,   ZTSFS,  DELZ,   ZDELZW, ZZBOTW, &
                                 ZFTEMP, ZFVAP,  ZRIB, &
                                 ISAND, &
                                 ZAILCG,ZAILCGS,ZFCANC,ZFCANCS,ZCO2CONC,ZCO2I1CG, &
                                 ZCO2I1CS,ZCO2I2CG,ZCO2I2CS,COSZS,XDIFFUS,ZSLAI, &
                                 ICC,    ctem_on, ZRMATCTEM,ZFCANCMX, L2MAX, &
                                 NOL2PFTS,  ZCFLUXCG,  ZCFLUXCS, &
                                 ZANCSVEG,ZANCGVEG, ZRMLCSVEG,ZRMLCGVEG, &
                                 TCSNOW, GSNOW, &
                                 ITC,    ITCG,   ITG,    N,      I,  J,  TRNCH,  KOUNT, &
                                 IC,IG,  IZREF,  ISLFD,  NLANDCS,NLANDGS,NLANDC, NLANDG, NLANDI, &
                                 NBS, ISNOALB, DAYL, DAYL_MAX, &
                                 ipeatland, ancsmoss, angsmoss, ancmoss, angmoss, &
                                 rmlcsmoss, rmlgsmoss, rmlcmoss, rmlgmoss, Cmossmas, dmoss, &
                                 iday, pdd)

! CCCCCCCCCCCCCCCContinue here
!print *,'class_main N:',N
!print *,'class_main FCS(1:N)   :',FCS(1:N)   
!print *,'class_main QMELTC(1:N):',QMELTC(1:N)
!print *,'class_main FGS(1:N)   :',FGS(1:N)   
!print *,'class_main QMELTG(1:N):',QMELTG(1:N)
!print *,'class_main TBARGS(1:N,1):',TBARGS(1:N,1)

         call waterBudgetDriver(THLIQ, THICE, TS, ZTCAN, ZRCAN, ZSCAN, & ! Formerly CLASSW
                                ZRUNOFF,ZTRUNOFF,XSNO,  ZTSNOW, ZRHOSNO,ZALBSNO, &
                                ZWSNOW, ZZPOND, ZTPOND, ZGROWTH,ZTBASE, GFLUX, &
                                ZPCFC,  ZPCLC,  ZPCPN,  ZPCPG,  ZQFCF,  ZQFCL, &
                                ZQFN,   ZQFG,   ZQFC,   ZHMFC,  ZHMFG,  ZHMFN, &
                                ZHTCC,  ZHTCS,  ZHTC,   ZROFC,  ZROFN,  ZROVG, &
                                ZWTRS,  ZWTRG,  ZOVRFLW,ZSUBFLW,ZBASFLW, &
                                ZTOVRFL,ZTSUBFL,ZTBASFL,EVAPO,  QFLUX,  RHOAIR, &
                                TBARC,  TBARG,  TBARCS, TBARGS, THLIQC, THLIQG, &
                                THICEC, THICEG, HCPC,   HCPG,   RPCP,   TRPCP, &
                                SPCP,   TSPCP,  PCPR,   TA,     RHOSNI, ZGGEO, &
                                FFC,    FG,     FCS,    FGS,    TPONDC, TPONDG, &
                                TPNDCS, TPNDGS, EVAPC,  EVAPCG, EVAPG,  EVAPCS, &
                                EVPCSG, EVAPGS, QFREZC, QFREZG, QMELTC, QMELTG, &
                                RAICAN, SNOCAN, RAICNS, SNOCNS, FSVF,   FSVFS, &
                                CWLCAP, CWFCAP, CWLCPS, CWFCPS, TCANO, &
                                TCANS,  CHCAP,  CHCAPS, CMASSC, CMASCS, ZSNOW, &
                                GZEROC, GZEROG, GZROCS, GZROGS, G12C,   G12G, &
                                G12CS,  G12GS,  G23C,   G23G,   G23CS,  G23GS, &
                                TSNOCS, TSNOGS, WSNOCS, WSNOGS, RHOSCS, RHOSGS, &
                                ZPLIMC, ZPLIMG, ZPLMCS, ZPLMGS, ZTSFS, &
                                TCTOPC, TCBOTC, TCTOPG, TCBOTG, FROOT, FROOTS, &
                                ZTHPOR, ZTHLRET,ZTHLMIN,ZBI,    ZPSISAT,ZGRKSAT, &
                                ZTHLRAT,ZTHFC,  ZXDRAIN,ZHCPS,  DELZ, &
                                ZDELZW, ZZBOTW, ZXSLOPE,ZGRKFAC,ZWFSURF,ZWFCINT, &
                                ISAND,  igdr, &
                                IWF,    N,      I,      J,      KOUNT, &
                                TRNCH,  IC,     IG,     IG+1,   IG+2, &
                                NLANDCS,NLANDGS,NLANDC, NLANDG, NLANDI)


      ! ToDo: Add call to energyWaterBalanceCheck (1,...)
      !> energyWaterBalanceCheck does the initial calculations for the energy and water balance checks
!      call energyWaterBalanceCheck(1, CTVSTP, CTSSTP, CT1STP, CT2STP, CT3STP, ...)  ! Formerly CLASSZ


!         if (isnoalb == 1) then  ! New routine from Joe (KW)
!            call fourBandDriver(J, ZWSNOW, -1, ZTSNOW, & ! in
!                                ZSNOW, TCSNOW, GSNOW, PCSNGAT, WSNOGAT, & ! in
!                                ALBSGAT, RHOSGAT, & ! in/out
!                                TZSGAT, REFGAT) ! out
!         end if


      endif ! if (.not.dotile.and.doprec)
!
      doprec=dotile
!
    end do CLASS_MAIN_LOOP



! The following calculations are made only on points covered by the current
! mosaic level. This prevents a bug from occuring when unitialized input
! leads to NAN results
!
    FINAL_CALC : DO I=1,N
!
!     ADJUST TERMS DEPENDING ON THE VERTICAL DIFFUSION SOLVING OPTION
!     (NECESSARY MODIFICATION TO USE COUPLED WITH GEM)
!     TFLUX = TFLUX - ZCFLUX * "TPOTA" = TFLUX - ZCFLUX * [Ta + Z*g/Cp]
!     QFLUX = QFLUX - ZCFLUX * QA
!
          IF (IMPFLX) THEN
             TFLUX(I) = TFLUX(I) - ZCFLUX(I)*(TA(I)             &
                      + GRAV/CPD * ( FCS(I)*(ZZUSL(I)-DISPS(I)) &
                      +              FFC(I)*(ZZUSL(I)-DISP (I)) &
                      + ( FGS(I) + FG(I) )*ZZUSL(I)))
             QFLUX(I) = QFLUX(I) - ZCFLUX(I)*QA(I)
!     OFFLINE PREVIOUS OPTION
          ELSE
             ZCFLUX(I) = 0.
          ENDIF
!
!   BM CALCULATIONS
!     IN  ISBA,  BM = VMOD * ( KARMAN/ FM )^2
!     IN CLASS, CDM = ( KARMAN/ FM )^2
!   ADD CALCULATION FOR BM OUTSIDE CLASS TO FACILITATE CODE MAINTENANCE
!   AND VERSION CHANGES
!
          VMOD(I) = SQRT( MAX( VAMIN,UA(I)*UA(I)+VA(I)*VA(I) ) )
          ZBM(I)  = VMOD(I) * ZCDM(I)
!
!   ADD CALCULATION FOR Z0M & ZOT OUTSIDE CLASS TO FACILITATE CODE MAINTENANCE
!   AND VERSION CHANGES
!   CALCULTATION OF ROUGHNESS LENGHTS AVERAGED OVER THE 4 SURFACE TYPES :
!   CANOPY OVER SNOW, SNOW OVER GROUND, CANOPY OVER GROUND & BARE GROUND
!
          Z0M(I) = FCS(I)*EXP(ZOMLCS(I)) + FGS(I)*EXP(ZOMLNS(I)) &
                 + FFC(I)*EXP(ZOMLNC(I)) +  FG(I)*EXP(ZOMLNG(I))
          Z0H(I) = FCS(I)*EXP(ZOELCS(I)) + FGS(I)*EXP(ZOELNS(I)) &
                 + FFC(I)*EXP(ZOELNC(I)) +  FG(I)*EXP(ZOELNG(I))
!
!
          ALVIS_SOL(I) = 0.5*(ZALVS(I)+ZALIR(I))
          FSOLUACC(I)  = FSOLUACC(I)+FLUSOL(I)*ALVIS_SOL(I)*DT
          FIRUACC(I)   = FIRUACC(I) +QLWAVG(I) * DT
          ZFL(I)       = FFC(I)*GZEROC(I)+FG(I)*GZEROG(I)+FCS(I)*GZROCS(I)+ &
                         FGS(I)*GZROGS(I)
!
    end do FINAL_CALC
!
!  end do DO_MOSAIC
!



  zTT2m(:)    = st(:)                     ! instantaneos

if (prints_L) then
print*,'class_main'
print*,'class_main Physics output'
print*,'class_main ALVIS_SOL      :',minval(ALVIS_SOL),maxval(ALVIS_SOL),sum(ALVIS_SOL)/(N)
print*,'class_main QFLUX          :',minval(QFLUX),maxval(QFLUX),sum(QFLUX)/(N)
print*,'class_main TFLUX          :',minval(TFLUX),maxval(TFLUX),sum(TFLUX)/(N)
print*,'class_main ZCFLUX         :',minval(ZCFLUX),maxval(ZCFLUX),sum(ZCFLUX)/(N)
print*,'class_main ZBM            :',minval(ZBM),maxval(ZBM),sum(ZBM)/(N)
print*,'class_main QSENS          :',minval(QSENS),maxval(QSENS),sum(QSENS)/(N)
print*,'class_main QEVAP          :',minval(QEVAP),maxval(QEVAP),sum(QEVAP)/(N)
print*,'class_main ZFRV           :',minval(ZFRV),maxval(ZFRV),sum(ZFRV)/(N)
print*,'class_main HBL            :',minval(HBL),maxval(HBL),sum(HBL)/(N)
print*,'class_main ZILMO          :',minval(ZILMO),maxval(ZILMO),sum(ZILMO)/(N)
print*,'class_main QS             :',minval(QS),maxval(QS),sum(QS)/(N)
print*,'class_main ZTSURF         :',minval(ZTSURF),maxval(ZTSURF),sum(ZTSURF)/(N)
print*,'class_main SQ             :',minval(SQ),maxval(SQ),sum(SQ)/(N)
print*,'class_main ST             :',minval(ST),maxval(ST),sum(ST)/(N)
print*,'class_main SU             :',minval(SU),maxval(SU),sum(SU)/(N)
print*,'class_main SV             :',minval(SV),maxval(SV),sum(SV)/(N)
print*,'class_main ZTSRAD         :',minval(ZTSRAD),maxval(ZTSRAD),sum(ZTSRAD)/(N)
print*,'class_main ZFL            :',minval(ZFL),maxval(ZFL),sum(ZFL)/(N)
print*,'class_main ZFTEMP         :',minval(ZFTEMP),maxval(ZFTEMP),sum(ZFTEMP)/(N)
print*,'class_main ZFVAP          :',minval(ZFVAP),maxval(ZFVAP),sum(ZFVAP)/(N)
endif


! FILL THE ARRAYS TO BE AGGREGATED LATER IN S/R AGREGE
  CALL FILLAGG ( BUS, BUSSIZ, PTSURF, PTSURFSIZ, INDX_SOIL, SURFLEN )
!
  DEALLOCATE( PAIDAT,HGTDAT,ACVDAT,ACIDAT, &
              TBARC ,TBARG ,TBARCS,TBARGS, &
              THLIQC,THLIQG,THICEC,THICEG, &
              HCPC  ,HCPG  ,FROOT ,FROOTS, GFLUX , &
              TCTOPC,TCBOTC,TCTOPG,TCBOTG, &
              ISAND, IORG,  ITERCT )
!
  RETURN
END
