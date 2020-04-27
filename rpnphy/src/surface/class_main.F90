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
  use class_configs
  use mu_jdate_mod, only: jdate_day_of_year, mu_js2ymdhms
  use sfc_options
  use sfcbus_mod
  use phy_options, only: MU_JDATE_HALFDAY
  use ctem_params, only : initpftpars, nlat, ilg, ican, ignd

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
!016       K. Winger  (Mar 2010) -Add RUNOFFTOT, aggregated surface runoff
!                                 and DRAINTOT, aggregated base drainage
!017       D. Deacu   (Mar 2010) -Set Z0ORO = f(Z0)
!018       K. Winger  (Nov 2001) -Change GGEO to ZGGEO and set it to value from namelist
!019       L. Duarte  (May 2010) -Add code to run CTEM in coupled mode
!020       K. Winger  (Nov 2012) -Pass FSNOW to output, renamed to ZFSNOW
!021       L. Duarte  (    2014) -Update for use with CLASS 3.6/CTEM 2.0
!022       K. Winger  (Aug 2016) -Initialize if soil just appeard due to melting glacier
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
  integer :: IC,ICP1,IPAI,IHGT,IALC,IALS,IALG,IPCP, ICC
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
           TRSNOW(N),   QLWAVG(N), &
           ALVS  (N),   ALIR  (N),   QSOL  (N), &
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

! Thickness and depth of each soil layer in meters
      real :: DELZ(IG),ZBOT(IG)
!     DATA  DELZ    /0.10,0.25,3.75/
!     DATA  ZBOT    /0.10,0.35,4.10/
!
!******************************************************
!
!
!  Allocatable work arrays (depending on IG=CLASS_IG and IC=CLASS_IC)
!     real :: PAIDAT(N,IC),HGTDAT(N,IC),ACVDAT(N,IC),ACIDAT(N,IC)
!     real :: TBARC (N,IG),TBARG (N,IG),TBARCS(N,IG),TBARGS(N,IG),
!    1     THLIQC(N,IG),THLIQG(N,IG),THICEC(N,IG),THICEG(N,IG),
!    2     HCPC  (N,IG),HCPG  (N,IG),FROOT (N,IG),GFLUX (N,IG),
!    3     TCTOPC(N,IG),TCBOTC(N,IG),TCTOPG(N,IG),TCBOTG(N,IG)
!     integer ::  ISAND(N,IG),IORG(N,IG)
!     integer ::  ITERCT(N,6,50)
!
  integer, DIMENSION(:,:,:), ALLOCATABLE :: ITERCT
  integer, DIMENSION(:,:),   ALLOCATABLE :: ISAND,IORG
  REAL, DIMENSION(:,:), ALLOCATABLE :: PAIDAT,HGTDAT,ACVDAT,ACIDAT
  REAL, DIMENSION(:,:), ALLOCATABLE :: TBARC ,TBARG ,TBARCS,TBARGS, &
                                       THLIQC,THLIQG,THICEC,THICEG, &
                                       HCPC  ,HCPG  ,FROOT ,GFLUX, &
                                       TCTOPC,TCBOTC,TCTOPG,TCBOTG
!
  real :: XDIFFUS(N)
  real :: zlyglfmas(N,CTEM_ICC)

!
  real,pointer,dimension(:)   :: FSOLUACC, FIRUACC, SU, SV, ST, SQ, ALVIS_SOL
  real,pointer,dimension(:)   :: EVAPO, QSENS, QEVAP, HBL, ZILMO, ZFRV, PS, QS
  real,pointer,dimension(:)   :: Z0H, Z0ORO, ZTSURF, ZTSRAD, UA, VA, TA, QA
  real,pointer,dimension(:)   :: TFLUX, QFLUX, ZCANG, FLUSOL, ZHUAIRCAN
  real,pointer,dimension(:)   :: ZDLAT, ZDLON, ZSDEPTH, ZZTSL, ZZUSL, QLWIN
  real,pointer,dimension(:)   :: ZTAIRCAN, ZTSNOW, ZTBASE, ZTPOND, ZZPOND
  real,pointer,dimension(:)   :: ZRHOSNO, ZRUNOFF, ZRUNOFFTOT, ZDRAINTOT
  real,pointer,dimension(:)   :: ZSCAN, ZINISOIL, XSNO, ZALBSNO, ZGROWTH
  real,pointer,dimension(:)   :: ZXDRAIN, ZXSLOPE, ZGRKFAC, ZWFSURF, ZWFCINT
  real,pointer,dimension(:)   :: ZCMAI, ZFSGV, ZFSGS, ZFSGG, ZFSNOW, ZFLGV
  real,pointer,dimension(:)   :: ZFLGS, ZFLGG, ZHFSC, ZHFSS, ZHFSG, ZHEVC
  real,pointer,dimension(:)   :: ZHEVS, ZHEVG, ZHMFC, ZHTCC, ZHTCS, ZPCFC, ZPCLC
  real,pointer,dimension(:)   :: ZPCPG, ZQFCF, ZQFCL, ZQFG, ZQFN
  real,pointer,dimension(:)   :: ZTBASFL, ZTOVRFL, ZTRUNOFF, ZTSUBFL, ZWSNOW
  real,pointer,dimension(:)   :: ZWTRC, ZWTRS, ZWTRG, ZROFC, ZROFN, ZROVG
  real,pointer,dimension(:)   :: ZOVRFLW, ZSUBFLW, ZBASFLW, ZQSWD, ZTCAN, ZRCAN
  real,pointer,dimension(:)   :: ZCFLUX, ZPCPN, ZHMFN, ZALGWET, ZALGDRY
  real,pointer,dimension(:)   :: FFC, FCS, FG, FGS, ZFL, ZRAINRATE, ZSNOWRATE
  real,pointer,dimension(:)   :: ZSNOW, TINDX, MOSFRAC

  real,pointer,dimension(:)   :: ZFTEMP, ZFVAP, ZRIB, ZCDH, ZCDM, ZBM

  real,pointer,dimension(:)   :: ZSLP, ZARE, ZLEG, ZANIS, ZEXCW, ZLBEDR, ZWT, ZWTNEW

  real,pointer,dimension(:)   :: ztt2m

  real,pointer,dimension(:,:) :: THLIQ, THICE, TS, ZFCANMX, ZSAND, ZCLAY
  real,pointer,dimension(:,:) :: ZORGM, ZDELZW, ZZBOTW, ZTHPOR, ZTHLMIN
  real,pointer,dimension(:,:) :: ZTHLRET, ZPSISAT, ZBI, ZPSIWLT, ZHCPS
  real,pointer,dimension(:,:) :: ZTCS, ZTSFS, ZTHFC
  real,pointer,dimension(:,:) :: ZGRKSAT, ZGRKTLD, ZTHLRAT, ZHTC
  real,pointer,dimension(:,:) :: ZZOLN, ZALVSC, ZALIRC, ZPAIMAX, ZPAIMIN
  real,pointer,dimension(:,:) :: ZCWGTMX, ZZRTMAX, ZRSMIN, ZQA50, ZVPDA, ZVPDB
  real,pointer,dimension(:,:) :: ZPSIGA, ZPSIGB, ZQFC, ZHMFG

  real,pointer,dimension(:,:) :: ZFCANCMX
  real,pointer,dimension(:,:) :: ZAILCG, ZAILCB, ZAILC, ZZOLNC, ZSLAI, ZBMASVEG
  real,pointer,dimension(:,:) :: ZCMASVEGC, ZVEGHGHT, ZROOTDPTH, ZALVSCTM, ZALIRCTM 
  real,pointer,dimension(:,:) :: ZPAIC, ZSLAIC
  real,pointer,dimension(:,:,:) :: ZRMATC, ZRMATCTEM

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
  integer :: ictemmod
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
  logical :: pftexist(n,ctem_icc)

! Fire-related variables
! Not used for now (LD)
  real    :: emit_co2gat(n,ctem_icc), &
             emit_cogat(n,ctem_icc), &
             emit_ch4gat(n,ctem_icc), &
             emit_nmhcgat(n,ctem_icc), &
             emit_h2gat(n,ctem_icc), &
             emit_noxgat(n,ctem_icc), &
             emit_n2ogat(n,ctem_icc), &
             emit_pm25gat(n,ctem_icc), &
             emit_tpmgat(n,ctem_icc), &
             emit_tcgat(n,ctem_icc), &
             emit_ocgat(n,ctem_icc), &
             emit_bcgat(n,ctem_icc), &
!             burnvegfgat(n,ctem_icc), &
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

  ICC  = CTEM_ICC
  IC   = CLASS_IC
  ICP1 = CLASS_IC+1
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
  if (ctem_on) then
    ictemmod=1
    call initpftpars(compete)
  else
    ictemmod=0
  endif

  mosaic=.false.  ! mosaic option for competition purposes
                       ! always set it to false, even when mosaic is used! (LD)
  dofire=.false.
  LNDUSEON=.false.

!
  ALLOCATE( PAIDAT(N,IC),HGTDAT(N,IC),ACVDAT(N,IC),ACIDAT(N,IC), &
            TBARC (N,IG),TBARG (N,IG),TBARCS(N,IG),TBARGS(N,IG), &
            THLIQC(N,IG),THLIQG(N,IG),THICEC(N,IG),THICEG(N,IG), &
            HCPC  (N,IG),HCPG  (N,IG),FROOT (N,IG),GFLUX (N,IG), &
            TCTOPC(N,IG),TCBOTC(N,IG),TCTOPG(N,IG),TCBOTG(N,IG), &
            ISAND (N,IG),IORG  (N,IG),ITERCT(N,6,50) )
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
  PS        (1:N)      => bus( x(PMOINS ,1,1) : )  !dynamic
  ZCLAY     (1:N,1:IG) => bus( x(CLAY   ,1,1) : )  !permanent
  ZORGM     (1:N,1:IG) => bus( x(ORGM   ,1,1) : )  !permanent
  ZCANG     (1:N)      => bus( x(CANG   ,1,1) : )  !volatile
  ZQSWD     (1:N)      => bus( x(QSWD   ,1,1) : )  !volatile, init in CLASS
  ZSAND     (1:N,1:IG) => bus( x(SAND   ,1,1) : )  !permanent
  TA        (1:N)      => bus( x(TMOINS ,1,nk) : )  !dynamic
  QA        (1:N)      => bus( x(HUMOINS,1,nk) : )  !dynamic, MOSAIC
  ZDLAT     (1:N)      => bus( x(DLAT   ,1,1) : )  !permanent
  ZDLON     (1:N)      => bus( x(DLON   ,1,1) : )  !permanent
  UA        (1:N)      => bus( x(UMOINS ,1,nk) : )  !dynamic
  VA        (1:N)      => bus( x(VMOINS ,1,nk) : )  !dynamic
  ZZTSL     (1:N)      => bus( x(ZTSL   ,1,1) : )  !volatile, calc in GEM
  ZZUSL     (1:N)      => bus( x(ZUSL   ,1,1) : )  !volatile, calc in GEM
  IF (RADSLOPE) THEN
     FLUSOL (1:N)      => bus( x(FLUSLOP,1,1) : )  !permanent
  ELSE
     FLUSOL (1:N)      => bus( x(FLUSOLIS,1,1) : )  !permanent
  ENDIF
  QLWIN     (1:N)      => bus( x(FDSI   ,1,1) : )  !permanent
  ZXDRAIN   (1:N)      => bus( x(XDRAIN ,1,1) : )  !permanent
  ZGRKFAC   (1:N)      => bus( x(GRKFAC ,1,1) : )
  ZWFSURF   (1:N)      => bus( x(WFSURF ,1,1) : )
  ZWFCINT   (1:N)      => bus( x(WFCINT ,1,1) : )
  ZRAINRATE (1:N)      => bus( x(RAINRATE,1,1) : )  !volatile, calc in GEM
  ZSNOWRATE (1:N)      => bus( x(SNOWRATE,1,1) : )  !volatile, calc in GEM
  ZSNOW     (1:N)      => bus( x(SNODP  ,1,indx_sfc ) : )  !permanent
!!!  ZRUNOFFTOT(1:N) => bus( x(RUNOFFTOT,1,indx_sfc) : )
!!!  ZDRAINTOT (1:N) => bus( x(DRAINTOT,1,indx_sfc) : )
!!!  ZINISOIL  (1:N) => bus( x(INISOIL_SL,1,1) : )



  igdr(i) = 1
  DO J=1,IG
     DO I=1,N
        ISAND(I,J) = NINT(ZSAND(I,J))
        IF (ISAND(I,J).GT.-3) IGDR(I) = J  ! Index of soil layer in which bedrock is encountered
        bterm(i,j) = 0.159*zclay(i,j)+2.91
     ENDDO
  ENDDO
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
     ZQSWD(I)  = 2*QSOL(I)*XDIFFUS(I)
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
    ZWFSURF   (1:N) => bus( x(WFSURF  ,1,1) : )

    Z0H       (1:N) => bus( x( Z0T    ,1,indx_sfc ) : )  !permanent
    Z0ORO     (1:N) => bus( x( Z0     ,1,indx_sfc ) : )  !permanent
    QFLUX     (1:N) => bus( x( ALFAQ  ,1,1 ) : )  !volatile, CLASST
    TFLUX     (1:N) => bus( x( ALFAT  ,1,1 ) : )  !volatile, CLASST
    ZCFLUX    (1:N) => bus( x( BT     ,1,1 ) : )  !volatile, TPREP
    ZBM       (1:N) => bus( x( BM     ,1,1 ) : )  !volatile, TPREP
    EVAPO     (1:N) => bus( x( WFLUX  ,1,1 ) : )  !volatile, CLASST
    FFC       (1:N) => bus( x( FCOVC  ,1,1 ) : )  !volatile, APREP
    QSENS     (1:N) => bus( x( FC     ,1,indx_sfc ) : )  !volatile, TPREP
    FCS       (1:N) => bus( x( FCOVCS ,1,1 ) : )  !volatile, TPREP
    FG        (1:N) => bus( x( FCOVG  ,1,1 ) : )  !volatile, TPREP
    FGS       (1:N) => bus( x( FCOVGS ,1,1 ) : )  !volatile, TPREP
    ZFRV      (1:N) => bus( x( FRV    ,1,indx_sfc ) : )  !permanent
    QEVAP     (1:N) => bus( x( FV     ,1,indx_sfc ) : )  !volatile, TPREP
    HBL       (1:N) => bus( x( HST    ,1,indx_sfc ) : )  !permanent
    ZILMO     (1:N) => bus( x( ILMO   ,1,indx_sfc ) : )  !permanent
    THICE     (1:N,1:IG) => bus( x(ISOIL,1,1 ) : )  !permanent
    MOSFRAC   (1:N) => bus( x( MOSFRACT ,1,1 ) : ) !permanent
    SQ        (1:N) => bus( x( QDIAG  ,1,1 ) : ) !permanent
    QS        (1:N) => bus( x( QSURF  ,1,indx_sfc ) : ) !permanent
    ST        (1:N) => bus( x( TDIAG  ,1,1 ) : ) !permanent
    TINDX     (1:N) => bus( x( TINDEX ,1,1 ) : ) !permanent
    TS        (1:N,1:IG) => bus( x( TSOIL  ,1,1 ) : ) !permanent
    SU        (1:N) => bus( x( UDIAG  ,1,1 ) : ) !permanent
    SV        (1:N) => bus( x( VDIAG  ,1,1 ) : ) !permanent
    THLIQ     (1:N,1:IG) => bus( x( WSOIL  ,1,1 ) : ) !permanent
    XSNO      (1:N) => bus( x( SNOMA  ,1,1 ) : ) !permanent
    ZALBSNO   (1:N) => bus( x( SNOAL  ,1,1 ) : ) !permanent
    ZALGDRY   (1:N) => bus( x( ALGDRY ,1,1 ) : ) !permanent
    ZALGWET   (1:N) => bus( x( ALGWET ,1,1 ) : ) !permanent
    ZBASFLW   (1:N) => bus( x( DRAIN  ,1,1 ) : ) !permanent
    ZBI       (1:N,1:IG) => bus( x( BBI    ,1,1 ) : ) !permanent
    ZCMAI     (1:N) => bus( x( CMAI   ,1,1 ) : )
    ZDELZW    (1:N,1:IG) => bus( x( DELZW  ,1,1 ) : ) !permanent
    ZFLGG     (1:N) => bus( x( FLGG   ,1,1 ) : ) !volatile, TPREP
    ZFLGS     (1:N) => bus( x( FLGS   ,1,1 ) : ) !volatile, TPREP
    ZFLGV     (1:N) => bus( x( FLGV   ,1,1 ) : ) !volatile, TPREP
    ZFSGG     (1:N) => bus( x( FSGG   ,1,1 ) : ) !volatile, TPREP
    ZFSGS     (1:N) => bus( x( FSGS   ,1,1 ) : ) !volatile, TPREP
    ZFSGV     (1:N) => bus( x( FSGV   ,1,1 ) : ) !volatile, TPREP
    ZFSNOW    (1:N) => bus( x( FSNOW  ,1,1 ) : ) !volatile
    ZGRKSAT   (1:N,1:IG) => bus( x( GRKSAT ,1,1 ) : ) !permanent
    ZGRKTLD   (1:N,1:IG) => bus( x( GRKTLD ,1,1 ) : ) !permanent
    ZGROWTH   (1:N) => bus( x( VEGGRO ,1,1 ) : ) !permanent
    ZHCPS     (1:N,1:IG) => bus( x( HCPS   ,1,1 ) : ) !permanent
    ZHEVC     (1:N) => bus( x( HEVC   ,1,1 ) : ) !volatile, TPREP
    ZHEVG     (1:N) => bus( x( HEVG   ,1,1 ) : ) !volatile, TPREP
    ZHEVS     (1:N) => bus( x( HEVS   ,1,1 ) : ) !volatile, TPREP
    ZHFSC     (1:N) => bus( x( HFSC   ,1,1 ) : ) !volatile, TPREP
    ZHFSG     (1:N) => bus( x( HFSG   ,1,1 ) : ) !volatile, TPREP
    ZHFSS     (1:N) => bus( x( HFSS   ,1,1 ) : ) !volatile, TPREP
    ZHMFC     (1:N) => bus( x( HMFC   ,1,1 ) : ) !volatile, TPREP
    ZHMFG     (1:N,1:IG) => bus( x( HMFG   ,1,1 ) : ) !volatile, WPREP
    ZHMFN     (1:N) => bus( x( HMFN   ,1,1 ) : ) !volatile, TPREP
    ZHTC      (1:N,1:IG) => bus( x( HTC    ,1,1 ) : ) !volatile, APREP
    ZHTCC     (1:N) => bus( x( HTCC   ,1,1 ) : ) !volatile, APREP
    ZHTCS     (1:N) => bus( x( HTCS   ,1,1 ) : ) !volatile, APREP
    ZOVRFLW   (1:N) => bus( x( OVERFL ,1,1 ) : ) !volatile, WPREP
    ZPCFC     (1:N) => bus( x( PCFC   ,1,1 ) : ) !volatile, WPREP
    ZPCLC     (1:N) => bus( x( PCLC   ,1,1 ) : ) !volatile, WPREP
    ZPCPG     (1:N) => bus( x( PCPG   ,1,1 ) : ) !volatile, WPREP
    ZPCPN     (1:N) => bus( x( PCFG   ,1,1 ) : ) !volatile, WPREP
    ZPSISAT   (1:N,1:IG) => bus( x( PSISAT ,1,1 ) : ) !permanent
    ZPSIWLT   (1:N,1:IG) => bus( x( PSIWLT ,1,1 ) : ) !permanent
    ZQFC      (1:N,1:IG) => bus( x( QFC    ,1,1 ) : ) !*volatile, WPREP
    ZQFCF     (1:N) => bus( x( QFCF   ,1,1 ) : ) !*volatile, WPREP
    ZQFCL     (1:N) => bus( x( QFCL   ,1,1 ) : ) !*volatile, WPREP
    ZQFG      (1:N) => bus( x( QFG    ,1,1 ) : ) !*volatile, WPREP
    ZQFN      (1:N) => bus( x( QFN    ,1,1 ) : ) !*volatile, WPREP
    ZRCAN     (1:N) => bus( x( WVEG   ,1,1 ) : ) !permanent
    ZRHOSNO   (1:N) => bus( x( SNODEN ,1,1 ) : ) !permanent
    ZROFC     (1:N) => bus( x( ROFC   ,1,1 ) : ) !*volatile, WPREP
    ZROFN     (1:N) => bus( x( ROFN   ,1,1 ) : ) !*volatile, WPREP
    ZROVG     (1:N) => bus( x( ROVG   ,1,1 ) : ) !*volatile, WPREP
    ZRUNOFF   (1:N) => bus( x( RUNOFF ,1,1 ) : ) !volatile, CLASSW
    ZSCAN     (1:N) => bus( x( IVEG   ,1,1 ) : ) !permanent
    ZSUBFLW   (1:N) => bus( x( SUBFLW ,1,1 ) : ) !volatile, WPREP
    ZTBASE    (1:N) => bus( x( TBASE  ,1,1 ) : ) !permanent
    ZTCAN     (1:N) => bus( x( TVEG   ,1,1 ) : ) !permanent
    ZTCS      (1:N,1:IG) => bus( x( TCS    ,1,1 ) : ) !permanent
    ZTHFC     (1:N,1:IG) => bus( x( THFC   ,1,1 ) : ) !permanent
    ZTHLMIN   (1:N,1:IG) => bus( x( THLMIN ,1,1 ) : ) !permanent
    ZTHLRAT   (1:N,1:IG) => bus( x( THLRAT ,1,1 ) : ) !permanent
    ZTHLRET   (1:N,1:IG) => bus( x( THLRET ,1,1 ) : ) !permanent
    ZTHPOR    (1:N,1:IG) => bus( x( THPOR  ,1,1 ) : ) !permanent
    ZTPOND    (1:N) => bus( x( TPOND  ,1,1 ) : ) !permanent
    ZTSNOW    (1:N) => bus( x( TSNO   ,1,1 ) : ) !permanent
    ZTSRAD    (1:N) => bus( x( TSRAD  ,1,1 ) : ) !permanent
    ZTSURF    (1:N) => bus( x( TSURF  ,1,1 ) : ) !permanent
    ZWSNOW    (1:N) => bus( x( WSNOW  ,1,1 ) : ) !permanent
    ZWTRC     (1:N) => bus( x( WTRC   ,1,1 ) : ) !*volatile, APREP
    ZWTRG     (1:N) => bus( x( WTRG   ,1,1 ) : ) !*volatile, APREP
    ZWTRS     (1:N) => bus( x( WTRS   ,1,1 ) : ) !*volatile, APREP
    ZXSLOPE   (1:N) => bus( x( XSLOPE ,1,1 ) : ) !permanent
    ZZBOTW    (1:N,1:IG) => bus( x( ZBOTW  ,1,1 ) : ) !permanent
    ZZPOND    (1:N) => bus( x( ZPOND  ,1,1 ) : ) !permanent
    FSOLUACC  (1:N) => bus( x( FSOLUPAF,1,1 ) : ) !permanent
    FIRUACC   (1:N) => bus( x( FIRUPAF,1,1 ) : ) !permanent
    ALVIS_SOL (1:N) => bus( x( ALVIS  ,1,1 ) : ) !permanent
    ZFL       (1:N) => bus( x( FL     ,1,1 ) : ) !volatile, CLASS
!
    ZSDEPTH   (1:N) => bus( x( SDEPTH ,1,1 ) : ) !permanent
    ZFCANMX   (1:N,1:ICP1) => bus( x( FCANMX ,1,1 ) : ) !permanent
    ZZOLN     (1:N,1:ICP1) => bus( x( ZOLN   ,1,1 ) : ) !permanent
    ZALVSC    (1:N,1:ICP1) => bus( x( ALVSC  ,1,1 ) : ) !permanent
    ZALIRC    (1:N,1:ICP1) => bus( x( ALIRC  ,1,1 ) : ) !permanent
    ZPAIMAX   (1:N,1:IC) => bus( x( LAIMAX ,1,1 ) : ) !permanent
    ZPAIMIN   (1:N,1:IC) => bus( x( LAIMIN ,1,1 ) : ) !permanent
    ZCWGTMX   (1:N,1:IC) => bus( x( VEGMA  ,1,1 ) : ) !permanent
    ZZRTMAX   (1:N,1:IC) => bus( x( ROOTDP ,1,1 ) : ) !permanent
    ZRSMIN    (1:N,1:IC) => bus( x( STOMR  ,1,1 ) : ) !permanent
    ZQA50     (1:N,1:IC) => bus( x( QA50   ,1,1 ) : ) !permanent
    ZVPDA     (1:N,1:IC) => bus( x( VPDA   ,1,1 ) : ) !permanent
    ZVPDB     (1:N,1:IC) => bus( x( VPDB   ,1,1 ) : ) !permanent
    ZPSIGA    (1:N,1:IC) => bus( x( PSIGA  ,1,1 ) : ) !permanent
    ZPSIGB    (1:N,1:IC) => bus( x( PSIGB  ,1,1 ) : ) !permanent
    ZHUAIRCAN (1:N) => bus( x( HUAIRCAN,1,1 ) : ) !permanent
    ZTAIRCAN  (1:N) => bus( x( TAIRCAN,1,1 ) : ) !permanent
    ZTBASFL   (1:N) => bus( x( TBASFL ,1,1 ) : ) !volatile, WPREP
    ZTOVRFL   (1:N) => bus( x( TOVRFL ,1,1 ) : ) !volatile, WPREP
    ZTSUBFL   (1:N) => bus( x( TSUBFL ,1,1 ) : ) !volatile, WPREP
    ZTRUNOFF  (1:N) => bus( x( TRUNOFF,1,1 ) : ) !volatile, WPREP
    ZTSFS     (1:N,1:4) => bus( x( TSURFSA,1,1 ) : ) !permanent
!
    ZFTEMP    (1:N) => bus( x( FTEMP  ,1,indx_sfc ) : ) !permanent
    ZFVAP     (1:N) => bus( x( FVAP   ,1,indx_sfc ) : ) !permanent
    ZRIB      (1:N) => bus( x( RIB    ,1,1 ) : ) !permanent
    ZCDH      (1:N) => bus( x( CDH    ,1,1 ) : ) !*volatile, TPREP
    ZCDM      (1:N) => bus( x( CDM    ,1,1 ) : ) !*volatile, TPREP

    ZSLP      (1:N) => bus( x( SLPGW    ,1,1 ) : )
    ZARE      (1:N) => bus( x( ARE      ,1,1 ) : )
    ZLEG      (1:N) => bus( x( LEGGW    ,1,1 ) : )
    ZANIS     (1:N) => bus( x( ANIS     ,1,1 ) : )
    ZEXCW     (1:N) => bus( x( EXCW     ,1,1 ) : )
    ZLBEDR    (1:N) => bus( x( LBEDR    ,1,1 ) : )
    ZWT       (1:N) => bus( x( TOTW     ,1,1 ) : )
    ZWTNEW    (1:N) => bus( x( WTNEW    ,1,1 ) : )
!
!
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
!*************************************************************
!
!
!    KOUNT_EQ_0 : IF(KOUNT.EQ.0 .or. maxval(zinisoil).gt.0.5) THEN
    KOUNT_EQ_0 : IF(kount0) THEN
!
        if (igwscheme.ne.1) then
           do i=1,n
!              ZEXCW(I)=0.0
!              ZWTNEW(I)=0.0
              do J=1,ig
!                 thpf(I,J)=0.0
                 RW(I,J)=0.0
              enddo
 
              zsdepth(i)=zsdepth(i)-MOD(nint(zsdepth(i)*100.0),10)/100.0
!              WT(i)=WT2
 
              jji=0
              do 3212 j=1,IG
                 IF (zsdepth(i).LT.(zbot(j)+0.01).AND.jji.EQ.0) THEN
                    jji=1
                    Ibed(i)=j
                 ELSEIF (jji.EQ.0) THEN 
                    Ibed(i)=j
                 ENDIF
3212          continue
           enddo
        endif

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
        call soilProperties(ZTHPOR, ZTHLRET, ZTHLMIN, ZBI, ZPSISAT, & ! Formerly CLASSB
                            ZGRKSAT, ZTHLRAT, ZHCPS, ZTCS, &
                            ZTHFC, ZPSIWLT, ZDELZW, ZZBOTW, ZALGWET, &
                            ZALGDRY, ZSAND, ZCLAY, ZORGM, DELZ, ZBOT, ZSDEPTH, &
                            ISAND, igdr, N, NMIM, 1, N, NMIM, IG, ictemmod, &
                            ZWT, maxd, igwscheme)
!


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



    igwscheme_eq_1 : if (igwscheme.eq.1) then
       do I=1,N
          Isat(i)=ig
          Isats(i)=ig
       enddo
    else
! *********************
       do I=1,N
          IF (kount0) THEN
             ZWTNEW(I)=ZWT(I)
          ELSE
!             IF (ABS(WT(I,N)-ZWTNEW(I)).GT..5) then
!                print *,(WT(I,N)-ZWTNEW(I))
!             endif
             ZWT(I)=ZWTNEW(I)
             Ibed(i)=nint(ZLBEDR(i))
          ENDIF
! *********************
          IF (ZWTNEW(I)/(.2).LE.25.0) THEN
             XWT(I)=ZSDEPTH(I)+25.-ZWTNEW(I)/(.2)
          ELSEIF (ZWTNEW(I)/(.2).GT.25.0) THEN
            do 322 J=Ibed(I),1,-1
               dumw=0.0
               sdlz=0.0
               do 323 K=J+1,Ibed(I)
                    thpf(I,K)= MAX(ZTHPOR(I,K)-THICE(I,K)-0.00001,0.0)
                    dumw=dumw+(DELZ(K)*thpf(I,K))
                    sdlz=sdlz+DELZ(K)
323            ENDDO
               IF(ZTHPOR(I,J)-THICE(I,J).GT.0.0)  THEN                   
                    RW(I,J)=(ZWTNEW(I)-25.0*.2-dumw) / &
                    (ZTHPOR(I,J)-THICE(I,J)-0.00001)  ! this finds water remained for the next (up) layer
               ELSE
                    RW(I,J)=RW(I,J+1)   ! this mean the same water is still available when layer frozen compeletly
                    if (j.eq.ig) then
                      write(6,*),'class.ftn: j=ig at ',i,trnch,' ; RW = ',RW(I,J+1)
                      RW(I,J)=0.0
                    endif
               ENDIF
               IF(RW(I,J).LT.0.0) EXIT 
322         ENDDO
            XWT(I)=ZSDEPTH(I)-RW(I,J+1)-sdlz+DELZ(J+1)   
          ENDIF
!****************
          IF (XWT(I).LT.0.0) XWT(I)=0.0
          XWT(I) = MAX(maxd+0.0001,XWT(I))
          jji=0
          do 321 j=1,IG
            IF(ZSDEPTH(I).LE.0.1) THEN
                jji=1
                Isat(I)=1
            ELSEIF (ZSDEPTH(I).LT.(zbot(j)+0.001).AND.jji.EQ.0) THEN
                jji=1
                Isat(I)=j-1
                IF (ZTHPOR(I,j).EQ.0.0) THEN 
                    Isat(I)=j-1
                    IF (Ibed(I).EQ.Isat(I)+1) Ibed(I)=Ibed(I)-1
                ENDIF
            ELSEIF (XWT(I).LE.zbot(j).AND.jji.EQ.0) THEN
                jji=1
                Isat(I)=j-1
            ELSEIF (jji.EQ.0) THEN
                Isat(I)=j
            ENDIF
321       continue
!   *******************
            Isats(I) = Isat(I)   
            Isat(I)  = Ibed(I)
!***********************
       enddo
    endif igwscheme_eq_1
!
!
    doprec = .false.
!
    CLASS_MAIN_LOOP : do k=1,n+1
!
!      dotile = (k.ne.n+1.and.((nmos.EQ.0).OR.(MOSFRAC(k).gt.0)))
      dotile = (k.ne.n+1)
      if (dotile.and.doprec) then
         j = k
      elseif (dotile.and..not.doprec) then
         j = k
         i = k
      elseif (.not.dotile.and.doprec) then
         CALL atmosphericVarsCalc(VPD, TADP, PADRY, RHOAIR, RHOSNI, & ! Formerly CLASSI
                                  RPCP, TRPCP, SPCP, TSPCP, &
                                  TA, QA, PCPR, RRATE, SRATE, PS, &
                                  IPCP, N, I, J)
!=====================RLi ====================================== \
!     4             ZRH )
!=====================RLi ====================================== /
!
         CALL  radiationDriver(FFC, FG, FCS, FGS, ALVSCN, ALIRCN, & ! Formerly CLASSA
                               ALVSG,  ALIRG,  ALVSCS, ALIRCS, ALVSSN, ALIRSN, &
                               ALVSGC, ALIRGC, ALVSSC, ALIRSC, TRVSCN, TRIRCN, &
                               TRVSCS, TRIRCS, FSVF,   FSVFS, &
                               RAICAN, RAICNS, SNOCAN, SNOCNS, FRAINC, FSNOWC, &
                               FRAICS, FSNOCS, DISP,   DISPS,  ZOMLNC, ZOMLCS, &
                               ZOELNC, ZOELCS, ZOMLNG, ZOMLNS, ZOELNG, ZOELNS, &
                               CHCAP,  CHCAPS, CMASSC, CMASCS, CWLCAP, CWFCAP, &
                               CWLCPS, CWFCPS, RC,     RCS,    RBCOEF, FROOT, &
                               ZPLIMC, ZPLIMG, ZPLMCS, ZPLMGS, TRSNOW, ZSNOW, &
                               ZWSNOW, ALVS,   ALIR,   ZHTCC,  ZHTCS,  ZHTC, &
                               ZWTRC,  ZWTRS,  ZWTRG,  ZCMAI,  ZFSNOW, &
                               ZFCANMX,ZZOLN,  ZALVSC, ZALIRC, ZPAIMAX,ZPAIMIN, &
                               ZCWGTMX,ZZRTMAX,ZRSMIN, ZQA50,  ZVPDA,  ZVPDB, &
                               ZPSIGA, ZPSIGB, PAIDAT, HGTDAT, ACVDAT, ACIDAT, &
                               ASVDAT, ASIDAT, AGVDAT, AGIDAT, ZALGWET,ZALGDRY, &
                               THLIQ,  THICE,  TS,     ZRCAN,  ZSCAN,  ZTCAN, &
                               ZGROWTH,XSNO,   ZTSNOW, ZRHOSNO,ZALBSNO,ZBLEND, &
                               Z0ORO,  SNOLIM, ZPLMG0, ZPLMS0, &
                               FCLOUD, TA,     VPD,    RHOAIR, COSZS, &
                               QSWINV, ZDLAT,  ZDLON,  RHOSNI, DELZ,   ZDELZW, &
                               ZZBOTW, ZTHPOR, ZTHLMIN,ZPSISAT,ZBI,    ZPSIWLT, &
!     N                         ZHCPS,  ISAND,  IDAY,   N,      I,      J,
                               ZHCPS,  ISAND, &
!     P                         ZFCANCMX,ICC,CTEM1,CTEM2,ZRMATC,ZZOLNC,ZPAIC,
!     Q                         ZAILC,ZAILCG,ZCMASVEGC,L2MAX, NOL2PFTS,ZSLAIC,
!     S                         ZAILCGS,  ZFCANCS, ZFCANC,
                               ZFCANCMX,ICC,ictemmod,ZRMATC,ZZOLNC,ZCMASVEGC, &
                               ZAILC,ZPAIC,L2MAX, NOL2PFTS,ZSLAIC, &
                               ZAILCG,   ZAILCGS,  ZFCANC, ZFCANCS, &
                               IDAY,   N,      I,      J, &
                               TRNCH,kount,  IC,     ICP1,   IG,     IDISP,  IZREF, &
                               IWF,    IPAI,   IHGT,   IALC,   IALS,   IALG, &
                               ZALVSCTM, ZALIRCTM)
!
!
         call energyBudgetDriver(TBARC, TBARG, TBARCS, TBARGS, THLIQC, THLIQG, & ! Formerly CLASST
                                 THICEC, THICEG, HCPC,   HCPG,   TCTOPC, TCBOTC, TCTOPG, TCBOTG, &
                                 GZEROC, GZEROG, GZROCS, GZROGS, G12C,   G12G,   G12CS,  G12GS, &
                                 G23C,   G23G,   G23CS,  G23GS,  QFREZC, QFREZG, QMELTC, QMELTG, &
                                 EVAPC,  EVAPCG, EVAPG,  EVAPCS, EVPCSG, EVAPGS, TCANO,  TCANS, &
                                 RAICAN, SNOCAN, RAICNS, SNOCNS, CHCAP,  CHCAPS, TPONDC, TPONDG, &
                                 TPNDCS, TPNDGS, TSNOCS, TSNOGS, WSNOCS, WSNOGS, RHOSCS, RHOSGS, &
                                 ITERCT, ZCDH,   ZCDM,   QSENS,  TFLUX,  QEVAP,  EVAPO,  QFLUX, &
!     8                           EVPPOT, ZCFLUX, EVAPB,  ZTSRAD, QS,     ZTSURF, ST,     SU,
                                 EVPPOT, ZCFLUX, EVAPB,  ZTSRAD, QS,     ZTSURF, &
!     9                           SV,     SQ,     ZFSGV,  ZFSGS,  ZFSGG,  ZFLGV,  ZFLGS,  ZFLGG,
                                 ST,     SU,     SV,     SQ,     srh, &
                                 ZFSGV,  ZFSGS,  ZFSGG,  ZFLGV,  ZFLGS,  ZFLGG, &
                                 ZHFSC,  ZHFSS,  ZHFSG,  ZHEVC,  ZHEVS,  ZHEVG,  ZHMFC,  ZHMFN, &
                                 ZHTCC,  ZHTCS,  ZHTC,   ZQFCF,  ZQFCL,  CDRAG,  WTABLE, ZILMO, &
                                 ZFRV,   HBL, ZTAIRCAN, ZHUAIRCAN, ZZUSL, ZZTSL, ZUN,    ZTN, &
                                 VPD,    TADP,   RHOAIR, QSWINV, QSWINI, QLWIN,  UA,     VA, &
                                 TA,     QA,     PADRY,  FFC,    FG,     FCS,    FGS,    RBCOEF, &
!     F                           FSVF,   FSVFS,  ALVSCN, ALIRCN, ALVSG,  ALIRG,
                                 FSVF,   FSVFS,  PS,     vmod,   ALVSCN, ALIRCN, ALVSG,  ALIRG, &
                                 ALVSCS, ALIRCS, ALVSSN, ALIRSN, ALVSGC, ALIRGC, ALVSSC, ALIRSC, &
!     H                           TRVSCN, TRIRCN, TRVSCS, TRIRCS, RC,     RCS,    ZWTRG,
                                 TRVSCN, TRIRCN, TRVSCS, TRIRCS, RC,     RCS,    ZWTRG,  qlwavg, &
                                 FRAINC, FSNOWC, FRAICS, FSNOCS, CMASSC, CMASCS, DISP,   DISPS, &
                                 ZOMLNC, ZOELNC, ZOMLNG, ZOELNG, ZOMLCS, ZOELCS, ZOMLNS, ZOELNS, &
                                 TS,     THLIQ,  THICE,  ZTPOND, ZZPOND, ZTBASE, ZTCAN,  ZTSNOW, &
                                 ZSNOW,  TRSNOW, ZRHOSNO,ZWSNOW, ZTHPOR, ZTHLRET,ZTHLMIN,ZTHFC, &
                                 ZDLAT,  PCPR,   ZHCPS,  ZTCS,   ZTSFS,  DELZ,   ZDELZW, ZZBOTW, &
                                 ZFTEMP, ZFVAP,  ZRIB, &
                                 ISAND,  Isat,   igwscheme, &
!===================== RLi =====================================\
                                 ZAILCG,ZAILCGS,ZFCANC,ZFCANCS,ZCO2CONC,ZCO2I1CG, &
                                 ZCO2I1CS,ZCO2I2CG,ZCO2I2CS,COSZS,XDIFFUS,ZSLAI, &
                                 ICC,    ictemmod, ZRMATCTEM,ZFCANCMX, L2MAX, &
                                 NOL2PFTS,  ZCFLUXCG,  ZCFLUXCS, &
!    ------------ CTEM INPUTS ABOVE THIS LINE, OUTPUTS BELOW -----------|
!     U                           ZANCSVEG,ZANCGVEG, ZRMLCSVEG,ZRMLCGVEG,ZCANRES,
                                 ZANCSVEG,ZANCGVEG, ZRMLCSVEG,ZRMLCGVEG, &
!     +                           fieldsm,wiltsm,
                                 bterm,zpsisat,zgrksat, &
                                 ITC,    ITCG,   ITG,    N,      I,  J,  TRNCH,  KOUNT, &
                                 IC,IG,  IZREF,  ISLFD,  NLANDCS,NLANDGS,NLANDC, NLANDG, NLANDI)
!     V                           ZRH,ZSAND,ZCLAY,ZORGM,PS )
!    Note: SAND, CLAY and ORGM are not modified by this subroutine
!
!===================== RLi =====================================/
!
!     Compute long-wave upward flux from surface temperature
!     using Stefan-Boltzmann law
!

         call waterBudgetDriver(THLIQ, THICE, TS, ZTCAN, ZRCAN, ZSCAN, & ! Formerly CLASSW
                                ZRUNOFF,ZTRUNOFF,XSNO,  ZTSNOW, ZRHOSNO,ZALBSNO, &
                                ZWSNOW, ZZPOND, ZTPOND, ZGROWTH,ZTBASE, GFLUX, &
                                ZPCFC,  ZPCLC,  ZPCPN,  ZPCPG,  ZQFCF,  ZQFCL, &
                                ZQFN,   ZQFG,   ZQFC,   ZHMFC,  ZHMFG,  ZHMFN, &
                                ZHTCC,  ZHTCS,  ZHTC,   ZROFC,  ZROFN,  ZROVG, &
                                ZWTRS,  ZWTRG,  ZOVRFLW,ZSUBFLW,ZBASFLW, &
                                ZTOVRFL,ZTSUBFL,ZTBASFL,EVAPO, &
                                TBARC,  TBARG,  TBARCS, TBARGS, THLIQC, THLIQG, &
                                THICEC, THICEG, HCPC,   HCPG,   RPCP,   TRPCP, &
                                SPCP,   TSPCP,  PCPR,   TA,     RHOSNI, ZGGEO, &
                                FFC,    FG,     FCS,    FGS,    TPONDC, TPONDG, &
                                TPNDCS, TPNDGS, EVAPC,  EVAPCG, EVAPG,  EVAPCS, &
                                EVPCSG, EVAPGS, QFREZC, QFREZG, QMELTC, QMELTG, &
                                RAICAN, SNOCAN, RAICNS, SNOCNS, FROOT,  FSVF, &
                                FSVFS,  CWLCAP, CWFCAP, CWLCPS, CWFCPS, TCANO, &
                                TCANS,  CHCAP,  CHCAPS, CMASSC, CMASCS, ZSNOW, &
                                GZEROC, GZEROG, GZROCS, GZROGS, G12C,   G12G, &
                                G12CS,  G12GS,  G23C,   G23G,   G23CS,  G23GS, &
                                TSNOCS, TSNOGS, WSNOCS, WSNOGS, RHOSCS, RHOSGS, &
                                ZPLIMC, ZPLIMG, ZPLMCS, ZPLMGS, ZTSFS, &
                                TCTOPC, TCBOTC, TCTOPG, TCBOTG, &
                                ZTHPOR, ZTHLRET,ZTHLMIN,ZBI,    ZPSISAT,ZGRKSAT, &
                                ZTHLRAT,ZTHFC,  ZXDRAIN,ZHCPS,  DELZ, &
                                ZDELZW, ZZBOTW, ZXSLOPE,ZGRKFAC,ZWFSURF,ZWFCINT, &
!     P                          ISAND,  IWF,    N,      I,      J,      KOUNT,
                                ISAND,  igdr, &
                                IWF,    N,      I,      J,      KOUNT, &
                                TRNCH,  IC,     IG,     IG+1,   IG+2, &
                                NLANDCS,NLANDGS,NLANDC, NLANDG, NLANDI, ZWT, &
                                QDIST, QINT, QINC,QING,QINGS,QINCS,QDISC,Isat, &
                                ZWTNEW,DMLIQT,Ibed,Isats,XWT,ZEXCW,ZSLP,ZARE,ZLEG, &
                                ZANIS, igwscheme)
!!!                                 ZDDENS, ZINTFLOW, ZSANI, INTERFLOW_L) !HUZIY ADDED ZDDENS and ZINTFLOW FOR INTERFLOW
         if (igwscheme.ne.1) then
            do ik=i,j
               ZLBEDR(IK)=Ibed(IK)
            enddo
         endif

      endif ! if (.not.dotile.and.doprec)
!
      doprec=dotile
!
    end do CLASS_MAIN_LOOP



!Huziy: debug start
    call check_if_isand_psisat_are_consistent(ISAND, ZPSISAT, N, IG, ok_status)
    IF (.not. ok_status) THEN
       print*,"ISAND and psisat are not"// &
              "consistent after CLASS_MAIN_LOOP,"// &
              "kount = ", KOUNT
    ENDIF
!Huziy: debug finish



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
          ALVIS_SOL(I) = 0.5*(ALVS(I)+ALIR(I))
          FSOLUACC(I)  = FSOLUACC(I)+FLUSOL(I)*ALVIS_SOL(I)*DT
          FIRUACC(I)   = FIRUACC(I) +QLWAVG(I) * DT
          ZFL(I)       = FFC(I)*GZEROC(I)+FG(I)*GZEROG(I)+FCS(I)*GZROCS(I)+ &
                         FGS(I)*GZROGS(I)
!
    end do FINAL_CALC
!
!  end do DO_MOSAIC
!


  ztsurf   (1:n) => bus( x(tsurf,1,1) : )
  st       (1:n) => bus( x(tdiag,1,1) : )
  ztt2m    (1:n) => bus( x(tdiagtyp,1,indx_sfc)  : )

  zTT2m(:)    = st(:)                     ! instantaneos



! FILL THE ARRAYS TO BE AGGREGATED LATER IN S/R AGREGE
  CALL FILLAGG ( BUS, BUSSIZ, PTSURF, PTSURFSIZ, INDX_SOIL, SURFLEN )
!
  DEALLOCATE( PAIDAT,HGTDAT,ACVDAT,ACIDAT, &
              TBARC ,TBARG ,TBARCS,TBARGS, &
              THLIQC,THLIQG,THICEC,THICEG, &
              HCPC  ,HCPG  ,FROOT ,GFLUX , &
              TCTOPC,TCBOTC,TCTOPG,TCBOTG, &
              ISAND, IORG,  ITERCT )
!
  RETURN
END
