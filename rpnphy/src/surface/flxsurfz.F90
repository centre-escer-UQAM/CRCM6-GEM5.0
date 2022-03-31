SUBROUTINE FLXSURFZ(CDM, CDH, CFLUX, RIB, FTEMP, FVAP, ILMO, &
                    UE, FCOR, TA , QA , ZU, ZT, VMOD,        &
                    TG , QG , H , Z0M , Z0T,                 &
                    LZZ0, LZZ0T, FM, FH,N,IL1,IL2,FG,ITER,JL )

  use sfclayer_mod, only: sl_sfclayer,SL_OK

  IMPLICIT NONE

  INTEGER N,IL1,IL2,ITER(N),JL
  REAL CDM(N),CDH(N),CFLUX(N),RIB(N),FCOR(N),ILMO(N)
  REAL FTEMP(N),FVAP(N),TA(N),QA(N),ZU(N),VMOD(N)
  REAL TG(N),QG(N),H(N),Z0M(N),UE(N),ZT(N)
  REAL Z0T(N),LZZ0(N),LZZ0T(N)
  REAL fm(N),fh(N)
  REAL FG(N)

!
!Author
!          K. Winger (Mar 2022)
!Revision
!Object
!          Interface to routine sl_sfclayer
!          Reorganize fields into temporary arrays with only valid points
!
!Arguments
!
!          - Output -
! CDM      transfer coefficient of momentum times transfer coefficient of momentum
! CDH      transfer coefficient of momentum times transfer coefficient of heat
! CFLUX    transfer coefficient of momentum times temperature times wind
! RIB      bulk Richardson number
! FTEMP    temperature flux
! FVAP     vapor flux
! ILMO     (1/length of Monin-Obukov)
! UE       friction velocity 
! H        height of the boundary layer
! FM       momentum stability function
! FH       heat stability function
! LZZ0     log ((zu+z0m)/z0m)
! LZZ0T    log ((zt+z0m)/z0t)
!          - Input -
! FCOR     Coriolis factor
! ZU       (hghtm_air) height of the lowest momentum level (m)
! ZT       (hghtt_air) height of the lowest thermodynamic level (m)
! TA       potential temperature at first predictive level above surface
! QA       specific humidity     "    "      "        "      "     "
! VMOD     wind speed            "    "      "        "      "     "
! TG       surface temperature
! QG       specific humidity at the surface
! Z0M      roughness length for momentum      flux calculations
! Z0T      roughness length for heat/moisture flux calculations
! N        horizontal dimension
!

  ! Temporary input fields
  real, dimension (:)  , allocatable :: tmp_TA   , tmp_QA , tmp_VMOD
  real, dimension (:)  , allocatable :: tmp_ZU   , tmp_ZT
  real, dimension (:)  , allocatable :: tmp_TG   , tmp_QG
  real, dimension (:)  , allocatable :: tmp_Z0M  , tmp_Z0T, tmp_FCOR
  ! Temporary output fields
  real, dimension (:)  , allocatable :: tmp_CMU  , tmp_CTU, tmp_RIB
  real, dimension (:)  , allocatable :: tmp_ILMO , tmp_UE , tmp_H
  real, dimension (:)  , allocatable :: tmp_FTEMP, tmp_FVAP
  real, dimension (:)  , allocatable :: tmp_LZZ0 , tmp_LZZ0T
  real, dimension (:)  , allocatable :: tmp_FM   , tmp_FH

  integer :: i, err, tmp_i, tmp_n, real_i(N)
  real    :: CM, VA

! -------------------------------------------------------------------------------------

  ! Determine number of points to be treated
  tmp_n = 0
  do i=IL1,IL2
     if (FG(i) > 0. .and. ITER(i) == 1) then
        tmp_n = tmp_n + 1
     endif
  enddo

  ! If there are no points to be treated => return
  if ( tmp_n == 0 ) return

  ! Allocate temporary input fields which will only contain points to be treated
  allocate (tmp_TA(tmp_n)   , tmp_QA(tmp_n) , tmp_VMOD(tmp_n))
  allocate (tmp_ZU(tmp_n)   , tmp_ZT(tmp_n))
  allocate (tmp_TG(tmp_n)   , tmp_QG(tmp_n))
  allocate (tmp_Z0M(tmp_n)  , tmp_Z0T(tmp_n),tmp_FCOR(tmp_n))
  ! Allocate temporary output fields which will only contain points to be treated
  allocate (tmp_CMU(tmp_n)  , tmp_CTU(tmp_n))
  allocate (tmp_RIB(tmp_n)  , tmp_ILMO(tmp_n))
  allocate (tmp_FTEMP(tmp_n), tmp_FVAP(tmp_n))
  allocate (tmp_UE(tmp_n)   , tmp_H(tmp_n))
  allocate (tmp_LZZ0(tmp_n) , tmp_LZZ0T(tmp_n))
  allocate (tmp_FM(tmp_n)   , tmp_FH(tmp_n))

  ! Fill temporary input fields with only points that are to be treated
  tmp_i = 0
  do i=IL1,IL2
     if (FG(i) > 0.0 .and. iter(i) == 1) then
        tmp_i           = tmp_i + 1
        real_i(tmp_i)   = i
        tmp_TA(tmp_i)   = TA(i)
        tmp_QA(tmp_i)   = QA(i)
        tmp_VMOD(tmp_i) = VMOD(i)
        tmp_ZU(tmp_i)   = ZU(i)
        tmp_ZT(tmp_i)   = ZT(i)
        tmp_TG(tmp_i)   = TG(i)
        tmp_QG(tmp_i)   = QG(i)
        tmp_Z0M(tmp_i)  = Z0M(i)
        tmp_Z0T(tmp_i)  = Z0T(i)
        tmp_FCOR(tmp_i) = FCOR(i)
     endif
  enddo

  ! Pass VMOD as VDIR and FCOR as LAT because both are mandatory but neither is needed
  err = sl_sfclayer(tmp_TA, tmp_QA, tmp_VMOD , tmp_VMOD , tmp_ZU  , tmp_ZT,         &
                    tmp_TG, tmp_QG, tmp_Z0M, tmp_Z0T, tmp_FCOR, tmp_FCOR,           &
                    coefm=tmp_CMU   , coeft=tmp_CTU  , rib=tmp_RIB , ilmo=tmp_ILMO, &
                    flux_t=tmp_FTEMP, flux_q=tmp_FVAP, ue=tmp_UE   , h=tmp_H,       &
                    lzz0m=tmp_LZZ0  , lzz0t=tmp_LZZ0T, stabm=tmp_FM, stabt=tmp_FH)

  if (err /= SL_OK) then
     call physeterror('flxsurfz', 'Error in sl_sfclayer')
     return
  endif

  ! Deallocate temporary input fields
  deallocate (tmp_TA, tmp_QA, tmp_VMOD, tmp_ZU , tmp_ZT)
  deallocate (tmp_TG, tmp_QG, tmp_Z0M , tmp_Z0T, tmp_FCOR)

  ! Put points in output fields in correct position
  ! and "calculate" CFLUX, product of surface drag coefficient and wind speed
  do i=tmp_n,1,-1
     CM = tmp_CMU(i) / tmp_UE(i)
     VA = tmp_CMU(i) / CM / CM      ! Wind speed
     CDM(real_i(i))   = tmp_CMU(i) / VA
     CDH(real_i(i))   = tmp_CTU(i) / VA
     RIB(real_i(i))   = tmp_RIB(i)
     FTEMP(real_i(i)) = tmp_FTEMP(i)
     FVAP(real_i(i))  = tmp_FVAP(i)
     ILMO(real_i(i))  = tmp_ILMO(i)
     UE(real_i(i))    = tmp_UE(i)
     H(real_i(i))     = tmp_H(i)
     LZZ0(real_i(i))  = tmp_LZZ0(i)
     LZZ0T(real_i(i)) = tmp_LZZ0T(i)
     FM(real_i(i))    = tmp_FM(i)
     FH(real_i(i))    = tmp_FH(i)

     ! Product of surface drag coefficient and wind speed
     CFLUX(real_i(i)) = max(tmp_CTU(i), 0.)
  enddo

  ! Deallocate temporary output fields
  deallocate (tmp_CMU  , tmp_CTU ,  tmp_RIB, tmp_ILMO)
  deallocate (tmp_FTEMP, tmp_FVAP , tmp_UE , tmp_H)
  deallocate (tmp_LZZ0 , tmp_LZZ0T, tmp_FM , tmp_FH)

end subroutine flxsurfz
