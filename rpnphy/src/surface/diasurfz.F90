SUBROUTINE DIASURFZ(UDIAG,VDIAG,TDIAG,QDIAG,N,  &
                    TA,QA,UA,VA,TG,QG,Z0M,Z0T,ILMO,ZUA,ZTA, &
                    H,UE,FTEMP,FVAP,ZUG,ZTG,LAT,FCOR,FG,IL1,IL2,JL)
!SUBROUTINE DIASURFZ(UDIAG,VDIAG,TDIAG,QDIAG,N,  &
!                    UA,VA,TG,QG,Z0,Z0T,ILMO,ZA, &
!                    H,UE,FTEMP,FVAP,ZU,ZT,LAT,FG,IL1,IL2,JL)

  use sfclayer_mod, only: sl_prelim, sl_sfclayer, SL_OK

  IMPLICIT NONE
  INTEGER N,JL,IL1,IL2
  REAL ZTG(N),ZUG(N)
  REAL UDIAG(N),VDIAG(N),TDIAG(N),QDIAG(N),ZUA(N),ZTA(N)
  REAL TA(N),QA(N),UA(N),VA(N)
  REAL TG(N),QG(N),UE(N),FTEMP(N),FVAP(N)
  REAL ILMO(N),Z0T(N),Z0M(N),H(N),FG(N)
  REAL LAT(N),FCOR(N)

!
!Author
!          K. Winger (Mar 2022)
!Revision
!Object
!          Interface to routine sl_sfclayer
!          - Reorganize fields into temporary arrays with only valid points
!          - Aggregate output fields
!
!Arguments
!
!          - Output -
! UDIAG    U component of the wind at Z=ZUG
! VDIAG    V component of the wind at Z=ZUG
! TDIAG    temperature in kelvins at Z=ZTG
! QDIAG    specific humidity at Z=ZTG
!
!          - Input -
! N        number of points to process
! UA       U component of wind at Z=ZUA
! VA       V component of wind at Z=ZUA
! TA       Lowest level potential temperature (K)
! QA       Lowest level specific humidity (kg/kg)
! TG       temperature at the surface (skin) (Z=0) in Kelvins
! QG       specific humidity at the surface (skin)
! PS       surface pressure at the surface
! ILMO     inverse of MONIN-OBUKHOV lenth
! H        height of boundary layer
! UE       friction velocity
! Z0M      roughness lenth for winds
! Z0T      roughness lenth for temperature and moisture
! FTEMP    (flux_t) temperature flux at surface
! FVAP     (flux_q) vapor flux at surface
! ZUA      (hghtm_air)  height of the lowest momentum level (m)
! ZTA      (hghtt_air)  height of the lowest thermodynamic level (m)
! ZUG      (hghtm_diag) height for computation of wind components
! ZTG      (hghtt_diag) height for computation of temperature and moisture
! LAT      Latitude (rad)
! FCOR     Coriolis factor (/s)
! FG       Fraction of surface type being studied 
!

  ! Temporary input fields
  real, dimension (:)  , allocatable :: tmp_TA , tmp_QA , tmp_UA , tmp_VA
  real, dimension (:)  , allocatable :: tmp_TG , tmp_QG
  real, dimension (:)  , allocatable :: tmp_ZUA, tmp_ZTA, tmp_ZUG, tmp_ZTG
  real, dimension (:)  , allocatable :: tmp_Z0M, tmp_Z0T, tmp_LAT, tmp_FCOR

  ! Temporary output fields
  real, dimension (:)  , allocatable :: tmp_TDIAG, tmp_QDIAG, tmp_UDIAG, tmp_VDIAG

  ! Work fields
  real, dimension (:)  , allocatable :: VMOD, VDIR, PS
  integer :: i, err, tmp_i, tmp_n, real_i(N)

  ! TDIAGLIM: Limit temperature inversion in lowest layer [get 'tdiaglim']
  logical, parameter :: TDIAGLIM_FALSE = .false.

! -------------------------------------------------------------------------------------

  ! Determine number of points to be treated
  tmp_n = 0
  do i=IL1,IL2
     if (FG(i) > 0. ) then
        tmp_n = tmp_n + 1
     endif
  enddo

  ! If there are no points to be treated => return
  if ( tmp_n == 0 ) return

  ! Allocate temporary input fields which will only contain points to be treated
  allocate (tmp_TA(tmp_n) , tmp_QA(tmp_n) , tmp_UA(tmp_n) , tmp_VA(tmp_n))
  allocate (tmp_TG(tmp_n) , tmp_QG(tmp_n))
  allocate (tmp_ZUA(tmp_n), tmp_ZTA(tmp_n), tmp_ZUG(tmp_n), tmp_ZTG(tmp_n))
  allocate (tmp_Z0M(tmp_n), tmp_Z0T(tmp_n), tmp_LAT(tmp_n), tmp_FCOR(tmp_n))

  ! Allocate temporary output fields which will only contain points to be treated
  allocate (tmp_TDIAG(tmp_n), tmp_QDIAG(tmp_n), tmp_UDIAG(tmp_n), tmp_VDIAG(tmp_n))

  ! Allocate work fields
  allocate (VMOD(tmp_n), VDIR(tmp_n), PS(tmp_n))

  ! Fill temporary input fields with only points that are to be treated
  tmp_i = 0
  do i=IL1,IL2
     if (FG(i) > 0.0) then
        tmp_i           = tmp_i + 1
        real_i(tmp_i)   = i
        tmp_TA(tmp_i)   = TA(i)
        tmp_QA(tmp_i)   = QA(i)
        tmp_UA(tmp_i)   = UA(i)
        tmp_VA(tmp_i)   = VA(i)
        tmp_TG(tmp_i)   = TG(i)
        tmp_QG(tmp_i)   = QG(i)
        tmp_ZUA(tmp_i)  = ZUA(i)
        tmp_ZTA(tmp_i)  = ZTA(i)
        tmp_ZUG(tmp_i)  = ZUG(i)
        tmp_ZTG(tmp_i)  = ZTG(i)
        tmp_Z0M(tmp_i)  = Z0M(i)
        tmp_Z0T(tmp_i)  = Z0T(i)
        tmp_LAT(tmp_i)  = LAT(i)
        tmp_FCOR(tmp_i) = FCOR(i)
     endif
  enddo


  ! Preliminaries
  ! Calculate wind modulus and direction
  PS = 100000.  ! Dummy field: mandatory but not needed
  err = sl_prelim(tmp_TA, tmp_QA, tmp_UA, tmp_VA, PS, tmp_ZUA, &
                  spd_air=VMOD, dir_air=VDIR, min_wind_speed=sqrt(2.5))
!                 spd_air=VMOD, dir_air=VDIR, min_wind_speed=sqrt(vamin))
  if (err /= SL_OK) then
     call physeterror('diasurfz', 'error returned by sl_prelim()')
     return
  endif

!print*,'N,tmp_n:',N,tmp_n
!print*,'tmp_TA:',tmp_TA
!print*,'tmp_QA:',tmp_QA
!print*,'VMOD:',VMOD
!print*,'VDIR:',VDIR
!print*,'tmp_ZUA:',tmp_ZUA
!print*,'tmp_ZTA:',tmp_ZTA
!print*,'tmp_TG:',tmp_TG
!print*,'tmp_QG:',tmp_QG
!print*,'tmp_Z0M:',tmp_Z0M
!print*,'tmp_Z0T:',tmp_Z0T
!print*,'tmp_LAT:',tmp_LAT
!print*,'tmp_FCOR:',tmp_FCOR
!print*,'tmp_ZUG:',tmp_ZUG
!print*,'tmp_ZTG:',tmp_ZTG

  ! Calculate diagnostic levels via sl_sfclayer
  err = sl_sfclayer(tmp_TA, tmp_QA, VMOD   , VDIR   , tmp_ZUA, tmp_ZTA,            &
                    tmp_TG, tmp_QG, tmp_Z0M, tmp_Z0T, tmp_LAT, tmp_FCOR,           &
                    hghtm_diag_row=tmp_ZUG, hghtt_diag_row=tmp_ZTG,                &
                    t_diag=tmp_TDIAG, q_diag=tmp_QDIAG,                            &
                    u_diag=tmp_UDIAG, v_diag=tmp_VDIAG,                            &
                    tdiaglim=TDIAGLIM_FALSE)
                    
  if (err /= SL_OK) then
     call physeterror('diasurfz', 'error returned by sl_sfclayer()')
     return
  endif

  ! Deallocate temporary input fields
  deallocate (tmp_TA , tmp_QA , tmp_UA , tmp_VA)
  deallocate (tmp_TG , tmp_QG)
  deallocate (tmp_ZUA, tmp_ZTA, tmp_ZUG, tmp_ZTG)
  deallocate (tmp_Z0M, tmp_Z0T, tmp_LAT, tmp_FCOR)

  ! Deallocate work fields
  deallocate (VMOD, VDIR)

  ! Put points in output fields in correct position
  ! and aggregate the different surface fractions at the same time 
  ! (There is one call to diasurfz per CLASS surface one fraction)
  do i=tmp_n,1,-1
     TDIAG(real_i(i)) = TDIAG(real_i(i)) + tmp_TDIAG(i)*FG(real_i(i))
     QDIAG(real_i(i)) = QDIAG(real_i(i)) + tmp_QDIAG(i)*FG(real_i(i))
     UDIAG(real_i(i)) = UDIAG(real_i(i)) + tmp_UDIAG(i)*FG(real_i(i))
     VDIAG(real_i(i)) = VDIAG(real_i(i)) + tmp_VDIAG(i)*FG(real_i(i))
  enddo

  ! Deallocate temporary output fields
  deallocate (tmp_TDIAG, tmp_QDIAG, tmp_UDIAG, tmp_VDIAG)

end subroutine diasurfz

