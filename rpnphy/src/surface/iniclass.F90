subroutine iniclass

  !@Object Initialize CLASS fields
  !@Author K. Winger
  !@Revisions
  ! 001 K. Winger (UQAM/ESCER) Jun 2020 - Initial version


  use sfc_options  , only : vamin, delt
  use sfclayer_mod , only : set_class_const
  use classicParams, only : runParamsFile, prepareGlobalParams, GROWYR

  implicit none
!  character(len=*), intent(in) :: F_path
!  integer, intent(in) :: F_myproc
!  integer :: F_istat
  logical, save :: first = .true.

if (first) then
  first = .false.
  print *,'In iniclass.F90'
  print *,'iniclass HELLOOOOOOOOOOOOOOOOOOOOOOOOOOO AAA'
endif

  runParamsFile = "../CLASS_input_table"

  ! Pass some constants from GEM to CLASS
  call set_class_const(delt,vamin)

  ! Set pure CLASS/CTEM constants
  call prepareGlobalParams

print *,'GROWYR:',GROWYR(1,1,:)

end subroutine iniclass
