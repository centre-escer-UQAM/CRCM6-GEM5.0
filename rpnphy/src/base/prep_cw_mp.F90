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
!** S/P PREP_CW

subroutine prep_cw_mp(f, fsiz, d, dsiz, v, vsiz, ficebl, kount, trnch, task, ni, nkm1)

   use phy_options
   use phybus

   implicit none

#include <arch_specific.hf>

  integer fsiz, dsiz, vsiz, ni, n, nkm1
  integer kount, trnch, task
  real, target :: f(fsiz), d(dsiz), v(vsiz)
  real ficebl(ni,nkm1)


!Author
!          D. Paquin-Ricard (June 2017)
!          P. Vaillancourt (July  2016)
!
!Object
!          When MP schemes are used; merge water contents from "implicit cloud sources"  for radiation
!
!Arguments
!
!          - Input -
! dsiz     dimension of d
! fsiz     dimension of f
! vsiz     dimension of v
! ficebl   fraction of ice
! kount    index of timestep
! trnch    number of the slice
! task     task number
! ni       horizontal dimension
! nk       vertical dimension
!
!          - Output -
!
!          - Input/Output -
! d        dynamic             bus
! f        permanent variables bus
! v        volatile (output)   bus
!***********************************************************************

  integer ik, i, k
  real, target,dimension(ni,nkm1) :: zero
  real, pointer, dimension(:,:) :: zfbl, zfdc, zfsc, zftot, zfxp, zfmp, ziwcimp, zlwcimp, &
       zqldi, zqlsc, zqsdi, zqssc, zqtbl

  zero   (1:ni,1:nkm1) =  0.0
!  zftot  (1:ni,1:nkm1) => f(ftot:)
!  zfxp   (1:ni,1:nkm1) => f(fxp:) !DPR-june2017 to calculate ftot, summation of fmp and fxp
!  zfmp   (1:ni,1:nkm1) => f(fmp:) !DPR-june2017 use fmp as the variable for sommation over implicit CF
!  zlwcimp(1:ni,1:nkm1) => f(lwcimp:)
!  ziwcimp(1:ni,1:nkm1) => f(iwcimp:)
  nullify(zftot);    if (ftot > 0)   zftot(1:ni,1:nkm1) => f(ftot:)
  nullify(zfxp);     if (fxp > 0)    zfxp(1:ni,1:nkm1) => f(fxp:)
  nullify(zfmp);     if (fmp > 0)    zfmp(1:ni,1:nkm1) => f(fmp:)
  nullify(zlwcimp);  if (lwcimp > 0) zlwcimp(1:ni,1:nkm1) => f(lwcimp:)
  nullify(ziwcimp);  if (iwcimp > 0) ziwcimp(1:ni,1:nkm1) => f(iwcimp:)

!
  if (any(convec == (/'KFC     ', 'KFC2    ', 'BECHTOLD'/))) then
     nullify(zqldi);    if (qldi > 0) zqldi  (1:ni,1:nkm1) => f(qldi:)
     nullify(zqsdi);    if (qsdi > 0) zqsdi  (1:ni,1:nkm1) => f(qsdi:)
     nullify(zfdc);     if (fdc > 0)  zfdc   (1:ni,1:nkm1) => f(fdc:)
  else
     zqldi => zero(1:ni,1:nkm1)
     zqsdi => zero(1:ni,1:nkm1)
     zfdc => zero(1:ni,1:nkm1)
  endif
!
  if (conv_shal /= 'NIL') then
     nullify(zqlsc);   if (qlsc > 0) zqlsc  (1:ni,1:nkm1) => v(qlsc:)
     nullify(zqssc);   if (qssc > 0) zqssc  (1:ni,1:nkm1) => v(qssc:)
     nullify(zfsc);    if (fsc > 0)  zfsc   (1:ni,1:nkm1) => f(fsc:)
  else
     zqlsc => zero(1:ni,1:nkm1)
     zqssc => zero(1:ni,1:nkm1)
     zfsc => zero(1:ni,1:nkm1)
  endif

  if(fluvert=='MOISTKE') then
     nullify(zfbl);   if (fbl > 0) zfbl   (1:ni,1:nkm1) => f(fbl:)
     nullify(zqtbl);  if (qtbl > 0) zqtbl  (1:ni,1:nkm1) => f(qtbl:)
  else
     zfbl => zero(1:ni,1:nkm1)
     zqtbl => zero(1:ni,1:nkm1)
  endif

!
! PV-avril2016 : simplified version for MP schemes only
!              - agregates implicit clouds only;
!              - assumes that moistke is an implicit source of clouds amongst others (eliminate choice between mtke and exp clouds)
!              - could choose a maximum overlap of clouds instead of random?
!              - condensates are agregated assuming max overlap while fractions assume random ???
!
!     Add the cloud water (liquid and solid) coming from PBL, shallow  and deep cumulus clouds
!     note that all condensates must be GRID-SCALE values (not in-cloud)

     if(fluvert=='MOISTKE') then
        do k=1,nkm1
           do i=1,ni
                 zlwcimp(i,k) =  zqtbl(i,k) * (1.0 - ficebl(i,k) )
                 ziwcimp(i,k) =  zqtbl(i,k) * ficebl(i,k)
           enddo
        enddo
     endif

     do k=1,nkm1
        do i=1,ni
           zlwcimp(i,k) = zlwcimp(i,k) + zqldi(i,k) + zqlsc(i,k)
           ziwcimp(i,k) = ziwcimp(i,k) + zqsdi(i,k) + zqssc(i,k)
        enddo
     enddo


!     Agregate implicit cloud fractions
!         FBL is for PBL clouds from moistke
!         FDC is for deep convection clouds
!         FSC is for the shallow convection clouds
!     DPR-june2017 FMP is the sum over the implicit, FXP: the explicit, FTOT: total
  do k=1,nkm1
     do i=1,ni
        zfmp(i,k) = min(  max( &                                      ! random overlap
             1. - (1.-zfbl(i,k))*(1.-zfdc(i,k))*(1.-zfsc(i,k)) , &
             0.   ),1.)
!        zfmp(i,k) = min(max( zfbl(i,k)+zfdc(i,k)+zfsc(i,k) ,0.) ,1.)  ! maximum overlap
        zftot(i,k)=min(max(zfmp(i,k)+zfxp(i,k),0.),1.) ! DPR-june2017 maximum overlap
!       zftot(i,k)=min(max( 1.-(1.-zfmp(i,k))*(1.-zfxp(i,k)) ,0.),1.) ! DPR-june2017random overlap

     enddo
  enddo

end subroutine prep_cw_mp
