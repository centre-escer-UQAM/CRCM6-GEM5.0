!-------------------------------------a LICENCE BEGIN -------------------------
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
!-------------------------------------- LICENCE END ---------------------------

subroutine cldoppro_MP(d, dsiz, f, fsiz, v, vsiz, &
     taucs, omcs, gcs, taucl, omcl, gcl, &
     liqwcin, icewcin, &
     liqwpin, icewpin, cldfrac, &
     tt, sig, ps, m, &
     lmx, nk, nkp, mpcat, kount)

   use phy_options
   use phybus

   implicit none

#include <arch_specific.hf>
#include <rmnlib_basics.hf>

   include "nbsnbl.cdk"

   integer :: dsiz, fsiz, vsiz
   integer :: lmx, m, nk, nkp, mpcat, kount

   real, target :: d(dsiz), f(fsiz), v(vsiz)
   real :: taucs(lmx,nk,nbs), omcs(lmx,nk,nbs), gcs(lmx,nk,nbs)
   real :: taucl(lmx,nk,nbl), omcl(lmx,nk,nbl), gcl(lmx,nk,nbl)
   real :: liqwcin(lmx,nk), icewcin(lmx,nk)
   real :: liqwpin(lmx,nk), icewpin(lmx,nk)
   real :: cldfrac(lmx,nk), tt(m,nk),sig(lmx,nk),ps(lmx)


   !@Authors
   !        D. Paquin-Ricard june 2017 continuing adaptation
   !        p. vaillancourt, mai 2016 (merged simplified/adapted prep_cw_rad ,cldroppro and diagno_cw_rad)
   !
   !@Object
   !        calculate cloud optical properties for ccc radiative transfer scheme for MP schemes
   !        calculate vertical integrals of liquid/ice hydrometeors
   !        output merged fields of liquid water content,  ice water content and cloud fraction
   !
   !@Arguments
   !          - output -
   ! taucs    cloud solar optical thickness
   ! omcs     cloud solar scattering albedo
   ! gcs      cloud solar asymmetry factor
   ! taucl    cloud longwave optical thickness
   ! omcl     cloud longwave scattering albedo
   ! gcl      cloud longwave asymmetry factor
   ! topthw   total integrated optical thickness of water in the visible
   ! topthi   total integrated optical thickness of ice in the visible
   ! ctp      cloud top pressure
   ! ctt      cloud top temperature
   ! ecc      effective cloud cover (nt)
   !          - input -
   ! liqwpin  liquid water path in g/m2
   ! icewpin  solid water path in g/m2
   ! liqwcin  in-cloud liquid water content
   !          in kg water/kg air
   ! icewcin  in-cloud ice water content
   !          in kg water/kg air
   ! cldfrac  layer cloud amount (0. to 1.) (lmx,nk)
   ! tt       layer temperature (k) (m,nk)
   ! sig      sigma levels (0. to 1.) (lmx,nk; local sigma)
   ! ps       surface pressure (n/m2) (lmx)
   ! mg       ground cover (ocean=0.0,land <= 1.)  (lmx)
   ! ml       fraction of lakes (0.-1.) (lmx)
   ! lmx      number of profiles to compute
   ! m        first dimension of temperature (usually lmx)
   ! nk       number of layers
   ! mpcat    number of ice categories to consider; depends on STCOND and p3_ncat
   ! kount    time step number
   !
   !***********************************************************************
   !
   external cldoppro_data

!#include "tdpack_const.hf"
include "thermoconsts.inc"
#include "cldop.cdk"
   include "phyinput.cdk"
   include "surface.cdk"
   include "nocld.cdk"

   real, dimension(lmx,nk) :: transmissint
   real, dimension(lmx,nk) :: trans_exp,trav2d
   logical, dimension(lmx) :: top
   real, dimension(m,nk) :: aird
   real, dimension(m,nk) :: rew
   real, dimension(m,nk) :: rei
   real, dimension(m,nk) :: rec_cdd
   real, dimension(m,nk) :: vs1
   real, dimension(lmx) :: trmin
   real, dimension(lmx) :: tmem
   real, dimension(lmx) :: tot,att, z_exp
   real, dimension(lmx,nkp,nkp) :: ff
   integer, dimension(lmx     ) :: ih
   integer, dimension(lmx     ) :: ib

   !*********************************************************

   integer, parameter :: IOPTREW = 1
   integer, parameter :: IOPTREI = 2

   real, save :: third = 0.3333333 !#TODO: make this a computed const?

   integer :: i, j, k, kind, ip, l
   real :: rec_grav, cut, tmp1
   real :: dp(lmx,nk), press
   logical :: strcld, nostrlwc, readfield_L
   logical :: nocloud(lmx,nk)
   real :: rew1, rew2, rew3, dg, dg2, dg3, tausw, omsw, gsw, tausi, omsi, gsi, y1, y2, y3
   real :: taulw, omlw, glw, tauli, omli, gli
   real :: tausimp(mpcat), omsimp(mpcat), gsimp(mpcat), taulimp(mpcat), omlimp(mpcat), glimp(mpcat)
   real :: tauswmp, omswmp, gswmp
   real :: taulwmp, omlwmp, glwmp
   real :: xnu, xnu2
   real :: lwpinmp(lmx,nk), cldfmp(lmx,nk), lwcinmp(lmx,nk)
   real :: iwcinmp(lmx,nk,mpcat), iwpinmp(lmx,nk,mpcat), effradi(lmx,nk,mpcat)
   real :: iwpinmps(lmx,nk), iwcinmps(lmx,nk)

   real, pointer, dimension(:,:) :: zqcmoins, zqimoins, zqnmoins, zqgmoins
   real, pointer, dimension(:,:) :: zi1qtm, zi2qtm, zi3qtm, zi4qtm
   real, pointer, dimension(:,:) :: zeffradc, zeffradi1 , zeffradi2 , zeffradi3 , zeffradi4
   real, pointer, dimension(:) :: ztopthw,ztopthi,ztcc,zecc,zeccl,zeccm,zecch,znt
   real, pointer, dimension(:) :: zmg,zml,zctp,zctt
   real, pointer, dimension(:,:) :: ziwcimp,zlwcimp,zhumoins,ztmoins,zpmoins,zsigw,zfxp,zfmp,zftot
   real, pointer, dimension(:)   :: ztlwp, ztiwp,ztlwpin,ztiwpin
   real, pointer, dimension(:,:) :: zlwcrad, ziwcrad, zcldrad, zqcplus, zgraupel, &
        zqiplus, zsnow, zqi_cat1, zqi_cat2, zqi_cat3, zqi_cat4

   nullify(zqcmoins);  if (qcmoins > 0) zqcmoins  (1:lmx,1:nk) => d(qcmoins:)

   nullify(zeffradc);  if (effradc > 0) zeffradc(1:lmx,1:nk) => f(effradc:)
   nullify(zeffradi1); if (effradi1 > 0) zeffradi1(1:lmx,1:nk) => f(effradi1:)
   nullify(zeffradi2); if (effradi2 > 0) zeffradi2(1:lmx,1:nk) => f(effradi2:)
   nullify(zeffradi3); if (effradi3 > 0) zeffradi3(1:lmx,1:nk) => f(effradi3:)
   nullify(zeffradi4); if (effradi4 > 0) zeffradi4(1:lmx,1:nk) => f(effradi4:)

   nullify(zqimoins); if (qimoins > 0) zqimoins(1:lmx,1:nk) => d(qimoins:)
   nullify(zqnmoins); if (qnmoins > 0) zqnmoins(1:lmx,1:nk) => d(qnmoins:)
   nullify(zqgmoins); if (qgmoins > 0) zqgmoins(1:lmx,1:nk) => d(qgmoins:)

   nullify(zi1qtm); if (i1qtmoins > 0) zi1qtm(1:lmx,1:nk) => d(i1qtmoins:)
   nullify(zi2qtm); if (i2qtmoins > 0) zi2qtm(1:lmx,1:nk) => d(i2qtmoins:)
   nullify(zi3qtm); if (i3qtmoins > 0) zi3qtm(1:lmx,1:nk) => d(i3qtmoins:)
   nullify(zi4qtm); if (i4qtmoins > 0) zi4qtm(1:lmx,1:nk) => d(i4qtmoins:)

   nullify(zhumoins); if (humoins > 0) zhumoins  (1:lmx,1:nk) => d(humoins:)
   nullify(ztmoins);  if (tmoins > 0)  ztmoins   (1:lmx,1:nk) => d(tmoins:)
   nullify(zpmoins);  if (pmoins > 0)  zpmoins   (1:lmx,1:nk) => f(pmoins:)
   nullify(zsigw);    if (sigw > 0)    zsigw     (1:lmx,1:nk) => d(sigw:)
   nullify(zftot);    if (ftot > 0)    zftot     (1:lmx,1:nk) => f(ftot:)
   nullify(zfmp);     if (fmp > 0)     zfmp      (1:lmx,1:nk) => f(fmp:)
   nullify(zfxp);     if (fxp > 0)     zfxp      (1:lmx,1:nk) => f(fxp:)
   nullify(ziwcimp);  if (iwcimp > 0)  ziwcimp   (1:lmx,1:nk) => f(iwcimp:)
   nullify(zlwcimp);  if (lwcimp > 0)  zlwcimp   (1:lmx,1:nk) => f(lwcimp:)
   nullify(zqcplus);  if (qcplus > 0)  zqcplus   (1:lmx,1:nk) => d(qcplus:)

   nullify(ztopthw);  if (topthw > 0) ztopthw  (1:lmx)       => f(topthw:)
   nullify(ztopthi);  if (topthi > 0) ztopthi  (1:lmx)       => f(topthi:)
   nullify(zecc);     if (ecc > 0)    zecc     (1:lmx)       => f(ecc:)
   nullify(zeccl);    if (eccl > 0)   zeccl    (1:lmx)       => f(eccl:)
   nullify(zeccm);    if (eccm > 0)   zeccm    (1:lmx)       => f(eccm:)
   nullify(zecch);    if (ecch > 0)   zecch    (1:lmx)       => f(ecch:)
   nullify(ztcc);     if (tcc > 0)    ztcc     (1:lmx)       => f(tcc:)
   nullify(znt);      if (nt > 0)     znt      (1:lmx)       => f(nt:)
   nullify(zmg);      if (mg > 0)     zmg      (1:lmx)       => f(mg:)
   nullify(zml);      if (ml > 0)     zml      (1:lmx)       => f(ml:)
   nullify(zctp);     if (ctp > 0)    zctp     (1:lmx)       => v(ctp:)
   nullify(zctt);     if (ctt > 0)    zctt     (1:lmx)       => v(ctt:)
   nullify(ztlwp);    if (tlwp > 0)   ztlwp    (1:lmx)       => f( tlwp:)
   nullify(ztiwp);    if (tiwp > 0)   ztiwp    (1:lmx)       => f( tiwp:)
   nullify(ztlwpin);  if (tlwpin > 0) ztlwpin  (1:lmx)      => f( tlwpin:)
   nullify(ztiwpin);  if (tiwpin > 0) ztiwpin  (1:lmx)      => f( tiwpin:)
   nullify(zlwcrad);  if (lwcrad > 0) zlwcrad  (1:lmx,1:nk)  => v( lwcrad:)
   nullify(ziwcrad);  if (iwcrad > 0) ziwcrad  (1:lmx,1:nk)  => v( iwcrad:)
   nullify(zcldrad);  if (cldrad > 0) zcldrad  (1:lmx,1:nk)  => v( cldrad:)

   rec_grav=1./grav
   nostrlwc = (climat.or.stratos)
   ztlwp    = 0.0
   ztiwp    = 0.0
   ztlwpin  = 0.0
   ztiwpin  = 0.0
   zlwcrad  = 0.0
   ziwcrad  = 0.0
   zcldrad  = 0.0
   iwcinmps = 0.0
   iwcinmp  = 0.0
   lwcinmp  = 0.0

   if (kount == 0) then
      if (inilwc) then
         if (stcond /= 'NIL' ) then

            !     initialiser le champ d'eau nuageuse ainsi que la
            !     fraction nuageuse pour l'appel a la radiation a kount=0
            !     seulement.
            !     ces valeurs seront remplacees par celles calculees dans
            !     les modules de condensation.
            call cldwin(zfmp,zlwcimp,ztmoins,zhumoins,zpmoins, &
                 trav2d,zsigw,lmx,nk,satuco)
         endif
      endif
   endif


   !     begin -code previously in prep_cw_rad
   !     Initialization of cloud (from prep_cw_rad)
   readfield_L = .false.
   if ( any(dyninread_list_s == 'qc').or.&
        any(phyinread_list_s(1:phyinread_n) == 'tr/mpqc:p') )then
      readfield_L = .true.
      do k=1,nk
         do i=1,lmx
            !zlwc(i,k) = zqcplus(i,k)
            lwcinmp(i,k) = max(zqcplus(i,k),0.)
         enddo
      enddo
   endif

   if ( (stcond(1:6)=='MP_MY2').and. &
        (any(dyninread_list_s == 'mpqi') .or. any(phyinread_list_s(1:phyinread_n) == 'tr/mpqi:p')) .and. &
        (any(dyninread_list_s == 'mpqs') .or. any(phyinread_list_s(1:phyinread_n) == 'tr/mpqs:p')) .and. &
        (any(dyninread_list_s == 'mpqg') .or. any(phyinread_list_s(1:phyinread_n) == 'tr/mpqg:p')) )then
      readfield_L = .true.
      zqiplus(1:lmx,1:nk) => d(qiplus:)
      zsnow(1:lmx,1:nk) => d(qnplus:)
      zgraupel(1:lmx,1:nk) => d(qgplus:)
      !ziwc = zqiplus + zsnow
      do k=1,nk
         do i=1,lmx
            iwcinmp(i,k,1) = max(zqiplus(i,k),0.)
            iwcinmp(i,k,2) = max(zsnow(i,k),0.)
            iwcinmp(i,k,3) = max(zgraupel(i,k),0.)
         enddo
      enddo

   elseif ( (stcond=='MP_P3') .and. &
        (any(dyninread_list_s == 'i1qt') .or. any(phyinread_list_s(1:phyinread_n) == 'tr/i1qt:p')) )then

      readfield_L = .true.
      zqi_cat1(1:lmx,1:nk) => d(i1qtplus:)
      !ziwc = zqi_cat1
      do k=1,nk
         do i=1,lmx
            iwcinmp(i,k,1) = max(zqi_cat1(i,k),0.)
         enddo
      enddo

      if ( (p3_ncat >= 2).and. &
           (any(dyninread_list_s == 'i2qt') .or. any(phyinread_list_s(1:phyinread_n) == 'tr/i2qt:p')) )then

         zqi_cat2(1:lmx,1:nk) => d(i2qtplus:)
         !ziwc = ziwc + zqi_cat2
         do k=1,nk
            do i=1,lmx
               iwcinmp(i,k,2) = max(zqi_cat2(i,k),0.)
            enddo
         enddo

         if ( (p3_ncat >= 3).and. &
              (any(dyninread_list_s == 'i3qt') .or. any(phyinread_list_s(1:phyinread_n) == 'tr/i3qt:p')) )then
            zqi_cat3(1:lmx,1:nk) => d(i3qtplus:)
            !ziwc = ziwc + zqi_cat3
            do k=1,nk
               do i=1,lmx
                  iwcinmp(i,k,3) = max(zqi_cat3(i,k),0.)
               enddo
            enddo

            if ( (p3_ncat >= 4).and. &
                 (any(dyninread_list_s == 'i4qt') .or. any(phyinread_list_s(1:phyinread_n) == 'tr/i4qt:p')) )then
               zqi_cat4(1:lmx,1:nk) => d(i4qtplus:)
               !ziwc = ziwc + zqi_cat4
               do k=1,nk
                  do i=1,lmx
                     iwcinmp(i,k,4) = max(zqi_cat4(i,k),0.)
                  enddo
               enddo
            endif ! p3_ncat>= 4

         endif ! p3_ncat>= 3

      endif ! p3_ncat>= 2

   endif  ! stcond == mp_p3

   if (readfield_L) then
      do k=1,nk
         do i=1,lmx
            do l=1,mpcat
               iwcinmps(i,k)=iwcinmps(i,k)+iwcinmp(i,k,l)
            enddo
            !if (zlwc(i,k)+ziwc(i,k) .gt. 1.e-6)then
            if (lwcinmp(i,k)+iwcinmps(i,k) .gt. 1.e-6)then
               !zftot(i,k)=1
               zfxp(i,k)=1
            else
               !zftot(i,k)=0
               zfxp(i,k)=0
            endif
            !zfmp(i,k)=zftot(i,k)
         enddo
      enddo
   endif


   !     assume implicit sources have been agregated into zlwc and ziwc in subroutine prep_cw_MP

   do k=1,nk
      do i=1,lmx
         cldfrac(i,k) = zfmp(i,k)
         cldfmp(i,k)  = zfxp(i,k)
      enddo
   enddo

   !...  "no stratospheric lwc" mode when CLIMAT or STRATOS = true
   !...  no clouds above TOPC or where HU < MINQ (see nocld.cdk for topc and minq)

   if (nostrlwc) then
      do k=1,nk
         do i=1,lmx
            press = sig(i,k)*ps(i)
            if (topc.gt.press .or. minq.ge.zhumoins(i,k) ) then
               cldfmp(i,k)   = 0.0
               cldfrac(i,k)  = 0.0
               zlwcimp (i,k) = 0.0 ! implicite
               ziwcimp (i,k) = 0.0 ! implicite
               lwcinmp (i,k) = 0.0 ! explicite
               do l=1,mpcat
                  iwcinmp (i,k,l) = 0.0 !explicite
               enddo
            endif
         enddo
      enddo
   endif


   do k=1,nk
      do i=1,lmx
         liqwcin(i,k)  = max(zlwcimp(i,k), 0.)
         icewcin(i,k)  = max(ziwcimp(i,k), 0.)
         lwcinmp(i,k)  = max(zqcmoins(i,k),0.)
      enddo
   enddo

   if (stcond(1:5)=='MP_P3') then
      do k=1,nk
         do i=1,lmx
            iwcinmp(i,k,1) = max(zi1qtm(i,k),   0.)
            effradi(i,k,1) = max(zeffradi1(i,k),0.)
         enddo
      enddo
      if(p3_ncat == 2) then
         do k=1,nk
            do i=1,lmx
               iwcinmp(i,k,2) = max(zi2qtm(i,k),   0.)
               effradi(i,k,2) = max(zeffradi2(i,k),0.)
            enddo
         enddo
      endif
      if(p3_ncat == 3) then
         do k=1,nk
            do i=1,lmx
               iwcinmp(i,k,3) = max(zi3qtm(i,k),   0.)
               effradi(i,k,3) = max(zeffradi3(i,k),0.)
            enddo
         enddo
      endif
      if(p3_ncat == 4) then
         do k=1,nk
            do i=1,lmx
               iwcinmp(i,k,4) = max(zi4qtm(i,k),   0.)
               effradi(i,k,4) = max(zeffradi4(i,k),0.)
            enddo
         enddo
      endif
   endif

   if (stcond(1:6)=='MP_MY2') then
      do k=1,nk
         do i=1,lmx
            iwcinmp(i,k,1) = max(zqimoins(i,k), 0.)
            effradi(i,k,1) = max(zeffradi1(i,k),0.)
            iwcinmp(i,k,2) = max(zqnmoins(i,k), 0.)
            effradi(i,k,2) = max(zeffradi2(i,k),0.)
            iwcinmp(i,k,3) = max(zqgmoins(i,k), 0.)
            effradi(i,k,3) = max(zeffradi3(i,k),0.)
         enddo
      enddo
   endif

   do k=1,nk
      do i=1,lmx
         !Normalize water contents to get in-cloud values
         if (cldfrac(i,k) .lt. 0.01) then
            liqwcin(i,k) = 0.
            icewcin(i,k) = 0.
            cldfrac(i,k) = 0.
         else
            tmp1 = 1./cldfrac(i,k)
            liqwcin(i,k)=liqwcin(i,k)*tmp1
            icewcin(i,k)=icewcin(i,k)*tmp1
         endif
         ! remove ice water content if effective radius .ge. 1e-4,
         ! since it will not be seen by radiation
         do l=1,mpcat
            if (effradi(i,k,l) .ge. 1.e-4)   then
               iwcinmp(i,k,l) = 0.
            endif
         enddo

         if (cldfmp(i,k) .lt. 0.01) then
            lwcinmp(i,k) = 0.
            do l=1,mpcat
               iwcinmp(i,k,l) = 0.
            enddo
            cldfmp(i,k) = 0.
         else
            tmp1 = 1./cldfmp(i,k)
            lwcinmp(i,k) = lwcinmp(i,k)*tmp1
            do l=1,mpcat
               iwcinmp(i,k,l) = iwcinmp(i,k,l)*tmp1
            enddo
         endif
      enddo
   enddo

   !     calculate in-cloud liquid and ice water paths in each layer (converted to g/m^2)

   do i=1,lmx
      dp(i,1) =0.5*(sig(i,1)+sig(i,2))*rec_grav*1000.
      dp(i,nk)=1.-0.5*(sig(i,nk)+sig(i,nk-1))*rec_grav*1000.
      dp(i,1) =max(dp(i,1)*ps(i),0.)
      dp(i,nk)=max(dp(i,nk)*ps(i),0.)
   end do
   do k=2,nk-1
      do i=1,lmx
         dp(i,k)=0.5*(sig(i,k+1)-sig(i,k-1))*rec_grav*1000.
         dp(i,k)=max(dp(i,k)*ps(i),0.)
      end do
   end do

   do k=1,nk
      do i=1,lmx
         icewpin(i,k)  = icewcin(i,k)*dp(i,k)
         liqwpin(i,k)  = liqwcin(i,k)*dp(i,k)
         lwpinmp(i,k)  = lwcinmp(i,k)*dp(i,k)
         iwpinmps(i,k) = 0.
         iwcinmps(i,k) = 0.
      enddo
   enddo

   do k=1,nk
      do i=1,lmx
         do l=1,mpcat
            iwpinmp(i,k,l)=iwcinmp(i,k,l)*dp(i,k)
            iwpinmps(i,k)=iwpinmps(i,k)+iwpinmp(i,k,l)
            iwcinmps(i,k)=iwcinmps(i,k)+iwcinmp(i,k,l)
         enddo
      end do
   end do

   do k=1,nk
      do i=1,lmx
         ztlwpin(i)= ztlwpin(i)+liqwpin(i,k)+lwpinmp(i,k)
         ztiwpin(i)= ztiwpin(i)+icewpin(i,k)+iwpinmps(i,k)
      enddo
   enddo

   !FIXME: coherence of thresholds on cloud fraction vs water path vs water content everywhere in code...
   cut = 0.001
   do k=1,nk
      do i=1,lmx
         if (lwpinmp(i,k).le. 0.001 .and. iwpinmps(i,k) .le. 0.001) then
            cldfmp(i,k) = 0.
         endif
         if (liqwpin(i,k).le.0.001.and.icewpin(i,k).le.0.001) then
            cldfrac(i,k) =0.
         endif
         nocloud(i,k)=(liqwpin(i,k)+lwpinmp(i,k)+icewpin(i,k)+iwpinmps(i,k) .le. cut)
         ztlwp(i)   = ztlwp(i)   + liqwpin(i,k)*cldfrac(i,k) + lwpinmp(i,k)*cldfmp(i,k)
         ztiwp(i)   = ztiwp(i)   + icewpin(i,k)*cldfrac(i,k) + iwpinmps(i,k)*cldfmp(i,k)
         zlwcrad(i,k) = liqwcin(i,k)*cldfrac(i,k) + lwcinmp(i,k)*cldfmp(i,k)
         ziwcrad(i,k) = icewcin(i,k)*cldfrac(i,k) + iwcinmps(i,k)*cldfmp(i,k)
         !cldfrac(i,k)=  min(cldfrac(i,k)+cldfmp(i,k), 1.0)
         cldfrac(i,k)=  max(min(cldfrac(i,k)+cldfmp(i,k), 1.0),0.0)
         if(nocloud(i,k)) then
            cldfrac(i,k)=0.0
         endif
         zcldrad(i,k) = cldfrac(i,k)
      end do
   end do

   !     conversion d'unites : tlwp et tiwp en kg/m2
   do i=1,lmx
      ztlwp(i)   = ztlwp(i) * 0.001
      ztiwp(i)   = ztiwp(i) * 0.001
   enddo

   !     end -code previously in prep_cw_rad
   !
   !     initialize output fields
   !
   do i =1, lmx
      ztopthw(i) =  0.0
      ztopthi(i) =  0.0
   end do

   ! begin-code : FOR EFFECTIVE RADII OF IMPLICIT CLOUDS (NON-mp SOURCES)
   !           For calculation of effective radius of water clouds,
   !           set cdd=cloud droplet concentration per cm^3 to 100 over oceans and 500 over land
   !...   ioptrew=1 : as in newrad - from h. barker, based on aircraft data, range 4-17 micron from slingo

   do k = 1, nk
      do i =1, lmx
         if (zmg(i) .le. 0.5 .and. zml(i) .le. 0.5)                  then
            rec_cdd(i,k) =  0.01
         else
            rec_cdd(i,k) =  0.002
         endif
         aird(i,k) =  sig(i,k) * ps(i) / ( tt(i,k) * rgasd )  !aird is air density in kg/m3
         vs1(i,k) = liqwcin(i,k) * aird(i,k) * rec_cdd(i,k)
      end do
   end do
   call vspown1(rew, vs1, third, nk * lmx)
   !
   do k = 1, nk
      do i = 1, lmx
         rew(i,k) = min ( max (4., 754.6 * rew(i,k)), 17.0)
      end do
   end do
   !
   !...   choice of effective radius for implicit ice clouds
   !...   ioptrei=1 : from cccma (icewcin must be in g/m3)
   !...   ioptrei=2 : constant in microns
   !
   if (ioptrei .eq. 1)                                            then
      do k = 1, nk
         do i = 1, lmx
            vs1(i,k) = 1000. * icewcin(i,k) * aird(i,k)
         enddo
      enddo
      call vspown1(rei, vs1, 0.216, nk * lmx)

      do k = 1, nk
         do i = 1, lmx
            rei(i,k) =  max (min (83.8*rei(i,k), 50.0), 20.0)
         end do
      end do
   elseif (ioptrei .eq. 2)                                       then
      rei = 30.
   endif
   ! end-code : FOR EFFECTIVE RADII OF IMPLICIT CLOUDS (NON-mp SOURCES)
   !
   !----------------------------------------------------------------------
   !     cloud radiative properties for radiation.
   !     taucs, omcs, gcs (taucl, omcl, gcl): optical depth, single
   !     scattering albedo, asymmetry factor for solar (infrared).
   !     rew: effective radius (in micrometer) for water cloud
   !     rei: effective radius (in micrometer) for ice cloud
   !     dg: geometry length for ice cloud
   !     liqwcin  (icewcin): liquid water (ice) content (in gram / m^3)
   !     liqwpin (icewpin): liquid water (ice) path length (in gram / m^2)
   !     cloud: cloud fraction
   !     parameterization for water cloud:
   !     dobbie, etc. 1999, jgr, 104, 2067-2079
   !     lindner, t. h. and j. li., 2000, j. clim., 13, 1797-1805.
   !     parameterization for ice cloud:
   !     fu 1996, j. clim., 9, 2223-2337.
   !     fu et al. 1998 j. clim., 11, 2223-2337.
   !     NOTE:  In two Fu papers, opt prop are tested for DG from 15-130microns or REI from 10-90 microns(using dg=1.5395xrei)
   !----------------------------------------------------------------------


   DO_NBS: do j = 1, nbs
      do k = 1, nk
         do i = 1, lmx
            IF_NOCLOUD: if (nocloud(i,k))                              then
               taucs(i,k,j) =  0.
               omcs(i,k,j)  =  0.
               gcs(i,k,j)   =  0.
            else
               rew2 =  rew(i,k) * rew(i,k)
               rew3 =  rew2 * rew(i,k)
               if (liqwpin(i,k) .gt. 0.001)                          then
                  tausw =  liqwpin(i,k) * &
                       (aws(1,j) + aws(2,j) / rew(i,k) + &
                       aws(3,j) / rew2 + aws(4,j) / rew3)
                  omsw  =  1.0 - (bws(1,j) + bws(2,j) * rew(i,k) + &
                       bws(3,j) * rew2 + bws(4,j) * rew3)
                  gsw   =  cws(1,j) + cws(2,j) * rew(i,k) + &
                       cws(3,j) * rew2 + cws(4,j) * rew3
               else
                  tausw =  0.
                  omsw  =  0.
                  gsw   =  0.
               endif


               if (icewpin(i,k) .gt. 0.001)                          then
                  dg   =  1.5396 * rei(i,k)
                  dg2  =  dg  * dg
                  dg3  =  dg2 * dg
                  tausi =  icewpin(i,k) * ( ais(1,j) + ais(2,j) / dg )
                  omsi  =  1.0 - (bis(1,j) + bis(2,j) * dg + &
                       bis(3,j) * dg2 + bis(4,j) * dg3)
                  gsi   =  cis(1,j) + cis(2,j) * dg + cis(3,j) * dg2 + &
                       cis(4,j) * dg3
               else
                  tausi =  0.
                  omsi  =  0.
                  gsi   =  0.
               endif


               if (lwpinmp(i,k) .gt. 0.001)                          then
                  rew1 =  zeffradc(i,k) * 1.e+6
                  if (kount.eq.0) rew1 = 10.    ![microns] assign value at step zero (in case QC is non-zero at initial conditions)
                  rew1 = min ( max (4., rew1), 40.0)
                  rew2 =  rew1*rew1
                  rew3 =  rew1*rew1*rew1
                  tauswmp =  lwpinmp(i,k) * &
                       (aws(1,j) + aws(2,j) / rew1 + &
                       aws(3,j) / rew2 + aws(4,j) / rew3)
                  omswmp  =  1.0 - (bws(1,j) + bws(2,j) * rew1 + &
                       bws(3,j) * rew2 + bws(4,j) * rew3)
                  gswmp   =  cws(1,j) + cws(2,j) * rew1 + &
                       cws(3,j) * rew2 + cws(4,j) * rew3
               else
                  tauswmp =  0.
                  omswmp =  0.
                  gswmp   =  0.
               endif
               !

               do l=1,mpcat
                  if (iwpinmp(i,k,l) .gt. 0.001 .and. effradi(i,k,l) .lt. 1.e-4)   then
                     dg   =  1.5396 * effradi(i,k,l)*1.e6
                     if (kount.eq.0) dg = 1.5396*50.   ![microns] assign value at step zero (in case "QI" is non-zero at initial conditions)
                     dg   =  min(max(dg,15.),110.)  ! max value is lower otherwise, model is crashing
                     dg2  =  dg  * dg
                     dg3  =  dg * dg *dg
                     tausimp(l) =  iwpinmp(i,k,l) * ( ais(1,j) + ais(2,j) / dg )
                     omsimp(l)  =  1.0 - (bis(1,j) + bis(2,j) * dg + &
                          bis(3,j) * dg2 + bis(4,j) * dg3)
                     gsimp(l)   =  cis(1,j) + cis(2,j) * dg + cis(3,j) * dg2 + &
                          cis(4,j) * dg3
                  else
                     tausimp(l) =  0.
                     omsimp(l)  =  0.
                     gsimp(l)   =  0.
                  endif
               enddo

               !PV agregate SW optical properties for liq-imp + ice-imp + liq-exp + ice-exp

               ! old recipe - see cldoppro3
               !              taucs(i,k,j)  =  tausw + tausi
               !                y1          =  omsw * tausw
               !                y2          =  omsi * tausi
               !                omcs(i,k,j) = (y1 + y2) / taucs(i,k,j)
               !                gcs (i,k,j) = (y1 * gsw + y2 * gsi) / (y1 + y2)

               y1  =  tausw + tausi
               y2  =  tausw*omsw + tausi*omsi
               y3  =  tausw*omsw*gsw + tausi*omsi*gsi

               y1  =  y1 + tauswmp
               y2  =  y2 + tauswmp*omswmp
               y3  =  y3 + tauswmp*omswmp*gswmp

               do l=1,mpcat
                  y1  =  y1 +  tausimp(l)
                  y2  =  y2 +  tausimp(l)*omsimp(l)
                  y3  =  y3 +  tausimp(l)*omsimp(l)*gsimp(l)
               enddo

               taucs(i,k,j)  =  y1

               if (y1 .gt. 0.0)                            then
                  omcs(i,k,j)          =  y2/y1
                  gcs (i,k,j)          =  y3/y2
               else
                  omcs(i,k,j) =  0.
                  gcs (i,k,j) =  0.
               endif
               !
               !     calculate the optical depth for water and ice cloud in visible
               !
               if (j .eq. 1)                                         then
                  ztopthw(i) =  ztopthw(i) + tausw + tauswmp
                  ztopthi(i) =  ztopthi(i) + tausi
                  do l=1,mpcat
                     ztopthi(i) =  ztopthi(i) + tausimp(l)
                  enddo
               endif

            endif IF_NOCLOUD

         enddo
      enddo
   enddo DO_NBS

   DO_NBL: do j = 1, nbl
      do k = 1, nk
         do i = 1, lmx
            IF_NOCLOUD2: if (nocloud(i,k))                              then
               taucl(i,k,j) =  0.
               omcl(i,k,j)  =  0.
               gcl(i,k,j)   =  0.
            else
               rew2 =  rew(i,k) * rew(i,k)
               rew3 =  rew2 * rew(i,k)
               !
               if (liqwpin(i,k) .gt. 0.001)                          then
                  taulw =  liqwpin(i,k) * (awl(1,j) + awl(2,j) * rew(i,k)+ &
                       awl(3,j) / rew(i,k) + awl(4,j) / rew2 + &
                       awl(5,j) / rew3)
                  omlw  =  1.0 - (bwl(1,j) + bwl(2,j) / rew(i,k) + &
                       bwl(3,j) * rew(i,k) + bwl(4,j) * rew2)
                  glw   =  cwl(1,j) + cwl(2,j) / rew(i,k) + &
                       cwl(3,j) * rew(i,k) + cwl(4,j) * rew2
               else
                  taulw =  0.
                  omlw  =  0.
                  glw   =  0.
               endif
               !----------------------------------------------------------------------
               !     since in fu etc. the param. is for absorptance, so need a factor
               !     icewpin(i,k) / tauli for single scattering albedo
               !----------------------------------------------------------------------
               !
               if (icewpin(i,k) .gt. 0.001)        then
                  dg   =  1.5396 * rei(i,k)
                  dg2  =  dg  * dg
                  dg3  =  dg2 * dg
                  tauli =  icewpin(i,k) * (ail(1,j) + ail(2,j) / dg + &
                       ail(3,j) / dg2)
                  omli  =  1.0 - (bil(1,j) / dg + bil(2,j) + &
                       bil(3,j) * dg + bil(4,j) * dg2) * &
                       icewpin(i,k) / tauli
                  gli   =  cil(1,j) + cil(2,j) * dg + cil(3,j) * dg2 + &
                       cil(4,j) * dg3
               else
                  tauli =  0.
                  omli  =  0.
                  gli   =  0.
               endif

               if (lwpinmp(i,k) .gt. 0.001)                          then
                  rew1 =  zeffradc(i,k) * 1.e6
                  rew1 = min(max(4., rew1), 40.0)  !TEST, JM
                  if (kount.eq.0) rew1 = 10.    !assign value at step zero (in case QC is non-zero at initial conditions)
                  rew2 =  rew1*rew1
                  rew3 =  rew1*rew1*rew1
                  taulwmp =  lwpinmp(i,k) * (awl(1,j) + awl(2,j) * rew1+ &
                       awl(3,j) / rew1 + awl(4,j) / rew2 + &
                       awl(5,j) / rew3)
                  if (taulwmp<0.) then
                     print*, '** calc, taulwp: ',taulwmp,rew1,zeffradc(i,k)
                     stop
                  endif


                  omlwmp  =  1.0 - (bwl(1,j) + bwl(2,j) / rew1 + &
                       bwl(3,j) * rew1 + bwl(4,j) * rew2)
                  glwmp   =  cwl(1,j) + cwl(2,j) / rew1 + &
                       cwl(3,j) * rew1 + cwl(4,j) * rew2
               else
                  taulwmp =  0.
                  omlwmp  =  0.
                  glwmp   =  0.
               endif
               !
               do l=1,mpcat
                  if (iwpinmp(i,k,l) .gt. 0.001 .and. effradi(i,k,l) .lt. 1.e-4)   then
                     dg   =  1.5396 * effradi(i,k,l)*1.e6
                     if (kount.eq.0) dg = 1.5396*50.     !assign value at step zero (in case "QI" is non-zero at initial conditions)
                     dg   =  min(max(dg,15.),110.) ! max value is lower to avoid model crashing!
                     dg2  =  dg  * dg
                     dg3  =  dg2 * dg
                     taulimp(l) =  iwpinmp(i,k,l) * (ail(1,j) + ail(2,j) / dg + &
                          ail(3,j) / dg2)
                     omlimp(l)  =  1.0 - (bil(1,j) / dg + bil(2,j) + &
                          bil(3,j) * dg + bil(4,j) * dg2) * &
                          iwpinmp(i,k,l) / taulimp(l)
                     glimp(l)   =  cil(1,j) + cil(2,j) * dg + cil(3,j) * dg2 + &
                          cil(4,j) * dg3
                  else
                     taulimp(l) =  0.
                     omlimp(l)  =  0.
                     glimp(l)   =  0.
                  endif
               enddo
               !
               !PV agregate LW optical properties for liq-imp + ice-imp + liq-exp + ice-exp
               !

               y1  =  taulw + tauli
               y2  =  taulw*omlw + tauli*omli
               y3  =  taulw*omlw*glw + tauli*omli*gli

               y1  =  y1 + taulwmp
               y2  =  y2 + taulwmp*omlwmp
               y3  =  y3 + taulwmp*omlwmp*glwmp

               do l=1,mpcat
                  y1  =  y1 +  taulimp(l)
                  y2  =  y2 +  taulimp(l)*omlimp(l)
                  y3  =  y3 +  taulimp(l)*omlimp(l)*glimp(l)
               enddo

               taucl(i,k,j)  =  y1

               if (taucl(i,k,j)<0.) then
                  print*, '*** taucl: ',taucl(i,k,j),taulw,tauli,taulwmp,taulimp(1)
                  stop
               endif


               if (y1 .gt. 0.0)                            then
                  omcl(i,k,j)          =  y2/y1
                  gcl (i,k,j)          =  y3/y2
               else
                  omcl(i,k,j) =  0.
                  gcl (i,k,j) =  0.
               endif

            endif IF_NOCLOUD2
         enddo
      enddo
   enddo DO_NBL

   !
   !     diagnostics: cloud top pressure (ctp) and temperature (ctt)
   !     using the cloud optical depth at window region (band 6) to
   !     calculate the emissivity
   !
   !     calcul des indices IH et IB pour nuages 2-D
   !     IH = niveau le plus pres de sigma=0.4
   !     IB = niveau le plus pres de sigma=0.7
   !
   do k = 1, nk
      do i = 1, lmx
         vs1(i,k) = - 1.64872 * taucl(i,k,6)
         if (sig(i,k).le.0.4) ih(i) = k
         if (sig(i,k).le.0.7) ib(i) = k

      enddo
   enddo

   call vsexp ( trans_exp, vs1, lmx*nk )

   do i = 1, lmx
      zctp (i)   = 110000.
      zctt (i)   = 310.
      top(i) = .true.
      transmissint(i,1) = 1. - cldfrac(i,1) * (1.-trans_exp(i,1) )
      if ( (1. - transmissint(i,1)) .gt. 0.99 .and. top(i) )    then
         zctp(i) = sig(i,1)*ps(i)
         zctt(i) = tt(i,1)
         top(i) = .false.
      end if
   end do
   !
   do k = 2, nk
      do i = 1, lmx
         transmissint(i,k)=transmissint(i,k-1) * (1. - cldfrac(i,k) * &
              (1.-trans_exp(i,k) ) )
         if ( (1. - transmissint(i,k)) .gt. 0.99 .and. top(i) )  then
            zctp(i) = sig(i,k)*ps(i)
            zctt(i) = tt(i,k)
            top(i) = .false.
         end if
      end do
   end do
   !
   !...  compute total, high cloud, middle cloud and low cloud effective cloud cover (nt) as in radir7
   !     using the cloud optical depth at window region (band 6) to
   !     calculate the emissivity

   do l=1,nkp-1
      do i=1,lmx
         ff(i,l,l)=1.
         tmem(i)=1.
         trmin(i)=1.
      enddo
      ip=l+1
      do k=ip,nkp
         kind=k-2
         kind=max0(kind,1)
         do i=1,lmx
            xnu=1.-cldfrac(i,k-1)*(1.-trans_exp(i,k-1) )
            if(cldfrac(i,kind).lt.0.01) then
               tmem(i)= ff(i,l,k-1)
               trmin(i)= xnu
            else
               trmin(i)=min(trmin(i),xnu)
            endif
            ff(i,l,k)= tmem(i) * trmin(i)
         enddo
      enddo
   enddo

   do  i=1,lmx
      zecc(i)=1.-ff(i,1,nkp)
      zecch(i) = 1. - ff(i, 1   ,IH(i))
      zeccm (i) = 1. - ff(i,IH(i),IB(i))
      zeccl (i) = 1. - ff(i,IB(i),nkp   )
   enddo

   !...  compute total true cloud cover using maximum-random cloud overlap assumption

   do i=1,lmx
      ff(i,1,1)=1.
      tmem(i)=1.
      trmin(i)=1.
   enddo
   do k=2,nkp
      kind=k-2
      kind=max0(kind,1)
      do i=1,lmx
         xnu=1.-cldfrac(i,k-1)
         if(cldfrac(i,kind).lt.0.01) then
            tmem(i)= ff(i,1,k-1)
            trmin(i)= xnu
         else
            trmin(i)=min(trmin(i),xnu)
         endif
         ff(i,1,k)= tmem(i) * trmin(i)
      enddo
   enddo

   do  i=1,lmx
      ztcc(i)=1.-ff(i,1,nkp)
      ! new NT formulation: TCC*(1-exp(-0.1* total cloud optical depth for visible))
      tot(i)=ztopthw(i)+ztopthi(i)
      ! att is a tuneable factor
      att(i)=-0.1*tot(i)
   enddo

   call vsexp (z_exp,att,lmx)

   do  i=1,lmx
      znt(i)=ztcc(i)*(1.-z_exp(i))
   enddo

   return
end subroutine cldoppro_MP

