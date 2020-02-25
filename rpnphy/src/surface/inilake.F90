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

subroutine inilake(ni, trnch)
   use sfc_options, only: schmlake, ltran0
   use sfcbus_mod
   use flake_parameters, only: H_Ice_max
   implicit none
#include <arch_specific.hf>
#include <rmnlib_basics.hf>

   integer ni, trnch

   !@Author Katja Winger (February 2017)
   !@Object Initialize the lake parameters
  !@Arguments
   !             - Input -
   ! NI          longueur d'une tranche horizontale

   include "sfcinput.cdk"

   integer :: i, k
   real, parameter :: TMELS = 273.15
   real, parameter :: TMELI = 273.05
   real, pointer, dimension(:)   :: zlakedepth, zlaketransp
   real, pointer, dimension(:)   :: zfrv_li, zfrv_lw, zlakect, zlakehice, zlakehml, zlakefice, zlaketbot, zlaketice, zlaketmnw, zlaketwml
   real, pointer, dimension(:,:) :: zlaketp, zvegf
   real, pointer, dimension(:)   :: ztwater, zglsea, zicedp, ztmice

#define MKPTR1D(NAME1,NAME2)     nullify(NAME1); if (vd%NAME2%i > 0 .and. associated(busptr(vd%NAME2%i)%ptr)) NAME1(1:ni)                                 => busptr(vd%NAME2%i)%ptr(:,trnch)
#define MKPTR2D(NAME1,NAME2)     nullify(NAME1); if (vd%NAME2%i > 0 .and. associated(busptr(vd%NAME2%i)%ptr)) NAME1(1:ni,1:vd%NAME2%mul*vd%NAME2%niveaux) => busptr(vd%NAME2%i)%ptr(:,trnch)
#define MKPTR2DL(NAME1,NAME2,NN) nullify(NAME1); if (vd%NAME2%i > 0 .and. associated(busptr(vd%NAME2%i)%ptr)) NAME1(1:ni,1:NN)                            => busptr(vd%NAME2%i)%ptr(:,trnch)

   MKPTR1D(zlakedepth,lakedepth)
   MKPTR1D(zlaketransp,laketransp)
   MKPTR1D(zfrv_li,frv_li)
   MKPTR1D(zfrv_lw,frv_lw)
   MKPTR1D(zlakect,lakect)
   MKPTR1D(zlakehice,lakehice)
   MKPTR1D(zlakehml,lakehml)
   MKPTR1D(zlakefice,lakefice)
   MKPTR1D(zlaketbot,laketbot)
   MKPTR1D(zlaketice,laketice)
   MKPTR1D(zlaketmnw,laketmnw)
   MKPTR1D(zlaketwml,laketwml)

   MKPTR1D(ztwater,twater)
   MKPTR1D(zglsea,glsea)
   MKPTR1D(zicedp,icedp)
   MKPTR1D(ztmice,tmice)

   MKPTR2D(zvegf,vegf)

!  ------

   ! Initial values (if not read from initial conditions): 
   ! no snow, no ice cover
   ! Surface water temp. = 10C
   ! 1M mixed layer
   ! Bottom temp. = 4C

   ! Lake depth
   if (any('lakedepth' == phyinread_list_s(1:phyinread_n)))  then
      do i=1,ni
         if (zlakedepth(i) < 30.) zlakedepth(i) = 100.    ! Here zlakedepth is still in decimeters!!!
      enddo
      zlakedepth = zlakedepth*0.1  ! Convert read zlakedepth from decimeters to meters
   ! From here on zlakedepth is in meters!
   elseif (any('vegf' == phyinread_list_s(1:phyinread_n))) then
      do i=1,ni
         if (zvegf(i,3) <= 0.5) then
           zlakedepth(i) = 10.
         else
           zlakedepth(i) = 60.
         endif
      end do
   else
      zlakedepth = 30.
   endif

   ! Lake water transparency
   if (all('laketransp' /= phyinread_list_s(1:phyinread_n))) &
      zlaketransp = ltran0

   ! Lake ice fraction
   if (all('lakefice' /= phyinread_list_s(1:phyinread_n))) then
      if (any('glsea' == phyinread_list_s(1:phyinread_n))) then
         zlakefice = zglsea
      else
         zlakefice = 0.0
      endif
   endif

   ! Lake ice thickness
   if (all('lakehice' /= phyinread_list_s(1:phyinread_n))) then
      do i=1,ni
         if (any('icedp' == phyinread_list_s(1:phyinread_n))) then
            zlakehice(i)   = min(zicedp(i),H_Ice_max)
         elseif (zlakefice(i) >= 0.0) then
            zlakehice(i)   = 0.1
         else
            zlakehice(i)   = 0.0
         endif
      enddo
   endif

   ! Lake ice surface temperature
   if (all('laketice' /= phyinread_list_s(1:phyinread_n))) then
      if (any('tmic' == phyinread_list_s(1:phyinread_n))) then
         zlaketice   = ztmice
      else
         zlaketice   = TMELS
      endif
   endif

   ! Lake mixed layer thickness
   if (all('lakehml'  /= phyinread_list_s(1:phyinread_n))) &
      zlakehml(i)    = 1.0

   if     (schmlake == 'FLAKE') then
      ! Limit lake depth to 60m
      do i=1,ni
         zlakedepth(i) = min(zlakedepth(i), 60.)
      enddo
      ! Lake ice surface friction velocity
      if (all('frv_li'   /= phyinread_list_s(1:phyinread_n))) &
         zfrv_li     = 0.01
      ! Lake water surface friction velocity
      if (all('frv_lw'   /= phyinread_list_s(1:phyinread_n))) &
         zfrv_lw     = 0.01
      ! Lake mixed layer temperature in FLake
      if (all('laketwml' /= phyinread_list_s(1:phyinread_n))) then
         if (any('twater' == phyinread_list_s(1:phyinread_n))) then
            zlaketwml   = ztwater
         else
            zlaketwml   = TMELS + 7.5 ! Summer: TMELS + 10.0; Winter: TMELS + 5.0
         endif
      endif
      ! Lake bottom temperature in FLake
      if (all('laketbot' /= phyinread_list_s(1:phyinread_n))) &
         zlaketbot   = TMELS +  4.0
      ! Lake water average temperature in FLake
      if (all('laketmnw' /= phyinread_list_s(1:phyinread_n))) &
         zlaketmnw   = TMELS +  5.0   ! Summer: TMELS +  5.0; Winter: TMELS + 4.5
      ! Lake shape factor (thermocline) in FLake
!print *,'inilake zlaketwml : ',zlaketwml(:)
!print *,'inilake zlaketmnw : ',zlaketmnw(:)
!print *,'inilake zlaketbot : ',zlaketbot(:)
!print *,'inilake zlakehml  : ',zlakehml(:)
!print *,'inilake zlakedepth: ',zlakedepth(:)
      if (all('lakect'   /= phyinread_list_s(1:phyinread_n))) then
         do i=1,ni
            if (abs(zlaketwml(i)-zlaketbot(i)) .gt. 0.1) then
               zlakect(i)     = min(0.8,(zlaketwml(i) - zlaketmnw(i))  &
                              / (zlaketwml(i) - zlaketbot(i))          &  ! Crashes when lake average temperature equals bottom temperature
                              / (1.-zlakehml(i)/zlakedepth(i)))
            else
              zlakect(i)     = 0.8
            endif
!print *,'inilake zlakect: ',zlakect(i)
         enddo
      endif
   endif


   ! Set lake fields to zero in points without a lake fraction
   if (any('vegf' == phyinread_list_s(1:phyinread_n))) then
      do i=1,ni
         if (zvegf(i,3) <= 0.0) then
            zlakedepth(i)  = 0.
            zlaketransp(i) = 0.
            zlakefice(i)  = 0.
            zlakehice(i)  = 0.
            zlakehice(i)  = 0.
            zfrv_li(i)    = 0.
            zfrv_lw(i)    = 0.
            zlakehml(i)   = 0.
            if     (schmlake == 'FLAKE') then
               zlaketwml(i) = 0.
               zlaketbot(i) = 0.
               zlaketmnw(i) = 0.
               zlakect(i)   = 0.
            endif
         endif
      enddo
   endif
   return
end subroutine inilake
