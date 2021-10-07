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

!/@*
subroutine inisurf4(kount, ni, nk, trnch)
   use sfc_options
   use sfcbus_mod
   use svs_configs
   use class_configs
   implicit none
#include <arch_specific.hf>
#include <rmnlib_basics.hf>
#include <msg.h>
   !@Object Transfer and initialize geophysical fields for the
   !        surface schemes
   !@Arguments
   !       - Input/Ouput -
   ! f        field for permanent physics variables
   ! fsiz     dimension of f
   ! e        field for entry variables
   ! esiz     dimension of e
   ! ni       horizontal dimension

   integer ni, nk, kount, trnch

   !@Author Stephane Belair (February 1999)
   !@Revisions
   ! 001 K. Winger (UQAM/ESCER)(Sep 2019 - Add lake fraction
   ! 002 K. Winger (UQAM/ESCER) Jun 2020 - Section for CLASS added
   !@NOTE: This subroutine expects snow depth in cm.
   !       The snow depth is converted in metre (in this s/r)
   !       when the 'entry variables' are transfered to the
   !       permanent variables.
   !*@/

#include "tdpack_const.hf"
   include "isbapar.cdk"
   include "sfcinput.cdk"

   real, parameter :: z0ice = 0.001
   real, parameter :: z0sea = 0.001

   real, save :: almin  = 0.50
   real, save :: tauf   = 0.24
   real, save :: tauday = 24.

   real    :: tempsum, tempclay, tempsand, land_frac
   integer :: i, k, nk1, l, tsoil_id,wsoil_id,isoil_id, sand_id,clay_id
   real*8  :: sum_poids_8

   real, pointer, dimension(:) :: &
        zdrainaf, zemisr, zemistg, zemistgen, zfvapliqaf, zglacier, zglsea, &
        zglsea0, zicedp, ziceline, zlhtg, zmg, zml, zresa, zresagr, zresavg, &
        zresasa, zresasv, zslop, zsnoal, zsnoalen, zsnoagen, zsnodpl, zsnoden, &
        zsnoma, zsnoro, zsnvden, zsnvdp, zsnvma, ztsrad, ztwater, zwveg, &
        zwsnow, zz0en, zz0veg, zz0tveg
   real, pointer, dimension(:,:) :: &
        zalvis, zclay, zclayen, zisoil, zrunofftotaf, zdraintotaf, zsand, zsanden, zsnodp, &
        ztglacier, ztmice, ztmoins, ztsoil, zvegf, zwsoil, zz0, zz0t
   ! pointers added for lake fraction
   real, pointer, dimension(:) :: &
        zlakefr
   ! pointers added for CLASS
   real, pointer, dimension(:)   :: zsdepth, zxdrain, zz0oro
   real, pointer, dimension(:,:) :: zrootdp, zorgm, zmcmai, zmexcw

   logical, save :: CLASS_nml_read  = .false.

#define MKPTR1D(NAME1,NAME2) nullify(NAME1); if (vd%NAME2%i > 0 .and. associated(busptr(vd%NAME2%i)%ptr)) NAME1(1:ni) => busptr(vd%NAME2%i)%ptr(:,trnch)
#define MKPTR2D(NAME1,NAME2) nullify(NAME1); if (vd%NAME2%i > 0 .and. associated(busptr(vd%NAME2%i)%ptr)) NAME1(1:ni,1:vd%NAME2%mul*vd%NAME2%niveaux) => busptr(vd%NAME2%i)%ptr(:,trnch)

   !-------------------------------------------------------------

   !  Nothing to process if no fields were read
   if (phyinread_n == 0) return

   MKPTR1D(zdrainaf,drainaf)
   MKPTR1D(zemisr,emisr)
   MKPTR1D(zemistg,emistg)
   MKPTR1D(zemistgen,emistgen)
   MKPTR1D(zfvapliqaf,fvapliqaf)
   MKPTR1D(zglacier,glacier)
   MKPTR1D(zglsea,glsea)
   MKPTR1D(zglsea0,glsea0)
   MKPTR1D(zicedp,icedp)
   MKPTR1D(ziceline,iceline)
   MKPTR1D(zlhtg,lhtg)
   MKPTR1D(zmg,mg)
   MKPTR1D(zml,ml)
   MKPTR1D(zresa,resa)
   MKPTR1D(zresagr,resagr)
   MKPTR1D(zresavg,resavg)
   MKPTR1D(zresasa,resasa)
   MKPTR1D(zresasv,resasv)
   MKPTR1D(zslop,slop)
   MKPTR1D(zsnoal,snoal)
   MKPTR1D(zsnoalen,snoalen)
   MKPTR1D(zsnoagen,snoagen)
   MKPTR1D(zsnoden,snoden)
   MKPTR1D(zsnodpl,snodpl) 
   MKPTR1D(zsnoma,snoma)
   MKPTR1D(zsnoro,snoro)
   MKPTR1D(zsnvden,snvden)
   MKPTR1D(zsnvdp,snvdp) 
   MKPTR1D(zsnvma,snvma)
   MKPTR1D(ztsrad,tsrad)
   MKPTR1D(ztwater,twater)
   MKPTR1D(zwveg,wveg)
   MKPTR1D(zwsnow,wsnow)
   MKPTR1D(zz0en,z0en)
   MKPTR1D(zz0veg,z0veg)
   MKPTR1D(zz0tveg,z0tveg)

   MKPTR2D(zalvis,alvis)
   MKPTR2D(zclay,clay)
   MKPTR2D(zclayen,clayen)
   MKPTR2D(zdraintotaf,draintotaf)
   MKPTR2D(zisoil,isoil)
   MKPTR2D(zrunofftotaf,runofftotaf)
   MKPTR2D(zsand,sand)
   MKPTR2D(zsanden,sanden)
   MKPTR2D(zsnodp,snodp)
   MKPTR2D(ztglacier,tglacier)
   MKPTR2D(ztmice,tmice)
   MKPTR2D(ztmoins,tmoins)
   MKPTR2D(ztsoil,tsoil)
   MKPTR2D(zvegf,vegf)
   MKPTR2D(zwsoil,wsoil)
   MKPTR2D(zz0,z0)
   MKPTR2D(zz0t,z0t)
   
   ! for lakes
   MKPTR1D(zlakefr,lakefr)

   ! for CLASS
   if (schmsol == 'CLASS') then
      MKPTR1D(zsdepth,sdepth)
      MKPTR1D(zxdrain,xdrain)
      MKPTR1D(zz0oro,z0oro)

      MKPTR2D(zrootdp,rootdp)
      MKPTR2D(zorgm,orgm)
      MKPTR2D(zmexcw,excw)
      MKPTR2D(zmcmai,cmai)
   endif


   ! Several treatments on geophysical fields valid for isba
   ! the water temperature (tm) is decreased for points where the
   ! filtering of mountains lead to an icrease of the water level
   ! (old subroutine modtmtp of gem's dynamic library)

   ! Other consistency tests ...

   if (any('snodp' == phyinread_list_s(1:phyinread_n))) then
!VDIR NODEP
      do k=1,nsurf
         do i=1,ni
            zsnodp(i,k) = max( 0., zsnodp(i,k))
         end do
      end do
   endif

   if (any('tglacier' == phyinread_list_s(1:phyinread_n))) then
!VDIR NODEP
      do i=1,ni
         ztglacier(i,1) = min( trpl, ztglacier(i,1))
         ztglacier(i,2) = min( trpl, ztglacier(i,2))
      end do
   endif

   ! From the "entry" to the "permanent" bus
   !
   !========================================================================
   !          for variables common to all surface schemes
   !========================================================================
   !
   !
!VDIR NODEP
   DO_I: do i=1,ni

      ! Make sure all surface fractions are above critmask/-lac and adjust MG to VF
      if (any('vegf' == phyinread_list_s(1:phyinread_n))) then

         ! Ocean fraction need to be twice as large as 'critmask' 
         ! so it does not completely disappear when sea ice appears
         if (zvegf(i,1) < critmask*2.) zvegf(i,1) = 0.

         ! Glacier fraction
         if (zvegf(i,2) < critmask   ) zvegf(i,2) = 0.

         ! Lake fraction
         if (zvegf(i,3) < critlac    ) then
            ! If there is an ocean fraction convert lake into ocean
            if (zvegf(i,1) > 0. ) zvegf(i,1) = zvegf(i,1) + zvegf(i,3)
            zvegf(i,3) = 0.
         endif

         ! Land fraction
         land_frac = 1. - zvegf(i,1) - zvegf(i,2) - zvegf(i,3)
         if (land_frac  < critmask   ) then
            land_frac = 0.
            zvegf(i,4:26)  = 0.
         endif

         ! Make sure the sum of all surface fractions is still 1.
         sum_poids_8 = land_frac + zvegf(i,1) + zvegf(i,2) + zvegf(i,3)
         sum_poids_8 = 1. / sum_poids_8
         do l = 1,26
           zvegf(i,l) = zvegf(i,l) * sum_poids_8
         enddo

         ! Recalculate MG
         zmg(i) = max(0., min(1., 1. - zvegf(i,1) - zvegf(i,3)))
!print *,'inisurf: i,zmg,zvegf(1),zvegf(3):',i,zmg(i),zvegf(i,1),zvegf(i,3)

      endif


      if (any('alvis' == phyinread_list_s(1:phyinread_n))) then
         nk1 = size(zalvis,2)
         zalvis(i,indx_soil   ) = zalvis(i,nk1)
         zalvis(i,indx_glacier) = zalvis(i,nk1)
         zalvis(i,indx_water  ) = zalvis(i,nk1)
         zalvis(i,indx_ice    ) = zalvis(i,nk1)
         zalvis(i,indx_agrege ) = zalvis(i,nk1)
         if (schmurb.ne.'NIL') then
            zalvis(i,indx_urb ) = zalvis(i,nk1)
         endif
         if (schmlake.ne.'NIL') then
            zalvis(i,indx_lake) = zalvis(i,nk1)
         endif
      endif

      if (kount == 0 .and. .not.any('emisr' == phyinread_list_s(1:phyinread_n))) then
         zemisr(i) = 1.
      endif

      !       --- snodp deja en metres
      if (any('snodp' == phyinread_list_s(1:phyinread_n))) then
         zsnodp(i,indx_water  ) = 0.0
      endif

      if (any('tsoil' == phyinread_list_s(1:phyinread_n)) .and. &
          all('tsrad' /= phyinread_list_s(1:phyinread_n))) then
         ztsrad(i) = ztsoil(i,1)
      endif
      if (any('z0en' == phyinread_list_s(1:phyinread_n))) then
         zz0 (i,indx_soil   ) = max(zz0en(i),z0min)
         zz0 (i,indx_glacier) = max(zz0en(i),Z0GLA)
         zz0 (i,indx_water  ) = z0sea
         zz0 (i,indx_ice    ) = z0ice
         zz0 (i,indx_agrege ) = max(zz0en(i),z0min)
         if (schmlake.ne.'NIL') then
            zz0 (i,indx_lake) = z0sea
         endif
      endif
      if (any('z0en' == phyinread_list_s(1:phyinread_n)) .and. &
          all('z0t'  /= phyinread_list_s(1:phyinread_n))) then
         zz0t(i,indx_soil   ) = max(zz0en(i),z0min)
         zz0t(i,indx_glacier) = max(zz0en(i),Z0GLA)
         zz0t(i,indx_water  ) = z0sea
         zz0t(i,indx_ice    ) = z0ice
         zz0t(i,indx_agrege ) = max(zz0en(i),z0min)
         if (schmlake.ne.'NIL') then
            zz0t(i,indx_lake) = z0sea
         endif
      endif
      if (any('z0veg' == phyinread_list_s(1:phyinread_n))) then
         zz0veg (i) = max(zz0veg(i),z0min)
         zz0tveg(i) = max(zz0veg(i),z0min)
      endif
      if (any('glsea0' == phyinread_list_s(1:phyinread_n))) then
         zglsea (i) = zglsea0(i)
      endif
      !       Mask for the lakes
      if (any('vegf' == phyinread_list_s(1:phyinread_n))) then
         zml(i) = zvegf(i,3)
         if (schmlake /= 'NIL') then
!           if (any('mgen' == phyinread_list_s(1:phyinread_n))) then
             zlakefr(i) = zml(i) !!! * zmg(i)
!           endif
         else
           zlakefr(i) = 0.
         endif
      endif
      if (kount == 0 .and. .not.icelac) ziceline(i) = 1.

      if(kount == 0) then
         !# total surface runoff
         zrunofftotaf(i,1:nsurf+1) = 0.0
         !# total drainage
         zdraintotaf(i,1:nsurf+1) = 0.0
         !# evaporation
         zfvapliqaf(i) = 0.0
      endif

   end do DO_I


   if (any('tmice' == phyinread_list_s(1:phyinread_n))) then
      do k=1,nl
         do i=1,ni
            ztmice(i,k) = min(tcdk, ztmice(i,k))
         end do
      end do
   endif


   ! Limit sea ice thickness (KW)
   if ( icemax.ge.0.0 .and. any('icedp' == phyinread_list_s(1:phyinread_n))) then
      do i=1,ni
         zicedp(i) = min(zicedp(i), icemax)
      end do
   endif


   ! Limit snow depth (KW)
   if ( snowmax.ge.0.0 .and. any('snodp' == phyinread_list_s(1:phyinread_n))) then 
      do k=1,nsurf+1
         do i=1,ni
            zsnodp(i,k) = min(zsnodp(i,k), snowmax)
         end do
      end do
   endif


   ! Set snow depth for points without specific surface fractions to 0. (KW)
   ! No snow over water
   zsnodp(:,indx_water) = 0.0
   do i=1,ni
     ! Glaciers
     if ( zvegf(i,2) .lt. critmask ) zsnodp(i,indx_glacier) = 0.0
     ! Either ocean, glacier or lake => no soil
     if ( zvegf(i,1) + zvegf(i,2) + zvegf(i,3) .ge. 1-critmask ) zsnodp(i,indx_soil   ) = 0.0
   end do


   !========================================================================
   !                             for lakes only
   !========================================================================

   call lacs4(climat, ni, trnch)

   !========================================================================
   !     Special cases

   if (any('icedp' == phyinread_list_s(1:phyinread_n))) then
!VDIR NODEP
      do i=1,ni
         !           no snow allowed in the absence of marine ice
         if (zicedp(i).lt.himin) then
            zsnodp(i,indx_ice) = 0.0
         endif
      end do
   endif

   !=========================================================================
   !                                      FOR ISBA ... FOR ISBA ... FOR ISBA
   !=========================================================================

   IF_ISBA: if (schmsol == 'ISBA') then

      if (kount == 0 .and. all('resa' /= phyinread_list_s(1:phyinread_n))) zresa(1:ni) = 50.

      ! Special operations for the snow variables
      !
      ! Careful here about the units:
      ! "snoro" is the relative density of snow, 
      !         i.e., rho_ice / rho_water (no units)
      ! "snoma" is the snow water equivalent in mm (i.e., kg / m2)
      ! "snoal" is the snow albedo determined from the snow age
      !
      ! Note that "snoag" is in hours ... (tauday also)

!VDIR NODEP
      do i=1,ni
         if (any('snoro' == phyinread_list_s(1:phyinread_n))) then
            zsnoro(i) = max(100.,zsnoro(i)) / rauw
         endif
         if (any('snoro' == phyinread_list_s(1:phyinread_n)) .or. &
              any('snodp' == phyinread_list_s(1:phyinread_n))) then
            zsnoma(i) = rauw * zsnoro(i) * zsnodp(i,indx_soil)
         endif
      end do

      ! For the albedo, there are two possibilities:
      !
      ! 1) if switch "snoalb_anl" is true, then the "i6"
      !    record in the starting standard file (snoalen) contains the snow albedo
      !
      ! 2) if switch "snoalb_anl" is false, then we use the snow age (snoagen) 
      !    to derive the snow albedo
      !
      IF_SNO_ALB: if (snoalb_anl) then

         if (any('snoalen' == phyinread_list_s(1:phyinread_n))) then
            do i=1,ni
               zsnoal(i)  =  zsnoalen(i)
            end do
         endif

      else

         ! snow albedo is determined from the snow age according to two different
         ! expressions depending if the snow pack is melting or not

         if (any('snoagen' == phyinread_list_s(1:phyinread_n)) .or. &
              any('snoalen' == phyinread_list_s(1:phyinread_n))) then
!VDIR NODEP
            do i=1,ni
               if (ztmoins(i,nk).lt.trpl) then
                  zsnoal(i)  = ansmax - todry*zsnoagen(i)/tauday
               else
                  zsnoal(i)  = (ansmax-almin) * &
                       exp( -tauf*zsnoagen(i)/tauday ) &
                       + almin
               end if
               zsnoal(i)  = max( zsnoal(i) , almin )
               zsnoal(i)  = min( zsnoal(i) , ansmax )
            end do
         endif

      end if IF_SNO_ALB

      !  Initialize the parameters that depend on vegetation
      if (any('vegf' == phyinread_list_s(1:phyinread_n)) .or. &
           (kntveg > 0 .and. mod(kount,kntveg) == 0)) then
         call inicover2(kount, ni, trnch)
      endif

      ! Sand and clay fractions of the soil are taken as simple averages 
      ! of the first 3 layers

!VDIR NODEP
      do i=1,ni
         if (any('sand' == phyinread_list_s(1:phyinread_n))) then
            zsand(i,1) = (zsand(i,1) + zsand(i,2) + zsand(i,3)) / 3.
         endif
         if (any('clay' == phyinread_list_s(1:phyinread_n))) then
            zclay(i,1) = (zclay(i,1) + zclay(i,2) + zclay(i,3)) / 3.
         endif
      end do

      ! Make sure the entry fields are coherent ...

      call coherence3(ni, trnch)

      ! Initialize the soil characteristics using the soil texture

      if (any('clay' == phyinread_list_s(1:phyinread_n)) .or. &
           any('sand' == phyinread_list_s(1:phyinread_n))) then
         call inisoili2(ni, trnch)
      endif

   end if IF_ISBA
!=========================================================================
!                                      FOR SVS  ... FOR SVS  ... FOR SVS 
!=========================================================================
!
!
   IF_SVS: IF (schmsol.EQ.'SVS') THEN
!
!VDIR NODEP
         DO i=1,ni
            if (kount == 0) then
               ! calculate snow mass
               zsnoma(i)  = zsnoden(i) * zsnodpl(i)
               zsnvma(i)  = zsnvden(i) * zsnvdp(i)


               zdrainaf(i)        = 0.0
               if ( read_emis ) &
                    zemistg(i)         = zemistgen(i)
               zresagr(i)         = 100.
               zresavg(i)         = 50.
               zresasa(i)         = 100.
               zresasv(i)         = 100.               
               ! DDeacu: Ensure that slope is positive and set its minimum value             
               if ( zmg(i).gt.critmask ) then
                  zslop(i)  = min ( max( abs( zslop(i) ) , 5.e-03 ) , 1.0 ) 
               else
                  zslop(i)  = 0.0
               endif
               
            endif      

      END DO
!
!
!                          Initialize the parameters that depend
!                          on vegetation
!
   
      if (any('vegf' == phyinread_list_s(1:phyinread_n))) then
         call inicover_svs(0, ni, trnch)
      endif
!
!
!
!
!                           Sand and clay fractions 
!
!VDIR NODEP
      kount_zero: if ( kount == 0 ) then
         soil_data: if ( soiltext == "GSDE" .or. soiltext == "SLC" &
              .or. soiltext == "SOILGRIDS" ) then 
            DO k=1,nl_stp
               DO i=1,ni
                  watmask1: if (zmg(i).lt.critmask) then
                     ! OVER WATER...
                     zsand  (i,k)    = 0.0
                     zclay  (i,k)    = 0.0
                  else
                     ! OVER LAND
                     
                     if (zsanden(i,k)+zclayen(i,k).lt.critexture) then
                        !                If no sand and clay component
                        !                attribute to these points characteristics
                        !                of typical loamy soils
                        zsand(i,k) = 35.
                        zclay(i,k) = 35.
                     else 
                        !                 Minimum of 1% of sand and clay 
                        zsand(i,k) =  max( zsanden(i,k) , 1.0) 
                        
                        zclay(i,k) =  max( zclayen(i,k) , 1.0)
                        
                        if ( zsand(i,k)+zclay(i,k).gt.100 ) then
                           ! reduce sand & clay  percentage proportionally 
                           tempsum= zsand(i,k) + zclay(i,k)
                           zsand(i,k) = zsand(i,k)/tempsum * 100.
                           zclay(i,k) = zclay(i,k)/tempsum * 100.
                        endif
                     endif
                  endif watmask1
                  
               enddo
            enddo
            ! read in texture, do coherence check 
            ! initialize soil characteristics 
            call inisoili_svs( ni, trnch )
         endif soil_data

         ! Make sure the entry fields are coherent ...

         call coherence3(ni, trnch)


      endif kount_zero ! kount =0.0

     END IF IF_SVS



   !========================================================================
   !                             for TEB only
   !========================================================================

   ! Note that TEB variables do not support reading for kount>0:  phyincread_list_s
   !  would need to be processed within initown() to implement this support.
   if (kount == 0 .and. schmurb == 'TEB') &
        call initown2(ni, nk, trnch)


   !========================================================================
   !                       for interactive lakes only
   !========================================================================

   IF_LAKES: if (schmlake /= 'NIL' .and. kount == 0) then

      call inilake(ni,trnch)

   endif IF_LAKES


   !========================================================================
   !                       for CLASS/CTEM only
   !========================================================================
!print *,'inisurf: trnch =',trnch
!print *,'inisurf: zsand(1,:):',zsand(1,:)
!print *,'inisurf: zsand(:, 1):',zsand(:, 1)
!print *,'inisurf: zsand(:,16):',zsand(:,16)
!print *,'inisurf: zclay(:, 1):',zclay(:, 1)
!print *,'inisurf: zclay(:,16):',zclay(:,16)
   IF_CLASS: if (schmsol == 'CLASS') then

      !  Initialize the parameters that depend on vegetation
      if (any('vegf' == phyinread_list_s(1:phyinread_n)) .or. &
           (kntveg > 0 .and. mod(kount,kntveg) == 0)) then
         call inicover2(kount, ni, trnch)
      endif

      ! Make sure number of soil levels read is correct
      ! -----------------------------------------------

      if (any('tsoil' == phyinread_list_s(1:phyinread_n)) ) then
      
         ! Find sand & clay id
         tsoil_id = 0
         wsoil_id = 0
         isoil_id = 0
         do k=1,phyinread_n
            if (phyinread_list_s(k) == 'tsoil') tsoil_id = k
            if (phyinread_list_s(k) == 'wsoil') wsoil_id = k
            if (phyinread_list_s(k) == 'isoil') isoil_id = k
            if (tsoil_id /= 0 .and. wsoil_id /= 0 .and. isoil_id /= 0) exit  ! All IDs found -> exit loop
         end do
         if (tsoil_id == 0 .or. wsoil_id == 0 .or. isoil_id == 0) then
            call msg(MSG_ERROR,'(inisurf) I0, I1, and/or I2 not read')
            return
         endif
   
         ! Make sure number of soil levels read is correct
         if (phyinread_list_nk(tsoil_id) /= class_ig .or. &
             phyinread_list_nk(wsoil_id) /= class_ig .or. &
             phyinread_list_nk(isoil_id) /= class_ig ) then
            call msg(MSG_ERROR,'(inisurf) I0, I1, and/or I2: wrong number of levels read')
            return
         endif
      endif

      ! Copy read sand & clay into all levels
      ! -------------------------------------
      if (any('sand' == phyinread_list_s(1:phyinread_n)) .and. &
          any('clay' == phyinread_list_s(1:phyinread_n)) ) then

         ! Find sand & clay id
         sand_id = 0
         clay_id = 0
         do k=1,phyinread_n
            if (phyinread_list_s(k) == 'sand') sand_id = k
            if (phyinread_list_s(k) == 'clay') clay_id = k
            if (sand_id /= 0 .and. clay_id /= 0) exit
         end do
!print *,'inisurf: sand_id,clay_id:',sand_id,clay_id
!print *,'inisurf: sand_nk,clay_nk:',phyinread_list_nk(sand_id),phyinread_list_nk(clay_id)
         if (sand_id == 0 .or. clay_id == 0) then
            call msg(MSG_ERROR,'(inisurf) sand and/or clay not read')
            return
         endif
   
         ! Copy read lowest read sand level into all levels below
         ! If number of levels read is 1, field is in lowest level, class_ig ...
         if (phyinread_list_nk(sand_id) == 1) then
            do i=1,ni
               zsand(i,1                           :class_ig-1) = zsand(i,class_ig)
            end do
         ! ... otherwise the field is in the top n levels
         else
            do i=1,ni
               zsand(i,phyinread_list_nk(sand_id)+1:class_ig  ) = zsand(i,phyinread_list_nk(sand_id))
            end do
         endif
   
         ! Copy read lowest read clay level into all levels below
         ! If number of levels read is 1, field is in lowest level, class_ig ...
         if (phyinread_list_nk(clay_id) == 1) then
            do i=1,ni
               zclay(i,1                           :class_ig-1) = zclay(i,class_ig)
            end do
         ! ... otherwise the field is in the top n levels
         else
            do i=1,ni
               zclay(i,phyinread_list_nk(clay_id)+1:class_ig  ) = zclay(i,phyinread_list_nk(clay_id))
            end do
         endif
!print *,'inisurf: zsand(:, 1):',zsand(:, 1)
!print *,'inisurf: zsand(:, 2):',zsand(:, 2)
!print *,'inisurf: zsand(:, 3):',zsand(:, 3)
!print *,'inisurf: zsand(:, 4):',zsand(:, 4)
!print *,'inisurf: zsand(:, 5):',zsand(:, 5)
!print *,'inisurf: zsand(:, 6):',zsand(:, 6)
!print *,'inisurf: zsand(:, 7):',zsand(:, 7)
!print *,'inisurf: zsand(:,15):',zsand(:,15)
!print *,'inisurf: zsand(:,16):',zsand(:,16)
!print *,'inisurf: zclay(:, 1):',zclay(:, 1)
!print *,'inisurf: zclay(:, 2):',zclay(:, 2)
!print *,'inisurf: zclay(:, 3):',zclay(:, 3)
!print *,'inisurf: zclay(:,15):',zclay(:,15)
!print *,'inisurf: zclay(:,16):',zclay(:,16)

         ! Initialize organic matter to zero
         zorgm  = 0.
         zmcmai = 0.
         zmexcw = 0.


         ! Adjust sand & clay values
         do k=1,class_ig
            do i=1,ni

               if (zmg(i).lt.critmask) then
                  ! OVER WATER...
                  zsand  (i,k)    = 0.0
                  zclay  (i,k)    = 0.0
               else
                  ! OVER LAND
                  if (zsand(i,k)+zclay(i,k).lt.critexture) then
                     !                If no sand and clay component
                     !                attribute to these points characteristics
                     !                of typical loamy soils
                     zsand(i,k) = 35.
                     zclay(i,k) = 35.
                  else
                     !                 Minimum of 1% of sand and clay 
                     zsand(i,k) =  max( zsand(i,k) , 1.0)
                     zclay(i,k) =  max( zclay(i,k) , 1.0)

                     ! If the sum of sand + clay is greater than 100% ...
                     if ( zsand(i,k)+zclay(i,k).gt.100 ) then
                        ! ... reduce sand & clay  percentage proportionally 
                        tempsum= zsand(i,k) + zclay(i,k)
                        zsand(i,k) = zsand(i,k)/tempsum * 100.
                        zclay(i,k) = zclay(i,k)/tempsum * 100.
                     endif
                  endif
               endif

            enddo
         enddo

         do i=1,ni

            ! orographic roughness length
            if (any('z0oro' == phyinread_list_s(1:phyinread_n))) then
               ! Re-initialize Z0M and Z0M to Z0oro
               if (any('z0en' == phyinread_list_s(1:phyinread_n))) then
                  zz0 (i,indx_soil   ) = max(zz0oro(i),z0min)
                  zz0t(i,indx_soil   ) = max(zz0en(i),z0min)
               endif

               zz0oro (i) = max(zz0oro(i),z0min)
            else
               zz0oro (i) = max(zz0en(i),z0min)
            endif
!zz0oro (i) = z0min

            ! Set minimum soil moisture to 0.04
            do k=1,class_ig
               zwsoil(i,k)= max(0.04,zwsoil(i,k))
            enddo

            ! If soil is frozen set liquid water content to 0.04
            do k=1,class_ig
               if (ztsoil(i,k).lt.273.15) zwsoil(i,k)=0.04
            enddo

            ! Set drainage index for water flow at bottom of soil profile
            zxdrain(i) = 1.0 - zvegf(i,23)

            ! Keep depth to bedrock between 0.38 and 3.0 m
            zsdepth(i) = max(0.38, zsdepth(i))
            zsdepth(i) = min(3.00, zsdepth(i))

            ! Make sure root depth does not go beyond depth to bedrock
            do k=1,class_ic
              zrootdp(i,k) = min(zrootdp(i,k), zsdepth(i))
            enddo

         enddo

      endif ! if sand/clay got read

      ! Read CLASS namelist 'CLASS_input_table'
      if (.not. CLASS_nml_read) then
         if (trnch==1) then
            call iniclass
            CLASS_nml_read = .true.
         endif
      endif

      ! Make sure the entry fields are coherent ...
      call coherence3(ni, trnch)

   endif IF_CLASS

   return
end subroutine inisurf4
