!-------------------------------------- LICENCE BEGIN -------------------------
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

!/@*
function sfc_nml2(F_namelist) result(F_istat)
   use str_mod, only: str_concat, str_toreal
   use sfc_options
   use sfcbus_mod
   use cpl_itf, only: cpl_nml
   use sfclayer_mod, only: sl_put, SL_OK
   implicit none
#include <arch_specific.hf>
   !@objective Set defaults values and read surface namelist
   !@Revisions
   ! 001 K. Winger (ESCER/UQAM) - Add 'schmlake' 
   !@Arguments
   ! F_namelist    File name containing the namelists to read
   character(len=*), intent(in) :: F_namelist
   !@return
   integer :: F_istat
   !*@/

#include <msg.h>
#include <rmnlib_basics.hf>
#include <WhiteBoard.hf>
#include <clib_interface_mu.hf>
   include "tebcst.cdk"

   integer, external :: msg_getUnit

   integer, parameter :: SFC_NML_ERR = RMN_ERR
   integer, parameter :: SFC_NML_OK  = RMN_OK + 1
   integer, parameter :: CPL_NML_OK  = 1

   logical :: offline
   integer :: err, unout
   !-------------------------------------------------------------------
   F_istat = SFC_NML_ERR

   err = sfc_options_init()
   if (.not.RMN_IS_OK(err)) return
   err = sfc_nml_read(F_namelist)
   if (.not.RMN_IS_OK(err)) return
   err = sfc_nml_check()
   if (.not.RMN_IS_OK(err)) return
   call sfc_nml_print()
   err = sfc_nml_post_init()
   if (.not.RMN_IS_OK(err)) return

   unout = msg_getUnit(MSG_INFO)
   err   = cpl_nml(F_namelist, unout)
   if (.not.RMN_IS_OK(err)) then
      call msg(MSG_ERROR,'(sfc_nml) Probleme in cpl_nml')
      return
   endif
   cplocn = (err ==  CPL_NML_OK)
   if (cplocn) err = cpl_nml('print', unout)

   F_istat = SFC_NML_OK
   !-------------------------------------------------------------------
   return


contains
 

   function sfc_nml_read(m_namelist) result(m_istat)
      implicit none
      character(len=*), intent(in) :: m_namelist
      integer :: m_istat
      character(len=1024) :: msg_S, namelist_S, name_S
      integer :: istat, unf
      !----------------------------------------------------------------
      m_istat = RMN_ERR
      name_S = '(sfc_nml)'
      namelist_S = '&surface_cfgs'

      unf = 0
      istat = clib_isreadok(m_namelist)
      if (RMN_IS_OK(istat)) istat = fnom(unf, m_namelist, 'SEQ+OLD', 0)
      if (.not.RMN_IS_OK(istat)) then
         write(msg_S,'(a,a,a,a)') trim(name_S), &
              ' Using default config, ', &
              'Namelist file Not found/readable: ', trim(m_namelist)
         call msg(MSG_INFO, msg_S)
         m_istat = RMN_OK
         return
      endif

      read(unf, nml=surface_cfgs, iostat=istat)
      if (istat == 0) then     !# Read ok
         m_istat = RMN_OK + 1
      else if (istat < 0) then !# EOF, nml not found
         write(msg_S,'(a,a,a,a,a)') trim(name_S), &
              ' Using default config, Namelist ',&
              trim(namelist_S), ' not available in file: ', trim(m_namelist)
         call msg(MSG_INFO, msg_S)
         m_istat = RMN_OK
      else !# Error
         write(msg_S,'(a,a,a,a,a)') trim(name_S), &
              ' Namelist', trim(namelist_S), &
              ' invalid in file: ', trim(m_namelist)
         call msg(MSG_ERROR, msg_S)
      endif
      istat = fclos(unf)
      !----------------------------------------------------------------
      return
   end function sfc_nml_read


   subroutine sfc_nml_print()
      implicit none
      integer, external :: msg_getUnit
      integer :: unout
      !----------------------------------------------------------------
      unout = msg_getUnit(MSG_INFO)
      if (unout > 0) then
         write(unout, nml=surface_cfgs)
      endif
      !----------------------------------------------------------------
      return
   end subroutine sfc_nml_print


   function sfc_nml_check() result(m_istat)
      use sfc_options, only : nl_svs_default, dp_svs_default
      use svs_configs, only : nl_svs,dl_svs
      implicit none
      integer :: m_istat

      integer :: istat,k,nk
      character(len=1024) :: msg_S
      !----------------------------------------------------------------
      m_istat = RMN_ERR

      istat = clib_toupper(diusst)
      istat = clib_toupper(ice_emiss)
      istat = clib_toupper(schmsol)
      istat = clib_toupper(schmurb)
      istat = clib_toupper(schmlake)
      istat = clib_toupper(sl_func_stab)
      istat = clib_toupper(sl_func_unstab)
      istat = clib_toupper(snow_emiss)
      istat = clib_toupper(isba_soil_emiss)
      istat = clib_toupper(soiltext)
      istat = clib_toupper(water_emiss)
      istat = clib_toupper(z0mtype)
      istat = clib_toupper(z0ttype)
      istat = clib_toupper(z0tevol)

      if (.not.any(schmsol == SCHMSOL_OPT)) then
         call str_concat(msg_S, SCHMSOL_OPT,', ')
         call msg(MSG_ERROR,'(sfc_nml_check) schmsol = '//trim(schmsol)//' : Should be one of: '//trim(msg_S))
         return
      endif

      if (.not.any(schmurb == SCHMURB_OPT)) then
         call str_concat(msg_S, SCHMURB_OPT,', ')
         call msg(MSG_ERROR,'(sfc_nml_check) schmurb = '//trim(schmurb)//' : Should be one of: '//trim(msg_S))
         return
      endif

      if (.not.any(schmlake == SCHMLAKE_OPT)) then
         call str_concat(msg_S, SCHMLAKE_OPT,', ')
         call msg(MSG_ERROR,'(sfc_nml_check) schmlake = '//trim(schmlake)//' : Should be one of: '//trim(msg_S))
         return
      endif

      if (.not.any(diusst == DIUSST_OPT)) then
         call str_concat(msg_S, DIUSST_OPT,', ')
         call msg(MSG_ERROR,'(sfc_nml_check) diusst = '//trim(diusst)//' : Should be one of: '//trim(msg_S))
         return
      endif

      if (.not.any(z0mtype == Z0MTYPE_OPT)) then
         call str_concat(msg_S, Z0MTYPE_OPT,', ')
         call msg(MSG_ERROR,'(sfc_nml_check) z0mtype = '//trim(z0mtype)//' : Should be one of: '//trim(msg_S))
         return
      endif

     if (.not.any(z0ttype == Z0TTYPE_OPT)) then
         call str_concat(msg_S, Z0TTYPE_OPT,', ')
         call msg(MSG_ERROR,'(sfc_nml_check) z0ttype = '//trim(z0ttype)//' : Should be one of: '//trim(msg_S))
         return
      endif

     if (.not.any(z0tevol == Z0TEVOL_OPT)) then
         call str_concat(msg_S, Z0TEVOL_OPT,', ')
         call msg(MSG_ERROR,'(sfc_nml_check) z0tevol = '//trim(z0tevol)//' : Should be one of: '//trim(msg_S))
         return
      endif

      if ( (min(z0tlat(1), z0tlat(2)) <  0.0)    .or. &
           (max(z0tlat(1), z0tlat(2)) > 90.0)    .or. &
           (z0tlat(1) > z0tlat(2))         )    then
         write(msg_S,*) z0tlat(1), z0tlat(2)
         call msg(MSG_ERROR, '(sfc_nml_check) zz0tlat = '//trim(msg_S)//' : should be in 0. to 90. and z0tlat(1) < z0tlat(2)')
         return
      endif

      if (.not.any(sl_func_stab == SL_FUNC_STAB_OPT)) then
         call str_concat(msg_S, SL_FUNC_STAB_OPT,', ')
         call msg(MSG_ERROR,'(sfc_nml_check) sl_func_stab = '//trim(sl_func_stab)//' : Should be one of: '//trim(msg_S))
         return
      endif
      if (sl_func_stab /= 'DELAGE97' .and. .not.sl_z0ref) then
         call msg(MSG_ERROR,'(sfc_nml_check) sl_z0ref must be .TRUE. for any class of stability functions (sl_func_stab) other than DELAGE97')
         return
      endif

      if (.not.any(sl_func_unstab == SL_FUNC_UNSTAB_OPT)) then
         call str_concat(msg_S, SL_FUNC_UNSTAB_OPT,', ')
         call msg(MSG_ERROR,'(sfc_nml_check) sl_func_unstab = '//trim(sl_func_unstab)//' : Should be one of: '//trim(msg_S))
         return
      endif
      if (sl_func_unstab /= 'DELAGE92' .and. .not.sl_z0ref) then
         call msg(MSG_ERROR,'(sfc_nml_check) sl_z0ref must be .TRUE. for any class of stability functions (sl_func_unstab) other than DELAGE92')
         return
      endif
 
      if (.not.RMN_IS_OK(str_toreal(ice_emiss_const,ice_emiss))) then
         call msg(MSG_ERROR,'(sfc_nml_check) ice_emiss = '//trim(ice_emiss)//' : Should be a floating point number between 0 and 1')
         return
      endif

      if (.not.RMN_IS_OK(str_toreal(snow_emiss_const,snow_emiss))) then
         call msg(MSG_ERROR,'(sfc_nml_check) snow_emiss = '//trim(snow_emiss)//' : Should be a floating point number between 0 and 1')
         return
      endif

      if (.not.RMN_IS_OK(str_toreal(isba_soil_emiss_const,isba_soil_emiss))) then
         if (.not.any(isba_soil_emiss == ISBA_SOIL_EMISS_OPT)) then
            call str_concat(msg_S, ISBA_SOIL_EMISS_OPT,', ')
            call msg(MSG_ERROR,'(sfc_nml_check) isba_soil_emiss = '//trim(isba_soil_emiss)//' : Should be a number or one of: '//trim(msg_S))
            return
         endif
      endif

      if (.not.RMN_IS_OK(str_toreal(water_emiss_const,water_emiss))) then
         call msg(MSG_ERROR,'(sfc_nml_check) water_emiss = '//trim(water_emiss)//' : Should be a floating point number between 0 and 1')
         return
      endif

      ! ------ CHECK SVS OPTIONS --------------
      IF_SVS: if (schmsol == 'SVS') then

         ! check number of SOIL LEVELS
         nk = 0
         do k=1,size(dp_svs)
            if (dp_svs(k) < 0.) exit
            nk = nk+1
         enddo

         ! KDP option for SVS now OBSOLETE, USERS SHOULD USE KHYD instead
         !--- MAKE MODEL ABORT IF USER SPECIFIES A VALUE WHEN USING SVS
         if (kdp /= -1) then
            write(msg_S,*) '(sfc_nml_check) KDP=',kdp,' should not be used with SVS. '//&
                 'This sfc config variable is OBSOLETE and should be removed from surface_cfgs. '//&
                 'To change last/deepest soil layer considered during the accumulation '//&
                 'of lateral flow and drainage, use surface_cfg variable KHYD instead'
            call msg(MSG_ERROR, msg_S)
            return
         endif

         ! USE DEFAULT NUMBER & DEPTH of LAYERS if NOT SPECIFIED
         ! FOR KHYD: the "last soil layer considered during the accumulation of lateral flow and drainage"
         ! (1) if using default layers, and khyd not specified then:
         ! set khyd = last SVS layer
         ! (2) if using default layer and khyd specified, use khyd
         ! (3) if users specifies soil layers, then forced user to specify khyd also.
         !#TODO: should keep only the nml check here, var setting should be move to sfc_nml_post_init
         if (nk == 0) then
            call msg(MSG_INFO,'(sfc_nml_check) SVS soil layer depth dp_svs'//&
                 ' NOT SPECIFIED by user, will use DEFAULT VALUE')
            ! default number of soil layers
            nk = nl_svs_default
            ! By default "last soil layer for accumulation of lateral flow and drainage" LEVEL = nk,
            if (khyd <= 1) then
               call msg(MSG_INFO,'(sfc_nml_check) Using default SVS configuration with '//&
                    '"last hydro drainage/lateral flow layer" level KHYD = last soil layer')
               khyd = nl_svs_default
            else
               call msg(MSG_INFO, '(sfc_nml_check) Using default SVS configuration for soil '//&
                    'layer but NOT for "last hydro drainage/lateral flow layer" KHYD')
            endif
            ! default values for depth of SVS values [in meters]
            do k=1,nk
               dp_svs(k) = dp_svs_default(k)
            enddo
         else if (nk < 3) then
            write(msg_S,*) '(sfc_nml_check) SVS requires a minimum of 3 soil layers. '//&
                 'Only ',nk, ' depths were specified in dp_svs=' , dp_svs(1:nk)
            call msg(MSG_ERROR, msg_S)
            return
         else if (nk >= 3 .and. khyd <= 1) then
            write(msg_S,*) '(sfc_nml_check)  Last hydro drainage/lateral flow layer level: '//&
                 'KHYD needs to be specified by USER because not using DEFAULT SVS set-up'
            call msg(MSG_ERROR, msg_S)
            return
         endif

         ! CHECK THAT DEPTH IS INCREASING and SAVE NUMBER OF LAYERS
         do k=1,nk-1
            if (dp_svs(k) >= dp_svs(k+1)) then
               write(msg_S,*) '(sfc_nml_check) SVS layer depths badly specified: '//&
                    'dp_svs(k)=',dp_svs(k),' (k=',k,') should be shallower than '//&
                    'dl_svs(k+1)=',dp_svs(k+1),' (k+1=',k+1,')'
               call msg(MSG_ERROR, msg_S)
               return
            endif
         enddo
         ! Save depth in dl_svs ... to be used by svs
         !#TODO: this is not a check per se, should be move to sfc_nml_post_init
         nl_svs = nk
         if (.not.allocated(dl_svs)) allocate(dl_svs(nl_svs))
         do k=1,nl_svs
            dl_svs(k) = dp_svs(k)
         enddo

         if (khyd > nl_svs) then
            write(msg_S, '(a,i3,a,i3)') '(sfc_nml_check) Last hydro drainage/lateral '//&
                 'flow layer in SVS khyd=',khyd,' cannot exceed number of SVS soil layers =', nl_svs
            call msg(MSG_ERROR, msg_S)
            return
         endif

         if (.not.any(soiltext == SOILTEXT_OPT)) then
            call str_concat(msg_S, SOILTEXT_OPT, ', ')
            call msg(MSG_ERROR, '(sfc_nml_check) soiltext = '//trim(soiltext)//&
                 ' : Should be one of: '//trim(msg_S))
            return
         endif

      endif IF_SVS

      ! ------ CHECK CLASS OPTIONS --------------
      if( schmsol == 'CLASS') then
      if ( class_ig < 3 ) then
         write(msg_S,*) class_ig
         call msg(MSG_ERROR, '(sfc_nml_check) class_ig = '//trim(msg_S)//' : should be at least 3')
         return
      endif

      if (any(schmsol_lev(1:class_ig) <= 0.0) .and. &
          .not.(class_ig == 3 .and.  all(schmsol_lev(1:class_ig) == 0.0))) then
         call msg(MSG_ERROR,'(sfc_nml_check) Invalid soil level definition: all levels should have a thickness greater than 0.')
         return
      endif

      endif



      m_istat = RMN_OK
      !----------------------------------------------------------------
      return
   end function sfc_nml_check


   function sfc_nml_post_init() result(m_istat)
      implicit none
      integer :: m_istat

#include "tdpack_const.hf"

      integer :: istat, options, iverb, idxmax
      !----------------------------------------------------------------
      m_istat = RMN_ERR
 
      !# nsurf nb of surface types (schemes)
      !# phybusinit / gesdict defines nrow = nagg = nsurf+1
      idxmax = max(indx_soil, indx_glacier, indx_water, indx_ice, indx_agrege)
      if (schmurb  /= 'NIL') idxmax = max(idxmax, indx_urb )
      if (schmlake /= 'NIL') idxmax = max(idxmax, indx_lake)
      nsurf = idxmax - 1

      !# Z0tlat to radians
      z0tlat(:) = z0tlat(:) * PI/180.

      !# Put some sfc options value on th WB for the main physic to access
      options = WB_REWRITE_NONE+WB_IS_LOCAL
      iverb   = wb_verbosity(WB_MSG_INFO)
      istat = WB_OK
      istat = min(wb_put('sfc/icelac', icelac, options), istat)
      istat = min(wb_put('sfc/nclass', nclass, options), istat)
      istat = min(wb_put('sfc/nicelvl', nl, options), istat)
      istat = min(wb_put('sfc/nsurf', nsurf, options), istat)
      istat = min(wb_put('sfc/schmsol', schmsol, options), istat)
      istat = min(wb_put('sfc/schmurb', schmurb, options), istat)
      istat = min(wb_put('sfc/schmlake', schmlake, options), istat)
      istat = min(wb_put('sfc/snoalb_anl', snoalb_anl, options), istat)
      iverb = wb_verbosity(iverb)
      if (istat /= WB_OK) then
         call msg(MSG_ERROR, '(sfc_nml) Problem in wb_put')
         return
      endif

      !# Offline special case
      iverb = wb_verbosity(WB_MSG_FATAL)
      istat = wb_get('itf_phy/OFFLINE', offline)
      if (istat /= WB_OK) offline = .false.
      istat = wb_verbosity(iverb)

      !# Surface Layer module init
      istat = SL_OK
      if (istat == SL_OK) istat = sl_put('beta', beta)
      if (istat == SL_OK) istat = sl_put('rineutral',sl_rineutral)
      if (istat == SL_OK) istat = sl_put('tdiaglim',tdiaglim)
      if (istat == SL_OK) istat = sl_put('sl_stabfunc_stab',sl_func_stab)
      if (istat == SL_OK) istat = sl_put('sl_stabfunc_unstab',sl_func_unstab)
      if (istat == SL_OK) istat = sl_put('z0ref',sl_z0ref)
      if (istat /= SL_OK) then
         call msg(MSG_ERROR,'(sfc_nml) Problem initializing the surface layer module')
         return
      endif

      m_istat = RMN_OK
      !----------------------------------------------------------------
      return
   end function sfc_nml_post_init

end function sfc_nml2
