!---------------------------------- LICENCE BEGIN -------------------------------
! GEM - Library of kernel routines for the GEM numerical atmospheric model
! Copyright (C) 1990-2010 - Division de Recherche en Prevision Numerique
!                       Environnement Canada
! This library is free software; you can redistribute it and/or modify it
! under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, version 2.1 of the License. This library is
! distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
! without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
! PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
! You should have received a copy of the GNU Lesser General Public License
! along with this library; if not, write to the Free Software Foundation, Inc.,
! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
!---------------------------------- LICENCE END ---------------------------------
!
!**s/r dcmip_supercell - Setup for supercell (DCMIP 2016)

      subroutine dcmip_supercell (F_u,F_v,F_w,F_t,F_zd,F_s,F_topo,F_q,F_pert,F_thbase,Mminx,Mmaxx,Mminy,Mmaxy,Nk)

      use canonical
      use supercell

      use geomh

      use glb_ld
      use cstv
      use lun
      use ver
      use gmm_itf_mod
      use ptopo
      implicit none

      integer Mminx,Mmaxx,Mminy,Mmaxy,Nk,F_pert
      real F_u     (Mminx:Mmaxx,Mminy:Mmaxy,Nk), &
           F_v     (Mminx:Mmaxx,Mminy:Mmaxy,Nk), &
           F_w     (Mminx:Mmaxx,Mminy:Mmaxy,Nk), &
           F_t     (Mminx:Mmaxx,Mminy:Mmaxy,Nk), &
           F_zd    (Mminx:Mmaxx,Mminy:Mmaxy,Nk), &
           F_s     (Mminx:Mmaxx,Mminy:Mmaxy),    &
           F_topo  (Mminx:Mmaxx,Mminy:Mmaxy),    &
           F_q     (Mminx:Mmaxx,Mminy:Mmaxy,Nk), &
           F_thbase(Mminx:Mmaxx,Mminy:Mmaxy,Nk)

      !object
      !==================================
      !   Setup for supercell (DCMIP 2016)
      !===================================


      !-------------------------------------------------------------------------------

      integer i,j,k,istat

      real(8) x_a_8,y_a_8,utt_8,vtt_8,s_8(2,2),rlon_8

      real(8)  :: lon,     & ! Longitude (radians)
                  lat,     & ! Latitude (radians)
                  z          ! Altitude (m)

      real(8)  :: p          ! Pressure  (Pa)

      integer  :: zcoords    ! 0 if p coordinates are specified
                             ! 1 if z coordinates are specified

      real(8)  :: u,       & ! Zonal wind (m s^-1)
                  v,       & ! Meridional wind (m s^-1)
               !!!w,       & ! Vertical Velocity (m s^-1)
                  t,       & ! Temperature (K)
                  tv,      & ! Virtual Temperature (K)
                  thetav,  & ! Virtual potential temperature (K)
               !!!phis,    & ! Surface Geopotential (m^2 s^-2)
                  ps,      & ! Surface Pressure (Pa)
                  rho,     & ! density (kg m^-3)
                  q,       & ! water vapor mixing ratio (kg/kg)
!                  dth,     & ! Theta-Theta_EQ
!                  dpr,     & ! Pressure-Pressure_EQ
                  thbaseX   ! Basic potential temperature
!                  dtt        ! TT-TT_EQ

      real(8), parameter:: ztop_8 = 20000.d0, & ! Model Top
                           zh_8   = 15000.d0    ! Altitude of the Rayleigh damped layer

      real(8) zblen_max_8,zblen_top_8,work1_8,work2_8,fact_8

      !-------------------------------------------------------------------------------

      if (Lun_out > 0) write (Lun_out,1000)

      !Factors in Rayleigh damped layer (Based on set_betav)
      !-----------------------------------------------------
      zblen_max_8 = Cstv_Zsrf_8 - zh_8/8780.2 !HREF = 8780.2 R T_ref/grav as in set_geom
      zblen_top_8 = Cstv_Ztop_8

      fact_8 = 1.0d0

      istat = gmm_get(gmmk_fcu_s,fcu)
      istat = gmm_get(gmmk_fcv_s,fcv)
      istat = gmm_get(gmmk_fcw_s,fcw)

      zcoords = 0  ! p coordinates are specified

      if (Lun_out>0) write(Lun_out,*) 'DCMIP_SUPERCELL INITIALIZATION: TEMPERATURE STARTED'

      !Initial conditions: T,ZD,W,Q,S,TOPO,THBASE
      !------------------------------------------
      do k = 1,Nk

         do j = 1,l_nj

            lat   = geomh_y_8(j)
            y_a_8 = geomh_y_8(j)

            if (Ptopo_couleur == 0) then

               do i = 1,l_ni

                  lon = geomh_x_8(i)

                  call supercell_test (lon,lat,p,z,zcoords,Ver_a_8%t(k),Ver_b_8%t(k), &
                                       Cstv_pref_8,u,v,t,tv,thetav,ps,rho,q,F_pert,thbaseX)

                  F_t   (i,j,k) = tv
                  F_q   (i,j,k) = q
                  F_s   (i,j)   = log(ps/Cstv_pref_8)
                  F_topo(i,j)   = 0. ! It is zero
                  F_zd  (i,j,k) = 0. ! It is zero
                  F_w   (i,j,k) = 0. ! It is zero

                  !Basic potential temperature
                  !---------------------------
                  F_thbase(i,j,k) = thbaseX

                  work1_8 = zblen_max_8-(Ver_a_8%t(k)+Ver_b_8%t(k)*log(ps/Cstv_pref_8))
                  work2_8 = zblen_max_8-zblen_top_8
                  work1_8 = min(1.d0,max(0.d0,work1_8/work2_8))

                  fcw(i,j,k) = work1_8*work1_8*min(1.d0,fact_8)

               end do

            else

               do i = 1,l_ni

                  x_a_8 = geomh_x_8(i) - acos(-1.D0)

                  call smat(s_8,rlon_8,lat,x_a_8,y_a_8)

                  lon = rlon_8 + acos(-1.D0)

                  call supercell_test (lon,lat,p,z,zcoords,Ver_a_8%t(k),Ver_b_8%t(k), &
                                       Cstv_pref_8,utt_8,vtt_8,t,tv,thetav,ps,rho,q,F_pert,thbaseX)

                  F_t   (i,j,k) = tv
                  F_q   (i,j,k) = q
                  F_s   (i,j)   = log(ps/Cstv_pref_8)
                  F_topo(i,j)   = 0. ! It is zero
                  F_zd  (i,j,k) = 0. ! It is zero
                  F_w   (i,j,k) = 0. ! It is zero

                  !Basic potential temperature
                  !---------------------------
                  F_thbase(i,j,k) = thbaseX

                  work1_8 = zblen_max_8-(Ver_a_8%t(k)+Ver_b_8%t(k)*log(ps/Cstv_pref_8))
                  work2_8 = zblen_max_8-zblen_top_8
                  work1_8 = min(1.d0,max(0.d0,work1_8/work2_8))

                  fcw(i,j,k) = work1_8*work1_8*min(1.d0,fact_8)

               end do

            end if

         end do

      end do

      if (Lun_out>0) write(Lun_out,*) 'DCMIP_SUPERCELL INITIALIZATION: TEMPERATURE COMPLETED'

      if (Lun_out>0) write(Lun_out,*) 'DCMIP_SUPERCELL INITIALIZATION:      U TRUE STARTED'

      !Initial conditions: U True
      !--------------------------
      do k = 1,Nk

         do j = 1,l_nj

            lat   = geomh_y_8(j)
            y_a_8 = geomh_y_8(j)

            if (Ptopo_couleur == 0) then

               do i = 1,l_niu

                  lon = geomh_xu_8(i)

                  call supercell_test (lon,lat,p,z,zcoords,Ver_a_8%m(k),Ver_b_8%m(k), &
                                       Cstv_pref_8,u,v,t,tv,thetav,ps,rho,q,F_pert,thbaseX)

                  F_u(i,j,k) = u

                  work1_8 = zblen_max_8-(Ver_a_8%m(k)+Ver_b_8%m(k)*log(ps/Cstv_pref_8))
                  work2_8 = zblen_max_8-zblen_top_8
                  work1_8 = min(1.d0,max(0.d0,work1_8/work2_8))

                  fcu(i,j,k) = work1_8*work1_8*min(1.d0,fact_8)

               end do

            else

               do i = 1,l_niu

                  x_a_8 = geomh_xu_8(i) - acos(-1.D0)

                  call smat(s_8,rlon_8,lat,x_a_8,y_a_8)

                  lon = rlon_8 + acos(-1.D0)

                  call supercell_test (lon,lat,p,z,zcoords,Ver_a_8%m(k),Ver_b_8%m(k), &
                                       Cstv_pref_8,utt_8,vtt_8,t,tv,thetav,ps,rho,q,F_pert,thbaseX)

                  u = s_8(1,1)*utt_8 + s_8(1,2)*vtt_8

                  F_u(i,j,k) = u

                  work1_8 = zblen_max_8-(Ver_a_8%m(k)+Ver_b_8%m(k)*log(ps/Cstv_pref_8))
                  work2_8 = zblen_max_8-zblen_top_8
                  work1_8 = min(1.d0,max(0.d0,work1_8/work2_8))

                  fcu(i,j,k) = work1_8*work1_8*min(1.d0,fact_8)

               end do

            end if

         end do

      end do

      if (Lun_out>0) write(Lun_out,*) 'DCMIP_SUPERCELL INITIALIZATION:      U TRUE COMPLETED'

      if (Lun_out>0) write(Lun_out,*) 'DCMIP_SUPERCELL INITIALIZATION:      V TRUE STARTED'

      !Initial conditions: V True
      !--------------------------
      do k = 1,Nk

         do j = 1,l_njv

            lat   = geomh_yv_8(j)
            y_a_8 = geomh_yv_8(j)

            if (Ptopo_couleur == 0) then

               do i = 1,l_ni

                  lon = geomh_x_8(i)

                  call supercell_test (lon,lat,p,z,zcoords,Ver_a_8%m(k),Ver_b_8%m(k), &
                                       Cstv_pref_8,u,v,t,tv,thetav,ps,rho,q,F_pert,thbaseX)

                  F_v(i,j,k) = v

                  work1_8 = zblen_max_8-(Ver_a_8%m(k)+Ver_b_8%m(k)*log(ps/Cstv_pref_8))
                  work2_8 = zblen_max_8-zblen_top_8
                  work1_8 = min(1.d0,max(0.d0,work1_8/work2_8))

                  fcv(i,j,k) = work1_8*work1_8*min(1.d0,fact_8)

               end do

            else

               do i = 1,l_ni

                  x_a_8 = geomh_x_8(i) - acos(-1.D0)

                  call smat(s_8,rlon_8,lat,x_a_8,y_a_8)

                  lon = rlon_8 + acos(-1.D0)

                  call supercell_test (lon,lat,p,z,zcoords,Ver_a_8%m(k),Ver_b_8%m(k), &
                                       Cstv_pref_8,utt_8,vtt_8,t,tv,thetav,ps,rho,q,F_pert,thbaseX)

                  v = s_8(2,1)*utt_8 + s_8(2,2)*vtt_8

                  F_v(i,j,k) = v

                  work1_8 = zblen_max_8-(Ver_a_8%m(k)+Ver_b_8%m(k)*log(ps/Cstv_pref_8))
                  work2_8 = zblen_max_8-zblen_top_8
                  work1_8 = min(1.d0,max(0.d0,work1_8/work2_8))

                  fcv(i,j,k) = work1_8*work1_8*min(1.d0,fact_8)

               end do

            end if

         end do

      end do

      if (Lun_out>0) write(Lun_out,*) 'DCMIP_SUPERCELL INITIALIZATION:      V TRUE COMPLETED'

      !-------------------------------------------------------------------------------

      return

 1000 format( &
      /,'USE INITIAL CONDITIONS FOR SUPERCELL: (S/R DCMIP_SUPERCELL)', &
      /,'===========================================================',/,/)

      end subroutine dcmip_supercell
