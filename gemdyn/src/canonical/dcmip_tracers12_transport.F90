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

!**s/r dcmip_tracers12_transport - Setup for 3D Hadley-like meridional circulation (DCMIP 2012)

      subroutine dcmip_tracers12_transport (F_u,F_v,F_zd,F_t,F_q,F_topo,F_s,F_q1, &
                                            Mminx,Mmaxx,Mminy,Mmaxy,Nk)

      use dcmip_2012_init_1_2_3
      use gem_options
      use geomh

      use glb_ld
      use cstv
      use lun
      use ver
      use gmm_itf_mod
      use ptopo
      implicit none

      integer Mminx,Mmaxx,Mminy,Mmaxy,Nk

      real F_u    (Mminx:Mmaxx,Mminy:Mmaxy,Nk), &
           F_v    (Mminx:Mmaxx,Mminy:Mmaxy,Nk), &
           F_zd   (Mminx:Mmaxx,Mminy:Mmaxy,Nk), &
           F_t    (Mminx:Mmaxx,Mminy:Mmaxy,Nk), &
           F_q    (Mminx:Mmaxx,Mminy:Mmaxy,Nk), &
           F_s    (Mminx:Mmaxx,Mminy:Mmaxy)   , &
           F_topo (Mminx:Mmaxx,Mminy:Mmaxy)   , &
           F_q1   (Mminx:Mmaxx,Mminy:Mmaxy,Nk)

      !object
      !===============================================================
      !   Setup for 3D Hadley-like meridional circulation (DCMIP 2012)
      !===============================================================


      !-----------------------------------------------------------------------

      integer i,j,k

      real(8) x_a_8,y_a_8,utt_8,vtt_8,s_8(2,2),rlon_8

      real(8)  :: &
                  lon,     & ! Longitude (radians)
                  lat,     & ! Latitude (radians)
                  z          ! Height (m)

      real(8)  :: p          ! Pressure  (Pa)

      integer  :: zcoords    ! 0 or 1 see below

      real(8)  :: &
                  u,       & ! Zonal wind (m s^-1)
                  v,       & ! Meridional wind (m s^-1)
                  w,       & ! Vertical Velocity (m s^-1)
                  zd,      & ! Zdot GEM
                  t,       & ! Temperature (K)
                  tv,      & ! Virtual Temperature (K)
                  phis,    & ! Surface Geopotential (m^2 s^-2)
                  ps,      & ! Surface Pressure (Pa)
                  rho,     & ! density (kg m^-3)
                  q,       & ! Specific Humidity (kg/kg)
                  q1,      & ! Tracer q1 (kg/kg)
                  time       ! Current time step

      ! if zcoords = 1, then we use z and output p
      ! if zcoords = 0, then we use p

      !-----------------------------------------------------------------------

      if (Lun_out > 0) write (Lun_out,1000)

      time = Lctl_step * Cstv_dt_8

      zcoords = 0

      !Initial conditions: T,ZD,Q,Q1,S,TOPO
      !------------------------------------
      do k = 1,Nk

         do j = 1,l_nj

            lat   = geomh_y_8(j)
            y_a_8 = geomh_y_8(j)

            if (Ptopo_couleur == 0) then

               do i = 1,l_ni

                  lon = geomh_x_8(i)

                  call test1_advection_hadley (lat,p,z,zcoords,Ver_a_8%t(k),Ver_b_8%t(k),Cstv_pref_8,u,v,w,zd,t,tv,phis,ps,rho,q,q1,time)

                  F_t   (i,j,k) = tv
                  F_q   (i,j,k) = q
                  F_s   (i,j)   = log(ps/Cstv_pref_8)
                  F_topo(i,j)   = phis
                  F_zd  (i,j,k) = zd

                  !Tracers
                  !-------
                  F_q1  (i,j,k) = q1

               end do

            else

               do i = 1,l_ni

                  x_a_8 = geomh_x_8(i) - acos(-1.D0)

                  call smat(s_8,rlon_8,lat,x_a_8,y_a_8)

                  lon = rlon_8 + acos(-1.D0)

                  call test1_advection_hadley (lat,p,z,zcoords,Ver_a_8%t(k),Ver_b_8%t(k),Cstv_pref_8,utt_8,vtt_8,w,zd,t,tv,phis,ps,rho,q,q1,time)

                  F_t   (i,j,k) = tv
                  F_q   (i,j,k) = q
                  F_s   (i,j)   = log(ps/Cstv_pref_8)
                  F_topo(i,j)   = phis
                  F_zd  (i,j,k) = zd

                  !Tracers
                  !-------
                  F_q1  (i,j,k) = q1

               end do

            end if

         end do

      end do

      !Initial conditions: U True
      !--------------------------
      do k = 1,Nk

         do j = 1,l_nj

            lat   = geomh_y_8(j)
            y_a_8 = geomh_y_8(j)

            if (Ptopo_couleur == 0) then

               do i = 1,l_niu

                  lon = geomh_xu_8(i)

                  call test1_advection_hadley (lat,p,z,zcoords,Ver_a_8%m(k),Ver_b_8%m(k),Cstv_pref_8,u,v,w,zd,t,tv,phis,ps,rho,q,q1,time)

                  F_u(i,j,k) = u

               end do

            else

               do i = 1,l_niu

                  x_a_8 = geomh_xu_8(i) - acos(-1.D0)

                  call smat(s_8,rlon_8,lat,x_a_8,y_a_8)

                  lon = rlon_8 + acos(-1.D0)

                  call test1_advection_hadley (lat,p,z,zcoords,Ver_a_8%m(k),Ver_b_8%m(k),Cstv_pref_8,utt_8,vtt_8,w,zd,t,tv,phis,ps,rho,q,q1,time)

                  u = s_8(1,1)*utt_8 + s_8(1,2)*vtt_8

                  F_u(i,j,k) = u

               end do

            end if

         end do

      end do

      !Initial conditions: V True
      !--------------------------
      do k = 1,Nk

         do j = 1,l_njv

            lat   = geomh_yv_8(j)
            y_a_8 = geomh_yv_8(j)

            if (Ptopo_couleur == 0) then

               do i = 1,l_ni

                  lon = geomh_x_8(i)

                  call test1_advection_hadley (lat,p,z,zcoords,Ver_a_8%m(k),Ver_b_8%m(k),Cstv_pref_8,u,v,w,zd,t,tv,phis,ps,rho,q,q1,time)

                  F_v(i,j,k) = v

               end do

            else

               do i = 1,l_ni

                  x_a_8 = geomh_x_8(i) - acos(-1.D0)

                  call smat(s_8,rlon_8,lat,x_a_8,y_a_8)

                  lon = rlon_8 + acos(-1.D0)

                  call test1_advection_hadley (lat,p,z,zcoords,Ver_a_8%m(k),Ver_b_8%m(k),Cstv_pref_8,utt_8,vtt_8,w,zd,t,tv,phis,ps,rho,q,q1,time)

                  v = s_8(2,1)*utt_8 + s_8(2,2)*vtt_8

                  F_v(i,j,k) = v

               end do

            end if

         end do

      end do

      !-----------------------------------------------------------------------

      return

 1000 format( &
      /,'PRESCRIBED CONDITIONS FOR 3D Hadley-like meridional circulation (DCMIP 2012 T12)',   &
      /,'================================================================================',/,/)

      end subroutine dcmip_tracers12_transport
