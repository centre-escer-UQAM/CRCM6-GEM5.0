!-------------------------------------- LICENCE BEGIN ------------------------
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
!-------------------------------------- LICENCE END --------------------------
!/@*
SUBROUTINE hines_print (flux_u, flux_v, drag_u, drag_v, alt, sigma_t,    &
       &                  sigma_alpha, v_alpha, m_alpha,                   &
       &                  iu_print, iv_print,                              &
       &                  ilprt1, ilprt2, levprt1, levprt2, naz,           &
       &                  nlons, nlevs, nazmth)
    implicit none
#include <arch_specific.hf>
#include <rmnlib_basics.hf>

    INTEGER  naz, ilprt1, ilprt2, levprt1, levprt2
    INTEGER  nlons, nlevs, nazmth
    INTEGER  iu_print, iv_print
    REAL*8  flux_u(nlons,nlevs), flux_v(nlons,nlevs)
    REAL*8  drag_u(nlons,nlevs), drag_v(nlons,nlevs)
    REAL*8  alt(nlons,nlevs), sigma_t(nlons,nlevs)
    REAL*8  sigma_alpha(nlons,nlevs,nazmth)
    REAL*8  v_alpha(nlons,nlevs,nazmth), m_alpha(nlons,nlevs,nazmth)
!@Author aug. 8/95 - c. mclandress
!@Object
!  Print out altitude profiles of various quantities from
!  hines' doppler spread gravity wave drag parameterization scheme.
!  (note: only for naz = 4 or 8). 
!@Arguments
!            - Input -
! iu_print   1 to print out values in east-west direction.
! iv_print   1 to print out values in north-south direction.
! ilprt1     first longitudinal index to print.
! ilprt2     last longitudinal index to print.
! levprt1    first altitude level to print.
! levprt2    last altitude level to print.
!*@/
    !
    !  internal variables.
    !
    INTEGER  n_east, n_west, n_north, n_south
    INTEGER  i, l
    !-----------------------------------------------------------------------
    !
    !  azimuthal indices of cardinal directions.
    !
    n_east = 1
    IF (naz.EQ.4)  THEN
       n_west  = 3       
       n_north = 2
       n_south = 4       
    ELSE IF (naz.EQ.8)  THEN
       n_west  = 5       
       n_north = 3
       n_south = 7       
    END IF
    !
    !  print out values for range of longitudes.
    !
    DO i = ilprt1,ilprt2
       !
       !  print east-west wind, sigmas, cutoff wavenumbers, flux and drag.
       !
       IF (iu_print.EQ.1)  THEN
          WRITE (RMN_STDOUT,*) 
          WRITE (RMN_STDOUT,'(a,i3)') 'hines gw (east-west) at longitude i =',i
          WRITE (RMN_STDOUT,6005) 
6005      FORMAT (15x,' u ',2x,'sig_e',2x,'sig_t',3x,'m_e',  &
               &            4x,'m_w',4x,'fluxu',5x,'gwdu')
          DO l = levprt1,levprt2
             WRITE (RMN_STDOUT,6701) alt(i,l)/1.e3, v_alpha(i,l,n_east),   &
                  &                          sigma_alpha(i,l,n_east), sigma_t(i,l),  &
                  &                          m_alpha(i,l,n_east)*1.e3,   &
                  &                          m_alpha(i,l,n_west)*1.e3,  &
                  &                          flux_u(i,l)*1.e5, drag_u(i,l)*24.*3600.
          END DO
6701      FORMAT (' z=',f7.2,1x,3f7.1,2f7.3,f9.4,f9.3)
       END IF
       !
       !  print north-south winds, sigmas, cutoff wavenumbers, flux and drag.
       !
       IF (iv_print.EQ.1)  THEN
          WRITE(RMN_STDOUT,*) 
          WRITE(RMN_STDOUT,'(a,i3)') 'hines gw (north-south) at longitude i =',i
          WRITE(RMN_STDOUT,6006) 
6006      FORMAT (15x,' v ',2x,'sig_n',2x,'sig_t',3x,'m_n',   &
               &            4x,'m_s',4x,'fluxv',5x,'gwdv')
          DO l = levprt1,levprt2
             WRITE (RMN_STDOUT,6701) alt(i,l)/1.e3, v_alpha(i,l,n_north),    &
                  &                          sigma_alpha(i,l,n_north), sigma_t(i,l),   &
                  &                          m_alpha(i,l,n_north)*1.e3,    &
                  &                          m_alpha(i,l,n_south)*1.e3,   &
                  &                          flux_v(i,l)*1.e5, drag_v(i,l)*24.*3600.
          END DO
       END IF
       !
    END DO
    !
    RETURN
    !-----------------------------------------------------------------------
  END SUBROUTINE hines_print
