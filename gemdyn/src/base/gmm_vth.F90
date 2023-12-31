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
module gmm_vth
   implicit none
   public
   save

!
!______________________________________________________________________
!                                                                      |
!  GMM variables at TIME th (t0-dt/2)                                  |
!______________________________________________________________________|
!                    |                                                 |
! NAME               | DESCRIPTION                                     |
!--------------------|-------------------------------------------------|
! zdth               | Zdot: generalized vertical velocity             |
!--------------------|-------------------------------------------------|
!  xth               | upstream x position                             |
!  yth               | upstream y position                             |
!  zth               | upstream z position                             |
!----------------------------------------------------------------------
!

      real, pointer, dimension (:,:,:) :: zdth  => null()
      real, pointer, dimension (:    ) :: xth   => null()
      real, pointer, dimension (:    ) :: yth   => null()
      real, pointer, dimension (:    ) :: zth   => null()

      integer, parameter :: MAXNAMELENGTH = 32

      character(len=MAXNAMELENGTH) ::  gmmk_zdth_s, &
                          gmmk_xth_s , gmmk_yth_s , gmmk_zth_s

end module gmm_vth
