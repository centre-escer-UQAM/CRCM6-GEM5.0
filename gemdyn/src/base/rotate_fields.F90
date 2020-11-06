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
!!! rotate_vectors - rotates vector components

module rotate_fields

   use tdpack
   use glb_ld
   use grid_options
   use geomh

   implicit none

   private

   public rotate_vectors, rotate_vectors_pset, rotate_vectors_punset

!author Bernard Dugas - RPN (Jan->Nov 2008)

!revision
! v3.3.2  - Dugas B.    - Code inspired from ZONGINI for local tiles
! v4.812  - Dugas B. (ESCER/UQAM) - Adaptation to GEM 4.8. This means...
!                      - Convertion to a F90 module
!                      - Some re-organisation to minimize operations
!                      - Use the Geomg_(x/y) pair rather than the no-longer
!                        available Geomn_(lon/lat)gs one and change the
!                        correponsding call to LLACAR
!         - Dugas B. (ESCER/UQAM) - Comment the OpenMP directives
!                      - Save the polar rotation angle 'rot' and use it to
!                        determine wether (UU,VV) needs to be re-calculated
!         - Winger K. (ESCER/UQAM) - Correction to first call to llacar_8
!         - Dugas B.  (ESCER/UQAM) - Modify the scope of lonngs/latgs 
! v5.0.0  - Dugas B.  (ESCER/UQAM) - Adaptation to GEM 5.0. This essentially
!                        means converting the cdk includes to 'using' the
!                        corresponding F90 modules  
!
!object
!      This package calculates the local rotation angles between
!      a rotated frame of reference and the normal geographical one.
!      This set of angles can be used to convert model wind vectors
!      to ordinary zonal and meridional winds.

!      Thus, the (UM,VM) model winds vectors convert to (UG,VG)
!      geographical winds via the following transformation:

!      /                         \ /    \   /    \
!      |  cos( phi )  sin( phi ) | | UM | = | UG |
!      | -sin( phi )  cos( phi ) | | VM |   | VG |
!      \                         / \    /   \    /

   character(len=3), save :: initialize_S = 'yes'
   integer,  save :: l_ni0=-1,l_nj0=-1
   real(8),  save :: rot=0.0_8

   logical,  save :: pset_L = .false.
   integer,  save :: pni,pnj,poffi,poffj

   real(8),  dimension (:,:), pointer,save :: cosphi,sinphi

   contains
!  ---------------------------------------------------------------
   subroutine rotate_vectors( uu,vv,lminx,lmaxx,lminy,lmaxy,F_nk, DO_IT )

      implicit none

      logical :: DO_IT
      integer :: lminx,lmaxx,lminy,lmaxy,F_nk
      real    :: uu(lminx:lmaxx,lminy:lmaxy,F_nk)
      real    :: vv(lminx:lmaxx,lminy:lmaxy,F_nk)

!arguments
! uu,vv      I/O    - input (uz,vz) and output (ug,vg)

!implicits
!#include "dcst.cdk"
!#include "glb_ld.cdk"
!#include "grd.cdk"
!#include "geomg.cdk"

!!
!     Work space for coordinate calculations

      real(8),  parameter :: polonlat=90.0_8

      real(8),  dimension (:), allocatable :: lon, lat, longs, latgs
      real(8),  dimension (:), allocatable :: lonr, latr
      real(8),  dimension (:), allocatable :: cart, carot

      real(8)   ri_8(3,3),costeta,sinteta
      real(8)   pcart(3),pcarot(3)
      real(8)   polon,polat

      integer   i,j,k,l, ij
      real(8)   mp(3,2),mt(2,3)
      real(8)   pis2_8,rad2deg_8, rmp, hold

!     ---------------------------------------------------------------

      pis2_8 = Pi_8/2.0 ; rad2deg_8 = 90.0/pis2_8

      INITIAL_PHI_CALCULATION : &
      if (initialize_S == 'yes' .or. &
     .not. (l_ni0 == l_ni .and. l_nj0 == l_nj) ) then

!!$omp critical (PHI_CALCULATION)

       ! Ce 'GOTO 100' est REQUIS en mode OMP. !
       !   Autrement, il n'est jamais active   !
         if (initialize_S == 'no' .and. &
             l_ni0 == l_ni .and. l_nj0 == l_nj) goto 100

         l_ni0 = l_ni ; l_nj0 = l_nj

       ! print *,"lminx,lmaxx,  lminy,lmaxy ",lminx,lmaxx,lminy,lmaxy
       ! print *,"l_minx,l_maxx,l_i0,l_ni,G_ni ",l_miny,l_maxy,l_j0,l_nj,G_ni
       ! print *,"l_miny,l_maxy,l_j0,l_nj,G_nj ",l_minx,l_maxx,l_i0,l_ni,G_nj

         allocate( cosphi(l_ni,l_nj),sinphi(l_ni,l_nj) )
         cosphi = 1.0 ; sinphi = 0.0

!     ---------------------------------------------------------------

         GRD_ROULE_BLOCK : if (Grd_roule) then

!     Calcul de l'angle entre les poles des deux grilles

            do i=1,3
               do j=1,3
                  ri_8(i,j) = Grd_rot_8(j,i)
               end do
            end do

            call llacar_8( pcart, polonlat,polonlat, 1,1 )
            call mxma8( ri_8,1,3,pcart,1,3,pcarot,1,3, 3,3,1)
            call cartall_8( polon, polat, pcarot, 1 )

            rot  = 90. - polat

         endif GRD_ROULE_BLOCK

         ROTATED_GEOMETRY : if ( abs( rot )-0.001 > 0._8 ) then
            
!     Allouer les variables temporaire

            allocate( lon  (l_ni*l_nj)  , lat  (l_ni*l_nj)   )
            allocate( lonr (l_ni*l_nj)  , latr (l_ni*l_nj)   )
            allocate( longs(lminx:lmaxx), latgs(lminy:lmaxy) )
            allocate( cart (3*l_ni*l_nj), carot(3*l_ni*l_nj) )

            longs(1:l_ni) = geomh_x_8(1:l_ni) * rad2deg_8
            latgs(1:l_nj) = geomh_y_8(1:l_nj) * rad2deg_8

!     Calcul des latitudes et longitudes de la
!     grille tournee dans le cadre non-tourne

            call llacar_8( cart, longs(1),latgs(1), l_ni,l_nj )

            call mxma8( ri_8,1,3,cart,1,3,carot,1,3, 3,3, l_ni*l_nj )
            call cartall_8( lon, lat, carot, l_ni*l_nj)

            do i=1,l_ni*l_nj
               lon(i)  = mod(lon(i) + 360.0,360.0)
               lonr(i) = lon(i)*Pi_8/180.
               latr(i) = lat(i)*Pi_8/180.
            end do

            do j=1,l_nj 
            do i=1,l_ni

               ij = (j-1)*l_ni+i

!     Calcul de la matrice de rotation des vents
!     ------------------------------------------

!     Definir les composantes requises de M' 
!     [ ou M':(u,v)geo --> (dx/dt,dy/dt,dz/dt)geo ]
               mp(1,1) = -sin( lonr(ij) )
               mp(2,1) =  cos( lonr(ij) )
               mp(3,1) =  0.0
!!!            mp(1,2) = -sin( latr(ij) )*cos( lonr(ij) )
!!!            mp(2,2) = -sin( latr(ij) )*sin( lonr(ij) )
!!!            mp(3,2) =  cos( latr(ij) )

!     Definir les composantes de MT, la transposee de M
!     [ ou M:(u,v)mod --> (dx/dt,dy/dt,dz/dt)mod ]
               mt(1,1) = -Geomh_sx_8(i)
               mt(1,2) =  Geomh_cx_8(i)
               mt(1,3) =  0.0
               mt(2,1) = -Geomh_sy_8(j)*Geomh_cx_8(i)
               mt(2,2) = -Geomh_sy_8(j)*Geomh_sx_8(i)
               mt(2,3) =  Geomh_cy_8(j)
!
!     Calculer la premiere colonne du produit MT RT M' = TT
!     [ ou R:(repere modele) --> (repere geographique) ] 
               sinteta = 0.0
               costeta = 0.0
!
!     On ne calcule donc que -TT(1,1) (= sin(theta)) et
!     TT(2,1) (= cos(theta))
!
               do k=1,3
                  rmp     = 0.0
                  do l=1,3
                     rmp  = rmp+Grd_rot_8(k,l)*mp(l,1)
                  enddo
                  sinteta = sinteta - mt(1,k)*rmp
                  costeta = costeta + mt(2,k)*rmp
               enddo

!     Angle de rotation phi = theta + PI/2

               cosphi(i,j) = -sinteta
               sinphi(i,j) =  costeta

!     Trouver theta a partir de sin(theta) et cos(theta)
!!!            if ( costeta .ne. 0.0 ) then
!!!               theta(i,j) = atan( sinteta/costeta )
!!!            else if ( sinteta .gt. 0.0 ) then
!!!               theta(i,j) = pis2
!!!            else if ( sinteta .lt. 0.0 ) then
!!!               theta(i,j) = -pis2
!!!            endif

!     theta est defini dans l'interval [ -pi , +pi ]
!!!            if ( costeta .lt. 0.0 ) then
!!!               if ( sinteta .ge. 0.0 ) then
!!!                  theta(i,j) = theta(i,j) + Pi_8
!!!               else
!!!                  theta(i,j) = theta(i,j) - Pi_8
!!!               endif
!!!            endif

!!!            phi(i,j) = theta(i,j) + pis2

            enddo
            enddo

!     Desallouer les variables temporaire

            if (Grd_roule) then
               deallocate( lon  , lat   )
               deallocate( lonr , latr  )
               deallocate( longs, latgs )
               deallocate( cart , carot )
            endif

         endif ROTATED_GEOMETRY

         initialize_S = 'no'

  100    continue

!!$omp end critical (PHI_CALCULATION)

      endif INITIAL_PHI_CALCULATION

      if ( DO_IT .and. ( abs( rot )-0.001 > 0._8 ) ) then

         if (.not.pset_L) then
            poffi = 0 ; poffj = 0 ; pni = l_ni ; pnj = l_nj
         endif

         do k=1,F_nk
            do j=1+poffj,pnj+poffj
               do i=1+poffi,pni+poffi
                  hold      =  cosphi(i,j)*uu(i,j,k)+sinphi(i,j)*vv(i,j,k)
                  vv(i,j,k) = -sinphi(i,j)*uu(i,j,k)+cosphi(i,j)*vv(i,j,k)
                  uu(i,j,k) =  hold
               enddo
            enddo
         enddo

      endif

      return
   end subroutine rotate_vectors
!  ---------------------------------------------------------------
   subroutine rotate_vectors_pset( F_ni0,F_nj0,F_offi,F_offj )

         ! Define offsets for physic's grid

         integer :: F_ni0,F_nj0,F_offi,F_offj

         pset_L = .true.
         pni = F_ni0 ; pnj = F_nj0 ; poffi = F_offi ; poffj = F_offj

      return

   end subroutine rotate_vectors_pset
!  ---------------------------------------------------------------
   subroutine rotate_vectors_punset( )

         ! Un-Define phusic grid's offsets

         pset_L = .false.

      return
   end subroutine rotate_vectors_punset
!  ---------------------------------------------------------------
end module rotate_fields
