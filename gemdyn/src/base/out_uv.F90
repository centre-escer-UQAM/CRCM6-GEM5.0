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

!**s/r out_uv - output winds

      subroutine out_uv (levset, set)
      use vertical_interpolation, only: vertint2
      use vGrid_Descriptors, only: vgrid_descriptor,vgd_get,vgd_free,VGD_OK,VGD_ERROR
      use vgrid_wb, only: vgrid_wb_get
      use gmm_vt1
      use gmm_pw
      use gem_options
      use glb_ld
      use out_mod
      use out3
      use levels
      use outp
      use outd
      use gmm_itf_mod
      use rotate_fields, only: rotate_vectors
      implicit none
#include <arch_specific.hf>

      integer levset, set


      type(vgrid_descriptor) :: vcoord
      logical :: write_diag_lev,near_sfc_L
      integer ii,i,j,k,istat,kind,nko
      integer i0,in,j0,jn,pnuu,pnvv,pnuv,psum,pnu,pnv
      integer, dimension(:), allocatable :: indo
      integer, dimension(:), pointer     :: ip1m
      real uu(l_minx:l_maxx,l_miny:l_maxy,G_nk+1),&
           vv(l_minx:l_maxx,l_miny:l_maxy,G_nk+1)
      real, dimension(:    ), allocatable::prprlvl,rf
      real, dimension(:    ), pointer, save :: hybm  => null()
      real, dimension(:,:  ), pointer :: udiag,vdiag => null()
      real, dimension(:,:,:), allocatable:: uv_pres,uu_pres,vv_pres,cible
      real, dimension(:,:,:), allocatable:: w1,w2
      real hybm_gnk2(1)
      integer ind0(1)
!
!-------------------------------------------------------------------
!
      pnuu=0 ; pnvv=0 ; pnuv=0 ; pnu=0 ; pnv=0

      do ii=1,Outd_var_max(set)
         if (Outd_var_S(ii,set) == 'UU')then
            pnuu=ii
         endif
         if (Outd_var_S(ii,set) == 'VV')then
            pnvv=ii
         endif
         if (Outd_var_S(ii,set) == 'UV')then
            pnuv=ii
         endif
         if (Outd_var_S(ii,set) == 'U' )then
            pnu =ii
         endif
         if (Outd_var_S(ii,set) == 'V' )then
            pnv =ii
         endif
      enddo

      psum=pnuu+pnuv+pnvv+pnu+pnv
      if (psum == 0)return

      i0 = 1 ; in = l_ni
      j0 = 1 ; jn = l_nj

      nullify (pw_uu_plus, pw_vv_plus, udiag, vdiag)
      istat = gmm_get(gmmk_pw_uu_plus_s, pw_uu_plus)
      istat = gmm_get(gmmk_pw_vv_plus_s, pw_vv_plus)
      istat = gmm_get(gmmk_diag_uu_s   , udiag     )
      istat = gmm_get(gmmk_diag_vv_s   , vdiag     )

      if (pnu+pnv > 0) then
         allocate ( w1(l_minx:l_maxx,l_miny:l_maxy,G_nk+1) )
         allocate ( w2(l_minx:l_maxx,l_miny:l_maxy,G_nk+1) )
      endif

      if (Level_typ_S(levset) == 'M') then  ! Output on model levels

         kind=Level_kind_ip1
!        Setup the indexing for output
         allocate (indo(G_nk+1))
         call out_slev2 ( Level(1,levset), Level_max(levset), &
                          Level_momentum,indo,nko,near_sfc_L)

         write_diag_lev= near_sfc_L .and. out3_sfcdiag_L

!        Retreieve vertical coordinate description
         if ( .not. associated(hybm) ) then
            nullify(ip1m,hybm)
            istat = vgrid_wb_get('ref-m',vcoord,ip1m)
            deallocate(ip1m); nullify(ip1m)
            if (vgd_get(vcoord,'VCDM - vertical coordinate (m)',hybm) /= VGD_OK) istat = VGD_ERROR
            istat = vgd_free(vcoord)
         endif
         hybm_gnk2(1)=hybm(G_nk+2)
         ind0(1)=1

         if ( (pnuu /= 0) .or. (pnvv /= 0) ) then
            call out_fstecr3(pw_uu_plus,l_minx,l_maxx,l_miny,l_maxy,hybm,&
                   'UU  ',Outd_convmult(pnuu,set),Outd_convadd(pnuu,set),&
                   kind,-1,G_nk, indo, nko,Outd_nbit(pnuu,set),.false. )
            call out_fstecr3(pw_vv_plus,l_minx,l_maxx,l_miny,l_maxy,hybm,&
                   'VV  ',Outd_convmult(pnvv,set),Outd_convadd(pnvv,set),&
                   kind,-1,G_nk, indo, nko,Outd_nbit(pnvv,set),.false. )
            if (write_diag_lev) then
               call out_fstecr3(udiag, l_minx,l_maxx,l_miny,l_maxy,&
                                hybm_gnk2, 'UU', Outd_convmult(pnuu,set),&
                                Outd_convadd(pnuu,set),Level_kind_diag,-1,1,&
                                ind0,1,Outd_nbit(pnuu,set),.false. )
               call out_fstecr3(vdiag, l_minx,l_maxx,l_miny,l_maxy,&
                                hybm_gnk2, 'VV', Outd_convmult(pnuu,set),&
                                Outd_convadd(pnuu,set),Level_kind_diag,-1,1,&
                                ind0,1,Outd_nbit(pnuu,set),.false. )
            endif
         endif

         if (pnuv /= 0) then
            do k = 1, G_nk
               do j = j0, jn
               do i = i0, in
                  uu(i,j,k) = sqrt(pw_uu_plus(i,j,k)*pw_uu_plus(i,j,k)+ &
                                   pw_vv_plus(i,j,k)*pw_vv_plus(i,j,k))
               enddo
               enddo
            enddo
            call out_fstecr3(uu,l_minx,l_maxx,l_miny,l_maxy,hybm      ,&
                 'UV  ',Outd_convmult(pnuv,set),Outd_convadd(pnuv,set),&
                 kind,-1,G_nk, indo, nko, Outd_nbit(pnuv,set),.false. )
            if (write_diag_lev) then
               do j = j0, jn
               do i = i0, in
                  uu(i,j,1) = sqrt(udiag(i,j)*udiag(i,j)+ &
                                   vdiag(i,j)*vdiag(i,j))
               enddo
               enddo
               call out_fstecr3(uu,l_minx,l_maxx,l_miny,l_maxy, hybm_gnk2,&
                       'UV  ',Outd_convmult(pnuv,set),Outd_convadd(pnuv,set),&
                       Level_kind_diag,-1,1,ind0,1,Outd_nbit(pnuv,set),.false. )
            endif
         endif

         if (pnu+pnv > 0) then

            w1(:,:,1:G_nk+1) = uu(:,:,1:G_nk+1)
            w2(:,:,1:G_nk+1) = vv(:,:,1:G_nk+1)

            call rotate_vectors( w1,w2, l_minx,l_maxx,l_miny,l_maxy,G_nk+1, .true. )

            if (pnu > 0) then
               call out_fstecr3(w1,l_minx,l_maxx,l_miny,l_maxy,hybm,      &
                   'U   ',Outd_convmult(pnu,set),Outd_convadd(pnu,set),   &
                   kind,-1,G_nk, indo, nko,Outd_nbit(pnu,set),.false. )
               if (write_diag_lev)                                        &
                  call out_fstecr3(w1(l_minx,l_miny,G_nk+1),              &
                     l_minx,l_maxx,l_miny,l_maxy, hybm_gnk2,              &
                    'U   ',Outd_convmult(pnu,set),Outd_convadd(pnu,set),  &
                    Level_kind_diag,-1,1,ind0,1, Outd_nbit(pnu,set),.false. )
            endif

            if (pnv > 0) then
               call out_fstecr3(w2,l_minx,l_maxx,l_miny,l_maxy,hybm,      &
                  'V   ',Outd_convmult(pnv,set),Outd_convadd(pnv,set),    &
                  kind,-1,G_nk, indo, nko,Outd_nbit(pnv,set),.false. )
               if (write_diag_lev)                                        &
                  call out_fstecr3(w2(l_minx,l_miny,G_nk+1),              &
                     l_minx,l_maxx,l_miny,l_maxy, hybm(G_nk+2),           &
                     'V   ',Outd_convmult(pnv,set),Outd_convadd(pnv,set), &
                     Level_kind_diag,-1,1,ind0,1, Outd_nbit(pnv,set),.false. )
            endif

         endif

         deallocate(indo)

      else   ! Output on pressure levels

         nullify (pw_log_pm)
         istat= gmm_get(gmmk_pw_log_pm_s, pw_log_pm)

!        Set kind to 2 for pressure output
         kind=2
!        Setup the indexing for output
         nko=Level_max(levset)
         allocate ( indo(nko), rf(nko) , prprlvl(nko), &
                    cible(l_minx:l_maxx,l_miny:l_maxy,nko) )
         do i = 1, nko
            indo(i)=i
            rf(i)= Level(i,levset)
            prprlvl(i) = rf(i) * 100.0
            cible(:,:,i) = log(prprlvl(i))
         enddo

         allocate(uu_pres(l_minx:l_maxx,l_miny:l_maxy,nko   ))
         allocate(vv_pres(l_minx:l_maxx,l_miny:l_maxy,nko   ))
         uu(:,:,1:G_nk) = pw_uu_plus(:,:,1:G_nk)
         vv(:,:,1:G_nk) = pw_vv_plus(:,:,1:G_nk)
         if (out3_sfcdiag_L) then
            uu(:,:,G_nk+1) = udiag
            vv(:,:,G_nk+1) = vdiag
         else
            uu(:,:,G_nk+1) = pw_uu_plus(:,:,G_nk)
            vv(:,:,G_nk+1) = pw_vv_plus(:,:,G_nk)
         endif

         if (nko > G_nk+1 .and. pnu+pnv > 0)then
            deallocate( w1,w2 )
            allocate ( w1(l_minx:l_maxx,l_miny:l_maxy,nko) )
            allocate ( w2(l_minx:l_maxx,l_miny:l_maxy,nko) )
         endif

!        Vertical interpolation

         call vertint2 ( uu_pres, cible, nko, uu, pw_log_pm, G_nk+1,&
                         l_minx,l_maxx,l_miny,l_maxy, 1,l_ni,1,l_nj,&
                         inttype=Out3_vinterp_type_S )
         call vertint2 ( vv_pres, cible, nko, vv, pw_log_pm, G_nk+1,&
                         l_minx,l_maxx,l_miny,l_maxy, 1,l_ni,1,l_nj,&
                         inttype=Out3_vinterp_type_S )

         if (pnuv /= 0) then
            allocate(uv_pres(l_minx:l_maxx,l_miny:l_maxy,nko   ))
            do k =  1, nko
            do j = j0, jn
            do i = i0, in
              uv_pres(i,j,k) = sqrt(uu_pres(i,j,k)*uu_pres(i,j,k)+ &
                                    vv_pres(i,j,k)*vv_pres(i,j,k))
            enddo
            enddo
            enddo
            if (Outd_filtpass(pnuv,set) > 0) &
               call filter2( uv_pres,Outd_filtpass(pnuv,set),&
                             Outd_filtcoef(pnuv,set), &
                             l_minx,l_maxx,l_miny,l_maxy,nko)
            call out_fstecr3(uv_pres,l_minx,l_maxx,l_miny,l_maxy,rf   ,&
                 'UV  ',Outd_convmult(pnuv,set),Outd_convadd(pnuv,set),&
                 kind,-1,nko, indo, nko, Outd_nbit(pnuv,set),.false. )
            deallocate (uv_pres)
         endif

         if ( (pnuu /= 0) .or. (pnvv /= 0) )then
            if (Outd_filtpass(pnuu,set) > 0) &
                call filter2( uu_pres,Outd_filtpass(pnuu,set),&
                              Outd_filtcoef(pnuu,set), &
                              l_minx,l_maxx,l_miny,l_maxy,nko)
            call out_fstecr3(uu_pres,l_minx,l_maxx,l_miny,l_maxy,rf   ,&
                 'UU  ',Outd_convmult(pnuu,set),Outd_convadd(pnuu,set),&
                 kind,-1,nko, indo, nko, Outd_nbit(pnuu,set),.false. )
            if (Outd_filtpass(pnvv,set) > 0) &
                 call filter2( vv_pres,Outd_filtpass(pnvv,set),&
                               Outd_filtcoef(pnvv,set), &
                               l_minx,l_maxx,l_miny,l_maxy,nko)
            call out_fstecr3(vv_pres,l_minx,l_maxx,l_miny,l_maxy,rf   ,&
                 'VV  ',Outd_convmult(pnvv,set),Outd_convadd(pnvv,set),&
                 kind,-1,nko, indo, nko, Outd_nbit(pnvv,set),.false. )
         endif

         if (pnu+pnv > 0) then

            w1(:,:,1:nko) = uu_pres(:,:,1:nko)
            w2(:,:,1:nko) = vv_pres(:,:,1:nko)

            call rotate_vectors( w1,w2, l_minx,l_maxx,l_miny,l_maxy,nko, .true. )

            if (pnu > 0) then
               if (Outd_filtpass(pnu,set) > 0)                           &
                  call filter2( w1,Outd_filtpass(pnu,set),               &
                                Outd_filtcoef(pnu,set),                  &
                                l_minx,l_maxx,l_miny,l_maxy,nko )
                  call out_fstecr3( w1,l_minx,l_maxx,l_miny,l_maxy,rf,   &
                    'U   ',Outd_convmult(pnu,set),Outd_convadd(pnu,set), &
                    kind,-1,nko, indo, nko, Outd_nbit(pnu,set),.false. )
            endif

            if (pnv > 0) then
               if (Outd_filtpass(pnv,set) > 0)                           &
                  call filter2( w2,Outd_filtpass(pnv,set),               &
                               Outd_filtcoef(pnv,set),                   &
                               l_minx,l_maxx,l_miny,l_maxy,nko )
                  call out_fstecr3( w2,l_minx,l_maxx,l_miny,l_maxy,rf,   &
                    'V   ',Outd_convmult(pnv,set),Outd_convadd(pnv,set), &
                    kind,-1,nko, indo, nko, Outd_nbit(pnv,set),.false. )
            endif

         endif

         deallocate(indo,rf,prprlvl,uu_pres,vv_pres,cible)
         if (pnu+pnv > 0) deallocate( w1,w2 )

      endif
!
!-------------------------------------------------------------------
!
      return
      end subroutine out_uv
