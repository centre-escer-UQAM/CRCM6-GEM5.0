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
!***s/r yyg_initblenu - to initialize communication pattern for U field
!
      Subroutine yyg_initblenu()
      use tdpack
      use glb_ld
      use glb_pil
      use ptopo
      use yyg_blnu
      implicit none
#include <arch_specific.hf>
!
!author
!     Abdessamad Qaddouri/V.Lee - October 2009
!  PLEASE consult Abdessamad or Vivian before modifying this routine.
!
!revision
!  v4.8  V.Lee - Correction for limiting range in point found on other grid
!

      integer i,j,kk,ii,jj,ki,ksend,krecv
      integer imx1,imx2
      integer imy1,imy2
      integer xmin,xmax,xmaxu,ymin,ymax,ymaxv
      integer adr
      integer, dimension (:), pointer :: recv_len, send_len
      real*8  xx_8(G_niu,G_nj),yy_8(G_niu,G_nj)
      real*8  xg_8(1-G_ni:2*G_ni),yg_8(1-G_nj:2*G_nj)
      real*8  xgu_8(1-G_ni:2*G_ni-1),ygv_8(1-G_nj:2*G_nj-1)
      real*8  s(2,2),h1,h2
      real*8  x_d,y_d,x_a,y_a
      real*8 TWO_8
      parameter( TWO_8   = 2.0d0 )
!
!     Localise could get point way outside of the actual grid in search
!     So extend all global arrays: xg_8,yg_8, xgu_8,ygv_8

      do i=1,G_ni
         xg_8(i) = G_xg_8(i)
      end do
      do j=1,G_nj
         yg_8(j) = G_yg_8(j)
      enddo

      do i=-G_ni+1,0
         xg_8(i) = xg_8(i+G_ni) - TWO_8*pi_8
      end do
      do i=G_ni+1,2*G_ni
         xg_8(i) = xg_8(i-G_ni) + TWO_8*pi_8
      end do

      yg_8( 0    ) = -(yg_8(1) + pi_8)
      yg_8(-1    ) = -TWO_8*pi_8 -  &
           (yg_8(0)+yg_8(1)+yg_8(2))
      yg_8(G_nj+1) =  pi_8 - yg_8(G_nj)
      yg_8(G_nj+2) =  TWO_8*pi_8 - &
           (yg_8(G_nj+1)+yg_8(G_nj)+yg_8(G_nj-1))
      do j=-2,-G_nj+1,-1
         yg_8(j) = 1.01*yg_8(j+1)
      end do
      do j=G_nj+3,2*G_nj
         yg_8(j) = 1.01*yg_8(j-1)
      end do

      do i=1-G_ni,2*G_ni-1
      xgu_8(i)=0.5D0 *(xg_8(i+1)+xg_8(i))
      enddo
      do j=1-G_nj,2*G_nj-1
      ygv_8(j)= 0.5D0*(yg_8(j+1)+yg_8(j))
      enddo
!
      do j=1,G_nj
      do i=1,G_niu
         xx_8(i,j)=xgu_8(i)
         yy_8(i,j)=yg_8(j)
      enddo
      enddo
      xmin=1-G_ni
      xmax=2*G_ni
      xmaxu=2*G_ni-1
      ymin=1-G_nj
      ymax=2*G_nj
      ymaxv=2*G_nj-1

!Delta xg, yg is not identical between xg(i) and xg(i+1)
!h1, h2 used in this routine is ok as it is a close estimate for
!creating YY pattern exchange and it works on the global tile

      h1=xg_8(2)-xg_8(1)
      h2=yg_8(2)-yg_8(1)
!
! And allocate temp vectors needed for counting for each processor
!
      allocate (recv_len (Ptopo_numproc))
      allocate (send_len (Ptopo_numproc))
      recv_len (:)=0
      send_len (:)=0

!
! FIRST PASS is to find the number of processor to tag for
! communication and the number of items to send and receive for each
! processor before allocating the vectors
!
!
!
      do j=1+glb_pil_s, G_nj-glb_pil_n
      do i=1+glb_pil_w, G_niu-glb_pil_e
!        U vector
         x_d=xx_8(i,j)-acos(-1.D0)
         y_d=yy_8(i,j)
         call smat(s,x_a,y_a,x_d,y_d)
         x_a=x_a+(acos(-1.D0))

         call localise_blend(imx1,imy1,x_a,y_a, &
                          xgu_8,yg_8,xmin,xmaxu,ymin,ymax,h1,h2)
         call localise_blend(imx2,imy2,x_a,y_a, &
                          xg_8,ygv_8,xmin,xmax,ymin,ymaxv,h1,h2)


! check if this point can be found in the other grid
! It is important to do this check before min-max
!   (Imx,Imy )could be zero or negatif or 1<(Imx,Imy )<(G_ni,G_nj)

         if (imx1 > 1+glb_pil_w .and. imx1 < G_niu-glb_pil_e .and. &
             imy1 > 1+glb_pil_s .and. imy1 < G_nj -glb_pil_n  .and. &
             imx2 > 1+glb_pil_w .and. imx2 < G_ni -glb_pil_e .and. &
             imy2 > 1+glb_pil_s .and. imy2 < G_njv-glb_pil_n) then

             imx1 = min(max(imx1-1,glb_pil_w+1),G_niu-glb_pil_e-3)
             imy1 = min(max(imy1-1,glb_pil_s+1),G_nj -glb_pil_n-3)
             imx2 = min(max(imx2-1,glb_pil_w+1),G_ni -glb_pil_e-3)
             imy2 = min(max(imy2-1,glb_pil_s+1),G_njv-glb_pil_n-3)

!

! check to collect from who
             if (i >= l_i0.and.i <= l_i0+l_ni-1 .and. &
                 j >= l_j0.and.j <= l_j0+l_nj-1      ) then
                 do kk=1,Ptopo_numproc
                    if (max(imx1,imx2) >= Ptopo_gindx(1,kk).and. &
                        max(imx1,imx2) <= Ptopo_gindx(2,kk).and. &
                        max(imy1,imy2) >= Ptopo_gindx(3,kk).and. &
                        max(imy1,imy2) <= Ptopo_gindx(4,kk) ) then
                        recv_len(kk)=recv_len(kk)+1
                    endif
                 enddo

             endif

! check to send to who
             if (max(imx1,imx2) >= l_i0.and.         &
                 max(imx1,imx2) <= l_i0+l_ni-1 .and. &
                 max(imy1,imy2) >= l_j0.and.         &
                 max(imy1,imy2) <= l_j0+l_nj-1) then
                 do kk=1,Ptopo_numproc
                    if (i >= Ptopo_gindx(1,kk).and.&
                        i <= Ptopo_gindx(2,kk).and.&
                        j >= Ptopo_gindx(3,kk).and.&
                        j <= Ptopo_gindx(4,kk)     )then
                        send_len(kk)=send_len(kk)+1
                    endif
                 enddo
             endif
         endif
      enddo
      enddo
!
!
!
! Obtain sum of elements to send and receive for each processor
! and the total memory needed to store and receive for each processor
!
     Bln_usend_all=0
     Bln_urecv_all=0
     Bln_usendmaxproc=0
     Bln_urecvmaxproc=0

     do kk=1,Ptopo_numproc
        Bln_usend_all=send_len(kk)+Bln_usend_all
        Bln_urecv_all=recv_len(kk)+Bln_urecv_all

        if (send_len(kk) > 0) Bln_usendmaxproc=Bln_usendmaxproc+1
        if (recv_len(kk) > 0) Bln_urecvmaxproc=Bln_urecvmaxproc+1
     enddo

!
!     print *,'Allocate common vectors'
      allocate (Bln_urecvproc(Bln_urecvmaxproc))
      allocate (Bln_urecv_len(Bln_urecvmaxproc))
      allocate (Bln_urecv_adr(Bln_urecvmaxproc))

      allocate (Bln_usendproc(Bln_usendmaxproc))
      allocate (Bln_usend_len(Bln_usendmaxproc))
      allocate (Bln_usend_adr(Bln_usendmaxproc))
      Bln_urecvproc(:) = 0
      Bln_urecv_len(:) = 0
      Bln_urecv_adr(:) = 0
      Bln_usendproc(:) = 0
      Bln_usend_len(:) = 0
      Bln_usend_adr(:) = 0

!    print*,'Bln_usendmaxproc=',Bln_usendmaxproc,'recvmaxproc=',Bln_urecvmaxproc

     ksend=0
     krecv=0
     Bln_usend_all=0
     Bln_urecv_all=0
!
! Fill the lengths and addresses for selected processors to communicate
!
     do kk=1,Ptopo_numproc
        if (send_len(kk) > 0) then
            ksend=ksend+1
            Bln_usendproc(ksend)=kk
            Bln_usend_len(ksend)=send_len(kk)

            Bln_usend_adr(ksend)= Bln_usend_all
            Bln_usend_all= Bln_usend_all + Bln_usend_len(ksend)
        endif
        if (recv_len(kk) > 0) then
            krecv=krecv+1
            Bln_urecvproc(krecv)=kk
            Bln_urecv_len(krecv)=recv_len(kk)

            Bln_urecv_adr(krecv)= Bln_urecv_all
            Bln_urecv_all= Bln_urecv_all + Bln_urecv_len(krecv)
        endif

     enddo
!    print *,'krecv=',krecv,'Bln_urecvmaxproc=',Bln_urecvmaxproc
!    print *,'ksend=',ksend,'Bln_usendmaxproc=',Bln_usendmaxproc

!     print *,'Bln_urecv_all=',Bln_urecv_all, 'Bln_usend_all=',Bln_usend_all

!
! Now allocate the vectors needed for sending and receiving each processor
!
      if (Bln_urecv_all > 0) then
          allocate (Bln_urecv_i(Bln_urecv_all))
          allocate (Bln_urecv_j(Bln_urecv_all))
          Bln_urecv_i(:) = 0
          Bln_urecv_j(:) = 0
      endif

      if (Bln_usend_all > 0) then
          allocate (Bln_usend_imx1(Bln_usend_all))
          allocate (Bln_usend_imy1(Bln_usend_all))
          allocate (Bln_usend_imx2(Bln_usend_all))
          allocate (Bln_usend_imy2(Bln_usend_all))
          allocate (Bln_usend_xxr(Bln_usend_all))
          allocate (Bln_usend_yyr(Bln_usend_all))
          allocate (Bln_usend_s1(Bln_usend_all))
          allocate (Bln_usend_s2(Bln_usend_all))
          Bln_usend_imx1(:) = 0
          Bln_usend_imy1(:) = 0
          Bln_usend_imx2(:) = 0
          Bln_usend_imy2(:) = 0
          Bln_usend_xxr(:) = 0.0
          Bln_usend_yyr(:) = 0.0
          Bln_usend_s1(:) = 0.0
          Bln_usend_s2(:) = 0.0
      endif
!

      recv_len(:)=0
      send_len(:)=0
!
! SECOND PASS is to initialize the vectors with information for communication
!
!
      do j=1+glb_pil_s, G_nj-glb_pil_n
      do i=1+glb_pil_w, G_niu-glb_pil_e
!        U vector
         x_d=xx_8(i,j)-acos(-1.D0)
         y_d=yy_8(i,j)
         call smat(s,x_a,y_a,x_d,y_d)
         x_a=x_a+(acos(-1.D0))
         call localise_blend(imx1,imy1,x_a,y_a, &
                          xgu_8,yg_8,xmin,xmaxu,ymin,ymax,h1,h2)
         call localise_blend(imx2,imy2,x_a,y_a, &
                          xg_8,ygv_8,xmin,xmax,ymin,ymaxv,h1,h2)

! check if this point can be found in the other grid
! It is important to do this check before min-max
!   (Imx,Imy )could be zero or negatif or 1<(Imx,Imy )<(G_ni,G_nj)

         if (imx1 > 1+glb_pil_w .and. imx1 < G_niu-glb_pil_e .and. &
             imy1 > 1+glb_pil_s .and. imy1 < G_nj -glb_pil_n  .and. &
             imx2 > 1+glb_pil_w .and. imx2 < G_ni -glb_pil_e .and. &
             imy2 > 1+glb_pil_s .and. imy2 < G_njv-glb_pil_n) then

             imx1 = min(max(imx1-1,glb_pil_w+1),G_niu-glb_pil_e-3)
             imy1 = min(max(imy1-1,glb_pil_s+1),G_nj -glb_pil_n-3)
             imx2 = min(max(imx2-1,glb_pil_w+1),G_ni -glb_pil_e-3)
             imy2 = min(max(imy2-1,glb_pil_s+1),G_njv-glb_pil_n-3)

!

! check to collect from who
             if (i >= l_i0.and.i <= l_i0+l_ni-1 .and. &
                 j >= l_j0.and.j <= l_j0+l_nj-1      ) then
                 do kk=1,Bln_urecvmaxproc
                    ki=Bln_urecvproc(kk)
                    if (max(imx1,imx2) >= Ptopo_gindx(1,ki).and. &
                        max(imx1,imx2) <= Ptopo_gindx(2,ki).and. &
                        max(imy1,imy2) >= Ptopo_gindx(3,ki).and. &
                        max(imy1,imy2) <= Ptopo_gindx(4,ki) ) then
                        recv_len(kk)=recv_len(kk)+1
                        adr=Bln_urecv_adr(kk)+recv_len(kk)
                        ii=i-l_i0+1
                        jj=j-l_j0+1
                        Bln_urecv_i(adr)=ii
                        Bln_urecv_j(adr)=jj
                    endif
                 enddo
             endif

! check to send to who
             if (max(imx1,imx2) >= l_i0.and.         &
                 max(imx1,imx2) <= l_i0+l_ni-1 .and. &
                 max(imy1,imy2) >= l_j0.and.         &
                 max(imy1,imy2) <= l_j0+l_nj-1) then
                 do kk=1,Bln_usendmaxproc
                    ki=Bln_usendproc(kk)
                    if (i >= Ptopo_gindx(1,ki).and. &
                        i <= Ptopo_gindx(2,ki).and. &
                        j >= Ptopo_gindx(3,ki).and. &
                        j <= Ptopo_gindx(4,ki)     )then
                        send_len(kk)=send_len(kk)+1
                        adr=Bln_usend_adr(kk)+send_len(kk)
                        Bln_usend_imx1(adr)=imx1-l_i0+1
                        Bln_usend_imy1(adr)=imy1-l_j0+1
                        Bln_usend_imx2(adr)=imx2-l_i0+1
                        Bln_usend_imy2(adr)=imy2-l_j0+1
                        Bln_usend_xxr(adr)=x_a
                        Bln_usend_yyr(adr)=y_a
                        Bln_usend_s1(adr)=s(1,1)
                        Bln_usend_s2(adr)=s(1,2)
                    endif
                 enddo
             endif
         endif
      enddo
      enddo

!
!
!Check receive lengths from each processor
!     do ki=1,Bln_urecvmaxproc
!        kk=Bln_urecvproc(ki)
!        if (Ptopo_couleur == 0) then
!            kkproc = kk+Ptopo_numproc-1
!        else
!            kkproc = kk -1
!        endif
!    write(6,1000) 'Bln_urecvw_len',kkproc,Bln_urecvw_len(kk),Bln_urecvw_adr(kk)
!    write(6,1000) 'Bln_urecve_len',kkproc,Bln_urecve_len(kk),Bln_urecve_adr(kk)
!    write(6,1000) 'Bln_urecvs_len',kkproc,Bln_urecvs_len(kk),Bln_urecvs_adr(kk)
!    write(6,1000) 'Bln_urecvn_len',kkproc,Bln_urecvn_len(kk),Bln_urecvn_adr(kk)
!   enddo

!Check send lengths to each processor

!     do ki=1,Bln_usendmaxproc
!        kk=Bln_usendproc(ki)
!        if (Ptopo_couleur == 0) then
!            kkproc = kk+Ptopo_numproc-1
!        else
!            kkproc = kk -1
!        endif
! write(6,1000) 'Bln_usendw_len',kkproc,Bln_usendw_len(kk),Bln_usendw_adr(kk)
! write(6,1000) 'Bln_usende_len',kkproc,Bln_usende_len(kk),Bln_usende_adr(kk)
! write(6,1000) 'Bln_usends_len',kkproc,Bln_usends_len(kk),Bln_usends_adr(kk)
! write(6,1000) 'Bln_usendn_len',kkproc,Bln_usendn_len(kk),Bln_usendn_adr(kk)
!     enddo
      deallocate (recv_len,send_len)

 1000 format(a15,i3,'=',i5,'bytes, addr=',i5)
 1001 format(a15,i3,'=',i4,'bytes   i:', i3,' j:',i3)

!
      return
      end

