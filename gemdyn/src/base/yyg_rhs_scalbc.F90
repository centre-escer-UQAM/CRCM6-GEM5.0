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
!*s/r yyg_rhs_scalbc - to interpolate and exchange RHS data
!


      Subroutine yyg_rhs_scalbc(rhs_dst1,sol_src1,minx,maxx,miny,maxy, &
                                NK,iter,Linfini)
      use gem_options
      use geomh
      use glb_ld
      use glb_pil
      use sol
      use ptopo
      use yyg_rhs
      implicit none
#include <arch_specific.hf>

!
!author
!     Abdessamad Qaddouri/V.Lee  - October 2009
!

      integer minx,maxx,miny,maxy,NK
      real*8  sol_src1(minx:maxx,miny:maxy,Nk)
      real*8  rhs_dst1(minx:maxx,miny:maxy,Nk)
      real*8  tab_src_8(l_minx:l_maxx,l_miny:l_maxy,Nk)
      real    tab_src(l_minx:l_maxx,l_miny:l_maxy,Nk)
      real    linfini,L1,L2,L3(2),L4(2),infini(Nk)
      integer ierr,i,j,k,kk,kk_proc,m,mm,iter,adr
      real, dimension (:,:), allocatable :: recv_pil,send_pil
!     integer status(MPI_STATUS_SIZE)
!     integer stat(MPI_STATUS_SIZE,Ptopo_numproc)
      integer status
      integer request(Ptopo_numproc*2)
      real*8, dimension (:),allocatable ::  send_Rhsx_8
      integer tag1,tag2,recvlen,sendlen,ireq
      tag2=14
      tag1=13


      tab_src=0.
      linfini=0.0
      do k=1,NK
         infini(k)=0.0
      do j=1,l_nj
      do i=1,l_ni
         tab_src(i,j,k)=real(sol_src1(i,j,k))
      enddo
      enddo
      enddo

      sendlen=0
      recvlen=0
      ireq=0
      do kk=1,Rhsx_sendmaxproc
         sendlen=max(sendlen,Rhsx_send_len(kk))
      enddo
      do kk=1,Rhsx_recvmaxproc
         recvlen=max(recvlen,Rhsx_recv_len(kk))
      enddo

      call rpn_comm_xch_halo(tab_src, l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,NK, &
                  G_halox,G_haloy,G_periodx,G_periody,l_ni,0 )

      tab_src_8(:,:,:)=dble(tab_src(:,:,:))

!     print *,'sendlen=',sendlen,' recvlen=',recvle
      if (sendlen > 0) then
          allocate(send_pil(sendlen*NK,Rhsx_sendmaxproc))
          allocate(send_Rhsx_8(sendlen*NK))

      endif
      if (recvlen > 0) then
          allocate(recv_pil(recvlen*NK,Rhsx_recvmaxproc))
      endif

!
      do 100 kk=1,Rhsx_sendmaxproc
!
!        For each processor (in other colour)
         if (Ptopo_couleur == 0) then
             kk_proc = Rhsx_sendproc(kk)+Ptopo_numproc-1
         else
             kk_proc = Rhsx_sendproc(kk)-1
         endif

!        prepare to send to other colour processor
         if (Rhsx_send_len(kk) > 0) then
!            prepare something to send

             adr=Rhsx_send_adr(kk)+1
             call int_cub_lag8(send_Rhsx_8,tab_src_8,              &
                             Rhsx_send_imx(adr:adr+Rhsx_send_len(kk)), &
                             Rhsx_send_imy(adr:adr+Rhsx_send_len(kk)), &
                             geomh_x_8,geomh_y_8,l_minx,            &
                             l_maxx,l_miny,l_maxy, NK,              &
                             Rhsx_send_xxr(adr:adr+Rhsx_send_len(kk)), &
                             Rhsx_send_yyr(adr:adr+Rhsx_send_len(kk)), &
                             Rhsx_send_len(kk))

             mm=0
             do m=1,Rhsx_send_len(kk)
                do k=1,NK
                   mm=mm+1
                   infini(k)=max(infini(k),abs(real(send_Rhsx_8(mm))-Sol_rhs(mm,1,kk)))
                   Sol_rhs(mm,1,kk)=real(send_Rhsx_8(mm))
                   send_pil(mm,KK)=real(send_Rhsx_8(mm)*Rhsx_send_sten(adr+m-1))
                enddo
             enddo

             linfini=infini(1)
             do k=2,NK
                linfini=max(linfini,infini(k))
             enddo

             ireq = ireq+1
!            print *,'RHS_scalbc','sending',Rhsx_send_len(kk)*NK,' to ',kk_proc
!            call MPI_ISend (send_pil(1,KK),Rhsx_send_len(kk)*NK,MPI_REAL, &
!                                        kk_proc,tag2+Ptopo_world_myproc, &
!                                        MPI_COMM_WORLD,request(ireq),ierr)
             call RPN_COMM_ISend (send_pil(1,KK),Rhsx_send_len(kk)*NK,&
                                  'MPI_REAL', kk_proc,tag2+Ptopo_world_myproc, &
                                  'MULTIGRID',request(ireq),ierr)
         endif
 100 continue
!
!
!        check to receive from other colour processor
!
      do 101 kk=1,Rhsx_recvmaxproc
!        For each processor (in other colour)

         if (Ptopo_couleur == 0) then
             kk_proc = Rhsx_recvproc(kk)+Ptopo_numproc-1
         else
             kk_proc = Rhsx_recvproc(kk)-1
         endif
         if (Rhsx_recv_len(kk) > 0) then
!            detect something to receive

             ireq = ireq+1
!            print *,'RHS_scalbc','receiving ',Rhsx_recv_len(kk)*NK,' from ',kk_proc
!            call MPI_IRecv(recv_pil(1,KK),Rhsx_recv_len(kk)*NK,MPI_REAL, &
!                   kk_proc,tag2+kk_proc,MPI_COMM_WORLD,request(ireq),ierr)
             call RPN_COMM_IRecv(recv_pil(1,KK),Rhsx_recv_len(kk)*NK,'MPI_REAL',&
                    kk_proc,tag2+kk_proc,'MULTIGRID',request(ireq),ierr)
         endif

 101  continue

!Wait for all done sending and receiving
!     call mpi_waitall(ireq,request,stat,ierr)
      call RPN_COMM_waitall_nostat(ireq,request,ierr)

      l2=0
      do kk=1,Rhsx_sendmaxproc
         do i=1,sendlen*NK
            l2 =max(l2,abs(Sol_rhs(i,1,kk)))
         enddo
      enddo
!     call mpi_allreduce(L2,L1,1,MPI_REAL,MPI_MAX,Ptopo_intracomm,ierr)
      call RPN_COMM_allreduce(L2,L1,1,"MPI_REAL","MPI_MAX","grid",ierr)
      L3(2)=L1
      L2   =L1



      if (sendlen > 0) deallocate(send_pil,send_Rhsx_8)

!     Check the precision
!     call mpi_allreduce(Linfini,L1,1,MPI_REAL,MPI_MAX,Ptopo_intracomm,ierr)
      call RPN_COMM_allreduce(Linfini,L1,1,"MPI_REAL","MPI_MAX","grid",ierr)
      L3(1)=L1
      L4=0.
      Linfini=L1

      if (Ptopo_myproc == 0.and.Ptopo_couleur == 0) then
!         call MPI_Send(L3,2,MPI_REAL,1,tag1,Ptopo_intercomm,ierr)
!         call MPI_Recv(L4,2,MPI_REAL,1,tag2,Ptopo_intercomm,status,ierr)
          call RPN_COMM_Send(L3,2,'MPI_REAL',1,tag1,'GRIDPEERS',ierr)
          call RPN_COMM_Recv(L4,2,'MPI_REAL',1,tag2,'GRIDPEERS',status,ierr)
      endif
      if (Ptopo_myproc == 0.and.Ptopo_couleur == 1) then
!         call MPI_Recv(L4,2,MPI_REAL,0,tag1,Ptopo_intercomm,status,ierr)
!         call MPI_Send(L3,2,MPI_REAL,0,tag2,Ptopo_intercomm,ierr)
          call RPN_COMM_Recv(L4,2,'MPI_REAL',0,tag1,'GRIDPEERS',status,ierr)
          call RPN_COMM_Send(L3,2,'MPI_REAL',0,tag2,'GRIDPEERS',ierr)
      endif
      Linfini=max(Linfini,L4(1))
      L2=max(L2,L4(2))
      if (iter > 1) linfini=linfini/L2
!     call mpi_bcast(linfini,1,mpi_real,0,Ptopo_intracomm,ierr)
      call RPN_COMM_bcast(linfini,1,'MPI_REAL',0,'GRID',ierr)
!      print *,'YYG_rhs:iter=',iter,linfini

! Now fill my results if I have received something

      if (recvlen > 0) then

          do 200 kk=1,Rhsx_recvmaxproc
             mm=0
             do m=1,Rhsx_recv_len(kk)
             adr=Rhsx_recv_adr(kk)+m
             do k=1,NK
                mm=mm+1
                rhs_dst1(Rhsx_recv_i(adr),Rhsx_recv_j(adr),k)= &
                rhs_dst1(Rhsx_recv_i(adr),Rhsx_recv_j(adr),k)-recv_pil(mm,KK)
             enddo
          enddo
 200  continue
      deallocate(recv_pil)

      endif
!
!
      return
      end

