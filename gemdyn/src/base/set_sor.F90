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

!**s/r set_sor - initialization of all control parameters for output

      subroutine set_sor ()
      use iso_c_binding
      use timestr_mod
      implicit none
#include <arch_specific.hf>

!author
!     V. Lee    - rpn - July  2004 (from dynout2 v3_12)
!
!revision
! v4_80 - Desgagne M.       - major re-factorization of output

#include <WhiteBoard.hf>
#include <clib_interface_mu.hf>
#include "lun.cdk"
#include "ptopo.cdk"
#include "out.cdk"
#include "out3.cdk"
#include "grd.cdk"
#include "level.cdk"
#include "timestep.cdk"
#include "schm.cdk"
#include "step.cdk"
#include "lctl.cdk"
#include "outd.cdk"
#include "outp.cdk"
#include "outc.cdk"
#include "hgc.cdk"
#include "cstv.cdk"
#include "path.cdk"
      include "rpn_comm.inc"

      integer,external :: srequet

      integer,parameter :: NBUS = 3
      character(len=9) :: BUS_LIST_S(NBUS) = &
                  (/'Entry    ', 'Permanent', 'Volatile '/)

      character*5 blank_S
      character*8 unit_S
      character*256 fn
      character(len=32) :: varname_S,outname_S,bus0_S
      logical iela
      integer pnerror,i,idx,k,j,levset,ierr,cnt
      integer, dimension (4,2) :: ixg, ixgall
      integer istat,options,rsti
      real, dimension (4,2) :: rot, rotall
!
!-------------------------------------------------------------------
!
      if (Lun_out.gt.0) write(Lun_out,5200)

      pnerror    = 0

      fn = trim(Path_outcfg_S)
      inquire (FILE=fn,EXIST=iela)
      if (iela) then
         pnerror = pnerror +  srequet()
      else
         pnerror = pnerror + 1
      endif
      
      if  (pnerror .gt. 0) then
         if (Lun_out.gt.0) write(Lun_out,5000) pnerror
      endif
      
      Out3_accavg_L = .false.
      Out3_ndigits = max(3,Out3_ndigits)
      select case(Out3_unit_s(3:3))
      case('STE') 
         unit_S = 'TIMESTEPS'
         Out3_ndigits = max(6,Out3_ndigits)
      case('SEC')
         unit_S = 'SECONDS'
         Out3_ndigits = max(6,Out3_ndigits)
      case('MIN')
         unit_S = 'MINUTES'
      case('HOU')
         unit_S = 'HOURS'
      case('DAY')
         unit_S = 'DAYS'
      case('MON')
         unit_S = 'MONTHS'
      case default
         unit_S = 'HOURS'
      end select
      if (Lun_out.gt.0) write(Lun_out,3000)unit_S

      ! Transfer filter and xnbit info to requested variables
      do k=1, Outd_sets
         do j=1,Outd_var_max(k)
            do i=1,Out3_filtpass_max
               if (Outd_varnm_S(j,k) .eq. Out3_filt_S(i)) then
                  Outd_filtpass(j,k) = Out3_filtpass(i)
                  Outd_filtcoef(j,k) = Out3_filtcoef(i)
               endif
            enddo
            do i=1,Out3_xnbits_max
               if (Outd_varnm_S(j,k) .eq. Out3_xnbits_S(i)) then
                  Outd_nbit(j,k) = Out3_xnbits(i)
               endif
            enddo
         enddo
      enddo
      do k=1, Outp_sets
         do j=1,Outp_var_max(k)
            do i=1,Out3_filtpass_max
               if (Outp_varnm_S(j,k) .eq. Out3_filt_S(i)) then
                  Outp_filtpass(j,k) = Out3_filtpass(i)
                  Outp_filtcoef(j,k) = Out3_filtcoef(i)
               endif
            enddo
            do i=1,Out3_xnbits_max
               if (Outp_varnm_S(j,k) .eq. Out3_xnbits_S(i)) then
                  Outp_nbit(j,k) = Out3_xnbits(i)
               endif
            enddo
         enddo
         if (Outp_avg_L(k).or.Outp_accum_L(k)) Out3_accavg_L = .true.
      enddo
      do k=1, Outc_sets
         do j=1,Outc_var_max(k)
            do i=1,Out3_filtpass_max
               if (Outc_varnm_S(j,k) .eq. Out3_filt_S(i)) then
                  Outc_filtpass(j,k) = Out3_filtpass(i)
                  Outc_filtcoef(j,k) = Out3_filtcoef(i)
               endif
            enddo
            do i=1,Out3_xnbits_max
               if (Outc_varnm_S(j,k) .eq. Out3_xnbits_S(i)) then
                  Outc_nbit(j,k) = Out3_xnbits(i)
               endif
            enddo
         enddo
      enddo

! Print table of dynamic variables requested for output

      if (Lun_out.gt.0) then
         write(Lun_out,900)
         write(Lun_out,1006)
         write(Lun_out,901)
         cnt=0
         do j=1,Outd_sets
            do i=1,Outd_var_max(j)
               write(Lun_out,1007) Outd_var_S(i,j)(1:4),&
                 Outd_varnm_S(i,j)(1:16),Outd_nbit(i,j), &
                 Outd_filtpass(i,j),Outd_filtcoef(i,j),Level_typ_s(Outd_lev(j))
               cnt = cnt + 1
            enddo
         enddo
         write(Lun_out,1004) cnt 
         write(Lun_out,1006)
         write(Lun_out,2001)
      endif

      ixg = 0 ; ixgall = 0 ; rot = 0. ; rotall = 0.
      if (Ptopo_myproc.eq.0) then
         ixg(1,Ptopo_couleur+1)= Hgc_ig1ro
         ixg(2,Ptopo_couleur+1)= Hgc_ig2ro
         ixg(3,Ptopo_couleur+1)= Hgc_ig3ro
         ixg(4,Ptopo_couleur+1)= Hgc_ig4ro
         rot(1,Ptopo_couleur+1)= Grd_xlat1
         rot(2,Ptopo_couleur+1)= Grd_xlon1
         rot(3,Ptopo_couleur+1)= Grd_xlat2
         rot(4,Ptopo_couleur+1)= Grd_xlon2
      endif
      
      call RPN_COMM_allreduce ( ixg  , ixgall  ,       8,&
           RPN_COMM_INTEGER,"MPI_SUM",RPN_COMM_MULTIGRID,ierr )
      call RPN_COMM_allreduce ( rot  , rotall  ,       8,&
           RPN_COMM_REAL   ,"MPI_SUM",RPN_COMM_MULTIGRID,ierr )

      ierr = timestr2step(rsti,Fcst_rstrt_S,Cstv_dt_8)

      Out_unf= 0 ; Out_laststep_S = ' '
      Out_ixg(1:4) = ixgall(:,1)  ;  Out_ixg(5:8) = ixgall(:,2)
      Out_rot(1:4) = rotall(:,1)  ;  Out_rot(5:8) = rotall(:,2)
      Out_etik_S   = Out3_etik_s ; Out_gridtyp_S= 'E'
      Out_endstepno= min(Step_total+Step_initial, Lctl_step+rsti)
      Out_endstepno= Out_endstepno - Step_initial !in step_kount space
      Out_deet     = int(Cstv_dt_8)

      options = WB_REWRITE_NONE+WB_IS_LOCAL
      istat= wb_put('model/Output/etik', Out3_etik_S, options)      

  900 format(/'+',35('-'),'+',17('-'),'+',5('-'),'+'/'| DYNAMIC VARIABLES REQUESTED FOR OUTPUT              |',5x,'|')
  901 format('|',1x,'OUTPUT',1x,'|',2x,'   OUTCFG   ',2x,'|',2x,' BITS  |','FILTPASS|FILTCOEF| LEV |')
  902 format('|',1x,'OUTPUT',1x,'|',2x,'PHYSIC NAME ',2x,'|',2x,' BITS  |','FILTPASS|FILTCOEF| LEV |')
 1002 format('|',5X,a9,' Bus ',40x, '|')
 1004 format('|',5X,'Number of elements: ',i4,40x, '|')
 1006 format('+',8('-'),'+',16('-'),'+',9('-'),'+',8('-'),'+',8('-'),'+',5('-'))
 1007 format('|',2x,a4,2x,'|',a16,'|',i5,'    |',i8,'|',f8.3,'|',a4,' |')
 2001 format('* Note: NO filter is applied to 3D fields on M levels')
 3000 format(/,'SET_SOR - OUTPUT FILES will be in ',A8)
 5000 format( &
           ' TOTAL NUMBER OF WARNINGS ENCOUNTERED IN', &
           ' DIRECTIVE SETS: ', I5)
 5200 format(/,'INITIALIZATION OF OUTPUT PRODUCTS (S/R SET_SOR)', &
             /,'===============================================')
!
!-------------------------------------------------------------------
!
      return
      end