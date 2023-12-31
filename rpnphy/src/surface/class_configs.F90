#if defined(DOC)
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
!Revisions
! 001      L. Duarte   (Dec 2008) - CLASS_IG is now a user-defined option
!                                   and is declared in options.cdk
! 002      K. Winger   (Jun 2020) - Turned into module
#endif

module class_configs
!     parameters of land scheme CLASS
!     integer CLASS_IC, CLASS_IG
      integer CLASS_IC
      integer CTEM_ICC
!
!     number of vegetation classes that have a special treatment
      parameter (CLASS_IC=4)
!
!     number of vegetation classes (PFTs) in CTEM
      parameter (CTEM_ICC=9)
!
!     other parameters used for the correspondence between CTEM and CLASS PFTs:
!     number of CTEM PFTs within each CLASS PFT
      INTEGER, PARAMETER :: NOL2PFT(CLASS_IC)=(/2,3,2,2/)

!     index of the first CTEM PFT within each CLASS PFT
      INTEGER, PARAMETER :: FIRSTPFT(CLASS_IC)=(/1,3,6,8/)

end module
