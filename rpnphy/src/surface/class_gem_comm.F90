!S/R class_GEM_comm  -  Communication subroutine (previously called "CLASSD")
subroutine class_gem_comm (nom,valeur)
!
  use classicParams, only: DELT, VAMIN, &
                           SPHAIR, GRAV, VKC, RGAS, RGASV, SBC, TFREZ, &
                           BETA, FACTN, HMIN, PI
!
      character *(*) nom
      character      nomc*8
      real           valeur

!
!Author
!          B. Bilodeau (June 2004)
!
!Revisions
! 001      B. Dugas    (Aug 2007) - Simplify the interface
! 002      B. Dugas    (Jan 2009) - Add support for VAMIN
! 003      K. Winger   (May 2020) - Renamed to class_GEM_comm
!                                 - Adapted to CLASSIC
!
!Object
!          Communication subroutine used to define
!          certain REAL constants in the CLASS common
!          blocks (mainly in PHYCON and CLASSD2)
!
!Arguments
!
!          - Input -
! NOM      name of the option to be treated
! VALEUR   value of the constant
!
!Notes
!
!Implicites
!

!
!     conversion de minuscules a majuscules
      call low2up(nom, nomc)

print *,'Transferred to CLASS: ',nomc, ' = ', valeur

!
      if      (nomc.eq.'ANGMAX')                  THEN
!
            angmax = valeur
!
      else if (nomc.eq.'AS')                      THEN
!
            as = valeur
!
      else if (nomc.eq.'ASX')                     THEN
!
            asx = valeur
!
      else if (nomc.eq.'BETA')                    THEN
!
            beta = valeur
!
      else if (nomc.eq.'BS')                      THEN
!
            bs = valeur
!
      else if (nomc.eq.'VKC')                     THEN
!
            vkc = valeur
!
      else if (nomc.eq.'CI')                      THEN
!
            ci = valeur
!
      else if (nomc.eq.'CPD')                     THEN
!
            cpd = valeur
!
      else if (nomc.eq.'DELT')                    THEN
!
            delt = valeur
!
      else if (nomc.eq.'DELTA')                   THEN
!
            delta = valeur
!
      else if (nomc.eq.'FACTN')                   THEN
!
            factn = valeur
!
      else if (nomc.eq.'GRAV')                    THEN
!
            grav = valeur
!
      else if (nomc.eq.'HMIN')                    THEN
!
            hmin = valeur
!
      else if (nomc.eq.'PI')                      THEN
!
            pi = valeur
!
      else if (nomc.eq.'RGAS')                    THEN
!
            rgas = valeur
!
      else if (nomc.eq.'RGASV')                   THEN
!
            rgasv = valeur
!
      else if (nomc.eq.'SBC')                     THEN
!
            sbc = valeur
!
      else if (nomc.eq.'SPHAIR')                  THEN
!
            sphair = valeur
!
      else if (nomc.eq.'TFREZ')                   THEN
!
            tfrez = valeur
!
      else if (nomc.eq.'VAMIN')                   THEN
!
            vamin = valeur
!
      else                                                 
            write(6,1010) nomc
            call errorHandler('class_gem_comm',-1)
      endif
!
!
1010   FORMAT ( ' *****************************************', &
              / ' *****************************************', &
              / ' *                                       *', &
              / ' ***** ABORT ***** ABORT ***** ABORT *****', &
              / ' *                                       *', &
              / ' *  ', A8, 'IS AN INVALID OPTION         *', &
              / ' *     OF SUBROUTINE CLASSD              *', &
              / ' *                                       *', &
              / ' *****************************************', &
              / ' *****************************************')
!
!
  RETURN
END
