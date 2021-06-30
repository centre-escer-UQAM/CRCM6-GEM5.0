***S/R CLASSD  -  Communication subroutine
      subroutine CLASSD (nom,valeur)
*
      character *(*) nom
      character      nomc*8
      real           valeur

*
*Author
*          B. Bilodeau (June 2004)
*
*Revisions
* 001      B. Dugas    (Aug 2007) - Simplify the interface
* 002      B. Dugas    (Jan 2009) - Add support for VAMIN
*
*Object
*          Communication subroutine used to define
*          certain REAL constants in the CLASS common
*          blocks (mainly in PHYCON and CLASSD2)
*
*Arguments
*
*          - Input -
* NOM      name of the option to be treated
* VALEUR   value of the constant
*
*Notes
*
*Implicites
*
**
*
      REAL            DELT,TFREZ
      COMMON /CLASS1/ DELT,TFREZ
      REAL            RGAS,RGASV,GRAV,SBC,VKC,CT,VMIN
      COMMON /CLASS2/ RGAS,RGASV,GRAV,SBC,VKC,CT,VMIN
      REAL            TCW,TCICE,TCSAND,TCCLAY,TCOM,TCDRYS
      COMMON /CLASS3/ TCW,TCICE,TCSAND,TCCLAY,TCOM,TCDRYS
      REAL            RHOSOL,RHOOM
      COMMON /CLASS3/ RHOSOL,RHOOM
      REAL            HCPW,HCPICE,HCPSOL,HCPOM,HCPSND,HCPCLY
      COMMON /CLASS4/ HCPW,HCPICE,HCPSOL,HCPOM,HCPSND,HCPCLY
      REAL            SPHW,SPHICE,SPHVEG,SPHAIR,RHOW,RHOICE
      COMMON /CLASS4/ SPHW,SPHICE,SPHVEG,SPHAIR,RHOW,RHOICE
      REAL            TCGLAC,CLHMLT,CLHVAP
      COMMON /CLASS4/ TCGLAC,CLHMLT,CLHVAP
      REAL            PI,GROWYR        ,ZOLNG,ZOLNS,ZOLNI
      COMMON /CLASS6/ PI,GROWYR(18,4,2),ZOLNG,ZOLNS,ZOLNI
      REAL            ZORAT   ,ZORATG
      COMMON /CLASS6/ ZORAT(4),ZORATG
      REAL            DELTA,CGRAV,CKARM,CPD
      COMMON /PHYCON/DELTA,CGRAV,CKARM,CPD
      REAL            AS,ASX,CI,BS,BETA,FACTN,HMIN,ANGMAX
      COMMON /CLASSD2/AS,ASX,CI,BS,BETA,FACTN,HMIN,ANGMAX
      REAL            VAMIN
      COMMON /CLASSD3/VAMIN
      integer          ifrsoil
      common /switches/ifrsoil
*
*     conversion de minuscules a majuscules
      call low2up(nom, nomc)
*
      if      (nomc.eq.'ANGMAX')                  THEN
*
            angmax = valeur
*
      else if (nomc.eq.'AS')                      THEN
*
            as = valeur
*
      else if (nomc.eq.'ASX')                     THEN
*
            asx = valeur
*
      else if (nomc.eq.'BETA')                    THEN
*
            beta = valeur
*
      else if (nomc.eq.'BS')                      THEN
*
            bs = valeur
*
      else if (nomc.eq.'CGRAV')                   THEN
*
            cgrav = valeur
*
      else if (nomc.eq.'CKARM')                   THEN
*
            ckarm = valeur
*
      else if (nomc.eq.'CI')                      THEN
*
            ci = valeur
*
      else if (nomc.eq.'CPD')                     THEN
*
            cpd = valeur
*
      else if (nomc.eq.'DELT')                    THEN
*
            delt = valeur
*
      else if (nomc.eq.'DELTA')                   THEN
*
            delta = valeur
*
      else if (nomc.eq.'FACTN')                   THEN
*
            factn = valeur
*
      else if (nomc.eq.'GRAV')                    THEN
*
            grav = valeur
*
      else if (nomc.eq.'HMIN')                    THEN
*
            hmin = valeur
*
      else if (nomc.eq.'IFRSOIL')                 THEN
*
            ifrsoil = valeur
*
      else if (nomc.eq.'PI')                      THEN
*
            pi = valeur
*
      else if (nomc.eq.'RGAS')                    THEN
*
            rgas = valeur
*
      else if (nomc.eq.'RGASV')                   THEN
*
            rgasv = valeur
*
      else if (nomc.eq.'SBC')                     THEN
*
            sbc = valeur
*
      else if (nomc.eq.'SPHAIR')                  THEN
*
            sphair = valeur
*
      else if (nomc.eq.'TFREZ')                   THEN
*
            tfrez = valeur
*
      else if (nomc.eq.'VAMIN')                   THEN
*
            vamin = valeur
*
      else                                                 
            write(6,1010) nomc
            call xit('CLASSD',-1)
      endif
*
*
1010   FORMAT ( ' *****************************************',
     +        / ' *****************************************',
     +        / ' *                                       *',
     +        / ' ***** ABORT ***** ABORT ***** ABORT *****',
     +        / ' *                                       *',
     +        / ' *  ', A8, 'IS AN INVALID OPTION         *',
     +        / ' *     OF SUBROUTINE CLASSD              *',
     +        / ' *                                       *',
     +        / ' *****************************************',
     +        / ' *****************************************')
*
*
      RETURN
      END
