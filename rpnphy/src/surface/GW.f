      SUBROUTINE GW(THLIQ,THICE,DT,THPOR,
     1              BI,PSISAT,GRKSAT,DELZW,
     2              IGRD,IG22,IGP1,IGP2,ILG,IL1,IL2,JL,N,IGDR,
     3              ZBOTW,THPORF,I,GRKSATF,THLMIN,THLMAX,IWT,QDIS1,
     4              QIN3,WT1,WTnew,DMLIQ,igrinfl,KQIN,Ibed,DMLIQJ,
     5              Isat,xsi,ZWT,TBARW,EXCW,ARE,LEG,SLP,SANI,igwscheme)  

      IMPLICIT NONE

      integer   igwscheme
      REAL      FFF, CMIC, SMPFZ, WH_ZWT,A1,A2,A3,WTold,
     1          WH,QIN3,WT(ILG),WA,WS,WTSUB,ZD,beta
      REAL      RSBMAX, ZLLS, ZLLS_IWT
      REAL      FCRMAX
      REAL      TIMEAN, ZWT(ILG),ZWT1, QDIS1, ROUS,
     1          S_NODE, XS
      INTEGER   IGRD  (ILG),    IGDR  (ILG), I, K, IZ, ICK,
     1          IG,IGP1,ILG,IWT(ILG),J,Isat(ILG),IG22,jindex,
     2          igrinfl, KQIN(ILG), Ibed(ILG), IL1,IL2,JL,N, IGP2,
     3          jji,JM,JK
      REAL      THPORF(ILG,IG22), DELZW(ILG,IG22),THPOR(ILG,IG22),
     1          THLIQ(ILG,IG22),THICE(ILG,IG22),THLMIN(ILG,IG22),
     2          THLMAX(ILG,IG22), GRKSAT(ILG,IG22),TBARW(ILG,IG22)
      REAL      BI(ILG,IG22), PSISAT(ILG,IG22),EXCW(ILG)
      REAL      ZBOTW (ILG,IG22), GRKSATF(ILG,IG22)
      REAL      SMC(IG22) ,ZNODE(IG22),  DT(ILG)
      REAL      GRKAQ(ILG),DTHGW,DMLIQ,DMLIQJ(ILG),cumdelz
      REAL      WT1(ILG),WTnew(ILG),H2,Waq,cumdelzw(IG22),MLIQJ,MLIQk   
      REAL      MLIQ(IG22),CAP,QINold,icefrac(ILG,IG22),icefracsum,MLIQI           
      REAL      MLIQold,maxd,xsi(ILG,IG22), GRKL(IG22)
      REAL      sdlz,dumw,RW(IG22),SANI(ILG),
     3          ARE(ILG), SLP(ILG),    LEG(ILG)
      ICK=igwscheme-1 ! Switch for CKGW=1 and CKGWL=2: 
      FFF=2.50    ! decay factor
      beta=8.0    !
      RSBMAX=5.   ! the maximum sub-surface runoff when the grid cell average water table depth is zero. This must be calibrated for the model later 
      FCRMAX=0.0  ! impermeable fraction of frozen soil 
      TIMEAN=10.2 ! Grid-cell mean topographic index
      ROUS =0.2   ! Specific yeild
      CMIC=0.40   ! Micropore content (0-1) 0-close to free drainage
      maxd=0.30
      IF (DT(I).LT.1800. .AND. igrinfl.EQ.1)  THEN
          WT(I)=WTnew(I)*1000.0
      ELSE   
          WT(I)=WT1(I)*1000.0
      ENDIF    ! m----> mm 
!************************************************
!     Initial calculations
!************************************************     
      
      IG=IGDR(I)

      do 10 J=1,IG22
         xsi(I,J)=0.0
         SMC(J)=THLIQ(I,J)  ! Total soil moisture content, mm 
         MLIQ(J)=THLIQ(I,J)*DELZW(I,J)*1.E3  ! Liquid water content at each layer, mm
         GRKL(J)=0.0
         IF (THPORF(I,J).GT.0.0) GRKL(J)=MIN(GRKSATF(I,J)*(min(1.0,
     1      (THLIQ(I,J))/THPORF(I,J)))**(2.*BI(I,J)+3.),GRKSATF(I,J))
!           IF (THPORF(I,J).GT.0.0) GRKL(J)=MIN(GRKSATF(I,J)*(min(1.0,
!           1      (THLIQ(I,J)+THICE(I,J))/
!           1      THPORF(I,J)))**(2.*BI(I,J)+3.),GRKSATF(I,J))
c           IF (THPORF(I,J).GT.0.0) GRKL(J)=GRKSATF(I,J)
10    enddo

      ZNODE(1)=DELZW(I,1)/2.0           ! Layer depth at center point of layer, m
      cumdelzw(1)=DELZW(I,1)
      DO 21 J=2,IG22
         cumdelzw(J)=cumdelzw(J-1)+DELZW(I,J)
21    enddo

      do 20 J=2,IG22
         ZNODE(J) = cumdelzw(J-1)+0.5*DELZW(I,J)  ! m
20    enddo
!*********************************************************
!  the cumulative depth of the layers up to soil layer bottom
!  that is the bedrock or groundwater depth
!********************************************************* 
      ZLLS=0.0
      icefracsum =0.0
      do k=1,IGDR(I)      ! atleast for 5 layers as min level for WT is 1.5 m
          ZLLS=ZLLS+DELZW(I,k)    ! m
C         icefrac(I,k)=min(1.0,THICE(I,k)/(THICE(I,k)+THLIQ(I,k)))
C
c         IF(THPOR(I,k).GT.0.0.and.k.lt.Ibed(I)) THEN
c         IF (Ibed(I).ge.3) THEN
c            icefrac(I,k)=min(1.0,THICE(I,k)/THPOR(I,3))
c         ELSE
c            icefrac(I,k)=min(1.0,THICE(I,k)/THPOR(I,MIN(Ibed(I),k)))
c         ENDIF
          IF((THICE(I,k)+THLIQ(I,k)).GT.0.0) THEN
              icefrac(I,k)=min(1.0,THICE(I,k)/(THICE(I,k)+THLIQ(I,k)))
          ELSE
              icefrac(I,k)=0.0
          ENDIF
              icefracsum = icefracsum + icefrac(I,k) * DELZW(I,k)  
      enddo
      ZLLS=-ZLLS

!*********************************

       IWT(I)=Isat(I)+1
!                         mustbe careful with this as Isat must return you the level of water table and 
!			not the position for aquifer. So make sure that Isat is the place for the water table (layer number)
!***********************************************
! Finding WA at the previous stage as a function of 
! WT at the previuos stage, and before any change in 
! WT
!***********************************************
      IF (IWT(I).LE.IG) THEN
         WS=0.0
            DO J = IWT(I),IG
               WS=WS + THPORF(I,J) * DELZW(I,J)*1.E3  !mm
            ENDDO
         WA=WT(I)-WS                                 !mm
      ELSE
         WA=WT(I)     ! this also be activated for the case IG<5
      ENDIF
!************************************************  
!     Finding DISCHARGE OUT AND IN 
!************************************************
      FCRMAX = max(0.0,exp(-3.0*(1.0-(icefracsum/ABS(ZLLS))))-exp(-3.0))
      QDIS1=(1.0-FCRMAX)*RSBMAX*EXP(-TIMEAN)*EXP(-FFF*(ZWT(I)))
      IF(THICE(I,Isat(I)).GE.THLIQ(I,Isat(I))) QDIS1=0.0       
      S_NODE = MIN(1.0,SMC(Isat(I))/THPORF(I,Isat(I))) 
      S_NODE = MAX(S_NODE,REAL(0.01,KIND=8)) 
      SMPFZ = -PSISAT(I,Isat(I))*1000.*S_NODE**(-BI(I,Isat(I))) 
     1        *(1.0+beta*THICE(I,Isat(I)))**2.0 ! m --> mm 
      SMPFZ = MAX(-80000.0,CMIC*SMPFZ) 
      WH_ZWT = - ZWT(I) * 1.E3 !(mm)

      IF(DELZW(I,Isat(I)).LE.0.00005.AND.Isat(I).GE.2) THEN
          GRKAQ(I)=GRKL(Isat(I)-1)*(1.0-EXP(-FFF*(ZWT(I)-ZNODE(Isat(I)-
     1        1))))/(FFF*(ZWT(I)-ZNODE(Isat(I)-1)))                      ! m/s  
          WH = SMPFZ - ZNODE(Isat(I)-1)*1.E3    ! (mm)    
          QIN3=-GRKAQ(I)*(WH_ZWT-WH)/
     1    ((ZWT(I)-ZNODE(Isat(I)-1))*1.E3)*1.E3 !mm/s          
      ELSE
          GRKAQ(I)=GRKL(Isat(I))*(1.0-EXP(-FFF*(ZWT(I)-ZNODE(Isat(I))))
     1         )/(FFF*(ZWT(I)-ZNODE(Isat(I))))                           ! m/s
          WH = SMPFZ - ZNODE(Isat(I))*1.E3    ! (mm)   
c$$$          IF (ZWT.le.ZNODE(Isat(I))) then
c$$$              ZWT=ZWT
c$$$              print *, 'check Isat and ibed'
c$$$          endif
          QIN3=-GRKAQ(I)*(WH_ZWT-WH)/
     1    ((ZWT(I)-ZNODE(Isat(I)))*1.E3)*1.E3
      ENDIF

      IF(THICE(I,Isat(I)).GE.THLIQ(I,Isat(I))) QIN3=0.0
      QIN3=MAX(-10.0/DT(I),MIN(5./DT(I),QIN3))          !mm/s
!************************************************
! Water storage in the aquifer + saturated soil
!************************************************
      WTold=MAX(WT(I),5000.)
      WT(I)=WT(I)+(QIN3-QDIS1)*DT(I)!(mm)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  if the water table is within the aquifer
      H2=0.0
      Waq=0.0
      MLIQJ=0.0
      MLIQk=0.0
      DMLIQJ(I)=0.0
      WTSUB=0.0
      XS=0.0 
      IF(IWT(I).EQ.IG+1) THEN 
         WA=WA+(QIN3-QDIS1)*DT(I) !(mm) 
         WT(I)=WA 
         ZWT1= (-ZLLS+ 25.) - WA/1000./ROUS!(m) 

! 
!    ZWT correction for maxd
!
         IF(ZWT1.LE.maxd)  THEN      
            H2=maxd-ZWT1
            Waq=H2*0.2           
            IF (THLIQ(I,Isat(I)).GT.THICE(I,Isat(I))) THEN
!               QDIS1=QDIS1+((Waq)/DT(I))*1.0e3
               EXCW(I)=EXCW(I)+((Waq)/DT(I))*1.0e3
               WA=MAX(WA-Waq*1000.0,0.0)                ! mm
               WT(I)=WA 
            ELSEIF(QIN3.GT.0.0.AND.QIN3.GE.((Waq)*1.0e3/DT(I)))  THEN
               QIN3=QIN3-((Waq)/DT(I))*1.0e3
               WA=MAX(WA-Waq*1000.0,0.0)                ! mm
               WT(I)=WA                        ! mm     
            ELSEIF(QIN3.GT.0.0.AND.QIN3.LE.((Waq)*1.0e3/DT(I))) THEN
               QIN3=0.0
               WA=MAX(WA-Waq*1000.0,0.0)    ! here assumed that we can reduced it, although max reduction is QIN3 
               WT(I)=WA            ! mm
            ENDIF
            ZWT1= maxd     !m
         ENDIF
! 
!    END of ZWT correction for maxd
!
         MLIQ(IG) = MLIQ(IG) - QIN3 * DT(I)! [mm]
         MLIQ(IG) = MLIQ(IG) + MAX(0.,(WA - WTold)) ! mm

         if (ig.gt.1) then
         MLIQJ=MLIQ(IG-1)+MAX((MLIQ(IG)-THLMAX(I,IG)*
     1         DELZW(I,IG)*1.E3),0.0)
         IF (MLIQ(IG).GT.(THLMAX(I,IG)*DELZW(I,IG)*1.E3)
     1       .AND.MLIQJ.LT.(THLMAX(I,IG-1)*DELZW(I,IG-1)*1.E3))
     2        THEN
             MLIQ(IG-1)=MLIQ(IG-1)+(MLIQ(IG)-THLMAX(I,IG)*
     1                  DELZW(I,IG)*1.E3)
             DMLIQJ(I)=(MLIQ(IG)-THLMAX(I,IG)*DELZW(I,IG)*1.E3)/1000.
             MLIQ(IG)=THLMAX(I,IG)*DELZW(I,IG)*1.E3
         ENDIF

         MLIQk=MLIQ(IG-1)-MAX((-MLIQ(IG)+THLMIN(I,IG)*
     1         DELZW(I,IG)*1.E3),0.0)         
         IF (MLIQ(IG).LT.(THLMIN(I,IG)*DELZW(I,IG)*1.E3).AND.
     1       MLIQk.GT.(THLMIN(I,IG-1)*DELZW(I,IG-1)*1.E3)
     2       ) THEN
             MLIQ(IG-1)=MLIQ(IG-1)+(MLIQ(IG)-THLMIN(I,IG)*
     1                  DELZW(I,IG)*1.E3)
             DMLIQJ(I)=(MLIQ(IG)-THLMIN(I,IG)*DELZW(I,IG)*1.E3)/1000.
             MLIQ(IG)=THLMIN(I,IG)*DELZW(I,IG)*1.E3
         ENDIF
         endif
         DMLIQ=MAX(0.,(WA - WTold))/1000.0
         WA= MIN(WA, 5000.)
         KQIN(I)=1
      ELSE
! if water table is within in bottom soil layer
        IF (IWT(I).EQ.IG) THEN 
            do 322 JM=Ibed(I),1,-1
               dumw=0.0
               sdlz=0.0
               do 323 JK=JM+1,Ibed(I)
                  dumw=dumw+(DELZW(I,JK)*(THPORF(I,JK)
     1            ))
                  sdlz=sdlz+DELZW(I,JK)
323            ENDDO
               IF(THPORF(I,JM).GT.0.0) THEN
                    RW(JM)=(WT(I)/1000.0-25.0*.2-dumw)
     1                    /(THPORF(I,JM))          ! this finds water remained for the next (up) layer
               ELSE
                    RW(JM)=RW(JM+1)   ! this mean the same water is still available when layer frozen compeletly
                    if (jm.eq.ig22) then
                      write(6,*),'GW.f: j=ig at ',i,
     &                           ' ; RW = ',RW(JM+1)
                      RW(jm)=0.0
                    endif
               ENDIF
               IF(RW(JM).LT.0.0) EXIT
322         ENDDO
            IF(JM.EQ.Ibed(I).AND.RW(JM).LT.0.0) THEN 
              ZWT1=-ZLLS-RW(JM)
            ELSE
              ZWT1=-ZLLS+(-RW(JM+1)-sdlz+DELZW(I,JM+1))
            ENDIF
!*******************************************
!    ZWT correction for maxd
!
           IF(ZWT1.LE.maxd.AND.ZBOTW(I,IGDR(I)).LT.maxd.AND.
     1                        ZBOTW(I,IGDR(I)).GT.ZWT1) THEN
             Waq=(maxd-ZBOTW(I,IGDR(I)))*0.2
             WS=WT(I)-WA
             H2=ZBOTW(I,IGDR(I))-ZWT1
                 EXCW(I)=EXCW(I)+((WS+Waq)*1.0e3/DT(I))
                 WT(I)=WT(I)-WS-Waq
                 IF (WT(I).LT.5000.0) THEN
                     WA=WT(I)
                 ELSE
                     WA=5000.0
                 ENDIF
             ZWT1=maxd
           ELSEIF(ZWT1.LE.maxd.AND.ZBOTW(I,IGDR(I)).GT.maxd) THEN
             H2=maxd-ZWT1
             Waq=H2*0.2
             IF (THLIQ(I,Isat(I)).GT.THICE(I,Isat(I))) THEN
                EXCW(I)=EXCW(I)+((Waq)/DT(I))*1.0e3
                WT(I)=MAX(WT(I)-Waq*1000.0,0.0)
                IF(WT(I).LT.5000.0) THEN
                   WA=WT(I)
                ELSE
                   WA=5000.0
                ENDIF
             ELSEIF(QIN3.GT.0.0.AND.QIN3.GE.((Waq)*1.0e3/DT(I)))  THEN
                QIN3=QIN3-((Waq)/DT(I))*1.0e3
                WT(I)=MAX(WT(I)-Waq*1000.0,0.0)
                IF(WT(I).LT.5000.0) THEN
                   WA=WT(I)
                ELSE
                   WA=5000.0
                ENDIF
             ELSEIF(QIN3.GT.0.0.AND.QIN3.LE.((Waq)*1.0e3/DT(I))) THEN
                QIN3=0.0
                WT(I)=MAX(WT(I)-Waq*1000.0,0.0)
                IF(WT(I).LT.5000.0) THEN
                   WA=WT(I)
                ELSE
                   WA=5000.0
                ENDIF
             ENDIF
             ZWT1= maxd     !m        
           ENDIF
! 
!    END of ZWT correction for maxd
! 
          DMLIQ=0.0 ! CLM3.5 mondel assume that liq water in soil works independently from 
                       ! water table for this case
           KQIN(I)=0

! if water table is within in soil layer
        ELSE
           WS=0.
           DO IZ = IWT(I)+1,IG
              WS=WS + THPORF(I,IZ) * DELZW(I,IZ)*1.E3   !(mm)
           ENDDO

           DMLIQ=0.0
           KQIN(I)=0
           ZLLS_IWT=0.0
           do k=1,IWT(I)+1 
               ZLLS_IWT=ZLLS_IWT+DELZW(I,k)               ! (m)
           enddo   

           ZLLS_IWT=-ZLLS_IWT
            do 422 JM=Ibed(I),1,-1
               dumw=0.0
               sdlz=0.0
               do 423 JK=JM+1,Ibed(I)
                  dumw=dumw+DELZW(I,JK)*MAX(THPORF(I,JK)
     1            ,0.0)
                  sdlz=sdlz+DELZW(I,JK)
423            ENDDO
               IF(THPORF(I,JM).GT.0.0) THEN
                    RW(JM)=(WT(I)/1000.0-25.0*.2-dumw)
     1                    /(THPORF(I,JM))          ! this finds water remained for the next (up) layer
               ELSE
                    RW(JM)=RW(JM+1)   ! this mean the same water is still available when layer frozen compeletly
                    if (jm.eq.ig22) then
                      write(6,*),'GW.f(2): j=ig at ',i,
     &                           ' ; RW = ',RW(jm+1)
                      RW(jm)=0.0
                    endif
               ENDIF
               IF(RW(JM).LT.0.0) EXIT
422         ENDDO
            IF(JM.EQ.Ibed(I).AND.RW(JM).LT.0.0) THEN
              ZWT1=-ZLLS-RW(JM)
            ELSE
              ZWT1=-ZLLS+(-RW(JM+1)-sdlz+DELZW(I,JM+1))
            ENDIF

!  **************************************************
!    ZWT correction for maxd
!
           IF(ZWT1.LE.maxd.AND.ZBOTW(I,IGDR(I)).LT.maxd.AND.
     1                        ZBOTW(I,IGDR(I)).GT.ZWT1) THEN
             Waq=(maxd-ZBOTW(I,IGDR(I)))*0.2
             WS=WT(I)-WA
             H2=ZBOTW(I,IGDR(I))-ZWT1
             IF (THLIQ(I,Isat(I)).GT.THICE(I,Isat(I))) THEN
                EXCW(I)=EXCW(I)+((WS+Waq)*1.0e3/DT(I))
                WT(I)=WT(I)-WS-Waq
                IF(WT(I).LT.5000.0) THEN
                   WA=WT(I)
                ELSE
                   WA=5000.0
                ENDIF
             ELSEIF(QIN3.GT.0.0.AND.QIN3.GE.((WS+Waq)*1.0e3/DT(I))) THEN
                QIN3=QIN3-((WS+Waq)*1.0e3/DT(I))
                WT(I)=WT(I)-WS-Waq
                IF(WT(I).LT.5000.0) THEN
                   WA=WT(I)
                ELSE
                   WA=5000.0
                ENDIF
             ELSEIF(QIN3.GT.0.0.AND.QIN3.LE.((WS+Waq)*1.0e3/DT(I))) THEN
                QIN3=0.0
                WT(I)=WT(I)-WS-Waq
                IF(WT(I).LT.5000.0) THEN
                   WA=WT(I)
                ELSE
                   WA=5000.0
                ENDIF
             ENDIF
             ZWT1=maxd
           ELSEIF(ZWT1.LE.maxd.AND.ZBOTW(I,IGDR(I)).GT.maxd) THEN
             H2=maxd-ZWT1
             Waq=H2*0.2
             IF (THLIQ(I,Isat(I)).GT.THICE(I,Isat(I))) THEN
                EXCW(I)=EXCW(I)+((Waq)/DT(I))*1.0e3
                WT(I)=MAX(WT(I)-Waq*1000.0,0.0)
                IF(WT(I).LT.5000.0) THEN
                   WA=WT(I)
                ELSE
                   WA=5000.0
                ENDIF
             ELSEIF(QIN3.GT.0.0.AND.QIN3.GE.((Waq)*1.0e3/DT(I)))  THEN
                QIN3=QIN3-((Waq)/DT(I))*1.0e3
                WT(I)=MAX(WT(I)-Waq*1000.0,0.0)
                IF(WT(I).LT.5000.0) THEN
                   WA=WT(I)
                ELSE
                   WA=5000.0
                ENDIF
             ELSEIF(QIN3.GT.0.0.AND.QIN3.LE.((Waq)*1.0e3/DT(I))) THEN
                QIN3=0.0
                WT(I)=MAX(WT(I)-Waq*1000.0,0.0)
                IF(WT(I).LT.5000.0) THEN
                   WA=WT(I)
                ELSE
                   WA=5000.0
                ENDIF
             ENDIF
            ZWT1= maxd     !m        
           ENDIF
!    End of ZWT correction
!
        ENDIF

        CAP=5000.-MAX(0.0,(maxd-ZBOTW(I,IGDR(I)))*0.2*1000.)
        DO J=3,IG
            CAP=CAP+max(THPORF(I,J)*DELZW(I,J),0.0)*1000.
        ENDDO
        WT(I)=max(WT(I),0.0)
        IF(WT(I).GT.CAP) THEN
           WT(I)=min(WT(I),CAP)
      ENDIF

c$$$        WTSUB = 0. 
c$$$        DO J = IG,IG
c$$$            WTSUB = WTSUB + GRKSATF(I,J)*DELZW(I,J)
c$$$        END DO
            WTSUB =GRKSATF(I,IG)*DELZW(I,IG)

!***********************************************************
!!  subsurface runoff from saturated layer
!*********************************************************** 
         IF(ICK.EQ.2) THEN
             call baseflw(THLIQ,THICE,
     2       IG22,ILG,N,IGDR,
     3       I,IWT,DT,DELZW,THLMIN,
     4       ZNODE,FFF,ARE,LEG,SLP,SANI,MLIQ,xsi,QDIS1,Isat,GRKSAT)
         ELSE
             DO J = IG, IG
                MLIQ(J) = MLIQ(J) - QDIS1*DT(I)*GRKSATF(I,J)*DELZW(I,J)
     1                    /WTSUB
                xsi(I,J)= QDIS1*DT(I)*GRKSATF(I,J)*DELZW(I,J)/WTSUB
             END DO
         ENDIF

      END IF

!************************************************  
!     Updating soil moisture content 
!************************************************
       DO J = IWT(I), IG-1
       
          IF (MLIQ(J) .LT.THLMIN(I,J)*DELZW(I,J)*1.E3 ) THEN 
             XS = (THLMIN(I,J))*DELZW(I,J)*1.E3-MLIQ(J)   !(mm)
          ELSE 
             XS = 0.0
          END IF 
       
          MLIQ(J) = MLIQ(J)+XS                          !(mm)
          WT(I)= WT(I)-XS                     ! (mm)
          IF (WT(I).LE.5000.0) WA= WA-XS 
          IF (KQIN(I).EQ.0) THEN
             QDIS1=MAX(QDIS1-XS/DT(I),0.0)
             xsi(I,J)=MAX(0.0,xsi(I,J)-XS)
          ENDIF
       END DO
        
c      Updating and checking the bottom layer 
       IF (MLIQ(IG) .LT. THLMIN(I,IG)*DELZW(I,IG)*1.E3) THEN
          IF(QIN3.GT.0.0.AND.KQIN(I).GT.0) THEN
             XS = (THLMIN(I,IG)+0.00001)*DELZW(I,IG)*1.E3-MLIQ(IG)    !(mm)
             QIN3=MAX(0.0,QIN3-XS/DT(I))
          ELSEIF(KQIN(I).EQ.0) THEN
              XS = (THLMIN(I,IG)+0.00001)*DELZW(I,IG)*1.E3-MLIQ(IG)    !(mm)
              QDIS1=MAX(0.0,QDIS1-XS/DT(I))
              xsi(I,IG)=MAX(0.0,xsi(I,IG)-XS)
          ENDIF
       ELSE
          XS = 0. 
       END IF

       MLIQ(IG)= MLIQ(IG) + XS             !(mm)
       WT(I)= WT(I)-XS                     ! (mm)
       IF (WT(I).LE.5000.0) WA= WA-XS

       IF (MLIQ(IG) .GT. (THLMAX(I,IG))*DELZW(I,IG)*1.E3) THEN
          XS = (THLMAX(I,IG))*DELZW(I,IG)*1.E3-MLIQ(IG)    !(mm)
          IF(QIN3.LT.0.0.AND.KQIN(I).EQ.1) THEN
             QIN3=(QIN3-XS/DT(I))
          ENDIF
          If (KQIN(I).EQ.0)  THEN
               QDIS1=MAX(QDIS1-XS/DT(I),0.0)
               xsi(I,IG)=xsi(I,IG)-XS/DT(I)
          ENDIF
       ELSE
          XS = 0.
       END IF
       MLIQ(IG)= MLIQ(IG) + XS             !(mm)
       WT(I)= WT(I)-XS                     ! (mm)
       
       DO J = 1,IG
         THLIQ(I,J)= MLIQ(J) /(DELZW(I,J)*1.E3)   !(mm)/mm
       END DO
       DTHGW=DMLIQ/(DELZW(I,IG)*1.E3)       
       WTnew(I)=WT(I)*1.0E-3       
      RETURN                                                                      
      END
         



