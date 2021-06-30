      SUBROUTINE baseflw(THLIQ,THICE,
     2        IG22,ILG,N,IGDR,
     3        I,IWT,DT,DELZW,THLMIN,
     4        ZNODE,FFF,ARE,LEG,SLP,SANI1,MLIQ,xsi,qbt,Isat,GRKSAT)  

      IMPLICIT NONE

      REAL   FFF
      INTEGER   IGDR  (ILG), I, K, IZ,
     1          IG,ILG,IWT(ILG),J,IG22,Isat(ILG),
     2          N
      REAL      DELZW(ILG,IG22),DT(ILG),THLMIN(ILG,IG22),
     1          THLIQ(ILG,IG22),THICE(ILG,IG22)
      REAL      GRKSAT(ILG,IG22)
      REAL      ZNODE(IG22),MLIQ(IG22),xsi(ILG,IG22)
      REAL      ANIS, qb(IG22),MLIQO,
     1          T(IG22), Zc, Fliq(IG22),qbt
      REAL      ARE(ILG), SLP(ILG),LEG(ILG)
      REAL      SANI1(ILG),qbtold
c      FFF=2.5
      IF(SANI1(I).GT.0) then
c         ANIS=50000.0
         ANIS=SANI1(I)
      ELSE
         ANIS=0.0
      ENDIF
      Zc=1.0
      MLIQO=0.0
      qbtold=qbt
      qbt=0.0
!      IF (ANIS.GT.0.0.AND.SLP(I).EQ.0.0.AND.LEG(I).EQ.0.0.AND.ARE(I)
      IF (ANIS.GE.0.0.AND.SLP(I).EQ.0.0.AND.LEG(I).EQ.0.0.AND.ARE(I)
     1    .EQ.0.0) THEN                
             SLP(I)=0.004045
             ARE(I)=2296.003
             LEG(I)=43418.28
      ENDIF

      DO 10 J=1,IG22                                                               
          qb(J)=0.0
          IF (J.LE.IGDR(I)) THEN
            Fliq(J)=1.0-THICE(I,J)/(THICE(I,J)+THLIQ(I,J))
          ELSE
            T(J)=0.0
            Fliq(J)=0.0
          ENDIF
          IF (IWT(I).LE.IGDR(I)) THEN
            IF  (J.LT.IWT(I)) THEN
                 T(J)=0.0
            ELSEIF (J.EQ.IWT(I)) THEN
                 T(J)=Fliq(J)*ANIS*GRKSAT(I,J)*exp(-FFF*(ZNODE(J)-
     1           Zc))/FFF*(exp(FFF*(ZNODE(J)-ZNODE(J-1))/2)-1.0)
            ELSEIF (J.GT.IWT(I).and.J.LE.IGDR(I)) THEN
                 T(J)=Fliq(J)*ANIS*GRKSAT(I,J)*exp(-FFF*(ZNODE(J)-
     1           Zc))/FFF*(exp(FFF*(ZNODE(J)-ZNODE(J-1)))-1.0)
            ENDIF
          ELSE
             T(J)=0.0
          ENDIF
          qb(J)=min((T(J)*tan((SLP(I)))*LEG(I))/(ARE(I)*1.E6)*1.E3,
     1               max(THLIQ(I,J)-THICE(I,J)-THLMIN(I,J),0.0)*
     2                DELZW(I,J)*1.E3/DT(I))
          if (qb(J).LT.0.0) print *,qb(J)
          IF(THICE(I,Isat(I)).GE.THLIQ(I,Isat(I))) qb(J)=0.0
          qbt=qbt+qb(J)
          IF (J.GE.IWT(I).and.J.LE.IGDR(I)) THEN
            MLIQO=MLIQ(J)
            MLIQ(J)=max(MLIQ(J)-qb(J)*DT(I),THLMIN(I,J)*DELZW(I,J)
     1              *1.0e3)
            xsi(I,J)=max(0.0,MLIQO-MLIQ(J))
          ENDIF
c      print *,qb(J),xsi(I,J)
10     CONTINUE 
c      if (qbt.GT.qbtold) then
c          print *, qbt, qbtold
c      endif
      RETURN                                                                      
      END


