      SUBROUTINE REONIN( NYOBJT, NUOBJT, XYZNO, MNOND0,  ONDEC0 )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++012
C BUT :    REMPLIR LES 2 COMPOSANTES DE L'ONDE COMPLEXE INITIALE
C -----    C-A-D A L'INSTANT TEMPS INITIAL PARTIE REELLE ET IMAGINAIRE
C          DANS LE CAS OU ELLE EST CONSTANTE EN TOUS LES NOEUDS
C          OU DEFINIE PAR UNE FONCTION UTILISATEUR EN UN NOEUD
C          DE COORDONNEES XNO,YNO,ZNO A L'INSTANT TEMPS
C
C ENTREES :
C ---------
C NYOBJT : NUMERO DU TYPE D'OBJET (1:POINT, 2:LIGNE, ..., 5:OBJET )
C NUOBJT : NUMERO DE L'OBJET DANS SON LEXIQUE
C XYZNO  : LES 3 COORDONNEES DU NOEUD
C MNOND0 : ADRESSE MCN DU TABLEAU 'ONDEINIT'
C
C SORTIE :
C --------
C ONDEC0 : TABLEAU (NDIM) DES NDIM COMPOSANTES L'ONDE
C          AU NOEUD XNO, YNO, ZNO A L'INSTANT TEMPS
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++012
C AUTEUR : ALAIN PERRONNET Texas A & M University at QATAR  Fevrier 2011
C23456---------------------------------------------------------------012
      include"./incl/donela.inc"
      include"./incl/a___ondeinit.inc"
      include"./incl/ctemps.inc"
C
      include"./incl/pp.inc"
      COMMON             MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (RMCN(1),MCN(1))
      DOUBLE PRECISION  XYZNO(3), ONDEC0(1:2), PARAM(7)
C
C     TYPE DES DONNEES L'ONDE OU GRADIENT DU DEPLACEMENT AU NOEUD
      LTOND0 = MCN( MNOND0 + WTOND0 )
C
C     REMPLISSAGE SELON LE TYPE
      IF( LTOND0 .EQ. 1 ) THEN
C
C        ONDE CONSTANTE
C        ==============
         DO I = 1, 2
            ONDEC0(I) = RMCN( MNOND0 + WAOND0 - 1 + I )
         ENDDO
C
      ELSE IF( LTOND0 .EQ. -1 ) THEN
C
C        FONCTION UTILISATEUR
C        ====================
         PARAM(1) = TEMPS
         PARAM(2) = XYZNO(1)
         PARAM(3) = XYZNO(2)
         PARAM(4) = XYZNO(3)
         PARAM(5) = NYOBJT
         PARAM(6) = NUOBJT
         DO I = 1, 2
C           LE NUMERO DE LA PARTIE RELLE ou IMAGINAIRE de L'ONDE
            PARAM(7) = I
            CALL FONVAL( MCN(MNOND0+WFOND0), 7, PARAM,
     %                   NCODEV, ONDEC0(I) )
         ENDDO
C
      ENDIF
C
      RETURN
      END
