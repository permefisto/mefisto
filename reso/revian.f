      SUBROUTINE REVIAN( NYOBJT, NUOBJT, XPT, YPT, ZPT, MNVIAN,
     %                   VITANG )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++012
C BUT :    REMPLIR LE VECTEUR VITANG AVEC LA VALEUR DES 3 COMPOSANTES
C -----    DE LA VITESSE ANGULAIRE LUE DANS a___vitesseangulaire
C
C ENTREES:
C --------
C NYOBJT : NUMERO DU TYPE OBJET (1:POINT, 2:LIGNE, ... )
C NUOBJT : NUMERO DE L'OBJET DANS SON LEXIQUE
C XPT,YPT,ZPT : LES 3 COORDONNEES DU POINT
C MNVIAN : ADRESSE MCN DU TABLEAU a___vitesseangulaire
C
C SORTIE :
C --------
C VITANG : VECTEUR(3) DE LA VITESSE ANGULAIRE
C          REMPLI PARTIELLEMENT DE 1 A NBVIAN OU TOTALEMENT
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++012
C AUTEUR: ALAIN PERRONNET LJLL UPMC & Saint Pierre du Perray  Avril 2010
C23456---------------------------------------------------------------012
      include"./incl/donflu.inc"
      include"./incl/a___vitesseangulaire.inc"
      include"./incl/ctemps.inc"
      DOUBLE PRECISION XPT, YPT, ZPT, XYZ(7)
      DOUBLE PRECISION VITANG(*)
      include"./incl/pp.inc"
      COMMON           MCN(MOTMCN)
      REAL             RMCN(1)
      EQUIVALENCE     (MCN(1),RMCN(1))
C
C     NOMBRE COMPOSANTES EN UN NOEUD DE LA VITESSE ANGULAIRE
      NBVIAN = MCN( MNVIAN + WBVIAN )
C
C     REMPLISSAGE SELON LE TYPE DE CONDITION INITIALE DE LA VITANG
      IF( MCN( MNVIAN + WTVIAN ) .EQ. 1 ) THEN
C
C        NBVIAN COMPOSANTES CONSTANTES DE LA VITESSE ANGULAIRE
C        =====================================================
         DO K=1,NBVIAN
            VITANG(K) = RMCN( MNVIAN + WAVIAN -1 + K )
         ENDDO
C
      ELSE
C
C        FONCTION UTILISATEUR VITANG(t,x,y,z,ntyobj,nuobj,nocomp)
C        ====================
         XYZ(1) = TEMPS
         XYZ(2) = XPT
         XYZ(3) = YPT
         XYZ(4) = ZPT
         XYZ(5) = NYOBJT
         XYZ(6) = NUOBJT
         DO K=1,NBVIAN
            XYZ(7) = K
            CALL FONVAL( MCN(MNVIAN+WFVIAN), 7, XYZ,
     %                   NCODEV, VITANG(K) )
         ENDDO
C
      ENDIF
C
      RETURN
      END
