      SUBROUTINE REBLVI( NYOBJT, NUOBJT, XPI,YPI,ZPI, MNBLVI,
     %                   NBCOBV, BLVITE )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++012
C BUT :    REMPLIR LE TABLEAU BLVITE DES VITESSES BLOQUEES DU FLUIDE
C -----
C
C ENTREES:
C --------
C NYOBJT : NUMERO DU TYPE OBJET (1:POINT, 2:LIGNE, ... )
C NUOBJT : NUMERO DE L'OBJET DANS SON LEXIQUE
C XPI,YPI,ZPI : LES 3 COORDONNEES DU POINT D'INTEGRATION NUMERIQUE
C MNBLVI : ADRESSE MCN DU TABLEAU 'BLVITESSE'
C
C SORTIE :
C --------
C NBCOBV : LE NOMBRE DE COMPOSANTES DE LA VITESSE A FIXER
C BLVITE : LE TABLEAU (NBCOBV) DES VITESSES EXERCEES SUR L'OBJET
C
C ATTENTION: SEULES NBCOBV COMPOSANTES DES DEPLACEMENTS PARMI
C ---------- LES NDIM POSSIBLES SONT INITIALISEES!
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++012
C AUTEUR : BENHAMADOUCHE SOFIANE                            Janvier 2000
C AJOUT  : ALAIN PERRONNET LJLL UPMC & St Pierre du Perray Decembre 2012
C....6...............................................................012
      include"./incl/a___blvitesse.inc"
      include"./incl/ctemps.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL             RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
C
      DOUBLE PRECISION XPI,YPI,ZPI,XYZ(7)
      DOUBLE PRECISION BLVITE(3)
      EQUIVALENCE     (MCN(1),RMCN(1))
C
C     LE NOMBRE DE COMPOSANTES DE VITESSES A DEFINIR
      NBCOBV = MCN( MNBLVI + WBCOBV )
C
C     REMPLISSAGE SELON LE TYPE DES DONNEES DES VITESSES
      MN = MNBLVI + WUCOBV + NBCOBV
C
C     type of the imposed velocity
      LTBLVI = MCN( MNBLVI + WTBLVI )

ccc   variable LTBLVI 'type vitesse fixee' entier 
ccc   ( 1 : 'vitesse fixee constante'
ccc   ,-1 : 'fonction vitesse fixee'
ccc   , 0 : 'vitesse fixee a supprimer'
ccc   , 2 : 'vitesse convectee au temps tn' ) ;

      IF( LTBLVI .EQ. 1 ) THEN
C
C        VALEURS CONSTANTES POUR CET OBJET
C        =================================
         DO J=1,NBCOBV
C           LA VALEUR DE BLOCAGE DE LA VITESSE EST EMPILEE
            BLVITE( J ) = RMCN( MN - 1 + J )
         ENDDO
C
      ELSE IF( LTBLVI .EQ. 2 ) THEN
C
C        VITESSE CONVECTEE AU TEMPS tn POUR CET OBJET
C        ============================================
         DO J=1,NBCOBV
C           LA VALEUR CODE DE BLOCAGE DE VITESSE CONVECTEE EST EMPILEE
            BLVITE( J ) = 1D222
         ENDDO
C
      ELSE
C
C        FONCTION UTILISATEUR
C        ====================
         XYZ(1) = TEMPS
         XYZ(2) = XPI
         XYZ(3) = YPI
         XYZ(4) = ZPI
         XYZ(5) = NYOBJT
         XYZ(6) = NUOBJT
         DO J=1,NBCOBV
C           LE NUMERO DE LA COMPOSANTE
            XYZ(7) = MCN( MNBLVI + WUCOBV + J - 1 )
            CALL FONVAL( MCN(MN), 7, XYZ,
     %                   NCODEV, BLVITE(J) )
         ENDDO
      ENDIF
C
      RETURN
      END
