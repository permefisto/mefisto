       SUBROUTINE ITEML3( XYZ , NOMOBJ , NUOBJT )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : SAUVEGARDER DANS LE TABLEAU DES ITEMS-LIGNES UNE LIGNE VISIBLE
C ----- SUR L'ECRAN A PARTIR DE SES COORDONNEES OBJET EN 3D
C       EN SORTIE LA PLUME DE TRACE EST POSITIONNEE SUR L'ITEM
C
C ENTREES :
C ---------
C XYZ    : COORDONNEES 3D OBJET DE L'ITEM-LIGNE
C NOMOBJ : NOM DE LA LIGNE DANS SON LEXIQUE
C NUOBJT : NUMERO DE LA LIGNE DANS SON LEXIQUE
C
C SORTIE  :
C ---------
C LE TABLEAU MNITEL EST AUGMENTE EVENTUELLEMENT DE L'ITEM
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC PARIS       AVRIL 1990
C.......................................................................
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      REAL              XYZ(3),XYZA(3)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      CHARACTER*(*)     NOMOBJ
      CHARACTER*48      CARSYM
C
C     COORDONNEES ECRAN DE L'ITEM
      CALL XYZAXO( XYZ, XYZA )
      NX = NUPXEX( XYZA(1) )
      NY = NUPXEY( XYZA(2) )
C
      IF( NX .GE. 0 .AND. NX .LE. LAPXFE  .AND.
     %    NY .GE. 0 .AND. NY .LE. LHPXFE       ) THEN
C
C         L'ITEM VISIBLE EST MIS DANS LE TABLEAU MNITEL
          IF( MCN(MNITEL+2) .GE. MCN(MNITEL+1) ) THEN
C            TABLEAU TROP PETIT:LA TAILLE DU TABLEAU EST AUGMENTEE
             CALL ITEMAU( MNITEL )
          ENDIF
C
C         LE NOMBRE D'ITEM EST AUGMENT DE 1
          MCN(MNITEL+2) = MCN(MNITEL+2) + 1
C
C         L'ADRESSE MCN DE L'ITEM
          MNI = MNITEL + MCN(MNITEL) * MCN(MNITEL+2)
C
C         LES COORDONNEES IMAGE DE CET ITEM , SON CODE
          MCN( MNI     ) = NX
          MCN( MNI + 1 ) = NY
C
C         NUMERO DE L'OBJET DANS SON LEXIQUE
          MCN( MNI + 2 ) = NUOBJT
C
C        TRACE EFFECTIF DE L'ITEM SELON LE SYMBOLE '#'
C        TRACE CONDITIONNEL DU NOM DE LA LIGNE
         IF( LPLIGN .NE. 0 ) THEN
            CARSYM = '#' // NOMOBJ
C           LE NOMBRE DE CARACTERES DE CARSYM
            NB = NUDCNB( CARSYM )
C           TRACE DE LA CHAINE = SYMBOLE + EVENTUEL NOM
            CALL SYMBOLE2D( NCOLIG, XYZA(1), XYZA(2), CARSYM(1:NB) )
         ENDIF
      ENDIF
      END
