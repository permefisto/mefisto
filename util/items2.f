       SUBROUTINE ITEMS2( XY , NOMOBJ , NUOBJT )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : SAUVEGARDER DANS LE TABLEAU DES ITEMS-SURFACES UNE SURFACE VISIBLE
C ----- SUR L'ECRAN A PARTIR DE SES COORDONNEES OBJET EN 2D
C       EN SORTIE LA PLUME DE TRACE EST POSITIONNEE SUR L'ITEM
C
C ENTREES :
C ---------
C XY     : COORDONNEES 2D OBJET DE L'ITEM-SURFACE
C NOMOBJ : NOM DE LA SURFACE DANS SON LEXIQUE
C NUOBJT : NUMERO DE LA SURFACE DANS SON LEXIQUE
C
C SORTIE  :
C ---------
C LE TABLEAU MNITES EST AUGMENTE EVENTUELLEMENT DE L'ITEM
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC PARIS         MAI 1994
C.......................................................................
      PARAMETER        ( NUTYOB = 3 )
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      REAL              XY(2)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      CHARACTER*(*)     NOMOBJ
      CHARACTER*48      CARSYM
C
C     COORDONNEES ECRAN DE L'ITEM
      NX = NUPXEX( XY(1) )
      NY = NUPXEY( XY(2) )
C
      IF( NX .GE. 0 .AND. NX .LE. LAPXFE .AND.
     %    NY .GE. 0 .AND. NY .LE. LHPXFE ) THEN
C
C        L'ITEM VISIBLE EST MIS DANS LE TABLEAU MNITES
         IF( MCN(MNITES+2) .GE. MCN(MNITES+1) ) THEN
C           TABLEAU TROP PETIT:LA TAILLE DU TABLEAU EST AUGMENTEE
            CALL ITEMAU( MNITES )
         ENDIF
C
C        LE NOMBRE D'ITEM EST AUGMENT DE 1
         MCN(MNITES+2) = MCN(MNITES+2) + 1
C
C        L'ADRESSE MCN DE L'ITEM
         MNI = MNITES + MCN(MNITES) * MCN(MNITES+2)
C
C        LES COORDONNEES IMAGE DE CET ITEM , SON CODE
         MCN( MNI     ) = NX
         MCN( MNI + 1 ) = NY
C
C        NUMERO DE L'OBJET DANS SON LEXIQUE
         MCN( MNI + 2 ) = NUOBJT
C
C        TRACE EFFECTIF DE L'ITEM SELON LE SYMBOLE '.'
C        TRACE CONDITIONNEL DU NOM DE LA LIGNE
         IF( LPSURF .NE. 0 ) THEN
            CARSYM = '.' // NOMOBJ
C           LE NOMBRE DE CARACTERES DE CARSYM
            NB = NUDCNB( NOMOBJ ) + 1
C           TRACE DE LA CHAINE = SYMBOLE + EVENTUEL NOM
            CALL SYMBOLE2D( NCOSUR, XY(1), XY(2), CARSYM(1:NB) )
         ENDIF
      ENDIF
      END
