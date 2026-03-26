       SUBROUTINE ITEMV3( XYZ , NOMOBJ , NUOBJT )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : SAUVEGARDER DANS LE TABLEAU DES ITEMS-VOLUMES UN VOLUME VISIBLE
C ----- SUR L'ECRAN A PARTIR DE SES COORDONNEES OBJET EN 3D
C       EN SORTIE LA PLUME DE TRACE EST POSITIONNEE SUR L'ITEM
C
C ENTREES :
C ---------
C XYZ    : COORDONNEES 3D OBJET DE L'ITEM-VOLUME
C NOMOBJ : NOM DU VOLUME DANS SON LEXIQUE
C NUOBJT : NUMERO DU VOLUME DANS SON LEXIQUE
C
C SORTIE  :
C ---------
C LE TABLEAU MNITEV EST AUGMENTE EVENTUELLEMENT DE L'ITEM
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC PARIS         MAI 1994
C.......................................................................
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      REAL              XYZ(3)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      CHARACTER*(*)     NOMOBJ
      CHARACTER*48      CARSYM
      REAL              XYZA(3)
C
C     COORDONNEES ECRAN DE L'ITEM
      CALL XYZAXO( XYZ, XYZA )
      NX = NUPXEX( XYZA(1) )
      NY = NUPXEY( XYZA(2) )
C
      IF( NX .GE. 0 .AND. NX .LE. LAPXFE  .AND.
     %    NY .GE. 0 .AND. NY .LE. LHPXFE       ) THEN
C
C        L'ITEM VISIBLE EST MIS DANS LE TABLEAU MNITEV
         IF( MCN(MNITEV+2) .GE. MCN(MNITEV+1) ) THEN
C           TABLEAU TROP PETIT:LA TAILLE DU TABLEAU EST AUGMENTEE
            CALL ITEMAU( MNITEV )
         ENDIF
C
C        LE NOMBRE D'ITEM EST AUGMENT DE 1
         MCN(MNITEV+2) = MCN(MNITEV+2) + 1
C        L'ADRESSE MCN DE L'ITEM
         MNI = MNITEV + MCN(MNITEV) * MCN(MNITEV+2)
C
C        LES COORDONNEES IMAGE DE CET ITEM , SON CODE
         MCN( MNI     ) = NX
         MCN( MNI + 1 ) = NY
C
C        NUMERO DE L'OBJET DANS SON LEXIQUE
         MCN( MNI + 2 ) = NUOBJT
C
C        TRACE EFFECTIF DE L'ITEM SELON LE SYMBOLE '*'
C        TRACE CONDITIONNEL DU NOM DU POINT
         IF( LPVOLU .NE. 0 ) THEN
            CARSYM = '*' // NOMOBJ
C           LE NOMBRE DE CARACTERES DE CARSYM
            NB = NUDCNB( NOMOBJ ) + 1
C           TRACE DE LA CHAINE = SYMBOLE + EVENTUEL NOM
            CALL SYMBOLE2D( NCOVOL, XYZA(1), XYZA(2), CARSYM(1:NB) )
         ENDIF
      ENDIF
      END
