       SUBROUTINE ITEMS3( XYZ , NOMOBJ , NUOBJT )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : SAUVEGARDER DANS LE TABLEAU DES ITEMS-SURFACES UNE SURFACE VISIBLE
C ----- SUR L'ECRAN A PARTIR DE SES COORDONNEES OBJET EN 3D
C       EN SORTIE LA PLUME DE TRACE EST POSITIONNEE SUR L'ITEM

C ENTREES :
C ---------
C XYZ    : COORDONNEES 3D OBJET DE L'ITEM-SURFACE
C NOMOBJ : NOM DE LA SURFACE DANS SON LEXIQUE
C NUOBJT : NUMERO DE LA SURFACE DANS SON LEXIQUE

C SORTIE  :
C ---------
C LE TABLEAU MNITES EST AUGMENTE EVENTUELLEMENT DE L'ITEM
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC PARIS         MAI 1994
C.......................................................................
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      REAL              XYZ(3),XYZA(3)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      CHARACTER*(*)     NOMOBJ
      CHARACTER*48      CARSYM

C     COORDONNEES ECRAN DE L'ITEM
      CALL XYZAXO( XYZ, XYZA )
      NX = NUPXEX( XYZA(1) )
      NY = NUPXEY( XYZA(2) )

      IF( NX .GE. 0 .AND. NX .LE. LAPXFE  .AND.
     %    NY .GE. 0 .AND. NY .LE. LHPXFE       ) THEN

C        L'ITEM VISIBLE EST MIS DANS LE TABLEAU MNITES
         IF( MCN(MNITES+2) .GE. MCN(MNITES+1) ) THEN
C           TABLEAU TROP PETIT:LA TAILLE DU TABLEAU EST AUGMENTEE
            CALL ITEMAU( MNITES )
         ENDIF

C        LE NOMBRE D'ITEM EST AUGMENTE DE 1
         MCN(MNITES+2) = MCN(MNITES+2) + 1

C        L'ADRESSE MCN DE L'ITEM
         MNI = MNITES + MCN(MNITES) * MCN(MNITES+2)

C        LES COORDONNEES IMAGE DE CET ITEM, SON CODE
         MCN( MNI     ) = NX
         MCN( MNI + 1 ) = NY

C        NUMERO DE L'OBJET DANS SON LEXIQUE
         MCN( MNI + 2 ) = NUOBJT

C        TRACE EFFECTIF DE L'ITEM SELON LE SYMBOLE '.'
C        TRACE CONDITIONNEL DU NOM DE LA LIGNE
         IF( LPSURF .NE. 0 ) THEN
            CARSYM = '.' // NOMOBJ
C           LE NOMBRE DE CARACTERES DE CARSYM
            NB = NUDCNB( NOMOBJ ) + 1
C           TRACE DE LA CHAINE = SYMBOLE + EVENTUEL NOM
            CALL SYMBOLE2D( NCOSUR, XYZA(1), XYZA(2), CARSYM(1:NB) )
         ENDIF
      ENDIF

      RETURN
      END
