       SUBROUTINE ITEMP2( XY , NOMOBJ , NUOBJT )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : SAUVEGARDER DANS LE TABLEAU DES ITEMS-POINTS UN POINT VISIBLE
C ----- SUR L'ECRAN A PARTIR DE SES COORDONNEES OBJET EN 2D
C       EN SORTIE LA PLUME DE TRACE EST POSITIONNEE SUR L'ITEM

C ENTREES:
C --------
C XY     : COORDONNEES 2D OBJET DE L'ITEM-POINT
C NOMOBJ : NOM DU POINT DANS SON LEXIQUE
C NUOBJT : NUMERO DU POINT DANS SON LEXIQUE
C
C SORTIE :
C --------
C LE TABLEAU MNITEP EST AUGMENTE EVENTUELLEMENT DE L'ITEM
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC PARIS         MAI 1994
C.......................................................................
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      REAL              XY(2)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      CHARACTER*(*)     NOMOBJ
      CHARACTER*48      CARSYM

C     COORDONNEES ECRAN DE L'ITEM
      NX = NUPXEX( XY(1) )
      NY = NUPXEY( XY(2) )

      IF( NX .GE. 0 .AND. NX .LE. LAPXFE .AND.
     %    NY .GE. 0 .AND. NY .LE. LHPXFE ) THEN

C        L'ITEM VISIBLE EST MIS DANS LE TABLEAU MNITEP
         IF( MCN(MNITEP+2) .GE. MCN(MNITEP+1) ) THEN
C           TABLEAU TROP PETIT:LA TAILLE DU TABLEAU EST AUGMENTEE
            CALL ITEMAU( MNITEP )
         ENDIF

C        LE NOMBRE D'ITEM EST AUGMENTE DE 1
         MCN(MNITEP+2) = MCN(MNITEP+2) + 1
C        L'ADRESSE MCN DE L'ITEM
         MNI    = MNITEP + MCN(MNITEP) * MCN(MNITEP+2)
C        LES COORDONNEES IMAGE DE CET ITEM , SON CODE
         MCN( MNI     ) = NX
         MCN( MNI + 1 ) = NY
C        NUMERO DE L'OBJET DANS SON LEXIQUE
         MCN( MNI + 2 ) = NUOBJT

C        TRACE EFFECTIF DE L'ITEM SELON LE SYMBOLE '+'
         CARSYM = '+'
         NB     = 1

C        TRACE CONDITIONNEL DU NOM DU POINT
         IF( LPPOIN .NE. 0 ) THEN
            CARSYM = '+' // NOMOBJ
C           LE NOMBRE DE CARACTERES DE CARSYM
            NB = NUDCNB( NOMOBJ ) + 1
         ENDIF

C        TRACE DE LA CHAINE = SYMBOLE + EVENTUEL NOM
         CALL SYMBOLE2D( NCOPOI, XY(1), XY(2), CARSYM(1:NB) )
      ENDIF

      RETURN
      END
