       SUBROUTINE ITEMST2( NOST, XYZSOM )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    SAUVEGARDER DANS LE TABLEAU DES ITEMS-SOMMETS UN SOMMET VISIBLE
C -----    SUR L'ECRAN A PARTIR DE SES COORDONNEES OBJET EN 2D

C ENTREES:
C --------
C NOST     : NUMERO DU SOMMET DANS LE TABLEAU XYZSOM DU MAILLAGE
C XYZSOM : XYZ DES SOMMETS DU MAILLAGE

C SORTIE :
C --------
C LE TABLEAU MNITST (cf mecoit.inc) EST AUGMENTE EVENTUELLEMENT DE L'ITEM
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  Saint Pierre du Perray            Avril 2020
C.......................................................................
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/pp.inc"
      COMMON         MCN(MOTMCN)
      REAL           XYZSOM(3,*)

C     COORDONNEES ECRAN DE L'ITEM
      NX = NUPXEX( XYZSOM(1,NOST) )
      NY = NUPXEY( XYZSOM(2,NOST) )

      IF( NX .GE. 0 .AND. NX .LE. LAPXFE   .AND.
     %    NY .GE. 0 .AND. NY .LE. LHPXFE ) THEN

C        L'ITEM SOMMET VISIBLE EST MIS DANS LE TABLEAU MNITST
C        ----------------------------------------------------
         IF( MCN(MNITST+2) .GE. MCN(MNITST+1) ) THEN
C           TABLEAU TROP PETIT:LA TAILLE DU TABLEAU EST AUGMENTEE
            CALL ITEMAU( MNITST )
         ENDIF

C        LE NOMBRE D'ITEMS EST AUGMENTE DE 1
         MCN(MNITST+2) = MCN(MNITST+2) + 1

C        L'ADRESSE MCN DE L'ITEM
         MNI = MNITST + MCN(MNITST) * MCN(MNITST+2)

C        LES COORDONNEES IMAGE DE CET ITEM, SON CODE
         MCN( MNI     ) = NX
         MCN( MNI + 1 ) = NY

C        NUMERO DE SOMMET DANS LE MAILLAGE XYZSOM
         MCN( MNI + 2 ) = NOST

C        TRACE EN ROUGE DE L'ITEM SOMMET '+NOST'
         CALL TRST2D( NCROUG, NOST, XYZSOM )

      ENDIF

      RETURN
      END
