      SUBROUTINE GAR2P1BP1( NBSOM, NUELEM, NONOEF, AG, NBNOVI, BG )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    ELIMINATION DE GAUSS SUR LES 2 COMPOSANTES DU SECOND MEMBRE
C -----    POUR LES DL BARYCENTRES DES TRIANGLES BREZZI-FORTIN
C
C ENTREES:
C --------
C NBSOM  : NOMBRE DE SOMMETS DE LA TRIANGULATION
C NUELEM : NUMERO DU TRIANGLE BREZZI-FORTIN A TRAITER
C NONOEF : NO DES 3 SOMMETS DU TRIANGLE NUELEM
C AG     : LA LIGNE DE AE DU 4-EME DEGRE DE LIBERTE AU BARYCENTRE
C
C MODIFIE:
C --------
C BG     : VECTEUR SECOND MEMBRE AVANT ET APRES GAUSS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC Saint Pierre du Perray Decembre 2008
C23456---------------------------------------------------------------012
      DOUBLE PRECISION  AG(4), BG(NBNOVI,2), PIVOT
      INTEGER           NONOEF(3)
C
C     TRIANGULATION DE GAUSS SUR LA LIGNE 4 DU VECTEUR
      DO K = 1, 2
         PIVOT = BG( NBSOM+NUELEM, K ) / AG(4)
         DO I = 1, 3
            N = NONOEF( I )
            BG(N,K) = BG(N,K) - AG(I) * PIVOT
         ENDDO
      ENDDO
C
      RETURN
      END
