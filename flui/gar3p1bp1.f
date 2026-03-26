      SUBROUTINE GAR3P1BP1( NBSOM, NUELEM, NONOEF, AG, NBNOVI, BG )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    ELIMINATION DE GAUSS SUR LES 3 COMPOSANTES DU SECOND MEMBRE
C -----    POUR LE DL BARYCENTRE DU TETRAEDRE DE BREZZI-FORTIN
C
C ENTREES:
C --------
C NBSOM  : NOMBRE DE SOMMETS DE LA TETRAEDRISATION
C NUELEM : NUMERO DU TETRAEDRE BREZZI-FORTIN A TRAITER
C NONOEF : NO DES 4 SOMMETS DU TETRAEDRE NUELEM
C AG     : LIGNE DE AE DU 5-EME DEGRE DE LIBERTE AU BARYCENTRE
C
C MODIFIE:
C --------
C BG     : VECTEUR SECOND MEMBRE AVANT ET APRES GAUSS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC Saint Pierre du Perray Decembre 2010
C23456---------------------------------------------------------------012
      DOUBLE PRECISION  AG(5), BG(NBNOVI,*), PIVOT
      INTEGER           NONOEF(4)
C
C     TRIANGULATION DE GAUSS SUR LA LIGNE 5 DU VECTEUR
      DO K = 1, 3
         PIVOT = BG( NBSOM+NUELEM, K ) / AG(5)
         DO I = 1, 4
            N = NONOEF( I )
            BG(N,K) = BG(N,K) - AG(I) * PIVOT
         ENDDO
      ENDDO
C
      RETURN
      END
