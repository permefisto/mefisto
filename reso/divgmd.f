      SUBROUTINE DIVGMD( NDIM, NBNOVI, VG, MATDIAG )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    DIVISER CHAQUE COLONNE DU VECTEUR VG PAR LA MATRICE DIAGONALE
C -----
C
C ENTREES:
C --------
C NDIM   : NOMBRE DE COLONNES DU VECTEUR GLOBAL VG
C NBNOVI : NOMBRE DE DL D'UNE COMPOSANTE DU VECTEUR VG
C MATDIAG: MATDIAG(NBNOVI) LA MATRICE DIAGONALE
C
C SORTIES:
C --------
C VG     : NDIM COMPOSANTES MODIFIEES DU VECTEUR GLOBAL
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC & St Pierre du Perray  Decembre 2010
C23456---------------------------------------------------------------012
      DOUBLE PRECISION  VG(NBNOVI,NDIM), MATDIAG(NBNOVI)
C
      DO NS=1,NBNOVI
         DO K=1,NDIM
            VG( NS, K ) = VG( NS, K ) / MATDIAG( NS )
         ENDDO
      ENDDO
C
      RETURN
      END
