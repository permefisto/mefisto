      SUBROUTINE UUMGRP( NDIM, NBNOVI, VG0, MDIAG, CGRADP,  VG1 )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   MDIAG VG1 = MDIAG VG0 + C GRAD P => VG1
C -----
C
C ENTREES:
C --------
C NDIM   : NOMBRE DE COLONNES DU VECTEUR GLOBAL VG
C NBNOVI : NOMBRE DE DL D'UNE COMPOSANTE DU VECTEUR VG
C VG0    : NDIM COMPOSANTES MODIFIEES DU VECTEUR GLOBAL
C MDIAG  : MDIAG(NBNOVI) LA MATRICE DIAGONALE
C
C SORTIES:
C --------
C VG1    : NDIM COMPOSANTES MODIFIEES DU VECTEUR GLOBAL
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC & St Pierre du Perray  Decembre 2010
C23456---------------------------------------------------------------012
      DOUBLE PRECISION  VG0(NBNOVI,NDIM), VG1(NBNOVI,NDIM),
     %                  CGRADP(NBNOVI,NDIM), MDIAG(NBNOVI)
C
      DO NS=1,NBNOVI
         DO K=1,NDIM
            VG1( NS, K ) = VG0( NS, K ) + CGRADP( NS, K ) / MDIAG( NS )
         ENDDO
      ENDDO
C
      RETURN
      END
