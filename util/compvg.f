       SUBROUTINE COMPVg( NTDL, NODLIB, VG )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT:     COMPRESSION PAR SUPPRESSION DES DL FIXES DEFINIS PAR NODLIB
C ----
C ENTREES:
C --------
C NTDL   : NOMBRE TOTAL DE DEGRES DE LIBERTE
C NODLIB : NODLIB(I) = NUMERO DU DEGRE DE LIBERTE S'IL EST LIBRE
C                     -INDICE DANS LA LISTE DES DL BLOQUES S'IL EST BLOQUE
C VG     : VECTEUR GLOBAL DE 1 A NTDL
C
C SORTIE :
C --------
C Vg     : VECTEUR GLOBAL SANS LES COMPOSANTES DES DL FIXES
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET St PIERRE DU PERRAY & LJLL UPMC  JANVIER 2010
C23456---------------------------------------------------------------012
      DOUBLE PRECISION  VG(1:NTDL)
      INTEGER           NODLIB(NTDL)
C
      DO K = 1, NTDL
C
C        NUMERO DU DL DANS LA NUMEROTATION DES DEGRES DE LIBERTE NON FIXES
C        ATTENTION: LES DL LIBRES ONT UN NUMERO DL LIBRE NODLIB(K)<K
C                   AUX DL AVANT ET LA NUMEROTATION EST CROISSANTE
C                   CF L'ORDRE DE NUMEROTATION DANS redlib.f
         NDL = NODLIB(K)
         IF( NDL .GT. 0 ) THEN
            VG(NDL) = VG(K)
         ENDIF
C
      ENDDO
C
      RETURN
      END
