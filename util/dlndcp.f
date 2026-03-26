      SUBROUTINE DLNDCP( NBCOMP, NBNOEU, VND, VCP )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT:     REORDONNER LES DL DONNES PAR NOEUDS
C ----     EN LES DL PAR COMPOSANTES
C
C ENTREES:
C --------
C NBCOMP : NOMBRE DE COMPOSANTES
C NBNOEU : NOMBRE DE NOEUDS
C VND    : VALEUR DONNEE PAR NOEUDS
C
C SORTIES:
C --------
C VCP    : VALEUR DONNEE PAR COMPOSANTES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Texas A & M University at QATAR  Fevrier 2011
C23456---------------------------------------------------------------012
      DOUBLE PRECISION  VND(NBCOMP,NBNOEU), VCP(NBNOEU,NBCOMP)
C
      DO K=1,NBCOMP
         DO N=1,NBNOEU
            VCP( N, K ) = VND( K, N )
         ENDDO
      ENDDO
C
      RETURN
      END
