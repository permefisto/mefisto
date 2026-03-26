      SUBROUTINE DLCPND( NBCOMP, NBNOEU, VCP, VND )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT:     REORDONNER LES DL DONNES PAR COMPOSANTES
C ----     EN LES DL PAR NOEUDS
C
C ENTREES:
C --------
C NBCOMP : NOMBRE DE COMPOSANTES
C NBNOEU : NOMBRE DE NOEUDS
C VCP    : VALEUR DONNEE PAR COMPOSANTES
C
C SORTIES:
C --------
C VND    : VALEUR DONNEE PAR NOEUDS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Texas A & M University at QATAR  Fevrier 2011
C23456---------------------------------------------------------------012
      DOUBLE PRECISION  VCP(NBNOEU,NBCOMP), VND(NBCOMP,NBNOEU)
C
      DO K=1,NBCOMP
         DO N=1,NBNOEU
            VND( K, N ) = VCP( N, K )
         ENDDO
      ENDDO
C
      RETURN
      END
