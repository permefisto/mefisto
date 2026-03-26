      SUBROUTINE ASB0D( N, A, U,  V )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DE V = A * U
C -----    	
C
C ENTREES:
C --------
C N      : NOMBRE DE LIGNES DES VECTEURS ET DE LA MATRICE SYMETRIQUE A
C A      : MATRICE PLEINE SYMETRIQUE
C U      : VECTEUR
C
C SORTIE :
C --------
C V      : VECTEUR RESULTAT
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : MIGUEL FERNANDEZ VARELA ANALYSE NUMERIQUE UPMC   JANVIER 1998
C-----------------------------------------------------------------------
      DOUBLE PRECISION A(N*(N+1)/2), U(N), V(N)
C
C     LA BOUCLE SUR LES LIGNES	
      IG=0
      DO 100 I=1,N
C
          V(I)=0.D0
C
C         SOMMATION SUR LES COLONNES
          DO 10 J=1,I
             IG=IG+1
             V(I)=V(I)+A(IG)*U(J)
 10       CONTINUE
C
C         SOMMATION SUR LES LIGNES (PAR SYMETRIE)
          DO 20 J=I+1,N
             IA=J*(J-1)/2+I
             V(I)=V(I)+A(IA)*U(J)
 20       CONTINUE
C
 100  CONTINUE
      RETURN
      END
