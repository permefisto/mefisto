       SUBROUTINE tVgAeVg( N, NODL, VG, AE, RESULT )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT:    CALCUL DU PRODUIT TVECTEUR MATRICE VECTEUR MATRICIEL DANS LE
C ----    CAS OU LA MATRICE EST SYMETRIQUE STOCKEE SOUS FORME TRIANGULAIRE
C         RESULT = tVE AE VE REEL DOUBLE PRECISION
C
C ENTREES:
C --------
C N      : NOMBRE DE LIGNES DE V et DE LIGNES ET COLONNES DE AE
C NODL   : NUMERO GLOBAL DES N DEGRES DE LIBERTE LOCAUX
C VG     : VECTEUR GLOBAL DE 1 A NTDL
C          VE = VG( NODL() )
C AE     : MATRICE (N,N) SYMETRIQUE STOCKEE DE HAUT EN BAS ET
C          DE LA GAUCHE VERS LA DROITE.
C          SEULE LA PARTIE TRIANGULAIRE SUPERIEURE EST STOCKEE
C
C SORTIE :
C --------
C RESULT : = tVE AE VE REEL DOUBLE PRECISION
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET St PIERRE DU PERRAY & LJLL UPMC  JANVIER 2010
C23456---------------------------------------------------------------012
      DOUBLE PRECISION  VG(1:*), AE(N*(N+1)/2), RESULT
      INTEGER           NODL(N)
C
      K = 0
      RESULT = 0D0
      DO I = 1, N
         NODLI = NODL(I)
         DO J = 1, I
            K = K + 1
            RESULT = RESULT + VG( NODLI ) * AE(K) * VG( NODL(J) )
         ENDDO
         KT = K
         DO J = I+1, N
            KT = KT + J - 1
            RESULT = RESULT + VG( NODLI ) * AE(KT) * VG( NODL(J) )
         ENDDO
      ENDDO
C
      RETURN
      END
