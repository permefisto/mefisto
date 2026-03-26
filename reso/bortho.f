      SUBROUTINE BORTHO( N, NBVP, VECPR, BVECPR,  V )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    B ORTHOGONALISER LE VECTEUR V AUX NBVP VECTEURS PROPRES
C -----
C ENTREES:
C --------
C N      : NOMBRE DE COMPOSANTES DES VECTEURS
C NBVP   : NOMBRE DE VECTEURS PROPRES
C VECPR  : NBVP VECTEURS B-ORTHOGONAUX
C BVECPR : NBVP B * VECPR OU B EST LA MATRICE DE MASSE
C
C MODIFIE:
C --------
C V      : LE VECTEUR B-ORTHOGONAL AUX NBVP VECPR
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET UPMC ANALYSE NUMERIQUE PARIS     DECEMBRE 1998
C2345X7..............................................................012
      DOUBLE PRECISION V(N), VECPR(N,NBVP), BVECPR(N,NBVP), tVBVk
C
      DO K=1,NBVP
C
         tVBVk = 0D0
         DO I=1,N
            tVBVk = tVBVk + V(I) * BVECPR(I,K)
         ENDDO
C
         DO I=1,N
            V(I) = V(I) - tVBVk * VECPR(I,K)
         ENDDO
C
      ENDDO
C
      RETURN
      END
