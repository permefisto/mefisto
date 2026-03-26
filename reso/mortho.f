      SUBROUTINE MORTHO( N,      NBVP, VECPR, MGVECPR,
     %                   NCODSM, MUMG, MG,    V, W )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    MG ORTHONORMALISER LE VECTEUR V AUX NBVP VECTEURS PROPRES
C -----
C ENTREES:
C --------
C N      : NOMBRE DE COMPOSANTES DES VECTEURS
C NBVP   : NOMBRE DE VECTEURS PROPRES
C VECPR  : NBVP VECTEURS M-ORTHOGONAUX
C MGVECPR: NBVP MG * VECPR OU MG EST LA MATRICE DE MASSE
C
C NCODSM : CODE DE STOCKAGE DE LA MATRICE GLOBALE DE CAPACITE
C          1 SI MATRICE SYMETRIQUE
C MUMG   : TABLEAU DES POINTEURS SUR LES COEFFICIENTS DIAGONAUX DE MG
C MG     : MATRICE DE MASSE PROFIL
C
C MODIFIE:
C --------
C V      : LE VECTEUR MG-ORTHONORME AUX NBVP VECPR
C W      : LE VECTEUR MG V
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET TEXAS A & M UNIVERSITY at QATAR     MARS 2011
C2345X7..............................................................012
      DOUBLE PRECISION  MG(*), V(N), W(N), VECPR(N,NBVP),MGVECPR(N,NBVP)
      INTEGER           MUMG(0:N)
      DOUBLE PRECISION  tVMgVk
      INTRINSIC         SQRT
C
C     Mg-ORTHOGONALISATION
      DO K=1,NBVP
C
         tVMgVk = 0D0
         DO I=1,N
            tVMgVk = tVMgVk + V(I) * MGVECPR(I,K)
         ENDDO
C
         DO I=1,N
            V(I) = V(I) - tVMgVk * VECPR(I,K)
         ENDDO
C
      ENDDO
C
C     Wk = MG Vk
      CALL MAPRVE( 0, 1D0, N, NCODSM, MUMG, MG, V,  W )
C
C     MG ORTHONORMALISATION DE V et W
      tVMgVk = 0D0
      DO I=1,N
         tVMgVk = tVMgVk + V(I) * W(I)
      ENDDO
      tVMgVk = SQRT( tVMgVk )

      DO I=1,N
         V(I) = V(I) / tVMgVk
         W(I) = W(I) / tVMgVk
      ENDDO
C
      RETURN
      END
