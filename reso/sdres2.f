       SUBROUTINE SDRES2(RHO,U,Z,F,P,R,RL)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : MISE A JOUR DES VECTEURS ITERES DE L'ALGORITHME
C -----

C ENTREES :
C ---------
C RHO  : PARAMETRE DE DESCENTE
C U    : VECTEUR A MODIFIER
C Z    : VECTEUR A MODIFIER
C F    : VECTEUR A MODIFIER
C P    : VECTEUR DIRECTION
C R    : VECTEUR DIRECTION
C RL   : VECTEUR DIRECTION
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY     ANALYSE NUMERIQUE UPMC PARIS  MAI 1990
C23456---------------------------------------------------------------012
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      COMMON / SOUSDO / NORESO,NTDL,NDSM,NDIM,NPIMAX
      DOUBLE PRECISION  RHO,U(NTDL,NDSM),Z(NTDL,NDSM),F(NTDL,NDSM),
     &                      P(NTDL,NDSM),R(NTDL,NDSM),RL(NTDL,NDSM)

      DO I  = 1 , NTDL
         DO NS = 1 , NDSM
            U(I,NS) = U(I,NS) + RHO * Z(I,NS)
            F(I,NS) = F(I,NS) + RHO * P(I,NS)
            R(I,NS) = R(I,NS) + RHO * RL(I,NS)
         ENDDO
      ENDDO

      RETURN
      END
