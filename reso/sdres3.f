       SUBROUTINE SDRES3(RAPPORT,V1,V2)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCUL DE LA NOUVELLE DIRECTION DE DESCENTE
C -----

C ENTREES :
C ---------
C RAPPORT : PARAMETRE D'ORTHOGONALISATION
C V1,V2   : LES VECTEURS
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY     ANALYSE NUMERIQUE UPMC PARIS  MAI 1990
C23456---------------------------------------------------------------012
      COMMON / SOUSDO / NORESO,NTDL,NDSM,NDIM,NPIMAX
      DOUBLE PRECISION  RAPPORT,V1(NTDL,NDSM),V2(NTDL,NDSM)

      DO I = 1 , NTDL
         DO NS = 1 , NDSM
            V1(I,NS) = V2(I,NS) + RAPPORT * V1(I,NS)
         ENDDO
      ENDDO

      RETURN
      END
