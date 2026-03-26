      SUBROUTINE TESTL2( N, U, V, EPSR, EPSA,
     %                   RESULT, VALNORM )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   ON VEUT TESTER LA NORME RELATIVE DE LA DIFFERENCE ENTRE
C -----   DEUX VECTEURS U ET V DE N COMPOSANTES.
C         SI LA ||NORME DE U|| < 0.0001
C         ALORS PAS DE DIVISION DE VALNORM PAR LA NORME DE U
C
C ENTREES :
C ---------
C N        : NOMBRE DE COMPOSANTES DU VECTEUR
C U , V    : VECTEURS A TESTER
C EPSR,EPSA: LES SEUILS DE TEST EN ABSOLU ET RELATIF
C
C SORTIE  :
C ---------
C RESULT  : LE RESULTAT DU TEST
C           =1 SI ||U-V|| SOUS LE SEUIL
C           =0 SINON
C VALNORM : ||U-V|| EN RELATIF OU ABSOLU
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEURS: CHANRION OLIVIER & JEAN-MICHEL ROVARCH
C          DEA D'ANALYSE NUMERIQUE PARIS  JANVIER 1998
C ......................................................................
      IMPLICIT NONE
      INTEGER           N, RESULT, I
      DOUBLE PRECISION  U(N), V( N )
      DOUBLE PRECISION  EPSA, EPSR, NORMU, NORMDF, VALNORM
C
C     LA NORME L2 DU VECTEUR U ET LA NORME L2 DU VECTEUR U-V
      NORMU  = 0D0
      NORMDF = 0D0
      DO 10 I=1,N
         NORMU  = NORMU  + U( I ) ** 2
         NORMDF = NORMDF + ( U( I ) - V( I ) ) ** 2
   10 CONTINUE
C
C     TEST ABSOLU OU RELATIF EN FONCTION DE LA NORME DE U
C     ---------------------------------------------------
      IF( NORMU .LE. 1D-12 ) THEN
C
C        TEST EN ABSOLU
         IF ( NORMDF .LE. EPSA*EPSA ) THEN
            RESULT=1
         ELSE
            RESULT=0
         ENDIF
         VALNORM=NORMDF
C
      ELSE
C
C        TEST EN RELATIF
         VALNORM = NORMDF / NORMU
         IF ( VALNORM .LE. EPSR*EPSR ) THEN
            RESULT=1
         ELSE
            RESULT=0
         ENDIF
      ENDIF
C
      RETURN
      END
