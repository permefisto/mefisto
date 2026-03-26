      SUBROUTINE TGSIDE( XYZTG1 , XYZTG2 , LSIGNE )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    IDENTIFIER OU NOM 2 TANGENTES DEFINIES PAR 3 COORDONNEES
C -----
C
C ENTREES :
C ---------
C XYZTG1   : LA PREMIERE TANGENTE
C XYZTG2   : LA SECONDE  TANGENTE
C
C SORTIE :
C --------
C LSIGNE : +1 SI LES 2 TANGENTES SONT JUGEES IDENTIQUES
C          -1 SI LES 2 TANGENTES SONT JUGEES OPPOSEES
C           0 SI LES 2 TANGENTES SONT JUGEES DIFFERENTES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    DECEMBRE 1996
C23456---------------------------------------------------------------012
      COMMON / EPSSSS / EPZERO,EPSXYZ
      REAL              XYZTG1(3),XYZTG2(3)
C
      LSIGNE = 0
C
C     CARRE DE LA NORME DE LA TANGENTE
      TG1 = XYZTG1(1)**2 + XYZTG1(2)**2 + XYZTG1(3)**2
      TG2 = XYZTG2(1)**2 + XYZTG2(2)**2 + XYZTG2(3)**2
C
      IF( ABS(TG1-TG2) .GT. TG2*EPSXYZ ) RETURN
C
C     IDENTIFICATION POSSIBLE AU SIGNE PRES
      EPSTG = SQRT( TG2 ) * EPSXYZ
C
C     BOUCLE SUR LES 3 COMPOSANTES
      DO 20 K=1,3
         R1 = ABS( XYZTG1( K ) )
         R2 = ABS( XYZTG2( K ) )
         IF( ABS(R1-R2) .GT. EPSTG ) GOTO 30
C
C        LES VALEURS ABSOLUES SONT EGALES. LES SIGNES?
C        ATTENTION AUX VALEURS PRESQUE NULLES DUES
C        AUX ERREURS D'ARRONDIS
         IF( R1 .GT. EPSTG ) THEN
            IF( LSIGNE .EQ. 0 ) THEN
C              PREMIERE COMPARAISON
               IF( XYZTG1(K)*XYZTG2(K) .LT. 0 ) THEN
                  LSIGNE = -1
               ELSE
                  LSIGNE =  1
               ENDIF
            ELSE
C              AU MOINS SECONDE COMPARAISON
               IF( XYZTG1(K)*XYZTG2(K) .LT. 0 ) THEN
C                 LE SIGNE RECENSE DOIT ETRE -
C                 SINON PAS D'IDENTIFICATION
                  IF( LSIGNE .GT. 0 ) GOTO 30
               ELSE
C                 LE SIGNE RECENSE DOIT ETRE +
C                 SINON PAS D'IDENTIFICATION
                  IF( LSIGNE .LT. 0 ) GOTO 30
               ENDIF
            ENDIF
         ENDIF
 20   CONTINUE
C
C     LES 2 TANGENTES SONT IDENTIFIEES AU SIGNE LSIGNE PRES
      RETURN
C
C     LES TANGENTES SONT DIFFERENTES
 30   LSIGNE = 0
      END
