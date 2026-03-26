      SUBROUTINE ORDXYZ( XYZA , NPA , IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : REORDONNER LES 2 POINTS SELON LES X CROISSANTS
C ----- SI EGALITE SELON LES Y CROISSANTS
C       SI EGALITE SELON LES Z CROISSANTS
C       SI EGALITE DES 2 POINTS IERR = 1 EN SORTIE
C
C ENTREES ET SORTIES:
C -------------------
C XYZA   : COORDONNEES DES 2 POINTS DE R ** 3
C          XYZA(.,I) LE I-EME POINT AVEC I =1 OU 2
C NPA    : NUMERO DES 2 POINTS DANS UNE CERTAINE NUMEROTATION
C          EN SORTIE LE 1-ER A LE PLUS PETIT X ...
C
C SORTIE :
C --------
C IERR   : 0 SI PAS D'ERREUR
C          1 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : ALAIN PERRONNET ANALYSE NUMRIQUE PARIS  JANVIER 1986
C...............................................................................
      INTEGER NPA(2)
      REAL    XYZA(3,2)
C
      IERR = 0
      DO 10 I=1,3
         IF( XYZA(I,1) .NE. XYZA(I,2) ) THEN
            IF( XYZA(I,1) .GT. XYZA(I,2) ) GOTO 20
            RETURN
         ENDIF
 10   CONTINUE
C
C     ERREUR LES 2 POINTS SONT CONFONDUS
      IERR = 1
      RETURN
C
C     PERMUTATION DES 2 POINTS
 20   I      = NPA(1)
      NPA(1) = NPA(2)
      NPA(2) = I
      DO 30 I=1,3
         X         = XYZA(I,1)
         XYZA(I,1) = XYZA(I,2)
         XYZA(I,2) = X
 30   CONTINUE
      END
