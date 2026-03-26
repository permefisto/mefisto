      SUBROUTINE CL2VED ( NBVAR, ALFA1, B01, ALFA2, B02, B )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    B <= ALFA1 * B01 + ALFA2 * B02
C -----    B , B01 , B02 VECTEURS DE NBVAR VARIABLES
C          ALFA1 , ALFA2 REELS DOUBLE PRECISION
C
C ENTREES:
C --------
C NBVAR  : NOMBRE DE VARIABLES DES VECTEURS B, B01, B02
C ALFA1  : PREMIER REEL DOUBLE PRECISION DE LA COMBINAISON LINEAIRE
C B01    : PREMIER VECTEUR REEL DOUBLE PRECISION DE LA COMBINAISON LINEAIRE
C ALFA2  : SECOND  REEL DOUBLE PRECISION DE LA COMBINAISON LINEAIRE
C B02    : SECOND  VECTEUR REEL DOUBLE PRECISION DE LA COMBINAISON LINEAIRE
C
C SORTIE :
C --------
C B      : VECTEUR RESULTAT DE LA COMBINAISON LINEAIRE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS    NOVEMBRE 1979
C2345X7..............................................................012
      DOUBLE PRECISION  B01(NBVAR),B02(NBVAR),B(NBVAR),ALFA1,ALFA2
C
C     TESTS SUIVANT LES VALEURS DE ALFA1 ET ALFA2
C     POUR REDUIRE LE NOMBRE DES OPERATIONS
C     =============================================
      IF( ALFA1 .EQ. 0.D0 ) THEN
C
         IF( ALFA2 .EQ. 1.D0 ) THEN
C           ALFA1=0.D0 ET ALFA2=1.D0
            DO 2 I=1,NBVAR
               B(I) = B02(I)
 2          CONTINUE
C
         ELSE
C
C           ALFA1=0.D0 ET ALFA2 QUELCONQUE
            DO 4 I=1,NBVAR
               B(I) = ALFA2 * B02(I)
 4          CONTINUE
C
         ENDIF
C
      ELSE IF( ALFA2 .EQ. 0.D0 ) THEN
C
         IF( ALFA1 .EQ. 1.D0 ) THEN
C           ALFA1=1.D0 ET ALFA2=0.D0
            DO 6 I=1,NBVAR
               B(I) = B01(I)
 6          CONTINUE
C
         ELSE
C
C           ALFA1 QUELCONQUE ET ALFA2=0.D0
            DO 8 I=1,NBVAR
               B(I) = ALFA1 * B01(I)
 8          CONTINUE
C
         ENDIF
C
      ELSE IF( ALFA1 .EQ. 1.D0 ) THEN
C
C        ALFA1=1.D0
         IF( ALFA2 .EQ. 1.D0 ) THEN
C
C           ALFA1=1.D0,ALFA2=1.D0
            DO 10 I=1,NBVAR
               B(I) = B01(I) + B02(I)
 10         CONTINUE
C
         ELSE
C
C           ALFA1=1.D0 ET ALFA2<>1.D0
            DO 20 I=1,NBVAR
               B(I) = B01(I) + ALFA2 * B02(I)
 20         CONTINUE
C
         ENDIF
C
      ELSE
C
C        ALFA1<>1D0
         IF( ALFA2 .EQ. 1.D0 ) THEN
C
C           ALFA1<>1.D0 ET ALFA2=1.D0
            DO 30 I=1,NBVAR
               B(I) = ALFA1 * B01(I) + B02(I)
 30         CONTINUE
C
         ELSE
C
C           ALFA1<>1.D0 ET ALFA2<>1.D0
            DO 40 I=1,NBVAR
               B(I) = ALFA1 * B01(I) + ALFA2 * B02(I)
 40         CONTINUE
C
         ENDIF
C
      ENDIF
      RETURN
      END
