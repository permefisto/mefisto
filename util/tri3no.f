      SUBROUTINE TRI3NO( NAVANT , NAPRES )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRI CROISSANT DE 3 ENTIERS
C -----
C ENTREE :
C --------
C NAVANT : LES 3 ENTIERS A TRIER
C
C SORTIE :
C --------
C NAPRES : LES 3 ENTIERS APRES LE TRI
C
C REMARQUE : EN ENTREE NAVANT ET NAPRES PEUVENT ETRE CONFONDUS
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE UPMC PARIS    JANVIER 1987
C2345X7..............................................................012
      INTEGER  NAVANT(3), NAPRES(3)
C
C     PROTECTION DES 3 VALEURS D'ENTREE
      DO I=1,3
         NAPRES(I) = NAVANT(I)
      ENDDO
C
C     LE TRI
C
C     SI   NAPRES(1) > NAPRES(2) ILS SONT PERMUTES
      N = NAPRES(1)
      IF( N .GT. NAPRES(2) )  THEN
         NAPRES(1) = NAPRES(2)
         NAPRES(2) = N
      ENDIF
C     ICI NAPRES(1) < NAPRES(2)
C
C     SI  NAPRES(2) > NAPRES(3) ILS SONT PERMUTES
      N = NAPRES(2)
      IF( N .GT. NAPRES(3) ) THEN
         NAPRES(2) = NAPRES(3)
         NAPRES(3) = N
C        ICI NAPRES(2) < NAPRES(3)
C
C        SI  NAPRES(1) > NAPRES(2) ILS SONT PERMUTES
         N = NAPRES(1)
         IF( N .GT. NAPRES(2) )  THEN
            NAPRES(1) = NAPRES(2)
            NAPRES(2) = N
         ENDIF
      ENDIF
C
      RETURN
      END
