      SUBROUTINE TRI3N2( NAVAN1, NAVAN2, NAPRE1, NAPRE2 )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRI CROISSANT DES 3 ENTIERS  DE NAVAN1
C -----    ET MEMES PERMUTATIONS SUR NAVAN2
C
C ENTREE :
C --------
C NAVAN1 : LES 3 ENTIERS A TRIER
C NAVAN2 : LES 3 ENTIERS DEVANT SUBIR LES MEMES PERMUTATIONS
C
C SORTIE :
C --------
C NAPRE1 : LES 3 ENTIERS APRES LE TRI
C NAPRE2 : LES 3 ENTIERS AYANT SUBI LES MEMES PERMUTATIONS
C
C REMARQUE : EN ENTREE NAVAN1 ET NAPRE1 PEUVENT ETRE CONFONDUS
C            EN ENTREE NAVAN2 ET NAPRE2 PEUVENT ETRE CONFONDUS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC  SEPTEMBRE 1991
C....................................................................012
      INTEGER  NAVAN1(3),NAVAN2(3),NAPRE1(3),NAPRE2(3)
C
C     PROTECTION DES 3 VALEURS D'ENTREE
      DO 5 I=1,3
         NAPRE1(I) = NAVAN1(I)
         NAPRE2(I) = NAVAN2(I)
 5    CONTINUE
C
C     LE TRI
C
C     SI   NAPRE1(1) > NAPRE1(2) ILS SONT PERMUTES
      N = NAPRE1(1)
      IF( N .GT. NAPRE1(2) )  THEN
         NAPRE1(1) = NAPRE1(2)
         NAPRE1(2) = N
         N         = NAPRE2(1)
         NAPRE2(1) = NAPRE2(2)
         NAPRE2(2) = N
      ENDIF
C     ICI NAPRE1(1) < NAPRE1(2)
C
C     SI  NAPRE1(2) > NAPRE1(3) ILS SONT PERMUTES
      N = NAPRE1(2)
      IF( N .GT. NAPRE1(3) ) THEN
         NAPRE1(2) = NAPRE1(3)
         NAPRE1(3) = N
C        ICI NAPRE1(2) < NAPRE1(3)
         N         = NAPRE2(2)
         NAPRE2(2) = NAPRE2(3)
         NAPRE2(3) = N
C
C        SI  NAPRE1(1) > NAPRE1(2) ILS SONT PERMUTES
         N = NAPRE1(1)
         IF( N .GT. NAPRE1(2) )  THEN
            NAPRE1(1) = NAPRE1(2)
            NAPRE1(2) = N
            N         = NAPRE2(1)
            NAPRE2(1) = NAPRE2(2)
            NAPRE2(2) = N
         ENDIF
      ENDIF
      END
