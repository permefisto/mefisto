      SUBROUTINE TRI3SAR( NOSOMI, NOAREI, NOSOMF, NOAREF )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRI CROISSANT DU NUMERO DES 3 SOMMETS EN COORDONNANT
C -----    POUR CONSERVER LA NUMEROTATION DES ARETES

C ENTREE :
C --------
C NOSOMI : LES 3 SOMMETS A TRIER
C NOAREI : LES 3 ARETES A COORDONNER

C SORTIE :
C --------
C NOSOMF : LES 3 SOMMETS APRES TRI
C NOAREF : LES 3 ARETES COORDONNEES

C REMARQUE : EN ENTREE NOSOMI ET NOSOMF PEUVENT ETRE CONFONDUS
C            EN ENTREE NOAREI ET NOAREF PEUVENT ETRE CONFONDUS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC  SEPTEMBRE 1991
C....................................................................012
      INTEGER  NOSOMI(3),NOAREI(3),NOSOMF(3),NOAREF(3)

C     PROTECTION DES 3 VALEURS D'ENTREE
      DO I=1,3
         NOSOMF(I) = NOSOMI(I)
         NOAREF(I) = NOAREI(I)
      ENDDO

C     LE TRI

C     SI   NOSOMF(1) > NOSOMF(2) ILS SONT PERMUTES
      N = NOSOMF(1)
      IF( N .GT. NOSOMF(2) )  THEN
         NOSOMF(1) = NOSOMF(2)
         NOSOMF(2) = N
C        PERMUTATION DES ARETES 2 ET 3
         N         = NOAREF(2)
         NOAREF(2) = NOAREF(3)
         NOAREF(3) = N
      ENDIF

C     ICI NOSOMF(1) < NOSOMF(2)

C     SI  NOSOMF(2) > NOSOMF(3) ILS SONT PERMUTES
      N = NOSOMF(2)
      IF( N .GT. NOSOMF(3) ) THEN
         NOSOMF(2) = NOSOMF(3)
         NOSOMF(3) = N
C        ICI NOSOMF(2) < NOSOMF(3)
C        PERMUTATION DES ARETES 1 ET 3
         N         = NOAREF(1)
         NOAREF(1) = NOAREF(3)
         NOAREF(3) = N

C        SI  NOSOMF(1) > NOSOMF(2) ILS SONT PERMUTES
         N = NOSOMF(1)
         IF( N .GT. NOSOMF(2) )  THEN
            NOSOMF(1) = NOSOMF(2)
            NOSOMF(2) = N
C           PERMUTATION DES ARETES 2 ET 3
            N         = NOAREF(2)
            NOAREF(2) = NOAREF(3)
            NOAREF(3) = N
         ENDIF
      ENDIF

      RETURN
      END
