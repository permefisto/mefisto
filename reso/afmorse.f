      SUBROUTINE AFMORSE( NBLIGN, NTDL, LPLIGN, LPCOLO, A )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT:     AFFICHER LES COEFFICIENTS DE LA MATRICE MORSE
C ----
C ENTREES:
C --------
C NBLIGN : NOMBRE DE LIGNES A AFFICHER
C NTDL   : NOMBRE DE LIGNES ET COLONNES DE LA MATRICE MORSE
C LPLIGN : POINTEUR SUR LE DERNIER COEFFICIENT DE CHAQUE LIGNE
C LPCOLO : NUMERO DE COLONNE DES COEFFICIENTS STOCKES DE LA MATRICE MORSE
C A      : COEFFICIENTS DE LA MATRICE MORSE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET TIMS NTU TAIPEI TAIWAN          NOVEMBRE 2009
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
      INTEGER           LPLIGN(0:NTDL), LPCOLO(*)
      DOUBLE PRECISION  A(*)
C
10000 FORMAT(/,' Une MATRICE MORSE de',I11,' LIGNES'/1X,100(1H=))
10010 FORMAT( 5( '  A(',I5,',',I5,')=',D16.8 ) )
C
20000 FORMAT(/,' A CONDENSED MATRIX of ',I11,' LINES'/1X,100(1H=))
20020 FORMAT( 5( '  A(',I5,',',I5,')=',D16.8 ) )
C
      IF( LANGAG .EQ. 0 ) THEN
         WRITE (IMPRIM,10000 ) NTDL
         DO 10 I = 1, NBLIGN
            I1 = LPLIGN(I-1) + 1
            I2 = LPLIGN(I)
            WRITE(IMPRIM,10010) (I,LPCOLO(K),A(K), K=I1,I2)
 10      CONTINUE
      ELSE
         WRITE (IMPRIM,20000 ) NTDL
         DO 20 I = 1, NBLIGN
            I1 = LPLIGN(I-1) + 1
            I2 = LPLIGN(I)
            WRITE(IMPRIM,20020) (I,LPCOLO(K),A(K), K=I1,I2)
 20      CONTINUE
      ENDIF
C
      RETURN
      END
