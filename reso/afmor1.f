      SUBROUTINE AFMOR1( NTDL , LPLIGN , LPCOLO )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT: AFFICHER LA STRUCTURE DE LA MATRICE MORSE
C ----
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  PARAMETRES D'ENTREE :
C  -------------------
C NTDL   : ORDRE DE LA MATRICE
C LPLIGN :
C LPCOLO : LES POINTEURS DE LA MATRICE MORSE
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY  ANALYSE NUMERIQUE UPMC PARIS      MAI 1990
C23456---------------------------------------------------------------012
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      INTEGER           LPLIGN(NTDL+1),LPCOLO(1)
C
 100  FORMAT(//1X,80(1H=)/,' MATRICE MORSE DE RANG',I5)
 200  FORMAT(1X,80(1H-)/,' LIGNE',I8,'  ADRESSES',2I5/
     &       100(' COLONNES',20I5/))
C
      WRITE (IMPRIM,100 ) NTDL
       DO 1 I=1,NTDL
          I1 = LPLIGN(I)+1
          I2 = LPLIGN(I+1)
          WRITE (IMPRIM,200 ) I,I1,I2,
     &          (LPCOLO(K),K=I1,I2)
 1     CONTINUE
C
       RETURN
       END
