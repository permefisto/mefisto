      SUBROUTINE INTCOT( X, NBINTV, XI, INDICE )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RECHERCHER L'INDICE DE L'INTERVALLE CONTENANT X PARMI LES XI
C -----    LE PREMIER ETANT EN FAIT LE ZERO-EME
C ENTREES :
C ---------
C X      : L'ABSCISSE DU POINT
C NBINTV : LE NOMBRE D'INTERVALLES
C XI     : LES ABSCISSES DES EXTREMITES DES INTERVALLES
C
C SORTIES :
C ---------
C INDICE : INDICE DE 0 A NBINTV-1 DE L'INTERVALLE CONTENANT X
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS       SEPTEMBRE 1993
C234567..............................................................012
      REAL  XI(0:NBINTV)
C
      IF( X .LE. XI(1) ) THEN
         INDICE = 0
         RETURN
      ENDIF
      DO 10 INDICE=1,NBINTV-2
         IF( X .LE. XI(INDICE+1) ) RETURN
 10   CONTINUE
      INDICE = NBINTV - 1
      END
