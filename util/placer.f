      SUBROUTINE PLACER( N1 , N2 , A , G )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : PLACER A(N1) A SA PLACE DEFINITIVE M ENTRE N1 ET N2 DE A
C -----
C ENTREE :
C --------
C N1,N2  : BORNES DU TRI DANS A
C
C ENTREES ET SORTIES :
C --------------------
C A      : TABLEAU A TRIER SUR LUI-MEME
C
C SORTIE :
C --------
C G      : POSITION DEFINITIVE DE A(N1) INITIAL
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS       MARS 199O
C23456---------------------------------------------------------------012
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      INTEGER A(1:*)
      INTEGER D,G,AM,AUX
C     D,G POINTEUR SUR LE COEFFICIENT DROITE ET GAUCHE DANS A
      AM = A(N1)
      G = N1+1
      D = N2
C
 10   IF( G .LT. D ) THEN
         IF( AM .GT. A(G) ) THEN
            G = G + 1
         ELSE
            IF( AM .LE. A(D) ) THEN
               D = D - 1
            ELSE
               AUX  = A(G)
               A(G) = A(D)
               A(D) = AUX
               D    = D - 1
               G    = G + 1
            ENDIF
         ENDIF
         GOTO 10
      ENDIF
C
C     SI IN-DECREMENTATION INCORRECTE: RECTIFICATION
      IF( AM .LE. A(G) ) G = G - 1
C
C     AM EST POSITIONNE A SA PLACE DEFINITIVE PAR PERMUTATION AVEC A(G)
      A(N1) = A(G)
      A(G) = AM
      END
