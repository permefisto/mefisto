      SUBROUTINE TRIRAR( N , A , L , LAPILE , IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : TRIER PAR QUICK SORT LE TABLEAU A DE NOMBRE REELS
C
C ENTREE :
C --------
C N      : NOMBRE DE TERMES A TRIER DANS A
C
C ENTREES ET SORTIES :
C --------------------
C A      : TABLEAU A TRIER SUR LUI-MEME
C
C AUXILIAIRE :
C ------------
C L      : NOMBRE D'ENTIERS DE LAPILE
C LAPILE : PILE DES ARGUMENTS D'APPEL DU TRI RAPIDE
C
C SORTIE :
C --------
C IERR   : 0 SI PAS DE PROBLEME
C          >0 SI TRI NON REALISE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS       MARS 199O
C23456---------------------------------------------------------------012
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      INTEGER           LAPILE(1:2,1:L),G,D
      REAL              A(1:N),AM,AUX
C
C     1,N EST EMPILE
      LH = 1
      LAPILE(1,1) = 1
      LAPILE(2,1) = N
C
C     TANT QUE LA PILE EST NON VIDE FAIRE
 10   IF( LH .GT. 0 ) THEN
C
C        LECTURE DU HAUT DE LA PILE
         N1 = LAPILE(1,LH)
         N2 = LAPILE(2,LH)
C
C        DEPILER
         LH = LH - 1
C
C        TRAITEMENT N1 N2
         IF( N2 .GT. N1 ) THEN
C           G PLACE DEFINITIVE DANS A DE LA VALEUR A(N1)
C           A DROITE DE G A(I) >= A(G)
C             GAUCHE           <  A(G)
CCC            CALL PLACER( N1 , N2 , A , G )
C           D,G POINTEUR SUR LE COEFFICIENT DROITE ET GAUCHE DANS A
            AM = A(N1)
            G  = N1+1
            D  = N2
C
 20         IF( G .LT. D ) THEN
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
               GOTO 20
            ENDIF
C
C           SI IN-DECREMENTATION INCORRECTE: RECTIFICATION
            IF( AM .LE. A(G) ) G = G - 1
C
C           AM EST POSITIONNE A SA PLACE DEFINITIVE PAR PERMUTATION AVEC A(G)
            A(N1) = A(G)
            A(G) = AM
C
C           Y A T IL ASSEZ DE PLACE DANS LA PILE?
            IF( LH + 2 .GT. L ) THEN
               WRITE(IMPRIM,*) 'PILE SATUREE DANS LE TRI'
               IERR = 1
               RETURN
            ENDIF
C
C           EMPILER G+1,N2
            LH = LH + 1
            LAPILE(1,LH) = G+1
            LAPILE(2,LH) = N2
C
C           EMPILER N1,G-1
            LH = LH + 1
            LAPILE(1,LH) = N1
            LAPILE(2,LH) = G-1
         ENDIF
         GOTO 10
      ENDIF
      END
