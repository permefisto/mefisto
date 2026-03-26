      SUBROUTINE ANSETO( NS,    N1FEOC, NFETOI, XYZSOM, COSE2P,
     %                   NBNORM )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   CALCULER LE NOMBRE DE NORMALES AUX FACES AUTOUR DU SOMMET NS
C -----
C
C ENTREES:
C --------
C NS     : LE POINT EN COURS DE TRAITEMENT
C N1FEOC : NUMERO DE LA PREMIERE FACE DE NFETOI
C NFETOI : LES 3 SOMMETS DE CHAQUE FACE ET FACE SUIVANTE DE L'ETOILE
C XYZSOM : LES 3 COORDONNEES DES SOMMETS
C COSE2P : COSINUS DU PETIT ANGLE DE COPLANEARITE
C
C SORTIE :
C --------
C NBNORM : NOMBRE DE NORMALES AUX FACES DE SOMMET NS ( <=3 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC   DECEMBRE 1991
C....................................................................012
      COMMON / UNITES / LECTEU,IMPRIM,NUNITE(30)
      INTEGER           NFETOI(5,*)
      REAL              XYZSOM(3,*),NORMAL(3,3)
C
      NBNORM = 0
C
C     POSITION DANS NFETOI DE LA 1-ERE FACE DE L'ETOILE
      NF = N1FEOC
C
C     TANT QU'IL EXISTE UNE FACE DE L'ETOILE FAIRE
 10   IF( NF .GT. 0 ) THEN
C
C        CETTE FACE CONTIENT ELLE LE SOMMET NS ?
         DO 20 I=1,3
            IF( NFETOI(I,NF) .EQ. NS ) GOTO 30
 20      CONTINUE
         GOTO 90
C
C        FACE DE SOMMET NS  CALCUL DU VECTEUR NORMAL
 30      CALL NORFA3( XYZSOM(1,NFETOI(1,NF)),
     %                XYZSOM(1,NFETOI(2,NF)),
     %                XYZSOM(1,NFETOI(3,NF)), NORMAL(1,NBNORM+1),
     %                I )
         IF( I .NE. 0 ) GOTO 90
C
         IF( NBNORM .LE. 0 ) THEN
C           IL N'EXISTE PAS D'AUTRE NORMALE
            NBNORM = 1
         ELSE
C           CALCUL DU PRODUIT SCALAIRE
            N = 0
            DO 40 I=1,NBNORM
               COSA = PROSCR( NORMAL(1,I), NORMAL(1,NBNORM+1), 3 )
               COSA = ABS( COSA )
               IF( 1.0 - COSA .LE. 1.0 - COSE2P ) THEN
C                 NORMALE IDENTIFIEE
                  N = I
               ENDIF
 40         CONTINUE
C
            IF( N .EQ. 0 ) THEN
C              UNE NOUVELLE NORMALE
               NBNORM = NBNORM + 1
               IF( NBNORM .EQ. 3 ) RETURN
            ENDIF
         ENDIF
C
C        LA FACE SUIVANTE
 90      NF = NFETOI(5,NF)
         GOTO 10
      ENDIF
      END
