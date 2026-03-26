        SUBROUTINE VOEXP4(NBSA,NBSH,NPSOM,NSENS,N1,N2,N3)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERER LES CONSTANTES DE NUMEROTATION N1 N2 ET N3
C -----    DES SOMMETS D'UNE FACE D'UN PENTAEDRE  SUIVANT
C          LE NUMERO DU PREMIER ET LE SENS DE PARCOURS
C
C ENTREES:
C --------
C NBSA,NBSH  : NOMBRE DE SOMMETS DANS CHAQUE DIRECTION DE LA FACE
C NBX ,NBY   : MEME CHOSE POUR LE QUADRILATERE
C NPSOM      : NUMERO DU SOMMET COINCIDANT AVEC LE SOMMET 1 DE LA FACE
C NSENS      : SENS DE PARCOURS DE LA FACE
C
C SORTIES:
C --------
C N1 N2 N3 : LES CONSTANTES DE NUMEROTATION
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY      ANALYSE NUMERIQUE PARIS   DECEMBRE  1988
C.......................................................................
      IF ( NSENS .EQ. 1 ) THEN
         IF      ( NPSOM .EQ. 1 ) THEN
            N1 = - NBSH
            N2 = 1
            N3 = NBSH
         ELSE IF ( NPSOM .EQ. 2 ) THEN
            N1 = 1
            N2 = NBSA
            N3 = -1
         ELSE IF ( NPSOM .EQ. 3 ) THEN
            N1 = NBSH * ( NBSA + 1 ) + 1
            N2 = - 1
            N3 = - NBSH
         ELSE IF ( NPSOM .EQ. 4 ) THEN
            N1 = NBSA * NBSH
            N2 = - NBSA
            N3 = 1
         END IF
      ELSE
         IF      ( NPSOM .EQ. 1 ) THEN
            N1 = - NBSA
            N2 = NBSA
            N3 = 1
         ELSE IF ( NPSOM .EQ. 2 ) THEN
            N1 = 1
            N2 = - 1
            N3 = NBSH
         ELSE IF ( NPSOM .EQ. 3 ) THEN
            N1 = NBSA * ( NBSH + 1 ) + 1
            N2 = - NBSA
            N3 = - 1
         ELSE IF ( NPSOM .EQ. 4 ) THEN
            N1 = NBSH * NBSA
            N2 = 1
            N3 = - NBSH
         END IF
      END IF
      END
