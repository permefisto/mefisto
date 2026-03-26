        SUBROUTINE VOEXH4(NBS1,NBS2,NPSOM,NSENS,N1,N2,N3)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERER LES CONSTANTES DE NUMEROTATION N1 N2 ET N3
C -----    DES SOMMETS D'UNE FACE SUIVANT LE NUMERO DU PREMIER
C          ET LE SENS DE PARCOURS
C
C ENTREES:
C --------
C NBS1,NBS2  : NOMBRE DE SOMMETS DANS CHAQUE DIRECTION DE LA FACE
C NPSOM    : NUMERO DU SOMMET COINCIDANT AVEC LE SOMMET 1 DE LA FACE
C NSENS    : SENS DE PARCOURS DE LA FACE
C
C SORTIES:
C --------
C N1 N2 N3 : LES CONSTANTES DE NUMEROTATION
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY      ANALYSE NUMERIQUE PARIS   DECEMBRE  1988
C.......................................................................
C
      IF ( NSENS .EQ. 1 ) THEN
         IF      ( NPSOM .EQ. 1 ) THEN
            N1 = - NBS1
            N2 = 1
            N3 = NBS1
         ELSE IF ( NPSOM .EQ. 2 ) THEN
            N1 = 1
            N2 = NBS2
            N3 = -1
         ELSE IF ( NPSOM .EQ. 3 ) THEN
            N1 = NBS1 * ( NBS2 + 1 ) + 1
            N2 = - 1
            N3 = - NBS1
         ELSE IF ( NPSOM .EQ. 4 ) THEN
            N1 = NBS1 * NBS2
            N2 = - NBS2
            N3 = 1
         END IF
      ELSE
         IF      ( NPSOM .EQ. 1 ) THEN
            N1 = - NBS2
            N2 = NBS2
            N3 = 1
         ELSE IF ( NPSOM .EQ. 2 ) THEN
            N1 = 1
            N2 = - 1
            N3 = NBS1
         ELSE IF ( NPSOM .EQ. 3 ) THEN
            N1 = NBS2 * ( NBS1 + 1 ) + 1
            N2 = - NBS2
            N3 = - 1
         ELSE IF ( NPSOM .EQ. 4 ) THEN
            N1 = NBS1 * NBS2
            N2 = 1
            N3 = - NBS1
         END IF
      END IF
C
      RETURN
      END
