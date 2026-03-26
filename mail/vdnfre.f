      SUBROUTINE VDNFRE( NF0, NBFAPE, NOFAPE, LEFACO )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RETIRER LA FACE NF0 DE LEFACO DE LA LISTE NOFAPE
C -----    AINSI QUE TOUTES LES FACES VOISINES PARTAGEANT L'UNE
C          DES 3 ARETES DE LA FACE NF0

C ENTREES:
C --------
C NF0    : NUMERO DE LA FACE A RETIRER  DES FACES PERDUES
C NBFAPE : NOMBRE DE FACES PERDUES

C MODIFIE:
C --------
C NOFAPE : NUMERO DES FACES PERDUES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC        MAI 1993
C2345X7..............................................................012
      INTEGER           NOFAPE(1:NBFAPE), LEFACO(11,0:*)

      DO 10 I=1,NBFAPE
         NFP = NOFAPE(I)
         IF( NFP .LE. 0   ) GOTO 10
         IF( (NFP .EQ. NF0) .OR.
     %       (NFP .EQ. LEFACO(6,NF0)) .OR.
     %       (NFP .EQ. LEFACO(7,NF0)) .OR.
     %       (NFP .EQ. LEFACO(8,NF0)) ) THEN
C            FACE CENSEE ETRE RETROUVEE
             NOFAPE(I) = 0
         ENDIF
 10   ENDDO

      RETURN
      END
