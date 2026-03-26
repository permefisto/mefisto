       SUBROUTINE NOFAPH( NCOGEL, NOSOEF, FACE )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    FACE CONTIENT EN SORTIE LE NUMERO DES SOMMETS DE LA OU
C -----    LES 3 FACES QUADRANGULAIRES NE CONTENANT PAS LE PLUS PETIT
C          NUMERO DE SOMMET DE L'EF
C
C ENTREE :
C--------
C NCOGEL : CODE GEOMETRIQUE DE L'EF
C          5:TETRAEDRE, 6:PENTAEDRE, 7:HEXAEDRE, 9:PYRAMIDE A BASE CARREE
C NOSOEF : NUMEROS DES 5 OU 6 OU 8 SOMMETS DE L'EF
C
C SORTIE :
C --------
C FACE   : NUMERO DES SOMMETS > PLUS PETIT NUMERO DE SOMMET DE L'EF
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEURS: F.CHOUKROUN O.RICOU UPMC ANALYSE NUMERIQUE PARIS JANVIER 1990
C MODIFS : PERRONNET ALAIN Laboratoire J-L. LIONS UPMC PARIS   Mars 2007
C--------------------------------------------------------------------012
       INTEGER NOSOEF(8), FACE(4,3)
C
       IF( NCOGEL .EQ. 9 ) THEN
C
C         PYRAMIDE
          FACE(1,1) = NOSOEF(1)
          FACE(2,1) = NOSOEF(2)
          FACE(3,1) = NOSOEF(3)
          FACE(4,1) = NOSOEF(4)
C
       ELSE IF( NCOGEL .EQ. 6 ) THEN
C
C         PENTAEDRE
          FACE(1,1) = NOSOEF(2)
          FACE(2,1) = NOSOEF(3)
          FACE(3,1) = NOSOEF(6)
          FACE(4,1) = NOSOEF(5)
C
       ELSE
C
C         HEXAEDRE
          FACE(1,1) = NOSOEF(5)
          FACE(2,1) = NOSOEF(6)
          FACE(3,1) = NOSOEF(7)
          FACE(4,1) = NOSOEF(8)
C
          FACE(1,2) = NOSOEF(2)
          FACE(2,2) = NOSOEF(3)
          FACE(3,2) = NOSOEF(7)
          FACE(4,2) = NOSOEF(6)
C
          FACE(1,3) = NOSOEF(3)
          FACE(2,3) = NOSOEF(4)
          FACE(3,3) = NOSOEF(8)
          FACE(4,3) = NOSOEF(7)
C
       ENDIF
       END
