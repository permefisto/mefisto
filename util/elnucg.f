      SUBROUTINE ELNUCG( NUTYEL , NCOGEL )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   FOURNIR SELON LE NUMERO DU TYPE DE L'ELEMENT FINI
C -----   SON CODE GEOMETRIQUE
C
C PARAMETRES D ENTREE :
C ---------------------
C NUTYEL : NUMERO DU TYPE DE L'ELEMENT FINI
C
C PARAMETRES RESULTAT :
C ---------------------
C NCOGEL : CODE GEOMETRIQUE DE L ELEMENT FINI
C          1:NOEUD     2:SEGMENT  3:TRIANGLE  4:QUADRANGLE  5:TETRAEDRE
C          6:PENTAEDRE 7:HEXAEDRE 8:6-CUBE    9:PYRAMIDE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE UPMC PARIS    OCTOBRE 1990
C2345X7..............................................................012
      include"./incl/gsmenu.inc"
      CHARACTER*4    LISTEL(3,9), KNOMEL(2)
      DATA LISTEL / 'NOEU' , 'D   ' , '    ' ,
     &              'SEGM' , 'ENT ' , '    ' ,
     &              'TRIA' , 'NGLE' , '    ' ,
     &              'QUAD' , 'RANG' , 'LE  ' ,
     &              'TETR' , 'AEDR' , 'E   ' ,
     &              'PENT' , 'AEDR' , 'E   ' ,
     &              'HEXA' , 'EDRE' , '    ' ,
     &              '6CUB' , 'E   ' , '    ' ,
     &              'PYRA' , 'MIDE' , '    ' /
C
C     LE NOM DE L'ELEMENT FINI
      CALL ELNUNM( NUTYEL, KNOMEL )
C
C     RECHERCHE DU CODE GEOMETRIQUE
      DO 10 NCOGEL=1,9
         IF( KNOMEL(1) .EQ. LISTEL(1,NCOGEL) ) RETURN
   10 CONTINUE
C
C     ERREUR
C     ------
      NBLGRC(NRERR) = 1
      KERR(1) = 'ELNUCG: '//KNOMEL(1)// ' ' // KNOMEL(2) //' INCONNU'
      CALL LEREUR
      NCOGEL = 0
      RETURN
      END
