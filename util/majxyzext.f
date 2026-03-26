      SUBROUTINE MAJXYZEXT( NBCOOR, XYZEFS )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   MISE A JOUR DU CADRE MINMAX DES NBCOOR COORDONNEES DES SOMMETS
C ----    D'UN MAILLAGE D'EF D'EXTREMES XYZEFS
C
C ENTREES :
C ---------
C NBCOOR : LE NOMBRE DE COORDONNEES D'UN SOMMET DU MAILLAGE
C XYZEFS : LES COORDONNEES EXTREMES DES SOMMETS OU NOEUDS DU MAILLAGE
C
C SORTIES :
C ---------
C DANS LE COMMON /XYZEXT/
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE DU PERRAY  JUILLET 2013
C23456---------------------------------------------------------------012
      include"./incl/xyzext.inc"
      REAL      XYZEFS(6,2)

      DO K=1,NBCOOR
         COOEXT(K,1) = MIN( COOEXT(K,1), XYZEFS(K,1) )
         COOEXT(K,2) = MAX( COOEXT(K,2), XYZEFS(K,2) )
      ENDDO

      RETURN
      END

