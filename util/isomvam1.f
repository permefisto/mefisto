      SUBROUTINE ISOMVAM1( DFM1, XYZE, XYZS )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LES 3 COORDONNEES RESULTAT DE L'ISOMETRIE INVERSE
C -----    DEFINIE PAR LA MATRICE DFM1 et VECT
C
C ENTREES:
C --------
C DFM1   : MATRICE(3,4) DEFINISSANT L'ISOMETRIE INVERSE
C XYZE   : LES 3 COORDONNEES A TRANSFORMER
C
C SORTIE :
C --------
C XYZS   : LES 3 COORDONNEES RESULTANTES
C          XYZS(3) = DFM1(3,1:3) * (XYZE(3) - DFM1(3,4))
C
C ATTENTION: SI XYZE = XYZS EN ENTREE RISQUE D'ERREUR !
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY Novembre 2011
C2345X7..............................................................012
      DOUBLE PRECISION  DFM1(3,4), D
      REAL              XYZE(3), XYZS(3)
      INTRINSIC         REAL
C
      DO I=1,3
         D = 0D0
         DO J=1,3
            D = D + DFM1(I,J) * ( XYZE(J) - DFM1(I,4) )
         ENDDO
         XYZS(I) = REAL( D )
      ENDDO
C
      RETURN
      END
