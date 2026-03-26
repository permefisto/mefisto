      SUBROUTINE ISOMVA( DMATRI, XYZE, XYZS )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LES 3 COORDONNEES RESULTAT DE L'ISOMETRIE
C -----    DEFINIE PAR LA MATRICE DMATRI
C
C ENTREES:
C --------
C DMATRI : MATRICE(3,4) DEFINISSANT L'ISOMETRIE
C XYZE   : LES 3 COORDONNEES A TRANSFORMER
C
C SORTIE :
C --------
C XYZS   : LES 3 COORDONNEES RESULTANTES
C          XYZS(3) = DMATRI(3,1:3) * XYZE(3) + DMATRI(3,4)
C
C ATTENTION: SI XYZE = XYZS EN ENTREE RISQUE D'ERREUR !
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS  JANVIER 1986
C.....................................................................
      DOUBLE PRECISION  DMATRI(3,4), D
      REAL              XYZE(3), XYZS(3)
      INTRINSIC         REAL
C
      DO 20 I=1,3
         D = DMATRI(I,4)
         DO 10 J=1,3
            D = D + DMATRI(I,J) * XYZE(J)
 10      CONTINUE
         XYZS(I) = REAL( D )
 20   CONTINUE
C
      RETURN
      END
