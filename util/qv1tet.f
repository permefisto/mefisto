      SUBROUTINE QV1TET( PTXYZD, NOSOTE, VTE, QTE )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   CALCUL DE LA QUALITE ET DU VOLUME D'UN TETRAEDRE
C -----   DEFINI PAR SES 4 NUMEROS DE SOMMETS DANS LE TABLEAU PTXYZD

C ENTREES:
C --------
C PTXYZD : X Y Z DISTANCE SOUHAITEE DES SOMMETS DU MAILLAGE
C NOSOTE : NUMERO PTXYZD DES 4 SOMMETS DU TETRAEDRE

C SORTIES:
C --------
C VTE : LE VOLUME DU TETRAEDRE (>0 SI ORIENTE COMME UN REPERE)
C                               <0 SI TETRAEDRE DEGENERE OU MAL ORIENTE)
C QTE : QUALITE DU TETRAEDRE 
C       VTE<0    QTE=-1  SOMMET DE L'AUTRE COTE DE LA FACE OPPOSEE
C       VTE=0    QTE= 0  4 SOMMETS DANS UN MEME PLAN
C       VTE>0  0<QTE<=1  TETRAEDRE CORRECT
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  Saint Pierre du Perray              Mai 2020
C2345X7..............................................................012
      DOUBLE PRECISION  PTXYZD(4,*), VTE, ARMIN, ARMAX, SURFTR(4)
      REAL              QTE
      INTEGER           NOSOTE(4)

      CALL QUATETD( PTXYZD( 1, NOSOTE(1) ),
     %              PTXYZD( 1, NOSOTE(2) ),
     %              PTXYZD( 1, NOSOTE(3) ),
     %              PTXYZD( 1, NOSOTE(4) ),
     %              ARMIN, ARMAX, SURFTR, VTE, QTE )

      RETURN
      END
