      SUBROUTINE AJ1ESR( N, NBENT, ENTIER )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    AJOUTER L'ENTIER N DANS LE TABLEAU ENTIER(1:NBENT)
C -----    SANS REPETITION C-A-D
C          N EST AJOUTE EN NBENT SI NON RETROUVE
C
C ENTREES:
C --------
C N      : L'ENTIER A AJOUTER
C
C MODIFIE :
C ---------
C NBENT : NOMBRE D'ENTIERS STOCKES AVANT ET APRES DANS LE TABLEAU ENTIER
C
C SORTIES:
C --------
C ENTIER : LE TABLEAU DES ENTIERS STOCKES SANS REPETITION
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC   NOVEMBRE 1992
C2345X7..............................................................012
      INTEGER  ENTIER(1:*)
C
      DO 10 K=1,NBENT
         IF( N .EQ. ENTIER(K) ) RETURN
 10   CONTINUE
C
C     ENTIER NON RETROUVE => AJOUTE
      NBENT = NBENT + 1
      ENTIER( NBENT ) = N
      END
