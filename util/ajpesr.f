      SUBROUTINE AJPESR( NB, ENT, NBENT, ENTIER )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    AJOUTER LES NB ENTIERS ENT  DANS LE TABLEAU ENTIER(1:NBENT)
C -----    SANS REPETITION C-A-D
C          N EST AJOUTE EN NBENT SI NON RETROUVE
C
C ENTREES:
C --------
C NB     : LE NOMBRE D'ENTIERS A AJOUTER
C ENT    : LES NB ENTIERS A AJOUTER
C
C MODIFIE :
C ---------
C NBENT : NOMBRE D'ENTIERS STOCKES AVANT ET APRES DANS LE TABLEAU ENTIER
C
C SORTIES:
C --------
C ENTIER : LE TABLEAU DES ENTIERS STOCKES SANS REPETITION
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC   DECEMBRE 1992
C2345X7..............................................................012
      INTEGER  ENT(1:NB),ENTIER(1:*)
C
      DO 20 I=1,NB
         N = ENT( I )
         DO 10 K=1,NBENT
            IF( N .EQ. ENTIER(K) ) GOTO 20
 10      CONTINUE
C
C        ENTIER NON RETROUVE => AJOUTE
         NBENT = NBENT + 1
         ENTIER( NBENT ) = N
 20   CONTINUE
      END
