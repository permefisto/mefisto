         SUBROUTINE  VOEXT3(NBSA,POXY1,POXY2,POXY3,COSO)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    MAILLER UN TETRAEDRE DEFINI PAR SES FACES
C -----    PAR LA METHODE DE PERRONNET, PIERROT ET VAZEILLES
C ENTREES :
C ---------
C NBSA   : NOMBRE DE SOMMETS PAR ARETE
C
C ENTREES ET SORTIES :
C --------------------
C COSO   : COORDONNEES DES SOMMETS DU MAILLAGE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY ANALYSE NUMERIQUE UPMC PARIS MARS 1989
C23456---------------------------------------------------------------012
C
      REAL   COSO(3,*),POXY1(NBSA),POXY2(NBSA),POXY3(NBSA)
C
      NSPLAN = 10
      DO 1 N1 = 4 , NBSA - 1
C        ON APPLIQUE LA METHODE DU TRIANGLE ALGEBRIQUE
         CALL SUEXT2(N1,POXY1,POXY2,POXY3,COSO(1,NSPLAN+1))
         NSPLAN = NSPLAN + (N1 * N1 + N1 ) / 2
 1    CONTINUE
C
      RETURN
      END
