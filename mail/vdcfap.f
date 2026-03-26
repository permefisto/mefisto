      SUBROUTINE VDCFAP( NBFAPE, NOFAPE )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    DETRUIRE LES FACES RETROUVEES DU TABLEAU NOFAPE
C -----
C
C ENTREES ET SORTIES :
C --------------------
C NBFAPE : NOMBRE DE FACES PERDUES  AVANT ET APRES
C NOFAPE : NUMERO DES NBFAPE FACES PERDUES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS UPMC  SEPTEMBRE 1991
C....................................................................012
      INTEGER           NOFAPE(NBFAPE)
C
C     LES FACES NULLES DU TABLEAU NOFAPE SONT DETRUITES ET COMPRESSEES
      NBFARE = 0
      DO 10 I=1,NBFAPE
         IF( NOFAPE(I) .EQ. 0 ) THEN
            NBFARE = NBFARE + 1
         ELSE
            NOFAPE(I-NBFARE) = ABS( NOFAPE(I) )
         ENDIF
 10   CONTINUE
      NBFAPE = NBFAPE - NBFARE

      RETURN
      END
