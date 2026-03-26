      SUBROUTINE IMPTAD( NBLIGN, NBCOLO, VP )
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    AFFICHER UN TABLEAU REEL DOUBLE PRECISION (NBLIGN,NBCOLO)
C -----    LIGNE PAR LIGNE
C
C ENTREES:
C --------
C NBLIGN : LE NOMBRE DE COLONNES DU TABLEAU VP
C NBCOLO : LE NOMBRE DE LIGNES   DU TABLEAU VP
C VP     : LE TABLEAU A AFFICHER
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET      ANALYSE NUMERIQUE UPMC PARIS    AOUT 1998
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      DOUBLE PRECISION VP(NBLIGN,NBCOLO)
      COMMON /UNITES/  LECTEU,IMPRIM,NUNITE(30)
C
10001 FORMAT(130('*')/T30,'Les ',I4,' VECTEURS de ',I6,' COMPOSANTES
     +'/130('*'))
20001 FORMAT(130('*')/T30,'The ',I4,' VECTORS of',I6,' COMPOSANTES
     +'/130('*'))
10003 FORMAT('VECTEUR ',I4)
20003 FORMAT('VECTOR ',I4)
10004 FORMAT(10G15.7)
C
      IF( LANGAG .EQ. 0 ) THEN
         WRITE (IMPRIM,10001) NBLIGN, NBCOLO
         DO 2 I=1,NBLIGN
            WRITE (IMPRIM,10003) I
            WRITE (IMPRIM,10004) (VP(I,J),J=1,NBCOLO)
 2       CONTINUE
      ELSE
         WRITE (IMPRIM,20001) NBLIGN, NBCOLO
         DO 4 I=1,NBLIGN
            WRITE (IMPRIM,20003) I
            WRITE (IMPRIM,10004) (VP(I,J),J=1,NBCOLO)
 4       CONTINUE
      ENDIF
      END
