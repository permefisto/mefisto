      SUBROUTINE PILEOB( NOTYPE, NOOBJT, MNDFSO,
     %                   LHPILE, MXPILE, MNPILE )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : EMPILER LE TYPE, LE NUMERO D'UN OBJET, ADRESSE MCN 'DEFINITION'
C -----
C
C ENTREES:
C --------
C NOTYPE : NUMERO DU TYPE 1:POINT 2:LIGNE ... 6:ISOMETRIE
C NOOBJT : NUMERO DE L'OBJET A EMPILER
C MNDFSO : L'ADRESSE DANS MCN DU TABLEAU 'DEFINITION' DU SUPER-OBJET
C          DONT L'OBJET COURANT EST UN CONSTITUANT
C
C MODIFIES:
C ---------
C LHPILE : HAUTEUR ACTUELLE DE LA PILE
C MXPILE : NOMBRE  MOTS DU TABLEAU PILE
C MNPILE : ADRESSE MCN DU TABLEAU PILE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS            MAI 1987
C.....................................................................
      include"./incl/pp.inc"
      COMMON    MCN(MOTMCN)
      INTEGER   NOOBJT
C
C     ADRESSE 1-ER MOT A EMPILER
      MN = MNPILE + LHPILE
C
C     LA HAUTEUR DE LA PILE EST MISE A JOUR
      LHPILE = LHPILE + 3
C     LE TABLEAU EST-IL ASSEZ GRAND ?
      IF( LHPILE .GT. MXPILE ) THEN
C        LA PILE TROP COURTE EST ALLONGEE
         CALL TNMCAU( 'ENTIER', MXPILE, LHPILE, MXPILE, MNPILE )
         MXPILE = LHPILE
      ENDIF
C
C     LES NBOBJT OBJETS SONT EMPILES PAR LA FIN DU TABLEAU NOOBJT
C     LE TYPE DE L'OBJET
      MCN( MN ) = NOTYPE
ccc      print *,'pileob: mn=',mn,' notype=',notype
C     LE NUMERO DE L'OBJET
      MCN( MN + 1 ) = ABS( NOOBJT )
ccc      print *,'pileob: mn+1=',mn+1,' no objet=',noobjt, MCN( MN + 1 )
C     L'ADRESSE MCN DU TABLEAU 'DEFINITION' DU SUPER-OBJET
      MCN( MN + 2 ) = MNDFSO
ccc      print *,'pileob: mn+2=',mn+2,' mndfso=',MCN( MN + 2 )
C
      RETURN
      END
