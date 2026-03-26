      SUBROUTINE ELNMNU( NOMELT, NO )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    FOURNIR SELON LE NOM DU TYPE DE L ELEMENT (NOMELT(1),(2))
C -----    SON NO DANS LES DIVERS SOUS-PROGRAMMES UTILITAIRES
C
C PARAMETRES D ENTREE :
C ---------------------
C NOMELT : NOM SUR 2 MOTS DE 4 CARACTERES DU NOM DE L ELEMENT FINI
C          L INTERPOLATION A ETE PRISE EN COMPTE
C
C PARAMETRE RESULTAT :
C --------------------
C NO     : NO DU TYPE D'ELEMENT FINI DANS LES SP UTILITAIRES
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS      JANVIER 1991
C ......................................................................
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      CHARACTER*4       NOMELT(2)
      include"./incl/gsmenu.inc"
      include"./incl/nomele.inc"
C
C     RECHERCHE DE L ELEMENT
C     ----------------------
      DO 200 J=1,NBMXEL
         DO 150 I=1,2
            IF( NOMELT(I) .NE. NOMELE(I,J) ) GOTO 200
  150    CONTINUE
         NO = J
         RETURN
  200 CONTINUE
C
C     ERREUR ELEMENT FINI INCONNU
C     ---------------------------
      NBLGRC(NRERR) = 1
      KERR(1) ='ELNMNU:ELEMENT FINI INCONNU'
      KERR(2) = NOMELT(1)//' '//NOMELT(2)
      CALL LEREUR
      WRITE(IMPRIM,10200) NOMELT,NOMELE
10200 FORMAT(' ERREUR ELNMNU:ELEMENT FINI',2(1X,A4),' INCONNU PARMI'/
     &10(3X,A4,1X,A4))
      NO = 0
      RETURN
      END
