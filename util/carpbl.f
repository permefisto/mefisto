      SUBROUTINE CARPBL( NL , NC )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   PASSER AU PROCHAIN CARACTERE BLANC
C -----
C
C ENTREES ET SORTIES :
C --------------------
C NL NC : NUMERO DE LIGNE ET COLONNE DU
C         EN ENTREE CARACTERE AU DELA DUQUEL CHERCHER
C                   LE PREMIER CARACTERE BLANC
C         EN SORTIE PREMIER CARACTERE BLANC
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS  OCTOBRE 1989
C23456---------------------------------------------------------------012
      include"./incl/lu.inc"
      include"./incl/gsmenu.inc"
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
C
C     LE CARACTERE SUIVANT
 1    NC = NC + 1
      IF( NC .GT. NBCALI ) THEN
C        PASSAGE A LA LIGNE SI POSSIBLE
         IF( NL .GE. MXKLG ) THEN
            NBLGRC(NRERR) = 1
            KERR(1) = 'LU:DEBORDEMENT DU NOMBRE DE LIGNES PERMISES'
            CALL LEREUR
            NL = 0
            NC = 0
            RETURN
         ENDIF
C        LA NOUVELLE LIGNE
         NL = NL + 1
         NC = 0
         GOTO 1
      ENDIF
C
C     LE CARACTERE EST-IL UN NON BLANC ?
      IF( KLG(NL)(NC:NC) .NE. ' ' ) GOTO 1
      END
