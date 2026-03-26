      SUBROUTINE CARPNB( NL , NC )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : PASSER AU PREMIER CARACTERE SUIVANT NON BLANC
C -----
C
C ENTREES ET SORTIES :
C --------------------
C NL NC : NUMERO DE LIGNE ET COLONNE DU
C         EN ENTREE CARACTERE AU DELA DUQUEL CHERCHER
C                   LE PREMIER CARACTERE NON BLANC
C         EN SORTIE PREMIER CARACTERE NON BLANC
C         NL=0 SI DEBORDEMENT AU DELA DE LA LIGNE LHKLG DE KLG
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS  OCTOBRE 1989
C23456---------------------------------------------------------------012
      include"./incl/lu.inc"
C
C     LE CARACTERE SUIVANT
 1    NC = NC + 1
      IF( NC .GT. NBCALI ) THEN
C        PASSAGE A LA LIGNE SI POSSIBLE
         IF( NL .GE. LHKLG ) THEN
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
C     LE CARACTERE EST-IL UN BLANC ?
      IF( KLG(NL)(NC:NC) .EQ. ' ' ) GOTO 1
      END
