      SUBROUTINE UNSEBL( TEXTE , NUDCAR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :REECRIRE LE TEXTE AVEC AU PLUS UN BLANC LA OU IL Y EN A PLUSIEURS
C -----AUCUN BLANC EN DEBUT DE LIGNE
C
C ENTREE ET SORTIE :
C ------------------
C TEXTE  : LA CHAINE DE CARACTERES
C
C SORTIE :
C --------
C NUDCAR : POSITION DU DERNIER CARACTERE NON BLANC EN SORTIE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     OCTOBRE 1990
C23456---------------------------------------------------------------012
      CHARACTER*(*)  TEXTE
C
C     LE DERNIER CARACTERE NON BLANC
      NUDCAR  = NUDCNB( TEXTE )
C
C     N1 LE PREMIER CARACTERE NON BLANC
      N1 = 1
      DO 10 N=1,NUDCAR
         IF( TEXTE(N:N) .EQ. ' ' ) THEN
            N1 = N1 + 1
         ELSE
            GOTO 20
         ENDIF
 10   CONTINUE
C
 20   IF( N1 .GT. 1 ) THEN
         TEXTE  = TEXTE(N1:NUDCAR)
         NUDCAR = NUDCAR - N1 + 1
      ENDIF
C
      N1 = 1
C
 30   N  = INDEX( TEXTE(N1:NUDCAR) , '  ' )
      IF( N .GT. 0 ) THEN
C        LE SECOND CARACTERE BLANC EST SUPPRIME
         N2     = N1 + N - 2
         TEXTE(N2+1:NUDCAR) = TEXTE(N2+2:NUDCAR)
         N1     = N2 + 2
         NUDCAR = NUDCAR - 1
         GOTO 30
      ENDIF
C
C     MISE A BLANC DES CARACTERES TRANSLATES
      N = LEN( TEXTE )
      IF( NUDCAR .LT. N ) TEXTE(NUDCAR+1:N) = ' '
      END
