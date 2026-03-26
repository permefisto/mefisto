      FUNCTION NUDCNB( TEXTE )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCULER LE NUMERO DU DERNIER CARACTERE NON BLANC D'UNE
C ----- CHAINE DE CARACTERES
C
C ENTREES :
C ---------
C TEXTE  : LA CHAINE DE CARACTERES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS       MARS  1990
C23456---------------------------------------------------------------012
      CHARACTER*(*)  TEXTE
C
C     LA BOUCLE SUR LES CARACTERES EN PARTANT DES DERNIERS
      L = LEN( TEXTE )
      DO 5 NUDCNB=L,1,-1
C        RECHERCHE DU DERNIER CARACTERE NON BLANC
         IF( TEXTE(NUDCNB:NUDCNB) .NE. ' ' ) RETURN
 5    CONTINUE
      NUDCNB = 0
      END
