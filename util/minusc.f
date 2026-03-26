      SUBROUTINE MINUSC( CHAINE )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : TRANSFORMER LES EVENTUELLES LETTRES MAJUSCULES DE CHAINE
C ----- EN MINUSCULES EN PASSANT PAR LE CODAGE DE LA TABLE ASCII
C
C ENTREE ET SORTIE :
C ------------------
C CHAINE : LA CHAINE DE CARACTERES A TRAITER
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS    FEVRIER 1990
C23456---------------------------------------------------------------012
      CHARACTER*(*) CHAINE
      CHARACTER*1   LETMIN
C
C     BOUCLE SUR LES CARACTERES
      DO 10 I = 1 , LEN( CHAINE )
         CHAINE(I:I) = LETMIN( CHAINE(I:I) )
 10   CONTINUE
      END
