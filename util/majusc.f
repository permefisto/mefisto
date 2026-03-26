      SUBROUTINE MAJUSC( CHAINE )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : TRANSFORMER LES EVENTUELLES LETTRES MINUSCULES DE CHAINE
C ----- EN MAJUSCULES EN PASSANT PAR LE CODAGE DE LA TABLE ASCII
C
C ENTREE ET SORTIE :
C ------------------
C CHAINE : LA CHAINE DE CARACTERES A TRAITER
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN  ANALYSE NUMERIQUE UPMC PARIS    FEVRIER 1990
C23456---------------------------------------------------------------012
      CHARACTER*(*) CHAINE
      CHARACTER*1   LETMAJ
C
C     BOUCLE SUR LES CARACTERES
      DO 10 I = 1 , LEN( CHAINE )
         CHAINE(I:I) = LETMAJ( CHAINE(I:I) )
 10   CONTINUE
      END
