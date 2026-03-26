      INTEGER FUNCTION IND1BL( NOM )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    RETOURNER LE NUMERO DU 1ER BLANC SUIVI UNIQUEMENT DE BLANCS
C -----
C
C ENTREES :
C ---------
C NOM    : NOM SUIVI DE 0 OU PLUSIEURS BLANCS
C
C SORTIES :
C ---------
C IND1BL : NUMERO DU 1-ER BLANC SUIVI UNIQUEMENT DE BLANCS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS FEVRIER 1989
C23456---------------------------------------------------------------012
      CHARACTER*(*) NOM
C
C     LE NOMBRE DE CARACTERES DE NOM
      L = LEN( NOM )
C
C     REPERAGE DES BLANCS PAR LA FIN
      DO 100 IND1BL=L,1,-1
         IF( NOM(IND1BL:IND1BL) .NE. ' ' ) RETURN
 100  CONTINUE
      END
