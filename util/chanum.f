      SUBROUTINE CHANUM( KNO , NO )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    REALISER LA CONVERSION : CHAINE DE CARACTERES ==> ENTIER
C -----
C ENTREE :
C --------
C KNO    : CHAINE DE CARACTERES A CONVERTIR EN ENTIER
C
C SORTIE :
C --------
C NO     : ENTIER
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS  NOVEMBRE 1987
C.......................................................................
      INTEGER        NO
      CHARACTER*(*)  KNO
      CHARACTER*24   BUFFER
      CHARACTER*5    KFORMA
C
C     ELIMINATION DES BLANCS
      NC1 =0
C
 10   NC1 = NC1 + 1
      IF( KNO(NC1:NC1) .EQ. ' ' ) GOTO 10
C
C     RECHERCHE DU BLANC SUIVANT
      NCF = NC1
C
 20   NCF = NCF + 1
      IF( KNO(NCF:NCF) .NE. ' ' ) GOTO 20
C
C     LE NOMBRE S'ETEND DE NC1 A NCF SOIT NBC CHIFFRES
      NBC = NCF - NC1
C
C     CONVERSION DU NOMBRE DE CHIFFRES
      WRITE (BUFFER,'(I2)') NBC
C
C     LA SPECIFICATION DU FORMAT POUR CADRER JUSTE SANS BLANC
      KFORMA = '(I' // BUFFER(1:2) // ')'
C
C     CONVERSION DES CARACTERES EN NBC CHIFFRES
      BUFFER = KNO(NC1:NCF)
      READ (BUFFER,KFORMA) NO
      END
