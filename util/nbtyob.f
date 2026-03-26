      INTEGER FUNCTION NBTYOB( NOTY , NOOB )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     CODER DANS UN ENTIER LE NUMERO DE TYPE DE L'OBJET
C -----     ET SON NUMERO DANS LE LEXIQUE
C
C ENTREE :
C --------
C NOTY   : NUMERO DU TYPE DE L'OBJET OU SOUS OBJET (CF NTYOBS)
C           1:'POINT'   2:'LIGNE'  3:'SURFACE'  4:'VOLUME' 5:'OBJET'
C          -1:'SOMMET' -2:'ARETE' -3:'FACE'    -4:'CUBE'
C NOOB   : NUMERO DE L'OBJET
C
C SORTIE :
C --------
C NBTYOB : NUMERO CODE SOUS LA FORME
C          = SIGNE(NOTY) * ( NOOB + 1 000 000 * ABS(NOTY) )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS   AVRIL 1987
C.....................................................................
C
C     LE CODAGE DU NUMERO DE TYPE ET NUMERO DANS LE LEXIQUE
      IF( NOTY .LE. 6 ) THEN
C        TYPE POINT OU LIGNE OU SURFACE OU VOLUME OU OBJET
         NBTYOB = NOOB + NOTY * 1 000 000
      ELSE
C        TYPE SOMMET OU ARETE OU FACE OU CUBE
         NBTYOB = - ( NOTY * 1 000 000 + NOOB )
      ENDIF
      END
