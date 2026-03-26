      FUNCTION MXCANB( NLMENU , KMENU )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCULER LE NUMERO MAXIMAL DE COLONNE
C ----- POSITION DU DERNIER CARACTERE NON BLANC DES LIGNES DE KMENU
C
C ENTREES :
C ---------
C NLMENU : NOMBRE DE LIGNES DE KMENU
C KMENU  : LES NLMENU LIGNES DE CARACTERES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS       MARS  1990
C23456---------------------------------------------------------------012
      CHARACTER*(*)  KMENU(NLMENU)
C
      MXCANB = 0
C     LE NOMBRE DE CARACTERES DECLARES DE KMENU
      NBC    = LEN( KMENU(1) )
C
C     LA BOUCLE SUR LES LIGNES
      DO 10 I=1,NLMENU
C        LA BOUCLE SUR LES CARACTERES FINAUX DE LA LIGNE I
         DO 5 J=NBC,1,-1
C           RECHERCHE DU DERNIER CARACTERE NON BLANC DE LA LIGNE I
            IF( KMENU(I)(J:J) .NE. ' ' ) GOTO 8
 5       CONTINUE
         J = 0
 8       MXCANB = MAX( MXCANB , J )
 10   CONTINUE
      END
