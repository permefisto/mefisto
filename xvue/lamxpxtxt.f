      FUNCTION LAMXPXTXT( NLMENU , KMENU )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCULER LE NOMBRE MAXIMAL EN PIXELS D'UN TEXTE
C -----
C
C ENTREES :
C ---------
C NLMENU : NOMBRE DE LIGNES DE KMENU
C KMENU  : LES NLMENU LIGNES DE CARACTERES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     OCTOBRE 1996
C23456---------------------------------------------------------------012
      CHARACTER*(*)  KMENU(NLMENU)
C
      LAMXPXTXT = 0
C
C     LE NOMBRE DE CARACTERES DECLARES DE KMENU
      NBC = LEN( KMENU(1) )
C
C     LA BOUCLE SUR LES LIGNES
      DO 10 I=1,NLMENU
C        LA BOUCLE SUR LES CARACTERES FINAUX DE LA LIGNE I
         DO 5 J=NBC,1,-1
C           RECHERCHE DU DERNIER CARACTERE NON BLANC DE LA LIGNE I
            IF( KMENU(I)(J:J) .NE. ' ' ) GOTO 8
 5       CONTINUE
         J = 0
C        LE NOMBRE DE PIXELS EN LARGEUR DU TEXTE
 8       IF( J .GT. 0 ) THEN
            CALL XVNBPIXELTEXTE( KMENU(I)(1:J), J, NBPXLA, NBPXHA )
            LAMXPXTXT = MAX( LAMXPXTXT, NBPXLA )
         ENDIF
 10   CONTINUE
      END
