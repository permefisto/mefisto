      INTEGER FUNCTION NOFOPRES()
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  RETOURNER LE NO DE LA FONCTION 'PRESSION_EXACTE(t,x,y,z)' ou
C -----  'EXACT_PRESSURE' DE PLUS GRAND NUMERO
C        0 SI AUCUNE DE CES FONCTIONS N'EST TROUVEE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC & Saint Pierre du Perray   Mars 2010
C2345X7..............................................................012
      include"./incl/ntmnlt.inc"
C
      CALL LXNMNO( NTFONC, 'PRESSION_EXACTE', NOFOPRF, MNFONC )
      CALL LXNMNO( NTFONC, 'EXACT_PRESSURE',  NOFOPRA, MNFONC )
C
C     LA PLUS GRANDE EST CENSEE ETRE LA PLUS RECENTE
      NOFOPRES = MAX( NOFOPRF, NOFOPRA )
C
      RETURN
      END
