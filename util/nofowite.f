      INTEGER FUNCTION NOFOWITE()
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  RETOURNER LE NO DE LA FONCTION 'WITESSE_EXACTE(t,x,y,z,nocomp)'
C -----  ou 'EXACT_WELOCITY' DE PLUS GRAND NUMERO
C        0 SI AUCUNE DE CES FONCTIONS N'EST TROUVEE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC  Saint Pierre du Perray Juillet 2011
C2345X7..............................................................012
      include"./incl/ntmnlt.inc"
C
      CALL LXNMNO( NTFONC, 'WITESSE_EXACTE', NOFOVITF, MNFONC )
      CALL LXNMNO( NTFONC, 'EXACT_WELOCITY', NOFOVITA, MNFONC )
C
C     LA PLUS GRANDE EST CENSEE ETRE LA PLUS RECENTE
      NOFOWITE = MAX( NOFOVITF, NOFOVITA )
C
      RETURN
      END
