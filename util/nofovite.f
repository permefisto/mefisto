      INTEGER FUNCTION NOFOVITE()
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  RETOURNER LE NO DE LA FONCTION 'VITESSE_EXACTE(t,x,y,z,nocomp)'
C -----  ou 'EXACT_VELOCITY' DE PLUS GRAND NUMERO
C        0 SI AUCUNE DE CES FONCTIONS N'EST TROUVEE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC & Saint Pierre du Perray   Mars 2010
C2345X7..............................................................012
      include"./incl/ntmnlt.inc"
C
      CALL LXNMNO( NTFONC, 'VITESSE_EXACTE', NOFOVITF, MNFONC )
      CALL LXNMNO( NTFONC, 'EXACT_VELOCITY', NOFOVITA, MNFONC )
C
C     LA PLUS GRANDE EST CENSEE ETRE LA PLUS RECENTE
      NOFOVITE = MAX( NOFOVITF, NOFOVITA )
C
      RETURN
      END
