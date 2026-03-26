      INTEGER FUNCTION NOFOTEEX()
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  RETOURNER LE NO DE LA FONCTION 'TEMPERATURE_EXACTE' ou 'EXACT_TEMPERATU
C -----  DE PLUS GRAND NUMERO
C        0 SI AUCUNE DE CES FONCTIONS N'EST TROUVEE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET Laboratoire J-L. LIONS UPMC PARIS Juillet 2007
C2345X7..............................................................012
      include"./incl/ntmnlt.inc"
C
      CALL LXNMNO( NTFONC, 'TEMPERATURE_EXACTE', NOFOTEF, MNFONC )
      CALL LXNMNO( NTFONC, 'EXACT_TEMPERATURE',  NOFOTEA, MNFONC )
C
C     LA PLUS GRANDE EST CENSEE ETRE LA PLUS RECENTE
      NOFOTEEX = MAX( NOFOTEF, NOFOTEA )
C
      RETURN
      END
