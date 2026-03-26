      INTEGER FUNCTION NOFODEEX()
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  RETOURNER LE NO DE FONCTION 'DEPLACEMENT_EXACT(t,x,y,z,nocomp)'
C -----  ou 'EXACT_DISPLACEMENT'  DE PLUS GRAND NUMERO
C        0 SI AUCUNE DE CES FONCTIONS N'EST TROUVEE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET Laboratoire J-L. LIONS UPMC PARIS Juillet 2007
C2345X7..............................................................012
      include"./incl/ntmnlt.inc"
C
      CALL LXNMNO( NTFONC, 'DEPLACEMENT_EXACT',   NOFOTEF, MNFONC )
      CALL LXNMNO( NTFONC, 'EXACT_DISPLACEMENT',  NOFOTEA, MNFONC )
C
C     LA PLUS GRANDE EST CENSEE ETRE LA PLUS RECENTE
      NOFODEEX = MAX( NOFOTEF, NOFOTEA )
C
      RETURN
      END
