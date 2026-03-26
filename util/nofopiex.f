      INTEGER FUNCTION NOFOPIEX()
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  RETOURNER LE NO DE LA FONCTION
C -----  PARTIE_IMAGINAIRE_EXACTE(t,x,y,z) ou EXACT_IMAGINARY_PART(t,x,y,z)
C        DE PLUS GRAND NUMERO
C        0 SI AUCUNE DE CES FONCTIONS N'EST TROUVEE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Texas A & M University at QATAR  Fevrier 2011
C2345X7..............................................................012
      include"./incl/ntmnlt.inc"
C
      CALL LXNMNO( NTFONC, 'PARTIE_IMAGINAIRE_EXACTE', NOFO1, MNFONC )
      CALL LXNMNO( NTFONC, 'EXACT_IMAGINARY_PART',     NOFO2, MNFONC )
C
C     LA PLUS GRANDE EST CENSEE ETRE LA PLUS RECENTE
      NOFOPIEX = MAX( NOFO1, NOFO2 )
C
      RETURN
      END
