      INTEGER FUNCTION NOFOTIEL()
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  RETOURNER le NUMERO de la FONCTION 'TAILLE_IDEALE(x,y,z)'
C -----            ou 'EDGE_LENGTH(x,y,z)' de PLUS GRAND NUMERO dans
C                  le TMS des FONCTIONS CREEES par l'UTILISATEUR
C        =0 SI AUCUNE DE CES 2 FONCTIONS N'EST TROUVEE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET Laboratoire J-L. LIONS UPMC PARIS Juillet 2007
C MODIFS: ALAIN PERRONNET Saint Pierre du Perray                Mai 2020
C2345X7..............................................................012
      include"./incl/ntmnlt.inc"
      include"./incl/darete.inc"

      CALL LXNMNO( NTFONC, 'TAILLE_IDEALE', NOFOTIF, MNFONC )
      CALL LXNMNO( NTFONC, 'EDGE_LENGTH',   NOFOTIA, MNFONC )

C     LA PLUS GRANDE EST CENSEE ETRE LA PLUS RECENTE
      NOFOTIEL = MAX( NOFOTIF, NOFOTIA )

C     SAUVEGARDE DANS LE COMMON /DEARET /
      NOFOTI = NOFOTIEL

      RETURN
      END
