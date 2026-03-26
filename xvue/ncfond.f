       INTEGER FUNCTION NCFOND()
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : DEFINIR LA COULEUR INITIALE DU FOND DE FENETRE MEFISTO
C -----
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC  PARIS  SEPTEMBRE 2006
C2345X7..............................................................012
      include"./incl/trvari.inc"
C
C     LES 16 COULEURS POSSIBLES
C     POUR LE FOND CHOISIR UNE COULEUR CLAIRE ET DIFFERENTE DE CELLE
C     DU TRACE DES POINTS LIGNES SURFACES ...
C     NCNOIR, NCROUG, NCVERT, NCBLEU, NCCYAN, NCJAUN, NCMAGE, NCBLAN,
C     NCGRIS, NCGRIM, NCGRIC, NCBEIG, NCORAN, NCSAUM, NCROSE, NCTURQ
C
ccc      NCFOND = NCGRIM
      NCFOND = NCBEIG

      RETURN
      END
