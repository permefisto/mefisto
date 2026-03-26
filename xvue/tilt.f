      SUBROUTINE TILT
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  ARRET EN CAS DE PROBLEME POUR DETERMINER SA POSITION
C -----
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS         MAI 1994
C2345X7..............................................................012
      include"./incl/gsmenu.inc"
C
      NBLGRC(NRERR) = 3
      KERR(1) = 'TILT: ERREUR A DEBOGUER'
      KERR(2) = 'RELANCER APRES AVOIR CORRIGER L''ERREUR'
      KERR(3) = 'ATTENTION: APRES LE CLIC SOURIS => ARRET'
      CALL LEREUR
C
C     SAISIE D'UN POINT PAR CLIC DE LA SOURIS OU ENTREE D'UN CARACTERE
      CALL SAIPTC( NOTYEV, NX, NY, NOCHAR )
C
C     ARRET DEFINITIF
      STOP 'TILT du logiciel MEFISTO'
      END
