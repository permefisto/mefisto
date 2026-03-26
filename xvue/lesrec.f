      SUBROUTINE LESREC
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : CALCULER LA TAILLE, LA POSITION ET TRACER LES RECTANGLES
C ----  HISTORIQUE, MENU, INVITE, LIGNE DE SAISIE, LIGNES LUES
C       ET RESULTAT
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 1998
C23456---------------------------------------------------------------012
C
C     TRACE DU RESULTAT
      CALL LERESU
C
C     TRACE DU NOM DE L'OBJET
      CALL LHISTO
C
C     TRACE DU MENU
      CALL LEMENU
C
C     TRACE DES LIGNES LUES
      CALL LIGLUE
C
C     TRACE DE L'INVITE
      CALL LINVIT
C
C     TRACE LA LIGNE DE SAISIE
      CALL LIGSAI
      RETURN
      END
