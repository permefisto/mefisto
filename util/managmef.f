      SUBROUTINE MANAGMEF
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : MANAGEMENT GESTION de MEFISTO a l'AIDE d'UTILITAIRES sur
C ----- les TMS, FICHIERS, UNITES, ...
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : Alain PERRONNET  Saint PIERRE du PERRAY            AVRIL 2020
C2345X7..............................................................012

C     LECTURE DU MOT-CLE OU OPTION A EXECUTER PARMI TOUTES LES OPTIONS
C     CONTENUES DANS LE MENU DE NOM 'managmef'
C     ----------------------------------------------------------------
 1    CALL LIMTCL( 'managmef', NMTCL )

C     TRAITEMENT DE L'OPTION NMTCL
      IF( NMTCL .LT.  0 ) GOTO 9999

      GOTO( 10, 20, 30, 40, 9999 ), NMTCL

C     SUIVI_TMS
 10   CALL SUITMS
      GOTO 1

C     SUIVI DES FICHIERS DE LA MS
 20   CALL SUIFMS
      GOTO 1

C     UTILITAIRES DE GESTION DES UNITES DE LECTURE AFFICHAGE
 30   CALL SUIVES
      GOTO 1

C     DESTRUCTION DE TMS POINT, LIGNE, ..., OBJET, ...
 40   CALL TUER
      GOTO 1

C     Sortie
 9999 RETURN
      END


