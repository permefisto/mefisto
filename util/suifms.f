      SUBROUTINE SUIFMS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : SUIVI DES FICHIERS DE LA MEMOIRE SECONDAIRE
C -----
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS        MARS 1990
C2345X7..............................................................012
C
C     LECTURE DU MOT CLE A TRAITER
C     ----------------------------
 10   CALL LIMTCL( 'suivi_ms' , NMTCL )
      IF( NMTCL .LE. 0 ) GOTO 9000
C
C     AJOUT_UN_FICHIER
C     ================
      CALL AJFICH( NCVALS )
      IF( NCVALS .EQ. -1 ) GOTO 9000
C     FERMETURE DE LA MS POUR SE DONNER LA POSSIBILITE DE
C     DECLARER UN NOUVEAU FICHIER
C     CALCUL DE NMTCL POUR TRAITER UNE ERREUR SUR APOLLO
      NMTCL = - MAX(1,NMTCL)
      CALL ARRET( NMTCL )
      GOTO 10
C
 9000 RETURN
      END
