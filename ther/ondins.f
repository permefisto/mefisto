      SUBROUTINE ONDINS( KNOMOB, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LE DEPLACEMENT DU A DES ONDES INSTATIONNAIRES
C -----    SELON LES DIFFERENTES METHODES PROGRAMMEES
C          NOM DU SP avec INS POUR HOMOGENEITE AVEC THERMIQUE ET ELASTICITE
C
C ENTREE :
C --------
C KNOMOB : NOM DE L'OBJET D'ONDE INSTATIONNAIRE A TRAITER
C
C SORTIE :
C --------
C IERR   : 0 SI PAS D'ERREUR D'EXECUTION, NON NUL SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEURS: ALAIN PERRONNET  ANALYSE NUMERIQUE UPMC PARIS       MARS 1999
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/cnonlin.inc"
C
      CHARACTER*(*)  KNOMOB
C
C     LECTURE DES MOTS CLE
C     ====================
C     L'ANCIEN HISTORIQUE EST EFFACE
      CALL RECTEF( NRHIST )
      NBLGRC(NRHIST) = 1
      IF( LANGAG .EQ. 0 ) THEN
         KHIST(1) = 'OBJET: ' // KNOMOB
      ELSE
         KHIST(1) = 'OBJECT: ' // KNOMOB
      ENDIF
      CALL LHISTO
C
C     ONDE_INSTATIONNAIRE_LINEAIRE
C     COEFFICIENTS INDEPENDANTS DU TEMPS ET DE LA TEMPERATURE
      TESTNL = 0
      CALL ONDES( KNOMOB, IERR )
      GOTO 9999
c
ccc      CALL LIMTCL( 'resoonin', NMTCL )
ccc      IF( NMTCL .LE.  0 ) GOTO 9999
ccc      GOTO( 100, 200, 200 ), NMTCL
cccC
cccC     ONDE_INSTATIONNAIRE_LINEAIRE
cccC     COEFFICIENTS INDEPENDANTS DU TEMPS ET DE LA TEMPERATURE
ccc 100  TESTNL = 0
ccc      CALL ONDES( KNOMOB, IERR )
ccc      GOTO 9999
cccC
cccC     ONDE_INSTATIONNAIRE_NON_LINEAIRE
ccc 200  NBLGRC(NRERR) = 2
ccc      IF( LANGAG .EQ. 0 ) THEN
ccc         KERR(1) = 'OPTION NON PROGRAMMEE'
ccc         KERR(2) = 'OPTION COEFFICIENT INDEPENDANT du TEMPS FORCEE'
ccc      ELSE
ccc         KERR(1) = 'NOT PROGRAMMED OPTION'
ccc         KERR(2) = 'OPTION COEFFICIENT INDEPENDENT of TIME FORCED'
ccc      ENDIF
ccc      CALL LEREUR
ccc      GOTO 100
c
C
C     SORTIE
 9999 RETURN
      END
