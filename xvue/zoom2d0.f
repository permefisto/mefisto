      SUBROUTINE ZOOM2D0( NOTYEV )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    INITIALISATION POUR UN TRACE 2D AVEC TRANSLATION ZOOM
C -----
C SORTIE :
C --------
C NOTYEV :    NO DE L'EVENEMENT
C          =0 SI ABANDON DEMANDE
C          =1 SI EXECUTION SANS ERREUR
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR:PERRONNET ALAIN ANALYSE NUMERIQUE LJLL UPMC PARIS NOVEMBRE 2003
C2345X...............................................................012
      include"./incl/trvari.inc"

C     POUR EVITER PLUSIEURS SUCCESSIVES EXECUTIONS
      IF( NORBITE .NE. 0 ) GOTO 9000

cccC     EN ATTENTE D'EVENEMENT
ccc      CALL XVSOURIS( NOTYEV, NOBOUTON, NOPX0, NOPY0 )
cccC
ccc      IF( NOTYEV .EQ. 0 ) THEN
cccC
cccC        ABANDON DEMANDE
ccc         NOBOUTON = 0
ccc         RETURN
cccC
ccc      ENDIF

C     FENETRE MAXIMALE POUR VOIR L'OBJET
      CALL ISOFENETRE( XOBMIN, XOBMAX, YOBMIN, YOBMAX )

      AXOPTV(1) = ( XOBMIN + XOBMAX ) / 2
      AXOPTV(2) = ( YOBMIN + YOBMAX ) / 2
      AXOPTV(3) = 0

      NOETATBOUTON = 0
      NOPX0 = -1
      NOPY0 = -1

C     LA MEMOIRE PIXELS EST EFFACEE
      CALL EFFACEMEMPX

C     TEMOIN DE PASSAGE DANS ZOOM2D0
      NORBITE = 1
 9000 NOTYEV  = 1

      RETURN
      END
