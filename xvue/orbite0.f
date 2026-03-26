      SUBROUTINE ORBITE0( NOTYEV )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    INITIALISATION POUR UN TRACE 3D AVEC TRANSLATION ORBITE ZOOM
C -----
C SORTIE :
C --------
C NOTYEV : NO DE L'EVENEMENT =0 SI ABANDON DEMANDE, NON NUL SINON
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR:PERRONNET ALAIN ANALYSE NUMERIQUE LJLL UPMC PARIS NOVEMBRE 2003
C2345X...............................................................012
      include"./incl/trvari.inc"
      include"./incl/xyzext.inc"
C
C     POUR EVITER PLUSIEURS SUCCESSIVES EXECUTIONS
      IF( NORBITE .NE. 0 ) RETURN
C
ccc      print*,'orbite0 1: cooext X=',cooext(1,1),cooext(1,2),
ccc     %       '  Y=',cooext(2,1),cooext(2,2),
ccc     %       '  Z=',cooext(3,1),cooext(3,2)

      DIAMOB = SQRT( (COOEXT(1,2)-COOEXT(1,1))**2
     %             + (COOEXT(2,2)-COOEXT(2,1))**2
     %             + (COOEXT(3,2)-COOEXT(3,1))**2 )
C
cccC     EN ATTENTE D'EVENEMENT
ccc      CALL XVSOURIS( NOTYEV, NOBOUTON, NOPX, NOPY )
cccC
ccc      IF( NOTYEV .EQ. 0 ) THEN
cccC
cccC        ABANDON DEMANDE
ccc         NOBOUTON = 2
ccc         RETURN
cccC
ccc      ENDIF
C
C     LONGITUDE ET LATITUDE ACTUELLE EN RADIANS
      CALL LOLARA( AXOPTV(1), AXOPTV(2), AXOPTV(3),
     %             AXOEIL(1), AXOEIL(2), AXOEIL(3),
     %             RLONGITUDE0, RLATITUDE0, LEPOLE )
C
C     BORNES EN LARGEUR ET HAUTEUR 2D DE LA FENETRE DE PROJECTION
      AXLMIN = AXOLAR / 1000
      AXLMAX = AXOLAR * 3
      AXHMIN = AXOHAU / 1000
      AXHMAX = AXOHAU * 3
C
C     SAUVEGARDE DANS LE COMMON / TRVAR5 /
C     NOBOUTON: =1 a 3 NUMERO DU BOUTON ACTIF DE LA SOURIS
C               =0 SI AUCUN ACTIF
      NOBOUTON = 0
C     NOETATBOUTON: 1 SI LE BOUTON NON NUL ACTUEL EST ENFONCE NON RELACHE
C                   0 SI LE BOUTON NON NUL ACTUEL A ETE ENFONCE PUIS RELACHE
      NOETATBOUTON = 0
C     NOPX0, NOPY0: NUMERO EN X ET Y DU PIXELS
C     VALEUR INDIQUANT UNE NON INITIALISATION
      NOPX0 = -1
      NOPY0 = -1
C
C     TEMOIN DE PASSAGE DANS ORBITE0
C     NORBITE    : 0 AU DEPART
C                  1 APRES EXECUTION DE ZOOM[12]D0 ou ORBITE0
C                >=2 APRES EXECUTION D'AU MOINS UNE FOIS  ZOOM[12]D1 ou ORBITE1
C                  0 A L'ABANDON DE  ZOOM[12]D1 ou ORBITE1
      NORBITE = 1
      NOTYEV  = 1
C
cccC     LA MEMOIRE PIXELS EST EFFACEE
ccc      CALL EFFACEMEMPX
C
ccc      print*,'orbite0 2: cooext X=',cooext(1,1),cooext(1,2),
ccc     %       '  Y=',cooext(2,1),cooext(2,2),
ccc     %       '  Z=',cooext(3,1),cooext(3,2)
C
      RETURN
      END
