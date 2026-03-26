      SUBROUTINE TRFLUX( NOPROJ, NCAS0,  NCAS1,  NDIM,   KNOMOB, NTLXOB,
     %                   MODECO, NBTYEL, MNELEM, MNPOGE, NDPGST,
     %                   TIMES )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER LES FLECHES REPRESENTANT LES FLUX NORMAUX DE CHALEUR
C -----    AUX POINTS D'INTEGRATION NUMERIQUE DES
C          FACES en 3D, ARETES en 2D, SOMMETS en 1D des EF
C
C ENTREES :
C ---------
C NOPROJ : SI OBJET EN 6D
C          TYPE DE PROJECTION 0 CI-DESSOUS FIXE LES COORDONNEES A ZERO
C          -1 PAS DE PROJECTION TRAITEMENT en XYZ NORMAL
C           1 : 'X Y Z 0 0 0'
C           2 : 'X Y 0 U 0 0'
C           3 : 'X 0 0 U V 0'
C           4 : '0 0 0 U V W'
C NCAS0  : NUMERO DU PREMIER JEU DE SOLUTION A TRACER
C NCAS1  : NUMERO DU DERNIER JEU DE SOLUTION A TRACER
C NDIM   : DIMENSION DE L'ESPACE DE L'OBJET ( 1 OU 2 OU 3 )
C KNOMOB : NOM DE L'OBJET
C NTLXOB : NUMERO DU TMS LEXIQUE DE L'OBJET A TRAITER
C MODECO : MODE DE TRACE DES VECTEURS DU TMS D'ADRESSE MNDEPL
C          =1 CE SONT SOLUTIONS
C          =2 CE SONT DES MODES PROPRES
C NBTYEL : LE NOMBRE DE TYPES D'ELEMENTS DANS CETTE TOPOLOGIE
C MNTOPO : ADRESSE MCN DU TABLEAU TOPOLOGIE
C MNELEM : ADRESSE MCN DES TABLEAUX ELEMENTS DE CETTE TOPOLOGIE
C MNPOGE : ADRESSE MCN DU TABLEAU POINTS GEOMETRIQUES DE L'OBJET
C NDPGST : CODE TRAITEMENT DES SOMMETS POINTS NOEUDS DU MAILLAGE
C          0 : NOEUDS=POINTS=SOMMETS
C          1 : NOEUDS=POINTS#SOMMETS
C          2 : NOEUDS#POINTS=SOMMETS
C          3 : NOEUDS#POINTS#SOMMETS
C TIMES  : TEMPS DU CALCUL DES NBVECT VECTEURS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    Novembre 1994
C MODIFS : ALAIN PERRONNET Laboratoire JL LIONS UPMC Paris Decembre 2006
C MODIFS : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY Novembre 2010
C23456---------------------------------------------------------------012
      PARAMETER     ( LIGCON=0, LIGTIR=1 )
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/a_objet__topologie.inc"
      include"./incl/a___npef.inc"
      include"./incl/a___xyzpoint.inc"
      include"./incl/ponoel.inc"
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/a___fluxpt.inc"
      include"./incl/a___aretefr.inc"
      include"./incl/ctemps.inc"
      include"./incl/xyzext.inc"
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL             RMCN(1)
      EQUIVALENCE     (RMCN(1),MCN(1))
      REAL              FLUXMX,FLUXCM
      CHARACTER*(*)     KNOMOB
      REAL              TIMES(NCAS0:NCAS1)
C
      IERR   = 0
      NOPT   = 1
      CMPFLU = 0.0
      CMFLEC = 2.5
C     PAR DEFAUT LES FLECHES SONT MAGENTA
      NCOUFL = NCMAGE
      NOAXE  = 0
C
C     RECHERCHE DU MIN ET MAX DES COORDONNEES DES POINTS
C     RECHERCHE DU MAXIMUM EN VALEUR ABSOLUE DES FLUX NORMAUX
C     =======================================================
      CALL MXFLUX( NTLXOB, NBTYEL, MNELEM,  COOEXT, FLUXMX, IERR )
      IF( IERR .NE. 0 ) RETURN
      IF( FLUXMX .LE. 0 ) FLUXMX = MAX( FLUXMX, 1.0 )
      FLUXCM = FLUXMX
C
C     LECTURE DES DONNEES POUR DEFINIR LA TAILLE DES FLECHES
C     ======================================================
 100  CALL LIMTCL( 'tracflux', NMTCL )
      IF( NMTCL .LE.  0 ) GOTO 9000
C
      GOTO( 110, 120, 100, 100, 150, 160, 170, 100, 100, 100,
     %      100, 100, 100, 100, 400, 500 ) , NMTCL
C
C     NOMBRE DE CM POUR TRACER LA FLECHE MAXIMALE
C     -------------------------------------------
 110  NOPT   = 1
      NCVALS = 0
      CALL INVITE( 66 )
      CALL LIRRSP( NCVALS, CMFLEC )
      IF( NCVALS .EQ. -1 ) GOTO 100
C     CMFLEC : NOMBRE DE CM DU DE LA FLECHE MAXIMALE
      IF( CMFLEC .LT. 0. ) CMFLEC = -CMFLEC
      GOTO 100
C
C     1CM  VAUT EN FLUX
C     -----------------
 120  NOPT   = 2
      NCVALS = 0
      CALL INVITE( 101 )
      CALL LIRRSP(NCVALS,FLUXCM)
      IF( NCVALS .EQ. -1 ) GOTO 100
      IF( FLUXCM .EQ. 0. ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'FLUX NORMAL NUL A CORRIGER'
         ELSE
            KERR(1) = 'NULL NORMAL FLUX. TO BE CORRECTED'
         ENDIF
         CALL LEREUR
         GOTO 100
      ENDIF
      FLUXCM = ABS( FLUXCM )
      GOTO 100
C
C     COULEUR des ARETES du MAILLAGE
C     ------------------------------
 150  CALL LIMTCL( 'couleur0' , I )
      IF( I .EQ. -1 ) THEN
         GOTO 9000
      ELSE IF( I .EQ. -2 ) THEN
         NCOUAF = -2
      ELSE IF( I .EQ. 0 ) THEN
C        COULEUR NOIRE
         NCOUAF = 0
      ELSE
C        COULEUR RESERVEE
         NCOUAF = N1COEL + I
      ENDIF
      GOTO 100
C
C     TYPE du TRAIT des ARETES du MAILLAGE
C     ------------------------------------
 160  CALL LIMTCL( 'typtrait' , I )
      IF( I .EQ. -1 ) GOTO 100
      NTLAFR = I
      GOTO 100
C
C     COULEUR des FLECHES
C     -------------------
 170  CALL LIMTCL( 'couleurs' , I )
      IF( I .LT. 0 ) THEN
         GOTO 9000
      ELSE IF( I .EQ. 0 ) THEN
C        COULEUR NOIRE
         NCOUFL = 0
      ELSE
         NCOUFL = N1COEL + I
      ENDIF
      GOTO 100
C
C
C     TRACE DES FLECHES FLUX NORMAUX DE TEMPERATURE DANS TOUS LES EF
C     ==============================================================
 400  CALL TRFLEF( NOPROJ, KNOMOB, NTLXOB,
     %             NBTYEL, MNELEM, MNPOGE, NDPGST,
     %             NCAS0,  NCAS1,  NDIM,   NMTCL,  MODECO, NOPT,
     %             FLUXMX, FLUXCM, CMFLEC, CMPFLU, TIMES )
      GOTO 100
C
C     TRACE DES FLECHES FLUX NORMAUX DES EF 3D SECTIONNES PAR UN PLAN
C     POUR UNE SOLUTION NCAS
C     ===============================================================
 500  IF( NDIM .LE. 2 ) GOTO 400
C
C     MISE A JOUR DE L'ECHELLE DE FLUX
C     CMPFLU : NOMBRE DE CM POUR L'UNITE DE FLUX
      IF( NOPT .EQ. 1 ) THEN
         CMPFLU = CMFLEC / FLUXMX
      ELSE IF( NOPT .EQ. 2 ) THEN
         CMPFLU = 1. / FLUXCM
      ENDIF
      CALL TRFLEFSE( KNOMOB, NTLXOB, NOPROJ,
     %               NBTYEL, MNELEM, NDPGST,
     %               MCN(MNPOGE+WBCOOP), MCN(MNPOGE+WNBPOI),
     %               RMCN(MNPOGE+WYZPOI),
     %               MODECO, NCAS1, CMPFLU )
      GOTO 100
C
C     SORTIE SANS ERREUR
C     ==================
 9000 RETURN
      END
