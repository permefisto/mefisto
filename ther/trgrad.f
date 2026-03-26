      SUBROUTINE TRGRAD( NOPROJ, NCAS0,  NCAS1,  NDIM,  KNOMOB,  NTLXOB,
     %                   MODECO, NBTYEL, MNNPEF, MNPOGE, NDPGST,
     %                   TIMES )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER LES FLECHES REPRESENTANT LES GRADIENTS DE TEMPERATURE
C -----    AUX POINTS D'INTEGRATION NUMERIQUE DE TOUS LES EF
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
C          =1 CE SONT DEPLACEMENTS
C          =2 CE SONT DES MODES PROPRES
C NBTYEL : LE NOMBRE DE TYPES D'ELEMENTS DANS CETTE TOPOLOGIE
C MNTOPO : ADRESSE MCN DU TABLEAU TOPOLOGIE
C MNNPEF : ADRESSE MCN DU TABLEAU(NBTYEL) DES ADRESSES MCN DES TMS 'NPEF'
C MNPOGE : ADRESSE MCN DU TABLEAU POINTS GEOMETRIQUES DE L'OBJET
C NDPGST : CODE TRAITEMENT DES SOMMETS POINTS NOEUDS DU MAILLAGE
C          0 : NOEUDS=POINTS=SOMMETS
C          1 : NOEUDS=POINTS#SOMMETS
C          2 : NOEUDS#POINTS=SOMMETS
C          3 : NOEUDS#POINTS#SOMMETS
C TIMES  : TEMPS DU CALCUL DES NBVECT VECTEURS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 1994
C MODIFS : ALAIN PERRONNET Laboratoire JL LIONS UPMC Paris Decembre 2006
C MODIFS : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY Novembre 2010
C23456---------------------------------------------------------------012
      PARAMETER    ( LIGCON=0, LIGTIR=1 )
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/a_objet__topologie.inc"
      include"./incl/a___npef.inc"
      include"./incl/a___dtemperature.inc"
      include"./incl/a___xyzpoint.inc"
      include"./incl/ponoel.inc"
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/a___aretefr.inc"
      include"./incl/ctemps.inc"
      include"./incl/xyzext.inc"
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              GRADMX,GRADCM
      CHARACTER*(*)     KNOMOB
      REAL              TIMES(NCAS0:NCAS1)
C
      IERR   = 0
      NOPT   = 1
      CMPGRA = 0.0
      CMFLEC = 2.5
C     PAR DEFAUT LES FLECHES SONT ORANGES
      NCOUFL = NCORAN
      NOAXE  = 0
      NCAS   = NCAS1
C
C     RECHERCHE DU MIN ET MAX DES COORDONNEES DES POINTS
C     RECHERCHE DU MAXIMUM EN VALEUR ABSOLUE DES GRADIENTS
C     ====================================================
      CALL MXGRAD( NTLXOB, NBTYEL, MNNPEF,  COOEXT, GRADMX, IERR )
      IF( IERR .NE. 0 ) RETURN
      IF( GRADMX .LE. 0 ) GRADMX = 1.0
      GRADCM = GRADMX
C
C     LECTURE DES DONNEES POUR DEFINIR LA TAILLE DES FLECHES
C     ======================================================
 100  CALL LIMTCL( 'tracgrad', NMTCL )
      IF( NMTCL .LE.  0 ) GOTO 9000
C
      GOTO( 110, 120, 100, 100, 150, 160, 170, 100, 100, 100,
     %      100, 100, 100, 100, 300, 400 ) , NMTCL
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
C     1CM  VAUT EN GRADIENT
C     ---------------------
 120  NOPT   = 2
      NCVALS = 0
      CALL INVITE( 103 )
      CALL LIRRSP(NCVALS,GRADCM)
      IF( NCVALS .EQ. -1 ) GOTO 100
      IF( GRADCM .EQ. 0. ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'GRADIENT NUL a CORRIGER'
         ELSE
            KERR(1) = 'NULL GRADIENT to BE CORRECTED'
         ENDIF
         CALL LEREUR
         GOTO 100
      ENDIF
      GRADCM = ABS( GRADCM )
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
C     TRACE DES FLECHES GRADIENTS DANS TOUS LES EF
C     ============================================================
 300  CALL TRGREF( NOPROJ, KNOMOB, NTLXOB,
     %             NBTYEL, MNNPEF, MNPOGE, NDPGST,
     %             NCAS0,  NCAS1,  NDIM,   NMTCL,  MODECO, NOPT,
     %             GRADMX, GRADCM, CMFLEC, CMPGRA, TIMES )
      GOTO 100
C
C     TRACE DES FLECHES GRADIENTS DES EF 3D SECTIONNES PAR UN PLAN
C     POUR UNE SOLUTION NCAS
C     ============================================================
 400  IF( NDIM .LE. 2 ) GOTO 300
      CALL TRGREFSE( KNOMOB, NTLXOB, NOPROJ,
     %               NBTYEL, MNNPEF, NDPGST,
     %               MCN(MNPOGE+WBCOOP), MCN(MNPOGE+WNBPOI),
     %               MCN(MNPOGE+WYZPOI),
     %               MODECO, NCAS1, NOPT,
     %               GRADMX, GRADCM, CMFLEC, CMPGRA )
      GOTO 100
C
C     SORTIE SANS ERREUR
C     ==================
 9000 RETURN
      END
