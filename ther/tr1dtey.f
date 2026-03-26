      SUBROUTINE TR1DTEY( KNOMOB, MODECO, NDIM,
     %                    NBTYEL, MNELEM, MNPOGE, NDPGST,
     %                    NCAS0,  NCAS1,  NTDL,   TEMPER,
     %                    TMIN,   NOEMIN, NCAMIN, TMAX,  NOEMAX, NCAMAX,
     %                    TEMPSS )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  TRACE DE LA TEMPERATURE SOUS FORME D'UNE COURBE Y=TEMPERATURE(t,x)
C -----  AVEC ZONES DE COULEURS ENTRE L'AXE ET LA COURBE
C
C ENTREES:
C --------
C KNOMOB : NOM DE L'OBJET
C MODECO : MODE DE TRACE DES VECTEURS DU TMS D'ADRESSE MNDEPL
C          =1 CE SONT DES ISO-SOLUTIONS
C          =2 CE SONT DES MODES PROPRES
C          =3 CE SONT DES PRESSIONS P1 A PARTIR D'UNE INTERPOLATION
C             SOIT en 2D P2 SOIT P1+BULLE P3 OU en 3D P1 OU en 3D P2
C          =4 CE SONT DES ERREURS AUX NOEUDS D'UN MAILLAGE
C          =5 CE SONT DES SOLUTIONS
C NDIM   : DIMENSION DE L'OBJET (1 OU 2 OU 3)
C NBTYEL : LE NOMBRE DE TYPES D'ELEMENTS DANS CET OBJET
C MNELEM : ADRESSE MCN DES TABLEAUX ELEMENTS DE CET OBJET
C MNPOGE : ADRESSE MCN DU TABLEAU POINTS GEOMETRIQUES DE L'OBJET
C NDPGST : CODE TRAITEMENT DES SOMMETS POINTS NOEUDS DU MAILLAGE
C          0 : NOEUDS=POINTS=SOMMETS
C          1 : NOEUDS=POINTS#SOMMETS
C          2 : NOEUDS#POINTS=SOMMETS
C          3 : NOEUDS#POINTS#SOMMETS
C NCAS0  : NUMERO DU PREMIER JEU DE SOLUTION A TRACER
C NCAS1  : NUMERO DU DERNIER JEU DE SOLUTION A TRACER
C NTDL   : NOMBRE TOTAL DE DEGRES DE LIBERTE
C TEMPER : TEMPERATURE(NTDL,NCAS0:NCAS1) LES NCAS0:NCAS1 VECTEURS TEMPERATURE
C TMIN   : TEMPERATURE MINIMALE DU CAS NCAS
C NOEMIN : NUMERO DU NOEUD OU LA TEMPERATURE EST MINIMALE
C TMAX   : TEMPERATURE MAXIMALE DU CAS NCAS
C NOEMAX : NUMERO DU NOEUD OU LA TEMPERATURE EST MAXIMALE
C TEMPSS : TEMPS DES NCAS0:NCAS1 VECTEURS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR:ALAIN PERRONNET LJLL UPMC PARIS & St Pierre du Perray JUIN 2009
C MODIFS:ALAIN PERRONNET LJLL UPMC & St Pierre du Perray   Novembre 2010
C23456---------------------------------------------------------------012
      PARAMETER   ( LIGCON=0, LIGTIR=1 )
      IMPLICIT INTEGER (W)
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/a___npef.inc"
      include"./incl/a___vecteur.inc"
      include"./incl/a___xyzpoint.inc"
      include"./incl/a___xyznoeud.inc"
      include"./incl/ponoel.inc"
      include"./incl/donele.inc"
      include"./incl/trvari.inc"
      include"./incl/xyzext.inc"
      include"./incl/mecoit.inc"
      include"./incl/ctemps.inc"
C
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
      CHARACTER*(*)     KNOMOB
      DOUBLE PRECISION  TEMPER(NTDL,NCAS0:NCAS1)
C
      CHARACTER*100     KNOM
      REAL              TEMPSS(NCAS0:NCAS1)
      REAL              X(3), TEMP(3), COULN(3), XY(2,4), COULS(4)
      INTEGER           NONOEF(3)
C
      MOREE2 = MOTVAR(6)
C
C     NOMBRE DE COULEURS DISPONIBLES
      NBCOUL = NDCOUL - N1COUL + 1
C
C     LA PALETTE ARC EN CIEL
      CALL PALCDE(11)
C
C     LES ARETES SONT TRACEES AVEC UNE EPAISSEUR
      CALL XVEPAISSEUR( 1 )
C
C     NBPOI  NOMBRE DE POINTS DU MAILLAGE DE L'OBJET
      NBPOI  = MCN(MNPOGE+WNBPOI)
C
C     NBCOOR NOMBRE DE COORDONNEES D'UN POINT
      NBCOOR = MCN(MNPOGE+WBCOOP)
C
C     MIN ET MAX DES COORDONNEES DES POINTS DU MAILLAGE DE L'OBJET
      CALL MAJEXT( MNPOGE )
      CALL MIMXPT( NBCOOR, NBPOI, RMCN(MNPOGE+WYZPOI), COOEXT )
      XMIN = COOEXT(1,1)
      XMAX = COOEXT(1,2)
C
C     LA FENETRE DE TRACE EST COMPLETEE EN ORDONNEE
      COOEXT(2,1) = TMIN
      COOEXT(2,2) = TMAX
      COOEXT(3,1) = 0.0
      COOEXT(3,2) = 0.0
C
C     DEFINITION DES MARGES DU TRACE
      HX = ( COOEXT(1,2) - COOEXT(1,1) ) / 40
      TH = ( TMAX - TMIN ) / 20
      CALL FENETRE( COOEXT(1,1)-HX, COOEXT(1,2)+HX, TMIN-TH, TMAX+TH )
      TH2 = TH / 2
C     LE TYPE DE LA VISEE
      NDIMLI = 1
      NOTYVI = 1
      CALL XVEPAISSEUR( 2 )
C
C     -------------------
C     OPTIONS DE LA VISEE
C     -------------------
 10   CALL VISE1D( NMTCL )
      IF( NMTCL .LT. 0 ) RETURN
C     INITIALISATION DE L'ORBITE ZOOM DEPLACEMENT
      IF( LORBITE .NE. 0 ) THEN
         CALL ZOOM1D0( NOTYEV )
         IF( NOTYEV .EQ. 0 ) GOTO 10
      ENDIF
C
 20   DO 1000 NCAS=NCAS0, NCAS1
C
C        L'ECRAN PIXELS EST EFFACE (PAS DE SCINTILLEMENT)
         CALL EFFACEMEMPX
C
C        LE TEMPS DU VECTEUR SOLUTION NCAS
         TEMPS = TEMPSS( NCAS )
C
C        LE TRACE DES AXES 2D
         CALL TRAXE2
C
C        LA SIGNIFICATION DES AXES
         CALL TEXTE2D(NCNOIR, COOEXT(1,1),    TMAX,'Solution')
         CALL TEXTE2D(NCNOIR, COOEXT(1,2)-HX, TMIN, 'X' )
C
C        LA COULEUR DE LA TEMPERATURE TMIN
         COULMIN = N1COUL
C
C        BOUCLE SUR LES DIFFERENTS TYPES D'ELEMENTS FINIS DU MAILLAGE
C        ============================================================
         DO 90 NTEF = 0, NBTYEL-1
C
C           L'ADRESSE MCN DU DEBUT DU TABLEAU NPEF DE CE TYPE D'EF
            MNELE = MCN( MNELEM + NTEF )
C
C           LE NUMERO DU TYPE DES ELEMENTS FINIS
            NUTYEL = MCN( MNELE + WUTYEL )
C
C           LE NOMBRE DE TELS ELEMENTS
            NBELEM = MCN(MNELE + WBELEM )
C
C           LES CARACTERISTIQUES DE L'ELEMENT FINI
            CALL ELTYCA( NUTYEL )
C
C           L'ADRESSE MCN DU TABLEAU 'NPEF' POUR CE TYPE D'EF
            MNPGEL = MNELE + WUNDEL
            IF( NDPGST .GE. 2 ) THEN
               MNPGEL = MNPGEL + NBELEM * MCN(MNELE+WBNDEL)
            ENDIF
C
            DO 50 NUELEM = 1 , NBELEM
C
C              LE NUMERO DES NOEUDS DE L'EF
               CALL EFNOEU( MNELE, NUELEM,  NBNOEF, NONOEF )
C
C              TRACE DU SEGMENT SOUS FORME D'UN TRAPEZE COLORE
               CALL XVTYPETRAIT( NTLAPL )
C
C              LA COULEUR DES NBNOEF NOEUDS
               DO 30 I=1,NBNOEF
C
C                 LA COORDONNEE X DU NOEUD I
                  X(I) = RMCN(MNPOGE+WYZPOI+(NONOEF(I)-1)*NBCOOR)
C
C                 LA TEMPERATURE AU NOEUD I
                  TEMP(I) = REAL( TEMPER( NONOEF(I), NCAS ) )
C
C                 LA COULEUR AU NOEUD I
                  IF( TMIN .NE. TMAX ) THEN
                     COULN(I)=(TEMP(I)-TMIN)/(TMAX-TMIN)*NBCOUL+N1COUL
                  ELSE 
                     COULN(I)=NDCOUL
                  ENDIF
                  IF( COULN(I) .GT. NDCOUL ) COULN(I)=NDCOUL
                  IF( COULN(I) .LT. N1COUL ) COULN(I)=N1COUL
C
 30            CONTINUE
C
               IF( NBNOEF .EQ. 2 ) THEN
C
C                 NOEUDS = 2 SOMMETS
C                 LE SOMMET 1
                  XY(1,1) = X(1)
                  XY(2,1) = TMIN
                  COULS(1)=COULMIN
C                 LE SOMMET 2
                  XY(1,2) = X(2)
                  XY(2,2) = TMIN
                  COULS(2)=COULMIN
C                 LE SOMMET 3
                  XY(1,3) = X(2)
                  XY(2,3) = TEMP(2)
                  COULS(3) = COULN(2)
C                 LE SOMMET 4
                  XY(1,4) = X(1)
                  XY(2,4) = TEMP(1)
                  COULS(4) = COULN(1)
                  CALL  QUADCOUL2D( XY, COULS )
                  NS2 = 2
C
               ELSE
C
C                 NOEUDS = 2 SOMMETS SUIVI DU MILIEU DU SEGMENT
C
C                 NOEUDS = 2 SOMMETS + QUADRANGLE 13 et 32
C                 LE SOMMET 1 = 1-er SOMMET
                  XY(1,1) = X(1)
                  XY(2,1) = TMIN
                  COULS(1)=COULMIN
C                 LE SOMMET 2 = MILIEU
                  XY(1,2) = X(3)
                  XY(2,2) = TMIN
                  COULS(2)=COULMIN
C                 LE SOMMET 3 = MILIEU
                  XY(1,3) = X(3)
                  XY(2,3) = TEMP(3)
                  COULS(3) = COULN(3)
C                 LE SOMMET 4 = 1-er SOMMET
                  XY(1,4) = X(1)
                  XY(2,4) = TEMP(1)
                  COULS(4) = COULN(1)
                  CALL  QUADCOUL2D( XY, COULS )
C
C                 LE SOMMET 1 = MILIEU
                  XY(1,1) = X(3)
                  XY(2,1) = TMIN
                  COULS(1)=COULMIN
C                 LE SOMMET 2 = 2-eme SOMMET
                  XY(1,2) = X(2)
                  XY(2,2) = TMIN
                  COULS(2)=COULMIN
C                 LE SOMMET 3 = 2-eme SOMMET
                  XY(1,3) = X(2)
                  XY(2,3) = TEMP(2)
                  COULS(3) = COULN(2)
C                 LE SOMMET 4 = MILIEU
                  XY(1,4) = X(3)
                  XY(2,4) = TEMP(3)
                  COULS(4) = COULN(3)
                  CALL  QUADCOUL2D( XY, COULS )
                  NS2 = 3
C
               ENDIF
C
C              TRACE DU SYMBOLE DES 2 SOMMETS ET DE L'ARETE
               XY(2,  1) = TMIN - TH2
               XY(2,NS2) = TMIN - TH2
               CALL SYMBOLE2D( NCGRIS, XY(1,1),   XY(2,1),   'I' )
               CALL SYMBOLE2D( NCGRIS, XY(1,NS2), XY(2,NS2), 'I' )
               CALL TRAIT2D(NCGRIS,XY(1,1),XY(2,1), XY(1,NS2),XY(2,NS2))
C
 50         CONTINUE
C
 90      CONTINUE
C
C        MINIMUM et MAXIMUM sont TRACES
         IF( NCAMIN .EQ. NCAS ) THEN
            MN = MNPOGE+WYZPOI+3*NOEMIN-3
            CALL SYMBOLE2D( NCGRIM, RMCN(MN), RMCN(MN+1),'.m' )
         ENDIF
         IF( NCAMAX .EQ. NCAS ) THEN
            MN = MNPOGE+WYZPOI+3*NOEMAX-3
            CALL SYMBOLE2D( NCGRIM, RMCN(MN), RMCN(MN+1),'.M' )
         ENDIF
C
C        RETOUR AU TRACE CONTINU DES LIGNES
         CALL XVTYPETRAIT( LIGCON )
C
C        EFFACEMENT DE LA LEGENDE SUR POSTSCRIPT
         IF ( LASOPS .NE. 0 ) THEN
            IF ( LASOPS .EQ. 1 ) THEN
               LASOPS = -11
            ELSE
               IF ( LASOPS .EQ. 2 ) THEN
                  LASOPS = -12
               ELSE
                  LASOPS = 0
                  NBLGRC(NRERR) = 2
                  IF( LANGAG .EQ. 0 ) THEN
                     KERR(1) = 'TR1DTEY: MAUVAISE VALEUR de LASOPS'
                     KERR(2) = '         ARRET du TRACE POSTSCRIPT'
                  ELSE
                     KERR(1) = 'TR1DTEY: BAD VALUE of LASOPS'
                     KERR(2) = '         STOP of POSTSCRIPT DRAWING'
                  ENDIF
                  CALL LEREUR
                  GOTO 10
               ENDIF
            ENDIF
            CALL XVPOSTSCRIPT(LASOPS)
            LASOPS = - LASOPS
            CALL XVPOSTSCRIPT(LASOPS)
         ENDIF
C
C        LE TRACE DU TITRE FINAL
C        =======================
         KNOM = 'OBJET: ' // KNOMOB
         I    = NUDCNB( KNOM )
         CALL XVCOULEUR( NCNOIR )
         CALL XVTEXTE( KNOM(1:I), I, 50, 30 )
C
C        LE TRACE DE LA LEGENDE : COULEURS => VALEURS
         NBCOUL = NDCOUL - N1COUL
         NCPAS  = NBCOUL / 10
         TPAS   = (TMAX-TMIN) / 10
         T      = TMIN
C
C        TRACE DE 11 VALEURS
         NCOUL = N1COUL
         NX    = LAPXFE - 160
         NY    = LHPXFE - 50
         DO 600 I=0,10
            CALL XVCOULEUR( NCOUL )
            CALL XVRECTANGLE( NX, NY, 30, 10 )
            WRITE( KNOM(1:10), '(G10.3)' ) T
            CALL XVTEXTE( KNOM(1:10), 10, NX+40, NY+10 )
            NCOUL = NCOUL + NCPAS
            T     = T  + TPAS
            NY    = NY - 15
 600     CONTINUE
C
C        RETOUR AU TRACE NORMAL POUR POSTSCRIPT
         IF ( LASOPS.NE.0 ) THEN
            LASOPS = LASOPS - 10
            CALL XVPOSTSCRIPT(LASOPS)
         ENDIF
C
C        RETOUR AUX PARAMETRES INITIAUX
         CALL XVEPAISSEUR( 1 )
         CALL XVTYPETRAIT( LIGCON )
C
C        DEFINITION DU TITRE ET FIN DU TRACE
         CALL LETITR( NOPROJ, MODECO, NCAS, TEMPS, KNOM )
C
C        ATTENDRE POUR LIRE LE TRACE
         CALL ATTENDSEC( TEMP2TRAC )
C
 1000 CONTINUE
C
C     RETOUR POUR UNE NOUVELLE VISEE
      IF( LORBITE .NE. 0 ) THEN
         CALL ZOOM1D1( NOTYEV )
         IF( NOTYEV .EQ. 0 ) GOTO 10
         GOTO 20
      ENDIF
C
      CALL CLICSO
      GOTO 10
C
C     SORTIE DU TRACE de la TEMPERATURE
      NDIM = 1
      COOEXT(1,1) = XMIN
      COOEXT(1,2) = XMAX
      COOEXT(2,1) = 0
      COOEXT(2,2) = 0
      COOEXT(3,1) = 0
      COOEXT(3,2) = 0
      INIEXT = 1
      MOAXYZ = 0
      XYZAMPLI(1) = 1
      XYZAMPLI(2) = 1
      XYZAMPLI(3) = 1
      RETURN
      END
