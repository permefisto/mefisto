      SUBROUTINE TRER3D( NDIM,   KNOMOB, MODECO0,
     %                   NBTYEL, MNTOPO, MNELEM, MNPOGE, MNNOEU, NDPGST,
     %                   NCAS0,  NCAS1,  NTDL,   ERREUR, dptemp,
     %                   ERRMIN, NOEMIN, NCAMIN, ERRMAX, NOEMAX, NCAMAX,
     %                   TIMES )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER L'ERREUR |ERREUR EXACTE(Noeud)-ERREUR CALCULEE(Noeud)|
C -----    A PARTIR DU TABLEAU ERREUR(NTDL,NDSM)
C
C ENTREES:
C --------
C NDIM   : DIMENSION DE L'ESPACE DE L'OBJET (2 OU 3 OU 6)
C KNOMOB : NOM DE L'OBJET
C MODECO0: MODE DE TRACE DES VECTEURS
C          =1 CE SONT DES TEMPERATURES
C          =2 CE SONT DES VECTEURS PROPRES
C          =3 CE SONT DES PRESSIONS P1 A PARTIR D'UNE INTERPOLATION
C             SOIT P2 SOIT P1+BULLE P3
C          =4 CE SONT DES ERREURS AUX NOEUDS D'UN MAILLAGE
C          =8 CE SONT DES ERREURS AUX NOEUDS D'UN MAILLAGE
C             DU MODULE D'UNE ONDE COMPLEXE
C            (CALCUL DU MODULE DE L'ERREUR COMPLEXE A FAIRE)
C NBTYEL : LE NOMBRE DE TYPES D'ELEMENTS FINIS DANS CETTE TOPOLOGIE
C MNTOPO : ADRESSE MCN DU TABLEAU TOPOLOGIE
C MNELEM : ADRESSE MCN DES TABLEAUX ELEMENTS FINIS DE CETTE TOPOLOGIE
C MNPOGE : ADRESSE MCN DU TABLEAU POINTS GEOMETRIQUES DE L'OBJET
C MNNOEU : ADRESSE MCN DU TABLEAU NOEUDS D'INTERPOLATION DU MAILLAGE
C NDPGST : CODE TRAITEMENT DES SOMMETS POINTS NOEUDS DU MAILLAGE
C          0 : NOEUDS=POINTS=SOMMETS
C          1 : NOEUDS=POINTS#SOMMETS
C          2 : NOEUDS#POINTS=SOMMETS
C          3 : NOEUDS#POINTS#SOMMETS
C NCAS0  : NUMERO DU PREMIER CAS A TRACER
C NCAS1  : NUMERO DU DERNIER CAS A TRACER
C NTDL   : NOMBRE TOTAL DE DEGRES DE LIBERTE (1 PAR NOEUD)
C ERREUR : ERREUR(NTDL,NCAS0:NCAS1) AUX NOEUDS
C ERRMIN : ERREUR MINIMALE EN UN NOEUD
C NOEMIN : NUMERO DU NOEUD OU L'ERREUR EST MINIMALE
C ERRMAX : ERREUR MAXIMALE EN UN NOEUD
C NOEMAX : NUMERO DU NOEUD OU L'ERREUR EST MAXIMALE
C TIMES  : TEMPS DU CALCUL DES NCAS0:NCAS1 VECTEURS ERREUR
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St Pierre du Perray  Fevrier 2009
C MODIFS : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY Novembre 2010
C MODIFS : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY  Juillet 2011
C MODIFS : ALAIN PERRONNET Saint Pierre du Perray              Mars 2021
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/gsmenu.inc"
      include"./incl/a___npef.inc"
      include"./incl/a___nsef.inc"
      include"./incl/a___xyzpoint.inc"
      include"./incl/a___xyznoeud.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___vecteur.inc"
      include"./incl/ponoel.inc"
      include"./incl/donele.inc"
      include"./incl/ctemps.inc"
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/xyzext.inc"
C
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      CHARACTER*(*)     KNOMOB
      DOUBLE PRECISION  ERREUR(NTDL,NCAS0:NCAS1)
      REAL              TIMES(NCAS0:NCAS1)
      type typ_dptab
         DOUBLE PRECISION, dimension(:), pointer :: dptab
      end type typ_dptab
      type( typ_dptab ),dimension(0:1) :: dptemp

      IF( INTERA .LE. 0 ) RETURN
      IERR   = 0
      NOPROJ = 0
C
C     LECTURE DE L'OPTION DE TRACE DE L'ERREUR
C     ----------------------------------------
 10   CALL LIMTCL( 'tracerr3', NMTCL0 )
      IF( NMTCL0 .LE. 0 ) GOTO 9999
C
C     POUR EVITER DES PUISSANCES NEGATIVES FAIBLES DE 10
      ERRMIN = 0
C
C     PARTICULARITE DE L'ERREUR DU MODULE DE L'ONDE COMPLEXE
      IF( MODECO0 .EQ. 8 ) THEN
         MODECO = 4
      ELSE
         MODECO = MODECO0
      ENDIF
C
      GOTO ( 210, 220, 230, 240, 250, 260 ), NMTCL0
C
C     TRACE DES SURFACES ISOERREURS
C     =============================
 210  CALL TRISOT( NDIM,    KNOMOB,  MODECO,
     %             NBTYEL,  MNELEM,  MNPOGE, NDPGST,
     %             NCAS0,   NCAS1,   NTDL,
     %             0,       ERREUR,  dptemp,
     %             ERRMIN,  NOEMIN,  NCAMIN, ERRMAX, NOEMAX, NCAMAX,
     %             TIMES )
      GOTO 10
C
C     TRACE DES ZONES DE COULEURS ISOERREURS
C     ======================================
 220  CALL TRZONT( NOPROJ,  NDIM,    KNOMOB, MODECO,
     %             NBTYEL,  MNELEM,  MNPOGE, MNNOEU, NDPGST,
     %             NCAS0,   NCAS1,   NTDL,
     %             0,       ERREUR,  dptemp,
     %             ERRMIN,  NOEMIN,  NCAMIN, ERRMAX, NOEMAX, NCAMAX,
     %             TIMES )
      GOTO 10
C
C     TRACE DES ZONES DE COULEURS ISOERREURS PAR SECTIONS X ou Y ou Z=CTE
C     ===================================================================
 230  CALL TRPLSE( 0,       KNOMOB,  MODECO,
     %             NBTYEL,  MNELEM,  MNPOGE, NDPGST,
     %             NCAS0,   NCAS1,   NTDL,
     %             0,       ERREUR,  dptemp,
     %             ERRMIN,  NOEMIN,  NCAMIN, ERRMAX, NOEMAX, NCAMAX,
     %             TIMES )
      GOTO 10
C
C     TRACE DES PROFILS D'ERREURS EN COULEURS PAR SECTIONS X ou Y ou Z=CTE
C     ====================================================================
 240  CALL TRPLSE( 1,       KNOMOB,  MODECO,
     %             NBTYEL,  MNELEM,  MNPOGE,  NDPGST,
     %             NCAS0,   NCAS1,   NTDL,
     %             0,       ERREUR,  dptemp,
     %             ERRMIN,  NOEMIN,  NCAMIN,  ERRMAX, NOEMAX, NCAMAX,
     %             TIMES )
      GOTO 10
C
C     TRACE DE L'ERREUR LE LONG D'UNE DROITE DEFINIE PAR 2 POINTS
C     ===========================================================
 250  CALL TRLLDR( NOPROJ,  NDIM,    KNOMOB,  MODECO,
     %             NBTYEL,  MNTOPO,  MNELEM,  MNPOGE, MNNOEU, NDPGST,
     %             NCAS0,   NCAS1,   NTDL,
     %             0,       ERREUR,  dptemp,
     %             ERRMIN,  NOEMIN,  NCAMIN,  ERRMAX, NOEMAX, NCAMAX,
     %             TIMES )
      GOTO 10
C
C     TRACE de la COURBE ERREUR(Temps)
C     ================================
 260  CALL LXLXOU( NTOBJE, KNOMOB, NTLXOB, MNLXOB )
      IF( NTLXOB .GT. 0 ) THEN
         CALL TRNLSERR( NTLXOB )
      ENDIF
      GOTO 10
C
 9999 RETURN
      END
