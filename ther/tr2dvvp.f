      SUBROUTINE TR2DVVP( KNOMOB, NDIM,   ValPr,  NBDLFX, VADLFX,
     %                    MNDLIB, NBDLIB, VECTLIB,
     %                    NBNOEU, VECTPR, NBTYEL, MNNPEF,NDPGST, MNXYZN,
     %                    MXSOMM, MNSOLE, MNCOPO, MXPIL3, MNPIL3,
     %                    MNVALS, MNFBAS, MNXYZS )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  TRACER UN VECTEUR PROPRE EN 2D
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC Paris & St Pierre du Perray Mai 2014
C23456---------------------------------------------------------------012
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE      (MCN(1),DMCN(1))
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/langue.inc"
      include"./incl/ponoel.inc"
      include"./incl/a___npef.inc"
      include"./incl/a___xyznoeud.inc"

      CHARACTER*(*)     KNOMOB
      DOUBLE PRECISION  VADLFX(NBDLFX), VECTLIB(NBDLIB), VECTPR(NBNOEU),
     %                  ValPr, VMIN, VMAX
      CHARACTER*160     KNOM
      CHARACTER*16      TEXT

      IF( NDIM .GE. 3 ) RETURN
C     NDIM=3 A TRAITER

C     COPIER LE VECTEUR PROPRE VECTLIB(NBDLIB) CALCULE
C     DANS LE TABLEAU VECTPR(NTDL) EN AJOUTANT LES NBDLFX DL
C     FIXES DE VALEURS DANS VADLFX(NBDLFX)
      IF( MNDLIB .LE. 0 ) THEN
         MNDLB = 1
      ELSE
         MNDLB = MNDLIB
      ENDIF
      CALL CPSIVV( 1, NBNOEU, NBDLIB, MCN(MNDLB), NBDLFX,  VADLFX,
     %             VECTLIB, VECTPR )

C     MIN ET MAX DU VECTPR
      VMIN = VECTPR(1)
      VMAX = VECTPR(1)
      DO I=2,NBNOEU
         IF( VECTPR(I) .LT. VMIN ) VMIN=VECTPR(I)
         IF( VECTPR(I) .GT. VMAX ) VMAX=VECTPR(I)
      ENDDO

C     L'ECRAN PIXELS EST EFFACE (PAS DE SCINTILLEMENT)
      CALL EFFACEMEMPX

C     TRACE DES AXES 2D
      CALL TRAXE2
C
C     BOUCLE SUR LES DIFFERENTS TYPES D'ELEMENTS FINIS DU MAILLAGE
      DO I = 0, NBTYEL-1
C
C        L'ADRESSE MCN DU DEBUT DU TABLEAU NPEF DE CE TYPE D'EF
         MNELE = MCN( MNNPEF + I )
C
C        LE NUMERO DU TYPE DES ELEMENTS FINIS
         NUTYEL = MCN( MNELE + WUTYEL )
C
C        LE NOMBRE DE TELS ELEMENTS
         NBELEM = MCN( MNELE + WBELEM )
C
C        LES CARACTERISTIQUES DE L'ELEMENT FINI
         CALL ELTYCA( NUTYEL )
C
C        L'ADRESSE MCN DU TABLEAU 'NPEF' POUR CE TYPE D'EF
         MNPGEL = MNELE + WUNDEL
         IF( NDPGST .GE. 2 ) THEN
            MNPGEL = MNPGEL + MCN(MNELE+WBELEM)*MCN(MNELE+WBNDEL)
         ENDIF

C        TRACE 2D EFFECTIF DES ZONES ISO-VALEUR DU VECTEUR PROPRE
         NBNOEU = MCN(MNXYZN+WNBNOE)
         NBCOOR = MCN(MNXYZN+WBCOON)
         CALL TRZON2( REAL(VMIN), REAL(VMAX), NDIM,
     %                1,  1,  1, NBNOEU, VECTPR,
     %                NUTYEL, NBELEM,
     %                NBNOE,  MCN(MNELE+WUNDEL),
     %                NBPOE,  MCN(MNPGEL),
     %                NBCOOR, NBNOEU, MCN(MNXYZN+WYZNOE),
     %                NBNOEU, MCN(MNXYZN+WYZNOE),
     %                MXSOMM, MCN(MNSOLE), MCN(MNCOPO),
     %                MXPIL3, MCN(MNPIL3),
     %                MCN(MNVALS), MCN(MNFBAS), MCN(MNXYZS) )
      ENDDO
C
C     LE TRACE DE LA LEGENDE: COULEURS => VALEURS
      CALL LEGCOULSO( REAL(VMIN), REAL(VMAX) )

C     TRACE DE LA 2-EME LIGNE DU TITRE DU TRACE
      MODECO = 2
      CALL TIT2LG( KNOMOB, MODECO )
C
C     DEFINITION DE LA 3-EME LIGNE DU TITRE
      IF( LANGAG .EQ. 0 ) THEN
         KNOM = 'VECTEUR PROPRE associe a la VALEUR PROPRE= '
      ELSE
         KNOM = 'EigenVector associated to EigenValue= '
      ENDIF

C     LA VALEUR PROPRE
      WRITE( TEXT, '(G14.6)' ) ValPr
C     SUPPRESSION DES BLANCS DE DEBUT ET INTERMEDIAIRES
      CALL TEXTSB( TEXT, L )
      I = NUDCNB( KNOM )
      KNOM(I+1:I+L) = TEXT(1:L)

C     LE MIN DU VECTEUR EN UN NOEUD
      I = NUDCNB( KNOM )
      KNOM(I+1:I+11) = ' Min VecPr='
      WRITE( TEXT, '(G14.6)' ) VMIN
      CALL TEXTSB( TEXT, L )
      I = NUDCNB( KNOM )
      KNOM(I+1:I+L) = TEXT(1:L)

C     LE MAX DU VECTEUR EN UN NOEUD
      I = NUDCNB( KNOM )
      KNOM(I+1:I+11) = ' Max VecPr='
      WRITE( TEXT, '(G14.6)' ) VMAX
      CALL TEXTSB( TEXT, L )
      I = NUDCNB( KNOM )
      KNOM(I+1:I+L) = TEXT(1:L)

C     RETOUR AUX PARAMETRES INITIAUX
      CALL XVEPAISSEUR( 1 )
      CALL XVTYPETRAIT( LIGCON )
C
C     TRACE DU TITRE
      CALL TRFINS( KNOM )
      print *,'tr2dvvp.f: KNOM=', KNOM

      RETURN
      END
