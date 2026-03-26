      SUBROUTINE TRVMTR21( MISTRE, KNOMOB, NBTYEL, MNELEM, NDPGST,
     %                     MODECO, NCAS,   CONTMN, CONTMX,
     %                     NBST,   MNCRIT,
     %                     NBPOIT, XYZPOI )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER EN 2D LE CRITERE DES CONTRAINTES DE VON MISES TRESCA
C -----    PAR ZONES DE COULEURS SUR LES EF 2D TRIANGLE OU QUADRANGLE
C
C ENTREES :
C ---------
C MISTRE : 1 POUR CRITERE DE VON MISES
C          2 POUR CRITERE DE TRESCA
C KNOMOB : NOM DE L'OBJET
C NBTYEL : LE NOMBRE DE TYPES D'ELEMENTS DANS CETTE TOPOLOGIE
C MNELEM : ADRESSE MCN DES TABLEAUX ELEMENTS DE CETTE TOPOLOGIE
C NDPGST : CODE TRAITEMENT DES SOMMETS POINTS NOEUDS DU MAILLAGE
C          0 : NOEUDS=POINTS=SOMMETS
C          1 : NOEUDS=POINTS#SOMMETS
C          2 : NOEUDS#POINTS=SOMMETS
C          3 : NOEUDS#POINTS#SOMMETS
C MODECO : MODE DES VECTEURS SOLUTIONS
C          =1 CE SONT DEPLACEMENTS
C          =2 CE SONT DES MODES PROPRES
C NCAS   : NUMERO DU CAS A TRAITER
C CONTMN : CONTRAINTE MIN
C CONTMX : CONTRAINTE MAX
C NBST   : NOMBRE DE SOMMETS DES TYPES D'EF DU CRITERE
C MNCRIT : ADRESSE MCN DES TABLEAUX CRITERE(NBPIEX,NBELFI)
C NBPOIT : NOMBRE TOTAL DE POINTS DU MAILLAGE
C XYZPOI : COORDONNEES DES POINTS DU MAILLAGE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : Alain PERRONNET Laboratoire J-L. LIONS UPMC PARIS    MAI 2007
C23456---------------------------------------------------------------012
      PARAMETER     (LIGCON=0)
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/a_objet__topologie.inc"
      include"./incl/a___npef.inc"
      include"./incl/a___xyzpoint.inc"
      include"./incl/a___xyznoeud.inc"
      include"./incl/a___face.inc"
      include"./incl/a___contrainte.inc"
      include"./incl/ponoel.inc"
      include"./incl/donele.inc"
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/ctemps.inc"
      include"./incl/xyzext.inc"
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
      CHARACTER*(*)     KNOMOB
      CHARACTER*160     KNOM
      REAL              XYZPOI(3,NBPOIT)
      REAL              CONTMN, CONTMX
      INTEGER           NBST(4), MNCRIT(4)
C
C     OPTIONS DE LA VISEE 2D DU TRACE
C     -------------------------------
 10   CALL VISE2D( NMTCL )
      IF( NMTCL .LT. 0 ) GOTO 9000
C
C     INITIALISATION DE TRANSLATION ZOOM
      CALL ZOOM2D0( NOTYEV )
      IF( NOTYEV .EQ. 0 ) GOTO 10
      GOTO 20
C
 15   IF( LORBITE .NE. 0 ) THEN
         CALL ZOOM2D1( NOTYEV )
         IF( NOTYEV .EQ. 0 ) GOTO 10
      ENDIF
C
C     TRACE DES AXES 2D
C     -----------------
 20   CALL TRAXE2
C
C     TRACE 2D DES TABLEAUX CRITERE(NBPIEX,NBELFI)
C     --------------------------------------------
      DO 30 I = 0, NBTYEL-1
C
C        NOMBRE DE SOMMETS POUR LE CRITERE
         NBPIEX = NBST(I+1)
C
C        L'ADRESSE MCN DU DEBUT DU TABLEAU NPEF
         MNELE = MCN( MNELEM + I )
C
C        LE NUMERO DU TYPE DES ELEMENTS FINIS
         NUTYEL = MCN( MNELE + WUTYEL )
C
C        LE NOMBRE DE TELS ELEMENTS
         NBELEM = MCN(MNELE + WBELEM )
C
C        LES CARACTERISTIQUES DE L'ELEMENT FINI
         CALL ELTYCA( NUTYEL )
C
C        L'ADRESSE MCN DU TABLEAU 'NPEF' POUR CE TYPE D'EF
         MNPGEL = MNELE + WUNDEL
         IF( NDPGST .GE. 2 ) THEN
            MNPGEL = MNPGEL + MCN(MNELE+WBELEM) * MCN(MNELE+WBNDEL)
         ENDIF
C
C        TRACE 2D EFFECTIF DU CRITERE PAR ZONES DE COULEURS
C        PASSAGE PAR UN SP AFIN DE BENEFICIER DES INDICES FORTRAN EN CLAIR
         CALL TRVMTR22( CONTMN, CONTMX,
     %                  NBPIEX, RMCN(MNCRIT(I+1)),
     %                  NBELEM, NBPOE,  MCN(MNPGEL),
     %                  NBPOIT, XYZPOI )
 30   CONTINUE
C
C     EFFACEMENT DE LA LEGENDE SUR POSTSCRIPT
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
               KERR(1) = 'TRVMTR21: MAUVAISE VALEUR DE LASOPS'
               KERR(2) = '          ARRET du TRACE POSTSCRIPT'
            ELSE
               KERR(1) = 'TRVMTR21: BAD VALUE of LASOPS'
               KERR(2) = '          STOP of POSTSCRIPT DRAWING'
            ENDIF
            CALL LEREUR
          ENDIF
        ENDIF
        CALL XVPOSTSCRIPT( LASOPS )
        LASOPS = - LASOPS
        CALL XVPOSTSCRIPT( LASOPS )
      ENDIF
C
C     LE TRACE DU TITRE FINAL
C     =======================
      IF( LANGAG .EQ. 0 ) THEN
         KNOM = 'OBJET: ' // KNOMOB
      ELSE
         KNOM = 'OBJECT: ' // KNOMOB
      ENDIF
      I    = NUDCNB( KNOM )
      CALL XVCOULEUR( NCGRIS )
      CALL XVTEXTE( KNOM(1:I), I, 50, 30 )
C
C     LE TRACE DE LA LEGENDE : COULEURS => VALEURS
      NBCOUL = NDCOUL - N1COUL
      NCPAS  = NBCOUL / 10
      TPAS   = (CONTMX-CONTMN) / 10
      T      = CONTMN
C
C     TRACE DE 11 VALEURS
      NCOUL = N1COUL
      NX    = LAPXFE - 170
      NY    = LHPXFE - 30
      DO 60 I=0,10
         CALL XVCOULEUR( NCOUL )
         CALL XVRECTANGLE( NX, NY, 30, 10 )
         WRITE( KNOM(1:10), '(G10.3)' ) T
         CALL XVTEXTE( KNOM(1:10), 10, NX+40, NY+10 )
         NCOUL = NCOUL + NCPAS
         T     = T  + TPAS
         NY    = NY - 15
 60   CONTINUE
C
C     RETOUR AU TRACE NORMAL POUR POSTSCRIPT
      IF ( LASOPS .NE. 0 ) THEN
        LASOPS = LASOPS - 10
        CALL XVPOSTSCRIPT( LASOPS )
      ENDIF
C
C     FIN DU TRACE
      IAVTIT = 1
      WRITE( KNOM(1:4), '(I4)' ) NCAS
      WRITE( KNOM(5:19), '(G15.6)' ) TEMPS
      IF( LANGAG .EQ. 0 ) THEN
         IF( MISTRE .EQ. 1 ) THEN
            IF( MODECO .EQ. 1 ) THEN
               CALL TRFINS( 'CRITERE de PLASTICITE de Von MISES: Cas'
     %                     //KNOM(1:4) // ' au temps ' // KNOM(5:19) )
            ELSE
               CALL TRFINS( 'CRITERE de PLASTICITE de Von MISES: '
     %  //'Frequence Propre '// KNOM(1:4) //': ' // KNOM(5:19) // ' Hz')
            ENDIF
         ELSE
            IF( MODECO .EQ. 1 ) THEN
               CALL TRFINS( 'CRITERE de PLASTICITE de TRESCA: Cas'
     %                      //KNOM(1:4) // ' au temps ' // KNOM(5:19) )
            ELSE
               CALL TRFINS( 'CRITERE de PLASTICITE de TRESCA: '//
     %    'Frequence Propre '// KNOM(1:4) //': ' // KNOM(5:19) // ' Hz')
            ENDIF
         ENDIF
      ELSE
         IF( MISTRE .EQ. 1 ) THEN
            IF( MODECO .EQ. 1 ) THEN
               CALL TRFINS( 'Von MISES''s PLASTICITY Criterion: Case'
     %                      //KNOM(1:4) // ' at time ' // KNOM(5:19) )
            ELSE
               CALL TRFINS( 'Von MISES''s PLASTICITY Criterion: '//
     %    'EigenFrequency '// KNOM(1:4) //': ' // KNOM(5:19) // ' Hz')
            ENDIF
         ELSE
            IF( MODECO .EQ. 1 ) THEN
               CALL TRFINS( 'TRESCA''s PLASTICITY Criterion: Case'
     %                     //KNOM(1:4) // ' at time ' // KNOM(5:19) )
            ELSE
               CALL TRFINS( 'TRESCA''s PLASTICITY Criterion: '//
     %    'EigenFrequency '// KNOM(1:4) //': ' // KNOM(5:19) // ' Hz')
            ENDIF
         ENDIF
      ENDIF
C
C     RETOUR POUR UNE NOUVELLE VISEE
      IF( LORBITE .NE. 0 ) GOTO 15
      CALL CLICSO
      GOTO 10
C
C     RETOUR AUX PARAMETRES INITIAUX
 9000 CALL XVEPAISSEUR( 1 )
      CALL XVTYPETRAIT( LIGCON )
      RETURN
      END
