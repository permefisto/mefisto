      SUBROUTINE TRCO2D( KNOMOB, NBNOEU, XYZNOE,
     &                   MNNPEF, NBVECT, VITESX, VITESY, VITMAX,
     &                   TEMPSV, NCAS0,  NCAS1 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    TRACER LA FONCTION COURANT DES DIFFERENTS CAS DE VITESSE
C -----
C ENTREES:
C --------
C KNOMOB : NOM DE L'OBJET FLUIDE
C NBNOEU : NOMBRE DE NOEUDS SUPPORT DE LA VITESSE
C XYZNOE : 3 COORDONNEES DES NOEUDS VITESSE DU FLUIDE
C MNNPEF : ADRESSE DE L'ADRESSE DU TABLEAU NPEF"TYPE EF (TRIA 2P1D ou 2P2C)
C NBVECT : NOMBRE DE CAS TRAITES
C VITESX : VITESSE EN X AUX NOEUDS
C VITESY : VITESSE EN Y AUX NOEUDS
C VITMAX : NORME MAXIMALE DE LA VITESSE
C TEMPSV : NBVECT TEMPS DES NBVECT VECTEURS VITESSE
C NCAS0  : NUMERO DU PREMIER CAS A TRACER
C NCAS1  : NUMERO DU DERNIER CAS A TRACER
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE DU PERRAY  Janvier 2011
C23456---------------------------------------------------------------012
      PARAMETER        (LIGCON=0, LIGTIR=1 )
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      include"./incl/langue.inc"
      include"./incl/gsmenu.inc"
      include"./incl/trvari.inc"
      include"./incl/mecoit.inc"
      include"./incl/ponoel.inc"
      include"./incl/a___npef.inc"
      include"./incl/ctemps.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
C
      INTRINSIC         SQRT
      REAL              XYZNOE(3,NBNOEU), TEMPSV(NBVECT)
      DOUBLE PRECISION  VITESX(NBNOEU,NBVECT), VITESY(NBNOEU,NBVECT)
      DOUBLE PRECISION  VITNORM, VITMAX
      INTEGER           NONOEF(6)
      CHARACTER*(*)     KNOMOB
      CHARACTER*80      KNOM
      INTRINSIC         REAL, INT
C
C     LA PALETTE 2 : LES COULEURS VARIENT RAPIDEMENT ET S'ASSOMBRISSENT
C     LA PALETTE 11: ARC EN CIEL
      CALL PALCDE(11)
C
C     LE NOMBRE-1 DE COULEURS DISPONIBLES
      NBCOUL = NDCOUL - N1COUL
C
C     OPTIONS DE LA VISEE POUR VOIR LES PARTICULES
C     ============================================
 10   CALL VISE2D( NMTCL )
      IF( NMTCL .LT. 0 ) GOTO 9000
C
C     INITIALISATION DE TRANSLATION ZOOM
      IF( LORBITE .NE. 0 ) THEN
         CALL ZOOM2D0( NOTYEV )
         IF( NOTYEV .EQ. 0 ) GOTO 10
      ENDIF
C
C     ======================================================================
C     BOUCLE SUR LES CAS A TRACER
C     ======================================================================
 20   DO 100 NCAS = NCAS0, NCAS1
C
      CALL EFFACEMEMPX
C
C     TRACE DES AXES 2D
C     -----------------
      CALL TRAXE2
C
C     -----------------------
C     TRACE DES ARETES DES EF
C     -----------------------
      CALL XVEPAISSEUR( 1 )
      CALL XVTYPETRAIT( NTLAFR )
C
C     MNELE : ADRESSE DU TABLEAU NPEF"TYPE EF (TRIA 2P1D ou 2P2C)
      MNELE  = MCN( MNNPEF )
C
C     NOMBRE D'ELEMENTS FINIS DE CE TYPE
      NBELFI = MCN( MNELE + WBELEM )
C
C     LE NUMERO DU TYPE DE L'ELEMENT FINI (TRIA 2P1D ou 2P2C)
      NUTYEL = MCN( MNELE + WUTYEL )
C
C     NUTYEF NO DU TYPE DU TRIANGLE ELEMENT FINI FLUIDE
C     NUTYEF = 1 SI TYPE BREZZI-FORTIN avec TRIA 2P1D
C            = 2 SI TYPE TAYLOR-HOOD   avec TRIA 2P2C
C
      DO NUELEM=1,NBELFI
C        LE NUMERO DES NOEUDS DE L'EF NUELEM
         CALL EFNOEU( MNELE, NUELEM, NBNOEF, NONOEF )
         K = 3
         DO I=1,3
            NS1 = NONOEF(K)
            NS2 = NONOEF(I)
            CALL TRAIT2D( NCOUAF, XYZNOE(1,NS1) , XYZNOE(2,NS1),
     %                            XYZNOE(1,NS2) , XYZNOE(2,NS2) )
            K = I
         ENDDO
      ENDDO
C
C     ------------------
C     TRACE DES VITESSES
C     ------------------
      CALL XVEPAISSEUR( NEPFLE )
      CALL XVTYPETRAIT( LIGCON )
C
      DO 80 NS1=1,NBNOEU
C
C        LES 2 COORDONNEES DU NOEUD
         X = XYZNOE(1,NS1)
         Y = XYZNOE(2,NS1)
C
C        NORME DE LA VITESSE
         VITNORM = SQRT( VITESX(NS1,NCAS)**2
     %                 + VITESY(NS1,NCAS)**2 )
C
C        LONGUEUR CM DE LA FLECHE DE LA VITESSE
         COXF  = REAL( VITESX( NS1, NCAS ) * CMPCON )
         COYF  = REAL( VITESY( NS1, NCAS ) * CMPCON )
         CONCM = REAL( VITNORM * CMPCON )
C
C        COULEUR DU BOIS DE LA FLECHE
         NOCOUL = INT( VITNORM / VITMAX * NBCOUL + N1COUL )
C
C        LE TRACE DE LA FLECHE DE LA VITESSE
         CALL T2FLEC( NOCOUL, X, Y, CONCM, COXF, COYF )
C
 80   CONTINUE
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
               KERR(1) = 'TRCO2D: MAUVAISE VALEUR de LASOPS'
               KERR(2) = 'ARRET du TRACE POSTSCRIPT'
            ELSE
               KERR(1) = 'TRCO2D: BAD VALUE of LASOPS'
               KERR(2) = 'STOP of the POSTSCRIPT DRAWING'
            ENDIF
            CALL LEREUR
          ENDIF
        ENDIF
        CALL XVPOSTSCRIPT(LASOPS)
        LASOPS = - LASOPS
        CALL XVPOSTSCRIPT(LASOPS)
      ENDIF
C
C     LE TRACE DU TITRE FINAL
C     =======================
C     LE TYPE DES CARACTERES
      IF( LANGAG .EQ. 0 ) THEN
         KNOM = 'OBJET: ' // KNOMOB
      ELSE
         KNOM = 'OBJECT: ' // KNOMOB
      ENDIF
      I = NUDCNB( KNOM )
      CALL XVCOULEUR( NCNOIR )
      CALL XVTEXTE( KNOM(1:I), I, 50, 30 )
C
C     RETOUR AUX PARAMETRES INITIAUX
      CALL XVEPAISSEUR( 1 )
      CALL XVTYPETRAIT( LIGCON )
C
C     CONSTRUCTION DU TITRE ET TRACE
      IF( LANGAG .EQ. 0 ) THEN
         KNOM = 'VITESSE       au TEMPS                              CM
     %par UNITE de VITESSE'
      ELSE
         KNOM = 'VELOCITY      at TIME                               CM
     %per VELOCITY UNIT'
      ENDIF
      WRITE( KNOM(10:13),   '(I4)'    ) NCAS
      WRITE( KNOM(22:36), '(G15.7)' ) TEMPSV(NCAS)
C     CMPCON : NOMBRE DE CM PAR UNITE DE VITESSE
      WRITE( KNOM(38:52), '(G15.7)' ) CMPCON
      CALL TRFINS( KNOM )
C
C     ATTENDRE POUR LIRE LE TRACE
      CALL ATTENDSEC( TEMP2TRAC )
C
C     FIN DE LA BOUCLE SUR LES CAS
 100  CONTINUE
C
C     RETOUR POUR UNE NOUVELLE VISEE
C     ------------------------------
      IF( LORBITE .NE. 0 ) THEN
         IF( NCAS0 .EQ. NCAS1 ) THEN
C           ZOOM  BOUTON ENFONCE et DEPLACE
            CALL ZOOM2D1( NOTYEV )
         ELSE
C           ZOOM  BOUTON ENFONCE et DEPLACE et RELACHE
            CALL ZOOM2D3( NOTYEV )
         ENDIF
         IF( NOTYEV .EQ. 0 ) GOTO 10
         GOTO 20
      ELSE
         CALL CLICSO
      ENDIF
      GOTO 10
C
C     RETOUR AU TRACE NORMAL POUR POSTSCRIPT
 9000 IF( LASOPS.NE.0 ) THEN
        LASOPS = LASOPS -10
        CALL XVPOSTSCRIPT(LASOPS)
      ENDIF
C
      RETURN
      END
