      SUBROUTINE TRTABLE( NBXYTAB,   XTABLE,
     %                    NC1YTABLE, N1YTABLE, YTABLE,
     %                    NMFVIDEO,  NCGRAF,
     %                    TXAxeX,    TXAxeY,
     %                    TITRE1,    TITRE2, TITRE3 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : TRACER LA COURBE (XTABLE(n),YTABLE(NC1YTABLE,n)) n=1,...,NBXYTAB
C -----
C ENTREES:
C --------
C NBXYTAB : NOMBRE DE POINTS A TRACER
C XTABLE  : ABSCISSE DES POINTS A TRACER
C NC1YTABLE: NUMERO DU 1-ER INDICE de YTABLE A TRACER A CHAQUE TEMPS
C N1YTABLE: NOMBRE DE VARIABLES DU 1-ER INDICE DU TABLEAU YTABLE
C YTABLE  : ORDONNEE DES POINTS A TRACER TABLEAU (N1YTABLE,NBXYTAB)
C NMFVIDEO: NOM DU FICHIER VIDEO a ECRIRE
C NCGRAF  : NO DE LA COULEUR DU GRAPHE
C TXAxeX  : TEXTE de l'AXE X
C TXAxeY  : TEXTE de l'AXE Y
C TITRE1  : 1-ERE LIGNE DU TITRE
C TITRE2  : 2-EME LIGNE DU TITRE
C TITRE3  : 3-EME LIGNE DU TITRE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  TEXAS A & M UNIVERSITY at QATAR Fevrier 2012
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/xvfontes.inc"
      include"./incl/trvari.inc"
      include"./incl/traaxe.inc"
      include"./incl/mecoit.inc"
      include"./incl/nmproj.inc"
      REAL             XTABLE(NBXYTAB)
      DOUBLE PRECISION YTABLE(N1YTABLE,NBXYTAB), DYMIN, DYMAX, DYEC, DYM
      CHARACTER*(*)    TXAxeX, TXAxeY, NMFVIDEO,
     %                 TITRE1, TITRE2, TITRE3
      INTRINSIC        REAL
C
ccc      print *,'TRTABLE ',nmfvideo
ccc      do m=1,nbxytab
ccc         print *,'xtable(',m,')=',xtable(m),
ccc     %           ' ytable=',(ytable(n,m),n=1,n1ytable)
ccc      enddo
C
C     LE CADRE DE LA FENETRE 2D DE TRACE
      XMIN = XTABLE(1)
      XMAX = XTABLE(1)
      DYMIN = YTABLE(NC1YTABLE,1)
      DYMAX = YTABLE(NC1YTABLE,1)
      DO M = 2, NBXYTAB
         XM = XTABLE(M)
         IF( XM .LT. XMIN ) XMIN = XM
         IF( XM .GT. XMAX ) XMAX = XM
         DYM = YTABLE(NC1YTABLE,M)
         IF( DYM .LT. DYMIN ) DYMIN = DYM
         IF( DYM .GT. DYMAX ) DYMAX = DYM
      ENDDO

      PRINT *,'TRTABLE: XMIN=',XMIN,'  XMAX=',XMAX,
     %                ' DYMIN=',DYMIN,' DYMAX=',DYMAX
C
C     TRAITEMENT DU CAS MIN=MAX ET TRES PETITES VALEURS
      IF( ABS(XMAX-XMIN) .LT. 1E-20 ) THEN
         XEC = XMAX / 10
         IF( ABS(XEC) .LT. 1E-20 ) XEC = 1.
         XMIN = XMIN - XEC
         XMAX = XMAX + XEC
      ENDIF
C
      IF( ABS(DYMAX-DYMIN) .LT. 1E-20 ) THEN
         DYEC = DYMAX / 10
         IF( ABS(DYEC) .LT. 1D-20 ) DYEC = 1D0
         DYMIN = DYMIN - DYEC
         DYMAX = DYMAX + DYEC
      ENDIF

      YMIN = REAL( DYMIN )
      YMAX = REAL( DYMAX )
C
C     ECARTEMENT
      XEC = ABS( XMAX - XMIN ) / 10.
      YEC = ABS( YMAX - YMIN ) / 10.

      print *,'TRTABLE: XMIN=',XMIN,'  XMAX=',XMAX,
     %               '  YMIN=',YMIN,'  YMAX=',YMAX,
     %         ' XEC=',XEC, ' YEC=',YEC
C
cccC     MISE SUR FICHIER.eps du TRACE
ccc      CALL xvinitierps( 1 )
C
      CALL EFFACEMEMPX
C
C     LE FOND EST BLANC
      NCOF0 = NCOFON
      CALL XVFOND( NCBLAN )
C
C     DEFINITION DE LA FENETRE DE VISION
      CALL FENETRE( XMIN-XEC, XMAX+XEC, YMIN-YEC, YMAX+YEC )
C
C     TRACE DE LA 1-ERE LIGNE EN HAUT A GAUCHE DE LA FENETRE DE TRACE
C     AVEC 'MEFISTO', NomUtilisateur, DATE, ...
      CALL TIT1LG
C
C     LE TRACE DES AXES 2D
      NETAXE = 0
      CALL TRAXE2
C
C     LA SIGNIFICATION DES AXES
      CALL CHOIXFONTE( 20 )
      L = NUDCNB(TXAxeX)
      CALL TEXTE2D( NCROUG, XMAX-XEC, YMIN-YEC*0.5, TXAxeX(1:L) )
C
      L = NUDCNB(TXAxeY)
      CALL TEXTE2D( NCVERT, XMIN-XEC, YMAX-YEC*0.45, TXAxeY(1:L) )
C
C     LE TRACE DES YTABLE(NC1YTABLE,m) EN FONCTION des XTABLE(m)
      XM0 = XTABLE(1)
      DYM = YTABLE(NC1YTABLE,1)
      IF( ABS(DYM) .LT. 1D-20 ) DYM = 0D0
      YM0 = REAL( DYM )
      DO M = 1, NBXYTAB
C
C        LE POINT A TRACER ( XTABLE(M), YTABLE(NC1YTABLE,M) )
          XM = XTABLE(M)
         DYM = YTABLE(NC1YTABLE,M)
         IF( ABS(DYM) .LT. 1D-27 ) DYM = 0D0
         YM = REAL( DYM )
CCCC        LE NO
CCC         CALL ENTIER2D(  NCGRAF, XM, YM, ITABLE(M) )
         CALL TRAIT2D( NCGRAF, XM0, YM0, XM, YM )
         CALL SYMBOLE2D( NCGRAF, XM, YM, 'x' )
C
         XM0 = XM
         YM0 = YM
C
      ENDDO

C     LE TITRE DU GRAPHIQUE
ccc      CALL CHOIXFONTE( 25 )
      L = NUDCNB( TITRE1 )
      IF( L .GT. 1 ) THEN
         CALL TEXTE2D( NCGRAF, XMIN+XEC, YMAX, TITRE1(1:L) )
      ENDIF

      L = NUDCNB( TITRE2 )
      IF( L .GT. 1 ) THEN
         CALL TEXTE2D( NCGRAF, XMIN+XEC, YMAX-YEC/2, TITRE2(1:L) )
      ENDIF

      L = NUDCNB( TITRE3 )
      IF( L .GT. 1 ) THEN
         CALL TEXTE2D( NCGRAF, XMIN+XEC, YMAX-YEC, TITRE3(1:L) )
      ENDIF

C     COPIE DE MEMPX DANS FENETRE
      CALL MEMPXFENETRE

C     POUR VIDER LE BUFFER DE X11
      CALL XVVOIR

      IF( LVIDEO .NE. 0 ) THEN

C     A PARTIR DE LA PIXMAP de la FENETRE X11 ACTUELLE CREATION DU FICHIER.xwd
      L = NUDCNB( NMFVIDEO )
      print *
      print *, 'xwd -xy -name Mefisto -out '//NMFVIDEO(1:L)//'.xwd',
     %        ' EST EXECUTE'
      CALL SYSTEM('xwd -xy -name Mefisto -out '//NMFVIDEO(1:L)//'.xwd')

C     CONVERSION DU FICHIER.xwd en le FICHIER.jpg
      LL = NUDCNB( NMPROJ )
      print *,    'convert ' // NMFVIDEO(1:L) // '.xwd' // ' '
     %                       // NMPROJ(1:LL)  // '_' // NMFVIDEO(1:L) //
     %                     '.jpg  EST EXECUTE'
      CALL SYSTEM('convert ' // NMFVIDEO(1:L) // '.xwd' // ' '
     %                       // NMPROJ(1:LL)  // '_' // NMFVIDEO(1:L) //
     %                     '.jpg')

C     DESTRUCTION DU FICHIER .xwd
      CALL SYSTEM( 'rm -Rf ' // NMFVIDEO(1:L) // '.xwd' )

cccC        MISE SUR FICHIER NMFVIDEO.eps du TRACE
cccC        ATTENTION PASSAGE PAR VARIABLE OBLIGATOIRE
ccc         NMAUX = NMFVIDEO
ccc         L = NUDCNB( NMAUX )
ccc         CALL xvsauverps( NMAUX(1:L), L )
cccC        CONVERSION DU FICHIER.eps EN FICHIER.jpg   NE MARCHE PAS
ccc         print *,     'convert ' //  NMAUX(1:L)//'.eps ' //
ccc     %                 NMAUX(1:L)//'.jpg     EST EXECUTE'
cccc        CALL SYSTEM( 'convert ' // NMAUX(1:L)//'.eps ' //
ccc     %                 NMAUX(1:L)//'.jpg ' )

      ENDIF

C     POUR ATTENDRE UN CLIC SOURIS ET  LIRE LE GRAPHIQUE
      CALL CHOIXFONTE( NPHFCO )
      CALL CLICSO

C     RETOUR AUX PARAMETRES INITIAUX
      CALL XVEPAISSEUR( 1 )
C     LE FOND REPREND SA COULEUR INITIALE
      CALL XVFOND( NCOF0 )

      RETURN
      END
