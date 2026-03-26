      SUBROUTINE TRTABLES( NBPT, XTABLES, YTABLES,
     %                     NMFVIDEO, NBMETH, METHOD, COULMETH,
     %                     TXAxeX, TXAxeY, TITRE1, TITRE2, TITRE3 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     TRACER LES SYMBOLES (XTABLES(n),YTABLES(n)) n=1,...,NBPT
C -----     DEFINIS PAR LE TABLEAU METHOD des METHODES

C ENTREES:
C --------
C NBPT    : NOMBRE DE POINTS A TRACER
C XTABLES : ABSCISSE DES POINTS A TRACER
C YTABLES : ORDONNEE DES POINTS A TRACER
C NMFVIDEO: NOM DU FICHIER VIDEO a ECRIRE
C NBMETH  : NOMBRE DE METHODES UTILISEES
C METHOD  : Caractere A TRACER suivi du NOM DE LA METHODE cf LEGENDE
C COULMETH: NUMERO DE LA COULEUR ET DE LA METHODE POUR CHAQUE POINT
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
      CHARACTER*48   METHOD(NBMETH)
      INTEGER        COULMETH(NBPT)
      REAL           XTABLES(NBPT), YTABLES(NBPT)
      CHARACTER*(*)  TXAxeX, TXAxeY, NMFVIDEO,
     %               TITRE1, TITRE2, TITRE3

C     LE CADRE DE LA FENETRE 2D DU TRACE DES POINTS=SYMBOLES
      XMIN = XTABLES(1)
      XMAX = XTABLES(1)
      YMIN = YTABLES(1)
      YMAX = YTABLES(1)
      DO M = 1, NBPT
         XM = XTABLES(M)
         IF( XM .LT. XMIN ) XMIN = XM
         IF( XM .GT. XMAX ) XMAX = XM
         YM = YTABLES(M)
         IF( YM .LT. YMIN ) YMIN = YM
         IF( YM .GT. YMAX ) YMAX = YM
      ENDDO

C     TRAITEMENT DU CAS MIN=MAX
      IF( XMIN .GE. XMAX ) THEN
         XEC = XMAX / 10
         IF( XEC .EQ. 0 ) XEC = 1
         XMIN = XMIN - XEC
         XMAX = XMAX + XEC
      ENDIF

      IF( YMIN .GE. YMAX ) THEN
         YEC = YMAX / 10
         IF( YEC .EQ. 0 ) YEC = 1
         YMIN = YMIN - YEC
         YMAX = YMAX + YEC
      ENDIF

C     ECARTEMENT
      XEC = ( XMAX - XMIN ) / 10
      YEC = ( YMAX - YMIN ) / 10

ccc      print *,'TRTABLES: XMIN=',XMIN,'  XMAX=',XMAX,
ccc     %                '  YMIN=',YMIN,'  YMAX=',YMAX

cccC     MISE SUR FICHIER.eps du TRACE
ccc      CALL xvinitierps( 1 )

      CALL EFFACEMEMPX

C     DEFINITION DE LA FENETRE DE VISION
      CALL FENETRE( XMIN-XEC, XMAX+XEC, YMIN-YEC, YMAX+YEC )

C     TRACE DE LA 1-ERE LIGNE EN HAUT A GAUCHE DE LA FENETRE DE TRACE
C     AVEC 'MEFISTO', NomUtilisateur, DATE, ...
      CALL TIT1LG

C     LE TRACE DES AXES 2D
      NETAXE = 0
      CALL TRAXE2

C     LA SIGNIFICATION DES AXES
      CALL CHOIXFONTE( 20 )
      L = NUDCNB(TXAxeX)
      CALL TEXTE2D( NCROUG, XMAX-XEC, YMIN-YEC/2, TXAxeX(1:L) )
C
      L = NUDCNB(TXAxeY)
      CALL TEXTE2D( NCVERT, XMIN-XEC, YMAX-YEC*1.3, TXAxeY(1:L) )

C     LE TRACE DES YTABLES(m) EN FONCTION des XTABLES(m)
c     1-ER POINT
      XM0 = XTABLES(1)
      YM0 = YTABLES(1)

      DO M = 1, NBPT

C        LE POINT A TRACER ( XTABLES(M), YTABLES(M) )
         XM = XTABLES(M)
         YM = YTABLES(M)
ccc         CALL TRAIT2D( NCBLAN, XM0, YM0, XM, YM )
         NOCOU = COULMETH(M)
C        POUR EVITER UN TRACE AVEC LA COULEUR DU FOND
         IF( NOCOU .GE. NCOFON ) NOCOU=NOCOU+1
         CALL SYMBOLE2D( NOCOU, XM,  YM,
     %                   METHOD( COULMETH(M) )(1:1) )

         XM0 = XM
         YM0 = YM

      ENDDO

C     LE TITRE DU GRAPHIQUE
ccc      CALL CHOIXFONTE( 25 )
      L = NUDCNB( TITRE1 )
      IF( L .GT. 1 ) THEN
         CALL TEXTE2D( NCBLAN, XMIN+XEC, YMAX+YEC/2, TITRE1(1:L) )
      ENDIF
C
      L = NUDCNB( TITRE2 )
      IF( L .GT. 1 ) THEN
         CALL TEXTE2D( NCBLAN, XMIN+XEC, YMAX, TITRE2(1:L) )
      ENDIF

      L = NUDCNB( TITRE3 )
      IF( L .GT. 1 ) THEN
         CALL TEXTE2D( NCBLAN, XMIN+XEC, YMAX-YEC/2, TITRE3(1:L) )
      ENDIF

C     TRACE SYMBOLE => METHODE
      XEC = ( XMAX - XMIN ) / 20
      YEC = ( YMAX - YMIN ) / 20
      XM = (XMIN+XMAX)/2 + XEC * 2
      YM = YMIN + (YMAX-YMIN)/3 + YEC * (NBMETH+2)
      DO M = 1, NBMETH
         L = NUDCNB( METHOD(M) )
         NOCOU = M
C        POUR EVITER UN TRACE AVEC LA COULEUR DU FOND
         IF( NOCOU .GE. NCOFON ) NOCOU=NOCOU+2
         CALL SYMBOLE2D( NOCOU, XM, YM, METHOD(M)(1:L) )
         YM = YM - YEC
      ENDDO

C     COPIE DE MEMPX DANS FENETRE
      CALL MEMPXFENETRE

C     POUR VIDER LE BUFFER DE X11
      CALL XVVOIR

      IF( LVIDEO .NE. 0 ) THEN

C     A PARTIR DE LA PIXMAP de la FENETRE X11 ACTUELLE CREATION DU FICHIER.xwd
      L  = NUDCNB( NMFVIDEO )
      LL = NUDCNB( NMPROJ )
      print *
      print *, 'xwd -xy -name Mefisto -out '//NMFVIDEO(1:L)//'.xwd',
     %        ' EST EXECUTE'
      CALL SYSTEM('xwd -xy -name Mefisto -out '//NMFVIDEO(1:L)//'.xwd')

C     CONVERSION DU FICHIER.xwd en le FICHIER.jpg
      print *,    'convert ' // NMFVIDEO(1:L) // '.xwd' // ' ' //
     %  NMPROJ(1:LL)  // '_' // NMFVIDEO(1:L) // '.jpg  EST EXECUTE'
      CALL SYSTEM('convert ' // NMFVIDEO(1:L) // '.xwd' // ' ' //
     %  NMPROJ(1:LL)  // '_' // NMFVIDEO(1:L) // '.jpg' )

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

      RETURN
      END
