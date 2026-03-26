      SUBROUTINE XVINIT
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : RECUPERER LES INFORMATIONS ET INITIALISER COULEURS ...
C ----- POUR FAIRE DES TRACES AVEC LA LIBRAIRIE XVUE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN ANALYSE NUMERIQUE UPMC  PARIS        MAI 1994
C2345X7..............................................................012
      include"./incl/langue.inc"
      include"./incl/trvari.inc"
      include"./incl/xvfontes.inc"
      include"./incl/xvpalette.inc"
      include"./incl/traaxe.inc"
      include"./incl/xyzext.inc"
      COMMON /UNITES/ LECTEU,IMPRIM,INTERA,NUNIT(29)
C
C     OUVERTURE XVUE
C     ==============
CCC      WRITE(IMPRIM,*)'Librairie Graphique XVUE'
CCCC     LA FENETRE EST DEJA OUVERTE PAR XTINIT => PAS D'APPEL DE XVINITGRAPHIQU
CCC      CALL XVINITGRAPHIQUE
C
C     INITIALISATION DES CARACTERISTIQUES DE L'ECRAN
C     ==============================================
C     LARGEUR ET HAUTEUR DE L'ECRAN EN MM
      CALL XVMMECRAN( LAMMEC, LHMMEC )
C     LARGEUR ET HAUTEUR EN PIXELS DE L'ECRAN
      CALL XVPXECRAN( LAPXEC, LHPXEC )
10000 FORMAT(1X,A,T22,' LARGEUR=',I4,'  HAUTEUR=',I4)
10001 FORMAT(1X,A,T20,'WIDTH=',I4,'  HEIGHT=',I4)
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10000) 'PIXELS de l''ECRAN',LAPXEC,LHPXEC
         WRITE(IMPRIM,10000) 'MM     de l''ECRAN',LAMMEC,LHMMEC
      ELSE
         WRITE(IMPRIM,10001) 'Screen PIXELS',LAPXEC,LHPXEC
         WRITE(IMPRIM,10001) 'Screen MM    ',LAMMEC,LHMMEC
      ENDIF
      CXMMPX = FLOAT(LAPXEC)/FLOAT(LAMMEC)
      CYMMPX = FLOAT(LHPXEC)/FLOAT(LHMMEC)
C
C     ------------------------------------------------------------------
C     -------------- PARTIE DEPENDANTE de l' ECRAN ---------------------
C     ------------------------------------------------------------------
C     DEFINITION DE LA FENETRE XVUE EN PIXELS ET OUVERTURE
C     ====================================================
C     TAILLE DE LA FENETRE GRAPHIQUE EN PIXELS A OUVRIR
      IF( LHPXEC .LE. 600 ) THEN
C        ECRAN 800x600
         LAPXFE = LAPXEC - 24
         LHPXFE = LHPXEC - 24
      ELSE IF( LHPXEC .LE. 800 ) THEN
C        ECRAN 1024x768
         LAPXFE = LAPXEC - 40
         LHPXFE = LHPXEC - 60
      ELSE IF( LHPXEC .LE. 900 ) THEN
C        VERSION MAC G4 1280x854 de DEFINITION D'ECRAN
         LAPXFE = LAPXEC - 96
         LHPXFE = LHPXEC - 48
      ELSE IF( LHPXEC .LE. 1080 ) THEN
C        ECRAN 1920 x 1080 ou ECRAN 1280 x 1024 OU ECRAN 1650 x 1050
         LAPXFE = LAPXEC - 108
         LHPXFE = LHPXEC - 100
      ELSE
C        ECRAN 1920 x 1200 ou ECRAN 1600 x 1200 OU PLUS
         LAPXFE = LAPXEC - 120
         LHPXFE = LHPXEC - 100
      ENDIF
C     LIMITATION A UN SEUL ECRAN
      IF( LAPXFE .GT. 1820 ) LAPXFE = 1820
      IF( LHPXFE .GT. 1000 ) LHPXFE = 1000
C
cccc     pour les traces TAMU
ccc      lapxfe = 1200
ccc      lhpxfe =  900
cccc
cccc     pour projecteur TIMS  Taipei Automne 2009
ccc      LAPXFE = 1000
ccc      LHPXFE = 700
c
ccc     pour faire des images.jpg
cc      LAPXFE = 1024
cc      LHPXFE =  768
C
C     AINSI X11 VA LIRE LE NOMBRE DE PIXELS ECRAN => MARCHE TOUJOURS
      CALL XVINFO( LAPXFE,LHPXFE,MAXFONTS,
     %             N1CORE,NDCORE,N1COEL,NDCOEL,N1COUL,NDCOUL,NBCOLO,
     %             NMFONT,NBCAFO,NBFONT,NCVISU )
C
C     QUELQUES DIMENSIONS ET CARACTERISTIQUES DE LA FENETRE XVUE
C     NOMBRE DE PIXELS EN LARGEUR ET HAUTEUR DE LA FENETRE ACTUELLE
      CALL XVPXFENETRE(LAPXFE,LHPXFE)
C     DIMENSION EN MM DE LA FENETRE OUVERTE
      XMMMIN=0.
C     LARGEUR EN MM DE LA FENETRE OUVERTE
      XMMMAX=LAPXFE/CXMMPX
      YMMMIN=0.
C     HAUTEUR EN MM DE LA FENETRE OUVERTE
      YMMMAX=LHPXFE/CXMMPX
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10000)'PIXELS de la FENETRE',LAPXFE,LHPXFE
         WRITE(IMPRIM,10000)'MM     de la FENETRE',NINT(XMMMAX)
     %                                            ,NINT(YMMMAX)
         WRITE(IMPRIM,*)'En X CONVERSION MM=>PIXELS =',CXMMPX
         WRITE(IMPRIM,*)'En Y CONVERSION MM=>PIXELS =',CYMMPX
      ELSE
         WRITE(IMPRIM,10001)'Window PIXELS',LAPXFE,LHPXFE
         WRITE(IMPRIM,10001)'Window MM    ',NINT(XMMMAX),NINT(YMMMAX)
         WRITE(IMPRIM,*)'X-CONVERSION MM=>PIXELS =',CXMMPX
         WRITE(IMPRIM,*)'Y-CONVERSION MM=>PIXELS =',CYMMPX
      ENDIF
C
C     EPAISSEUR EN M.M. DES TRAITS
      EPAIS = 0.5
C     EPAISSEUR EN PIXELS
      NBEPAI=MAX( NINT(CXMMPX*EPAIS), 1 )
C     NO DU TYPE DE TRAIT : CONTINU PAR DEFAUT
      NOTYTR = 0
C
C     GENERATION DANS LA PARTIE FORTRAN DE LA PALETTE INITIALE DES COULEURS
C     =====================================================================
10002 FORMAT(1X,A,T26,I4,' a ',I4)
10003 FORMAT(1X,A,T23,I4,' to ',I4)
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10002)'COULEURS RESERVEES   de', N1CORE, NDCORE
         WRITE(IMPRIM,10002)'COULEURS DISPONIBLES de', N1COUL, NDCOUL
      ELSE
         WRITE(IMPRIM,10003)'RESERVED COLORS from', N1CORE,NDCORE
         WRITE(IMPRIM,10003)'ALLOWED  COLORS from', N1COUL,NDCOUL
      ENDIF
C
      IF( NBCOLO .EQ. 2 ) THEN
C        NOIR ET BLANC
C        LA COULEUR NOIRE
         PROUGE(0)=0.
         PVERT (0)=0.
         PBLEU (0)=0.
C        LA COULEUR BLANCHE
         PROUGE(1)=1.
         PVERT (1)=1.
         PBLEU (1)=1.
      ELSE
C        CHARGEMENT DANS LA PARTIE FORTRAN DE LA PALETTE DES COULEURS INITIALES
C        DE X STOCKEES DANS LES TABLEAUX RED GREEN BLUE DE LA PARTIE ECRITE EN C
         CALL XVRECUPRGBDEC( NBCOLO, PROUGE, PVERT, PBLEU )
      ENDIF
C
C     LE NOMBRE DE COULEURS DE L'ECRAN
      IF( NDCOUL .NE. N1COUL ) THEN
C        LES 8 COULEURS ELEMENTAIRES ET LEUR NUMERO HP 700
C        LA COULEUR NOIRE EST LA COULEUR DE FOND
CCC         NCNOIR = 0 + N1COEL
         NCNOIR = 0
C        LA COULEUR ROUGE EST LA PREMIERE COULEUR ELEMENTAIRE
         NCROUG = 1 + N1COEL
         NCVERT = 2 + N1COEL
         NCBLEU = 3 + N1COEL
         NCCYAN = 4 + N1COEL
         NCJAUN = 5 + N1COEL
         NCMAGE = 6 + N1COEL
         NCBLAN = 7 + N1COEL
         NCGRIS = 8 + N1COEL
         NCGRIM = 9 + N1COEL
         NCGRIC =10 + N1COEL
         NCBEIG =11 + N1COEL
         NCORAN =12 + N1COEL
         NCSAUM =13 + N1COEL
         NCROSE =14 + N1COEL
         NCTURQ =15 + N1COEL
      ELSE
C        TRAITEMENT EN CAS D'UN ECRAN NOIR ET BLANC
         NCNOIR = 0
         NCROUG = 1
         NCVERT = 1
         NCBLEU = 1
         NCCYAN = 1
         NCJAUN = 1
         NCMAGE = 1
         NCBLAN = 1
         NCGRIS = 1
         NCGRIM = 1
         NCGRIC = 1
         NCBEIG = 1
         NCORAN = 1
         NCSAUM = 1
         NCROSE = 1
         NCTURQ = 1
      ENDIF
C
C     LE FOND DE L'ECRAN EST NOIR
      CALL EFFACE
C
C     CHARGEMENT DANS LA PARTIE FORTRAN DES FONTES DISPONIBLES DANS X
C     ===============================================================
      IF( NBFONT .GT. 0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,*) 'RECUPERATION de ', NBFONT,
     %                      ' FONTES de CARACTERES X11'
         ELSE
            WRITE(IMPRIM,*) 'RECOVERY of ', NBFONT,
     %                      ' X11-CHARACTER FONTS'
         ENDIF
      ELSE
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,*) 'ERREUR: PAS de FONTES de CARACTERES X11'
            WRITE(IMPRIM,*) 'REVOIR L''INSTALLATION DE X11'
         ELSE
            WRITE(IMPRIM,*) 'ERROR: NO X11-CHARACTER FONTS'
            WRITE(IMPRIM,*) 'MODIFY THE INSTALLATION of X11'
         ENDIF
         STOP 1
      ENDIF

      IF( NBFONT .LE. 0 .OR. NBFONT .GT. 512 ) NBFONT=512
      WRITE(IMPRIM,*)
C
C     CARACTERISTIQUES DES FONTES DISPONIBLES
      DO 20 I = NBFONT-1, 0, -1
C
cccC        AFFICHAGE DU NOM DE LA FONTE I
ccc         WRITE(IMPRIM,10004) I, NMFONT(I)(1:NBCAFO(I))
ccc10004    FORMAT('xvinit: FONTE ',I3,' de NOM=',A)
C
C        LA FONTE EST ELLE UNE FONTE USUELLE?
         K = INDEX( NMFONT(I), 'x' )
         IF( K .LE. 1 .OR. K .GT. 3 ) GOTO 20
C
         IF( INDEX( NMFONT(I), '5x7' ) .GT. 0 ) THEN
            NUFOHPX(7) = I
            GOTO 20
         ENDIF
C
         IF( INDEX( NMFONT(I), '5x8' ) .GT. 0 ) THEN
            NUFOHPX(8) = I
            GOTO 20
         ENDIF
C
         IF( INDEX( NMFONT(I), '6x9' ) .GT. 0 ) THEN
            NUFOHPX(9) = I
            GOTO 20
         ENDIF
C
         IF( INDEX( NMFONT(I), '6x10' ) .GT. 0 ) THEN
            NUFOHPX(10) = I
            GOTO 20
         ENDIF
C
         IF( INDEX( NMFONT(I), '6x12' ) .GT. 0 ) THEN
            NUFOHPX(12) = I
            GOTO 20
         ENDIF
C
         IF( INDEX( NMFONT(I), '6x13' ) .GT. 0 ) THEN
            NUFOHPX(13) = I
            GOTO 20
         ENDIF
C
         IF( INDEX( NMFONT(I), '7x13' ) .GT. 0 ) THEN
            NUFOHPX(13) = I
            GOTO 20
         ENDIF
C
         IF( INDEX( NMFONT(I), '7x14' ) .GT. 0 ) THEN
            NUFOHPX(14) = I
            GOTO 20
         ENDIF
C
         IF( INDEX( NMFONT(I), '8x13' ) .GT. 0 ) THEN
            NUFOHPX(13) = I
            GOTO 20
         ENDIF
C
         IF( INDEX( NMFONT(I), '8x16' ) .GT. 0 ) THEN
            NUFOHPX(16) = I
            GOTO 20
         ENDIF
C
         IF( INDEX( NMFONT(I), '9x15' ) .GT. 0 ) THEN
            NUFOHPX(15) = I
            GOTO 20
         ENDIF
C
         IF( INDEX( NMFONT(I), '10x20' ) .GT. 0 ) THEN
            NUFOHPX(20) = I
            GOTO 20
         ENDIF
C
         IF( INDEX( NMFONT(I), '12x24' ) .GT. 0 ) THEN
            NUFOHPX(24) = I
            GOTO 20
         ENDIF
C
 20   CONTINUE
C
C     COMPLETION DU TABLEAU NUFOHPX 
C     NUFOHPX( LHPX ) = NUMERO DE LA FONTE X11 DE HAUTEUR PIXEL LHPX
      NUFOHPX( 1) = NUFOHPX( 7)
      NUFOHPX( 2) = NUFOHPX( 7)
      NUFOHPX( 3) = NUFOHPX( 7)
      NUFOHPX( 4) = NUFOHPX( 7)
      NUFOHPX( 5) = NUFOHPX( 7)
      NUFOHPX( 6) = NUFOHPX( 7)
      NUFOHPX(11) = NUFOHPX(12)
      NUFOHPX(17) = NUFOHPX(16)
      NUFOHPX(18) = NUFOHPX(16)
      NUFOHPX(19) = NUFOHPX(20)
      NUFOHPX(21) = NUFOHPX(20)
      NUFOHPX(22) = NUFOHPX(20)
      NUFOHPX(23) = NUFOHPX(24)
      DO J = 25, MAXHPX
         NUFOHPX(J) = NUFOHPX(24)
      ENDDO
      NOPOCA = 0
      NOFONT = 0
C
ccc      DO J = 1, 32
cccC        NUMERO DE LA FONTE DE J PIXELS EN HAUTEUR
ccc         I = NUFOHPX(J)
ccc         IF( LANGAG .EQ. 0 ) THEN
ccc            WRITE(IMPRIM,10020) J, I, NMFONT(I)(1:NBCAFO(I))
ccc         ELSE
ccc            WRITE(IMPRIM,20020) J, I, NMFONT(I)(1:NBCAFO(I))
ccc         ENDIF
ccc      ENDDO
ccc10020 FORMAT('xvinit: Fonte choisie pour une HAUTEUR PIXELS=',I3,
ccc     %       '  No X11 de la Fonte=',I3,
ccc     %       '  NOM de la Fonte=',A)
ccc20020 FORMAT('xvinit: Chosen Font for a PIXEL HIGH=',I3,
ccc     %       '  X11 Font Number=',I3,
ccc     %       '  Font Name=',A)
C
C     ------------------------------------------------------------------
C     -------------- PARTIE DEPENDANTE de l' ECRAN ---------------------
C     ------------------------------------------------------------------
C     CHOIX DE LA FONTE COURANTE TEL QUE LE MENU TRACMAIL TIENNE DANS LA FENETRE
C     INITIALISATION DE LA VARIABLE NPHFCO EN FONCTION DE LHPXFE
      CALL LAHAFO
      CALL CHOIXFONTE( NPHFCO )
ccc
cccC     IMPRESSION DE LA FONTE CHARGEE NOFONT cf incl/xvfontes.inc
ccc      WRITE(IMPRIM,*)
ccc      WRITE(IMPRIM,*) 'FONTE CHARGEE ',NOFONT,'  ',
ccc     %                 NMFONT(NOFONT)(1:NBCAFO(NOFONT))
C
C     EPAISSEUR ET STYLE DES TRAITS A TRACER EN X
      CALL XVEPAISSEUR( NBEPAI )
      CALL XVTYPETRAIT( NOTYTR )
C
C     VALEURS PAR DEFAUT DES VARIABLES DE $MEFISTO/incl/xyzext.inc
      INIEXT = 0
      MOAXYZ = 0
      XYZAMPLI(1) = 1
      XYZAMPLI(2) = 1
      XYZAMPLI(3) = 1
      XYZAMPLI(4) = 1
C
C     PAS ENCORE DE TRACE DES AXES
      NETAXE = 0
C     PAS ENCORE DE ZOOM ou ORBITE
      NORBITE = 0
C
      RETURN
      END
