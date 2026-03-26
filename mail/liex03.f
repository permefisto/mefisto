      SUBROUTINE LIEX03( NTLXLI, LADEFI, RADEFI,
     %                   NTARLI, MNARLI, NTSOLI, MNSOLI, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERER LES ARETES DE LA LIGNE ARC DE CERCLE
C -----    OPTION 3; DEFINIE PAR 3 POINTS DE R**3 (TYPE 3) OU
C          OPTION 4; DEFINIE PAR 2 POINTS + UN RAYON + 1 POINT DEFINITION DU PLA
C
C ENTREES:
C --------
C LADEFI : TABLEAU ENTIER DE DEFINITION DE LA LIGNE
C          CF ~/td/d/a_ligne__definition
C RADEFI : TABLEAU REEL   DE DEFINITION DE LA LIGNE
C          CES 2 TABLEAUX LADEFI ET RADEFI ONT MEME ADRESSE A L'APPEL
C NTLXLI : NUMERO DU TABLEAU TS DU LEXIQUE DE LA LIGNE
C
C SORTIES:
C --------
C NTARLI : NUMERO      DU TMS 'NSEF' DES NUMEROS DES ARETES DE LA LIGNE
C MNARLI : ADRESSE MCN DU TMS 'NSEF' DES NUMEROS DES ARETES DE LA LIGNE
C          CF ~/td/d/a___nsef
C NTSOLI : NUMERO      DU TMS 'XYZSOMMET' DE LA LIGNE
C MNSOLI : ADRESSE MCN DU TMS 'XYZSOMMET' DE LA LIGNE
C          CF ~/td/d/a___xyzsommet
C IERR   : 0 SI PAS D'ERREUR
C          1 SI NOMBRE DE SOMMETS INCORRECT <2
C          2 POINT INITIAL OU FINAL NON INITIALISE
C          3 POINT INITIAL ET FINAL CONFONDUS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC  PARIS  SEPTEMBRE 1996
C2345X7..............................................................012
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/langue.inc"
      include"./incl/a_ligne__definition.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/gsmenu.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
C
      INTEGER           LADEFI(0:*)
      REAL              RADEFI(0:*)
C
      REAL              XYZ1(3), XYZ2(3), XYZ3(3), XYZ(3), XYZ4(3),
     %                  XY1(2),  XY2(2),  XY3(2), XY4(2), TG(3),
     %                  CENTR2(3), CENTRE(3)
      DOUBLE PRECISION  D2D3(3,3), TETA2D, DBLE, DEUXPI, ANGTOT, ARC,
     %                  ANGLE, ANGLE0, RAYON, A, B, C, D, H, A2, A3, HHH
      DOUBLE PRECISION  DTAILL, DTAIL2, DTAIL1, XYZD(3)
C
C     LE TYPE DE LA LIGNE
      NUTYLI = LADEFI( WUTYLI )
      IF( NUTYLI .NE. 3 .AND. NUTYLI .NE. 4 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'LIEX03: TYPE DE LIGNE INCORRECT'
         ELSE
            KERR(1) = 'LIEX03: INCORRECT TYPE of LINE'
         ENDIF
         CALL LEREUR
         IERR = 1
         GOTO 9999
      ENDIF
C
C     LE NOMBRE DES ARETES ET DES SOMMETS DE LA LIGNE
      NBARLI = LADEFI( WBARLI )
C
C     LA RAISON DE LA PROGRESSION GEOMETRIQUE
      RAIGEO = RADEFI( WAIGEO )
C
C     LE NUMERO UTILISATEUR DES POINTS EXTREMITES DE LA LIGNE
      NUPTIN = LADEFI( WUPTIN )
      NUPTFI = LADEFI( WUPTFI )
C
C     LE NUMERO DU POINT ENTRE LES EXTREMITES OU DE DEFINITION DU PLAN
      IF( NUTYLI .EQ. 3 ) THEN
         NUPTEN = LADEFI( WUPTEN )
      ELSE
         NUPTEN = LADEFI( WUPTAR )
      ENDIF
C
C     EXISTENCE OU NON DE LA FONCTION 'TAILLE_IDEALE(XYZ)' DES ARETES
C     ===============================================================
C     ICI LA CARTE EST SUPPOSEE ISOTROPE
      NOFOTI = NOFOTIEL()
C     NOFOTI>0 SI CETTE FONCTION EXISTE
      IF( NOFOTI .GT. 0 ) THEN
C        LA FONCTION 'TAILLE_IDEALE(X,Y,Z)' EXISTE
         IF( NBARLI .LT. 2 ) THEN
            NBARLI = 1
            RAIGEO = 1.0
         ENDIF
      ENDIF
      NBSOLI = NBARLI + 1
C
C     DONNEES CORRECTES?
      IF( NBSOLI .LT. 2 ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(2),'(I12)') NBSOLI-1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = KERR(2)(1:12) // 'NOMBRE INCORRECT D''ARETES'
         ELSE
            KERR(1) = KERR(2)(1:12) // 'INCORRECT EDGES NUMBER'
         ENDIF
         CALL LEREUR
         IERR = 2
         GOTO 9999
      ELSE IF( NUPTIN .LE. 0 .OR. NUPTFI .LE. 0 .OR.
     %         NUPTEN .LE. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = '1 DES 3 POINTS DE L''ARC N''EXISTE PAS'
         ELSE
            KERR(1) = '1 of 3 POINTS of the CIRCLE DOES NOT EXIST'
         ENDIF
         CALL LEREUR
         IERR = 3
         GOTO 9999
      ENDIF
C
C     COORDONNEES DU POINT INITIAL
      CALL LXNLOU( NTPOIN, NUPTIN, NTPOI, MN )
      CALL LXTSOU( NTPOI , 'XYZSOMMET',  NTSOM, MNSOM )
      MN = MNSOM + WYZSOM
      XYZ1( 1 ) = RMCN( MN )
      XYZ1( 2 ) = RMCN( MN + 1 )
      XYZ1( 3 ) = RMCN( MN + 2 )
C
C     COORDONNEES DU POINT FINAL
      CALL LXNLOU( NTPOIN, NUPTFI, NTPOI, MN )
      CALL LXTSOU( NTPOI , 'XYZSOMMET',  NTSOM, MNSOM )
      MN = MNSOM + WYZSOM
      XYZ2( 1 ) = RMCN( MN )
      XYZ2( 2 ) = RMCN( MN + 1 )
      XYZ2( 3 ) = RMCN( MN + 2 )
C
C     COORDONNEES DU POINT INTERMEDIAIRE
      CALL LXNLOU( NTPOIN, NUPTEN, NTPOI, MN )
      CALL LXTSOU( NTPOI , 'XYZSOMMET',  NTSOM, MNSOM )
      MN = MNSOM + WYZSOM
      XYZ3( 1 ) = RMCN( MN )
      XYZ3( 2 ) = RMCN( MN + 1 )
      XYZ3( 3 ) = RMCN( MN + 2 )
C
C     DISTANCE ENTRE CES 2 POINTS
      H = DIST2P( XYZ1, XYZ2 )
      IF( H .LE. 0. ) THEN
         NBLGRC(NRERR) = 2
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'POINTS INITIAL et FINAL CONFONDUS'
            KERR(2) = 'PAS DE CERCLE SOLUTION'
         ELSE
            KERR(1) = 'INITIAL and FINAL POINTS are SAME'
            KERR(2) = 'NO SOLUTION CIRCLE'
         ENDIF
         CALL LEREUR
         IERR = 3
         GOTO 9999
      ENDIF
C
C     LE RAYON DU CERCLE POUR LIGNE DE TYPE 4
      RAYON = 0
      IF( NUTYLI .EQ. 4 ) THEN
         RAYON = RADEFI( WAYARC )
         IF( RAYON-H*0.5 .LE. 1E-3*H ) THEN
            NBLGRC(NRERR) = 2
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'RAYON DU CERCLE <= DISTANCE(P1,P2)/2'
               KERR(2) = 'PAS DE CERCLE SOLUTION'
            ELSE
               KERR(1) = 'CIRCLE RADIUS <= DISTANCE(P1,P2)/2'
               KERR(2) = 'NO SOLUTION CIRCLE'
            ENDIF
            CALL LEREUR
            IERR = 4
            GOTO 9999
         ENDIF
      ENDIF
C
      H = DIST2P( XYZ1, XYZ3 )
      IF( H .LE. 0. ) THEN
         NBLGRC(NRERR) = 2
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'POINTS INITIAL ET INTERMEDIAIRE CONFONDUS'
            KERR(2) = 'PAS DE CERCLE SOLUTION'
         ELSE
            KERR(1) = 'INITIAL and INTERMEDIATE POINTS are SAME'
            KERR(2) = 'NO SOLUTION CIRCLE'
         ENDIF
         CALL LEREUR
         IERR = 3
         GOTO 9999
      ENDIF
C
      H = DIST2P( XYZ3, XYZ2 )
      IF( H .LE. 0. ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'POINTS FINAL ET INTERMEDIAIRE CONFONDUS'
         ELSE
            KERR(1) = 'FINAL and INTERMEDIATE POINTS are SAME'
         ENDIF
         CALL LEREUR
         IERR = 3
         GOTO 9999
      ENDIF
C
C     PASSAGE AUX COORDONNEES DES 3 POINTS DANS LEUR PLAN
C     ---------------------------------------------------
      CALL DF3D2D( XYZ1, XYZ3, XYZ2, D2D3, IERR )
      IF( IERR .GT. 0 ) THEN
         NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = '3 POINTS DU CERCLE COLINEAIRES'
         ELSE
            KERR(1) = '3 POINTS of CIRCLE ON A STRAIGHT LINE'
         ENDIF
         CALL LEREUR
         IERR = 4
         GOTO 9999
      ENDIF
C     LE POINT INITIAL A POUR COORDONNEES 2D 0. 0.
      XY1(1) = 0.
      XY1(2) = 0.
      CALL CH3D2D( XYZ1, D2D3, XYZ2, XY2 )
      CALL CH3D2D( XYZ1, D2D3, XYZ3, XY3 )
C
C     LES COORDONNEES 2D DU CENTRE DU CERCLE ET RAYON
C     -----------------------------------------------
      IF( NUTYLI .EQ. 3 ) THEN
C
C        CERCLE DEFINI PAR 3 POINTS DE R**3
         CALL CENCER( XY1, XY2, XY3, CENTR2 )
         RAYON = SQRT( CENTR2( 3 ) )
C
      ELSE
C
C        CERCLE DEFINI PAR 2 POINTS + RAYON + POINT DU PLAN
         IF( ABS(XY2(2)) .LE. 1E-4*RAYON ) THEN
C
C           LES 2 POINTS SONT HORIZONTAUX
            CENTR2(1) = XY2(1) * 0.5
            A2 = (XY3(1)-CENTR2(1))**2 + (XY3(2)-RAYON)**2
            A3 = (XY3(1)-CENTR2(1))**2 + (XY3(2)+RAYON)**2
            H  = (XY2(1)-XY1(1))**2 + (XY2(2)-XY1(2))**2
            C  = SQRT( RAYON**2 - H * 0.25 )
            IF( A2 .LT. A3 ) THEN
               CENTR2(2) =  REAL( C )
            ELSE
               CENTR2(2) = REAL( -C )
            ENDIF
C
         ELSE
C
C           LES 2 POINTS NE SONT PAS HORIZONTAUX
            H  = ( XY1(2) - XY2(2) )
            A2 = ( XY2(1) - XY1(1) ) / H
            A3 = ( H * (XY1(2)+XY2(2))+(XY1(1)-XY2(1))*(XY1(1)+XY2(1)) )
     %         / (2 * H)
            A = 1 + A2 * A2
            B = A2 * (A3-XY1(2)) - XY1(1)
            C = XY1(1)**2 + (A3-XY1(2))**2 - RAYON**2
            D = B * B - A * C
            IF( B .GT. 0 ) THEN
               CENTR2(1) = REAL( ( - B - SQRT(D) ) / A )
            ELSE
               CENTR2(1) = REAL( ( - B + SQRT(D) ) / A )
            ENDIF
C           LA SECONDE RACINE
            CENTR2(2) = REAL( C / A * CENTR2(1) )
C           CHOIX DU CENTRE COMME LE PLUS PROCHE DE P3
            A = A2 * CENTR2(1) + A3
            B = A2 * CENTR2(2) + A3
C
            A2 = (XY3(1)-CENTR2(1))**2 + (XY3(2)-A)**2
            A3 = (XY3(1)-CENTR2(2))**2 + (XY3(2)-B)**2
            IF( A2 .LT. A3 ) THEN
C              L'ORDONNEE DU CENTRE DU CERCLE
               CENTR2(2) =  REAL( A )
            ELSE
C              L'ABSCISSE ET L'ORDONNEE DU CENTRE DU CERCLE
               CENTR2(1) = CENTR2(2)
               CENTR2(2) = REAL( B )
            ENDIF
C
         ENDIF
      ENDIF
C
C     LES COORDONNEES DU CENTRE EN 3D
C     -------------------------------
      CALL CH2D3D( XYZ1, D2D3, CENTR2, CENTRE )
C
C     PASSAGE AUX COORDONNEES DEFINIES PAR LES 3 POINTS: CENTRE, P1, P3
      CALL DF3D2D( CENTRE, XYZ1, XYZ3, D2D3, IERR )
      IF( IERR .GT. 0 ) THEN
C        LE CENTRE ET P1 P3 SONT ALIGNES => ESSAI CENTRE, P1, P2
         CALL DF3D2D( CENTRE, XYZ1, XYZ2, D2D3, IERR )
      ENDIF
C
C     LES COORDONNEES DANS CE SYSTEME D'AXES CENTRE P1 P2 OU P3
      CALL CH3D2D( CENTRE, D2D3, XYZ2, XY2 )
C
C     CALCUL DE L'ANGLE CENTRE P1,CENTRE P2
      ANGTOT = TETA2D( DBLE( XY2(1) ), DBLE( XY2(2) ) )
C
C     2 FOIS PI
      DEUXPI = ATAN( 1D0 ) * 8D0
      LESIGN = 1
C
      IF( NUTYLI .EQ. 3 ) THEN
C        CALCUL DE L'ANGLE CENTRE P1,CENTRE P3
         CALL CH3D2D( CENTRE, D2D3, XYZ3, XY3 )
         A3 = TETA2D( DBLE( XY3(1) ), DBLE( XY3(2) ) )
         IF( ANGTOT .LT. A3 ) THEN
C           LE PAS DE L'ANGLE VA ETRE NEGATIF
            ANGTOT = ANGTOT - DEUXPI
            LESIGN = -1
         ENDIF
      ELSE
C        RECHERCHE DU PLUS PETIT DES ANGLES
         A      = DEUXPI - ANGTOT
         ANGTOT = MIN( ANGTOT, A )
      ENDIF
C
C     PROGRESSION GEOMETRIQUE OU NON ?
      IF( ABS(RAIGEO-1.0) .LT. 1.0E-3 ) THEN
C
C        ANGLES CONSTANTS
         A      = ANGTOT / NBARLI
         RAIGEO = 1.0
         NBTGS  = NBSOLI
C
      ELSE
C
C        ANGLES EN PROGRESSION GEOMETRIQUE. LE PREMIER ANGLE A
         A = ANGTOT * ( 1. - RAIGEO ) / ( 1. - RAIGEO ** NBARLI )
         NBTGS = 2 * NBARLI
      ENDIF
C
C     CALCUL DE LA LONGUEUR DE L'ARC TOTAL
      ARC = ABS( RAYON * ANGTOT )
C
C     EXISTENCE OU NON DE LA FONCTION 'TAILLE_IDEALE(XYZ)' DES ARETES
C     ===============================================================
      IF( NOFOTI .GT. 0 ) THEN
C
C        LA FONCTION 'TAILLE_IDEALE' EXISTE
C        DECLARATION DU TABLEAU DES TAILLES DES ARETES
         MXTAIL = 512
         CALL TNMCDC( 'REEL', MXTAIL, MNTAIL )
C        CALCUL DU NOMBRE D'ARCS DE CERCLE POUR DECLARER LES TABLEAUX
         NBARLI  = 1
         HHH     = ARC
         ANGLE   = 0
         XYZD(1) = XYZ1(1)
         XYZD(2) = XYZ1(2)
         XYZD(3) = XYZ1(3)
C
C        LES 3 PARAMETRES D'APPEL DE LA FONCTION 'TAILLE_IDEALE'
C        AU SOMMET XYZ DU SEGMENT SOMMET1 SOMMET2
 5       CALL FONVAL( NOFOTI, 3, XYZD,  NCODEV, DTAILL )
         IF( NCODEV .EQ. 0 ) THEN
C
C           NCODEV  : 0 DTAILL N'EST PAS INITIALISEE EN SORTIE
C                     1 DTAILL   EST     INITIALISEE EN SORTIE
            NBLGRC(NRERR) = 2
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = 'TAILLE_IDEALE(x,y,z) INCALCULABLE AU POINT'
            ELSE
               KERR(1) = 'EDGE_LENGTH(x,y,z) NOT COMPUTABLE at POINT'
            ENDIF
            KERR(2) = ' '
            WRITE(KERR(2)(1:15), '(E15.7)') XYZD(1)
            WRITE(KERR(2)(18:32),'(E15.7)') XYZD(2)
            WRITE(KERR(2)(35:49),'(E15.7)') XYZD(3)
            CALL LEREUR
            IERR = 1
            GOTO 9999
C
         ELSE
C
C           TAILLE AUTOUR DU SOMMET EST INITIALISEE
C           ESSAI DE CREER UN ARC DE LONGUEUR CETTE TAILLE
            DTAILL = ABS( DTAILL )
            IF( DTAILL .LT. HHH*0.65 ) THEN
C              CREATION D'UN POINT INTERMEDIAIRE => UNE ARETE DE PLUS
               IF( NBARLI .GE. MXTAIL-1 ) THEN
C                 SATURATION DU TABLEAU TAIL => IL EST AUGMENTE
                  CALL TNMCAU( 'REEL', MXTAIL, MXTAIL+512,
     %                         NBARLI-1, MNTAIL )
                  MXTAIL = MXTAIL + 512
               ENDIF
               NBARLI = NBARLI + 1
               ANGLE0 = ANGLE
               IPAS   = 0
C
C              LES COORDONNEES 2D DU POINT FINAL DE L'ARC
 7             ANGLE  = ANGLE0 + LESIGN * DTAILL / RAYON
               XY2(1) = REAL( RAYON * COS( ANGLE ) )
               XY2(2) = REAL( RAYON * SIN( ANGLE ) )
C              RETOUR AUX COORDONNEES 3D
               CALL CH2D3D( CENTRE, D2D3, XY2, XYZ3 )
               IF( IPAS .EQ. 0 ) THEN
C
C                 LA TAILLE IDEALE EN LA SECONDE EXTREMITE DE L'ARC
                  XYZD(1) = XYZ3(1)
                  XYZD(2) = XYZ3(2)
                  XYZD(3) = XYZ3(3)
                  CALL FONVAL( NOFOTI, 3, XYZD,  NCODEV, DTAIL2 )
                  DTAIL2 = ABS( DTAIL2 )
C
C                 LA TAILLE SOUHAITEE AU MILIEU DE L'ARC
                  A  = ANGLE0 + LESIGN * DTAILL / RAYON / 2
                  XY4(1) = REAL( RAYON * COS( A ) )
                  XY4(2) = REAL( RAYON * SIN( A ) )
C                 RETOUR AUX COORDONNEES 3D
                  CALL CH2D3D( CENTRE, D2D3, XY4, XYZ4 )
                  XYZD(1) = XYZ4(1)
                  XYZD(2) = XYZ4(2)
                  XYZD(3) = XYZ4(3)
                  CALL FONVAL( NOFOTI, 3, XYZD,  NCODEV, DTAIL1 )
                  DTAIL1 = ABS( DTAIL1 )
C
C                 STRATEGIE POUR LE CALCUL DE LA TAILLE FINALE
                  IF( DTAIL1 .LT. DTAILL .AND. DTAIL1 .LT. DTAIL2 ) THEN
C                    LA TAILLE AU MILIEU EST INFERIEURE A CELLES DES EXTREMITES
C                    CALCUL A NOUVEAU AVEC CETTE NOUVELLE TAILLE CAR MINIMUM LOC
                     DTAILL = DTAIL1
                     GOTO 7
                  ENDIF
C                 LA TAILLE MINIMALE EST CELLE D'UNE DES 2 EXTREMITES
C                 LA PLUS FAIBLE EST GARDEE
                  IF( DTAIL2 .LT. DTAILL ) THEN
C                    PAS DE NOUVEAU PASSAGE
                     IPAS   = 1
                     DTAILL = DTAIL2
                     GOTO 7
                  ENDIF
               ENDIF
C
C              LA TAILLE DE L'ARETE NBARLI-1
               RMCN( MNTAIL-2+NBARLI ) = REAL( LESIGN * DTAILL )
C              LA LONGUEUR RESTANTE
               HHH = HHH - DTAILL
C              LE SOMMET SUR LE CERCLE
               XYZD(1) = XYZ3(1)
               XYZD(2) = XYZ3(2)
               XYZD(3) = XYZ3(3)
               GOTO 5
            ENDIF
C
C           LA TAILLE DE L'ARETE NBARLI
            RMCN( MNTAIL-1+NBARLI ) = REAL( RAYON * ( ANGTOT-ANGLE ) )
         ENDIF
C
         IF( NBARLI .LE. 0 ) THEN
            NBLGRC(NRERR) = 1
            WRITE(KERR(2),'(I12)') NBARLI
            IF( LANGAG .EQ. 0 ) THEN
               KERR(1) = KERR(2)(1:12)//' => NOMBRE INCORRECT D''ARETES'
            ELSE
               KERR(1) = KERR(2)(1:12)//' => INCORRECT EDGES NUMBER'
            ENDIF
            CALL LEREUR
            IERR = 2
            GOTO 9999
         ENDIF
C
C        LE NOMBRE DE SOMMETS ET DE TANGENTES DE LA LIGNE
         NBSOLI = NBARLI + 1
         NBTGS  = 2 * NBARLI
      ENDIF
C
C     CONSTRUCTION DU TABLEAU 'XYZSOMMET'
C     -----------------------------------
      CALL LXTNDC( NTLXLI, 'XYZSOMMET', 'ENTIER',
     %                      WYZSOM + 3*(NBSOLI + NBTGS) )
      CALL LXTSOU( NTLXLI, 'XYZSOMMET',  NTSOLI, MNSOLI )
C
C     LE NOMBRE DE SOMMETS
      MCN( MNSOLI + WNBSOM ) = NBSOLI
C
C     LE NOMBRE DE TANGENTES
      MCN( MNSOLI + WNBTGS ) = NBTGS
C
C     LE PREMIER POINT DE L'ARC DE CERCLE = LE POINT INITIAL
C     ===================================
C     ADRESSE DU DEBUT DES COORDONNEES DU 1-ER SOMMET DE LA LIGNE
      MNS  = MNSOLI + WYZSOM
      MNTG = MNS + 3 * NBSOLI
C
C     AJOUT DES COORDONNEES DU POINT INITIAL EN (RAYON,0) DANS R**2
      RMCN( MNS     ) = XYZ1( 1 )
      RMCN( MNS + 1 ) = XYZ1( 2 )
      RMCN( MNS + 2 ) = XYZ1( 3 )
      MNS = MNS + 3
C
C     LES COORDONNEES DES SOMMETS ET TANGENTES DE L'ARC
C     =================================================
      IF( NOFOTI .EQ. 0 ) THEN
C
C        LA FONCTION 'TAILLE_IDEALE(X,Y,Z)' N' EXISTE PAS
C        ************************************************
C        LE POINT AU BOUT DE LA TANGENTE PLACEE AU POINT INITIAL DE L'ARC
         XY1(1) = REAL( RAYON )
         XY1(2) = REAL( RAYON * ANGTOT )
C
C        RETOUR AUX COORDONNEES 3D (TRANSFORMATION AFFINE)
         CALL CH2D3D( CENTRE, D2D3, XY1, XYZ )
C        LES COMPOSANTES DE LA TG = TRANSFORMEE R2 R3(P1+TG) -TRANSFORMEE R2 R3(
         TG(1) = XYZ(1) - XYZ1(1)
         TG(2) = XYZ(2) - XYZ1(2)
         TG(3) = XYZ(3) - XYZ1(3)
C
         IF( NBTGS .EQ. NBSOLI ) THEN
C           ARETES DE MEME LONGUEUR . CALCUL DE LA TANGENTE AU POINT INITIAL
C           CORRECTION DE LA DERIVEE EN MULTIPLIANT PAR H
            H = A / ANGTOT
            RMCN( MNTG     ) = REAL( TG(1) * H )
            RMCN( MNTG + 1 ) = REAL( TG(2) * H )
            RMCN( MNTG + 2 ) = REAL( TG(3) * H )
            MNTG = MNTG + 3
         ENDIF
C
         ANGLE0 = 0
         A3     = 1
C
         DO 20 I=1,NBARLI
C
C           L'ANGLE / OX DU POINT FINAL DE L'ARETE I
            ANGLE = ANGLE0 + A * A3
C
C           CALCUL DE LA LONGUEUR DU SOUS ARC I
            H  = A * A3 / ANGTOT
C
            IF( NBTGS .NE. NBSOLI ) THEN
C              PROGRESSION GEOMETRIQUE => 2 TGS PAR ARETE
C              LA TG AU SOMMET 1 DE L'ARETE A PARTIR DE LA TG
C              DU SOMMET 2 DE L'ARETE PRECEDENTE
C              CORRECTION DE LA DERIVEE EN MULTIPLIANT PAR H
               RMCN( MNTG     ) = REAL( TG(1) * H )
               RMCN( MNTG + 1 ) = REAL( TG(2) * H )
               RMCN( MNTG + 2 ) = REAL( TG(3) * H )
               MNTG = MNTG + 3
            ENDIF
C
C           LES COORDONNEES 2D DU POINT FINAL DE L'ARETE I
            XY2(1) = REAL( RAYON * COS( ANGLE ) )
            XY2(2) = REAL( RAYON * SIN( ANGLE ) )
C           RETOUR AUX COORDONNEES 3D
            CALL CH2D3D( CENTRE, D2D3, XY2, XYZ3 )
C           RANGEMENT DES 3 COORDONNEES
            RMCN( MNS     ) = XYZ3( 1 )
            RMCN( MNS + 1 ) = XYZ3( 2 )
            RMCN( MNS + 2 ) = XYZ3( 3 )
            MNS = MNS + 3
C
C           LE POINT AU BOUT DE LA TANGENTE PLACEE AU POINT INITIAL DE L'ARC
C           AVEC LA CORRECTION PAR H LA LONGUEUR DU SOUS INTERVALLE DU PARAMETRE
            XY3(1) = REAL( XY2(1) - RAYON * ANGTOT * SIN(ANGLE) )
            XY3(2) = REAL( XY2(2) + RAYON * ANGTOT * COS(ANGLE) )
C
C           RETOUR AUX COORDONNEES 3D (TRANSFORMATION AFFINE)
            CALL CH2D3D( CENTRE, D2D3, XY3, XYZ )
C           LES COMPOSANTES DE LA TG = TRANSFORMEE R2 R3(P1+TG) -TRANSFORMEE R2
            TG(1) = XYZ(1) - XYZ3(1)
            TG(2) = XYZ(2) - XYZ3(2)
            TG(3) = XYZ(3) - XYZ3(3)
C
C           STOCKAGE DE LA TANGENTE AU SOMMET 2 CORRIGEE PAR  * H
            RMCN( MNTG     ) = REAL( TG(1) * H )
            RMCN( MNTG + 1 ) = REAL( TG(2) * H )
            RMCN( MNTG + 2 ) = REAL( TG(3) * H )
            MNTG = MNTG + 3
C
C           LE POINT FINAL DE L'ARETE DEVIENT LE POINT INITIAL DE L'ARETE SUIVAN
            XYZ1(1) = XYZ3(1)
            XYZ1(2) = XYZ3(2)
            XYZ1(3) = XYZ3(3)
            ANGLE0  = ANGLE
C
C           LE PRODUIT DES RAISONS GEOMETRIQUES SELON LE RANG
            A3 = A3 * RAIGEO
 20      CONTINUE
C
      ELSE
C
C        LA FONCTION 'TAILLE_IDEALE(X,Y,Z)' EXISTE
C        *****************************************
C        LES AUTRES SOMMETS INTERMEDIAIRES
         ANGLE = 0
         DO 25 I=1,NBARLI
C
C           LA TAILLE DE L'ARETE I
            DTAILL = RMCN(MNTAIL-1+I)
C
C           L'ANGLE DE L'ARC NOUVELLE ARETE
            A = DTAILL / RAYON
C
C           LES COORDONNEES 2D DU POINT INITIAL DE L'ARC
            XY2(1) = REAL( RAYON * COS( ANGLE ) )
            XY2(2) = REAL( RAYON * SIN( ANGLE ) )
C
C           LE POINT AU BOUT DE LA TANGENTE PLACEE AU POINT INITIAL DE L'ARC
            XY3(1) = REAL( XY2(1) - DTAILL * SIN(ANGLE) )
            XY3(2) = REAL( XY2(2) + DTAILL * COS(ANGLE) )
C
C           RETOUR AUX COORDONNEES 3D (TRANSFORMATION AFFINE)
            CALL CH2D3D( CENTRE, D2D3, XY3, XYZ )
C
C           LES COMPOSANTES DE LA TG = TRANSFORMEE R2 R3(P1+TG) -TRANSFORMEE R2
            RMCN( MNTG     ) = XYZ(1) - RMCN(MNS-3)
            RMCN( MNTG + 1 ) = XYZ(2) - RMCN(MNS-2)
            RMCN( MNTG + 2 ) = XYZ(3) - RMCN(MNS-1)
            MNTG = MNTG + 3
C
C           L'ANGLE FINAL DU SOMMET DE L'ARC AVEC OX EN 2D
            ANGLE = ANGLE + A
C
C           LES COORDONNEES 2D DU POINT FINAL DE L'ARC DEJA CALCULEES DANS XY1
C           XY1: LES COORDONNEES 2D DU POINT FINAL DE L'ARC
            XY1(1) = REAL( RAYON * COS( ANGLE ) )
            XY1(2) = REAL( RAYON * SIN( ANGLE ) )
C
C           RETOUR AUX COORDONNEES 3D
            CALL CH2D3D( CENTRE, D2D3, XY1, RMCN(MNS) )
C
C           LE POINT AU BOUT DE LA TANGENTE PLACEE AU POINT FINAL DE L'ARC
            XY3(1) = REAL( XY1(1) - DTAILL * SIN(ANGLE) )
            XY3(2) = REAL( XY1(2) + DTAILL * COS(ANGLE) )
C
C           RETOUR AUX COORDONNEES 3D (TRANSFORMATION AFFINE)
            CALL CH2D3D( CENTRE, D2D3, XY3, XYZ )
C           LES COMPOSANTES DE LA TG = TRANSFORMEE R2 R3(P1+TG) -TRANSFORMEE R2
            RMCN( MNTG     ) = XYZ(1) - RMCN(MNS)
            RMCN( MNTG + 1 ) = XYZ(2) - RMCN(MNS+1)
            RMCN( MNTG + 2 ) = XYZ(3) - RMCN(MNS+2)
            MNTG = MNTG + 3
C
C           LE SOMMET SUR LE CERCLE
            XYZD(1) = RMCN(MNS)
            XYZD(2) = RMCN(MNS+1)
            XYZD(3) = RMCN(MNS+2)
            MNS  = MNS + 3
 25      CONTINUE
C
C        DESTRUCTION DU TABLEAU DES TAILLES
         CALL TNMCDS( 'REEL', MXTAIL, MNTAIL )
      ENDIF
C
C     LE DERNIER POINT DE L'ARC DE CERCLE
C     ECRASEMENT DES 3 DERNIERES COORDONNEES PAR CELLES EXACTES
      MNS = MNSOLI + WYZSOM + 3*NBSOLI - 3
      RMCN( MNS     ) = XYZ2( 1 )
      RMCN( MNS + 1 ) = XYZ2( 2 )
      RMCN( MNS + 2 ) = XYZ2( 3 )
C
C     AJOUT DE LA DATE
      CALL ECDATE( MCN(MNSOLI) )
C
C     AJOUT DU NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNSOLI + WBCOOR ) = 3
      MCN( MNSOLI + MOTVAR(6) ) = NONMTD( '~>>>XYZSOMMET' )
C
C     CONSTRUCTION DU TABLEAU 'NSEF' LIGNE STRUCTURE
C     ----------------------------------------------
      CALL LXTNDC( NTLXLI, 'NSEF', 'ENTIER', 1 + WBARSE + 4*NBARLI )
      CALL LXTSOU( NTLXLI, 'NSEF',  NTARLI, MNARLI )
C
C     LE TYPE DE L'OBJET : ICI LIGNE
      MCN( MNARLI + WUTYOB ) = 2
C     LE TYPE NON-FERME DE FERMETURE DU MAILLAGE
      MCN( MNARLI + WUTFMA ) = 0
C     LE NOMBRE DE SOMMETS PAR EF
      MCN( MNARLI + WBSOEF ) = 2
C     LE NOMBRE DE TANGENTES PAR EF
      MCN( MNARLI + WBTGEF ) = 2
C     LE NOMBRE D'ARETES DU SEGMENT STRUCTURE
      MCN( MNARLI + WBEFOB ) = NBARLI
C     LE NOMBRE D'EF AVEC TANGENTES
      MCN( MNARLI + WBEFTG ) = NBARLI
C     LE NOMBRE D'EF AVEC POINTEUR = NBEFTG = NBARLI
      MCN( MNARLI + WBEFAP ) = NBARLI
C     LE TYPE DU MAILLAGE : ICI SEGMENT STRUCTURE
      MCN( MNARLI + WUTYMA ) = 2
C     LE NOMBRE D'ARETES DU SEGMENT STRUCTURE
      MCN( MNARLI + WBARSE ) = NBARLI
C
C     LE NUMERO DE L'EF AVEC TANGENTES
      MNS = MNARLI + WBARSE
      DO 30 I=1,NBARLI
C        LE POINTEUR SUR LES EF A TG : L'ARETE I EST L'EF I AVEC 2 TGS
         MCN( MNS + I ) = I
C        LE NUMERO GEOMETRIQUE EST ICI LE CERCLE : 1
         MCN( MNS + NBARLI + I ) = 1
 30   CONTINUE
      MNS = MNS + 2 * NBARLI + 1
C
C     LE NUMERO DES TGS DE CHAQUE ARETE
      IF( NBTGS .NE. NBSOLI ) THEN
C        PROGRESSION GEOMETRIQUE => 2 TGS PAR ARETE
         K = 1
         DO 60 I=1,NBARLI
C           NUMERO DE LA TANGENTE AU SOMMET INITIAL DE L'ARETE I
            MCN(MNS  ) = K
C           -NUMERO DE LA TANGENTE AU SOMMET FINAL  DE L'ARETE I
            MCN(MNS+1) = -(K+1)
            MNS = MNS + 2
            K   = K + 2
 60      CONTINUE
      ELSE
C        UNE TANGENTE EN CHAQUE SOMMET
         DO 70 I=1,NBARLI
C           NUMERO DE LA TANGENTE AU SOMMET INITIAL DE L'ARETE I
            MCN(MNS  ) = I
C           -NUMERO DE LA TANGENTE AU SOMMET FINAL  DE L'ARETE I
            MCN(MNS+1) = -(I+1)
            MNS = MNS + 2
 70      CONTINUE
      ENDIF
C
C     AJOUT DE LA DATE
      CALL ECDATE( MCN(MNARLI) )
C
C     AJOUT DU NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNARLI + MOTVAR(6) ) = NONMTD( '~>>>NSEF' )
      IERR = 0
C
C     ERREUR
C     ======
 9999 RETURN
      END
