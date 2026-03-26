      SUBROUTINE LIEX08( NTLXLI, LADEFI, RADEFI,
     %                   NTARLI, MNARLI, NTSOLI, MNSOLI, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERER LES ARETES D'UN CERCLE de R3 ( TYPE DE LIGNE 8 )
C -----    GEOMETRIE C1
C
C ENTREES:
C --------
C LADEFI : TABLEAU ENTIER DE DEFINITION DE LA LIGNE
C RADEFI : TABLEAU REEL   DE DEFINITION DE LA LIGNE
C          CF ~/td/d/a_ligne__definition
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
C          > 0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC  PARIS  SEPTEMBRE 1996
C2345X7..............................................................012
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      include"./incl/gsmenu.inc"
      include"./incl/a_ligne__definition.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/ntmnlt.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
C
      INTEGER           LADEFI(0:*)
      REAL              RADEFI(0:*)
C
      INTEGER           NUPT(1:3)
      EQUIVALENCE      (NUPT(1),NUPT1), (NUPT(2),NUPT2), (NUPT(3),NUPT3)
      REAL              XYZ(3,3), TG(3),
     %                  XY(2,3), XY1(2), XYZ1(3),
     %                  CENRAY(3), CENTRE(3)
      DOUBLE PRECISION  PI, PI2, D2D3(3,3), XYZD(3),
     %                  DTAILL, DTAIL1, DTAIL2, DEUPIR, ANGLE, ANGLE0,
     %                  A, HHH, ANGDER
      DATA              PI/ 3.14159265358979312D0 /
C
C     LE TYPE DE LA LIGNE
      N = LADEFI( WUTYLI )
      IF( N .NE. 8 ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:4),'(I4)') N
         KERR(1) =  'LIEX08:TYPE INCORRECT '//KERR(MXLGER)(1:4)
         CALL LEREUR
         IERR = 1
         GOTO 9999
      ENDIF
C
C     LE NOMBRE D'ARETES ET DE SOMMETS DE LA LIGNE
      NBARLI = LADEFI( WBARLI )
C
C     EXISTENCE OU NON DE LA FONCTION 'TAILLE_IDEALE(XYZ)' DES ARETES
C     ===============================================================
C     ICI LA CARTE EST SUPPOSEE ISOTROPE
      NOFOTI = NOFOTIEL()
C     NOFOTI>0 SI CETTE FONCTION EXISTE
      IF( NOFOTI .GT. 0 ) THEN
C        LA FONCTION 'TAILLE_IDEALE(X,Y,Z)' EXISTE
         IF( NBARLI .LT. 2 ) THEN
            NBARLI = 3
         ENDIF
      ENDIF
C
C     DONNEES CORRECTES?
      IF( NBARLI .LT. 3 ) THEN
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:8),'(I8)') NBARLI
         KERR(1) = 'NOMBRE INCORRECT <2 D''ARETES '//KERR(MXLGER)(1:8)
         CALL LEREUR
         IERR = 2
         GOTO 9999
      ENDIF

      RAYONC = 0
      NS1CER = 0
      NUPLCT = 0

C     LE NUMERO DU TYPE DU CERCLE
      NUTYCI = LADEFI( WUTYCI )

      IF( NUTYCI .EQ. 1 ) THEN
C        CERCLE DEFINI PAR 3 POINTS
         NUPT1 = LADEFI( WUPT1C )
         NUPT2 = LADEFI( WUPT2C )
         NUPT3 = LADEFI( WUPT3C )
      ELSE IF( NUTYCI .EQ. 2 ) THEN
C        CERCLE DEFINI PAR LE CENTRE,UN POINT DU CERCLE,UN POINT DU PLAN
         NUPT1 = LADEFI( WUPTCE )
         NUPT2 = LADEFI( WUPTSC )
         NUPT3 = LADEFI( WUPDPL )
      ELSE IF( NUTYCI .EQ. 3 ) THEN
C        CERCLE DEFINI PAR LE CENTRE,UN RAYON et un PLAN X ou Y ou Z=Cte
         NUPT1  = LADEFI( WUPTCE )
         RAYONC = RADEFI( WAYDCI )
         NUPLCT = LADEFI( WUPLCT )
         GOTO 34
      ELSE
         NBLGRC(NRERR) = 1
         WRITE(KERR(MXLGER)(1:2),'(I2)') NUTYCI
         KERR(1) =  'TYPE DE CERCLE INCORRECT '//KERR(MXLGER)(1:2)
         CALL LEREUR
         IERR = 3
         GOTO 9999
      ENDIF

C     COORDONNEES ET DISTANCE ENTRE CES 3 POINTS
      ANGDER = 0
      DO 20 I=1,3
C        LES 3 COORDONNEES DU POINT I
         IF( NUPT(I) .LE. 0 ) THEN
            NBLGRC(NRERR) = 2
            WRITE(KERR(MXLGER)(1:5),'(I5)') I
            KERR(1) = 'LE POINT'//KERR(MXLGER)(1:5)
            KERR(2) = 'DE DEFINITION DU CERCLE EST INCONNU'
            CALL LEREUR
            IERR = 4
            GOTO 9999
         ENDIF
C
         CALL LXNLOU( NTPOIN, NUPT(I), NTPOI, MN )
         CALL LXTSOU( NTPOI , 'XYZSOMMET',  NTSOM, MNSOM )
         MN = MNSOM + WYZSOM
         XYZ( 1, I ) = RMCN( MN )
         XYZ( 2, I ) = RMCN( MN + 1 )
         XYZ( 3, I ) = RMCN( MN + 2 )
 20   CONTINUE

C     TEST DE LA COHERENCE DES DONNEES
      DO 30 I=1,3
         IF( I .NE. 3 ) THEN
            I1 = I + 1
         ELSE
            I1 = 1
         ENDIF
         A = DIST2P( XYZ(1,I), XYZ(1,I1) )
         IF( A .LE. 0. ) THEN
            NBLGRC(NRERR) = 1
            KERR(1) =  '2 POINTS CONFONDUS'
            CALL LEREUR
            IERR = 4
            GOTO 9999
         ENDIF
 30   CONTINUE

      IF( NUTYCI .EQ. 1 ) THEN

C        CERCLE DEFINI PAR 3 POINTS SUR LE CERCLE
C        ========================================
C        LE 1-ER POINT SUR LE CERCLE
         NS1CER = 1
C        PASSAGE AUX COORDONNEES 2D DES 3 POINTS DANS LEUR PLAN
         CALL DF3D2D( XYZ(1,1), XYZ(1,2), XYZ(1,3), D2D3, IERR )
         IF( IERR .GT. 0 ) THEN
            NBLGRC(NRERR) = 1
            KERR(1) =  '3 POINTS DU CERCLE ALIGNES'
            CALL LEREUR
            IERR = 5
            GOTO 9999
         ENDIF
C
C        LE 1-ER POINT SERT D'ORIGINE
         XY(1,1) = 0.
         XY(2,1) = 0.
         DO 32 I=2,3
C           LES COORDONNEES 2D DU POINT I
            CALL CH3D2D( XYZ(1,1), D2D3, XYZ(1,I), XY(1,I) )
 32      CONTINUE

C        LES COORDONNEES 2D DU CENTRE DU CERCLE ET CARRE DU RAYON
         CALL CENCER( XY(1,1), XY(1,2), XY(1,3), CENRAY )
         RAYONC = SQRT( CENRAY( 3 ) )

C        LES COORDONNEES DU CENTRE EN 3D
         CALL CH2D3D( XYZ(1,1), D2D3, CENRAY, CENTRE )
C
      ELSE IF( NUTYCI .EQ. 2 ) THEN

C        CERCLE DEFINI PAR SON CENTRE, UN POINT DU CERCLE
C                          ET  UN POINT DU PLAN
C        =================================================
         RAYONC    = DIST2P( XYZ(1,1), XYZ(1,2) )
         CENTRE(1) = XYZ(1,1)
         CENTRE(2) = XYZ(2,1)
         CENTRE(3) = XYZ(3,1)
         NS1CER    = 2

      ENDIF


 34   IF( NUTYCI .EQ. 3 ) THEN

C        CERCLE DEFINI PAR SON CENTRE, SON RAYON
C                          ET LE PLAN X ou Y ou Z=Cte
C        ============================================
         IF( RAYONC .LE. 0. ) THEN
            NBLGRC(NRERR) = 1
            KERR(1) = 'RAYON DU CERCLE <=0'
            CALL LEREUR
            IERR = 4
            GOTO 9999
         ENDIF

         CALL LXNLOU( NTPOIN, NUPT1, NTPOI, MN )
         CALL LXTSOU( NTPOI, 'XYZSOMMET',  NTSOM, MNSOM )
         MN = MNSOM + WYZSOM
         XYZ( 1, 1 ) = RMCN( MN )
         XYZ( 2, 1 ) = RMCN( MN + 1 )
         XYZ( 3, 1 ) = RMCN( MN + 2 )

         IF( NUPLCT .EQ. 1 ) THEN
C           cercle dans le plan YZ avec X=centre(1)=Cte
            XYZ( 1, 2 ) = XYZ( 1, 1 )
            XYZ( 2, 2 ) = XYZ( 2, 1 ) + RAYONC
            XYZ( 3, 2 ) = XYZ( 3, 1 )

            XYZ( 1, 3 ) = XYZ( 1, 1 )
            XYZ( 2, 3 ) = XYZ( 2, 1 )
            XYZ( 3, 3 ) = XYZ( 3, 1 ) + RAYONC

         ELSE IF( NUPLCT .EQ. 2 ) THEN
C           cercle dans le plan ZX avec Y=centre(2)=Cte
            XYZ( 1, 2 ) = XYZ( 1, 1 )
            XYZ( 2, 2 ) = XYZ( 2, 1 )
            XYZ( 3, 2 ) = XYZ( 3, 1 ) + RAYONC

            XYZ( 1, 3 ) = XYZ( 1, 1 ) + RAYONC
            XYZ( 2, 3 ) = XYZ( 2, 1 )
            XYZ( 3, 3 ) = XYZ( 3, 1 )

         ELSE IF( NUPLCT .EQ. 3 ) THEN
C           cercle dans le plan XY avec Z=centre(3)=Cte
            XYZ( 1, 2 ) = XYZ( 1, 1 ) + RAYONC
            XYZ( 2, 2 ) = XYZ( 2, 1 )
            XYZ( 3, 2 ) = XYZ( 3, 1 )

            XYZ( 1, 3 ) = XYZ( 1, 1 )
            XYZ( 2, 3 ) = XYZ( 2, 1 ) + RAYONC
            XYZ( 3, 3 ) = XYZ( 3, 1 )

         ENDIF

         CENTRE(1) = XYZ(1,1)
         CENTRE(2) = XYZ(2,1)
         CENTRE(3) = XYZ(3,1)
         NS1CER    = 2

      ENDIF
C
C     PASSAGE AUX COORDONNEES CENTRE P2, CENTRE P3
 35   CALL DF3D2D( CENTRE, XYZ(1,NS1CER), XYZ(1,3), D2D3, IERR )
      IF( IERR .GT. 0 ) THEN
C        LE CENTRE ET P_NS1CER P3 SONT ALIGNES
         IF( NS1CER .EQ. 1 ) THEN
            NS1CER = 2
            GOTO 35
         ENDIF
         NBLGRC(NRERR) = 1
         KERR(1) = ' CENTRE ET LES 2 POINTS ALIGNES'
         CALL LEREUR
         IERR = 5
         GOTO 9999
      ENDIF
C
      PI2    = PI + PI
      DEUPIR = PI2 * RAYONC
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
         NBARLI  = 0
         HHH     = DEUPIR
         ANGLE   = 0.0
         XY(1,1) = RAYONC
         XY(2,1) = 0.0
C        RETOUR AUX COORDONNEES 3D
         CALL CH2D3D( CENTRE, D2D3, XY(1,1),  XYZ(1,1) )
         XYZD(1) = XYZ(1,1)
         XYZD(2) = XYZ(2,1)
         XYZD(3) = XYZ(3,1)
C
C        LA 'TAILLE_IDEALE' AU SOMMET XYZD
 38      CALL FONVAL( NOFOTI, 3, XYZD,  NCODEV, DTAILL )
         IF( NCODEV .EQ. 0 ) THEN
C
C           NCODEV  : 0 DTAILL N'EST PAS INITIALISEE EN SORTIE
C                     1 DTAILL   EST     INITIALISEE EN SORTIE
            NBLGRC(NRERR) = 2
            KERR(1) = 'FONCTION TAILLE_IDEALE INCALCULABLE AU POINT'
            KERR(2) = ' '
            WRITE(KERR(2)(1:15), '(E15.7)') XYZD(1)
            WRITE(KERR(2)(18:32),'(E15.7)') XYZD(2)
            WRITE(KERR(2)(35:49),'(E15.7)') XYZD(3)
            CALL LEREUR
            IERR = 1
            GOTO 9999
C
         ELSE
C           TAILLE AUTOUR DU SOMMET EST INITIALISEE
C           ESSAI DE CREER UN ARC DE LONGUEUR CETTE TAILLE
            DTAILL = ABS( DTAILL )
            IF( DTAILL .LT. HHH*0.65 ) THEN
C              CREATION D'UN POINT INTERMEDIAIRE => UNE ARETE DE PLUS
               IF( NBARLI .GE. MXTAIL ) THEN
C                 SATURATION DU TABLEAU TAIL => IL EST AUGMENTE
                  CALL TNMCAU( 'REEL', MXTAIL, MXTAIL+512,
     %                         NBARLI, MNTAIL )
                  MXTAIL = MXTAIL + 512
               ENDIF
               NBARLI = NBARLI + 1
               ANGLE0 = ANGLE
               IPAS   = 0
C
C              LES COORDONNEES 2D DU POINT FINAL DE L'ARC
 7             ANGLE   = ANGLE0 + DTAILL / RAYONC
               XY(1,1) = REAL( RAYONC * COS( ANGLE ) )
               XY(2,1) = REAL( RAYONC * SIN( ANGLE ) )
C              RETOUR AUX COORDONNEES 3D
               CALL CH2D3D( CENTRE, D2D3, XY(1,1), XYZ(1,3) )
               IF( IPAS .EQ. 0 ) THEN
C
C                 LA TAILLE IDEALE EN LA SECONDE EXTREMITE DE L'ARC
                  XYZD(1) = XYZ(1,3)
                  XYZD(2) = XYZ(2,3)
                  XYZD(3) = XYZ(3,3)
                  CALL FONVAL( NOFOTI, 3, XYZD,  NCODEV, DTAIL2 )
                  DTAIL2 = ABS( DTAIL2 )
C
C                 LA TAILLE SOUHAITEE AU MILIEU DE L'ARC
                  A      = ANGLE0 + DTAILL / RAYONC / 2
                  XY1(1) = REAL( RAYONC * COS( A ) )
                  XY1(2) = REAL( RAYONC * SIN( A ) )
C                 RETOUR AUX COORDONNEES 3D
                  CALL CH2D3D( CENTRE, D2D3, XY1, XYZ1 )
                  XYZD(1) = XYZ1(1)
                  XYZD(2) = XYZ1(2)
                  XYZD(3) = XYZ1(3)
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
C              LA TAILLE DE L'ARETE NBARLI
               RMCN(MNTAIL-1+NBARLI) = REAL( DTAILL )
C              LA LONGUEUR RESTANTE DU CERCLE A MAILLER
               HHH = HHH - DTAILL
C              LE SOMMET SUR LE CERCLE
               XYZD(1) = XYZ(1,3)
               XYZD(2) = XYZ(2,3)
               XYZD(3) = XYZ(3,3)
               GOTO 38
            ENDIF
         ENDIF
C
C        L'ANGLE DE LA DERNIERE ARETE
         ANGDER = HHH / RAYONC
C        LA TAILLE DE LA DERNIERE ARETE
         RMCN(MNTAIL+NBARLI) = REAL( HHH )
         NBARLI = NBARLI + 1
C
         IF( NBARLI .LE. 1 ) THEN
            NBLGRC(NRERR) = 1
            KERR(1) =  'NOMBRE INCORRECT D''ARETES'
            CALL LEREUR
            IERR = 2
            GOTO 9999
         ENDIF
C
C        2 TG PAR SOMMET CAR ANGLES DIFFERENTS
         NBTGS = 2 * NBARLI
C
      ELSE
C        1 TG PAR SOMMET CAR ANGLES EGAUX
         NBTGS = NBARLI
      ENDIF
C
C     CONSTRUCTION DU TABLEAU 'XYZSOMMET'
C     -----------------------------------
      CALL LXTNDC( NTLXLI, 'XYZSOMMET', 'ENTIER',
     %             WYZSOM + 3*NBARLI + 3*NBTGS )
      CALL LXTSOU( NTLXLI, 'XYZSOMMET',  NTSOLI , MNSOLI )
C
C     LE NOMBRE DE SOMMETS
      MCN( MNSOLI + WNBSOM ) = NBARLI
C     LE NOMBRE DE TANGENTES
      MCN( MNSOLI + WNBTGS ) = NBTGS
C
C     ADRESSE DU DEBUT DES COORDONNEES DU 1-ER SOMMET DE LA LIGNE
      MNS  = MNSOLI + WYZSOM
C     ADRESSE DU DEBUT DES 3 COMPOSANTES DE LA 1-ERE TANGENTE
      MNTG = MNS + 3 * NBARLI
C
C     LES COORDONNEES DES NBARLI POINTS ET TANGENTES DU CERCLE
C     ========================================================
      IF( NOFOTI .LE. 0 ) THEN
C
C        PAS DE FONCTION TAILLE_IDEALE
C        *****************************
C        ANGLE CONSTANT ENTRE LES SOMMETS
         A     = PI2 / NBARLI
         ANGLE = 0.
         DO 40 I=1,NBARLI
C
C           LES COORDONNEES 2D DU POINT I
            XY(1,1) = REAL( RAYONC * COS( ANGLE ) )
            XY(2,1) = REAL( RAYONC * SIN( ANGLE ) )
C           RETOUR AUX COORDONNEES 3D
            CALL CH2D3D( CENTRE, D2D3, XY(1,1), XYZ(1,1) )
C           RANGEMENT DES 3 COORDONNEES
            RMCN( MNS     ) = XYZ( 1, 1 )
            RMCN( MNS + 1 ) = XYZ( 2, 1 )
            RMCN( MNS + 2 ) = XYZ( 3, 1 )
            MNS = MNS + 3
C
C           LA TANGENTE AU POINT FINAL DE L'ARETE I
C           L'INTERVALLE DU PARAMETRE T EST [0,1] POUR TOUT LE CERCLE
C           LE PARAMETRE T DE L'ARETE I EST DANS [(i-1)/NBARLI, i/NBARLI]
C           LE POINT AU BOUT DE LA TANGENTE PLACEE AU POINT INITIAL DE L'ARC
            XY(1,2) = REAL( XY(1,1) - A * RAYONC * SIN(ANGLE) )
            XY(2,2) = REAL( XY(2,1) + A * RAYONC * COS(ANGLE) )
C           RETOUR AUX COORDONNEES 3D (TRANSFORMATION AFFINE)
            CALL CH2D3D( CENTRE, D2D3, XY(1,2), XYZ(1,2) )
C           LES COMPOSANTES DE LA TG = TRANSFORMEE R2 R3(P1+TG)-TRANSFORMEE R2 R
            RMCN( MNTG     ) = XYZ(1,2) - XYZ(1,1)
            RMCN( MNTG + 1 ) = XYZ(2,2) - XYZ(2,1)
            RMCN( MNTG + 2 ) = XYZ(3,2) - XYZ(3,1)
            MNTG  = MNTG + 3
C
C           L'ANGLE OX POINT I
            ANGLE = ANGLE + A
 40      CONTINUE
C
      ELSE
C
C        LA FONCTION 'TAILLE_IDEALE(X,Y,Z)' EXISTE
C        *****************************************
         ANGLE = 0.0
         DO 45 I=1,NBARLI
C
C           LES COORDONNEES 2D DU PREMIER SOMMET DE L'ARETE I
            XY(1,1) = REAL( RAYONC * COS(ANGLE) )
            XY(2,1) = REAL( RAYONC * SIN(ANGLE) )
C           RETOUR AUX COORDONNEES 3D
            CALL CH2D3D( CENTRE, D2D3, XY(1,1),  XYZ(1,1) )
C           RANGEMENT DES 3 COORDONNEES
            RMCN( MNS     ) = XYZ( 1, 1 )
            RMCN( MNS + 1 ) = XYZ( 2, 1 )
            RMCN( MNS + 2 ) = XYZ( 3, 1 )
            MNS = MNS + 3
C
C           LA TAILLE DE L'ARETE I
            DTAILL = RMCN(MNTAIL-1+I)
C
C           L'ANGLE DE L'ARC = NOUVELLE ARETE
            A = DTAILL / RAYONC
C
C           LE POINT AU BOUT DE LA TANGENTE PLACEE AU POINT INITIAL DE L'ARC
            XY(1,2) = REAL( XY(1,1) - DTAILL * SIN(ANGLE) )
            XY(2,2) = REAL( XY(2,1) + DTAILL * COS(ANGLE) )
C           RETOUR AUX COORDONNEES 3D (TRANSFORMATION AFFINE)
            CALL CH2D3D( CENTRE, D2D3, XY(1,2), XYZ(1,2) )
C           LES COMPOSANTES DE LA TG = TRANSFORMEE R2 R3(P1+TG)-TRANSFORMEE R2 R
            TG(1) = XYZ(1,2) - XYZ(1,1)
            TG(2) = XYZ(2,2) - XYZ(2,1)
            TG(3) = XYZ(3,2) - XYZ(3,1)
C
C           LA TANGENTE DE L'ARC PRECEDENT
            IF( I .NE. 1 ) THEN
C              SOMMET QUELCONQUE
               MN   = MNTG
               MNTG = MNTG+3
            ELSE
C              SOMMET 1 DE TG 1 ET NBTGS
               MN   = MNSOLI + WYZSOM + 3*NBARLI + 3*NBTGS - 3
            ENDIF
C
C           LA TG DE L'ARC PRECEDENT
            RMCN( MN     ) = REAL( -TG(1) * ANGDER / A )
            RMCN( MN + 1 ) = REAL( -TG(2) * ANGDER / A )
            RMCN( MN + 2 ) = REAL( -TG(3) * ANGDER / A )
C
C           LA PREMIERE TANGENTE DE L'ARC
            RMCN( MNTG     ) = TG(1)
            RMCN( MNTG + 1 ) = TG(2)
            RMCN( MNTG + 2 ) = TG(3)
            MNTG = MNTG + 3
C
C           L'ANGLE FINAL DU SOMMET DE L'ARC AVEC OX EN 2D
            ANGLE  = ANGLE + A
            ANGDER = A
 45      CONTINUE
C
C        DESTRUCTION DU TABLEAU DES TAILLES
         CALL TNMCDS( 'REEL', MXTAIL, MNTAIL )
      ENDIF
C
C     AJOUT DE LA DATE
      CALL ECDATE( MCN(MNSOLI) )
C
C     AJOUT DU NUMERO DU TABLEAU DESCRIPTEUR
      MCN( MNSOLI + WBCOOR ) = 3
      MCN( MNSOLI + MOTVAR(6) ) = NONMTD( '~>>>XYZSOMMET' )
C
C     CONSTRUCTION DU TABLEAU 'NSEF' LIGNE NON STRUCTUREE
C     CAR LE DERNIER SOMMET EST AUSSI LE PREMIER
C     ---------------------------------------------------------
      CALL LXTNDC( NTLXLI, 'NSEF', 'ENTIER',
     %             WUSOEF+6*NBARLI )
      CALL LXTSOU( NTLXLI, 'NSEF',  NTARLI , MNARLI )
C
C     LE TYPE DE L'OBJET : ICI LIGNE
      MCN( MNARLI + WUTYOB ) = 2
C     LE TYPE FERME DE FERMETURE DU MAILLAGE
      MCN( MNARLI + WUTFMA ) = 1
C     LE NOMBRE DE SOMMETS PAR EF =2
      MCN( MNARLI + WBSOEF ) = 2
C     LES 2 TANGENTES DE CHAQUE EF (ARETE) SONT STOCKEES
      MCN( MNARLI + WBTGEF ) = 2
C     LE NOMBRE D'EF (ARETES) DU SEGMENT NON STRUCTURE
      MCN( MNARLI + WBEFOB ) = NBARLI
C     2 TANGENTES PAR ARETE
      MCN( MNARLI + WBEFTG ) = NBARLI
C     TOUTE ARETE A UN POINTEUR
      MCN( MNARLI + WBEFAP ) = NBARLI
C     LE TYPE DU MAILLAGE : ICI SEGMENT NON STRUCTURE
      MCN( MNARLI + WUTYMA ) = 0
C
C     LE NUMERO DES 2 SOMMETS DES NBARLI ARETES
      MN = MNARLI + WUSOEF
      DO 58 I=1,NBARLI
         MCN(MN  ) = I
         MCN(MN+1) = I + 1
         MN = MN + 2
 58   CONTINUE
C     LE DERNIER SOMMET DU CERCLE EST LE PREMIER
      MCN(MN-1) = 1
C
C     LE POINTEUR SUR CHAQUE EF A TG
      DO 60 I=1,NBARLI
         MCN(MN) = I
         MN = MN + 1
 60   CONTINUE
C
C     LE CODE CERCLE 1 DE CHAQUE ARETE
      DO 70 I=1,NBARLI
         MCN(MN) = 1
         MN = MN + 1
 70   CONTINUE
C
C     LE NUMERO DES 2 TANGENTES DES NBARLI ARETES
      IF( NOFOTI .LE. 0 ) THEN
C        ANGLE CONSTANT
         DO 80 I=1,NBARLI
            MCN(MN  ) = I
            MCN(MN+1) = -(I+1)
            MN = MN + 2
 80      CONTINUE
C        LA DERNIERE TANGENTE DU CERCLE EST LA PREMIERE
         MCN(MN-1) = -1
      ELSE
C        TAILLE_IDEALE PRESENTE
         DO 90 I=1,NBTGS
            MCN(MN) = I
            MN = MN + 1
 90      CONTINUE
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
