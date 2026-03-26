      SUBROUTINE QUNSAL( NTLXSU, NUCOTE, NBTGS,  NBSOCT, MNSOCT,
     %                   NBARLI, MNXYTG, MNNTGL,
     %                   NTFASU, MNFASU, NTSOFA, MNSOFA, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERER LE MAILLAGE D'UN QUADRANGLE ALGEBRIQUE STRUCTURE
C -----    OU NON STRUCTURE AVEC OU SANS TANGENTES
C
C ENTREES:
C --------
C NTLXSU : NUMERO DU TABLEAU TS DU LEXIQUE DU QUADRANGLE
C NUCOTE : NUCOTE(3) NUMERO DE 1 A 4 DU COTE 3 PARMI LES 4 LIGNES...
C          SI NUCOTE(3)=-4 LE COTE 3 EST LA LIGNE 4 A PARCOURIR
C                          EN SENS INVERSE DE SON RANGEMENT...
C          POUR LE SENS C1:S1S2 C2:S2S3 C3:S3S4 C4:S4S1
C
C                S4----<-------S3
C                |             |
C                |             |
C               \/             /\
C                |             |
C                |             |
C                S1---->-------S2
C
C NBTGS  : =0 SI PAS DE TANGENTES STOCKEES POUR LES 4 LIGNES COTES
C          >0 SINON
C NBSOCT : NOMBRE DE SOMMETS DES 4 COTES DU QUADRANGLE
C MNSOCT : ADRESSE MCN DU TABLEAU XYZSOMMET DES 4 LIGNES
C NBARLI : NOMBRE D'ARETES DE LA LIGNE  (0 SI PROBLEME RENCONTRE)
C MNXYTG : ADRESSE MCN DU TABLEAU DES 3 COMPOSANTES DES TANGENTES
C MNNTGL : ADRESSE MCN DU TABLEAU DES 2 NUMEROS DES TANGENTES DES
C          NBARLI ARETES DE LA LIGNE    (0 SI PROBLEME RENCONTRE)
C          TABLEAU A DETRUIRE EN FIN D'UTILISATION
C
C SORTIES:
C --------
C NTFASU : NUMERO      DU TMS 'NSEF' DES NUMEROS DES SOMMETS DES EF
C MNFASU : ADRESSE MCN DU TMS 'NSEF' DES NUMEROS DES SOMMETS DES EF
C          CF '~/TD/D/A___NSEF'
C NTSOFA : NUMERO      DU TMS 'XYZSOMMET' DE LA SURFACE
C MNSOFA : ADRESSE MCN DU TMS 'XYZSOMMET' DE LA SURFACE
C          CF '~/TD/D/A___XYZSOMMET'
C IERR   : 0 SI PAS D'ERREUR
C        > 0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS         OCTOBRE 1996
C234567..............................................................012
      include"./incl/ntmnlt.inc"
      include"./incl/gsmenu.inc"
      include"./incl/a___xyzsommet.inc"
      include"./incl/a___nsef.inc"
      include"./incl/a_surface__definition.inc"
      PARAMETER        (MXETRI=192)
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL             RMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1))
C     LA LONGUEUR DES ARETES DE LA LIGNE ENVELOPPE DANS R**3
      COMMON / PEL1R3 / PERMT3
C
      INTEGER           LADEFI(0:WULFTR)
      REAL              RADEFI(0:WULFTR)
      EQUIVALENCE      (LADEFI(0),RADEFI(0))
      COMMON /S09S01/   MNCOCU, NBSOCQ(4), NUCOTQ(4), MNSOCQ(4),
     %                  XYZ4ST(3,4), XYZL(3,4), XY(4), ST(2,4), UNMOT
      INTEGER           NUCOTE(4),NBSOCT(4),NBSTCT(4),MNSOCT(4)
      INTEGER           NBARLI(4),MNXYST(4),
     %                  MNXYTG(4),MNNTGL(4)
      REAL              RLONGC(4)
      INTEGER           MNCUCT(4)
      CHARACTER*24      KNOMLI
      DATA              KNOMLI / 'CONTOUR_QUADRANGLE_PLAM ' /
C
C     NECESSAIRE POUR NE PAS FAIRE PLUSIEURS FOIS LA DECLARATION DU TABLEAU
      MNCOCU = 0
      NBCOOR = 3
C
C     FORMATION DU CONTOUR FERME DU QUADRANGLE PLAN DONT LES LONGUEURS
C     DES 4 COTES SONT EGALES A CELLES DES COTES DU QUADRANGLE COURBE
C     ================================================================
      CALL LISTQU( KNOMLI, NUCOTE, NBSOCT, MNSOCT,
     %             XYZ4ST, RLONGC,
     %             NBS,    MNXY,   IERR )
      IF( IERR .NE. 0 ) RETURN
C     MNXY   : ADRESSE MCN DE LA PREMIERE COORDONNEE DU PREMIER SOMMET
C              DE LA LIGNE FORMANT LE CONTOUR DU QUADRANGLE PLAN DONT LA
C              LONGUEUR DES COTES EST CELLE DES COTES DU QUADRANGLE COURBE
C     LE SENS DE PARCOURS DES COTES N'A PAS ETE MODIFIE
C
C     PERIMETRE DU QUADRANGLE DANS R**3
      PERMT3 = RLONGC(1) + RLONGC(2) + RLONGC(3) + RLONGC(4)
C
C     SIMULATION DES DONNEES DU TABLEAU LADEFI POUR LA CONSTRUCTION
C     DE LA TRIANGULATION
C     =============================================================
C     variable NBLFTR nombre de lignes fermees contour de la surface
      LADEFI(WBLFTR) = 1
C     variable NBPTIT nombre de points internes futurs sommets
      LADEFI(WBPTIT) = 0
C     RECHERCHE DU NUMERO DE LA LIGNE CONTOUR
C     tableau  NULFTR(1..NBLFTR)nom des lignes fermees(enveloppe en premier)
      CALL LXNMNO( NTLIGN , KNOMLI , I, MN )
      LADEFI(WULFTR) = I
C
C     LES 2 COORDONNEES DES 4 SOMMETS DU QUADRANGLE PLAN SELON LA
C     LONGUEUR DES 4 COTES DU QUADRANGLE TRANSFINI
C     ===========================================================
      ST(1,1) = 0.0
      ST(2,1) = 0.0
      ST(1,2) = RLONGC(1)
      ST(2,2) = 0.0
      CALL SOQU4L( RLONGC, ST(1,3), ST(2,3), ST(1,4), ST(2,4), IERR )
C
C     TRIANGULATION DU QUADRANGLE PLAN FORME SELON LES LONGUEURS DES COTES
C     ====================================================================
C     LES 3 APPLICATIONS EN JEU
C     Q1  : CARRE => Q   Q1(U,V)=(R,S)
C     T   : CARRE => QUADRANGLE TRANSFINI  T(U,V)=(X,Y,Z)
C     F   : Q     => QUADRANGLE TRANSFINI  F(R,S)=(X,Y,Z)=T(Q1**(-1))
C
C     EXISTENCE OU NON DE LA FONCTION 'TAILLE_IDEALE' DES ARETES AUTOUR DU POINT
C     ICI LA CARTE EST SUPPOSEE ISOTROPE
      NOFOTI = NOFOTIEL()
C     NOFOTI>0 SI CETTE FONCTION EXISTE
C
C     QUELQUES COPIES DANS LE COMMON /S09S01/ PUISQUE LA FONCTION
C     TAILLE_IDEALE EXISTE
      DO 5 I=1,4
         NBSOCQ(I) = NBSOCT(I)
         NUCOTQ(I) = NUCOTE(I)
         MNSOCQ(I) = MNSOCT(I)
 5    CONTINUE
C
C     LA COORDONNEE CURVILIGNE HOMOGENE DE CHAQUE SOMMET DE CHACUN DES
C     4 COTES DU QUADRANGLE PLAN CONSTRUIT SELON LES LONGUEURS
      CALL PACTQT( 4, MNXY, NUCOTE, NBSOCT, RLONGC, MNCOCU )
C
C     variable ARETMX taille max des aretes des triangles equilateraux
      RADEFI(WRETMX)=(RLONGC(1)+RLONGC(2)+RLONGC(3)+RLONGC(4))/NBS*1.3
C
C     EXECUTION DE LA TRIANGULATION DU QUADRANGLE DEFINI PAR SES 4 SOMMETS
      CALL SUEX09( 1, NTLXSU, LADEFI, RADEFI,
     %             NTFASU, MNFASU, NTSOFA, MNSOFA, IERR )
      IF( IERR .NE. 0 ) RETURN
C     LE NOMBRE DE SOMMETS DE LA TRIANGULATION DU QUADRANGLE PLAN
      NBSOM = MCN( MNSOFA + WNBSOM )
C
C     TRANSFORMATION DU QUADRANGLE TRIANGULE PLAN EN LE CARRE DE REFERENCE
C     ====================================================================
      MN = MNSOFA + WYZSOM
      DO 10 I=1, NBSOM
         X = RMCN( MN )
         Y = RMCN( MN + 1 )
C        INVERSION D'UNE TRANSFORMATION Q1**2
         CALL FQ1INV( X, Y, ST, RMCN(MN), RMCN(MN+1), IERR )
         IF( IERR .NE. 0 ) RETURN
         MN = MN + 3
 10   CONTINUE
C
C     TRANSPORT DE CE QUADRANGLE DE REFERENCE TRIANGULE SUR SA SURFACE GAUCHE
C     =======================================================================
      NBS4CT = NBSOCT(1)+NBSOCT(2)+NBSOCT(3)+NBSOCT(4)
      IF( NBTGS .LE. 0 ) THEN
C
C        -------------------------------------------------
C        PAS DE TANGENTES A STOCKER DANS XYZSOMMET ET NSEF
C        -------------------------------------------------
C
         IF( MNCOCU .EQ. 0 ) THEN
C           LE TABLEAU D'ADRESSE MNCOCU N'EXISTE PAS ENCORE => IL EST CONSTRUIT
C           LA COORDONNEE CURVILIGNE HOMOGENE DE CHAQUE SOMMET DE CHACUN DES
C           4 COTES DU QUADRANGLE PLAN CONSTRUIT SELON LES LONGUEURS
C           ===================================================================
            CALL PACTQT( 4, MNXY, NUCOTE, NBSOCT, RLONGC, MNCOCU )
         ENDIF
C
C        TRANSPORT DE CE CARRE UNITE TRIANGULE SUR SA SURFACE GAUCHE
C        ===========================================================
         MN = MNSOFA + WYZSOM
         DO 20 J=1, NBSOM
C           LES COORDONNEES DU POINT DANS LE CARRE UNITE
            U = RMCN( MN )
            V = RMCN( MN + 1 )
C           T : CARRE => QUADRANGLE TRANSFINI  T(U,V)=(X,Y,Z)
C           LES COORDONNEES DU POINT DU QUADRANGLE TRANSFINI DANS R**3
            CALL XYZTRT( U, V, RMCN(MN) )
            MN = MN + 3
 20      CONTINUE
C
C        DESTRUCTION DU TABLEAU
         CALL TNMCDS( 'REEL', NBS4CT, MNCOCU )
C
      ELSE
C
C        ---------------------------------------------
C        EXISTENCE DE TANGENTES DANS XYZSOMMET ET NSEF
C        ---------------------------------------------
C
C        TRAITEMENT DU TMS 'NSEF' DU QUADRANGLE COURBE FINAL
C        ===================================================
C        PROTECTION DU NUMERO DES "4" SOMMETS DE CHAQUE TRIANGLE
C        DU CARRE ET DU QUADRANGLE
         NBEF   = MCN(MNFASU+WBEFOB)
         NBTGEF = 8
         CALL TNMCDC( 'ENTIER', 4*NBEF, MNNSTR )
         CALL TRTATA( MCN(MNFASU+WUSOEF), MCN(MNNSTR), 4*NBEF )
C
C        DESTRUCTION DU TMS 'NSEF' DE LA SURFACE
         CALL LXTSDS( NTLXSU, 'NSEF' )
C
C        CONSTRUCTION DU TABLEAU 'NSEF' DE CETTE SURFACE
         CALL LXTNDC( NTLXSU , 'NSEF' , 'ENTIER' ,
     %                WUSOEF+NBEF*( 6 + NBTGEF ))
         CALL LXTSOU( NTLXSU , 'NSEF' ,  NTFASU , MNFASU )
C        TYPE DE L'OBJET : SURFACE
         MCN ( MNFASU + WUTYOB ) = 3
C        SURFACE NON FERMEE
         MCN ( MNFASU + WUTFMA ) = 0
C        LE NOMBRE DE SOMMETS PAR FACE
         MCN ( MNFASU + WBSOEF ) = 4
C        LE NOMBRE DE TANGENTES STOCKEES PAR EF : SURFACE C1
         MCN ( MNFASU + WBTGEF ) = NBTGEF
C        LE NOMBRE D'EF DE LA SURFACE
         MCN ( MNFASU + WBEFOB ) = NBEF
C        LE NOMBRE D'EF AVEC TANGENTES DE LA SURFACE
         MCN ( MNFASU + WBEFTG ) = NBEF
C        LE NOMBRE D'EF AVEC POINTEUR SUR LES EF  A TG DE LA SURFACE
         MCN ( MNFASU + WBEFAP ) = NBEF
C        NUMERO DU TYPE DU MAILLAGE : NON STRUCTURE
         MCN ( MNFASU + WUTYMA ) = 0
C
C        COPIE DES 4 NUMEROS DES SOMMETS DES TRIANGLES
         CALL TRTATA( MCN(MNNSTR), MCN(MNFASU+WUSOEF), 4*NBEF )
C
C        LE TABLEAU DES POINTEURS ET DES TANGENTES
         MN = MNFASU + WUSOEF + 4 * NBEF - 1
         DO 81 N=1,NBEF
            MN = MN + 1
            MCN(MN) = N
C           CODE GEOMETRIQUE 'INTERPOLATION TRANSFINIE'
            MCN(MN+NBEF) = 16
 81      CONTINUE
         MN = MN + NBEF
         J  = 0
         DO 85 N=1,NBEF
C           LES 6 TANGENTES DU TRIANGLE N
            DO 83 I=1,6
               J  = J + 1
               MCN(MN+I) = J
 83         CONTINUE
C           LES 2 TANGENTES NULLES DU "4-EME" SOMMET DU TRIANGLE
            MCN(MN+7) = 0
            MCN(MN+8) = 0
            MN = MN + 8
 85      CONTINUE
C
C        LA DATE DE CREATION
         CALL ECDATE( MCN(MNFASU) )
C        LE NUMERO DU TABLEAU DESCRIPTEUR
         MCN( MNFASU + MOTVAR(6) ) = NONMTD( '~>>>NSEF' )
C
C        DESTRUCTION DU TABLEAU AUXILIAIRE
         CALL TNMCDS( 'ENTIER', 4*NBEF, MNNSTR )
C
C        TRAITEMENT DU TMS 'XYZSOMMET' DU QUADRANGLE COURBE FINAL
C        ========================================================
C
C        PROTECTION DES 2 COORDONNEES DES SOMMETS DU CARRE UNITE
         CALL TNMCDC( 'REEL', 2*NBSOM, MNSTCA )
         MNL0 = MNSOFA + WYZSOM
         MNL1 = MNSTCA
         DO 86 N=1,NBSOM
            RMCN(MNL1  ) = RMCN(MNL0  )
            RMCN(MNL1+1) = RMCN(MNL0+1)
            MNL1 = MNL1 + 2
            MNL0 = MNL0 + 3
 86      CONTINUE
C
C        ATTENTION: LA LOGIQUE DE RANGEMENT DES SOMMETS FRONTALIERS
C        DANS LE SENS DIRECT CHANGE A NOUVEAU EN REVENANT AU
C        RANGEMENT PAR COTES INITIAL
C        COTE 1 S1->S2  COTE 2 S2->S3  COTE 3 S4->S3 COTE 4 S1->S4
C
C                S4---->-------S3
C                |             |
C                |             |
C                ^             ^
C                |             |
C                |             |
C                S1---->-------S2
C
         NUCOTE(3) = -NUCOTE(3)
         NUCOTE(4) = -NUCOTE(4)
C
C        CALCUL DES COORDONNEES CURVILIGNES HOMOGENES DES LIGNES
C        DES 4 COTES DU QUADRANGLE COURBE
C        -------------------------------------------------------
         CALL CRCUCT( 4, NUCOTE, MNSOCT, RLONGC, MNCUCT )
C
C        MISE A JOUR DES COORDONNEES CURVILIGNES HOMOGENES DES
C        SOMMETS DES COTES DU QUADRANGLE PLAN
C        -----------------------------------------------------
C        LE SOMMET (0,0)
         MNC = MNCUCT(1)
         MN  = MNSTCA
         RMCN(MN  ) = 0
         RMCN(MN+1) = 0
         MN  = MN  + 2
C
C        LE COTE 1
C        LE NUMERO DE LA LIGNE DU COTE 1
         N = ABS( NUCOTE(1) )
C        L'ABSCISSE CURVILIGNE HOMOGENE DU COTE 1
         MNC = MNCUCT(1) + 1
         DO 96 I=2,NBSOCT(N)-1
            RMCN(MN  ) = RMCN(MNC)
            RMCN(MN+1) = 0
            MNC = MNC + 1
            MN  = MN  + 2
 96      CONTINUE
C
C        LE SOMMET (1,0)
         RMCN(MN  ) = 1
         RMCN(MN+1) = 0
         MN  = MN  + 2
C
C        LE COTE 2
C        LE NUMERO DE LA LIGNE DU COTE 2
         N = ABS( NUCOTE(2) )
C        L'ABSCISSE CURVILIGNE HOMOGENE DU COTE 2
         MNC = MNCUCT(2) + 1
         DO 97 I=2,NBSOCT(N)-1
            RMCN(MN  ) = 1
            RMCN(MN+1) = RMCN(MNC)
            MNC = MNC + 1
            MN  = MN  + 2
 97      CONTINUE
C
C        LE SOMMET (1,1)
         RMCN(MN  ) = 1
         RMCN(MN+1) = 1
         MN  = MN  + 2
C
C        LE COTE 3
C        LE NUMERO DE LA LIGNE DU COTE 3
         N = ABS( NUCOTE(3) )
C        L'ABSCISSE CURVILIGNE HOMOGENE DU COTE 3
         MNC = MNCUCT(3) + NBSOCT(N) - 2
         DO 98 I=2,NBSOCT(N)-1
            RMCN(MN  ) = RMCN(MNC)
            RMCN(MN+1) = 1
            MNC = MNC - 1
            MN  = MN  + 2
 98      CONTINUE
C
C        LE SOMMET (0,1)
         RMCN(MN  ) = 0
         RMCN(MN+1) = 1
         MN  = MN + 2
C
C        LE COTE 4
C        LE NUMERO DE LA LIGNE DU COTE 4
         N = ABS( NUCOTE(4) )
C        L'ABSCISSE CURVILIGNE HOMOGENE DU COTE 4
         MNC = MNCUCT(4) + NBSOCT(N) - 2
         DO 99 I=2,NBSOCT(N)-1
            RMCN(MN  ) = 0
            RMCN(MN+1) = RMCN(MNC)
            MNC = MNC - 1
            MN  = MN  + 2
 99      CONTINUE
C
C        RESERVATION DES COORDONNEES DES 2 DERIVEES EN CHAQUE SOMMET
         CALL TNMCDC( 'REEL', 6*NBSOM, MN2DER )
C
C        LE NOMBRE DE TANGENTES STOCKEES
         NBTGS = 6 * NBEF
C
C        DESTRUCTION DU TMS 'XYZSOMMET' DE LA SURFACE
         CALL LXTSDS( NTLXSU, 'XYZSOMMET' )
C
C        CONSTRUCTION DU TABLEAU 'XYZSOMMET' DE CETTE SURFACE
         CALL LXTNDC( NTLXSU, 'XYZSOMMET', 'MOTS',
     %                WYZSOM + 3 * (NBSOM+NBTGS) )
         CALL LXTSOU( NTLXSU, 'XYZSOMMET',  NTSOFA, MNSOFA )
C        LE NOMBRE DE COORDONNEES DES SOMMETS
         MCN( MNSOFA + WBCOOR ) = NBCOOR
C        LE NOMBRE DE SOMMETS
         MCN( MNSOFA + WNBSOM ) = NBSOM
C        LE NOMBRE DE TANGENTES
         MCN( MNSOFA + WNBTGS ) = NBTGS
C
C        L'ADRESSE DES SOMMETS DES 4 LIGNES=COTES DU QUADRANGLE COURBE
         CALL CRTGCT( 4, NUCOTE, MNSOCT, MNXYST, NBARLI, MNNTGL )
C
C        CALCUL DES 3 COORDONNEES DES SOMMETS DU QUADRANGLE COURBE
C        ET DES 3 COMPOSANTES DES 2 DERIVEES /Du ET /Dv
         MNL1 = MNSTCA
         MNL0 = MNSOFA + WYZSOM
         MNC  = MN2DER
         N1   = ABS( NUCOTE(1) )
         N2   = ABS( NUCOTE(2) )
         N3   = ABS( NUCOTE(3) )
         N4   = ABS( NUCOTE(4) )
         NBSTCT(1) = NBSOCT( N1 )
         NBSTCT(2) = NBSOCT( N2 )
         NBSTCT(3) = NBSOCT( N3 )
         NBSTCT(4) = NBSOCT( N4 )
C
         DO 180 N=1,NBSOM
C           LES 2 COORDONNEES DU POINT DANS LE CARRE UNITE PERMETTENT LE
C           CALCUL DES 3 COORDONNEES DU SOMMET DU QUADRANGLE PAR L'
C           INTERPOLATION TRANSFINIE DE COOK GORDON ET HALL
            XM = RMCN( MNL1 )
            YM = RMCN( MNL1 + 1 )
            CALL QU1STG( XM,     YM,   NBSTCT,
     S                   RMCN(MNCUCT(1)), RMCN(MNCUCT(2)),
     S                   RMCN(MNCUCT(3)), RMCN(MNCUCT(4)),
     S                   MNXYST(N1), MNXYST(N2), MNXYST(N3), MNXYST(N4),
     S                   MNXYTG(N1), MNXYTG(N2), MNXYTG(N3), MNXYTG(N4),
     S                   MNNTGL(N1), MNNTGL(N2), MNNTGL(N3), MNNTGL(N4),
     S                   RMCN(MNL0), RMCN(MNC) )
C           PASSAGE AU SOMMET SUIVANT
            MNL1 = MNL1 + 2
            MNL0 = MNL0 + 3
            MNC  = MNC  + 6
 180     CONTINUE
C
C        LA TRANSFORMATION DES TANGENTES CALCULEES EN CHAQUE SOMMET
C        DE LA TRIANGULATION DU CARRE UNITE EN LES 6 TANGENTES SUR
C        CHAQUE TRIANGLE DU MAILLAGE DU CARRE UNITE
C        ----------------------------------------------------------
         CALL QTDETG( MNSTCA, MN2DER, MNSOFA, MNFASU )
C
C        LA DATE DE CREATION
         CALL ECDATE( MCN(MNSOFA) )
C        LE NUMERO DU TABLEAU DESCRIPTEUR
         MCN( MNSOFA + MOTVAR(6) ) = NONMTD( '~>>>XYZSOMMET' )
C
C        DESTRUCTION DES TABLEAUX AUXILIAIRES
C        ------------------------------------
         CALL TNMCDS( 'REEL', NBS4CT,  MNCUCT(1) )
         CALL TNMCDS( 'REEL', 2*NBSOM, MNSTCA )
         CALL TNMCDS( 'REEL', 6*NBSOM, MN2DER )
         DO 190 N=1,4
            CALL TNMCDS( 'REEL', 3*NBSOCT(N), MNXYST(N) )
C           LES EVENTUELS TABLEAUX CREES DANS L'EXECUTION DU SP TGARLI
            IF( MNNTGL(N) .GT. 0 ) THEN
               CALL TNMCDS( 'ENTIER' , 3*NBARLI(N), MNNTGL(N) )
            ENDIF
 190     CONTINUE
      ENDIF
C
C     DESTRUCTION DE LA LIGNE CONTOUR DU QUADRANGLE PLAN
      CALL LXLXDS( NTLIGN, KNOMLI )
      END
