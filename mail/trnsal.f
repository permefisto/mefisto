        SUBROUTINE TRNSAL( NTLXSU, NUCOTE, NBTGS,
     %                     NBSOCT, MNSOCT,
     %                     NBARLI, MNXYTG, MNNTGL,
     %                     NTFASU, MNFASU, NTSOFA, MNSOFA, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    GENERER LE MAILLAGE D'UN TRIANGLE ALGEBRIQUE
C -----    NON STRUCTURE AVEC OU SANS TANGENTES
C          PAR INTERPOLATION TRANSFINIE CF CRAS A. PERRONNET 1/1998
C
C ENTREES:
C --------
C NTLXSU : NUMERO DU TABLEAU TS DU LEXIQUE DU TRIANGLE
C NUCOTE : NUCOTE(3) NUMERO DE 1 A 3 DU COTE 3 PARMI LES 3 LIGNES...
C          SI NUCOTE(3)=-2 LE COTE 3 EST LA LIGNE 2 A PARCOURIR
C                          EN SENS INVERSE DE SON RANGEMENT...
C          POUR LE SENS C1:S1S2 C2:S2S3 C3:S3S1
C
C                S3
C                |  \
C                |    \
C          C3   \/     /\ C2
C                |         \
C                |            \
C                S1---->-------S2
C                      C1
C
C NBTGS  : =0 SI PAS DE TANGENTES STOCKEES POUR LES 3 LIGNES COTES
C          >0 SINON
C NBSOCT : NOMBRE DE SOMMETS DES 3 COTES DU TRIANGLE
C MNSOCT : ADRESSE MCN DU TABLEAU XYZSOMMET DES 3 LIGNES
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
C          CF '~TD/D/A___NSEF'
C NTSOFA : NUMERO      DU TMS 'XYZSOMMET' DE LA SURFACE
C MNSOFA : ADRESSE MCN DU TMS 'XYZSOMMET' DE LA SURFACE
C          CF '~TD/D/A___XYZSOMMET'
C IERR   : 0 SI PAS D'ERREUR
C        > 0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS         OCTOBRE 1997
C234567..............................................................012
      IMPLICIT INTEGER (W)
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
      INTEGER           LADEFI(0:WULFTR)
      REAL              RADEFI(0:WULFTR)
      EQUIVALENCE      (LADEFI(0),RADEFI(0))
C
      COMMON /S09S01/   MNCOCU, NBSOCQ(4), NUCOTQ(4), MNSOCQ(4),
     %                  XYZ3ST(3,4), XYZL(3,4), XY(4), ST(2,4), UNMOT
C     LA LONGUEUR DES ARETES DE LA LIGNE ENVELOPPE DANS R**3
      COMMON / PEL1R3 / PERMT3
C
      INTEGER           NUCOTE(3),NBSOCT(3),NBSTCT(3),MNSOCT(3)
      INTEGER           NBARLI(3),MNXYST(3),
     %                  MNXYTG(3),MNNTGL(3),
     %                  MNCUCT(3)
      REAL              RLONGC(3)
      CHARACTER*24      KNOMLI
      DATA              KNOMLI / 'CONTOUR_TRIANGLE_PLAM   ' /
C
C     NECESSAIRE POUR NE PAS FAIRE PLUSIEURS FOIS LA DECLARATION DU TABLEAU
      MNCOCU = 0
C
C     FORMATION DU CONTOUR FERME DU TRIANGLE PLAN DONT LES LONGUEURS
C     DES 3 COTES SONT EGALES A CELLES DES COTES DU TRIANGLE COURBE
C     ==============================================================
      CALL LISTTR( KNOMLI, NUCOTE, NBSOCT, MNSOCT,
     %             XYZ3ST, RLONGC,
     %             NBS,    MNXY,   IERR )
      IF( IERR .NE. 0 ) RETURN
C
C     PERIMETRE DU TRIANGLE DANS R**3
      PERMT3 = RLONGC(1) + RLONGC(2) + RLONGC(3)
C
C     SIMULATION DES DONNEES DU TABLEAU LADEFI POUR LA CONSTRUCTION
C     DE LA TRIANGULATION
C     =============================================================
C     variable ARETMX 'taille max des aretes des triangles equilateraux'
      RADEFI(WRETMX) = (RLONGC(1) + RLONGC(2) + RLONGC(3) ) / NBS * 1.3
C     variable NBLFTR 'nombre de lignes fermees contour de la surface'
      LADEFI(WBLFTR) = 1
C     variable NBPTIT 'nombre de points internes futurs sommets'
      LADEFI(WBPTIT) = 0
C     RECHERCHE DU NUMERO DE LA LIGNE CONTOUR
C     tableau  NULFTR(1..NBLFTR)'nom des lignes fermees(enveloppe en premier)'
      CALL LXNMNO( NTLIGN , KNOMLI , I, MN )
      LADEFI(WULFTR) = I
C
C     CALCUL DES 2 COORDONNEES X3 Y3  DU 3-EME SOMMET DU TRIANGLE PLAN
C     ================================================================
      CALL SOTR3L( RLONGC, X3, Y3, IERR )
      ST(1,1) = 0.0
      ST(2,1) = 0.0
      ST(1,2) = RLONGC(1)
      ST(2,2) = 0.0
      ST(1,3) = X3
      ST(2,3) = Y3
C
C     TRIANGULATION DU TRIANGLE PLAN FORME SELON LES LONGUEURS DES COTES
C     ==================================================================
C     LES 3 APPLICATIONS EN JEU
C     P1  : TRIANGLE RECTANGLE UNITE => TRIANGLE PLAN       P1(U,V)=(R,S)
C     T   : TRIANGLE RECTANGLE UNITE => TRIANGLE ALGEBRIQUE  T(U,V)=(X,Y,Z)
C     F   : TRIANGLE PLAN  => TRIANGLE ALGEBRIQUE  F(R,S)=(X,Y,Z)=T(P1**(-1))
C
C     QUELQUES COPIES DANS LE COMMON /S09S01/
      DO 5 I=1,3
         NBSOCQ(I) = NBSOCT(I)
         NUCOTQ(I) = NUCOTE(I)
         MNSOCQ(I) = MNSOCT(I)
 5    CONTINUE
C
C     LA COORDONNEE CURVILIGNE HOMOGENE DE CHAQUE SOMMET DE CHACUN DES
C     3 COTES DU TRIANGLE PLAN EST CONSTRUITE SELON LES LONGUEURS DES COTES
      CALL PACTQT( 3, MNXY, NUCOTE, NBSOCT, RLONGC, MNCOCU )
C
C     EXECUTION DE LA TRIANGULATION DU TRIANGLE DEFINI PAR SES 3 SOMMETS
      CALL SUEX09( 6, NTLXSU, LADEFI, RADEFI,
     %             NTFASU, MNFASU, NTSOFA, MNSOFA, IERR )
      IF( IERR .NE. 0 ) RETURN
C     LE NOMBRE DE SOMMETS DE LA TRIANGULATION DU TRIANGLE PLAN (NON UNITE)
      NBSOM = MCN( MNSOFA + WNBSOM )
C
C     TRANSFORMATION DU TRIANGLE TRIANGULE PLAN EN LE TRIANGLE UNITE
C     ==============================================================
      MN = MNSOFA + WYZSOM
      DO 10 I=1, NBSOM
         X = RMCN( MN )
         Y = RMCN( MN + 1 ) / Y3
C        APPLICATION INVERSE DE P1: TRIANGLE RECTANGLE UNITE => TRIANGLE PLAN
C        P1(U,V)=(R,S)  P1**(-1)(R,S)=(U,V)
         RMCN( MN     ) = ( X - Y * X3 ) / RLONGC(1)
         RMCN( MN + 1 ) = Y
         MN = MN + 3
 10   CONTINUE
C
      IF( NBTGS .LE. 0 ) THEN
C
C        -------------------------------------------------
C        PAS DE TANGENTES A STOCKER DANS XYZSOMMET ET NSEF
C        -------------------------------------------------
C
C        TRANSPORT DE CE TRIANGLE UNITE TRIANGULE SUR SA SURFACE GAUCHE
C        ==============================================================
         MN = MNSOFA + WYZSOM
         DO 20 J=1, NBSOM
C           LES COORDONNEES DU POINT DANS LE TRIANGLE UNITE
            U = RMCN( MN )
            V = RMCN( MN + 1 )
C           T: TRIANGLE RECTANGLE UNITE => TRIANGLE ALGEBRIQUE T(U,V)=(X,Y,Z)
C           LES COORDONNEES DU POINT DU TRIANGLE TRANSFINI DANS R**3
            CALL XYZTTT( U, V, RMCN(MN) )
            MN = MN + 3
 20      CONTINUE
C
C        DESTRUCTION DU TABLEAU
         CALL TNMCDS( 'REEL', NBSOCT(1)+NBSOCT(2)+NBSOCT(3), MNCOCU )
C
      ELSE
C
C        ---------------------------------------------
C        EXISTENCE DE TANGENTES DANS XYZSOMMET ET NSEF
C        ---------------------------------------------
C
C        TRAITEMENT DU TMS 'NSEF' DU TRIANGLE COURBE FINAL
C        =================================================
C        PROTECTION DU NUMERO DES "4" SOMMETS DE CHAQUE TRIANGLE
C        DU TRIANGLE UNITE
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
C        TRAITEMENT DU TMS 'XYZSOMMET' DU TRIANGLE COURBE FINAL
C        ======================================================
C
C        PROTECTION DES 2 COORDONNEES DES SOMMETS DU TRIANGLE UNITE
         CALL TNMCDC( 'REEL', 2*NBSOM, MNSTTU )
         MNL0 = MNSOFA + WYZSOM
         MNL1 = MNSTTU
         DO 86 N=1,NBSOM
            RMCN(MNL1  ) = RMCN(MNL0  )
            RMCN(MNL1+1) = RMCN(MNL0+1)
            MNL1 = MNL1 + 2
            MNL0 = MNL0 + 3
 86      CONTINUE
C
C        CALCUL DES COORDONNEES CURVILIGNES HOMOGENES DES LIGNES
C        DES 3 COTES DU TRIANGLE COURBE
C        -------------------------------------------------------
         CALL CRCUCT( 3, NUCOTE, MNSOCT, RLONGC, MNCUCT )
C
C        SENS C1:S1S2 C2:S2S3 C3:S3S1
C
C                S3
C                |  \
C                |    \
C          C3   \/     /\ C2
C                |         \
C                |            \
C                S1---->-------S2
C                      C1
C
C        MISE A JOUR DES COORDONNEES CURVILIGNES HOMOGENES DES
C        SOMMETS DES COTES DU TRIANGLE UNITE MAILL'E
C        -----------------------------------------------------
C        LE SOMMET (0,0)
         MNC = MNCUCT(1)
         MN  = MNSTTU
         RMCN(MN  ) = 0
         RMCN(MN+1) = 0
         MN  = MN  + 2
C
C        LE COTE 1
C        LE NUMERO DE LA LIGNE DU COTE 1 S1->S2
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
C        LE NUMERO DE LA LIGNE DU COTE 2 S2->S3
         N = ABS( NUCOTE(2) )
C        L'ABSCISSE CURVILIGNE HOMOGENE DU COTE 2
         MNC = MNCUCT(2) + 1
         DO 97 I=2,NBSOCT(N)-1
            RMCN(MN  ) = 1 - RMCN(MNC)
            RMCN(MN+1) = RMCN(MNC)
            MNC = MNC + 1
            MN  = MN  + 2
 97      CONTINUE
C
C        LE SOMMET (0,1)
         RMCN(MN  ) = 0
         RMCN(MN+1) = 1
         MN  = MN  + 2
C
C        LE COTE 3
C        LE NUMERO DE LA LIGNE DU COTE 3 S3->S1
         N = ABS( NUCOTE(3) )
C        L'ABSCISSE CURVILIGNE HOMOGENE DU COTE 3
         MNC = MNCUCT(3) + 1
         DO 98 I=2,NBSOCT(N)-1
            RMCN(MN  ) = 0
            RMCN(MN+1) = 1 - RMCN(MNC)
            MNC = MNC + 1
            MN  = MN  + 2
 98      CONTINUE
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
C        LE NOMBRE DE COORDONNEES D'UN SOMMET
         MCN( MNSOFA + WBCOOR) = 3
C        LE NOMBRE DE SOMMETS
         MCN( MNSOFA + WNBSOM) = NBSOM
C        LE NOMBRE DE TANGENTES
         MCN( MNSOFA + WNBTGS) = NBTGS
C
C        L'ADRESSE DES SOMMETS DES 3 LIGNES=COTES DU TRIANGLE COURBE
         CALL CRTGCT( 3, NUCOTE, MNSOCT, MNXYST, NBARLI, MNNTGL )
C
C        CALCUL DES 3 COORDONNEES DES SOMMETS DU TRIANGLE COURBE
C        ET DES 3 COMPOSANTES DES 2 DERIVEES /Du ET /Dv
         MNL1 = MNSTTU
         MNL0 = MNSOFA + WYZSOM
         MNC  = MN2DER
         N1   = ABS( NUCOTE(1) )
         N2   = ABS( NUCOTE(2) )
         N3   = ABS( NUCOTE(3) )
         NBSTCT(1) = NBSOCT( N1 )
         NBSTCT(2) = NBSOCT( N2 )
         NBSTCT(3) = NBSOCT( N3 )
C
         DO 180 N=1,NBSOM
C           LES 2 COORDONNEES DU SOMMET DANS LE TRIANGLE UNITE
C           PERMETTENT LE CALCUL DES 3 COORDONNEES DU SOMMET DU
C           TRIANGLE COURBE PAR L'INTERPOLATION TRANSFINIE A. PERRONNET 1979
            XM = RMCN( MNL1 )
            YM = RMCN( MNL1 + 1 )
            CALL TR1STG( XM, YM, NBSTCT,
     S                   RMCN(MNCUCT(1)), RMCN(MNCUCT(2)),
     S                   RMCN(MNCUCT(3)),
     S                   MNXYST(N1), MNXYST(N2), MNXYST(N3),
     S                   MNXYTG(N1), MNXYTG(N2), MNXYTG(N3),
     S                   MNNTGL(N1), MNNTGL(N2), MNNTGL(N3),
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
         CALL QTDETG( MNSTTU, MN2DER, MNSOFA, MNFASU )
C
C        LA DATE DE CREATION
         CALL ECDATE( MCN(MNSOFA) )
C        LE NUMERO DU TABLEAU DESCRIPTEUR
         MCN( MNSOFA + MOTVAR(6) ) = NONMTD( '~>>>XYZSOMMET' )
C
C        DESTRUCTION DES TABLEAUX AUXILIAIRES
C        ------------------------------------
         CALL TNMCDS( 'REEL', NBSOCT(1)+NBSOCT(2)+NBSOCT(3), MNCUCT(1) )
         CALL TNMCDS( 'REEL', 2*NBSOM, MNSTTU )
         CALL TNMCDS( 'REEL', 6*NBSOM, MN2DER )
         DO 190 N=1,3
            CALL TNMCDS( 'REEL', 3*NBSOCT(N), MNXYST(N) )
C           LES EVENTUELS TABLEAUX CREES DANS L'EXECUTION DU SP TGARLI
            IF( MNNTGL(N) .GT. 0 ) THEN
               CALL TNMCDS( 'ENTIER' , 3*NBARLI(N), MNNTGL(N) )
            ENDIF
 190     CONTINUE
      ENDIF
C
C     DESTRUCTION DE LA LIGNE CONTOUR DU TRIANGLE PLAN
      CALL LXLXDS( NTLIGN, KNOMLI )
      END
