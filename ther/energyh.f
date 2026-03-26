      SUBROUTINE ENERGYH
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : VERSION 2D avec TRIANGULATIONS SUBDIVISES 1TRIANGLE=>4SOUS-TRIANGLES
C ----- TRACER LA COURBE Log( | KEh / PEh - 1 | )
C       CALCULES PAR SSPACE en FONCTION de -Log(h)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC Paris & Veulettes 76   Juillet 2008
C23456---------------------------------------------------------------012
      include"./incl/xvfontes.inc"
      include"./incl/trvari.inc"
      COMMON / UNITES / LECTEU , IMPRIM , INTERA , NUNITE(29)
      CHARACTER*120 TITRE
C
C     LES DONNEES CALCULEES PAR SSPACE
      PARAMETER (NBI2=10,NBH=4)
      INTEGER    NBTRIA(NBH), NBVERT(NBH), NBNODE(NBH)
      REAL       H(NBH)
      REAL       VALEUR(4,NBI2,NBH)
C     VALEUR(1,noEV,h) EIGENVALUE
C     VALEUR(2,noEV,h) KINETIC ENERGY
C     VALEUR(3,noEV,h) POTENTIAL ENERGY
C     VALEUR(4,noEV,h) KE/PE
C
C     TABLEAUX AUXILIAIRES
      REAL       PT(2,NBH), LNH0, LNH1, V0, V1
      DOUBLE PRECISION  DROITE(3)
C
C     INITIALISATIONS DES DONNEES
      DATA H      /  0.25,   0.5,   1.0,   2.0 /
      DATA NBVERT /  33789,  8487,  2142,  546 /
      DATA NBNODE / 134841, 33789,  8487, 2142 /
      DATA NBTRIA /  67264, 16816,  4204, 1051 /
C
      DATA VALEUR /
     %  1.00000763,   0.499996196,   0.500011464,   0.999969464,
     %  2.00003028,   0.999984719,   1.00004556,   0.999939164,
     %  2.00003076,   0.999984813,   1.00004588,   0.999938937,
     %  3.00007629,   1.50005506,   1.50002114,   1.00002261,
     %  3.0000875,   1.50013013,   1.49995729,   1.00011523,
     %  3.0000875,   1.50049381,   1.49959375,   1.0006002,
     %  4.00021744,   1.99962089,   2.00059673,   0.999512226,
     %  4.00035429,   1.99946869,   2.0008855,   0.999291912,
     %  4.00150967,   1.99885098,   2.00265853,   0.998098752,
     %  4.00196028,   1.99805519,   2.00390507,   0.997080755,
     %  1.00011992,   0.499941687,   0.500178232,   0.99952708,
     %  2.00046873,   0.999775767,   1.00069303,   0.999083373,
     %  2.00047636,   0.999773655,   1.0007027,   0.999071612,
     %  3.00116134,   1.4994478,   1.50171364,   0.998491169,
     %  3.00116754,   1.49941272,   1.50175474,   0.998440483,
     %  3.00152063,   1.49943779,   1.50208278,   0.998239122,
     %  4.00239372,   1.99871025,   2.00368351,   0.997517939,
     %  4.00248766,   1.99795908,   2.00452844,   0.99672274,
     %  4.00327969,   1.99769005,   2.00558957,   0.996061249,
     %  4.0040946,   1.9970164,   2.0070784,   0.994986745,
     %  1.00176942,   0.499186208,   0.502583213,   0.993240911,
     %  2.00665593,   0.997312805,   1.00934312,   0.98808105,
     %  2.00674438,   0.997407394,   1.00933701,   0.988180741,
     %  3.01611233,   1.4926739,   1.52343846,   0.97980584,
     %  3.01641631,   1.49224487,   1.52417152,   0.979053112,
     %  3.02115655,   1.4897471,   1.53140934,   0.972794836,
     %  4.02944136,   1.99404051,   2.03540099,   0.97967944,
     %  4.03109884,   1.98746459,   2.04363413,   0.972514873,
     %  4.03936911,   1.99535113,   2.04401798,   0.976190599,
     %  4.04021645,   1.99015771,   2.05005864,   0.97078087,
     %  1.03274202,   0.488351341,   0.544390661,   0.897060467,
     %  2.0415864,   0.974435021,   1.06715139,   0.913117886,
     %  2.04370689,   0.973666932,   1.07003989,   0.909935173,
     %  3.22876883,   1.45529762,   1.77347126,   0.82059273,
     %  3.23586559,   1.4592943,   1.77657121,   0.821410532,
     %  3.24639463,   1.55982084,   1.68657388,   0.924845843,
     %  4.23735762,   1.80866795,   2.42868969,   0.744709363,
     %  4.32709026,   1.96747866,   2.35961163,   0.833814614,
     %  4.38504505,   1.79847063,   2.58657437,   0.695309848,
     %  4.39006281,   1.80336577,   2.58669704,   0.6971693 /
C
C     LA FENETRE DE TRACE
      HMIN =  1E20
      HMAX = -1E20
      RMIN =  1E20
      RMAX = -1E20
      DO 3 K = 1, NBH
         HH = -Log( H(K) )
         HMIN = MIN( HMIN, HH )
         HMAX = MAX( HMAX, HH )
         DO 2 J = 1, NBI2
            V1 = Log( ABS( VALEUR(4,J,K) - 1.0 ) )
            RMIN = MIN( RMIN, V1 )
            RMAX = MAX( RMAX, V1 )
 2       CONTINUE
 3    CONTINUE
      print *,'HMIN=',HMIN,' HMAX=',HMAX,' RMIN=',RMIN,'  RMAX=',RMAX
      HH = ( HMAX - HMIN ) / 15
      RH = ( RMAX - RMIN ) / 15
C
C     MISE SUR FICHIER.eps du TRACE
      CALL xvinitierps( 1 )
      CALL EFFACE
      CALL FENETRE( HMIN-HH, HMAX+HH, RMIN-RH, RMAX+2*RH )
C
C     LA SIGNIFICATION DES AXES
      CALL CHOIXFONTE( 20 )
      CALL TEXTE2D( NCNOIR, HMIN, RMAX+RH/3, 'Log(| KEh/PEh - 1 |)')
      CALL TEXTE2D(NCNOIR, HMAX, RMIN, '-Log(h)' )
C
C     LE TRACE DES AXES 2D
      CALL TRAXE2
C
C     LE TRACE DES ( -Log(hi), Log( | KEh / PEh -1 | ) ) (-Log(hi))
      DO 20 J = 1, NBI2
C
C        CHOIX DE LA COULEUR DE TRACE
         NC = J
         IF( NC .GT. 9 ) NC = NC + 2
         print *
C
C        LE POINT A TRACER ( -Log(H(1)), Log( ABS( VALEUR(4,J,1) - 1.0 ) ) )
         LNH0 = -Log(H(1))
         V0  =  Log( ABS( VALEUR(4,J,1) - 1.0 ) )
C        LE NO DE LA VALEUR PROPRE A GAUCHE
         CALL ENTIER2D(  NC, HMAX+HH/7, V0, J )
         CALL SYMBOLE2D( NC, LNH0, V0, '*' )
C
C        DROITE PAR MOINDRES CARRES
         PT(1,1) = LNH0
         PT(2,1) = V0
C
C        BOUCLE SUR LA TAILLE DES ARETES H
         DO 10 K = 2, NBH
C
C           LE POINT A TRACER
            LNH1 = -Log(H(K))
            V1  =  Log( ABS( VALEUR(4,J,K) - 1.0 ) )
            CALL XVEPAISSEUR( 0 )
            CALL TRAIT2D(   NC, LNH0, V0,  LNH1, V1 )
            CALL SYMBOLE2D( NC, LNH1, V1, '*' )
C
C           DROITE PAR MOINDRES CARRES
            PT(1,K) = LNH1
            PT(2,K) = V1
C
            print *,'SLOPE j=',j,' k=',K,'=',
     %              (V1-V0)/(-Log(H(K))+Log(H(K-1)))
C
            LNH0 = LNH1
            V0  = V1
 10      CONTINUE
C
C        LE NO DE LA VALEUR PROPRE A DROITE
         CALL ENTIER2D( NC, HMIN-HH/2, V1, J )
C
C        DROITE PAR MOINDRES CARRES
         CALL DRDIMI( NBH, PT, DROITE, IERR )
C
         print *
         print *,'EIGV',J,': Least square LINE(',J,')=',DROITE
         print *,'Y=',-DROITE(1)/DROITE(2),' * X +',-DROITE(3)/DROITE(2)
         print *,'Least square LINE SLOPE',J,'=',-DROITE(1)/DROITE(2)
C
C        TRACE DE LA DROITE
         LNH0 = -Log(H(1))
         V0  = REAL( - ( DROITE(1) * LNH0 + DROITE(3) ) / DROITE(2) )
         LNH1 = -Log(H(NBH))
         V1  = REAL( - ( DROITE(1) * LNH1 + DROITE(3) ) / DROITE(2) )
         CALL XVEPAISSEUR( 2 )
         CALL TRAIT2D( NC, LNH0, V0,  LNH1, V1 )
C
C        LE NO DE LA DROITE A GAUCHE
         CALL ENTIER2D(  NC, LNH0, V0, J )
C        LE NO DE LA DROITE A DROITE
         CALL ENTIER2D( NC, LNH1, V1, J )
C
 20   CONTINUE
C
C     LE TITRE DU GRAPHIQUE
      CALL CHOIXFONTE( 25 )
      TITRE = '(-1/2)[Laplacian - (x^2 + y^2)] psi(x,y)= E*psi(x,y): Log
     %( | KEh/PEh - 1 | ) (-Log(h))'
      L = NUDCNB( TITRE )
      CALL TEXTE2D( NCROUG, HMIN+HH, RMAX+RH, TITRE(1:L) )
C
C     COPIE DE MEMPX DANS FENETRE
      CALL MEMPXFENETRE
C
C     POUR VIDER LE BUFFER DE X11
      CALL XVVOIR
C
C     MISE SUR FICHIER energyh.eps du TRACE
C     ATTENTION PASSAGE PAR VARIABLE OBLIGATOIRE
      NBC = 8
      CALL xvsauverps( 'energy2d', NBC )
      print *, 'NBC=',NBC
C
C     POUR ATTENDRE UN CLIC SOURIS ET  LIRE LE GRAPHIQUE
      CALL CHOIXFONTE( NPHFCO )
      CALL CLICSO
C
      CALL XVEPAISSEUR( 1 )
      RETURN
      END
