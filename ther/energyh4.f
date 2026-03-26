      SUBROUTINE ENERGYH4
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : VERSION 2D INITIALE SANS SUBDIVISIONS INTERNES DES EF
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
      DATA NBTRIA /  25076, 13139,  3641, 1051 /
      DATA NBVERT /  25076,  6649,  1861,  546 /
      DATA NBNODE /  25076, 26436,  7362, 2142 /
      DATA H      /  0.25, 0.5, 1.0, 2.0 /
C
      DATA VALEUR /
     %  1.00000465,   0.499997705,   0.500006911,   0.999981588,
     %  2.00001836,   0.999990819,   1.0000276,   0.999963224,
     %  2.00001836,   0.999990372,   1.00002804,   0.999962328,
     %  3.00004721,   1.50001014,   1.50003716,   0.999981989,
     %  3.0000484,   1.49924263,   1.50080583,   0.998958426,
     %  3.00005746,   1.49917066,   1.50088688,   0.998856534,
     %  4.00015831,   1.99315227,   2.00700599,   0.993097321,
     %  4.00033665,   1.99854203,   2.00179454,   0.998375202,
     %  4.00045824,   1.99687483,   2.00358355,   0.99665164,
     %  4.00158834,   1.99225396,   2.00933448,   0.991499415,
     %  1.00007868,   0.499960885,   0.500117791,   0.999686261,
     %  2.00029373,   0.999863247,   1.00043044,   0.999433052,
     %  2.00029516,   0.99986162,   1.00043343,   0.999428435,
     %  3.00076103,   1.49978651,   1.50097452,   0.99920851,
     %  3.00076985,   1.49963271,   1.50113714,   0.998997809,
     %  3.00100493,   1.49946529,   1.50153969,   0.99861849,
     %  4.00156307,   1.9989167,   2.00264634,   0.998137643,
     %  4.00188828,   1.99855977,   2.00332863,   0.997619531,
     %  4.00218821,   1.99846815,   2.00371993,   0.997378987,
     %  4.00299025,   1.99807184,   2.00491854,   0.996585047,
     %  1.00127208,   0.499416693,   0.501855382,   0.995140654,
     %  2.00487232,   0.997881528,   1.00699077,   0.990953992,
     %  2.00487423,   0.997883654,   1.00699046,   0.990956414,
     %  3.0117681,   1.49508774,   1.51668047,   0.985763167,
     %  3.0117712,   1.49494232,   1.51682876,   0.985570923,
     %  3.01508021,   1.4936835,   1.52139676,   0.981784332,
     %  4.0233717,   1.99481353,   2.02855819,   0.9833652,
     %  4.02358389,   1.99154349,   2.03204038,   0.980070825,
     %  4.03061199,   1.99134567,   2.03926634,   0.976501026,
     %  4.03201818,   1.98467896,   2.0473391,   0.969394354,
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
