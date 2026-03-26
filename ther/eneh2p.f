      SUBROUTINE ENEH2P
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : TRACER LES n COURBES Valeur Propre j = ENERGIEj(DistanceNeutrons)
C ----- D'UNE MOLECULE ION H2+ CALCULEES PAR SSPACE ou ITEINV
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Texas A & M University           Juillet 2005
C23456---------------------------------------------------------------012
      include"./incl/xvfontes.inc"
      include"./incl/trvari.inc"
      COMMON / UNITES / LECTEU , IMPRIM , INTERA , NUNITE(29)
      CHARACTER*80 TITRE
C
C     LES DONNEES CALCULEES PAR SSPACE ou ITEINV
      PARAMETER (NBENERG=11, NBDIST=13)
      REAL DIST(1:NBDIST)
      REAL ENERGY(1:NBENERG,1:NBDIST)
C
C     R=2*a
      DATA DIST  / 0.5, 1.0, 1.5, 1.8, 1.9, 2.0, 2.1, 2.2, 2.5,
     %             3.0, 4.0, 5.0, 8.0 /
      DATA ENERGY /
C     a=0.25
     % 0.2681706,   1.489415,    1.514104,     1.515010,    1.541900,
     % 1.784446,    1.7908931,   1.791726,     1.792254,    1.793537,
     % 1.794461,
C     a=0.50
     %-0.4494045,   0.4375051,   0.5276636,    0.5277818,   0.5787531,
     % 0.7692395,   0.7795283,   0.7806371,    0.7846013,   0.7851790,
     % 0.7853666,
C     a=0.75
     %-0.5797805,   0.04613118,  0.2160768,    0.2161276,   0.2790970,
     % 0.4179051,   0.4395814,   0.4437041,    0.4438240,   0.4519893,
     % 0.4520600,
C     a=0.9
     %-0.5973673,  -0.09344075,  0.1185861,    0.1186291,   0.1852763,
     % 0.3036364,   0.3269028,   0.3312954,    0.3315788,   0.3425681,
     % 0.3430987,
C     a=0.95
     %-0.5991067,  -0.1304783,   0.09388066,   0.09395196,  0.1613674,
     % 0.2731545,   0.2944282,   0.3017263,    0.3019121,   0.3144867,
     % 0.3144948,
C     a=1.0
     %-0.5994680,  -0.1636896,   0.07200654,   0.07207861,  0.1400811,
     % 0.2464110,   0.2676205,   0.2749847,    0.2766325,   0.2887060,
     % 0.2887495,
C     a=1.05
     %-0.5988125,  -0.1936356,   0.05270388,   0.05273167,  0.1211901,
     % 0.2220926,   0.2408125,   0.2506632,    0.2507473,   0.2656731,
     % 0.2657570,
C     a=1.1
     %-0.5973241,  -0.2204958,   0.03545368,   0.03548487,  0.1042922,
     % 0.2007734,   0.2170453,   0.2285543,    0.2298886,   0.2448030,
     % 0.2448608,
C     a=1.25
     %-0.5897233,  -0.2864977,  -0.006200359, -0.006173547, 0.06290934,
     % 0.1461017,   0.1568498,   0.1729460,    0.1739267,   0.1927107,
     % 0.1927815,
C     a=1.5
     %-0.5721911,  -0.3604319,  -0.05262957,  -0.05258019,  0.01548461,
     % 0.07747667,  0.08202978,  0.1047535,    0.1050122,   0.1304376,
     % 0.1339177,
C     a=2.0
     %-0.5370280,  -0.4330526,  -0.1004894,   -0.1003958,  -0.03756450,
     %-0.03241620,  0.007233994, 0.01979835,   0.01989723,  0.05625650,
     % 0.05688631,
C     a=2.5
     %-0.5100554,  -0.4592074,  -0.1210866,   -0.1209990,  -0.1013557,
     %-0.06366080, -0.03313263, -0.02949523,  -0.02938099,  0.01533443,
     %0.01543320,
C     a=4.0
     %-0.4656399,  -0.4600512,  -0.1672910,   -0.1340905,  -0.1338437,
     %-0.09360786, -0.09317230, -0.09292010,  -0.08274873, -0.03719096,
     %-0.03552936 /
C
C     LA FENETRE DE TRACE
CCC      CALL xvinitierps( 1 )
      CALL EFFACE
      CALL FENETRE( 0.0, 8.5,  -1.0, 2.0 )
      CALL XVEPAISSEUR( 2 )
      CALL CHOIXFONTE( 20 )
C
C     LA SIGNIFICATION DES AXES
      CALL TEXTE2D( NCNOIR,  0.05, 1.85, 'ENERGY' )
      CALL TEXTE2D( NCNOIR,  7.0, -0.9,  'R' )
      CALL CHOIXFONTE( 20 )
C
C     LE TRACE DES AXES 2D
      CALL TRAXE2
C
C     LE TRACE DES ( DISTi, ENERGYj(DISTi) )
      DO 20 J = 1, NBENERG
C        CHOIX DE LA COULEUR DE TRACE
         NC = J
         IF( NC .GT. 9 ) NC = NC + 2
         DO 10 I = 1, NBDIST
           IF( I .GT. 1 ) THEN
              CALL TRAIT2D( NC, DIST(I-1), ENERGY(J,I-1),
     %                          DIST(I  ), ENERGY(J,I  ) )
           ENDIF
           CALL SYMBOLE2D( NC,  DIST(I), ENERGY(J,I), '*' )
 10      CONTINUE
         CALL ENTIER2D( NC, DIST(J), ENERGY(J,J), J )
 20   CONTINUE
C
      CALL CHOIXFONTE( 35 )
      TITRE = 'H2+  Hydrogen Molecule Ion  11 energies'
      L = NUDCNB( TITRE )
      CALL TEXTE2D( NCROUG, 2.0, 1.8, TITRE(1:L) )
C
C     COPIE DE MEMPX DANS FENETRE
      CALL MEMPXFENETRE
C
C     POUR VIDER LE BUFFER DE X11
      CALL XVVOIR
C
CCC      CALL xvsauverps( 'energyH2plus', 12 )
C
C     POUR ATTENDRE UN CLIC SOURIS ET  LIRE LE GRAPHIQUE
      CALL CHOIXFONTE( NPHFCO )
      CALL CLICSO
C
CCCC     LECTURE D'UN ENTIER POUR LIRE LE GRAPHIQUE
CCC      READ( LECTEU, * ) I
C
CCC      CALL EFFACE
      RETURN
      END
