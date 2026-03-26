      SUBROUTINE ENERGYH0
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : TRACER LA COURBE Log( | KEh / PEh -(-1/2) | )
C ----- D'UN ATOME D'HYDROGENE CALCULES PAR SSPACE en FONCTION de h
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC Paris, Veulettes sur mer  Juin 2008
C23456---------------------------------------------------------------012
      include"./incl/xvfontes.inc"
      include"./incl/trvari.inc"
      COMMON / UNITES / LECTEU , IMPRIM , INTERA , NUNITE(29)
      CHARACTER*80 TITRE
C
C     LES DONNEES CALCULEES PAR SSPACE ou ITEINV
      PARAMETER (NBI2=14,NBH=3)
      INTEGER    NBHEXA(NBH), NBVERT(NBH), NBNODE(NBH)
      REAL       H(NBH)
      REAL       VALEUR(4,NBI2,NBH)
C     VALEUR(1,noEV,h) EIGENVALUE
C     VALEUR(2,noEV,h) KINETIC ENERGY
C     VALEUR(3,noEV,h) POTENTIAL ENERGY
C     VALEUR(4,noEV,h) KE/PE
      DATA       NBHEXA / 1024, 3240, 8192 /
      DATA       NBVERT / 1105, 3395, 8449 /
      DATA       NBNODE / 4305, 13377, 33473 /
      DATA       H      / 0.5, 0.3333333, 0.25 /
C
      DATA VALEUR /
     %     -0.5003565,       0.4998528,     -1.000209 ,  -0.4997482,
     %     -0.1248764,       0.1246383,     -0.2495147,  -0.4995228,
     %     -0.1248461,       0.1248282,     -0.2496742,  -0.4999641,
     %     -0.1248459,       0.1247335,     -0.2495794,  -0.4997748,
     %     -0.1248342,       0.1243868,     -0.2492210,  -0.4991023,
     %     -0.5545023E-01,   0.5544995E-01, -0.1109007,  -0.4999965,
     %     -0.5545023E-01,   0.5544941E-01, -0.1108997,  -0.4999962,
     %     -0.5540648E-01,   0.5581231E-01, -0.1112202,  -0.5018179,
     %     -0.5540648E-01,   0.5545690E-01, -0.1108634,  -0.5002273,
     %     -0.5540648E-01,   0.5544236E-01, -0.1108489,  -0.5001617,
     %     -0.5539725E-01,   0.5539760E-01, -0.1107958,  -0.4999973,
     %     -0.5539725E-01,   0.5539677E-01, -0.1107943,  -0.4999968,
     %     -0.5539725E-01,   0.5539654E-01, -0.1107938,  -0.4999966,
     %     -0.5538400E-01,   0.5560043E-01, -0.1109856,  -0.5009696,
     %
     %     -0.5003258,       0.5009562,     -1.001282,   -0.5003148,
     %     -0.1249831,       0.1250565,     -0.2500396,  -0.5001467,
     %     -0.1249591,       0.1249119,     -0.2498710,  -0.4999057,
     %     -0.1249589,       0.1248823,     -0.2498412,  -0.4998466,
     %     -0.1249513,       0.1248393,     -0.2497906,  -0.4997757,
     %     -0.5553090E-01,   0.5553120E-01, -0.1110623,  -0.5000003,
     %     -0.5553090E-01,   0.5553096E-01, -0.1110619,  -0.5000000,
     %     -0.5551843E-01,   0.5551862E-01, -0.1110372,  -0.5000001,
     %     -0.5551843E-01,   0.5551856E-01, -0.1110371,  -0.5000003,
     %     -0.5551843E-01,   0.5551851E-01, -0.1110370,  -0.5000000,
     %     -0.5550931E-01,   0.5564265E-01, -0.1111521,  -0.5005991,
     %     -0.5550931E-01,   0.5557434E-01, -0.1110837,  -0.5002923,
     %     -0.5550931E-01,   0.5559106E-01, -0.1111004,  -0.5003676,
     %     -0.5550422E-01,   0.5554189E-01, -0.1110461,  -0.5001694,
     %
     %     -0.5004582,       0.5011260,     -1.001584,   -0.5003334,
     %     -0.1250165,       0.1250674,     -0.2500840,  -0.5001018,
     %     -0.1249898,       0.1249537,     -0.2499434,  -0.4999278,
     %     -0.1249889,       0.1249715,     -0.2499605,  -0.4999652,
     %     -0.1249341,       0.1247715,     -0.2497056,  -0.4996744,
     %     -0.5555876E-01,   0.5544956E-01, -0.1110084,  -0.4995077,
     %     -0.5554758E-01,   0.5554774E-01, -0.1110955,  -0.4999998,
     %     -0.5554758E-01,   0.5554780E-01, -0.1110955,  -0.5000004,
     %     -0.5554600E-01,   0.5554595E-01, -0.1110921,  -0.4999992,
     %     -0.5554600E-01,   0.5554608E-01, -0.1110921,  -0.5000001,
     %     -0.5554600E-01,   0.5554595E-01, -0.1110920,  -0.4999995,
     %     -0.5554503E-01,   0.5582486E-01, -0.1113711,  -0.5012510,
     %     -0.5554503E-01,   0.5565022E-01, -0.1111954,  -0.5004725,
     %     -0.5554503E-01,   0.5555649E-01, -0.1111015,  -0.5000515 /
C
C     LA FENETRE DE TRACE
      RMIN =  1E20
      RMAX = -1E20
      DO 3 K = 1, NBH
         DO 2 J = 1, NBI2
ccc            V1 = Log( ABS( VALEUR(4,J,K) + 0.5 ) )
            V1 = ABS( VALEUR(4,J,K) + 0.5 )
            print *, 'H(',K,')=', H(K), ' V(',J,')=',V1
            RMIN = MIN( RMIN, V1 )
            RMAX = MAX( RMAX, V1 )
 2       CONTINUE
 3    CONTINUE
      print *,'RMIN=',RMIN,'  RMAX=',RMAX
      RH = ( RMAX - RMIN ) / 20
CCC      CALL xvinitierps( 1 )
      CALL EFFACE
      CALL FENETRE( h(NBH)-0.1, h(1)+0.1,  RMIN-RH, RMAX+RH )
      CALL XVEPAISSEUR( 2 )
      CALL CHOIXFONTE( 20 )
C
C     LA SIGNIFICATION DES AXES
ccc      CALL TEXTE2D( NCNOIR, h(NBH),RMAX+RH/2,'Log(| KEh/PEh -(-1/2) |)')
      CALL TEXTE2D(NCNOIR, h(NBH)-0.08, RMAX+RH/2,'| KEh/PEh -(-1/2) |')
      CALL TEXTE2D(NCNOIR, h(1)  +0.05, RMIN, 'h-mesh' )
      CALL CHOIXFONTE( 20 )
C
C     LE TRACE DES AXES 2D
      CALL TRAXE2
C
C     LE TRACE DES ( hi, Log( | KEh / PEh -(-1/2) | ) ) (hi)
      DO 20 J = 1, NBI2
C        CHOIX DE LA COULEUR DE TRACE
         NC = J
         IF( NC .GT. 9 ) NC = NC + 2
C
ccc         V0 = Log( ABS( VALEUR(4,J,1) + 0.5 ) )
         V0 = ABS( VALEUR(4,J,1) + 0.5 )
         DO 10 K = 2, NBH
ccc            V1 = Log( ABS( VALEUR(4,J,K) + 0.5 ) )
            V1 = ABS( VALEUR(4,J,K) + 0.5 )
            CALL TRAIT2D( NC, H(K-1), V0,  H(K), V1 )
            CALL SYMBOLE2D( NC,  H(K), V1, '*' )
            V0 = V1
 10      CONTINUE
         CALL ENTIER2D( NC, H(NBH), V1, J )
 20   CONTINUE
C
C     LE TITRE DU GRAPHIQUE
      CALL CHOIXFONTE( 35 )
ccc      TITRE = 'H Hydrogen Atom:  Log( | KEh / PEh -(-1/2) | ) (h)'
      TITRE = 'Hydrogen Atom: | KEh / PEh + 1/2 | (h)'
      L = NUDCNB( TITRE )
      CALL TEXTE2D( NCROUG, h(NBH), RMAX, TITRE(1:L) )
C
C     COPIE DE MEMPX DANS FENETRE
      CALL MEMPXFENETRE
C
C     POUR VIDER LE BUFFER DE X11
      CALL XVVOIR
C
CCC      CALL xvsauverps( 'energyH', NBI2 )
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
