      SUBROUTINE ERR6VP
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : TRACER LES ERREURS sur les VALEURS PROPRES sur [0,1]**6
C ----- en fonction des maillages nM
C       ERREUR(1..7,nM) |VPie-VPic| / VPie
C       avec VPie i-eme Valeur Propre exacte
C       avec VPic i-eme Valeur Propre calculee
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL Paris & St Pierre du Perray Fevrier 2009
C23456---------------------------------------------------------------012
      include"./incl/xvfontes.inc"
      include"./incl/trvari.inc"
      COMMON / UNITES / LECTEU , IMPRIM , INTERA , NUNITE(29)
      CHARACTER*120     TITRE
      DOUBLE PRECISION  PI2
C
C     LES DONNEES CALCULEES PAR THERMICER
      PARAMETER (NBM=3)
      INTEGER    NB6CUB(NBM), NBNODE(NBM)
      REAL       VALPRO(7,NBM)
      REAL       VALPRE(7)
      REAL       ERRELA(7,NBM)
C
C     INITIALISATIONS DES DONNEES
      DATA NBNODE/ 4096, 15625, 46656 /
      DATA NB6CUB/  729,  4096, 15625 /
C
      DATA VALPRO/ 59.4, 77.4, 77.4, 77.4, 95.4, 95.4, 95.4,
     %    59.2777, 81.9051, 81.9051, 81.9051, 96.891, 96.8912, 96.925,
     %    59.2426, 84.2426, 84.2426, 84.2426, 93.9324, 93.9345, 93.948/
C
C     VALEURS PROPRES EXACTES
      PI2 = (ATAN(1D0) * 4D0) ** 2
      VALPRE(1) = REAL( 6D0 * PI2 )
      VALPRE(2) = REAL( 9D0 * PI2 )
      VALPRE(3) = REAL( 9D0 * PI2 )
      VALPRE(4) = REAL( 9D0 * PI2 )
      VALPRE(5) = REAL( 9D0 * PI2 )
      VALPRE(6) = REAL( 9D0 * PI2 )
      VALPRE(7) = REAL( 9D0 * PI2 )
C
C     CALCUL DES ERREURS RELATIVES en %
      DO 2 M = 1, NBM
         DO 1 I=1,7
            ERRELA(I,M)=100 * ABS( VALPRE(I) - VALPRO(I,M) ) / VALPRE(I)
 1       CONTINUE
 2    CONTINUE
C
C     AFFICHAGE DES ERREURS RELATIVES
      DO 7 I=1, 7
         WRITE(IMPRIM,10001) I,
     %   (M, VALPRE(I), VALPRO(I,M), ERRELA(I,M),M=1,NBM)
 7    CONTINUE
10001 FORMAT(/'VALEUR PROPRE ',I1/
     %(' Maillage ',I1,' Valeur Exacte=',G15.6,
     % ' Valeur Calculee=',G15.6,' Erreur Relative=',F6.2,' %'))
C
C     LE CADRE DE LA FENETRE DE TRACE
      print *
      HMIN   =  1E20
      HMAX   = -1E20
      ERRMIN =  1E20
      ERRMAX = -1E20
      DO 4 M = 1, NBM
         HH = NBNODE( M )
         HMIN = MIN( HMIN, HH )
         HMAX = MAX( HMAX, HH )
         DO 3 I=1,7
            ER = ERRELA(I,M)
            ERRMIN = MIN( ERRMIN, ER )
            ERRMAX = MAX( ERRMAX, ER )
 3       CONTINUE
 4    CONTINUE
      print *,'HMIN=',HMIN,' HMAX=',HMAX,' ERRMIN=',ERRMIN,
     %        '  ERRMAX=',ERRMAX
      HH = ( HMAX - HMIN ) / 8
      RH = ( ERRMAX - ERRMIN ) / 8
C
C     MISE SUR FICHIER.eps du TRACE
      CALL xvinitierps( 1 )
      CALL EFFACE
      CALL FENETRE( HMIN-HH, HMAX+HH, ERRMIN-RH, ERRMAX+RH )
C
C     LA SIGNIFICATION DES AXES
      CALL CHOIXFONTE( 20 )
      CALL TEXTE2D( NCVERT, HMIN-HH/2, ERRMAX+RH/5,
     %             'Relative Error %')
      CALL TEXTE2D( NCMAGE, HMAX-HH/2, ERRMIN-RH/2, 'Vertices Nb' )
C
C     LE TRACE DES AXES 2D
      CALL TRAXE2
C
C     LE TRACE des ERREURS RELATIVES  en fonction de Nb de Sommets
      DO 20 M = 1, NBM
         NC = 0
         DO 10 I=1,7,3
C
C        CHOIX DE LA COULEUR DE TRACE
         NC = NC + 1
C
C        LE POINT A TRACER ( NBNODE( M ), ERRELA(I,M) )
         H  = NBNODE( M )
         ER = ERRELA(I,M)
C        LE NOMBRE DE SOMMETS DU MAILLAGE A DROITE
         CALL ENTIER2D(  NC, H, ER, NBNODE( M ) )
C        LE NO DE LA VALEUR PROPRE A GAUCHE
         CALL ENTIER2D(  NCROUG, H-HH/9, ER, I )
C        LE SYMBOLE
         CALL SYMBOLE2D( NCORAN, H, ER, '*' )
C
 10      CONTINUE
 20   CONTINUE
C
C     LE TITRE DU GRAPHIQUE
      CALL CHOIXFONTE( 25 )
      TITRE = '- Delta u = Ek u  on Omega=[0,1]**6, u=0 on Gamma '
      L = NUDCNB( TITRE )
      CALL TEXTE2D( NCROUG, HMIN+2*HH, ERRMAX+RH/2, TITRE(1:L) )
C
      TITRE='exact k eigv(u,v,w,x,y,z)= sin(n1 PI u) ... sin(n6 PI z) '
      L = NUDCNB( TITRE )
      CALL TEXTE2D( NCROUG, HMIN+2*HH, ERRMAX, TITRE(1:L) )
C
      TITRE='Ek = PI**2 ( n1**2 + ... + n6**2 )  ni integer>0'
      L = NUDCNB( TITRE )
      CALL TEXTE2D( NCROUG, HMIN+2*HH, ERRMAX-Rh/2, TITRE(1:L) )
C
      TITRE='Relative ERROR(1..7,nM) |Eigv ie-Eigv ic| / Eigv ie '
      L = NUDCNB( TITRE )
      CALL TEXTE2D( NCROUG, HMIN+2*HH, ERRMAX-Rh, TITRE(1:L) )
C
C     COPIE DE MEMPX DANS FENETRE
      CALL MEMPXFENETRE
C
C     POUR VIDER LE BUFFER DE X11
      CALL XVVOIR
C
C     MISE SUR FICHIER 6cuberr.eps du TRACE
C     ATTENTION PASSAGE PAR VARIABLE OBLIGATOIRE
      NBC = 7
      CALL xvsauverps( '6eigerr', NBC )
      print *, 'NBC=',NBC
C
C     POUR ATTENDRE UN CLIC SOURIS ET  LIRE LE GRAPHIQUE
      CALL CHOIXFONTE( NPHFCO )
      CALL CLICSO
C
      CALL XVEPAISSEUR( 1 )
      RETURN
      END
