      SUBROUTINE NLSESIKG( PENALI, NDIM,   MNXYZN,
     %                     NBTYEL, MNNPEF, NDPGST,
     %                     MNTPOB, MNTAUX, NUMIOB, NUMAOB, MNDOEL,
     %                     MNTHER, MNTAEL, MNX,
     %                     NORESO, NCODSK, NBRDKG, MNLPLI, MNLPCO, KG,
     %                     NBPTAF, IERR  )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    NLSE avec TESTNL=6:  CALCUL DE LA MATRICE GLOBALE nxn
C -----    [KG] = [Rho/dt - ALPHA LAPLACIEN + N(V0**2+W0**2)]
C          DES ELEMENTS FINIS LAGRANGE 2D OU 3D DE DEGRE 1 OU 2
C          
C          NLSE PROBLEM: COMPUTATION OF THE GLOBAL MATRIX nxn
C          FROM ISOPARAMETRIC LAGRANGE FINITE ELEMENT 2D or 3D
C          WITH POLYNOMIALS OF DEGREE 1 or 2 ON A REFERENCE FINITE ELEMENT
C
C ATTENTION: CETTE PROGRAMMATION SUPPOSE QUE
C  TEMPS, TEMPSINI, TEMPSFIN, PasTemps SONT INITIALISES incl/ctemps.inc
C  LA DENSITE DE MASSE Rho EST INDEPENDANTE   DE T, V, W
C  LA CONDUCTIVITE    Alfa EST INDEPENDANTE   DE T, V, W
C  LE COEFFICIENT D'ECHANGE g EST INDEPENDANT DE T, V, W
C  LES DL IMPOSES POUR V et W PAR DIRICHLET SONT AUX MEMES NOEUDS
C
C ENTREES:
C --------
C PENALI : COEFFICIENT DE PENALISATION DES TEMPERATURES FIXEES
C          ICI PENALI VAUT 1D30 POUR LE PRENDRE EN COMPTE
C          ou 0D0 SI ELLES SONT FIXEES DIRECTEMENT
C NDIM   : DIMENSION DE L'ESPACE DES COORDONNEES 1 ou 2 ou 3 ou 6
C          (SI AXISYMETRIE NDIM=2 X => R>=0 et Y=>Z et Z=0)
C MNXYZN : ADRESSE MCN DE TMS XYZNOEUD DE L'OBJET
C NBTYEL : NOMBRE DE TYPES D'EF DE L'OBJET
C MNNPEF : ADRESSE MCN DU TMC DES ADRESSES MCN DES TABLEAUX NPEF"TYPE EF
C NDPGST : CODE DE STOCKAGE DES SOMMETS POINTS NOEUDS
C MNTPOB : ADRESSE MCN DES TABLEAUX POLYNOMES DE BASE DES TYPES D'EF
C MNTAUX : ADRESSE MCN DES TABLEAUX AUXILAIRES
C NUMIOB : NUMERO MINIMAL DES OBJETS
C NUMAOB : NUMERO MAXIMAL DES OBJETS
C MNDOEL : ADRESSE MCN DES TMS DES DONNEES NLSE DE L'OBJET
C MNTHER : ADRESSE MCN DU TABLEAU DES DONNEES NLSE INTERNES
C          AU MOINS NPIMAX(=128) REELS DOUBLE PRECISION
C MNTAEL : ADRESSE MCN DES TABLEAUX ELEMENTAIRES
C MNX    : ADRESSE MCN DES NBCOOR COORDONNEES DES POINTS DE L'EF COURANT
C MNNODL : ADRESSE MCN DU TABLEAU DES NUMEROS DES NOEUDS DE L'EF COURANT
C          DANS $MEFISTO/incl/cthet.inc INITIALISE dans call thed1t.f
C
C NORESO : 1 STOCKAGE PROFIL DE LA MATRICE KG
C          2 STOCKAGE MORSE  DE LA MATRICE KG
C NCODSK : CODE DE STOCKAGE DE LA MATRICE GLOBALE DE CONDUCTIVITE
C          1 CAR MATRICE SYMETRIQUE
C NBRDKG : NOMBRE DE REELS DOUBLE PRECISION DE LA MATRICE GLOBALE KG
C MNLPLI : ADRESSE MCN DU TABLEAU POINTEUR DES COEFFICIENTS DIAGONAUX
C          DE LA MATRICE PROFIL OU MORSE DE K
C MNLPCO : ADRESSE MCN DU TABLEAU DES NUMEROS DE COLONNES DE LA MATRICE MORSE
C KG     : MATRICE GLOBALE DE NLSE SANS LE TERME NON LINEAIRE
C
C SORTIES:
C --------
C NBPTAF : NOMBRE DE POINTS PAR ARETE OU LE FLUX EST CALCULE
C IERR   : 0 SI PAS D'ERREUR, 1 SI EF NON PROGRAMME, 2 SI EF DEGENERE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & Saint Pierre du Perray  Mars 2014
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/donthe.inc"
      include"./incl/a___npef.inc"
      include"./incl/a___xyznoeud.inc"
      include"./incl/a___vecteur.inc"
      include"./incl/ponoel.inc"
      include"./incl/gsmenu.inc"
      include"./incl/cthet.inc"
C
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1),DMCN(1))
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
C
      INTEGER           NOOBSF(6), NOOBLA(12), NOOBPS(8)
      INTEGER           NUMIOB(4), NUMAOB(4),  MNDOEL(4)
      DOUBLE PRECISION  KG(NBRDKG)
C
      DOUBLE PRECISION  PENALI, D2PI, DELTA
      CHARACTER*4       NOMELE(2)
C
C     LE NOMBRE DE MOTS D'UNE VARIABLE REEL2
      MOREE2 = MOTVAR(6)
      D2PI   = ATAN( 1D0 ) * 8D0
      IERR   = 0
      MNPOLQ = 0
      MNPOLA = 0
      MNPOL  = 0
      MNPOIQ = 0
      MNPOID = 0
      MNPOIA = 0
      MNPDEL = 0
      MNF2   = 0
      MNDPOQ = 0
      MNDPOL = 0
      MNDPOA = 0
C
C     INITIALISATION DE LA MATRICE GLOBALE DU SYSTEME
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10002)
      ELSE
         WRITE(IMPRIM,20002)
      ENDIF
10002 FORMAT(/'CONSTRUCTION de la MATRICE GLOBALE KG=[Rho/dt -ALPHA LAPL
     %ACIEN + N(V0**2+W0**2)]')
20002 FORMAT(/'CONSTRUCTION of KG=[Rho/dt -ALPHA LAPLACIEN + N(V0**2+W0*
     %*2)] GLOBAL MATRIX ')
C
C     INITIALISATION A ZERO DE LA MATRICE GLOBALE
C     -------------------------------------------
      CALL AZEROD( NBRDKG, KG  )
C
C     LA GENERATION DES TABLEAUX ELEMENTAIRES ET DE LA MATRICE
C     ========================================================
C     LE NOMBRE DE POINTS D'INTEGRATION PAR ARETE EN 2D
      NBPTAF = 0
C
C     LA BOUCLE SUR LES TYPES D'ELEMENTS FINIS
      DO 500 NOTYEL = 1 , NBTYEL
C
C        L'ADRESSE DU TABLEAU NPEF"TYPE EF
         MNELE = MCN( MNNPEF - 1 + NOTYEL )
C
C        LE NUMERO DU TYPE DE L'ELEMENT FINI
         NUTYEL = MCN( MNELE + WUTYEL )
         IF( NUTYEL .LE. 4 ) THEN
C           EF AXISYMETRIQUE
            NOAXIS = 1
         ELSE
            NOAXIS = 0
         ENDIF
C
C        LE NOMBRE D'ELEMENTS FINIS DE CE TYPE
         NBELEM = MCN( MNELE + WBELEM )
C
C        L'ADRESSE MCN DES NOEUDS ET POINTS GEOMETRIQUES DES EF DE CE TYPE
         MNNDEL = MNELE + WUNDEL
         MNPGEL = MNNDEL
         IF( NDPGST .GE. 2 ) THEN
C           POINTS DIFFERENTS DES NOEUDS
            MNPGEL = MNPGEL + NBELEM * MCN(MNELE+WBNDEL)
         ENDIF
C
C        LES CARACTERISTIQUES DE L'ELEMENT FINI
         CALL ELNUNM( NUTYEL, NOMELE )
         CALL ELTYCA( NUTYEL )
C
C        LES ADRESSAGES POUR RECUPERER LES INFORMATIONS
C        DU TABLEAU DES POLY(POINTS INTEGRATION), ...
C        ----------------------------------------------
C        SELON LE TYPE DE L'ELEMENT FINI
         GOTO( 10,10,10,10, 1, 1,1, 1, 1, 1,
     %          1, 1,13, 1,10,10,1,10,13,10,
     %         10,10,10,10, 1, 1,1,10,10,10,
     %         10,10,10, 1, 1), NUTYEL
C
C        ERREUR
 1       NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'ERREUR NLSESIKG: TYPE EF '// NOMELE(1)
     %           // NOMELE(2) //' NON PROGRAMME'
         ELSE
            KERR(1) = 'ERROR NLSESIKG: FE TYPE '// NOMELE(1)
     %           // NOMELE(2) //' NOT PROGRAMMED'
         ENDIF
         CALL LEREUR
         IERR = 1
         GOTO 9999
C
C        RECHERCHE DU TABLEAU DE POBA ET PARTAGE EN P ET DP
 10      L = MNTPOB + (NOTYEL-1) * MXPOBA
C
C        LES VALEURS DES POLYNOMES DE L'EF DE DIMENSION - 1
C        EN 2D SEGMENT DE REFERENCE, EN 3D SURFACE DE REFERENCE
C        ......................................................
         IF( NDIM .EQ. 1 ) GOTO 12
C
         IA     = MCN( L  )
C        DIMENSION DE L ESPACE
C        NDIMA  = MCN( IA )
C        NOMBRE DE POLYNOMES D INTERPOLATION
         NBPOLA = MCN( IA + 1 )
C        NOMBRE DE POINTS D INTEGRATION NUMERIQUE
         NPIA   = MCN( IA + 2 )
         NBPTAF = MAX( NBPTAF, NPIA )
C        ADRESSE DES POIDS DES POINTS DE LA FORMULE D INTEGRATION
         MNPOIA = IA + 8
C        ADRESSE DES VALEURS DES POLYNOMES AUX POINTS D INTEGRATION
         MNPOLA = MNPOIA + MCN( IA + 3 ) * MOREE2 * NPIA
C        ADRESSE DES VALEURS DES DERIVEES DES POLYNOMES AUX POINTS D'INTEGRATION
         MNDPOA = MNPOLA + MCN( IA + 4 ) * MOREE2 * NBPOLA * NPIA
C
         IF( NBNSOM .EQ. 5 .OR. NBNSOM .EQ. 6 ) THEN
C           EF AVEC 2 TYPES DE FACES : PAR EXEMPLE LA PYRAMIDE et LE PENTAEDRE
C                                      POLA => TRIANGLE
C                                      POLQ => QUADRANGLE
C
C           LES POLYNOMES DE L'EF DE DIMENSION - 1 POUR LE SECOND TYPE DE
C           FACE=QUADRANGLE POSITIONNE EN 4-EME POSITION (CF SP ETTAEL)
C           ...........................................................
            IA     = MCN( L + 3  )
C           DIMENSION DE L ESPACE
C           NDIMQ  = MCN( IA )
C           NOMBRE DE POLYNOMES D INTERPOLATION
            NBPOLQ = MCN( IA + 1 )
C           NOMBRE DE POINTS D INTEGRATION NUMERIQUE
            NPIQ   = MCN( IA + 2 )
            NBPTAF = MAX( NBPTAF, NPIQ )
C           ADRESSE DES POIDS DES POINTS DE LA FORMULE D INTEGRATION
            MNPOIQ = IA + 8
C           ADRESSE DES VALEURS DES POLYNOMES AUX POINTS D INTEGRATION
            MNPOLQ = MNPOIQ + MCN( IA + 3 ) * MOREE2 * NPIQ
C           ADRESSE DES VALEURS DES DERIVEES DES POLYNOMES AUX POINTS D'INTEGRAT
            MNDPOQ = MNPOLQ + MCN( IA + 4) * MOREE2 * NBPOLQ * NPIQ
C
         ELSE
C
C           POLYNOMES DES FACES DE TYPE 2 = CEUX DE LA FACE DE TYPE 1
C           NDIMQ  = MCN( IA )
C           NOMBRE DE POLYNOMES D INTERPOLATION
            NBPOLQ = NBPOLA
C           NOMBRE DE POINTS D INTEGRATION NUMERIQUE
            NPIQ   = NPIA
C           ADRESSE DES POIDS DES POINTS DE LA FORMULE D INTEGRATION
            MNPOIQ = MNPOIA
C           ADRESSE DES VALEURS DES POLYNOMES AUX POINTS D INTEGRATION
            MNPOLQ = MNPOLA
C           ADRESSE DES VALEURS DES DERIVEES DES POLYNOMES AUX POINTS D'INTEGRAT
            MNDPOQ = MNDPOA
C
         ENDIF
C
C        LES VALEURS DES POLYNOMES DE L'EF DE DIMENSION TOTALE
C        EN 1D SEGMENT UNITE, EN 2D SURFACE DE REFERENCE,
C        EN 3D ou 6D VOLUME DE REFERENCE
C        .....................................................
         L      = L + 1
 12      IA     = MCN( L  )
C        DIMENSION DE L ESPACE
C        NDIM   = MCN( IA )
C        NOMBRE DE POLYNOMES D INTERPOLATION
         NBPOLY = MCN( IA + 1 )
C        NOMBRE DE POINTS D INTEGRATION NUMERIQUE
         NPI    = MCN( IA + 2 )
C        ADRESSE DES POIDS DES POINTS DE LA FORMULE D INTEGRATION
         MNPOID = IA + 8
C        ADRESSE DES VALEURS DES POLYNOMES AUX POINTS D INTEGRATION
         MNPOL  = MNPOID + MCN( IA + 3 ) * MOREE2 * NPI
C        ADRESSE DES VALEURS DES DERIVEES DES POLYNOMES AUX POINTS D'INTEGRATION
         MNDPOL = MNPOL + MCN( IA + 4 ) * MOREE2 * NBPOLY * NPI
C
C        LES TABLEAUX AUXILIAIRES
         MNF1   = MNTAUX
         MNF2   = MNF1   + MOREE2 * NPI
         MNPDEL = MNF1   + MOREE2 * NPI * NDIM
         MNDP   = MNPDEL + MOREE2 * NPI
         MNDFM1 = MNDP   + MOREE2 * NPI * NDIM * NBPOLY
C        AU TOTAL = MNDFM1 + MOREE2 * NPI * NDIM * NDIM
C        MOTAUX = MOREE2 * NPI * ( NDIM + 1 + NDIM * (NDIM+NBPOLY) )
         GOTO 30
C
C        TRIANGLE 2P1D  ET  TETRAEDRE 3P1D
C        =================================
 13      NBPOLY = NDIM + 1
         NPI    = 1
         NPIQ   = 1
         NPIA   = 1
         NBPTAF = MAX( NBPTAF, NPIA )
         MNF1   = MNTAUX
         MNDP   = MNF1 + MOREE2 * NDIM
         MNDFM1 = MNDP + MOREE2 * NDIM * NBPOLY
C        AU TOTAL = MNDFM1 + MOREE2 * NDIM * NDIM
C
C        LES TABLEAUX ELEMENTAIRES: KE, ME ET BE
C        =======================================
 30      MNCAEL = MNTAEL + NBPOLY * (NBPOLY+1) / 2  * MOREE2
C
C        LA BOUCLE SUR LES ELEMENTS FINIS DE CE TYPE NUTYEL
C        ==================================================
         MNDEF1 = MNDOEL(1)
         MNDEF2 = MNDOEL(2)
         MNDEF3 = MNDOEL(3)
         MNDEF4 = MNDOEL(4)
         NBJEUX = 1
         JEU    = 1
C
         DO 200 NUELEM = 1, NBELEM
C
C           LES NOEUDS DE L'ELEMENT FINI NUELEM
C           -----------------------------------
            CALL EFNOEU( MNELE, NUELEM, NBNDEL, MCN(MNNODL) )
C
C           LE NUMERO DE VOLUME  DE L'EF
C           LE NUMERO DE SURFACE DES FACES   DE L'EF
C           LE NUMERO DE LIGNE   DES ARETES  DE L'EF
C           LE NUMERO DE POINT   DES SOMMETS DE L'EF
C           ----------------------------------------
            CALL EFPLSV( MNELE , NUELEM,
     %                   NOVCEL, NOSFEL, NOLAEL, NOPSEL,
     %                   NOOBVC, NOOBSF, NOOBLA, NOOBPS )
C
C           LES COORDONNEES DES NBPOE NOEUDS DE L'ELEMENT FINI NUELEM
C           ---------------------------------------------------------
            CALL EFXYZP( NDIM, MNXYZN, NBELEM, NUELEM, MNPGEL, NBPOE,
     %                   RMCN(MNX) )
C
C           ===============================================================
C           LE CALCUL DES TABLEAUX AUXILIAIRES ET DES TABLEAUX ELEMENTAIRES
C           ===============================================================
            GOTO( 41, 41, 41, 41,  1,  1, 1,  1,  1,  1,
     %             1,  1, 42,  1, 41, 41, 1, 41, 49, 51,
     %            51, 51, 51, 51,  1,  1, 1,  1, 41,  1,
     %            51, 51,  1,  1  ), NUTYEL
C
C           ********************************************
C           2D OU AXISYMETRIQUE LAGRANGE ISOPARAMETRIQUE
C           ********************************************
C           LE CALCUL DES TABLEAUX AUXILIAIRES EN 2D OU AXISYMETRIE
C           -------------------------------------------------------
 41         CALL E12LAG( NBPOLY,     NPI,         MCN(MNPOID),
     %                   MCN(MNPOL), MCN(MNDPOL),
     %                   RMCN(MNX) , MCN(MNF1),   MCN(MNF2),
     %                   MCN(MNPDEL),MCN(MNDP),   MCN(MNDFM1) )
C
C           LA MATRICE ELEMENTAIRE DE CONDUCTIVITE EN 2D OU EN AXISYMETRIE
C           --------------------------------------------------------------
            CALL TR2LAG( D2PI, NOAXIS, RMCN(MNX), PENALI, NBJEUX, JEU,
     %                   NBNSOM, NOOBPS,
     %                   NUMIOB(1), NUMAOB(1), MCN(MNDEF1),
     %                   NBPOLA, NPIA,
     %                   MCN(MNPOIA),MCN(MNPOLA),MCN(MNDPOA),
     %                   NARET,NOOBLA,
     %                   NUMIOB(2),NUMAOB(2),MCN(MNDEF2),
     %                   NBPOLY,NPI,MCN(MNPOL),
     %                   NOOBSF,NUMIOB(3),NUMAOB(3),MCN(MNDEF3),
     %                   MCN(MNF1),MCN(MNF2),MCN(MNPDEL),MCN(MNDP),
     %                   MCN(MNTHER),MCN(MNTHER),MCN(MNTAEL),IERR )
            IF( IERR .NE. 0 ) GOTO 9999
            GOTO 70
C
C           ***************************
C           3D LAGRANGE ISOPARAMETRIQUE
C           ***************************
C           LE CALCUL DES TABLEAUX AUXILIAIRES EN 3D
C           ----------------------------------------
 51         CALL E13LAG( NBPOLY, NPI, MCN(MNPOID),
     %                   MCN(MNPOL),  MCN(MNDPOL),
     %                   RMCN(MNX) ,  MCN(MNF1),
     %                   MCN(MNPDEL), MCN(MNDP), MCN(MNDFM1) )
C
C           LA MATRICE ELEMENTAIRE DE CONDUCTIVITE EN 3D
C           --------------------------------------------
            CALL TR3LAG( NUTYEL, RMCN(MNX), PENALI, NBJEUX, JEU,
     %                   NBNSOM, NOOBPS,
     %                   NUMIOB(1),NUMAOB(1),MCN(MNDEF1),
     %                   NARET,NOOBLA,
     %                   NUMIOB(2),NUMAOB(2),MCN(MNDEF2),
     %                   NBPOLA, NPIA, MCN(MNPOIA),
     %                   MCN(MNPOLA),  MCN(MNDPOA),
     %                   NBPOLQ, NPIQ, MCN(MNPOIQ),
     %                   MCN(MNPOLQ),  MCN(MNDPOQ),
     %                   NFACE, NOOBSF,
     %                   NUMIOB(3), NUMAOB(3), MCN(MNDEF3),
     %                   NBPOLY, NPI, MCN(MNPOL),
     %                   NOOBVC,NUMIOB(4),NUMAOB(4),MCN(MNDEF4),
     %                   MCN(MNF1),   MCN(MNPDEL), MCN(MNDP),
     %                   MCN(MNTHER), MCN(MNTHER), MCN(MNTAEL) )
            GOTO 70
C
C           ************************************
C           2D TRIANGLE TRIA 2P1D LAGRANGE DROIT
C           ************************************
C           LA MATRICE ELEMENTAIRE DE CONDUCTIVITE EN 2D
C           --------------------------------------------
 42         CALL TR2P1D( RMCN(MNX),PENALI, NBJEUX, JEU,
     %                   NOOBPS,NUMIOB(1),NUMAOB(1),MCN(MNDEF1),
     %                   NOOBLA,NUMIOB(2),NUMAOB(2),MCN(MNDEF2),
     %                   NOOBSF,NUMIOB(3),NUMAOB(3),MCN(MNDEF3),
     %                   MCN(MNTHER), MCN(MNTAEL) )
            GOTO 70
C
C           *************************************
C           3D TETRAEDRE TETR 3P1D LAGRANGE DROIT
C           *************************************
 49         CALL E13P1D( RMCN(MNX), MCN(MNF1),
     %                   DELTA, MCN(MNDFM1), MCN(MNDP) )
C
C           LA MATRICE ELEMENTAIRE DE CONDUCTIVITE
C           --------------------------------------
            CALL TR3P1D( RMCN(MNX), DELTA,  MCN(MNDP), PENALI,
     %                   NBJEUX, JEU,
     %                   NOOBPS,NUMIOB(1),NUMAOB(1),MCN(MNDEF1),
     %                   NOOBLA,NUMIOB(2),NUMAOB(2),MCN(MNDEF2),
     %                   NOOBSF,NUMIOB(3),NUMAOB(3),MCN(MNDEF3),
     %                   NOOBVC,NUMIOB(4),NUMAOB(4),MCN(MNDEF4),
     %                   MCN(MNTHER), MCN(MNTAEL) )
C           GOTO 70
C
C           ASSEMBLAGE DU TABLEAU ELEMENTAIRE DE DL RANGES PAR BLOCS
C           DANS LE TABLEAU GLOBAL DE DL RANGES PAR NOEUDS
C           ========================================================
 70         IF( NORESO .EQ. 1 ) THEN
C
C              ASSEMBLAGE DE KE DANS LA MATRICE KG PROFIL GLOBALE
               CALL ASMEPC( NBPOLY, MCN(MNNODL),
     %                      NCODSK, MCN(MNTAEL), MCN(MNTAEL),
     %                      NCODSK, MCN(MNLPLI), KG      )
C
            ELSE IF( NORESO .GE. 2 ) THEN

C              ASSEMBLAGE DE KE DANS LA MATRICE KG MORSE GLOBALE
               CALL ASMEGC( NBPOLY, MCN(MNNODL),
     %                      NCODSK, MCN(MNTAEL), MCN(MNTAEL),
     %                      NCODSK, MCN(MNLPLI), MCN(MNLPCO), KG )
C
            ENDIF
C
C     FIN DE BOUCLE SUR LES EF DE TYPE NUTYEL
 200  CONTINUE
C
C     FIN DE BOUCLE SUR LES TYPES D'EF
 500  CONTINUE
C
 9999 RETURN
      END
