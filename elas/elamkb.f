      SUBROUTINE ELAMKB( IEMG,   IEKG,   IEBG,
     %                   PENALI, D2PI,   NDIM,
     %                   NBTYEL, MNNPEF, NDPGST,
     %                   MNTPOB, MXPOBA, MNTAUX,
     %                   MNXYZP, NUMIOB, NUMAOB, MNDOEL,
     %                   MNELAS, MNFORC, MNTAEL, MNX,   MNIP,   MNNODL,
     %                   MNTEMP, NTDLTE,
     %                   NORESO, MNLPLI, MNLPCO,
     %                   NBRDMG, MG,     NBRDKG, KG,    NTDLDE, MNBG,
     %                   NCODSM, NCODSK, NPIMAX )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DES MATRICES GLOBALES DE MASSE ET RAIDEUR ET
C -----    DU SECOND MEMBRE EN ELASTICITE INSTATIONAIRE
C          POUR DES ELEMENTS LAGRANGE 2D OU 3D DE DEGRE 1 OU 2
C          SUR UN TRIANGLE OU UN QUADRANGLE
C          OU UN TETRAEDRE OU UN PENTAEDRE OU UN HEXAEDRE
C
C ENTREES:
C --------
C IEMG   : 1 SI CALCUL DEMANDE DE LA MATRICE GLOBALE DE MASSE
C          0 SI PAS DE CALCUL  DE LA MATRICE GLOBALE DE MASSE
C IEKG   : 1 SI CALCUL DEMANDE DE LA MATRICE GLOBALE DE RAIDEUR
C          0 SI PAS DE CALCUL  DE LA MATRICE GLOBALE DE RAIDEUR
C IEBG   : 1 SI CALCUL DEMANDE DU SECOND MEMBRE GLOBAL
C          0 SI PAS DE CALCUL  DU SECOND MEMBRE GLOBAL
C
C PENALI : COEFFICIENT DE PENALISATION DE LA CONDITION DE DIRICHLET
C          NON 0D0 SI LE TMS "FIXATION" EST A PENALISER ET
C                  PENALI DOIT ETRE UNE GRANDE VALEUR REELLE DOUBLE
C          0D0 PAS DE PENALISATION DU TMS "FIXATION"
C D2PI   : 2 PI EN REEL DOUBLE PRECISION
C NDIM   : DIMENSION DE L'ESPACE DES COORDONNEES (SI AXISYMETRIE NDIM=2)
C          (X=>R>=0 Y=>Z Z=0)
C
C NBTYEL : NOMBRE DE TYPES D'EF DE L'OBJET
C MNNPEF : ADRESSE MCN DES TABLEAUX NPEF"TYPE EF
C NDPGST : CODE DE STOCKAGE DES SOMMETS POINTS NOEUDS
C MNTPOB : ADRESSE MCN DES TABLEAUX POLYNOMES DE BASE DES TYPES D'EF
C MXPOBA : NOMBRE MAXIMAL DE TABLEAUX POBA PAR TYPE D'EF
C MNTAUX : ADRESSE MCN DES TABLEAUX AUXILIAIRES
C
C MNXYZP : ADRESSE MCN DE TMS XYZSOMMET DE L'OBJET
C NUMIOB : NUMERO MINIMAL DES OBJETS
C NUMAOB : NUMERO MAXIMAL DES OBJETS
C MNDOEL : ADRESSE MCN DES DONNEES DE L'OBJET
C
C MNELAS : ADRESSE MCN DU TABLEAU TENSEUR D'ELASTICITE
C MNFORC : ADRESSE MCN DU TABLEAU DES FORCES EXERCEES EN UN POINT
C MNTAEL : ADRESSE MCN DES TABLEAUX ELEMENTAIRES
C MNX    : ADRESSE MCN DES COORDONNEES DES POINTS DE L'EF COURANT
C MNIP   : NUMERO DES DL PAR COMPOSANTES A PAR NOEUDS
C MNNODL : ADRESSE MCN DU TABLEAU DES NUMEROS DES NOEUDS DE L'EF COURANT
C
C MNTEMP : ADRESSE MCN DE LA TEMPERATURE DU PREMIER NOEUD
C          DU TABLEAU VECTEUR"TEMPERATURE A UTILISER
C NTDLTE : NOMBRE DE DL THERMIQUES DU MAILLAGE
C
C NORESO : 1 STOCKAGE PROFIL DES MATRICES
C          2 STOCKAGE MORSE  DES MATRICES
C MNLPLI : ADRESSE MCN DU TABLEAU POINTEUR DES COEFFICIENTS DIAGONAUX
C          DE LA MATRICE PROFIL OU MORSE DE K
C          SI EN SORTIE NCODSM/=0 ALORS M A MEME PROFIL QUE K
C                                 SINON M EST DIAGONALE
C MNLPCO : ADRESSE MCN DU TABLEAU DES NUMEROS DE COLONNES DE LA MATRICE MORSE
C
C NBRDMG : NOMBRE DE REELS DOUBLE PRECISION DE LA MATRICE GLOBALE MG
C          = NBRDKG SI MATRICE NON DIAGONALE, =NTDLDE SI MATRICE DIAGONALE
C MG     : LA MATRICE DE CAPACITE
C          0 SI IEMG=0
C NBRDKG : NOMBRE DE REELS DOUBLE PRECISION DE LA MATRICE GLOBALE KG
C KG     : LA MATRICE DE RAIDEUR
C          0 SI IEKG=0
C NTDLDE : NOMBRE DE DL DEPLACEMENT DU MAILLAGE
C MNBG   : ADRESSE MCN DU TABLEAU DU VECTEUR GLOBAL SECOND MEMBRE
C          0 SI IEBG=0
C
C NCODSM : CODE DE STOCKAGE DE LA MATRICE GLOBALE DE MASSE
C          0 SI MATRICE DIAGONALE
C          1 SI MATRICE SYMETRIQUE
C         -1 SI MATRICE NON SYMETRIQUE
C         NON INITIALISE SI IEMG=0
C NCODSK : CODE DE STOCKAGE DE LA MATRICE GLOBALE DE RAIDEUR
C          1 SI MATRICE SYMETRIQUE
C         -1 SI MATRICE NON SYMETRIQUE
C         NON INITIALISE SI IEKG=0
C
C SORTIE :
C --------
C NPIMAX : NOMBRE MAXIMAL DE POINTS D'INTEGRATION DES TYPES D'EF
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET   ANALYSE NUMERIQUE UPMC PARIS    MARS   1999
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      include"./incl/a___npef.inc"
      include"./incl/a___xyzpoint.inc"
      include"./incl/a___vecteur.inc"
      include"./incl/ponoel.inc"
      include"./incl/gsmenu.inc"
C
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
C
      include"./incl/pp.inc"
      COMMON            MCN (MOTMCN)
      REAL              RMCN(MOTMCN)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE      (MCN(1),RMCN(1),DMCN(1))
C
      INTEGER           NOOBSF(6),NOOBLA(12),NOOBPS(8)
      INTEGER           NUMIOB(4),NUMAOB(4),MNDOEL(4)
CCC      INTEGER           NOMTAB(5)
      DOUBLE PRECISION  MG(NBRDMG), KG(NBRDKG)
C
      DOUBLE PRECISION  D2PI, DELTA, PENALI
      CHARACTER*4       NOMELE(2)
C
      MNPOLQ = 0
      MNPOLA = 0
      MNPOL  = 0
      MNPOIQ = 0
      MNPOID = 0
      MNPOIA = 0
      MNDPOQ = 0
      MNDPOA = 0
      MNDPOL = 0
C
C     LE NOMBRE DE MOTS D'UNE VARIABLE REEL2
      MOREE2 = MOTVAR(6)
C
C     INITIALISATION DE LA MATRICE DE MASSE
      IF( IEMG .NE. 0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10001)
         ELSE
            WRITE(IMPRIM,20001)
         ENDIF
10001    FORMAT(' CONSTRUCTION de la MATRICE de MASSE   [M]')
20001    FORMAT(' CONSTRUCTION of the MASS MATRIX  [M]')
         CALL AZEROD( NBRDMG, MG )
      ENDIF
C
C     INITIALISATION DE LA MATRICE DE RAIDEUR
      IF( IEKG .NE. 0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10002)
         ELSE
            WRITE(IMPRIM,20002)
         ENDIF
10002    FORMAT(' CONSTRUCTION de la MATRICE de RAIDEUR [K]')
20002    FORMAT(' CONSTRUCTION of STIFFNESS MATRIX [K]')
         CALL AZEROD( NBRDKG, KG )
      ENDIF
C
C     INITIALISATION DU SECOND MEMBRE
      IF( IEBG .NE. 0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10003)
         ELSE
            WRITE(IMPRIM,20003)
         ENDIF
10003    FORMAT(' CONSTRUCTION du VECTEUR SECOND MEMBRE {b}')
20003    FORMAT(' CONSTRUCTION of SECOND MEMBER VECTOR {b}')
         CALL AZEROD( NTDLDE, DMCN( (MNBG+1)/2 ) )
      ENDIF
C
C     LA GENERATION DES TABLEAUX ELEMENTAIRES ET LES MATRICES DE
C     MASSE ET RAIDEUR SONT STOCKEES SOUS FORME PROFIL
C     ==========================================================
C     LE NOMBRE MAXIMAL DE POINTS D'INTEGRATION SUR UN ELEMENT
      NPIMAX = 0
      NDSM   = 1
C
C     LA BOUCLE SUR LES TYPES D'ELEMENTS FINIS
      DO 800 NOTYEL = 1, NBTYEL
C
C        L'ADRESSE DU TABLEAU NPEF"TYPE EF
         MNELE = MCN( MNNPEF - 1 + NOTYEL )
C        LE NUMERO DU TYPE DE L'ELEMENT FINI
         NUTYEL = MCN( MNELE + WUTYEL )
C
C        LE NOMBRE D'ELEMENTS FINIS DE CE TYPE
         NBELEM = MCN( MNELE + WBELEM )
C
C        L'ADRESSE MCN DES NOEUDS ET POINTS GEOMETRIQUES DES ELEMENTS FINIS
         MNNDEL = MNELE + WUNDEL
         MNPGEL = MNNDEL
         IF( NDPGST .GE. 2 ) THEN
            MNPGEL = MNPGEL + NBELEM * MCN(MNELE+WBNDEL)
         ENDIF
C
         IF( NUTYEL .LE. 4 ) THEN
            NOAXIS = 1
         ELSE
            NOAXIS = 0
         ENDIF
C
C        LES CARACTERISTIQUES DE CE TYPE D'ELEMENT FINI
C        ON TROUVE: NBPOE, NBNOE, NARET
         CALL ELNUNM( NUTYEL, NOMELE )
         CALL ELTYCA( NUTYEL )
C
C        LES ADRESSAGES POUR RECUPERER LES INFORMATIONS
C        DU TABLEAU DES POLY(POINTS INTEGRATION), ...
C        ----------------------------------------------
C        SELON LE TYPE DE L'ELEMENT FINI
         GOTO( 3, 3, 3, 3, 1, 1, 1, 1, 1, 1,
     %         1, 1, 3, 1, 3, 3, 1, 3, 2, 3,
     %         3, 3, 3, 3, 1, 1, 1, 1, 2, 1,
     %         3, 3, 1 ), NUTYEL
C
C        ERREUR
 1       NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'ERREUR: ELEMENT FINI '// NOMELE(1)
     %              // NOMELE(2) //' NON PROGRAMME'
         ELSE
            KERR(1) = 'ERROR: FINITE ELEMENT '// NOMELE(1)
     %              // NOMELE(2) //' NOT PROGRAMMED'
         ENDIF
         CALL LEREUR
         RETURN
C
C        TRIANGLE 2P1D ET TETRAEDRE 3P1D
C        ===============================
 2       NBPOLY = NDIM + 1
         NPI    = NBPOLY
         NPIMAX = MAX( NPIMAX, NPI )
         NBPOLQ = NDIM
         NPIA   = NDIM
         GOTO 9
C
C        LES AUTRES EF LAGRANGE ISOPARAMETRIQUES
C        =======================================
C        RECHERCHE DU TABLEAU DE POBA ET PARTAGE EN P ET DP
 3       L = MNTPOB + (NOTYEL-1) * MXPOBA
C
C        L'ELEMENT FINI : SURFACE DE REFERENCE
         MN     = MCN( L  )
C        DIMENSION DE L ESPACE
C        NDIMA  = MCN( MN )
C        NOMBRE DE POLYNOMES D INTERPOLATION
         NBPOLA = MCN( MN + 1 )
C        NOMBRE DE POINTS D INTEGRATION NUMERIQUE
         NPIA   = MCN( MN + 2 )
C        ADRESSE DES POIDS DES POINTS DE LA FORMULE D INTEGRATION
         MNPOIA = MN + 8
C        ADRESSE DES VALEURS DES POLYNOMES AUX POINTS D INTEGRATION
         MNPOLA = MNPOIA + MCN( MN + 3 ) * MOREE2 * NPIA
C        ADRESSE DES VALEURS DES DERIVEES DES POLYNOMES AUX POINTS D I
         MNDPOA = MNPOLA + MCN(MN + 4) * MOREE2 * NBPOLA * NPIA
C
         IF( NBNSOM .EQ. 5 .OR. NBNSOM .EQ. 6 ) THEN
C           EF AVEC 2 TYPES DE FACES : PAR EXEMPLE PYRAMIDE et PENTAEDRE
C                                      POLA => TRIANGLE
C                                      POLQ => QUADRANGLE
C
C           LES POLYNOMES DE L'EF DE DIMENSION - 1 POUR LE SECOND TYPE DE
C           FACE=QUADRANGLE POSITIONNE EN 4-EME POSITION (CF SP EETAEL)
            MN     = MCN( L + 2  )
C           DIMENSION DE L ESPACE
C           NDIMQ  = MCN( MN )
C           NOMBRE DE POLYNOMES D INTERPOLATION
            NBPOLQ = MCN( MN + 1 )
C           NOMBRE DE POINTS D INTEGRATION NUMERIQUE
            NPIQ   = MCN( MN + 2 )
            NBPTAF = MAX( NBPTAF, NPIQ )
C           ADRESSE DES POIDS DES POINTS DE LA FORMULE D INTEGRATION
            MNPOIQ = MN + 8
C           ADRESSE DES VALEURS DES POLYNOMES AUX POINTS D INTEGRATION
            MNPOLQ = MNPOIQ + MCN( MN + 3 ) * MOREE2 * NPIQ
C           ADRESSE DES VALEURS DES DERIVEES DES POLYNOMES AUX POINTS D'INTEGRAT
            MNDPOQ = MNPOLQ + MCN( MN + 4) * MOREE2 * NBPOLQ * NPIQ
C
         ELSE
C
C           POLYNOMES DES FACES DE TYPE 2 = CEUX DE LA FACE DE TYPE 1
C           NDIMQ  = MCN( MN )
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
C        L'ELEMENT FINI : VOLUME DE REFERENCE
         L      = L + 1
         MN     = MCN( L  )
C        DIMENSION DE L ESPACE
C        NDIM   = MCN( MN )
C        NOMBRE DE POLYNOMES D INTERPOLATION
         NBPOLY = MCN( MN + 1 )
C        NOMBRE DE POINTS D INTEGRATION NUMERIQUE
         NPI    = MCN( MN + 2 )
         NPIMAX = MAX( NPIMAX, NPI )
C        ADRESSE DES POIDS DES POINTS DE LA FORMULE D INTEGRATION
         MNPOID = MN + 8
C        ADRESSE DES VALEURS DES POLYNOMES AUX POINTS D INTEGRATION
         MNPOL  = MNPOID + MCN( MN + 3 ) * MOREE2 * NPI
C        ADRESSE DES VALEURS DES DERIVEES DES POLYNOMES AUX POINTS D I
         MNDPOL = MNPOL + MCN(MN + 4) * MOREE2 * NBPOLY * NPI
C
C        LE NOMBRE DE DEGRES DE LIBERTE (DL) DE L'EF
 9       NBDL   = NBPOLY * NDIM
C
C        LES TABLEAUX AUXILIAIRES
         MNF1   = MNTAUX
         MNF2   = MNF1   + MOREE2 * NPI
         MNPDEL = MNF1   + MOREE2 * NPI * NDIM
         MNDP   = MNPDEL + MOREE2 * NPI
         MNDFM1 = MNDP   + MOREE2 * NDIM * NBPOLY * NPI
C        AU TOTAL = MNDFM1 + MOREE2 * NDIM * NDIM * NPI
CCC         MMM1 = MNDFM1 + MOREE2 * NDIM * NDIM * NPI
C
         IF( NOAXIS .EQ. 1 ) THEN
C           LA MATRICE ELEMENTAIRE DE RIGIDITE EN AXISYMETRIE
            MNG1 = MNDFM1
            MNG2 = MNG1 + MOREE2 * NBPOLY * NBPOLY
            MNG3 = MNG2 + MOREE2 * NBPOLY
C           AU TOTAL = MNG3 + MOREE2 * 2 * NBPOLY
C           NOMBRE DE VARIABLES DOUBLE PRECISION AUXILIAIRES
         ELSE
C           LA MATRICE ELEMENTAIRE DE RIGIDITE EN 2D ET 3D
            MNG1 = MNDFM1
            MNG2 = MNG1 + MOREE2 * NBPOLY * (NBPOLY+1) / 2
            MNG3 = MNG2 + MOREE2 * NDIM * NBPOLY
C           AU TOTAL  = MNG3 + MOREE2 * NBPOLY * NBPOLY
C           NOMBRE DE VARIABLES DOUBLE PRECISION AUXILIAIRES
         ENDIF
CCC         MMM2 = MNG3 + MOREE2 * NBPOLY * NBPOLY
C
C        LES SECONDS MEMBRES ELEMENTAIRES
         MNGS1 = MNDFM1
         MNGS2 = MNGS1 + MOREE2 * NBPOLY * NDIM
         MNGS3 = MNGS2 + MOREE2 * NBPOLY * NDIM * NBPOLY
C        AU TOTAL = MNGS3 + MOREE2 * NBPOLQ * (NBPOLQ+1) / 2
CCC         MMM3 = MNGS3 + MOREE2 * NBPOLQ * (NBPOLQ+1) / 2
CCC         WRITE(IMPRIM,*) MAX(MMM1,MMM2,MMM3)-MNTAUX,
CCC     %   ' MOTS AUXILIAIRES NECESSAIRES POUR EF ',NOMELE
CCCC        LES CARACTERISTIQUES DES TABLEAUX DE L'ELEMENT FINI
CCCC        AU DESSOUS A SUPPRIMER SI CORRECT
CCC         CALL EETAEL( NUTYEL,NDSM,NBDL,NCODEM,NTPOBA,NOMTAB,MOTAUX )
CCC         WRITE(IMPRIM,*) MOTAUX*MOREE2,' MOTS AUXILIAIRES DECLARES'
C
C        LA LISTE DES NUMEROS DES DL PAR NOEUDS
C        ======================================
         CALL CALCIP( NDIM, NBPOLY, MCN(MNIP) )
C
C        LA BOUCLE SUR LES ELEMENTS FINIS DE CE TYPE NUTYEL
C        ==================================================
         DO 790 NUELEM = 1, NBELEM
C
C           LES NOEUDS DE L'ELEMENT FINI
C           LE NUMERO DE VOLUME  DE L'EF
C           LE NUMERO DE SURFACE DES FACES   DE L'EF
C           LE NUMERO DE LIGNE   DES ARETES  DE L'EF
C           LE NUMERO DE POINT   DES SOMMETS DE L'EF
C           ----------------------------------------
            CALL EFPLSV( MNELE , NUELEM ,
     %                   NOVCEL, NOSFEL, NOLAEL, NOPSEL ,
     %                   NOOBVC, NOOBSF, NOOBLA, NOOBPS )
C
C           LES COORDONNEES DES POINTS DE L'ELEMENT FINI
C           --------------------------------------------
            DO I=1,NBPOE
C              LE NUMERO DU I-EME POINT DE L'EF NUELEM
               N    = MCN( MNPGEL-1 + NUELEM + NBELEM*(I-1) )
               MNCE = MNXYZP + WYZPOI + (N-1) * 3
               L    = MNX - 1 + I
               RMCN( L         ) = RMCN( MNCE   )
               RMCN( L + NBPOE ) = RMCN( MNCE+1 )
               IF( NDIM .EQ. 3 ) THEN
                   RMCN( L + 2*NBPOE ) = RMCN( MNCE+2 )
               ENDIF
            ENDDO
C
C           LE NUMERO GLOBAL DES DEGRES DE LIBERTE DE L'ELEMENT FINI
C           --------------------------------------------------------
            L = 0
            DO 35 I=1,NBNOE
C              LE NUMERO DU I-EME NOEUD DE L'ELEMENT FINI NUELEM
               N = MCN( MNNDEL-1 + NUELEM + NBELEM*(I-1) )
               N = ( N - 1 ) * NDIM
               DO 30 J=1,NDIM
                  MCN( MNNODL + L ) = N + J
                  L = L + 1
 30            CONTINUE
 35         CONTINUE
C
C           ===============================================
C           LA MATRICE ELEMENTAIRE DE RAIDEUR EN ELASTICITE
C           SELON LE TYPE DE L'ELEMENT FINI
C           ===============================================
            IF( IEKG .EQ. 0 ) GOTO 100
            GOTO(  41,  41, 41,  41,9000,9000,9000,9000,9000,9000,
     %           9000,9000, 43,9000,  44,  44,9000,  44,  49,  50,
     %             50,  50, 50,  50,9000,9000,9000,9000,  43,9000,
     %             50,  50 ),NUTYEL
C
C           LE CALCUL DES TABLEAUX AUXILIAIRES EN AXISYMETRIE
C           -------------------------------------------------
 41         CALL E1LAXI( D2PI,NBPOLY,NPI,MCN(MNPOID),
     %                   MCN(MNPOL),MCN(MNDPOL),
     %                   RMCN(MNX),
     %                   MCN(MNF1),MCN(MNF2),
     %                   MCN(MNPDEL),MCN(MNDP),MCN(MNDFM1) )
            IF( NBPOLY .LE. 0 ) THEN
C               X=R<0 => ERREUR
                IERR = 7
                GOTO 9999
            ENDIF
C
C           LA MATRICE ELEMENTAIRE DE RIGIDITE EN AXISYMETRIE
C           -------------------------------------------------
            CALL ERLAXI(D2PI,   PENALI,    NBNSOM,    NARET, 
     %                  NOOBPS, NUMIOB(1), NUMAOB(1), MCN(MNDOEL(1)), 
     %                  NOOBLA, NUMIOB(2), NUMAOB(2), MCN(MNDOEL(2)), 
     %                  NOOBSF, NUMIOB(3), NUMAOB(3), MCN(MNDOEL(3)), 
     %                  NBPOLA, NBPOLY,    NPI, MCN(MNPOL ), 
     %                  MCN(MNIP),   MCN(MNF1), MCN(MNF2), 
     %                  MCN(MNPDEL), MCN(MNDP), 
     %                  MCN(MNG1),   MCN(MNG2), MCN(MNG3), 
     %                  MCN(MNELAS), MCN(MNTAEL) )
            GOTO 60
C
C           LA MATRICE ELEMENTAIRE DE RIGIDITE TRIA 2P1D
C           --------------------------------------------
 43         CALL ER2P1D( RMCN(MNX),PENALI,
     %                   NOOBPS,NUMIOB(1),NUMAOB(1),MCN(MNDOEL(1)),
     %                   NOOBLA,NUMIOB(2),NUMAOB(2),MCN(MNDOEL(2)),
     %                   NOOBSF,NUMIOB(3),NUMAOB(3),MCN(MNDOEL(3)),
     %                   MCN(MNELAS),MCN(MNTAEL))
            GOTO 60
C
C           LE CALCUL DES TABLEAUX AUXILIAIRES EN ELASTICITE 2D
C           ---------------------------------------------------
 44         CALL E12LAG( NBPOLY,      NPI,        MCN(MNPOID),
     %                   MCN(MNPOL),  MCN(MNDPOL),
     %                   RMCN(MNX),
     %                   MCN(MNF1),   MCN(MNF2),
     %                   MCN(MNPDEL), MCN(MNDP),  MCN(MNDFM1) )
C
C           LA MATRICE ELEMENTAIRE DE RIGIDITE EN ELASTICITE 2D
C           ---------------------------------------------------
            CALL ER2LAG(D2PI,   PENALI,    NBNSOM,    NARET, 
     %                  NOOBPS, NUMIOB(1), NUMAOB(1), MCN(MNDOEL(1)),
     %                  NOOBLA, NUMIOB(2), NUMAOB(2), MCN(MNDOEL(2)),
     %                  NOOBSF, NUMIOB(3), NUMAOB(3), MCN(MNDOEL(3)),
     %                  NBPOLA, NBPOLY, NPI , MCN(MNPOL), 
     %                  MCN(MNIP),   MCN(MNF1), MCN(MNF2), 
     %                  MCN(MNPDEL), MCN(MNDP), 
     %                  MCN(MNG1),   MCN(MNG2), MCN(MNG3), 
     %                  MCN(MNELAS), MCN(MNTAEL) )
            GOTO 60
C
C           LE CALCUL DES TABLEAUX AUXILIAIRES EN ELASTICITE 3D TETRAEDRE 3P1D
C           ------------------------------------------------------------------
 49         CALL E13P1D( RMCN(MNX), MCN(MNF1),
     %                   DELTA, MCN(MNDFM1), MCN(MNDP) )
C
C           LA MATRICE ELEMENTAIRE DE RIGIDITE EN ELASTICITE 3D TETRAEDRE 3P1D
C           ------------------------------------------------------------------
            CALL ER3P1D( RMCN(MNX), PENALI,
     %                   NOOBPS, NUMIOB(1), NUMAOB(1), MCN(MNDOEL(1)),
     %                   NOOBLA, NUMIOB(2), NUMAOB(2), MCN(MNDOEL(2)),
     %                   NOOBSF, NUMIOB(3), NUMAOB(3), MCN(MNDOEL(3)),
     %                   NOOBVC, NUMIOB(4), NUMAOB(4), MCN(MNDOEL(4)),
     %                   MCN(MNELAS),
     %                   DELTA,  MCN(MNDP), MCN(MNIP),
     %                   MCN(MNTAEL) )
            GOTO 60
C
C           LES TABLEAUX AUXILIAIRES ELASTICITE 3D LAGRANGE ISOPARAMETRIQUE
C           ---------------------------------------------------------------
 50         CALL E13LAG( NBPOLY,      NPI,         MCN(MNPOID),
     %                   MCN(MNPOL),  MCN(MNDPOL),
     %                   RMCN(MNX),   MCN(MNF1),
     %                   MCN(MNPDEL), MCN(MNDP),   MCN(MNDFM1) )
C
C           LA MATRICE ELEMENTAIRE DE RIGIDITE EN ELASTICITE 3D
C           ---------------------------------------------------
            CALL ER3LAG( PENALI,
     %                   NOOBPS, NUMIOB(1), NUMAOB(1), MCN(MNDOEL(1)),
     %                   NOOBLA, NUMIOB(2), NUMAOB(2), MCN(MNDOEL(2)),
     %                   NOOBSF, NUMIOB(3), NUMAOB(3), MCN(MNDOEL(3)),
     %                   NOOBVC, NUMIOB(4), NUMAOB(4), MCN(MNDOEL(4)),
     %                   NUTYEL, NBPOLY,    NPI,       MCN(MNPOL),
     %                   MCN(MNG1), MCN(MNG2), MCN(MNG3),
     %                   MCN(MNIP), MCN(MNF1), MCN(MNPDEL), MCN(MNDP),
     %                   MCN(MNELAS),
     %                   MCN(MNTAEL) )
C
C           ====================================================================
C           ASSEMBLAGE DE LA MATRICE ELEMENTAIRE DE RIGIDITE DANS LA MATRICE GLO
C           ====================================================================
 60         IF ( NORESO .EQ. 1 ) THEN
C
C              STOCKAGE PROFIL DES MATRICES
C              ----------------------------
               CALL ASMEPC( NBDL,   MCN(MNNODL),
     %                      NCODSK, MCN(MNTAEL), MCN(MNTAEL),
     %                      NCODSK, MCN(MNLPLI), KG  )
C
            ELSE IF( NORESO .EQ. 2 ) THEN
C
C              STOCKAGE MORSE DES MATRICES POUR LE GC
C              --------------------------------------
               CALL ASMEGC( NBDL,   MCN(MNNODL),
     %                      NCODSK, MCN(MNTAEL), MCN(MNTAEL),
     %                      NCODSK, MCN(MNLPLI), MCN(MNLPCO), KG )
            ENDIF
C
C           =============================================
C           LA MATRICE ELEMENTAIRE DE MASSE EN ELASTICITE
C           SELON LE TYPE DE L'ELEMENT FINI
C           =============================================
 100        IF( IEMG .EQ. 0 ) GOTO 200
            GOTO(  144, 144, 144,  144, 9000,9000,9000,9000,9000,9000,
     %            9000,9000, 142, 9000,  144, 144,9000, 144, 149, 150,
     %             150, 150, 150,  150, 9000,9000,9000,9000, 142,9000,
     %             150, 150 ),NUTYEL
C
C           LA MATRICE ELEMENTAIRE DE MASSE EN ELASTICITE TRIA 2P1D
C           -------------------------------------------------------
 142        CALL EM2P1D( RMCN(MNX),
     %                   NOOBSF(1),NUMIOB(3),NUMAOB(3),MCN(MNDOEL(3)),
     %                   NCODEM, MCN(MNTAEL) )
CCC
CCCC           TRANSFORMATION MATRICE DIAGONALE->MATRICE SYMETRIQUE
CCC            CALL DIASYM( 6, MCN(MNTAEL), MCN(MNTAEL) )
CCC            NCODEM = 1
            GOTO 160
CCCC
CCCC           LA MATRICE ELEMENTAIRE DE MASSE EN ELASTICITE QUAD 2Q1C NUTYEL=16
CCCC           -------------------------------------------------------
CCC 143        CALL EM2Q1C( RMCN(MNX),
CCC     %                   NOOBSF(1),NUMIOB(3),NUMAOB(3),MCN(MNDOEL(3)),
CCC     %                   NCODEM, MCN(MNTAEL) )
CCCC           TRANSFORMATION MATRICE DIAGONALE->MATRICE SYMETRIQUE
CCC            CALL DIASYM( 8, MCN(MNTAEL), MCN(MNTAEL) )
CCC            NCODEM = 1
CCC            GOTO 160
C
C           LA MATRICE ELEMENTAIRE DE MASSE EN ELASTICITE AXISYMETRIQUE OU 2D
C           -----------------------------------------------------------------
 144        CALL E12LAG( NBPOLY,      NPI,         MCN(MNPOID),
     %                   MCN(MNPOL),  MCN(MNDPOL),
     %                   RMCN(MNX),
     %                   MCN(MNF1),   MCN(MNF2),
     %                   MCN(MNPDEL), MCN(MNDP),   MCN(MNDFM1) )
C
            CALL EM2LAG( D2PI,   NOAXIS,
     %                   NOOBSF, NUMIOB(3), NUMAOB(3), MCN(MNDOEL(3)),
     %                   NBPOLY, NPI,       MCN(MNPOL),
     %                   MCN(MNF1), MCN(MNF2), MCN(MNPDEL), MCN(MNELAS),
     %                   NCODEM, MCN(MNTAEL) )
            GOTO 160
C
C           LA MATRICE ELEMENTAIRE DE MASSE EN ELASTICITE 3D TETRAEDRE 3P1D
C           ---------------------------------------------------------------
 149        CALL EM3P1D( RMCN(MNX), DELTA,
     %                   NOOBVC, NUMIOB(4), NUMAOB(4), MCN(MNDOEL(4)),
     %                   NCODEM, MCN(MNTAEL) )
CCC
CCCC           TRANSFORMATION MATRICE DIAGONALE->MATRICE SYMETRIQUE
CCC            CALL DIASYM( 12, MCN(MNTAEL), MCN(MNTAEL) )
CCC            NCODEM = 1
            GOTO 160
C
C           LA MATRICE ELEMENTAIRE DE MASSE EN ELASTICITE 3D
C           ------------------------------------------------
 150        CALL E13LAG( NBPOLY,      NPI,        MCN(MNPOID),
     %                   MCN(MNPOL),  MCN(MNDPOL),
     %                   RMCN(MNX),   MCN(MNF1),
     %                   MCN(MNPDEL), MCN(MNDP),  MCN(MNDFM1) )
C
            CALL EM3LAG( NOOBVC, NUMIOB(4), NUMAOB(4), MCN(MNDOEL(4)),
     %                   NBPOLY, NPI, MCN(MNPOL),
     %                   MCN(MNF1), MCN(MNPDEL), MCN(MNELAS),
     %                   NCODEM,    MCN(MNTAEL) )
C
C           ====================================================================
C           ASSEMBLAGE DE LA MATRICE ELEMENTAIRE DE MASSE DANS LA MATRICE GLOBAL
C           ====================================================================
 160        IF ( NORESO .EQ. 1 ) THEN
C
C              STOCKAGE PROFIL DES MATRICES
C              ----------------------------
               CALL ASMEPC( NBDL,   MCN(MNNODL),
     %                      NCODEM, MCN(MNTAEL), MCN(MNTAEL),
     %                      NCODSM, MCN(MNLPLI), MG  )
C
            ELSE IF( NORESO .EQ. 2 ) THEN
C
C              STOCKAGE MORSE DES MATRICES POUR LE GC
C              --------------------------------------
               CALL ASMEGC( NBDL,   MCN(MNNODL),
     %                      NCODEM, MCN(MNTAEL), MCN(MNTAEL),
     %                      NCODSM, MCN(MNLPLI), MCN(MNLPCO), MG )
            ENDIF
C
C           ===============================================
C           LA MATRICE ELEMENTAIRE DE RAIDEUR EN ELASTICITE
C           SELON LE TYPE DE L'ELEMENT FINI
C           ===============================================
 200        IF( IEBG .EQ. 0 ) GOTO 790
            GOTO( 620, 620, 620, 620,9000,9000,9000,9000,9000,9000,
     %           9000,9000, 630,9000, 630, 630,9000, 630, 639, 650,
     %            650, 650, 650, 650,9000,9000,9000,9000, 613,9000,
     %            650, 650 ),NUTYEL
C
C           LES SECONDS MEMBRES ELEMENTAIRES EN AXISYMETRIE
C           -----------------------------------------------
 613        CALL ES2P1D(RMCN(MNX),PENALI,
     %                  NOOBPS,NUMIOB(1),NUMAOB(1),MCN(MNDOEL(1)),
     %                  NOOBLA,NUMIOB(2),NUMAOB(2),MCN(MNDOEL(2)),
     %                  NOOBSF,NUMIOB(3),NUMAOB(3),MCN(MNDOEL(3)),
     %                  NUELEM,NBELEM,MCN(MNNDEL),
     %                  MNTEMP,NTDLTE,MCN(MNTEMP),
     %                  MCN(MNTAEL) )
            GOTO 660
C
C           LE CALCUL DES TABLEAUX AUXILIAIRES EN AXISYMETRIE
C           -------------------------------------------------
 620        CALL E1LAXI( D2PI,NBPOLY,NPI,MCN(MNPOID),
     %                   MCN(MNPOL),MCN(MNDPOL),
     %                   RMCN(MNX),
     %                   MCN(MNF1),MCN(MNF2),
     %                   MCN(MNPDEL),MCN(MNDP),MCN(MNDFM1) )
            IF( NBPOLY .LE. 0 ) THEN
C               X=R<0 => ERREUR
                IERR = 7
                RETURN
            ENDIF
C
C           LES SECONDS MEMBRES ELEMENTAIRES EN AXISYMETRIE
C           -----------------------------------------------
            CALL ESLAXI(D2PI,RMCN(MNX),PENALI,NBNSOM,NARET,NDSM,
     %                  NOOBPS,NUMIOB(1),NUMAOB(1),MCN(MNDOEL(1)),
     %                  NOOBLA,NUMIOB(2),NUMAOB(2),MCN(MNDOEL(2)),
     %                  NOOBSF,NUMIOB(3),NUMAOB(3),MCN(MNDOEL(3)),
     %                  NBPOLA,NPIA,MCN(MNPOIA),MCN(MNPOLA),MCN(MNDPOA),
     %                  NBPOLY,NPI ,MCN(MNPOL ),
     %                  NUELEM,NBELEM,MCN(MNNDEL),
     %                  MNTEMP,NTDLTE,MCN(MNTEMP),
     %                  MCN(MNELAS),MCN(MNFORC),MCN(MNFORC),
     %                  MCN(MNGS1),MCN(MNGS2), MCN(MNF1),MCN(MNF2),
     %                  MCN(MNPDEL),MCN(MNDP),MCN(MNIP),
     %                  MCN(MNTAEL) )
            GOTO 660
C
C           LE CALCUL DES TABLEAUX AUXILIAIRES EN ELASTICITE 2D
C           ---------------------------------------------------
 630        CALL E12LAG( NBPOLY,      NPI,         MCN(MNPOID),
     %                   MCN(MNPOL),  MCN(MNDPOL),
     %                   RMCN(MNX),
     %                   MCN(MNF1),   MCN(MNF2),
     %                   MCN(MNPDEL), MCN(MNDP),   MCN(MNDFM1) )
C
C           LES SECONDS MEMBRES ELEMENTAIRES EN ELASTICITE 2D
C           -------------------------------------------------
            CALL ES2LAG(RMCN(MNX),PENALI,NBNSOM,NARET,NDSM,
     %                  NOOBPS,NUMIOB(1),NUMAOB(1),MCN(MNDOEL(1)),
     %                  NOOBLA,NUMIOB(2),NUMAOB(2),MCN(MNDOEL(2)),
     %                  NOOBSF,NUMIOB(3),NUMAOB(3),MCN(MNDOEL(3)),
     %                  NBPOLA,NPIA,MCN(MNPOIA),MCN(MNPOLA),MCN(MNDPOA),
     %                  NBPOLY,NPI ,MCN(MNPOL ),
     %                  NUELEM,NBELEM,MCN(MNNDEL),
     %                  MNTEMP,NTDLTE,MCN(MNTEMP),
     %                  MCN(MNELAS),MCN(MNFORC),MCN(MNFORC),
     %                  MCN(MNGS1),MCN(MNGS2), MCN(MNF1),MCN(MNF2),
     %                  MCN(MNPDEL),MCN(MNDP),MCN(MNIP),
     %                  MCN(MNTAEL) )
            GOTO 660
C
C           LE CALCUL DES TABLEAUX AUXILIAIRES EN ELASTICITE 3D TETRAEDRE 3P1D
C           ------------------------------------------------------------------
 639        CALL E13P1D( RMCN(MNX), MCN(MNF1),
     %                   DELTA, MCN(MNDFM1), MCN(MNDP) )
C
C           LES SECONDS MEMBRES ELEMENTAIRES EN ELASTICITE 3D TETRAEDRE 3P1D
C           ----------------------------------------------------------------
            CALL ES3P1D( RMCN(MNX), PENALI,
     %                   NOOBPS, NUMIOB(1), NUMAOB(1), MCN(MNDOEL(1)),
     %                   NOOBLA, NUMIOB(2), NUMAOB(2), MCN(MNDOEL(2)),
     %                   NOOBSF, NUMIOB(3), NUMAOB(3), MCN(MNDOEL(3)),
     %                   NOOBVC, NUMIOB(4), NUMAOB(4), MCN(MNDOEL(4)),
     %                   MCN(MNELAS), MCN(MNFORC),
     %                   NUELEM, NBELEM, MCN(MNNDEL),
     %                   MNTEMP, NTDLTE, MCN(MNTEMP),
     %                   DELTA,  MCN(MNDP),
     %                   MCN(MNTAEL) )
            GOTO 660
C
C           LES TABLEAUX AUXILIAIRES ELASTICITE 3D LAGRANGE ISOPARAMETRIQUE
C           ---------------------------------------------------------------
 650        CALL E13LAG( NBPOLY, NPI, MCN(MNPOID),
     %                   MCN(MNPOL),  MCN(MNDPOL),
     %                   RMCN(MNX),   MCN(MNF1),
     %                   MCN(MNPDEL), MCN(MNDP), MCN(MNDFM1) )
C
C           LES SECONDS MEMBRES ELEMENTAIRES EN 3D
C           --------------------------------------
            CALL ES3LAG( RMCN(MNX), PENALI, NDSM,
     %                   NOOBPS, NUMIOB(1), NUMAOB(1), MCN(MNDOEL(1)),
     %                   NOOBLA, NUMIOB(2), NUMAOB(2), MCN(MNDOEL(2)),
     %                   NOOBSF, NUMIOB(3), NUMAOB(3), MCN(MNDOEL(3)),
     %                   NOOBVC, NUMIOB(4), NUMAOB(4), MCN(MNDOEL(4)),
     %                   NBPOLY, NPI,  MCN(MNPOL),
     %                   NBPOLA, NPIA, MCN(MNPOIA),
     %                   MCN(MNPOLA),  MCN(MNDPOA),
     %                   NBPOLQ, NPIQ, MCN(MNPOIQ),
     %                   MCN(MNPOLQ),  MCN(MNDPOQ),
     %                   MCN(MNELAS),  MCN(MNFORC), MCN(MNFORC),
     %                   MCN(MNGS1),   MCN(MNGS2),
     %                   NUTYEL, NUELEM, NBELEM,    MCN(MNNDEL),
     %                   MNTEMP, NTDLTE, MCN(MNTEMP),
     %                   MCN(MNF1), MCN(MNPDEL), MCN(MNDP),
     %                   MCN(MNTAEL) )
C
C           ASSEMBLAGE DU SECOND MEMBRE ELEMENTAIRE
C           =======================================
 660        CALL ASBEBG( NTDLDE, NDSM, NBDL, MCN(MNNODL), MCN(MNTAEL),
     %                   MCN(MNBG) )
C
 790     CONTINUE
C
 800  CONTINUE
      RETURN
C
C     ERREUR
 9000 NBLGRC(NRERR) = 1
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'ERREUR: ELEMENT FINI NON PROGRAMME'
      ELSE
         KERR(1) = 'ERROR: FINITE ELEMENT NOT PROGRAMMED'
      ENDIF
      CALL LEREUR
      RETURN
C
C     ERREUR R<0
 9999 NBLGRC(NRERR) = 1
      IF( LANGAG .EQ. 0 ) THEN
         KERR(1) = 'ERREUR: RAYON<0 en PROBLEME AXISYMETRIQUE'
      ELSE
         KERR(1) = 'ERROR: RADIUS<0 for AXISYMETRIC PROBLEM'
      ENDIF
      CALL LEREUR
      RETURN
      END
