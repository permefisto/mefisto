      SUBROUTINE THEMKB( NBJEUX, IEMG,   IEKG,   IEBG,   PENALI,
     %                   D2PI,   NDIM,   NTDL,   VITEGt,
     %                   NBTYEL, MNNPEF, NDPGST,
     %                   MNTPOB, MXPLBA, MNTAUX,
     %                   MNXYZP, NUMIOB, NUMAOB, MNDOEL,
     %                   MNTHER, MNTAEL, MNX,    MNNODL,
     %                   NORESO, MNLPLI, MNLPCO,
     %                   NBRDMG, MG,     NBRDKG, KG,   MNBG,
     %                   NCODSM, NCODSK, NBPTAF,
     %                   IERR  )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    EN THERMIQUE CALCUL DE LA MATRICE GLOBALE DE CAPACITE
C -----    ET/OU CONDUCTIVITE ET/OU DU SECOND MEMBRE GLOBAL
C          DES ELEMENTS FINIS LAGRANGE 1D OU 2D OU 3D DE DEGRE 1 OU 2
C
C          HEAT TRANSFER PROBLEM:
C          COMPUTATION OF THE CAPACITY GLOBAL MATRIX,
C                         THE CONDUCTIVITY GLOBAL MATRIX,
C                         THE SOURCE GLOBAL VECTOR
C          FROM ISOPARAMETRIC LAGRANGE FINITE ELEMENT 1D or 2D or 3D
C          WITH POLYNOMIALS OF DEGREE 1 or 2 ON A REFERENCE FINITE ELEMENT
C
C ENTREES:
C --------
C NBJEUX : NOMBRE DE JEUX DE DONNEES = DEGRE+1 DU POLYNOME P(Lambda)
C          >1 ICI POUR AVOIR UN POLYNOME DE DEGRE>=2
C IEMG   : =1 SI CALCUL DEMANDE DE LA MATRICE GLOBALE DE CAPACITE
C          =0 SI PAS DE CALCUL  DE LA MATRICE GLOBALE DE CAPACITE
C IEKG   : =1 SI CALCUL DEMANDE DE LA MATRICE GLOBALE DE CONDUCTIVITE
C          =0 SI PAS DE CALCUL  DE LA MATRICE GLOBALE DE CONDUCTIVITE
C IEBG   : =1 SI CALCUL DEMANDE DU SECOND MEMBRE GLOBAL
C          =0 SI PAS DE CALCUL  DU SECOND MEMBRE GLOBAL
C
C PENALI : COEFFICIENT DE PENALISATION DE LA CONDITION DE DIRICHLET
C          >0D0 SI LE TMS "CONTACT" EST A PENALISER ET PENALI DOIT ETRE GRAND
C          =0D0 PAS DE PENALISATION DU TMS "CONTACT"
C
C D2PI   : =2 PI EN REEL DOUBLE PRECISION
C NDIM   : DIMENSION DE L'ESPACE DES COORDONNEES 1 ou 2 ou 3 ou 6
C          (SI AXISYMETRIE NDIM=2 X => R>=0 et Y=>Z et Z=0)
C NTDL   : NOMBRE TOTAL DE DEGRES DE LIBERTE DU MAILLAGE DE L'OBJET
C VITEGt : VECTEUR(NTDL,NDIM) DES NDIM COMPOSANTES DE LA VITESSE
C          AUX NOEUDS DU MAILLAGE
C          CE TABLEAU EST UTILISE SEULEMENT EN CAS d'UN TRANSPORT
C          DE TEMPERATURE A CETTE VITESSE POUR LE SECOND MEMBRE
C          terme: -(VITEGt . Grad) Temperature

C NBTYEL : NOMBRE DE TYPES D'EF DE L'OBJET
C MNNPEF : ADRESSE MCN DU TMC DES ADRESSES MCN DES TABLEAUX NPEF"TYPE EF
C NDPGST : CODE DE STOCKAGE DES SOMMETS POINTS NOEUDS
C MNTPOB : ADRESSE MCN DES TABLEAUX POLYNOMES DE BASE DES TYPES D'EF
C MXPLBA : NOMBRE MAXIMAL DE TABLEAUX POBA PAR TYPE D'EF
C MNTAUX : ADRESSE MCN DES TABLEAUX AUXILAIRES
C MNXYZP : ADRESSE MCN DE TMS XYZSOMMET DE L'OBJET
C NUMIOB : NUMERO MINIMAL DES OBJETS
C NUMAOB : NUMERO MAXIMAL DES OBJETS
C MNDOEL : ADRESSE MCN DES TMS DES DONNEES THERMIQUES DE L'OBJET
C MNTHER : ADRESSE MCN DU TABLEAU DES DONNEES THERMIQUES INTERNES
C          AU MOINS NPIMAX(=128) REELS DOUBLE PRECISION
C MNTAEL : ADRESSE MCN DES TABLEAUX ELEMENTAIRES
C MNX    : ADRESSE MCN DES NBCOOR COORDONNEES DES POINTS DE L'EF COURANT
C MNNODL : ADRESSE MCN DU TABLEAU DES NUMEROS DES NOEUDS DE L'EF COURANT
C
C NORESO : 1 STOCKAGE PROFIL DES MATRICES
C          2 STOCKAGE MORSE  DES MATRICES
C          3 MULTI-PROCESSORS SOLVER of A x = b FROM A MORSE MATRIX A
C MNLPLI : ADRESSE MCN DU TABLEAU POINTEUR DES COEFFICIENTS DIAGONAUX
C          DE LA MATRICE PROFIL OU MORSE DE K
C          SI EN SORTIE NCODSM/=0 ALORS M A MEME PROFIL QUE K
C                                 SINON M EST DIAGONALE
C MNLPCO : ADRESSE MCN DU TABLEAU DES NUMEROS DE COLONNES DE LA MATRICE MORSE
C
C NBRDMG : NOMBRE DE REELS DOUBLE PRECISION DE LA MATRICE GLOBALE MG
C          = NBRDKG SI MATRICE NON DIAGONALE, =NTDL SI MATRICE DIAGONALE
C MG     : NBJEUX MATRICES GLOBALES DE CAPACITE
C          0 SI IEMG=0
C NBRDKG : NOMBRE DE REELS DOUBLE PRECISION DE LA MATRICE GLOBALE MG et KG
C KG     : NBJEUX MATRICES GLOBALES DE MASSE & CONDUCTIVITE
C          0 SI IEKG=0
C MNBG   : ADRESSE MCN DU TABLEAU DU VECTEUR GLOBAL SECOND MEMBRE
C          0 SI IEBG=0
C
C NCODSM : CODE DE STOCKAGE DE LA MATRICE GLOBALE DE CAPACITE
C          0 SI MATRICE DIAGONALE
C          1 SI MATRICE SYMETRIQUE
C         -1 SI MATRICE NON SYMETRIQUE
C         NON INITIALISE SI IEMG=0
C NCODSK : CODE DE STOCKAGE DE LA MATRICE GLOBALE DE CONDUCTIVITE
C          1 SI MATRICE SYMETRIQUE
C         -1 SI MATRICE NON SYMETRIQUE
C         NON INITIALISE SI IEKG=0
C
C SORTIE :
C --------
C NBPTAF : NOMBRE DE POINTS PAR ARETE OU LE FLUX EST CALCULE
C IERR   : 0 SI PAS D'ERREUR, 1 SI EF NON PROGRAMME, 2 SI EF DEGENERE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS        Aout 1998
C AUTEUR : ALAIN PERRONNET Texas A & M University           Juillet 2005
C AUTEUR : ALAIN PERRONNET TIMS NTU TAIPEI TAIWAN          Novembre 2009
C MODIFS : ALAIN PERRONNET Saint Pierre du Perray             Avril 2022
C23456---------------------------------------------------------------012
C$    use OMP_LIB
      include"./incl/langue.inc"
      include"./incl/donthe.inc"
      include"./incl/a___npef.inc"
      include"./incl/a___xyzpoint.inc"
      include"./incl/a___vecteur.inc"
      include"./incl/ponoel.inc"
      include"./incl/gsmenu.inc"
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
      DOUBLE PRECISION  MG(NBRDMG,NBJEUX), KG(NBRDKG,NBJEUX),
     %                  VITEGt(NTDL,NDIM)
C
      DOUBLE PRECISION  PENALI, D2PI, DELTA
      CHARACTER*4       NOMELE(2)
C
C     LE NOMBRE DE MOTS D'UNE VARIABLE REEL2
      MOREE2 = MOTVAR(6)
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
C     INITIALISATION DE LA MATRICE DE CAPACITE
      IF( IEMG .NE. 0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10001)
         ELSE
            WRITE(IMPRIM,20001)
         ENDIF
10001    FORMAT('THEMKB: CONSTRUCTION de la MATRICE GLOBALE [MG]')
20001    FORMAT('THEMKB: CONSTRUCTION of the GLOBAL MATRIX [MG]')
      ENDIF
C     INITIALISATION DE LA MATRICE DE CONDUCTIVITE
      IF( IEKG .NE. 0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10002)
         ELSE
            WRITE(IMPRIM,20002)
         ENDIF
10002    FORMAT('THEMKB: CONSTRUCTION de la MATRICE GLOBALE [KG]')
20002    FORMAT('THEMKB: CONSTRUCTION of the GLOBAL MATRIX [KG]')
      ENDIF
C     INITIALISATION DU SECOND MEMBRE
      IF( IEBG .NE. 0 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10003)
         ELSE
            WRITE(IMPRIM,20003)
         ENDIF
10003    FORMAT('THEMKB: CONSTRUCTION du VECTEUR {BG} SECOND MEMBRE')
20003    FORMAT('THEMKB: CONSTRUCTION of the SECOND MEMBER VECTOR {BG}')
      ENDIF
C
C     INITIALISATION A ZERO DES MATRICES GLOBALES ET VECTEURS GLOBAUX
C     A CALCULER POUR LES NBJEUX
C     ---------------------------------------------------------------
      MOTBG = NTDL * MOREE2
      MNBGJ = MNBG - MOTBG
      DO JEU = 1, NBJEUX
         IF( IEMG .GT. 0 ) CALL AZEROD( NBRDMG, MG(1,JEU)  )
         IF( IEKG .GT. 0 ) CALL AZEROD( NBRDKG, KG(1,JEU)  )
         IF( IEBG .GT. 0 ) THEN
            MNBGJ = MNBGJ + MOTBG
            CALL AZEROD( NTDL, DMCN( (MNBGJ+1)/MOREE2 ) )
         ENDIF
      ENDDO
C
C     LA GENERATION DES TABLEAUX ELEMENTAIRES ET LES MATRICES DE
C     CAPACITE ET CONDUCTIVITE SONT STOCKEES SOUS FORME PROFIL
C     ==========================================================
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
         NBELEM =  MCN( MNELE + WBELEM )
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
            KERR(1) = 'ERREUR THEMKB: TYPE EF '// NOMELE(1)
     %           // NOMELE(2) //' NON PROGRAMME'
         ELSE
            KERR(1) = 'ERROR THEMKB: FE TYPE '// NOMELE(1)
     %           // NOMELE(2) //' NOT PROGRAMMED'
         ENDIF
         CALL LEREUR
         IERR = 1
         GOTO 9999
C
C        RECHERCHE DU TABLEAU DE POBA ET PARTAGE EN P ET DP
 10      L = MNTPOB + (NOTYEL-1) * MXPLBA
C
         IF( NUTYEL .EQ. 30 ) THEN
C           6CUBE 6Q1C N'A POUR L'INSTANT PAS DE FACE...
            NBNSOM = 0
            NARET  = 0
            NFACE  = 0
            GOTO 12
         ENDIF
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
ccc        MNPDEL = MNF1 + MOREE2 * NPI * NDIM + MOREE2 * NPI * (NDIM-2) modif l
         MNPDEL = MNF1   + MOREE2 * NPI * NDIM
         MNDP   = MNPDEL + MOREE2 * NPI
         MNDFM1 = MNDP   + MOREE2 * NPI * NDIM * NBPOLY
C        AU TOTAL = MNDFM1 + MOREE2 * NPI * NDIM * NDIM
C        MOTAUX = MOREE2 * NPI * ( NDIM + 1 + NDIM * (NDIM+NBPOLY) )
C
C        ADRESSE DU SECOND MEMBRE ELEMENTAIRE
C        DERRIERE LA MATRICE SYMETRIQUE ELEMENTAIRE DE CONDUCTIVITE
         NBDL = NBPOLY
         GOTO 30
C
C        TRIANGLE 2P1D  ET  TETRAEDRE 3P1D
C        =================================
 13      NBPOLY = NDIM + 1
         NBDL   = NBPOLY
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
 30      IF( IEMG .NE. 0 ) THEN
            MNCAEL = MNTAEL + NBDL * (NBDL+1) / 2  * MOREE2
         ELSE
            MNCAEL = MNTAEL
         ENDIF
         MNSE = MNCAEL + NBDL * (NBDL+1) / 2  * MOREE2
C
C        LA BOUCLE SUR LES ELEMENTS FINIS DE CE TYPE NUTYEL
C        ==================================================
         MNDEF1 = MNDOEL(1)
         MNDEF2 = MNDOEL(2)
         MNDEF3 = MNDOEL(3)
         MNDEF4 = MNDOEL(4)
C
         DO 200 NUELEM = 1, NBELEM
C
C           LES NOEUDS DE L'ELEMENT FINI NUELEM
C           -----------------------------------
            CALL EFNOEU( MNELE, NUELEM, NBNDEL, MCN(MNNODL) )
C
C           LES POINTS GEOMETRIQUES DE L'ELEMENT
C           ------------------------------------
CCC            CALL EFPOGE( MNELE, NUELEM, NBPGEF, MCN(MNPOEF) )
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
C           LES COORDONNEES DES NBPOE POINTS DE L'ELEMENT FINI NUELEM
C           ---------------------------------------------------------
            CALL EFXYZP( NDIM, MNXYZP, NBELEM, NUELEM, MNPGEL, NBPOE,
     %                   RMCN(MNX) )
C
C           BOUCLE SUR LES JEUX DE DONNEES THERMIQUES
C           =========================================
            MNBGJ = MNBG - MOTBG
C
            DO 100 JEU = 1, NBJEUX
C
C              ===============================================================
C              LE CALCUL DES TABLEAUX AUXILIAIRES ET DES TABLEAUX ELEMENTAIRES
C              ===============================================================
               GOTO( 41, 41, 41, 41,  1,  1, 1,  1,  1,  1,
     %                1,  1, 42,  1, 41, 41, 1, 41, 49, 51,
     %               51, 51, 51, 51,  1,  1, 1, 40, 41, 61,
     %               51, 51, 40,  1  ), NUTYEL
C
C           *********************************************
C           1D LAGRANGE F:SEGMENT UNITE -> SEGMENT EST P1
C           *********************************************
C           LE CALCUL DES TABLEAUX AUXILIAIRES EN 1D
C           ----------------------------------------
 40         IF( RMCN(MNX) .GT. RMCN(MNX+1) ) THEN
C              PERMUTATION DES 2 SOMMETS POUR AVOIR UN EF DE MESURE>0
               XS = RMCN(MNX)
               RMCN(MNX) = RMCN(MNX+1)
               RMCN(MNX+1 ) = XS
C
               NST               = MCN( MNNODL )
               MCN( MNNODL     ) = MCN( MNNODL + 1 )
               MCN( MNNODL + 1 ) = NST
C
               NST       = NOOBPS(1)
               NOOBPS(1) = NOOBPS(2)
               NOOBPS(2) = NST
            ENDIF
            CALL E11LAG( NBPOLY, NPI, MCN(MNPOID),
     %                   MCN(MNPOL),  MCN(MNDPOL),
     %                   RMCN(MNX) ,  MCN(MNF1),
     %                   MCN(MNPDEL), MCN(MNDP), DELTA )
C
C           LA MATRICE ELEMENTAIRE DE CAPACITE EN 1D
C           ----------------------------------------
            IF( IEMG .NE. 0 ) THEN
C              MATRICE DE CAPACITE SYMETRIQUE NON DIAGONALE
               CALL TM1LAG( NBJEUX, JEU, NBPOLY, NPI, MCN(MNPOL),
     %                      NOOBLA(1), NUMIOB(2), NUMAOB(2),
     %                      MCN(MNDEF2),
     %                      MCN(MNF1), MCN(MNPDEL),
     %                      MCN(MNTHER), NCODEM, MCN(MNCAEL))
            ENDIF
C
C           LA MATRICE ELEMENTAIRE DE CONDUCTIVITE EN 1D
C           --------------------------------------------
            IF( IEKG .NE. 0 ) THEN
               CALL TR1LAG( RMCN(MNX), PENALI, NBJEUX, JEU,
     %                      NOOBPS,
     %                      NUMIOB(1), NUMAOB(1), MCN(MNDEF1),
     %                      NOOBLA,
     %                      NUMIOB(2), NUMAOB(2), MCN(MNDEF2),
     %                      NBPOLY, NPI, MCN(MNPOL),
     %                      MCN(MNF1),   MCN(MNPDEL), MCN(MNDP),
     %                      MCN(MNTHER), MCN(MNTHER), MCN(MNTAEL) )
               IF( IERR .NE. 0 ) GOTO 9999
            ENDIF
C
C           LES SECONDS MEMBRES ELEMENTAIRES EN 1D
C           --------------------------------------
            IF( IEBG .NE. 0 ) THEN
               CALL TS1LAG( RMCN(MNX), PENALI, NBJEUX, JEU,
     %                      NOOBPS,NUMIOB(1),NUMAOB(1),MCN(MNDEF1),
     %                      NOOBLA,NUMIOB(2),NUMAOB(2),MCN(MNDEF2),
     %                      NBPOLY,NPI,MCN(MNPOL),
     %                      MCN(MNF1),MCN(MNPDEL),MCN(MNDP),
     %                      0, MCN(MNSE) )
            ENDIF
            GOTO 70
C
C           ********************************************
C           2D OU AXISYMETRIQUE LAGRANGE ISOPARAMETRIQUE
C           ********************************************
C           LE CALCUL DES TABLEAUX AUXILIAIRES EN 2D OU AXISYMETRIE
C           -------------------------------------------------------
 41         CALL E12LAG( NBPOLY,NPI,MCN(MNPOID),
     %                  MCN(MNPOL), MCN(MNDPOL),
     %                  RMCN(MNX) , MCN(MNF1), MCN(MNF2),
     %                  MCN(MNPDEL),MCN(MNDP),MCN(MNDFM1) )
C
C           LA MATRICE ELEMENTAIRE DE CAPACITE EN 2D OU EN AXISYMETRIE
C           ----------------------------------------------------------
            IF( IEMG .NE. 0 ) THEN
CCCC
CCC               IF( NUTYEL .EQ. 16 ) THEN
CCCC                 LA MATRICE ELEMENTAIRE DE CAPACITE EN 2D
CCCC                 QUAD 2Q1C => MATRICE DE CAPACITE DIAGONALE
CCC                  CALL TM2Q1C( RMCN(MNX), NBJEUX, JEU,
CCC     %                         NOOBSF(1),NUMIOB(3),NUMAOB(3),
CCC     %                         MCN(MNDEF3),
CCC     %                         NCODEM, MCN(MNCAEL) )
CCCC
CCC               ELSE
CCC
C                 EF 2D ET AXISYMETRIE SAUF TRIA 2P1D ET QUAD 2Q1C
C                 MATRICE DE CAPACITE SYMETRIQUE NON DIAGONALE
                  CALL TM2LAG( D2PI, NOAXIS, NBJEUX, JEU,
     %                         NBPOLY,NPI,MCN(MNPOL),
     %                         NOOBSF(1),NUMIOB(3),NUMAOB(3),
     %                         MCN(MNDEF3),
     %                         MCN(MNF1),MCN(MNF2),MCN(MNPDEL),
     %                         MCN(MNTHER),NCODEM,MCN(MNCAEL))
CCC               ENDIF
            ENDIF
C
C           LA MATRICE ELEMENTAIRE DE CONDUCTIVITE EN 2D OU EN AXISYMETRIE
C           --------------------------------------------------------------
            IF( IEKG .NE. 0 ) THEN
               CALL TR2LAG( D2PI, NOAXIS, RMCN(MNX), PENALI, NBJEUX,JEU,
     %                      NBNSOM, NOOBPS,
     %                      NUMIOB(1), NUMAOB(1), MCN(MNDEF1),
     %                      NBPOLA, NPIA,
     %                      MCN(MNPOIA),MCN(MNPOLA),MCN(MNDPOA),
     %                      NARET,NOOBLA,
     %                      NUMIOB(2),NUMAOB(2),MCN(MNDEF2),
     %                      NBPOLY,NPI,MCN(MNPOL),
     %                      NOOBSF,NUMIOB(3),NUMAOB(3),MCN(MNDEF3),
     %                      MCN(MNF1),MCN(MNF2),MCN(MNPDEL),MCN(MNDP),
     %                      MCN(MNTHER),MCN(MNTHER),MCN(MNTAEL),IERR)
               IF( IERR .NE. 0 ) GOTO 9999
            ENDIF
C
C           LES SECONDS MEMBRES ELEMENTAIRES EN 2D OU AXISYMETRIE
C           -----------------------------------------------------
            IF( IEBG .NE. 0 ) THEN
               CALL TS2LAG( D2PI, NOAXIS, RMCN(MNX), PENALI, NBJEUX,JEU,
     %                      NBNSOM, NOOBPS,
     %                      NUMIOB(1),NUMAOB(1),MCN(MNDEF1),
     %                      NBPOLA,NPIA,
     %                      MCN(MNPOIA),MCN(MNPOLA),MCN(MNDPOA),
     %                      NARET,NOOBLA,
     %                      NUMIOB(2),NUMAOB(2),MCN(MNDEF2),
     %                      NBPOLY,NPI,MCN(MNPOL),
     %                      NOOBSF,NUMIOB(3),NUMAOB(3),MCN(MNDEF3),
     %                      NTDL, VITEGt,
     %                      MCN(MNF1),MCN(MNF2),MCN(MNPDEL),MCN(MNDP),
     %                      0, MCN(MNSE) )
            ENDIF
            GOTO 70
C
C           ***************************
C           3D LAGRANGE ISOPARAMETRIQUE
C           ***************************
C           LE CALCUL DES TABLEAUX AUXILIAIRES EN 3D
C           ----------------------------------------
 51         CALL E13LAG ( NBPOLY, NPI, MCN(MNPOID),
     %                    MCN(MNPOL),  MCN(MNDPOL),
     %                    RMCN(MNX) ,  MCN(MNF1),
     %                    MCN(MNPDEL), MCN(MNDP), MCN(MNDFM1) )
C
C           LA MATRICE ELEMENTAIRE DE CAPACITE EN 3D
C           ----------------------------------------
            IF( IEMG .NE. 0 ) THEN
              CALL TM3LAG( NBJEUX, JEU,
     %                     NBPOLY, NPI, MCN(MNPOL),
     %                     NOOBVC, NUMIOB(4), NUMAOB(4),MCN(MNDEF4),
     %                     MCN(MNF1),
     %                     MCN(MNPDEL),
     %                     MCN(MNTHER), NCODEM, MCN(MNCAEL) )
            ENDIF
C
C           LA MATRICE ELEMENTAIRE DE CONDUCTIVITE EN 3D
C           --------------------------------------------
            IF( IEKG .NE. 0 ) THEN
               CALL TR3LAG( NUTYEL, RMCN(MNX), PENALI, NBJEUX, JEU,
     %                      NBNSOM, NOOBPS,
     %                      NUMIOB(1),NUMAOB(1),MCN(MNDEF1),
     %                      NARET,NOOBLA,
     %                      NUMIOB(2),NUMAOB(2),MCN(MNDEF2),
     %                      NBPOLA, NPIA, MCN(MNPOIA),
     %                      MCN(MNPOLA),  MCN(MNDPOA),
     %                      NBPOLQ, NPIQ, MCN(MNPOIQ),
     %                      MCN(MNPOLQ),  MCN(MNDPOQ),
     %                      NFACE, NOOBSF,
     %                      NUMIOB(3), NUMAOB(3), MCN(MNDEF3),
     %                      NBPOLY, NPI, MCN(MNPOL),
     %                      NOOBVC,NUMIOB(4),NUMAOB(4),MCN(MNDEF4),
     %                      MCN(MNF1),   MCN(MNPDEL), MCN(MNDP),
     %                      MCN(MNTHER), MCN(MNTHER), MCN(MNTAEL) )
            ENDIF
C
C           LES SECONDS MEMBRES ELEMENTAIRES EN 3D
C           --------------------------------------
            IF( IEBG .NE. 0 ) THEN
               CALL TS3LAG( NUTYEL, RMCN(MNX), PENALI, NBJEUX, JEU,
     %                      NBNSOM, NOOBPS,
     %                      NUMIOB(1),NUMAOB(1),MCN(MNDEF1),
     %                      NARET,NOOBLA,
     %                      NUMIOB(2),NUMAOB(2),MCN(MNDEF2),
     %                      NBPOLA, NPIA, MCN(MNPOIA),
     %                      MCN(MNPOLA), MCN(MNDPOA),
     %                      NBPOLQ, NPIQ, MCN(MNPOIQ),
     %                      MCN(MNPOLQ), MCN(MNDPOQ),
     %                      NFACE, NOOBSF,
     %                      NUMIOB(3), NUMAOB(3), MCN(MNDEF3),
     %                      NBPOLY, NPI, MCN(MNPOL),
     %                      NOOBVC, NUMIOB(4),NUMAOB(4),MCN(MNDEF4),
     %                      NTDL, VITEGt,
     %                      MCN(MNF1), MCN(MNPDEL), MCN(MNDP),
     %                      0, MCN(MNSE) )
            ENDIF
            GOTO 70
C
C           ************************************
C           2D TRIANGLE TRIA 2P1D LAGRANGE DROIT
C           ************************************
C           LA MATRICE DIAGONALE ELEMENTAIRE DE CAPACITE
C           --------------------------------------------
 42         IF( IEMG .NE. 0 ) THEN
               CALL TM2P1D( RMCN(MNX), NBJEUX, JEU,
     %                      NOOBSF(1),NUMIOB(3),NUMAOB(3),
     %                      MCN(MNDEF3),NCODEM,MCN(MNCAEL) )
CCCC              TRANSFORMATION MATRICE DIAGONALE EN MATRICE SYMETRIQUE SUR ELL
CCC               CALL DIASYM( 3, MCN(MNCAEL), MCN(MNCAEL) )
CCC               NCODEM = 1
            ENDIF
C
C           LA MATRICE ELEMENTAIRE DE CONDUCTIVITE EN 2D
C           --------------------------------------------
            IF( IEKG .NE. 0 ) THEN
               CALL TR2P1D( RMCN(MNX),PENALI, NBJEUX, JEU,
     %                      NOOBPS,NUMIOB(1),NUMAOB(1),MCN(MNDEF1),
     %                      NOOBLA,NUMIOB(2),NUMAOB(2),MCN(MNDEF2),
     %                      NOOBSF,NUMIOB(3),NUMAOB(3),MCN(MNDEF3),
     %                      MCN(MNTHER), MCN(MNTAEL) )
            ENDIF
C
C           LES SECONDS MEMBRES ELEMENTAIRES
C           --------------------------------
            IF( IEBG .NE. 0 ) THEN
               CALL TS2P1D( RMCN(MNX),PENALI, NBJEUX, JEU,
     %                      NOOBPS,NUMIOB(1),NUMAOB(1),MCN(MNDEF1),
     %                      NOOBLA,NUMIOB(2),NUMAOB(2),MCN(MNDEF2),
     %                      NOOBSF,NUMIOB(3),NUMAOB(3),MCN(MNDEF3),
     %                      0, MCN(MNSE) )
            ENDIF
            GOTO 70
C
C           *************************************
C           3D TETRAEDRE TETR 3P1D LAGRANGE DROIT
C           *************************************
 49         CALL  E13P1D( RMCN(MNX), MCN(MNF1),
     %                    DELTA, MCN(MNDFM1), MCN(MNDP) )
C
C           LA MATRICE DIAGONALE ELEMENTAIRE DE CAPACITE
C           RANGEE DANS UNE MATRICE SYMETRIQUE
C           --------------------------------------------
            IF( IEMG .NE. 0 ) THEN
C              ATTENTION: STOCKAGE SYMETRIQUE ET NON DIAGONAL
C                         POUR COMPATIBILITE AVEC PENTAEDRE ET HEXAEDRE DEGRE 1
C                         DE MATRICE DE CAPACITE NON DIAGONALE
               CALL TM3P1D( RMCN(MNX), DELTA, NBJEUX, JEU,
     %                      NOOBVC, NUMIOB(4),NUMAOB(4),MCN(MNDEF4),
     %                      NCODEM, MCN(MNCAEL) )
CCCC              TRANSFORMATION MATRICE DIAGONALE EN MATRICE SYMETRIQUE SUR ELL
CCC               CALL DIASYM( 4, MCN(MNCAEL), MCN(MNCAEL) )
CCC               NCODEM = 1
            ENDIF
C
C           LA MATRICE ELEMENTAIRE DE CONDUCTIVITE
C           --------------------------------------
            IF( IEKG .NE. 0 ) THEN
               CALL TR3P1D( RMCN(MNX), DELTA,  MCN(MNDP), PENALI,
     %                      NBJEUX, JEU,
     %                      NOOBPS,NUMIOB(1),NUMAOB(1),MCN(MNDEF1),
     %                      NOOBLA,NUMIOB(2),NUMAOB(2),MCN(MNDEF2),
     %                      NOOBSF,NUMIOB(3),NUMAOB(3),MCN(MNDEF3),
     %                      NOOBVC,NUMIOB(4),NUMAOB(4),MCN(MNDEF4),
     %                      MCN(MNTHER), MCN(MNTAEL) )
            ENDIF
C
C           LES SECONDS MEMBRES ELEMENTAIRES
C           --------------------------------
            IF( IEBG .NE. 0 ) THEN
               CALL TS3P1D( RMCN(MNX), DELTA,  MCN(MNDP), PENALI,
     %                      NBJEUX, JEU,
     %                      NOOBPS,NUMIOB(1),NUMAOB(1),MCN(MNDEF1),
     %                      NOOBLA,NUMIOB(2),NUMAOB(2),MCN(MNDEF2),
     %                      NOOBSF,NUMIOB(3),NUMAOB(3),MCN(MNDEF3),
     %                      NOOBVC,NUMIOB(4),NUMAOB(4),MCN(MNDEF4),
     %                      0, MCN(MNSE) )
            ENDIF
            GOTO 70
C
C           ***************************************
C           6D LAGRANGE ISOPARAMETRIQUE  6CUBE 6Q1C
C           ***************************************
C           LE CALCUL DES TABLEAUX AUXILIAIRES EN 6D
C           ----------------------------------------
 61         CALL E16Q1C( NBPOLY, NPI, MCN(MNPOID),
     %                   MCN(MNPOL),  MCN(MNDPOL),
     %                   RMCN(MNX) ,  MCN(MNF1),
     %                   MCN(MNPDEL), MCN(MNDP), MCN(MNDFM1), IERR )
C           SI EF DEGENERE => RETOUR
            IF( IERR .NE. 0 ) RETURN
C
C           LA MATRICE ELEMENTAIRE DE CAPACITE EN 6D
C           ----------------------------------------
            IF( IEMG .NE. 0 ) THEN
              CALL TM6Q1C( NBJEUX, JEU, NBPOLY, NPI, MCN(MNPOL),
     %                     NOOBVC, NUMIOB(4), NUMAOB(4),MCN(MNDEF4),
     %                     MCN(MNF1), MCN(MNPDEL),
     %                     MCN(MNTHER), NCODEM, MCN(MNCAEL) )
            ENDIF
C
C           LA MATRICE ELEMENTAIRE DE CONDUCTIVITE EN 6D
C           --------------------------------------------
            IF( IEKG .NE. 0 ) THEN
               CALL TR6Q1C( NBJEUX, JEU, NBPOLY, NPI, MCN(MNPOL),
     %                      NOOBVC,NUMIOB(4),NUMAOB(4),MCN(MNDEF4),
     %                      MCN(MNF1),   MCN(MNPDEL), MCN(MNDP),
     %                      MCN(MNTHER), MCN(MNTHER), MCN(MNTAEL) )
            ENDIF
C
C           LES SECONDS MEMBRES ELEMENTAIRES EN 6D
C           --------------------------------------
            IF( IEBG .NE. 0 ) THEN
               CALL TS6Q1C( NBJEUX, JEU, NBPOLY, NPI, MCN(MNPOL),
     %                      NOOBVC, NUMIOB(4),NUMAOB(4),MCN(MNDEF4),
     %                      MCN(MNF1), MCN(MNPDEL), MCN(MNDP),
     %                      MCN(MNSE) )
            ENDIF
C            GOTO 70
C
C           ASSEMBLAGE DES TABLEAUX ELEMENTAIRES DANS LES TABLEAUX GLOBAUX
C           ==============================================================
 70         NBDL = NBNOE
C
            IF ( NORESO .EQ. 1 ) THEN
C
C              STOCKAGE PROFIL DES MATRICES
C              ----------------------------
               IF( IEMG .NE. 0 ) THEN
C                 ASSEMBLAGE MATRICE PROFIL DE CAPACITE
                  CALL ASMEPC( NBDL,   MCN(MNNODL),
     %                         NCODEM, MCN(MNCAEL), MCN(MNCAEL),
     %                         NCODSM, MCN(MNLPLI), MG(1,JEU)  )
               ENDIF
C
               IF( IEKG .NE. 0 ) THEN
C                 ASSEMBLAGE MATRICE PROFIL DE CONDUCTIVITE
                  CALL ASMEPC( NBDL,   MCN(MNNODL),
     %                         NCODSK, MCN(MNTAEL), MCN(MNTAEL),
     %                         NCODSK, MCN(MNLPLI), KG(1,JEU)  )
               ENDIF
C
            ELSE IF( NORESO .GE. 2 ) THEN
C
C              STOCKAGE MORSE DES MATRICES POUR LE GC
C              --------------------------------------
               IF( IEMG .NE. 0 ) THEN
C                 ASSEMBLAGE MATRICE MORSE DE CAPACITE
                  CALL ASMEGC( NBDL,  MCN(MNNODL),
     %                         NCODEM,MCN(MNCAEL),MCN(MNCAEL),
     %                        NCODSM,MCN(MNLPLI),MCN(MNLPCO),MG(1,JEU) )
               ENDIF
C
               IF( IEKG .NE. 0 ) THEN
C                 ASSEMBLAGE MATRICE MORSE DE CONDUCTIVITE
                  CALL ASMEGC( NBDL,  MCN(MNNODL),
     %                         NCODSK,MCN(MNTAEL),MCN(MNTAEL),
     %                        NCODSK,MCN(MNLPLI),MCN(MNLPCO),KG(1,JEU) )
               ENDIF
            ENDIF
C
            IF( IEBG .NE. 0 ) THEN
C              ASSEMBLAGE SECOND MEMBRE DANS LE SECOND MEMBRE GLOBAL DU JEU
               MNBGJ = MNBGJ + MOTBG
               CALL AS1BEBG( NBDL, MCN(MNNODL), MCN(MNSE),  MCN(MNBGJ) )
            ENDIF
C
C        FIN BOUCLE SUR LES JEUX
 100     CONTINUE
C
C     FIN DE BOUCLE SUR LES EF DE TYPE NUTYEL
 200  CONTINUE
C
C     FIN DE BOUCLE SUR LES TYPES D'EF
 500  CONTINUE
C
ccc      GOTO 9999
C
cccC     ERREUR: EF DEGENERE RENCONTRE
ccc 9900 NBLGRC(NRERR) = 1
ccc      IF( LANGAG .EQ. 0 ) THEN
ccc         KERR(1) = 'ERREUR THEMKB: un EF '// NOMELE(1)
ccc     %        // NOMELE(2) // ' est DEGENERE'
ccc      ELSE
ccc         KERR(1) = 'ERROR THEMKB: 1 FE '// NOMELE(1)
ccc     %        // NOMELE(2) // ' is DEGENERATED'
ccc      ENDIF
ccc      CALL LEREUR
ccc      IERR = 2
C
C     FIN DE LA BOUCLE SUR LES TYPES D'EF
 9999 RETURN
      END
