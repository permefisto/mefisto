      SUBROUTINE GP3BGP1( DeltaT, NDIM,  MNXYZN, NBNOMA, NBTYEL, MNNPEF,
     %                    MNTAUX, NUMIOB, NUMAOB, MNDOEL,
     %                    MNX,    Utn,    Utm, 
     %                    BG,     IERR  )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  GPE: TESTNL=10 ASSEMBLAGE DANS BG DES 4 COEFFICIENTS DES
C ----- 2 POLYNOMES DU 3-EME DEGRE DES EF TRIANGLES P1 LAGRANGE 2D

C -  Beta Vn+1m+1**3 
C -  Beta Wn+1m * Vn+1m+1**2
C + {  Rho/dt -V(X) -Beta Wn+1m**2 +[ Alfa d/dx +OmegaZ y, Alfa d/dy -OmegaZ x] [D] } Vn+1m
C + { -Rho/dt -V(X) -Beta Wn+1m**2 +[ Alfa d/dx -OmegaZ y, Alfa d/dy +OmegaZ x] [D] } Wn+1m
C + Rho/dt (Wn-Vn) - Fr -Fi = 0

C - Beta Wn+1m+1**3 
C + Beta Vn+1m * Wn+1m+1**2
C + {  Rho/dt -V(X) -Beta Vn+1m**2 +[ Alfa d/dx +OmegaZ y, Alfa d/dy -OmegaZ x] [D] } Wn+1m
C - { -Rho/dt -V(X) -Beta Vn+1m**2 +[ Alfa d/dx -OmegaZ y, Alfa d/dy +OmegaZ x] [D] } Vn+1m
C - Rho/dt (Wn+Vn) + Fr -Fi = 0


C ENTREES:
C --------
C DeltaT : PAS DE TEMPS EN DOUBLE PRECISION
C NDIM   : DIMENSION DE L'ESPACE DES COORDONNEES 1 ou 2 ou 3 ou 6
C          (SI AXISYMETRIE NDIM=2 X => R>=0 et Y=>Z et Z=0)
C MNXYZN : ADRESSE MCN DE TMS XYZNOEUD DE L'OBJET
C NBNOMA : NOMBRE TOTAL DE NOEUDS DU MAILLAGE
C
C NBTYEL : NOMBRE DE TYPES D'EF DE L'OBJET
C MNNPEF : ADRESSE MCN DU TMC DES ADRESSES MCN DES TABLEAUX NPEF"TYPE EF
C NDPGST : CODE DE STOCKAGE DES SOMMETS POINTS NOEUDS
C MNTAUX : ADRESSE MCN DES TABLEAUX AUXILAIRES
C NUMIOB : NUMERO MINIMAL DES OBJETS
C NUMAOB : NUMERO MAXIMAL DES OBJETS
C MNDOEL : ADRESSE MCN DES TMS DES DONNEES NLSE DE L'OBJET
C MNX    : ADRESSE MCN DES NBCOOR COORDONNEES DES POINTS DE L'EF COURANT

C Utn    : U(tn) (NBNOMA,2) U A L'INTANT tn
C Utm    : U(tn+1,m) (NBNOMA,2)  U A L'INTANT tn+1 iteration m

C SORTIES:
C --------
C BG     : VECTEUR GLOBAL BG(NBNOMA,2) SECOND MEMBRE
C IERR   : 0 SI PAS D'ERREUR, 1 SI EF NON PROGRAMME, 2 SI EF DEGENERE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR :ALAIN PERRONNET LJLL UPMC & Saint Pierre du Perray   Juin 2014
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
      INTEGER           NONOEF(4), NOOBSF(6), NOOBLA(12), NOOBPS(8)
      INTEGER           NUMIOB(4), NUMAOB(4),  MNDOEL(4)
C
      DOUBLE PRECISION  Utn(NBNOMA,2), Utm(NBNOMA,2), BG(NBNOMA,0:3,2)   
      DOUBLE PRECISION  DeltaT, Rho, Omega(3), Alfa, Beta, Force(2)
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
C     INITIALISATION A ZERO DU VECTEUR GLOBAL BG
C     ------------------------------------------
      CALL AZEROD( NBNOMA*8, BG )
C
C     LA GENERATION DES TABLEAUX ELEMENTAIRES ET DE LA MATRICE
C     ========================================================
C     LE NOMBRE DE POINTS D'INTEGRATION PAR ARETE EN 2D
      NBPTAF = 0
C
C     LA BOUCLE SUR LES TYPES D'ELEMENTS FINIS
      DO 500 NOTYEL = 1, NBTYEL
C
C        L'ADRESSE DU TABLEAU NPEF"TYPE EF
         MNELE = MCN( MNNPEF - 1 + NOTYEL )
C
C        LE NUMERO DU TYPE DE L'ELEMENT FINI
         NUTYEL = MCN( MNELE + WUTYEL )
         CALL ELNUNM( NUTYEL, NOMELE )

C        SEUL L'ELEMENT FINI P1 LAGRANGE EST UTILISE
         IF( NUTYEL .NE. 13 .AND. NUTYEL .NE. 19 ) THEN
C            ELEMENT FINI NON P1 => RETOUR
C            ERREUR
             NBLGRC(NRERR) = 1
             IF( LANGAG .EQ. 0 ) THEN
                KERR(1) = 'ERREUR GP3BGP1: TYPE EF '// NOMELE(1)
     %           // NOMELE(2) //' NON PROGRAMME'
             ELSE
                KERR(1) = 'ERROR GP3BGP1: FE TYPE '// NOMELE(1)
     %           // NOMELE(2) //' NOT PROGRAMMED'
             ENDIF
             CALL LEREUR
             IERR = 1
            GOTO 9999
         ENDIF
C
C        LE NOMBRE D'ELEMENTS FINIS DE CE TYPE P1
         NBELEM =  MCN( MNELE + WBELEM )
C
C        L'ADRESSE MCN DES NOEUDS ET POINTS GEOMETRIQUES DES EF DE CE TYPE
         MNNDEL = MNELE + WUNDEL
         MNPGEL = MNNDEL
C
C        LES CARACTERISTIQUES DE L'ELEMENT FINI
         CALL ELTYCA( NUTYEL )
C
C        TRIANGLE 2P1D  ET  TETRAEDRE 3P1D
         NBPOLY = NDIM + 1
         NPI    = 1
         NPIQ   = 1
         NPIA   = 1
         NBPTAF = MAX( NBPTAF, NPIA )
         MNF1   = MNTAUX
         MNDP   = MNF1 + MOREE2 * NDIM
         MNDFM1 = MNDP + MOREE2 * NDIM * NBPOLY
C        AU TOTAL = MNDFM1 + MOREE2 * NDIM * NDIM

         NBDLEF = NBPOLY * 2
         MNDEF1 = MNDOEL(1)
         MNDEF2 = MNDOEL(2)
         MNDEF3 = MNDOEL(3)
         MNDEF4 = MNDOEL(4)
         NBJEUX = 1
         JEU    = 1
C
C        LES NOEUDS DE L'ELEMENT FINI NUELEM=1
C        -------------------------------------
         NUELEM = 1
         CALL EFNOEU( MNELE, NUELEM, NBNDEL, NONOEF )
C
C        LE NUMERO DE VOLUME  DE L'EF
C        LE NUMERO DE SURFACE DES FACES   DE L'EF
C        LE NUMERO DE LIGNE   DES ARETES  DE L'EF
C        LE NUMERO DE POINT   DES SOMMETS DE L'EF
C        ----------------------------------------
         CALL EFPLSV( MNELE , NUELEM,
     %                NOVCEL, NOSFEL, NOLAEL, NOPSEL,
     %                NOOBVC, NOOBSF, NOOBLA, NOOBPS )
C
C        LES COORDONNEES DES NBPOE POINTS DE L'ELEMENT FINI NUELEM
C        ---------------------------------------------------------
         CALL EFXYZP( NDIM, MNXYZN, NBELEM, NUELEM, MNPGEL, NBPOE,
     %                RMCN(MNX) )

C        RECUPERATION DES COEFFICIENTS SUPPOSES CONSTANTS DU PROBLEME DE
C        GROSS-PITAEVSKII POUR UNE TRIANGULATION P1 (2P1D)
C        SUR LE PREMIER TRIANGLE DU MAILLAGE
C        ---------------------------------------------------------------
         CALL REDOGP( NDIM,   NBPOLY,    RMCN(MNX), NBJEUX, JEU,
     %                NOOBSF, NUMIOB(3), NUMAOB(3), MCN(MNDEF3),
     %                Rho, Omega, Alfa, Beta, Force )

C        ============================================
C        BOUCLE SUR LES ELEMENTS FINIS P1 DU MAILLAGE
C        ============================================
         DO 200 NUELEM = 1, NBELEM
C
C           LES NOEUDS DE L'ELEMENT FINI NUELEM
C           -----------------------------------
            CALL EFNOEU( MNELE, NUELEM, NBNDEL, NONOEF )
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
            CALL EFXYZP( NDIM, MNXYZN, NBELEM, NUELEM, MNPGEL, NBPOE,
     %                   RMCN(MNX) )

C           ************************************
C           2D TRIANGLE TRIA 2P1D LAGRANGE DROIT
C           ************************************
            CALL  GP3BE2P1( DeltaT, Rho, Omega, Alfa, Beta, Force,
     %                      NONOEF, RMCN(MNX), NBNOMA, Utn, Utm,
     %                      BG )

C        FIN DE BOUCLE SUR LES EF DE TYPE NUTYEL
 200     CONTINUE

C     FIN DE BOUCLE SUR LES TYPES D'EF
 500  CONTINUE

C     AFFICHAGE DES PREMIERS POLYNOMES
 9999 DO K=1,NBNOMA
         DO J=1,2

            IF( ABS(BG(K,0,J)) .GT. 0D0 .AND.
     %          ABS(BG(K,0,J)) .LT. 1D-28 ) THEN
ccc               PRINT 19999,(K,I,J,BG(K,I,J),I=3,0,-1)
               BG(K,0,J) = 0D0
            ENDIF

            IF(  ABS(BG(K,2,J)) .GT. 0D0 .AND.
     %           ABS(BG(K,2,J)) .LT. 1D-28 ) THEN
ccc               PRINT 19999,(K,I,J,BG(K,I,J),I=3,0,-1)
               BG(K,2,J) = 0D0
            ENDIF

         ENDDO
      ENDDO

ccc19999 FORMAT( 4('  BG(',I6,',',I1,',',I1,')=',G15.6) )
      RETURN
      END
