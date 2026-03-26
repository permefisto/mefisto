      SUBROUTINE THNORM( EXPO,   NDIM,   NTDL,   NDSM,  TEMPER,
     %                   NBCOOR, MNX,    MNNODL, NBTYEL, MNNPEF, NDPGST,
     %                   MNTPOB, MXPOBA, MNTAUX, MNXYZP,
     %                   IERR  )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LA NORME L1 DES NDSM TEMPERATURES**EXPO
C -----
C          normp = Som       Som    Omegal  (û**EXPO)(bl) Delta(bl)
C                  e dans E l=1...,L
C          Vol   = Som       Som    Omegal Delta(bl)
C                  e dans E l=1...,L
C
C EXPO   : EXPOSANT DE LA SOLUTION DE NORME L1 A CALCULER
C NDIM   : DIMENSION DE L'ESPACE DES COORDONNEES 2 ou 3 ou 6
C          (SI AXISYMETRIE NDIM=2 X => R>=0 et Y=>Z et Z=0)
C NTDL   : NOMBRE TOTAL DE DEGRES DE LIBERTE DU MAILLAGE DE L'OBJET
C NDSM   : NOMBRE DE VECTEURS SOLUTION
C TEMPER : TEMPER(NTDL,NDSM) LES NDSM VECTEURS SOLUTION
C NBCOOR : NOMBRE DE COORDONNEES DES NOEUDS
C MNX    : ADRESSE MCN DES NDIM COORDONNEES   DES POINTS DE L'EF COURANT
C MNNODL : ADRESSE MCN DU TABLEAU DES NUMEROS DES NOEUDS DE L'EF COURANT
C NBTYEL : NOMBRE DE TYPES D'EF DE L'OBJET
C MNNPEF : ADRESSE MCN DU TMC DES ADRESSES MCN DES TABLEAUX NPEF"TYPE EF
C NDPGST : CODE DE STOCKAGE DES SOMMETS POINTS NOEUDS
C MNTPOB : ADRESSE MCN DES TABLEAUX POLYNOMES DE BASE DES TYPES D'EF
C MXPOBA : NOMBRE MAXIMAL DE TABLEAUX POBA PAR TYPE D'EF
C MNTAUX : ADRESSE MCN DES TABLEAUX AUXILAIRES
C MNXYZP : ADRESSE MCN DE TMS XYZSOMMET DE L'OBJET
C
C SORTIE :
C --------
C IERR   : 0 SI PAS D'ERREUR, 1 SI EF NON PROGRAMME, 2 SI EF DEGENERE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET  LJLL UPMC & SAINT PIERRE DU PERRAY   MAI 2009
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
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
      DOUBLE PRECISION  TEMPER(NTDL,NDSM), EXPO, VOLEF, VOLUME, S, P
      DOUBLE PRECISION  DELTA
      CHARACTER*4       NOMELE(2)
C
C     LE NOMBRE DE MOTS D'UNE VARIABLE REEL2
      MOREE2 = MOTVAR(6)
      IERR   = 0
      VOLUME = 0D0
      MNPOL  = 0
      MNPOID = 0
      MNPDEL = 0
      MNF2   = 0
      MNDPOL = 0
C
C     LES INTEGRALES DES NDSM TEMPERATURES
C     ====================================
      CALL TNMCDC( 'REEL2', NDSM, MNINTP )
      CALL AZEROD( NDSM, MCN(MNINTP) )
C
C     LA BOUCLE SUR LES TYPES D'ELEMENTS FINIS
C     ========================================
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
     %         10,10,10,10, 1, 1,1, 1,10,10,
     %         10,10,10, 1), NUTYEL
C
C        ERREUR
 1       NBLGRC(NRERR) = 1
         IF( LANGAG .EQ. 0 ) THEN
            KERR(1) = 'ERREUR THNORM: TYPE EF '// NOMELE(1)
     %           // NOMELE(2) //' NON PROGRAMME'
         ELSE
            KERR(1) = 'ERROR THNORM: FE TYPE '// NOMELE(1)
     %           // NOMELE(2) //' NOT PROGRAMMED'
         ENDIF
         CALL LEREUR
         IERR = 1
         GOTO 9999
C
C        RECHERCHE DU TABLEAU DE POBA ET PARTAGE EN P ET DP
 10      L = MNTPOB + (NOTYEL-1) * MXPOBA
C
         IF( NUTYEL .EQ. 30 ) THEN
C           6CUBE 6Q1C N'A POUR L'INSTANT PAS DE FACE...
            NBNSOM = 0
            NARET  = 0
            NFACE  = 0
            GOTO 12
         ENDIF
C
C        LES VALEURS DES POLYNOMES DE L'EF DE DIMENSION TOTALE
C        EN 2D SURFACE DE REFERENCE, EN 3D ou 6D VOLUME DE REFERENCE
C        ...........................................................
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
         MNF1   = MNTAUX
         MNDP   = MNF1 + MOREE2 * NDIM
         MNDFM1 = MNDP + MOREE2 * NDIM * NBPOLY
C        AU TOTAL = MNDFM1 + MOREE2 * NDIM * NDIM
C
C        LES TEMPERATURES SUR L'EF
C        =========================
 30      CALL TNMCDC( 'REEL2', NBPOLY*NDSM, MNTEEF )
C
C        LA BOUCLE SUR LES ELEMENTS FINIS DE CE TYPE NUTYEL
C        ==================================================
         DO 100 NUELEM = 1, NBELEM
C
C           LES NOEUDS DE L'ELEMENT FINI NUELEM
C           -----------------------------------
            CALL EFNOEU( MNELE, NUELEM, NBNDEL, MCN(MNNODL) )
C
C           LES TEMPERATURES AUX NOEUDS DE L'EF
C           -----------------------------------
            MN = (MNTEEF+1)/2
            DO 36 N=1,NDSM
               DO 35 I=1,NBPOLY
                  DMCN(MN) = TEMPER( MCN(MNNODL-1+I), N )
                  MN = MN + 1
 35            CONTINUE
 36         CONTINUE
C
C           LES COORDONNEES DES NBPOE POINTS DE L'ELEMENT FINI NUELEM
C           ---------------------------------------------------------
            CALL EFXYZP( NDIM, MNXYZP, NBELEM, NUELEM, MNPGEL, NBPOE,
     %                   RMCN(MNX) )
C
C           ===================================================================
C           LE CALCUL DES TABLEAUX AUXILIAIRES ET DE LA MATRICE DE CONDUCTIVITE
C           ===================================================================
            GOTO( 41, 41, 41, 41,  1,  1, 1,  1,  1,  1,
     %             1,  1, 42,  1, 41, 41, 1, 41, 49, 51,
     %            51, 51, 51, 51,  1,  1, 1,  1, 41, 61,
     %            51, 51, 40,  1  ), NUTYEL
C
C           ********************************************
C           1D LAGRANGE ISOPARAMETRIQUE
C           ********************************************
C           LE CALCUL DES TABLEAUX AUXILIAIRES EN 2D OU AXISYMETRIE
 40         CALL E11LAG( NBPOLY, NPI, MCN(MNPOID),
     %                   MCN(MNPOL),  MCN(MNDPOL),
     %                   RMCN(MNX) ,  MCN(MNF1),
     %                   MCN(MNPDEL), MCN(MNDP), DELTA )
C
C           LA NORME L1 DES NDSM TEMPERATURES**EXPO
            CALL TL1LAG( EXPO, NBPOLY, NPI,  MCN(MNPOL), MCN(MNPDEL),
     &                   NDSM, MCN(MNTEEF),
     &                   VOLEF, MCN(MNINTP) )
            GOTO 80
C
C           ********************************************
C           2D OU AXISYMETRIQUE LAGRANGE ISOPARAMETRIQUE
C           ********************************************
C           LE CALCUL DES TABLEAUX AUXILIAIRES EN 2D OU AXISYMETRIE
 41         CALL E12LAG( NBPOLY, NPI, MCN(MNPOID),
     %                   MCN(MNPOL),  MCN(MNDPOL),
     %                   RMCN(MNX) ,  MCN(MNF1), MCN(MNF2),
     %                   MCN(MNPDEL), MCN(MNDP), MCN(MNDFM1) )
C
C           LA NORME L1 DES NDSM TEMPERATURES**EXPO
            CALL TL1LAG( EXPO, NBPOLY, NPI,  MCN(MNPOL), MCN(MNPDEL),
     &                   NDSM, MCN(MNTEEF),
     &                   VOLEF, MCN(MNINTP) )
            GOTO 80
C
C           ***************************
C           3D LAGRANGE ISOPARAMETRIQUE
C           ***************************
C           LE CALCUL DES TABLEAUX AUXILIAIRES EN 3D
 51         CALL E13LAG( NBPOLY, NPI, MCN(MNPOID),
     %                   MCN(MNPOL),  MCN(MNDPOL),
     %                   RMCN(MNX) ,  MCN(MNF1),
     %                   MCN(MNPDEL), MCN(MNDP), MCN(MNDFM1) )
C
C           LA NORME L1 DES NDSM TEMPERATURES**EXPO
            CALL TL1LAG( EXPO, NBPOLY, NPI,  MCN(MNPOL), MCN(MNPDEL),
     &                   NDSM, MCN(MNTEEF),
     &                   VOLEF, MCN(MNINTP) )
            GOTO 80
C
C           ************************************
C           2D TRIANGLE TRIA 2P1D LAGRANGE DROIT
C           ************************************
C           LA NORME L1 DES NDSM TEMPERATURES**EXPO
 42         CALL TL12P1D( EXPO, RMCN(MNX), NDSM, MCN(MNTEEF),
     &                   VOLEF, MCN(MNINTP) )
            GOTO 80
C
C           *************************************
C           3D TETRAEDRE TETR 3P1D LAGRANGE DROIT
C           *************************************
 49         CALL TL13P1D( EXPO, RMCN(MNX), NDSM, MCN(MNTEEF),
     &                    VOLEF, MCN(MNINTP) )
            GOTO 80
C
C           ***************************************
C           6D LAGRANGE ISOPARAMETRIQUE  6CUBE 6Q1C
C           ***************************************
 61         CALL E16Q1C( NBPOLY, NPI, MCN(MNPOID),
     %                   MCN(MNPOL),  MCN(MNDPOL),
     %                   RMCN(MNX) ,  MCN(MNF1),
     %                   MCN(MNPDEL), MCN(MNDP), MCN(MNDFM1), IERR )
C           SI EF DEGENERE => RETOUR
            IF( IERR .NE. 0 ) RETURN
C
C           LA NORME L1 DES NDSM TEMPERATURES**EXPO
            CALL TL1LAG( EXPO, NBPOLY, NPI,  MCN(MNPOL), MCN(MNPDEL),
     &                   NDSM, MCN(MNTEEF),
     &                   VOLEF, MCN(MNINTP) )
C
 80         VOLUME = VOLUME + VOLEF
 100     CONTINUE
C
C        DESTRUCTION DES TEMPERATURES DE L'EF
         CALL TNMCDS( 'REEL2', NBPOLY*NDSM, MNTEEF )
 500  CONTINUE
C     FIN DE LA BOUCLE SUR LES TYPES D'EF
C
C     AFFICHAGE DE LA NORME L1 DES SOLUTIONS**EXPO
C     ============================================
      WRITE(IMPRIM,*)
      IF( NDIM .EQ. 2 ) THEN
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10502) VOLUME
         ELSE
            WRITE(IMPRIM,20502) VOLUME
         ENDIF
      ELSE
         IF( NBCOOR .EQ. 3 ) THEN
            IF( LANGAG .EQ. 0 ) THEN
               WRITE(IMPRIM,10503) VOLUME
            ELSE
               WRITE(IMPRIM,20503) VOLUME
            ENDIF
         ELSE IF( NBCOOR .EQ. 6 ) THEN
            IF( LANGAG .EQ. 0 ) THEN
               WRITE(IMPRIM,10506) VOLUME
            ELSE
               WRITE(IMPRIM,20506) VOLUME
            ENDIF
         ENDIF
      ENDIF
10502 FORMAT('SURFACE DE L''OBJET=',D13.6)
10503 FORMAT('VOLUME DE L''OBJET=',D13.6)
10506 FORMAT('HYPER VOLUME DE L''OBJET=',D13.6)
20502 FORMAT('SURFACE of the OBJECT=',D13.6)
20503 FORMAT('VOLUME of the OBJECT=',D13.6)
20506 FORMAT('HYPER VOLUME of the OBJECT=',D13.6)
C
      IF( LANGAG .EQ. 0 ) THEN
         WRITE(IMPRIM,10511) NDSM,EXPO
      ELSE
         WRITE(IMPRIM,20511) NDSM,EXPO
      ENDIF
10511 FORMAT('NORME L1 DES',I3,' SOLUTION(S) **',D13.6)
20511 FORMAT('L1 NORM of',I3,' SOLUTION(S) **',D13.6)
      MN = (MNINTP+1)/2
      p  = expo / (expo-2)
      DO 550 N=1,NDSM
         S = DMCN(MN)
         IF( LANGAG .EQ. 0 ) THEN
            WRITE(IMPRIM,10550) N,S,S/VOLUME,1/p,(S/VOLUME)**(1/p)
         ELSE
            WRITE(IMPRIM,20550) N,S,S/VOLUME,1/p,(S/VOLUME)**(1/p)
         ENDIF
         MN = MN + 1
 550  CONTINUE
10550 FORMAT('SOLUTION',I3,'  NORME=',D15.8,'  NORME/VOLUME=',D15.8,
     %'  (NORME/VOLUME)**',G13.6,'=',D15.8 )
20550 FORMAT('SOLUTION',I3,'  NORM=',D15.8,'  NORM/VOLUME=',D15.8,
     %'  (NORM/VOLUME)**',G13.6,'=',D15.8 )
      WRITE(IMPRIM,*)
C
C     DESTRUCTION DU TABLEAU DEVENU INUTILE
      CALL TNMCDS( 'REEL2', NDSM, MNINTP )
C
 9999 RETURN
      END
