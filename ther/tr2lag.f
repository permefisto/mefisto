      SUBROUTINE TR2LAG( D2PI,   NOAXIS, X,      PENALI, NBJEUX, JEU,
     %                   NBSOMT, NOOBPS, NUMIPO, NUMAPO, LTDEPO,
     %                   NBPOLA, NPIA,   POIDSA, POLYA,  DPOLYA,
     %                   NBCOTE, NOOBLA, NUMILI, NUMALI, LTDELI,
     %                   NBPOLY, NPI,    POLY,
     %                   NOOBSF, NUMISU, NUMASU, LTDESU,
     %                   F1,     F2,     POIDEL, DP,
     %                   COND,   COEFTE, CONDUC, IERR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DE LA MATRICE DE CONDUCTIVITE DES ELEMENTS FINIS 2D
C -----    LAGRANGE ISOPARAMETRIQUES SAUF TRIANGLE 2P1D DE DEGRE 1
C
C ENTREES:
C --------
C D2PI   : 2 FOIS PI
C NOAXIS : 1 SI PROBLEME AXISYMETRIQUE, 0 SINON
C X      : LES 2 COORDONNEES DES NBPOLY POINTS DE L'ELEMENT FINI
C PENALI : COEFFICIENT DE PENALISATION DE LA CONDITION DE DIRICHLET
C NBJEUX : NOMBRE DE JEUX DE DONNEES
C JEU    : NUMERO DU JEU  DE DONNEES POUR CE CALCUL DE LA MATRICE ELEMENTAIRE
C
C NBSOMT : NOMBRE DE SOMMETS DE L'EF
C NOOBPS : NUMERO DE POINT DES SOMMETS
C NUMIPO : NUMERO MINIMAL DES OBJETS POINTS
C NUMAPO : NUMERO MAXIMAL DES OBJETS POINTS
C LTDEPO : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES THERMIQUES DES POINTS
C
C NBPOLA : NOMBRE DE POLYNOMES DE BASE SUR UN COTE DE L ELEMENT
C NPIA   : NOMBRE DE POINTS D INTEGRATION SUR UN COTE
C POIDSA : POIDS DES POINTS D INTEGRATION SUR UN COTE
C POLYA  : VALEUR DES POLYNOMES DE BASE AUX POINTS D INTEGRATION DU COTE
C DPOLYA : VALEUR DU GRADIENT DE CES MEMES POLYNOMES AUX POINTS ...
C
C NBCOTE : NOMBRE DES COTES DE L ELEMENT FINI
C NOOBLA : NUMERO DES OBJETS LIGNES DES ARETES DE L'ELEMENT
C NUMILI : NUMERO MINIMAL DES OBJETS LIGNES
C NUMALI : NUMERO MAXIMAL DES OBJETS LIGNES
C LTDELI : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES THERMIQUES DES LIGNES
C
C NBPOLY : NOMBRE DE POLYNOMES DE L'ELEMENT SURFACE
C NPI    : NOMBRE DE POINTS D INTEGRATION NUMERIQUE SUR LA SURFACE
C POLY   : VALEUR DES POLYNOMES DE BASE AUX POINTS D'INTEGRATION
C          POLY(I,L)= P(I) (XL)
C
C NOOBSF : NUMERO DE L'OBJET SURFACE DE CET ELEMENT
C NUMISU : NUMERO MINIMAL DES OBJETS SURFACES
C NUMASU : NUMERO MAXIMAL DES OBJETS SURFACES
C LTDESU : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES CONDUCTIVITE DES SURFACES
C
C F1     : COORDONNEES XX DES NPI POINTS D INTEGRATION DE L ELEMENT
C F2     : COORDONNEES YY DES NPI POINTS D INTEGRATION DE L ELEMENT
C POIDEL : DELTA * POIDS(NPI) DES NPI POINTS D INTEGRATION
C DP     : GRADIENT DES POLYNOMES DE BASE SUR L ELEMENT COURANT
C
C SORTIES:
C --------
C COND   : TENSEUR SYMETRIQUE DE CONDUCTIVITE (6 COEFFICIENTS)
C COEFTE : COEFFICIENT DE LA TEMPERATURE AUX POINTS D'INTEGRATION
C CONDUC : MATRICE ELEMENTAIRE DE CONDUCTIVITE
C IERR   : 7 SI PB AXISYMETRIQUE AVEC UNE ABSCISSE X<=0 EN UN POINT
C            D'INTEGRATION NUMERIQUE
C          0 SINON
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS      Octobre 1990
C MODIFS: ALAIN PERRONNET LJLL UPMC Saint Pierre du Perray     Aout 2011
C MODIFS: ALAIN PERRONNET LJLL UPMC Saint Pierre du Perray Decembre 2013
C AJOUTS: ALAIN PERRONNET LJLL UPMC Saint Pierre du Perray     Mars 2014
C23456---------------------------------------------------------------012
      include"./incl/gsmenu.inc"
      include"./incl/donthe.inc"
      include"./incl/cthet.inc"
      include"./incl/cnonlin.inc"
      include"./incl/ctemps.inc"
C
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE      (MCN(1), RMCN(1), DMCN(1))
C
      DOUBLE PRECISION D2PI,POIDSA(NPIA),POLYA(NBPOLA,NPIA),
     %                 DPOLYA(NBPOLA,NPIA)
      DOUBLE PRECISION POLY(NBPOLY,NPI),
     %                 F1(NPI),F2(NPI),POIDEL(NPI),
     %                 DP(2,NBPOLY,NPI),CONDUC(*)
      DOUBLE PRECISION PENALI,ECHANG,COND(6),CONDP(2,2),COEFTE(NPI)
      DOUBLE PRECISION Rho
      REAL             X(NBPOLY,2)
      INTEGER          LTDEPO( 1:MXDOTH, 1:NBJEUX, NUMIPO:NUMAPO )
      INTEGER          LTDELI( 1:MXDOTH, 1:NBJEUX, NUMILI:NUMALI )
      INTEGER          LTDESU( 1:MXDOTH, 1:NBJEUX, NUMISU:NUMASU )
      INTEGER          NOOBPS( 1:NBSOMT )
      INTEGER          NOOBLA( 1:NBCOTE )
C
      DOUBLE PRECISION GL(2),DGL(2),DELTA,PROSCD,S, XYZPI(3)
      INTEGER          NOPOAR(3)

C     LE NOMBRE DE MOTS D'UNE VARIABLE REEL2
      MOREE2 = MOTVAR(6)
C
C     INITIALISATION A ZERO DE LA MATRICE DE CONDUCTIVITE ELEMENTAIRE
C     ---------------------------------------------------------------
      CALL AZEROD( NBPOLY*(NBPOLY+1)/2, CONDUC )
C
C     ==========================
C     CONTRIBUTION DE LA SURFACE
C     ==========================
      IF( TESTNL .GT. 0 ) CALL NLDATA0( NBPOLY )
C     RECUPERATION DE LA TEMPERATURE AUX NBPOLY DL DE L ELEMENT FINI ou
C     RECUPERATION DE L'ONDE COMPLEXE AU TEMPS ACTUEL et INITIAL
C
C     SI LA CONDUCTIVITE N'EST PAS DECLAREE, SAUT DU CALCUL DE LA CONDUCTIVITE
      MNLT = LTDESU(LPCOND,JEU,NOOBSF)
      IF( MNLT .EQ. 0 ) GOTO 40
C
C     CONTRIBUTION DE LA CONDUCTIVITE
C     -------------------------------
      DO 30 L=1,NPI
C
         IF( TESTNL .GT. 0 ) THEN
C           CALCUL DE LA TEMPERATURE AU POINT D'INTEGRATION L
            TEMPEL = PROSCD( POLY(1,L), MCN(MNTHDL), NBPOLY )
         ENDIF
C
C        RECHERCHE DU TENSEUR SYMETRIQUE DE CONDUCTIVITE AU POINT
C        D'INTEGRATION L
         XYZPI(1) = F1(L)
         XYZPI(2) = F2(L)
         XYZPI(3) = 0D0
         CALL RECOND( 3, NOOBSF, 3, XYZPI, MNLT, COND )
C        si Gross-Pitaevskii: COND = -1/2 [Identite]

         IF( TESTNL .EQ. 9 ) THEN
C           [KG] = [Rho/PasTemps -N(V0**2+W0**2) + Alfa LAPLACIEN]
            COND(1) = -COND(1)
            COND(2) = -COND(2)
            COND(3) = -COND(3)
            COND(4) = -COND(4)
         ENDIF
C
         IF( NOAXIS .NE. 0 ) THEN
C           EF AXISYMETRIQUE
            IF( F1(L) .LE. 0D0 ) THEN
               NBLGRC(NRERR) = 2
               KERR(1) = 'ERREUR: PB AXISYMETRIQUE AVEC EF'
               KERR(2) = 'D''ABSCISSE X=RAYON NEGATIVE OU NULLE'
               CALL LEREUR
               IERR = 7
               RETURN
            ENDIF
            DELTA = POIDEL(L) * D2PI * F1(L)
         ELSE
C           EF NON AXISYMETRIQUE
            DELTA = POIDEL(L)
         ENDIF
C
C        COND = CONDUCTIVITE * POIDS * DELTA
         CONDP(1,1) = COND(1) * DELTA
         CONDP(2,1) = COND(2) * DELTA
         CONDP(1,2) = COND(2) * DELTA
         CONDP(2,2) = COND(3) * DELTA
C
C        CONDUC = CONDUC + T(DP) * CONDP * (DP)
         M = 0
         DO 25 I=1,NBPOLY
            DO 20 J=1,I
               S = 0D0
               DO 10 K=1,2
                  S = S + DP(K,I,L) * ( CONDP(K,1) * DP(1,J,L)
     %                                + CONDP(K,2) * DP(2,J,L) )
 10            CONTINUE
               M = M + 1
               CONDUC(M) = CONDUC(M) + S
 20         CONTINUE
 25      CONTINUE
 30   CONTINUE
C
C     CONTRIBUTION DU COEFFICIENT DEVANT LA TEMPERATURE
C     -------------------------------------------------
 40   MNLT = LTDESU(LPCOET,JEU,NOOBSF)
      IF( MNLT .GT. 0 ) THEN
         DO 50 L=1,NPI
C
C           RECHERCHE DU COEFFICIENT DE LA TEMPERATURE AU POINT D'INTEGRATION L
            XYZPI(1) = F1(L)
            XYZPI(2) = F2(L)
            XYZPI(3) = 0D0

            IF( TESTNL .GT. 0 ) THEN
               CALL NLDATA1( NBPOLY, POLY(1,L) )
C              PB NON LINEAIRE: CALCUL DES VALEURS
C              TEMPEL: TEMPERATURE   ACTUELLE AU POINT D'INTEGRATION L ou
C              TEMPEL: PARTIE REELLE ACTUELLE AU POINT D'INTEGRATION L de L'ONDE
C              IF( TESTNL .GE. 6 ) THEN  ONDE NLSE COMPLEXE
C                ONDEPI: PARTIE IMAGINAIRE ACTUELLE AU POINT D'INTEGRATION L de L'ONDE
C                TEMPEL0:PARTIE REELLE     INITIALE AU POINT D'INTEGRATION L de L'ONDE
C                ONDEPI0:PARTIE IMAGINAIRE INITIALE AU POINT D'INTEGRATION L de L'ONDE
            ENDIF

            IF( TESTNL .LT. 5 ) THEN

C              COEFFICIENT DEVANT LA TEMPERATURE A L'INSTANT TEMPS
               CALL RECOET( 3, NOOBSF, 3, XYZPI, MNLT, COEFTE(L) )

            ELSE

C              NLSE COEFFICIENT N(V**2+W**2) DEVANT L'ONDE U = V +iW A L'INSTANT TEMPS
C              EXEMPLE: Gross-Pitaevskii => -V - Beta |V**2+W**2|
               CALL RENLSE( 3, NOOBSF, 3, XYZPI, TEMPS, TEMPEL, ONDEPI,
     %                      MNLT, COEFTE(L) )

               IF( TESTNL .EQ. 9 ) THEN
C                 [KG] = [Rho/PasTemps -N(V0**2+W0**2) + Alfa LAPLACIEN]
                  COEFTE(L) = -COEFTE(L)
               ENDIF

               IF( TESTNL .EQ. 6 .OR. TESTNL .EQ. 9 ) THEN
C                 GROSS-PITAEVSKII DEMANDE LA DENSITE DE MASSE Rho / DeltaT
                  MN = LTDESU(LPMAST,JEU,NOOBSF)
                  IF( MN .GT. 0 ) THEN
                     CALL REMASS( 3, NOOBSF, 3, XYZPI, MN, Rho )
                     COEFTE(L) = Rho/PasTemps + COEFTE(L)
                  ENDIF
               ENDIF

            ENDIF
C
            IF( NOAXIS .NE. 0 ) THEN
C              EF AXISYMETRIQUE
               DELTA = POIDEL(L) * D2PI * F1(L)
            ELSE
C              EF NON AXISYMETRIQUE
               DELTA = POIDEL(L)
            ENDIF
C
C           COEF TEMPERATURE = COEFTE * POIDS * DELTA
            COEFTE(L) = COEFTE(L) * DELTA
C
 50      CONTINUE
C
         M = 0
         DO 80 J=1,NBPOLY
            DO 70 I=1,J
               S = 0.D0
               DO 60 L=1,NPI
                  S = S +  COEFTE(L) * POLY(J,L) * POLY(I,L)
 60            CONTINUE
               M = M + 1
               CONDUC(M) = CONDUC(M) + S
 70         CONTINUE
 80      CONTINUE
      ENDIF
C
C     =======================
C     CONTRIBUTION DES ARETES
C     =======================
      DO 200 K=1,NBCOTE
C
C        LE NUMERO DE LIGNE DU COTE K
         NOOB = NOOBLA(K)
         IF( NOOB .GT. 0 ) THEN
C
C           LE COTE K EST SUR UNE LIGNE. SUPPORT D'ECHANGE OU CONTACT PENALISE?
            IECHAN = 0
            IF( LTDELI(LPECHA,JEU,NOOB) .GT. 0 ) IECHAN = 1
            IF( LTDELI(LPCONT,JEU,NOOB) .GT. 0 .AND.
     %          PENALI .NE. 0D0                ) IECHAN = 2
C
            IF( IECHAN .EQ. 0 ) GOTO 200
C
C           UN TABLEAU ECHANGE OU CONTACT PENALISE EXISTE POUR CETTE LIGNE
C           LE NUMERO DES POINTS DU COTE
            NOPOAR(1) = K
            IF( K .NE. NBCOTE ) THEN
               NOPOAR(2) = K+1
            ELSE
               NOPOAR(2) = 1
            ENDIF
C           LE NUMERO DU POINT MILIEU
            NOPOAR(3) = K + NBCOTE
C
            IF( TESTNL .GT. 0 ) THEN
C              PB NON LINEAIRE:
C              RECUPERATION DE LA TEMPERATURE AUX NBPOLA DL DU COTE K
               MN = (MNTHET-1)/MOREE2
               DO 120 I=1,NBPOLA
                  NP = NOPOAR(I)
                  DMCN((MNTHDL-1)/2+I) = DMCN( MN+MCN(MNNODL+NP-1) )
 120           CONTINUE
            ENDIF
C
            DO 190 L=1,NPIA
C
C              CALCUL DES COORDONNEES DU POINT D INTEGRATION ET DU JACOBIEN ET D
               CALL E22LAG( NBPOLY,NBPOLA,NOPOAR,
     %                      POLYA(1,L),DPOLYA(1,L),
     %                      X,GL,DGL,DELTA)
C              EN SORTIE GL=LES 2 COORDONNEES DU POINT D'INTEGRATION
C
               IF( TESTNL .GT. 0 ) THEN
C                 CALCUL DE LA TEMPERATURE AU POINT D'INTEGRATION L
                  TEMPEL=PROSCD( POLYA(1,L), MCN(MNTHDL), NBPOLA )
               ENDIF
C
               IF( IECHAN .EQ. 1 ) THEN
C
C                 CALCUL DU COEFFICIENT D'ECHANGE
                  XYZPI(1) = GL(1)
                  XYZPI(2) = GL(2)
                  XYZPI(3) = 0D0
                  CALL REECHA( 2,NOOB, 3, XYZPI,
     %                         LTDELI(LPECHA,JEU,NOOB), ECHANG )
               ELSE
C
C                 CALCUL DU COEFFICIENT D'ECHANGE = CONTACT PENALISE
                  ECHANG = PENALI
C
               ENDIF
C
C              CALCUL DE TRANSPOSEE(P(GL)) * ECHANG * (P(GL))
               DELTA = DELTA * POIDSA(L) * ECHANG
C
C              SI ELEMENT AXISYMETRIQUE DELTA * 2 * PI * R
               IF( NOAXIS .NE. 0 ) DELTA = DELTA * D2PI * GL(1)
C
               DO 180 I=1,NBPOLA
                  DO 170 J=1,I
C                    LE STOCKAGE SYMETRIQUE IMPOSE NO LIGNE>= NO COLONNE
                     IF( NOPOAR(I) .GE. NOPOAR(J) ) THEN
                         NI = NOPOAR(I)
                         NJ = NOPOAR(J)
                     ELSE
                         NI = NOPOAR(J)
                         NJ = NOPOAR(I)
                     ENDIF
                     NJ = NI * (NI - 1) / 2 + NJ
                     CONDUC(NJ) = CONDUC(NJ) + DELTA * POLYA(I,L)
     %                                               * POLYA(J,L)
 170              CONTINUE
 180           CONTINUE
 190        CONTINUE
         ENDIF
 200  CONTINUE
C
C     =====================================================
C     CONTRIBUTION DES SOMMETS A LA PENALISATION DU CONTACT
C     =====================================================
      IF( PENALI .NE. 0D0 ) THEN
         DO 300 K=1,NBSOMT
C
            NOOB = NOOBPS(K)
            IF( NOOB .GT. 0 ) THEN
C
C              LE SOMMET K EST UN POINT. EST IL SUPPORT D'UN CONTACT PENALISE?
               IF( LTDEPO(LPCONT,JEU,NOOB) .GT. 0 ) THEN
C
C                 OUI: UN TABLEAU CONTACT PENALISE EXISTE POUR CE POINT
C                 CALCUL DU COEFFICIENT D'ECHANGE = CONTACT PENALISE
                  DELTA = PENALI
C
C                 SI ELEMENT AXISYMETRIQUE DELTA * 2 * PI * R
                  IF( NOAXIS .NE. 0 ) DELTA = DELTA * D2PI * X(K,1)
C
C                 LES SOMMETS SONT NUMEROTES EN PREMIER
C                 PUIS, VIENNENT LES EVENTUELS MILIEUX DES COTES
C                 COEFFICIENT DIAGONAL NI DE LA MATRICE DE CONDUCTIVITE ELEMENTA
                  NJ = K * ( K + 1 ) / 2
                  CONDUC(NJ) = CONDUC(NJ) + DELTA
C
               ENDIF
            ENDIF
 300     CONTINUE
      ENDIF
C
      RETURN
      END
