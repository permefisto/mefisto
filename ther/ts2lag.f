      SUBROUTINE TS2LAG( D2PI,   NOAXIS, X,      PENALI, NBJEUX, JEU,
     %                   NBSOMT, NOOBPS, NUMIPO, NUMAPO, LTDEPO,
     %                   NBPOLA, NPIA,   POIDSA, POLYA,  DPOLYA,
     %                   NBCOTE, NOOBLA, NUMILI, NUMALI, LTDELI,
     %                   NBPOLY, NPI,    POLY,
     %                   NOOBSF, NUMISU, NUMASU, LTDESU,
     %                   NBNOVI, VITEGt,
     %                   F1,     F2,     POIDEL, DP,
     %                   NOPART, BE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DU SECOND MEMBRE ELEMENTAIRE DES EF AXISYMETRIQUES
C -----    ou 2D LAGRANGE ISOPARAMETRIQUES

C ENTREES:
C --------
C D2PI   : 2 FOIS PI
C NOAXIS : 1 SI PROBLEME AXISYMETRIQUE,  0 SINON
C X      : COORDONNEES RAYON ET COTE DES NBPOLY POINTS DE L'EF
C          OU X Y DES NBPOLY POINTS DE L'EF
C PENALI : COEFFICIENT DE PENALISATION DE LA CONDITION DE DIRICHLET
C NBJEUX : NOMBRE DE JEUX DE DONNEES
C JEU    : NUMERO DU JEU  DE DONNEES POUR CE CALCUL DE LA MATRICE ELEMENTAIRE

C NBSOMT : NOMBRE DE SOMMETS DE L'EF
C NOOBPS : NUMERO DE POINT DES SOMMETS
C NUMIPO : NUMERO MINIMAL DES OBJETS POINTS
C NUMAPO : NUMERO MAXIMAL DES OBJETS POINTS
C LTDEPO : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES THERMIQUES DES POINTS

C NBPOLA : NOMBRE DE POLYNOMES DE BASE SUR UN COTE DE L'EF
C NPIA   : NOMBRE DE POINTS D INTEGRATION SUR UN COTE
C POIDSA : POIDS DES POINTS D INTEGRATION SUR UN COTE
C POLYA  : VALEUR DES POLYNOMES DE BASE AUX POINTS D INTEGRATION COTE
C DPOLYA : VALEUR DU GRADIENT DE CES MEMES POLYNOMES AUX POINTS ...

C NBCOTE : NOMBRE DES COTES DE L'EF
C NOOBLA : NOUMERO DES OBJETS LIGNES DES ARETES DE L'ELEMENT
C SOURCE : FLUX NORMAL DE TEMPERATURE A LA PAROI
C NUMILI : NUMERO MINIMAL DES OBJETS LIGNES
C NUMALI : NUMERO MAXIMAL DES OBJETS LIGNES
C LTDELI : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES DE THERMIQUE
C          DES OBJETS LIGNES

C NBPOLY : NOMBRE DE POLYNOMES
C NPI    : NOMBRE DE POINTS D INTEGRATION NUMERIQUE SUR LE RECTANGLE
C POLY   : POLY(I, L) = PI(RL, ZL)

C NOOBSF : NUMERO DE L'OBJET SURFACE DE CET ELEMENT
C NUMISU : NUMERO MINIMAL DES OBJETS SURFACES
C NUMASU : NUMERO MAXIMAL DES OBJETS SURFACES
C LTDESU : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES DE THERMIQUE
C          DES OBJETS SURFACES

C F1     : RAYON R DES NPI POINTS D'INTEGRATION
C F2     : COTE  Z DES NPI POINTS D'INTEGRATION
C POIDEL : DELTA * POIDS(NPI) DES NPI POINTS D INTEGRATION
C DP     : DP(2, NBPOLY, NPI) GRADIENT AUX POINTS D 'INTEGRATION DES
C          FONCTIONS DE BASE ISOPARAMETRIQUES
C X      : COORDONNEES RAYON ET COTE DES NBPOLY POINTS DE L'EF
C NOPART : POUR NLSE SEULEMENT AU NIVEAU DE SOURCE=FORCE et CONTACT=FIXATION
C          1 SI PARTIE REELLE TRAITEE ou 2 SI PARTIE IMAGINAIRE TRAITEE
C          0 SI INACTIF (CAS THERMIQUE STANDARD D'UNE SOURCE)

C NBNOVI : NOMBRE DE NOEUDS VITESSE DU MAILLAGE
C VITEGt : VECTEUR(NBNOVI,2) DES 2 COMPOSANTES DE LA VITESSE EN TOUS
C          LES NOEUDS DU MAILLAGE au TEMPS TEMPS
C          CE TABLEAU EST UTILISE SEULEMENT EN CAS d'UN TRANSPORT
C          DE TEMPERATURE A CETTE VITESSE POUR LE SECOND MEMBRE
C          terme: -(VITEGt . Grad) Temperature

C SORTIE :
C --------
C BE     : BE(NBPOLY) LE SECOND MEMBRE ELEMENTAIRE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     OCTOBRE 1990
C MODIFS : ALAIN PERRONNET TEXAS A & M University at QATAR  FEVRIER 2011
C MODIFS : ALAIN PERRONNET Saint Pierre du Perray             Avril 2022
C23456---------------------------------------------------------------012
      include"./incl/donthe.inc"
      include"./incl/cthet.inc"
      include"./incl/cnonlin.inc"
      include"./incl/a___fixation.inc"
      include"./incl/a___vitessefluide.inc"

      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE      (MCN(1), RMCN(1), DMCN(1))

      DOUBLE PRECISION  D2PI, SOURCE(2), VN(2), VITEGt(NBNOVI,2)
      DOUBLE PRECISION  POIDSA(NPIA), POLYA(NBPOLA, NPIA),
     %                  DPOLYA(NBPOLA, NPIA),
     %                  POLY(NBPOLY, NPI),
     %                  POIDEL(NPI),
     %                  DP(2, NBPOLY, NPI), F1(NPI), F2(NPI), XYZ(3),
     %                  BE(NBPOLY),
     %                  VITEBL(2), VITEF(9,2),
     %                  PENALI
      REAL              X(NBPOLY, 2)
      INTEGER           NOOBPS( 1:NBSOMT )
      INTEGER           NOOBLA( 1:NBCOTE )
      INTEGER           LTDEPO( 1:MXDOTH, 1:NBJEUX, NUMIPO:NUMAPO )
      INTEGER           LTDELI( 1:MXDOTH, 1:NBJEUX, NUMILI:NUMALI )
      INTEGER           LTDESU( 1:MXDOTH, 1:NBJEUX, NUMISU:NUMASU )

      INTEGER           NOPOAR(3)
      DOUBLE PRECISION  D, GL(2), DGL(2), DELTA, PROSCD

C     MISE A ZERO DE BE LE VECTEUR ELEMENTAIRE
C     ----------------------------------------
      CALL AZEROD( NBPOLY, BE )

      IF( TESTNL .GE. 1 .AND. MNTHET .GT. 0 ) THEN
C        LES SOURCES DEPENDENT DE LA TEMPERATURE
C        RECUPERATION DE LA TEMPERATURE AUX NBPOLY DL DE L ELEMENT FINI
         MN  = (MNTHET-1)/2
         MNT = (MNTHDL-1)/2
         DO I=1,NBPOLY
            NONOEU = MCN( MNNODL + (I-1) )
            DMCN( MNT+I ) = DMCN( MN+NONOEU )
         ENDDO
      ENDIF
      TEMPEL = 0D0
      ONDEPI = 0D0

C     ====================================
C     CONTRIBUTION DES SOURCES SURFACIQUES
C     ====================================
      IF( LTDESU(LPSOUR,JEU,NOOBSF) .GT. 0 ) THEN

         MNTHDLD = (MNTHDL+1 ) / 2
         DO L=1,NPI

C           LA CONTRIBUTION DES SOURCES DE CHALEUR
            IF( TESTNL .GE. 1 .AND. MNTHET .GT. 0 ) THEN
C              CALCUL DE LA TEMPERATURE AU POINT D'INTEGRATION L
               TEMPEL = PROSCD( POLY(1,L), DMCN(MNTHDLD), NBPOLY )
            ENDIF

C           LA VALEUR DES SOURCES SURFACIQUES EN CE POINT D'INTEGRATION
            XYZ(1) = F1(L)
            XYZ(2) = F2(L)
            XYZ(3) = 0D0
            IF( TESTNL .LE. 5 ) THEN
               CALL RESOUR( 3, NOOBSF, 3, XYZ,
     %                      LTDESU(LPSOUR,JEU,NOOBSF), SOURCE )
            ELSE
               CALL REFORC( 3, NOOBSF, 2, XYZ(1),XYZ(2),XYZ(3),
     %                                    0D0,   0D0,   0D0,
     %                      LTDESU(LPSOUR,JEU,NOOBSF), SOURCE )
               SOURCE(1) = SOURCE(NOPART)
            ENDIF

            IF( NOAXIS .EQ. 0 ) THEN
C              EF NON AXISYMETRIQUE
               D = POIDEL(L)
            ELSE
C              EF AXISYMETRIQUE => * 2 PI R(PT INTEGRATION L)
               D = POIDEL(L) * D2PI * F1(L)
            ENDIF

            DO I=1,NBPOLY
               BE(I) = BE(I) + D * POLY(I,L) * SOURCE(1)
            ENDDO

         ENDDO
      ENDIF

C     LA CONTRIBUTION DU TRANSPORT: - VITESSE * GRADIENT TEMPERATURE
C     ACTUELLEMENT CE TERME EST SOUSTRAIT DU SECOND MEMBRE + POINT FIXE
C     -----------------------------------------------------------------
      MNVIFL = LTDESU(LPVIFL,JEU,NOOBSF)
      IF( MNVIFL .GT. 0 ) THEN

C        TYPE DE LA DONNEE DE LA VITESSE DU FLUIDE
         LTVIFL = MCN( MNVIFL + WTVIFL )

C        RECUPERATION DE LA TEMPERATURE AUX NBPOLY DL DE L ELEMENT FINI
         MN  = (MNTHET-1)/2
         MNT = (MNTHDL-1)/2
         DO I=1,NBPOLY
            NONOEU = MCN( MNNODL + (I-1) )
            DMCN( MNT+I ) = DMCN( MN+NONOEU )
         ENDDO

         IF( LTVIFL .LT. 2 ) THEN

C           DONNEE: VITESSE du FLUIDE CONSTANTE ou FONCTION UTILISATEUR
C           -----------------------------------------------------------
            MNTHDLD = (MNTHDL+1)/2
            DO L=1,NPI

C              LA VALEUR DE LA VITESSE DU FLUIDE AU POINT D'INTEGRATION L
               IF( TESTNL .GE. 1 .AND. MNTHET .GT. 0 ) THEN
C                 CALCUL DE LA TEMPERATURE AU POINT D'INTEGRATION L
                  TEMPEL=PROSCD( POLY(1,L), DMCN(MNTHDLD), NBPOLY )
               ENDIF

C              COORDONNEES DU POINT D'INTEGRATION L
               XYZ(1) = F1(L)
               XYZ(2) = F2(L)
               XYZ(3) = 0D0
               CALL REVIFL( 3, NOOBSF, 2, 3, XYZ,
     %                      LTDESU(LPVIFL,JEU,NOOBSF), VITEBL )

C              - t[P] [V] [DP] [TEMPERATURE EF]
               DO I=1,NBPOLY

                  D = 0D0
                  DO K=1,NBPOLY
                     D = D + ( VITEBL(1)*DP(1,K,L) +
     %                         VITEBL(2)*DP(2,K,L) ) * DMCN(MNT+K)
                  ENDDO

                  IF( NOAXIS .EQ. 0 ) THEN
C                    EF NON AXISYMETRIQUE
                     D = D * POIDEL(L)
                  ELSE
C                    EF AXISYMETRIQUE => * 2 * PI * R(PT INTEGRATION L)
                     D = D * POIDEL(L) * D2PI * F1(L)
                  ENDIF
C
                  BE(I) = BE(I) - POLY(I,L) * D
C
               ENDDO
            ENDDO

         ELSE IF( LTVIFL .EQ. 2 ) THEN

C           DONNEE du VECTEUR VITESSE du FLUIDE
C           -----------------------------------
C           CALCUL DE LA VITESSE AU POINT D'INTEGRATION L
C           RECUPERATION DE LA VITESSE AUX NBPOLY DL DE L ELEMENT FINI
            DO I=1,NBPOLY
               NONOEU = MCN( MNNODL+(I-1) )
               DO K=1,2
                  VITEF(I,K) = VITEGt( NONOEU, K )
               ENDDO
            ENDDO

            DO L=1,NPI

C              VITESSE DU FLUIDE AU POINT D'INTEGRATION L
               VITEBL(1) = PROSCD( POLY(1,L), VITEF(1,1), NBPOLY )
               VITEBL(2) = PROSCD( POLY(1,L), VITEF(1,2), NBPOLY )

C              - t[P] [V] [DP] [TEMPERATURE EF]
               DO I=1,NBPOLY

                  D = 0D0
                  DO K=1,NBPOLY
                     D = D + ( VITEBL(1)*DP(1,K,L) +
     %                         VITEBL(2)*DP(2,K,L) ) * DMCN(MNT+K)
                  ENDDO

                  IF( NOAXIS .EQ. 0 ) THEN
C                    EF NON AXISYMETRIQUE
                     D = D * POIDEL(L)
                  ELSE
C                    EF AXISYMETRIQUE => * 2 * PI * R(PT INTEGRATION L)
                     D = D * POIDEL(L) * D2PI * F1(L)
                  ENDIF

                  BE(I) = BE(I) - POLY(I,L) * D

               ENDDO

            ENDDO

         ENDIF

      ENDIF
C
C     =======================================
C     CONTRIBUTIONS DES SOURCES SUR LES COTES
C     =======================================
      DO 50 K=1,NBCOTE
C
C        LE NUMERO DE LIGNE DU COTE K
         NOOB = NOOBLA(K)
         IF( NOOB .GT. 0 ) THEN
C
C           LE COTE K EST SUR UNE LIGNE. SUPPORT DE FLUX OU CONTACT PENALISE?
            IECHAN = 0
            IF( LTDELI(LPSOUR,JEU,NOOB) .GT. 0  ) IECHAN = 1
            IF( LTDELI(LPCONT,JEU,NOOB) .GT. 0 .AND.
     %          PENALI .NE. 0D0                 ) IECHAN = 2
C
            IF( IECHAN .EQ. 0 ) GOTO 50
C
C           UN TABLEAU SOURCE EXISTE POUR CETTE LIGNE
C           CALCUL DE LA CONTRIBUTION DE L'ARETE K A BE
C           ...........................................
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
            IF( TESTNL .GE. 1 ) THEN
C              RECUPERATION DE LA TEMPERATURE AUX DL DE L'ARETE K
               MN  = (MNTHET-1) / 2
               MNT = (MNTHDL-1) / 2
               DO I=1,NBPOLA
                  DMCN(MNT+I)=DMCN(MN+MCN(MNNODL+(NOPOAR(I)-1)))
               ENDDO
            ENDIF
C
            MNTHDLD = (MNTHDL+1)/2
            DO L=1,NPIA
C
C              CALCUL DES COORDONNEES DU POINT D INTEGRATION L
C              ET DU JACOBIEN EN CE POINT
               CALL E22LAG ( NBPOLY, NBPOLA, NOPOAR,
     %                       POLYA(1,L), DPOLYA(1,L),
     %                       X, GL, DGL, DELTA )
C              EN SORTIE GL=LES 2 COORDONNEES DU POINT D'INTEGRATION
C
C              CALCUL DU FLUX NORMAL
               IF( TESTNL .GE. 1 .AND. MNTHET .GT. 0 ) THEN
C                 CALCUL DE LA TEMPERATURE AU POINT D'INTEGRATION L
                  TEMPEL=PROSCD( POLYA(1,L), DMCN(MNTHDLD), NBPOLA )
               ENDIF
C
               IF( IECHAN .EQ. 1 ) THEN
C
C                 FLUX NORMAL DE CHALEUR: CONDITION NEUMANN OU FOURIER
                  XYZ(1) = GL(1)
                  XYZ(2) = GL(2)
                  XYZ(3) = 0D0
C
                  IF( TESTNL .LE. 5 ) THEN
C                    SOURCE(1) REQUISE
                     CALL RESOUR( 2, NOOB, 3, XYZ,
     %                            LTDELI(LPSOUR,JEU,NOOB), SOURCE )
                  ELSE
C                    FORCE(2) REQUISE
C                    LE VECTEUR NORMAL UNITAIRE
                     VN(1) =  DGL(2) / DELTA
                     VN(2) = -DGL(1) / DELTA
                     CALL REFORC( 2, NOOB, 2, XYZ(1),XYZ(2),XYZ(3),
     %                                        VN(1), VN(2), 0D0,
     %                            LTDELI(LPSOUR,JEU,NOOB), SOURCE )
                     SOURCE(1) = SOURCE(NOPART)
                  ENDIF
C
               ELSE
C
                  XYZ(1) = GL(1)
                  XYZ(2) = GL(2)
                  XYZ(3) = 0D0
                  IF( TESTNL .LE. 5 ) THEN
C                    CONTACT PENALISE = PENALI x TEMPERATURE
                     CALL RECONT( 2, NOOB, 3, XYZ,
     %                            LTDELI(LPCONT,JEU,NOOB), SOURCE )
                     SOURCE(1) = SOURCE(1) * PENALI
                  ELSE
C                    FIXATION(2) PENALISEE
                     MN = LTDELI(LPCONT,JEU,NOOB)
                     CALL REFIXA( 2, NOOB, XYZ(1),XYZ(2),XYZ(3), MN,
     %                            NBCOFI, SOURCE )
C                    NBCOFI LE NOMBRE DE COMPOSANTES FIXEES
                     DO I = 1, NBCOFI
C                       LE NUMERO DE LA COMPOSANTE FIXEE
                        NU = MCN( MN + WUCOFI - 1 + I )
                        IF( NU .EQ. NOPART ) THEN
                           SOURCE(1) = SOURCE(I) * PENALI
                        ENDIF
                     ENDDO
                  ENDIF
C
               ENDIF
C
               IF( NOAXIS .EQ. 0 ) THEN
C                 EF NON AXISYMETRIQUE
                  D = DELTA * POIDSA(L)
               ELSE
C                 EF AXISYMETRIQUE
                  D = DELTA * POIDSA(L) * D2PI * GL(1)
               ENDIF
C
               DO I=1,NBPOLA
                  II = NOPOAR(I)
                  BE(II) = BE(II) + D * POLYA(I,L) * SOURCE(1)
               ENDDO
            ENDDO
         ENDIF
 50   ENDDO
C
C     ============================================
C     CONTRIBUTION DU CONTACT PENALISE AUX SOMMETS
C     ============================================
      IF( PENALI .NE. 0D0 ) THEN
         DO K=1,NBSOMT
C
C           NO DE POINT DU SOMMET K
            NOOB = NOOBPS(K)
            IF( NOOB .GT. 0 ) THEN
C
C              LE SOMMET K EST UN POINT
C              EST IL SUPPORT D'UN CONTACT PENALISE?
               IF( LTDEPO(LPCONT,JEU,NOOB) .GT. 0 ) THEN
C
C                 OUI: UN TABLEAU CONTACT PENALISE EXISTE POUR CE POINT
C                 CONTACT PENALISE = PENALI x TEMPERATURE
                  IF( TESTNL .GE. 1 .AND. MNTHET .GT. 0 ) THEN
C                    PB NON LINEAIRE:
C                    RECUPERATION DE LA TEMPERATURE AU SOMMET K
                     TEMPEL = DMCN( (MNTHET-1)/2 + MCN(MNNODL+K-1) )
                  ENDIF
                  XYZ(1) = X(K,1)
                  XYZ(2) = X(K,2)
                  XYZ(3) = 0D0
                  IF( TESTNL .LE. 5 ) THEN
C                    CONTACT PENALISE = SOURCE(1)
                     CALL RECONT( 1, NOOB, 3, XYZ,
     %                            LTDEPO(LPCONT,JEU,NOOB), SOURCE )
                     SOURCE(1) = SOURCE(1) * PENALI
                  ELSE
C                    FIXATION(2) PENALISEE
                     MN = LTDEPO(LPCONT,JEU,NOOB)
                     CALL REFIXA( 1, NOOB, XYZ(1),XYZ(2),XYZ(3), MN,
     %                            NBCOFI, SOURCE )
C                    NBCOFI LE NOMBRE DE COMPOSANTES FIXEES
                     DO I = 1, NBCOFI
C                       LE NUMERO DE LA COMPOSANTE FIXEE
                        NU = MCN( MN + WUCOFI - 1 + I )
                        IF( NU .EQ. NOPART ) THEN
                           SOURCE(1) = SOURCE(I) * PENALI
                        ENDIF
                     ENDDO
                  ENDIF
C
C                 SI ELEMENT AXISYMETRIQUE DELTA * 2 * PI * R
                  IF(NOAXIS .NE. 0) SOURCE(1)= SOURCE(1) * D2PI * X(K,1)
C
C                 LE COEFFICIENT DU SECOND MEMBRE EST IMPOSE
                  BE(K) = BE(K) + SOURCE(1)
C
               ENDIF
            ENDIF
         ENDDO
      ENDIF
C
      RETURN
      END
