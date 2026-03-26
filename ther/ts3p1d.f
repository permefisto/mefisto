      SUBROUTINE TS3P1D( X,      DELTA,  DP,     PENALI, NBJEUX, JEU,
     &                   NOOBPS, NUMIPO, NUMAPO, LTDEPO,
     &                   NOOBLA, NUMILI, NUMALI, LTDELI,
     &                   NOOBSF, NUMISU, NUMASU, LTDESU,
     &                   NOOBVO, NUMIVO, NUMAVO, LTDEVO,
     &                   NOPART, BE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DU SECOND MEMBRE D'UN TETRAEDRE 3P1D
C -----    INTEGRATION NUMERIQUE AUX SOMMETS DU TETRAEDRE ET DES FACES
C
C ENTREES:
C --------
C X      : COORDONNEES RAYON ET COTE DES 4 POINTS DE L'EF
C          OU X Y DES 4 POINTS DE L'EF
C DELTA  : JACOBIEN DE LA TRANSFORMATION EF REFERENCE -> EF
C DP     : GRADIENT DES POLYNOMES DE BASE AUX SOMMETS DU TETRAEDRE
C PENALI : COEFFICIENT DE PENALISATION DE LA CONDITION DE DIRICHLET
C NBJEUX : NOMBRE DE JEUX DE DONNEES
C JEU    : NUMERO DU JEU  DE DONNEES POUR CE CALCUL DE LA MATRICE ELEMENTAIRE
C
C NOOBPS : NUMERO DE POINT DES SOMMETS
C NUMIPO : NUMERO MINIMAL DES OBJETS POINTS
C NUMAPO : NUMERO MAXIMAL DES OBJETS POINTS
C LTDEPO : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES THERMIQUES DES POINTS
C
C NOOBLA : NUMERO DES OBJETS LIGNES DES ARETES DE L'EF
C NUMILI : NUMERO MINIMAL DES OBJETS LIGNES
C NUMALI : NUMERO MAXIMAL DES OBJETS LIGNES
C LTDELI : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES THERMIQUES
C          DES OBJETS LIGNES
C
C NOOBSF : NUMERO DES SURFACES DES FACES DE L'EF
C NUMISU : NUMERO MINIMAL DES SURFACES UTILISEES
C NUMASU : NUMERO MAXIMAL DES SURFACES UTILISEES
C LTDESU : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES DE THERMIQUE
C          DES OBJETS SURFACES
C
C NOOBVO : NUMERO DE L'OBJET SURFACE DE CET ELEMENT
C NUMIVO : NUMERO MINIMAL DES OBJETS SURFACES
C NUMAVO : NUMERO MAXIMAL DES OBJETS SURFACES
C LTDEVO : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES DE THERMIQUE
C          DES OBJETS SURFACES
C NOPART : POUR NLSE SEULEMENT AU NIVEAU DE SOURCE=FORCE et CONTACT=FIXATION
C          1 SI PARTIE REELLE TRAITEE ou 2 SI PARTIE IMAGINAIRE TRAITEE
C          0 SI INACTIF (CAS THERMIQUE STANDARD D'UNE SOURCE)
C
C SORTIE :
C --------
C BE     : BE(4) LE SECOND MEMBRE ELEMENTAIRE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 1994
C MODIF  : ALAIN PERRONNET TEXAS A & M University at QATAR  FEVRIER 2011
C23456---------------------------------------------------------------012
      include"./incl/donthe.inc"
      include"./incl/ponoel.inc"
      include"./incl/cthet.inc"
      include"./incl/cnonlin.inc"
      include"./incl/a___fixation.inc"
C
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE      (MCN(1), RMCN(1), DMCN(1))
C
      DOUBLE PRECISION  DP(3,4), PENALI
      DOUBLE PRECISION  BE(4), VITEFL(3,4), TEMPRE(4)
      REAL              X(4,3)
      INTEGER           NOOBPS(1:4),
     %                  NOOBLA(1:6),
     %                  NOOBSF(1:4)
      INTEGER           LTDEPO(1:MXDOTH, 1:NBJEUX, NUMIPO:NUMAPO)
      INTEGER           LTDELI(1:MXDOTH, 1:NBJEUX, NUMILI:NUMALI)
      INTEGER           LTDESU(1:MXDOTH, 1:NBJEUX, NUMISU:NUMASU)
      INTEGER           LTDEVO(1:MXDOTH, 1:NBJEUX, NUMIVO:NUMAVO)
      INTEGER           NONOFK(3)
      DOUBLE PRECISION  SOURCE(3),VN(3), DELTA, DELTAK, XYZ(3),
     %                  DGL(2,3), DGLN
C
C     MISE A ZERO DE BE LE VECTEUR ELEMENTAIRE
C     ---------------------------------------
      ONDEPI = 0D0
      BE(1) = 0D0
      BE(2) = 0D0
      BE(3) = 0D0
      BE(4) = 0D0
C
C     ===================================
C     CONTRIBUTION DES SOURCES VOLUMIQUES
C     ===================================
      IF( LTDEVO(LPSOUR,JEU,NOOBVO) .GT. 0 ) THEN
C
         DO 15 K=1,4
C
C           LE POINT D'INTEGRATION K EST LE SOMMET K
            IF( TESTNL .GE. 1 .AND. MNTHET .GT. 0 ) THEN
C              LA TEMPERATURE AU SOMMET K DU TETRAEDRE
               MN = (MNTHET-1)/2
               TEMPEL = DMCN( MN+MCN(MNNODL+K-1) )
            ELSE
               TEMPEL = 0D0
            ENDIF
C
C           LA VALEUR DES SOURCES VOLUMIQUES EN CE SOMMET K
            XYZ(1) = X(K,1)
            XYZ(2) = X(K,2)
            XYZ(3) = X(K,3)
            IF( TESTNL .LE. 5 ) THEN
               CALL RESOUR( 4, NOOBVO, 3, XYZ,
     %                      LTDEVO(LPSOUR,JEU,NOOBVO), SOURCE )
            ELSE
               CALL REFORC( 4, NOOBVO, 3, XYZ(1),XYZ(2),XYZ(3),
     %                                    0D0,   0D0,   0D0,
     %                      LTDESU(LPSOUR,JEU,NOOBSF), SOURCE )
               SOURCE(1) = SOURCE(NOPART)
            ENDIF
C
C           LA CONTRIBUTION DE CETTE SOURCE DE CHALEUR
            BE(K) = SOURCE(1) * DELTA / 24D0
 15      CONTINUE
C
      ENDIF
C
C     CONTRIBUTION EVENTUELLE DE - VITESSE * GRADIENT TEMPERATURE
C     ACTUELLEMENT CE TERME EST SOUSTRAIT DU SECOND MEMBRE + POINT FIXE
C     -----------------------------------------------------------------
      IF( LTDEVO(LPVIFL,JEU,NOOBVO) .GT. 0 ) THEN
C
C        CALCUL DE LA TEMPERATURE  AUX 4 SOMMETS DU TETRAEDRE
C        CALCUL DU VECTEUR VITESSE AUX 4 SOMMETS DU TETRAEDRE
         DO 25 K=1,4
            IF( TESTNL .GE. 1 .AND. MNTHET .GT. 0 ) THEN
               MN = (MNTHET-1)/2
               TEMPRE(K) = DMCN( MN+MCN(MNNODL+K-1) )
            ELSE
               TEMPRE(K) = 0D0
            ENDIF
            TEMPEL = TEMPRE(K)
C           LES 3 COORDONNEES DU SOMMET K DU TETRAEDRE
            XYZ(1) = X(K,1)
            XYZ(2) = X(K,2)
            XYZ(3) = X(K,3)
            CALL REVIFL( 4, NOOBVO, 3, 3, XYZ,
     %                   LTDEVO(LPVIFL,JEU,NOOBVO), VITEFL(1,K) )
 25      CONTINUE
C
C        CONTRIBUTION AU SECOND MEMBRE
         DO 30 K=1,4
               BE(K) = BE(K) - (
     %     (VITEFL(1,K)*DP(1,1)+VITEFL(2,K)*DP(2,1)+VITEFL(3,K)*DP(3,1))
     %      *TEMPRE(1)
     %    +(VITEFL(1,K)*DP(1,2)+VITEFL(2,K)*DP(2,2)+VITEFL(3,K)*DP(3,2))
     %      *TEMPRE(2)
     %    +(VITEFL(1,K)*DP(1,3)+VITEFL(2,K)*DP(2,3)+VITEFL(3,K)*DP(3,3))
     %      *TEMPRE(3)
     %    +(VITEFL(1,K)*DP(1,4)+VITEFL(2,K)*DP(2,4)+VITEFL(3,K)*DP(3,4))
     %      *TEMPRE(4)             ) * DELTA / 24D0
 30      CONTINUE
C
      ENDIF
C
C     ======================================
C     CONTRIBUTION DES SOURCES SUR LES FACES
C     ======================================
      DO 80 K=1,4
C
C        NO DE SURFACE DE LA FACE K
         NOOB = NOOBSF(K)
         IF( NOOB .GT. 0 ) THEN
C
C           LA FACE EST SUR UNE SURFACE UTILISATEUR
            IECHAN = 0
            IF( LTDESU(LPSOUR,JEU,NOOB) .GT. 0 ) IECHAN = 1
            IF( LTDESU(LPCONT,JEU,NOOB) .GT. 0  .AND.
     %          PENALI .NE. 0D0                ) IECHAN = 2
C
            IF( IECHAN .EQ. 0 ) GOTO 80
C
C           UN TABLEAU SOURCE OU CONTACT EXISTE POUR CETTE FACE K
C           CALCUL DE LA CONTRIBUTION DE LA FACE K A BE
C           .....................................................
C           RECHERCHE DU NUMERO DES POINTS=NOEUDS DE LA FACE K
            CALL ELNOFA( 19, K, NBNOFK, NONOFK )
C           NONOFK(I) LES NUMEROS DES NOEUDS DE LA FACE K
C                     CE SONT AUSSI LES POINTS D'INTEGRATION
C
C           RECHERCHE DU JACOBIEN DE G
            N = NONOFK(1)
            DGL(1,1) = X( NONOFK(2), 1 ) - X( N, 1 )
            DGL(2,1) = X( NONOFK(3), 1 ) - X( N, 1 )
            DGL(1,2) = X( NONOFK(2), 2 ) - X( N, 2 )
            DGL(2,2) = X( NONOFK(3), 2 ) - X( N, 2 )
            DGL(1,3) = X( NONOFK(2), 3 ) - X( N, 3 )
            DGL(2,3) = X( NONOFK(3), 3 ) - X( N, 3 )
            CALL JAR2R3( DGL, DELTAK )
C
C           LE VECTEUR NORMAL UNITAIRE A LA FACE
            CALL VECNOR( DGL, DGLN, VN )
C
C           CALCUL DU FLUX NORMAL
            DO 60 L=1,3
C
C              LE NUMERO N DANS LE TETRAEDRE DU SOMMET L DE LA FACE K
               N = NONOFK( L )
C
               IF( TESTNL .GE. 1 .AND. MNTHET .GT. 0 ) THEN
C                 LA TEMPERATURE AU SOMMET
                  MN = (MNTHET-1)/2
                  TEMPEL=DMCN( MN+MCN(MNNODL+N-1) )
               ELSE
                  TEMPEL = 0D0
               ENDIF
C
C              LES 3 COORDONNEES DU SOMMET L DE LA FACE K
               XYZ(1) = X( N, 1 )
               XYZ(2) = X( N, 2 )
               XYZ(3) = X( N, 3 )
C
               IF( IECHAN .EQ. 1 ) THEN
C
                  IF( TESTNL .LE. 5 ) THEN
C                    SOURCE(1) REQUISE
C                    CALCUL DU FLUX NORMAL DU SOMMET L DE LA FACE K
                     CALL RESOUR( 3, NOOB, 3, XYZ,
     %                            LTDESU(LPSOUR,JEU,NOOB), SOURCE )
                  ELSE
C                    FORCE(2) REQUISE
C                    LE VECTEUR NORMAL UNITAIRE EST UTILISE
                     CALL REFORC( 3, NOOB, 3, XYZ(1),XYZ(2),XYZ(3),
     %                                        VN(1), VN(2), VN(3),
     %                            LTDESU(LPSOUR,JEU,NOOB), SOURCE )
                     SOURCE(1) = SOURCE(NOPART)
                  ENDIF
C
               ELSE
C
                  IF( TESTNL .LE. 5 ) THEN
C
C                    CONTACT PENALISE = PENALI x TEMPERATURE
                     CALL RECONT( 3, NOOB, 3, XYZ,
     %                            LTDESU(LPCONT,JEU,NOOB), SOURCE )
                     SOURCE(1) = SOURCE(1) * PENALI
                  ELSE
C                    FIXATION(2) PENALISEE
                     MN =  LTDESU(LPCONT,JEU,NOOB)
                     CALL REFIXA( 3, NOOB, XYZ(1),XYZ(2),XYZ(3), MN,
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
C              SOMMATION AVEC LE VECTEUR ELEMENTAIRE
               BE(N) = BE(N) + SOURCE(1) * DELTAK / 6D0
 60         CONTINUE
         ENDIF
 80   CONTINUE
C
C     ===========================================
C     CONTRIBUTION DU CONTACT PENALISE AUX ARETES
C     ===========================================
      IF( PENALI .NE. 0D0 ) THEN
         DO 130 K=1,6
C
C           NO DE LIGNE DE L'ARETE K
            NOOB = NOOBLA(K)
            IF( NOOB .GT. 0 ) THEN
C
C              LE SOMMET K EST UN POINT. EST IL SUPPORT D'UN CONTACT PENALISE?
               IF( LTDELI(LPCONT,JEU,NOOB) .GT. 0 ) THEN
C
C                 OUI: UN TABLEAU CONTACT PENALISE EXISTE POUR CE POINT
C                 CONTACT PENALISE = PENALI x TEMPERATURE
C
C                 LE NUMERO DES 2 SOMMETS DE L'ARETE K
                  GOTO( 111, 112, 113, 114, 115, 116 ) , K
 111              N = 1
                  L = 2
                  GOTO 118
 112              N = 2
                  L = 3
                  GOTO 118
 113              N = 3
                  L = 1
                  GOTO 118
 114              N = 1
                  L = 4
                  GOTO 118
 115              N = 2
                  L = 4
                  GOTO 118
 116              N = 3
                  L = 4
C
 118              XYZ(1) = X(N,1)
                  XYZ(2) = X(N,2)
                  XYZ(3) = X(N,3)
                  IF( TESTNL .LE. 5 ) THEN
C                    CONTACT PENALISE = PENALI x TEMPERATURE
C                    LE PREMIER SOMMET DE L'ARETE
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
C                 LE COEFFICIENT DU SECOND MEMBRE
                  BE(N) = BE(N) + SOURCE(1)
C
C                 LE SECOND SOMMET DE L'ARETE
                  XYZ(1) = X(L,1)
                  XYZ(2) = X(L,2)
                  XYZ(3) = X(L,3)
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
C                 LE COEFFICIENT DU SECOND MEMBRE
                  BE(L) = BE(L) + SOURCE(1)
C
               ENDIF
            ENDIF
 130     CONTINUE
C
C        ============================================
C        CONTRIBUTION DU CONTACT PENALISE AUX SOMMETS
C        ============================================
         DO K=1,4
C
C           NO DE POINT DU SOMMET K
            NOOB = NOOBPS(K)
            IF( NOOB .GT. 0 ) THEN
C
C              LE SOMMET K EST UN POINT. EST IL SUPPORT D'UN CONTACT PENALISE?
               IF( LTDEPO(LPCONT,JEU,NOOB) .GT. 0 ) THEN
C
C                 OUI: UN TABLEAU CONTACT PENALISE EXISTE POUR CE POINT
C                 CONTACT PENALISE = PENALI x TEMPERATURE
                  IF( TESTNL .GE. 1 .AND. MNTHET .GT. 0 ) THEN
C                    PB NON LINEAIRE:
C                    RECUPERATION DE LA TEMPERATURE AU SOMMET K
                     TEMPEL = DMCN( (MNTHET-1)/2 + MCN(MNNODL+K-1) )
                  ELSE
                     TEMPEL = 0D0
                  ENDIF
                  XYZ(1) = X(K,1)
                  XYZ(2) = X(K,2)
                  XYZ(3) = X(K,3)
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
C                 LE COEFFICIENT DU SECOND MEMBRE
                  BE(K) = BE(K) + SOURCE(1)
C
               ENDIF
            ENDIF
         ENDDO
      ENDIF
C
      RETURN
      END
