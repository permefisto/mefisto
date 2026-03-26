      SUBROUTINE TS3LAG( NOEF,   X,      PENALI, NBJEUX, JEU,
     %                   NBSOMT, NOOBPS, NUMIPO, NUMAPO, LTDEPO,
     %                   NBCOTE, NOOBLA, NUMILI, NUMALI, LTDELI,
     %                   NBPOTR, NPITR,  POIDTR, POLYTR, DPOLTR,
     %                   NBPOQU, NPIQU,  POIDQU, POLYQU, DPOLQU,
     %                   NBFACE, NOOBSF, NUMISU, NUMASU, LTDESU,
     %                   NBPOLY, NPI,    POLY,
     %                   NOOBVO, NUMIVO, NUMAVO, LTDEVO,
     %                   NBNOVI, VITEGt,
     %                   F,      POIDEL, DP,
     %                   NOPART, BE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DU SECOND MEMBRE DES EF TRIDIMENSIONNELS
C -----    LAGRANGE ISOPARAMETRIQUES
C
C ENTREES:
C --------
C NOEF   : NUMERO DU TYPE DE L'EF
C X      : COORDONNEES X Y Z DES NBPOLY POINTS DE L'EF
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
C NBCOTE : NOMBRE DES COTES DE L ELEMENT FINI
C NOOBLA : NUMERO DES OBJETS LIGNES DES ARETES DE L'ELEMENT
C NUMILI : NUMERO MINIMAL DES OBJETS LIGNES
C NUMALI : NUMERO MAXIMAL DES OBJETS LIGNES
C LTDELI : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES THERMIQUES DES LIGNES
C
C NBPOTR : NOMBRE DE POLYNOMES DE BASE SUR UNE FACE TRIANGULAIRE DE L'EF
C NPITR  : NOMBRE DE POINTS D INTEGRATION SUR UNE FACE TRIANGULAIRE
C POIDTR : POIDS DES POINTS D INTEGRATION SUR UNE FACE TRIANGULAIRE
C POLYTR : VALEUR DES POLYNOMES DE BASE AUX POINTS D INTEGRATION FACE TR
C DPOLTR : VALEUR DU GRADIENT DE CES MEMES POLYNOMES AUX POINTS ...
C
C NBPOQU : NOMBRE DE POLYNOMES DE BASE SUR UNE FACE QUADRANGULAIRE DE L'EF
C NPIQU  : NOMBRE DE POINTS D INTEGRATION SUR UNE FACE QUADRANGULAIRE
C POIDQU : POIDS DES POINTS D INTEGRATION SUR UNE FACE QUADRANGULAIRE
C POLYQU : VALEUR DES POLYNOMES DE BASE AUX POINTS D INTEGRATION FACE QU
C DPOLQU : VALEUR DU GRADIENT DE CES MEMES POLYNOMES AUX POINTS ...
C
C NBFACE : NOMBRE DE FACES DE L'EF
C NOOBSF : NUMERO DES SURFACES DES FACES DE L'EF
C NUMISU : NUMERO MINIMAL DES SURFACES UTILISEES
C NUMASU : NUMERO MAXIMAL DES SURFACES UTILISEES
C LTDESU : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES DE THERMIQUE DES SURFACES
C
C NBPOLY : NOMBRE DE POLYNOMES
C NPI    : NOMBRE DE POINTS D INTEGRATION NUMERIQUE SUR LE RECTANGLE
C POLY   : VALEUR DES POLYNOMES AUX POINTS D'INTEGRATION DE L'EF
C
C NOOBVO : NUMERO DE L'OBJET SURFACE DE CET ELEMENT
C NUMIVO : NUMERO MINIMAL DES OBJETS SURFACES
C NUMAVO : NUMERO MAXIMAL DES OBJETS SURFACES
C LTDEVO : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES DE THERMIQUE
C          DES OBJETS SURFACES
C
C F      : XX YY ZZ DES NPI POINTS D'INTEGRATION
C POIDEL : DELTA * POIDS(NPI) DES NPI POINTS D INTEGRATI
C DP     : DP(3, NBPOLY, NPI) GRADIENT AUX POINTS D 'INTEGRATION DES
C          FONCTIONS DE BASE LAGRANGE ISOPARAMETRIQUES
C X      : XX YY ZZ DES NBPOLY POINTS DE L'EF
C NOPART : POUR NLSE SEULEMENT AU NIVEAU DE SOURCE=FORCE et CONTACT=FIXATION
C          1 SI PARTIE REELLE TRAITEE ou 2 SI PARTIE IMAGINAIRE TRAITEE
C          0 SI INACTIF (CAS THERMIQUE STANDARD D'UNE SOURCE)
C
C NBNOVI : NOMBRE DE NOEUDS VITESSE DU MAILLAGE
C VITEGt : VECTEUR(NBNOVI,3) DES 3 COMPOSANTES DE LA VITESSE EN TOUS
C          LES NOEUDS DU MAILLAGE au TEMPS TEMPS
C          CE TABLEAU EST UTILISE SEULEMENT EN CAS d'UN TRANSPORT
C          DE TEMPERATURE A CETTE VITESSE POUR LE SECOND MEMBRE
C          terme: -(VITEGt . Grad) Temperature

C SORTIE :
C --------
C BE     : BE(NBPOLY) LE SECOND MEMBRE ELEMENTAIRE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     OCTOBRE 1994
C MODIFS : ALAIN PERRONNET TEXAS A & M University at QATAR  FEVRIER 2011
C MODIFS : ALAIN PERRONNET Saint Pierre du Perray             Avril 2022
C23456---------------------------------------------------------------012
      include"./incl/donthe.inc"
      include"./incl/ponoel.inc"
      include"./incl/cthet.inc"
      include"./incl/cnonlin.inc"
      include"./incl/a___fixation.inc"
      include"./incl/a___vitessefluide.inc"
C
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE      (MCN(1), RMCN(1), DMCN(1))
C
      DOUBLE PRECISION  POIDTR(NPITR),
     %                  POLYTR(NBPOTR, NPITR),
     %                  DPOLTR(2,NBPOTR, NPITR)
      DOUBLE PRECISION  POIDQU(NPIQU),
     %                  POLYQU(NBPOQU, NPIQU),
     %                  DPOLQU(2,NBPOQU, NPIQU)
      DOUBLE PRECISION  SOURCE(3),
     %                  POLY(NBPOLY, NPI),
     %                  POIDEL(NPI),
     %                  DP(3, NBPOLY, NPI)
      DOUBLE PRECISION  VITEGt(NBNOVI,3),
     %                  F(NPI,3),
     %                  BE(NBPOLY),
     %                  VITEFL(3), VITEBL(3), VITEF(27,3),
     %                  PENALI,
     %                  PROSCD,
     %                  XYZ(3), D
      REAL              X( NBPOLY, 3 )
      INTEGER           NOOBSF(NBFACE), NOOBLA(NBCOTE), NOOBPS(NBSOMT)
      INTEGER           LTDEPO(1:MXDOTH, 1:NBJEUX, NUMIPO:NUMAPO)
      INTEGER           LTDELI(1:MXDOTH, 1:NBJEUX, NUMILI:NUMALI)
      INTEGER           LTDESU(1:MXDOTH, 1:NBJEUX, NUMISU:NUMASU)
      INTEGER           LTDEVO(1:MXDOTH, 1:NBJEUX, NUMIVO:NUMAVO)
      INTEGER           NONOFK(8)
C
C     MISE A ZERO DE BE LE VECTEUR ELEMENTAIRE
C     ----------------------------------------
      CALL AZEROD( NBPOLY , BE )
C
      IF( TESTNL .GE. 1 .AND. MNTHET .GT. 0 ) THEN
C        LES SOURCES DEPENDENT DE LA TEMPERATURE
C        RECUPERATION DE LA TEMPERATURE AUX NBPOLY DL DE L ELEMENT FINI
         MN  = (MNTHET-1)/2
         MNT = (MNTHDL-1)/2
         DO I=1,NBPOLY
            DMCN( MNT+I )=DMCN( MN+MCN(MNNODL+(I-1)) )
         ENDDO
      ENDIF
      TEMPEL = 0D0
      ONDEPI = 0D0

C     ===================================
C     CONTRIBUTION DES SOURCES VOLUMIQUES
C     ===================================
      IF( LTDEVO(LPSOUR,JEU,NOOBVO) .GT. 0 ) THEN

         MNTHDLD = (MNTHDL+1 ) / 2
         DO L=1,NPI
C
C           LA CONTRIBUTION DES SOURCES DE CHALEUR
            IF( TESTNL .GE. 1 .AND. MNTHET .GT. 0 ) THEN
C              CALCUL DE LA TEMPERATURE AU POINT D'INTEGRATION L
               TEMPEL=PROSCD( POLY(1,L), DMCN(MNTHDLD), NBPOLY )
            ENDIF
C
C           LA VALEUR DES SOURCES VOLUMIQUES EN CE POINT L
            XYZ(1) = F(L,1)
            XYZ(2) = F(L,2)
            XYZ(3) = F(L,3)
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
            DO I=1,NBPOLY
               BE(I) = BE(I) + POIDEL(L) * POLY(I,L) * SOURCE(1)
            ENDDO
C
         ENDDO
C
      ENDIF
C
C     LA CONTRIBUTION DU TRANSPORT: - VITESSE * GRADIENT TEMPERATURE
C     ACTUELLEMENT CE TERME EST SOUSTRAIT DU SECOND MEMBRE + POINT FIXE
C     -----------------------------------------------------------------
      MNVIFL = LTDESU(LPVIFL,JEU,NOOBVO)
      IF( MNVIFL .GT. 0 ) THEN

C        TYPE DE LA DONNEE DE LA VITESSE DU FLUIDE
         LTVIFL = MCN( MNVIFL + WTVIFL )
C
C        RECUPERATION DE LA TEMPERATURE AUX NBPOLY DL DE L ELEMENT FINI
         MN  = (MNTHET-1)/2
         MNT = (MNTHDL-1)/2
         DO I=1,NBPOLY
            DMCN( MNT+I ) = DMCN( MN+MCN(MNNODL+I-1) )
         ENDDO

         IF( LTVIFL .LT. 2 ) THEN

C           DONNEE: VITESSE du FLUIDE CONSTANTE ou FONCTION UTILISATEUR
C           -----------------------------------------------------------
            MNTHDLD = (MNTHDL+1)/2
            DO L=1,NPI
C
C              LA VALEUR DE LA VITESSE DU FLUIDE AU POINT D'INTEGRATION L
               IF( TESTNL .GE. 1 .AND. MNTHET .GT. 0 ) THEN
C                 CALCUL DE LA TEMPERATURE AU POINT D'INTEGRATION L
                  TEMPEL=PROSCD( POLY(1,L), DMCN(MNTHDLD), NBPOLY )
               ENDIF

               XYZ(1) = F(L,1)
               XYZ(2) = F(L,2)
               XYZ(3) = F(L,3)
               CALL REVIFL( 4,NOOBVO, 3, 3, XYZ,
     %                      LTDEVO(LPVIFL,JEU,NOOBVO), VITEFL )
C
C              - t[P] [V] [DP] [TEMPERATURE EF]
               DO I=1,NBPOLY
                  D = 0D0
                  DO K=1,NBPOLY
                     D = D + ( VITEFL(1)*DP(1,K,L) +
     %                         VITEFL(2)*DP(2,K,L) +
     %                         VITEFL(3)*DP(3,K,L) ) * DMCN(MNT+K)
                  ENDDO
                  BE(I) = BE(I) - POIDEL(L) * POLY(I,L) * D
               ENDDO
            ENDDO
C
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
               VITEBL(3) = PROSCD( POLY(1,L), VITEF(1,3), NBPOLY )

C              - t[P] [V] [DP] [TEMPERATURE EF]
               DO I=1,NBPOLY

                  D = 0D0
                  DO K=1,NBPOLY
                     D = D + ( VITEBL(1)*DP(1,K,L) +
     %                         VITEBL(2)*DP(2,K,L) +
     %                         VITEBL(3)*DP(3,K,L) ) * DMCN(MNT+K)
                  ENDDO
                  BE(I) = BE(I) - POLY(I,L) * D * POIDEL(L)
C
               ENDDO

            ENDDO

         ENDIF

      ENDIF
C
C     =======================================
C     CONTRIBUTIONS DES SOURCES SUR LES FACES
C     =======================================
      DO 50 K=1,NBFACE
C
C        NO DE SURFACE DE LA FACE K
         NOOB = NOOBSF(K)
         IF( NOOB .GT. 0 ) THEN
C
C           LA FACE EST SUR UNE SURFACE SUPPORT D'ECHANGE OU CONTACT PENALISE?
            IECHAN = 0
            IF( LTDESU(LPSOUR,JEU,NOOB) .GT. 0 ) IECHAN = 1
            IF( LTDESU(LPCONT,JEU,NOOB) .GT. 0  .AND.
     %          PENALI .NE. 0D0            ) IECHAN = 2
C
            IF( IECHAN .EQ. 0 ) GOTO 50
C
C           UN TABLEAU SOURCE OU CONTACT EXISTE POUR CETTE FACE K
C           CALCUL DE LA CONTRIBUTION DE LA FACE K A BE
C           .....................................................
C           RECHERCHE DU NUMERO DES POINTS=NOEUDS DE LA FACE K
            CALL ELNOFA( NOEF, K, NBNOFK, NONOFK )
C           NONOFK(I) LES NUMEROS DES NOEUDS DE LA FACE K
C
            IF( NBSOFA(K) .EQ. 3 ) THEN
C
C              FACE TRIANGULAIRE
               CALL T43LAG( NBPOLY, X,      NBNOFK, NONOFK,
     &                      POLYTR, DPOLTR, NPITR,  POIDTR,
     &                      NOOB,   NUMISU, NUMASU, NBJEUX, JEU, LTDESU,
     &                      IECHAN, PENALI,
     &                      NOPART, BE )
            ELSE
C
C              FACE QUADRANGULAIRE
               CALL T43LAG( NBPOLY, X,      NBNOFK, NONOFK,
     &                      POLYQU, DPOLQU, NPIQU,  POIDQU,
     &                      NOOB,   NUMISU, NUMASU, NBJEUX, JEU, LTDESU,
     &                      IECHAN, PENALI,
     &                      NOPART, BE )
C
            ENDIF
         ENDIF
 50   ENDDO
C
C     ====================================================
C     CONTRIBUTION DES ARETES A LA PENALISATION DU CONTACT
C     ====================================================
      IF( PENALI .NE. 0D0 ) THEN
         DO K=1,NBCOTE
C
C           NO DE LIGNE DE L'ARETE K
            NOOB = NOOBLA(K)
            IF( NOOB .GT. 0 ) THEN
C
C              LE COTE K EST SUR UNE LIGNE
C              EST IL SUPPORT D'UN CONTACT PENALISE
               MN = LTDELI(LPCONT,JEU,NOOB)
               IF( MN .GT. 0 ) THEN
C
C                 OUI: UN TABLEAU CONTACT PENALISE EXISTE POUR CE COTE
C                 RECHERCHE DES NUMEROS LOCAUX DES NOEUDS DE L'ARETE K
                  NONOFK(1) = NOSOAR(1,K)
                  NONOFK(2) = NOSOAR(2,K)
                  IF( NBNOAR(K) .GT. 0 ) NONOFK(3) = NONOAR(1,K)
C
                  DO I=1,2+NBNOAR(K)
C
C                    LE NUMERO DU I-EME NOEUD DE L'ARETE K
                     NI = NONOFK(I)
C
C                    IL EXISTE UN CONTACT SUR CETTE LIGNE
                     XYZ(1) = X( NI, 1 )
                     XYZ(2) = X( NI, 2 )
                     XYZ(3) = X( NI, 3 )
                     IF( TESTNL .LE. 5 ) THEN
C                       CONTACT PENALISE = PENALI x TEMPERATURE
C                       LE PREMIER SOMMET DE L'ARETE
                        CALL RECONT( 2, NOOB, 3, XYZ,  MN, SOURCE )
                        SOURCE(1) = SOURCE(1) * PENALI
                     ELSE
C                       FIXATION(2) PENALISEE
                        CALL REFIXA( 2, NOOB, XYZ(1),XYZ(2),XYZ(3), MN,
     %                               NBCOFI, SOURCE )
C                       NBCOFI LE NOMBRE DE COMPOSANTES FIXEES
                        DO M = 1, NBCOFI
C                          LE NUMERO DE LA COMPOSANTE FIXEE
                           NU = MCN( MN + WUCOFI - 1 + M )
                           IF( NU .EQ. NOPART ) THEN
                              SOURCE(1) = SOURCE(M) * PENALI
                           ENDIF
                        ENDDO
                     ENDIF

C                    CALCUL DU CONTACT PENALISE = PENALI x TEMPERATURE
                     BE(NI) = BE(NI) + PENALI * SOURCE(1)

                  ENDDO

               ENDIF
            ENDIF
         ENDDO

C        =====================================================
C        CONTRIBUTION DES SOMMETS A LA PENALISATION DU CONTACT
C        =====================================================
         DO K=1,NBSOMT
C
C           NO DE POINT DU SOMMET K
            NOOB = NOOBPS(K)
            IF( NOOB .GT. 0 ) THEN
C
C              LE SOMMET K EST UN POINT. EST IL SUPPORT D'UN CONTACT PENALISE?
               MN = LTDEPO( LPCONT, JEU, NOOB )
               IF( MN .GT. 0 ) THEN
C
C                 OUI: UN TABLEAU CONTACT PENALISE EXISTE POUR CE POINT
C                 LES SOMMETS SONT NUMEROTES EN PREMIER
C                 PUIS, VIENNENT LES EVENTUELS MILIEUX DES COTES
                  IF( TESTNL .GE. 1 .AND. MNTHET .GT. 0 ) THEN
C                    PB NON LINEAIRE:
C                    RECUPERATION DE LA TEMPERATURE AU SOMMET K
                     TEMPEL = DMCN( (MNTHET-1)/2 + MCN(MNNODL+K-1) )
                  ENDIF
C
                  XYZ(1) = X( K, 1 )
                  XYZ(2) = X( K, 2 )
                  XYZ(3) = X( K, 3 )
                  IF( TESTNL .LE. 5 ) THEN
C                    CONTACT PENALISE = SOURCE(1)
                     CALL RECONT( 1, NOOB, 3, XYZ, MN, SOURCE )
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
C                 CALCUL DU CONTACT PENALISE = PENALI x TEMPERATURE
                  BE(K) = BE(K) + PENALI * SOURCE(1)
C
               ENDIF
            ENDIF
         ENDDO
      ENDIF
C
      RETURN
      END
