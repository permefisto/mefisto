      SUBROUTINE TS2P1D( X,      PENALI, NBJEUX, JEU,
     %                   NOOBPS, NUMIPO, NUMAPO, LTDEPO,
     %                   NOOBLA, NUMILI, NUMALI, LTDELI,
     %                   NOOBSF, NUMISU, NUMASU, LTDESU,
     %                   NOPART, BE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DU SECOND MEMBRE DU TRIANGLE 2P1D LAGRANGE DE DEGRE 1
C -----    INTEGRATION NUMERIQUE AUX SOMMETS DU TRIANGLE ET DES ARETES
C
C ENTREES:
C --------
C X      : COORDONNEES X ET Y DES 3 SOMMETS DU TRIANGLE
C PENALI : COEFFICIENT DE PENALISATION DE LA CONDITION DE DIRICHLET
C NBJEUX : NOMBRE DE JEUX DE DONNEES
C JEU    : NUMERO DU JEU  DE DONNEES POUR CE CALCUL DE LA MATRICE ELEMENTAIRE
C
C NOOBPS : NUMERO DE POINT DES SOMMETS
C NUMIPO : NUMERO MINIMAL DES OBJETS POINTS
C NUMAPO : NUMERO MAXIMAL DES OBJETS POINTS
C LTDEPO : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES THERMIQUES DES POINTS
C
C NOOBLA : NOUMERO DE LA LIGNE DES 3 ARETES DU TRIANGLE
C SOURCE : FLUX NORMAL DE LA TEMPERATURE A LA PAROI
C NUMILI : NUMERO MINIMAL DES OBJETS LIGNES
C NUMALI : NUMERO MAXIMAL DES OBJETS LIGNES
C LTDELI : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES DE THERMIQUES DES LIGNES
C          DES OBJETS LIGNES
C
C NOOBSF : NUMERO DE L'OBJET SURFACE DE CET ELEMENT
C SOURCE : SOURCES DE SURFACE
C          TABLEAU REMPLI DANS CE SOUS PROGRAMME
C NUMISU : NUMERO MINIMAL DES OBJETS SURFACES
C NUMASU : NUMERO MAXIMAL DES OBJETS SURFACES
C LTDESU : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES DE THERMIQUE DES SURFACES
C          DES OBJETS SURFACES
C NOPART : POUR NLSE SEULEMENT AU NIVEAU DE SOURCE=FORCE et CONTACT=FIXATION
C          1 SI PARTIE REELLE TRAITEE ou 2 SI PARTIE IMAGINAIRE TRAITEE
C          0 SI INACTIF (CAS THERMIQUE STANDARD D'UNE SOURCE)
C
C SORTIE :
C --------
C BE     : BE(3) LE SECOND MEMBRE ELEMENTAIRE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 1994
C23456---------------------------------------------------------------012
      include"./incl/donthe.inc"
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
      DOUBLE PRECISION  PENALI,
     %                  SOURCE(2), VN(2),
     %                  BE(3)
      REAL              X(3, 2)
      INTEGER           NOOBPS(3)
      INTEGER           NOOBLA(3)
      INTEGER           LTDEPO(1:MXDOTH, 1:NBJEUX, NUMIPO:NUMAPO)
      INTEGER           LTDELI(1:MXDOTH, 1:NBJEUX, NUMILI:NUMALI)
      INTEGER           LTDESU(1:MXDOTH, 1:NBJEUX, NUMISU:NUMASU)
      DOUBLE PRECISION  X21, Y21, X13, Y13, X32, Y32, DELTA, XYZ(3)
      DOUBLE PRECISION  D, VITEFL(2,3), TEMPRE(3)
      INTEGER           K12(2)
      EQUIVALENCE      (K12(1),K1), (K12(2),K2)
C
C     MISE A ZERO DU VECTEUR SECOND MEMBRE
      ONDEPI = 0D0
      BE(1) = 0D0
      BE(2) = 0D0
      BE(3) = 0D0
C
C     ====================================
C     CONTRIBUTION DES SOURCES SURFACIQUES
C     ====================================
      IF( LTDESU(LPSOUR,JEU,NOOBSF) .GT. 0 ) THEN
C
         X21 = X(2,1) - X(1,1)
         X32 = X(3,1) - X(2,1)
         X13 = X(1,1) - X(3,1)
C
         Y21 = X(2,2) - X(1,2)
         Y32 = X(3,2) - X(2,2)
         Y13 = X(1,2) - X(3,2)
C
         DELTA = ABS( X13 * Y21 - X21 * Y13 ) / 6D0
C
         DO 5 K=1,3
C
C           LE POINT D'INTEGRATION K EST LE SOMMET K
            IF( TESTNL .GE. 1 .AND. MNTHET .GT. 0 ) THEN
               MN = (MNTHET-1)/2
               TEMPEL = DMCN( MN+MCN(MNNODL+K-1) )
            ELSE
               TEMPEL = 0D0
            ENDIF
C
C           LA VALEUR DES SOURCES SURFACIQUES EN CE POINT
            XYZ(1) = X(K,1)
            XYZ(2) = X(K,2)
            XYZ(3) = 0D0
C
            IF( TESTNL .LE. 5 ) THEN
C              LA VALEUR DES SOURCES SURFACIQUES EN CE SOMMET K
               CALL RESOUR( 3, NOOBSF, 3, XYZ,
     %                      LTDESU(LPSOUR,JEU,NOOBSF), SOURCE )
            ELSE
               CALL REFORC( 3, NOOBSF, 2, XYZ(1),XYZ(2),XYZ(3),
     %                                    0D0,   0D0,   0D0,
     %                      LTDESU(LPSOUR,JEU,NOOBSF), SOURCE )
               SOURCE(1) = SOURCE(NOPART)
            ENDIF
C
C           LA CONTRIBUTION DE CETTE SOURCE DE CHALEUR
            BE(K) = SOURCE(1) * DELTA
 5       CONTINUE
C
      ENDIF
C
C     CONTRIBUTION EVENTUELLE DE - VITESSE * GRADIENT TEMPERATURE
C     ACTUELLEMENT CE TERME EST SOUSTRAIT DU SECOND MEMBRE + POINT FIXE
C     -----------------------------------------------------------------
      IF( LTDESU(LPVIFL,JEU,NOOBSF) .GT. 0 ) THEN
C        CALCUL DE LA TEMPERATURE  AUX 3 SOMMETS DU TRIANGLE
C        CALCUL DU VECTEUR VITESSE AUX 3 SOMMETS DU TRIANGLE
C
         DO 15 K=1,3
            IF( TESTNL .GE. 1 .AND. MNTHET .GT. 0 ) THEN
               MN = (MNTHET-1)/2
               TEMPRE(K) = DMCN( MN+MCN(MNNODL+K-1) )
            ELSE
               TEMPRE(K) = 0D0
            ENDIF
            TEMPEL = TEMPRE(K)
C           LES 2 COORDONNEES DU SOMMET K DU TRIANGLE
            XYZ(1) = X(K,1)
            XYZ(2) = X(K,2)
            XYZ(3) = 0D0
            CALL REVIFL( 3, NOOBSF, 2, 3, XYZ,
     %                   LTDESU(LPVIFL,JEU,NOOBSF), VITEFL(1,K) )
 15      CONTINUE
C
C        CONTRIBUTION AU SECOND MEMBRE
         X21 = X(2,1) - X(1,1)
         X32 = X(3,1) - X(2,1)
         X13 = X(1,1) - X(3,1)
C
         Y21 = X(2,2) - X(1,2)
         Y32 = X(3,2) - X(2,2)
         Y13 = X(1,2) - X(3,2)
C
C        SOUSTRACTION de t[P] t{V} [DP] {Temp n-1}
         DO 30 K=1,3
            BE(K)= BE(K)
     %           -( (X32*VITEFL(2,K)-Y32*VITEFL(1,K)) * TEMPRE(1)
     %             +(X13*VITEFL(2,K)-Y13*VITEFL(1,K)) * TEMPRE(2)
     %             +(X21*VITEFL(2,K)-Y21*VITEFL(1,K)) * TEMPRE(3) )/6D0
 30      CONTINUE
C
      ENDIF
C
C     =======================================
C     CONTRIBUTIONS DES SOURCES SUR LES COTES
C     =======================================
      DO 50 K=1,3
C
C        NO DE LIGNE DE L'ARETE K
         NOOB = NOOBLA(K)
         IF( NOOB .GT. 0 ) THEN
C
C           LE COTE K EST SUR UNE LIGNE. SUPPORT DE FLUX OU CONTACT PENALISE?
            IECHAN = 0
            IF( LTDELI(LPSOUR,JEU,NOOB) .GT. 0      ) IECHAN = 1
            IF( LTDELI(LPCONT,JEU,NOOB) .GT. 0 .AND.
     %          PENALI .NE. 0D0                     ) IECHAN = 2
C
            IF( IECHAN .EQ. 0 ) GOTO 50
C
C           UN TABLEAU SOURCE OU CONTACT PENALISE EXISTE POUR CETTE LIGNE
C           CALCUL DE LA CONTRIBUTION DE L'ARETE K A BE
C           .............................................................
C           LE NUMERO DES POINTS DU COTE K DANS LE TRIANGLE
            K1 = K
            IF( K1 .NE. 3 ) THEN
               K2 = K1 + 1
            ELSE
               K2 = 1
            ENDIF
C
C           LE VECTEUR ORTHOGONAL A L'ARETE K
            VN(1) = X(K2,2) - X(K ,2)
            VN(2) = X(K ,1) - X(K2,1)
C
C           LE JACOBIEN
            D = SQRT( VN(1)**2 + VN(2)**2 )
C
C           LE VECTEUR NORMAL UNITAIRE
            VN(1) = VN(1) / D
            VN(2) = VN(2) / D
C
C           LA LONGUEUR DE L'ARETE K / 2
            DELTA = SQRT( (X(K2,1)-X(K,1)) ** 2
     %                  + (X(K2,2)-X(K,2)) ** 2 ) * 0.5D0
C
C           CALCUL DU FLUX NORMAL AU SOMMET L DE L'ARETE K
C           ..............................................
            DO L = 1, 2
C
C              NO DE 1 A 3 DU SOMMET L DE L'ARETE K
               NSLK = K12(L)
C
               IF( TESTNL .GE. 1 .AND. MNTHET .GT. 0 ) THEN
                  MN = (MNTHET-1)/2
                  TEMPEL = DMCN( MN+MCN(MNNODL+NSLK-1) )
               ELSE
                  TEMPEL = 0D0
               ENDIF
C
C              LES 3 COORDONNEES DU SOMMET L DE L'ARETE K
               XYZ(1) = X(NSLK,1)
               XYZ(2) = X(NSLK,2)
               XYZ(3) = 0D0
C
               IF( IECHAN .EQ. 1 ) THEN
C
                  IF( TESTNL .LE. 5 ) THEN
C                    SOURCE(1) REQUISE
                     CALL RESOUR( 2, NOOB, 3, XYZ,
     %                            LTDELI(LPSOUR,JEU,NOOB), SOURCE )
                  ELSE
C                    FORCE(2) REQUISE
                     CALL REFORC( 2, NOOB, 2, XYZ(1),XYZ(2),XYZ(3),
     %                                        VN(1), VN(2), 0D0,
     %                            LTDELI(LPSOUR,JEU,NOOB), SOURCE )
                     SOURCE(1) = SOURCE(NOPART)
                  ENDIF
C
               ELSE
C
C                 CONTACT PENALISE = PENALI x TEMPERATURE
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
C              LA CONTRIBUTION AU SECOND MEMBRE DU FLUX NORMAL (NEUMANN)
               BE(NSLK) = BE(NSLK) + SOURCE(1) * DELTA
            ENDDO
C
         ENDIF
 50   CONTINUE
C
C     ============================================
C     CONTRIBUTION DU CONTACT PENALISE AUX SOMMETS
C     ============================================
      IF( PENALI .NE. 0D0 ) THEN
         DO 80 K=1,3
C
C           NO DE POINT DU SOMMET K
            NOOB = NOOBPS(K)
            IF( NOOB .GT. 0 ) THEN
C
C              LE SOMMET K EST UN POINT. EST IL SUPPORT D'UN CONTACT PENALISE?
               MNLT = LTDEPO(LPCONT,JEU,NOOB)
               IF( MNLT .GT. 0 ) THEN
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
C
C                 3 COORDONNEES DU SOMMET K
                  XYZ(1) = X(K,1)
                  XYZ(2) = X(K,2)
                  XYZ(3) = 0D0
C
                  IF( TESTNL .LE. 5 ) THEN
C                    CONTACT PENALISE = SOURCE(1)
                     CALL RECONT( 1, NOOB, 3, XYZ, MNLT, SOURCE )
                     SOURCE(1) = SOURCE(1) * PENALI
                  ELSE
C                    FIXATION(2) PENALISEE
                     CALL REFIXA( 1, NOOB, XYZ(1),XYZ(2),XYZ(3), MNLT,
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
 80      CONTINUE
      ENDIF
C
      RETURN
      END
