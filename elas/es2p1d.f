      SUBROUTINE ES2P1D( X,      PENALI,
     %                   NOOBPS, NUMIPO, NUMAPO, LTDEPO,
     %                   NOOBLA, NUMILI, NUMALI, LTDELI,
     %                   NOOBSF, NUMISU, NUMASU, LTDESU,
     %                   NUELEM, NBELEM, NUNDEL,
     %                   MNTEMP, NTDLTE, TEMPER,
     %                   BE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DES SECONDS MEMBRES DU TRIANGLE 2P1D
C -----    LAGRANGE DE DEGRE 1 EN ELASTICITE 2D
C
C ENTREES:
C --------
C X      : 2 COORDONNEES DES 3 SOMMETS DU TRIANGLE
C PENALI : 1/EPSILON DE LA PENALISATION DU DEPLACEMENT IMPOSE
C
C NOOBPS : NUMERO DE POINT DES SOMMETS
C NUMIPO : NUMERO MINIMAL  DES OBJETS POINTS
C NUMAPO : NUMERO MAXIMAL  DES OBJETS POINTS
C LTDEPO : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES D'ELASTICITE
C          DES POINTS
C
C NOOBLA : NUMERO DES OBJETS LIGNES DES ARETES DE L'ELEMENT FINI
C NUMILI : NUMERO MINIMAL DES OBJETS LIGNES
C NUMALI : NUMERO MAXIMAL DES OBJETS LIGNES
C LTDELI : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES D'ELASTICITE
C          DES OBJETS LIGNES
C
C NOOBSF : NUMERO DE L'OBJET SURFACE DE CET ELEMENT FINI
C NUMISU : NUMERO MINIMAL DES OBJETS SURFACES
C NUMASU : NUMERO MAXIMAL DES OBJETS SURFACES
C LTDESU : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES D'ELASTICITE
C          DES OBJETS SURFACES
C
C NUELEM : NUMERO DE L'ELEMENT FINI
C NBELEM : NOMBRE D'ELEMENTS FINIS
C NUNDEL : NUMERO DES NOEUDS DES ELEMENTS FINIS
C
C MNTEMP : >0 CALCUL DEMANDE DES CONTRAINTES THERMIQUES,=<0 SINON
C NTDLTE : NOMBRE DE NOEUDS THERMIQUES OU DL THERMIQUES DU MAILLAGE
C TEMPER : TEMPERATURE AUX NOEUDS DU MAILLAGE TEMPER(NTDLTE)
C
C SORTIES:
C --------
C BE     : BE(6) LE SECOND MEMBRE ELEMENTAIRE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE UPMC PARIS        MAI 1998
C-----------------------------------------------------------------------
      include "./incl/donela.inc"
      include "./incl/a___fixation.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              X(3,2)
      DOUBLE PRECISION  PENALI,
     %                  FORCE(3),
     %                  ELAS(6),
     %                  AUX(3),
     %                  TEMPER(NTDLTE),
     %                  BE(6)
      INTEGER           NUNDEL(1:NBELEM,1:3)
      INTEGER           NOOBPS(3),NOOBLA(3)
      INTEGER           LTDEPO(1:MXDOEL,NUMIPO:NUMAPO)
      INTEGER           LTDELI(1:MXDOEL,NUMILI:NUMALI)
      INTEGER           LTDESU(1:MXDOEL,NUMISU:NUMASU)
      DOUBLE PRECISION  TEMPIN,DILATA
C
      INTEGER           NOSOAR(2)
      EQUIVALENCE      (NOSOAR(1),K),(NOSOAR(2),KK)
      DOUBLE PRECISION  D, DELTA, XD, YD, VN(2)
      DOUBLE PRECISION  X21, Y21, X31, Y31, X32, Y32
C
C     INITIALISATION DE BE
      CALL AZEROD( 6, BE )
C
      X21 = X(2,1) - X(1,1)
      X31 = X(3,1) - X(1,1)
      X32 = X(3,1) - X(2,1)
C
      Y21 = X(2,2) - X(1,2)
      Y31 = X(3,2) - X(1,2)
      Y32 = X(3,2) - X(2,2)
C
      DELTA = ABS( X21 * Y31 - X31 * Y21 )
C
C     CONTRIBUTION DES EFFORTS SURFACIQUES
C     ------------------------------------
      IF( LTDESU(LPFORC,NOOBSF) .GT. 0 ) THEN
         L = 1
         D = DELTA / 6D0
         DO 15 K=1,3
C           LA VALEUR DES EFFORTS SURFACIQUES EN CE SOMMET K
            XD = X(K,1)
            YD = X(K,2)
            CALL REFORC( 3, NOOBSF, 2, XD,YD,0D0, 0D0,0D0,0D0,
     %                   LTDESU(LPFORC,NOOBSF), FORCE )
            BE(L  ) = BE(L  ) + D * FORCE(1)
            BE(L+1) = BE(L+1) + D * FORCE(2)
            L = L + 2
   15    CONTINUE
      ENDIF
C
C     CONTRIBUTIONS DES EFFORTS SUR LES COTES
C     ---------------------------------------
      DO 50 K=1,3
         NL = NOOBLA(K)
         IF( NL .GT. 0 ) THEN
            IF( LTDELI(LPFORC,NL) .GT. 0 ) THEN
C
C              LE NUMERO DES SOMMETS DU COTE K
               IF( K .NE. 3 ) THEN
                  KK = K+1
               ELSE
                  KK = 1
               ENDIF
C
C              LE VECTEUR ORTHOGONAL AU VECTEUR TANGENT
               VN(1) = X(KK,2) - X(K ,2)
               VN(2) = X(K ,1) - X(KK,1)
C
C              LE JACOBIEN
               D = SQRT( VN(1)**2 + VN(2)**2 )
C
C              LE VECTEUR NORMAL UNITAIRE
               VN(1) = VN(1) / D
               VN(2) = VN(2) / D
C
C              LE PRODUIT DU JACOBIEN PAR LE POIDS
               D = D * 0.5D0
C
               DO 40 M=1,2
C
C                 LE NUMERO DU SOMMET M DE L'ARETE K
                  MM = NOSOAR(M)
C
C                 CALCUL DE LA PRESSION EN CE SOMMET MM
                  XD = X(MM,1)
                  YD = X(MM,2)
                  CALL REFORC( 2, NL, 2,
     %                         XD,YD,0D0, VN(1),VN(2),0D0,
     %                         LTDELI(LPFORC,NL), FORCE )
                  L = 2 * MM
                  BE(L-1) = BE(L-1) + D * FORCE(1)
                  BE(L  ) = BE(L  ) + D * FORCE(2)
C
 40            CONTINUE
            ENDIF
         ENDIF
 50   CONTINUE
C
C     CONTRIBUTION DES CONTRAINTES INITIALES
C     --------------------------------------
      IF( LTDESU(LPCOIN,NOOBSF) .GT. 0 ) THEN
C        MISE A ZERO DU VECTEUR SOMME DES CONTRAINTES INITIALES AUX 3 SOMMETS
         AUX(1) = 0D0
         AUX(2) = 0D0
         AUX(3) = 0D0
         DO 70 L=1,3
C           LA VALEUR DES CONTRAINTES INITIALES EN CE SOMMET L
            XD = X(L,1)
            YD = X(L,2)
            CALL RECOIN( 3,NOOBSF,3, XD,YD,0.D0,
     %                   LTDESU(LPCOIN,NOOBSF), FORCE )
            AUX(1) = AUX(1) + FORCE(1)
            AUX(2) = AUX(2) + FORCE(2)
            AUX(3) = AUX(3) + FORCE(3)
 70      CONTINUE
C
C        AJOUT AUX SECONDS MEMBRES ELEMENTAIRES
         BE(1) = BE(1) - (AUX(3) * X32 - AUX(1) * Y32) / 6D0
         BE(2) = BE(2) - (AUX(2) * X32 - AUX(3) * Y32) / 6D0
         BE(3) = BE(3) - (AUX(1) * Y31 - AUX(3) * X31) / 6D0
         BE(4) = BE(4) - (AUX(3) * Y31 - AUX(2) * X31) / 6D0
         BE(5) = BE(5) - (AUX(3) * X21 - AUX(1) * Y21) / 6D0
         BE(6) = BE(6) - (AUX(2) * X21 - AUX(3) * Y21) / 6D0
      ENDIF
C
C     CONTRIBUTIONS DES CONTRAINTES THERMIQUES
C     ----------------------------------------
      IF( MNTEMP .GT. 0 .AND.
     &    LTDESU(LPDILA,NOOBSF) .GT. 0 ) THEN
C
         AUX(1) = 0D0
         AUX(2) = 0D0
         AUX(3) = 0D0
C
         DO 120 L=1,3
C           LA VALEUR DES CONTRAINTES INITIALES EN CE SOMMET L
            XD = X(L,1)
            YD = X(L,2)
C
C           LE TENSEUR DE L'ELASTICITE EN CE POINT D'INTEGRATION
            CALL REELAS( 3,NOOBSF,2, XD,YD,0D0,
     &                   LTDESU(LPYOUN,NOOBSF), ELAS )
C
C           LE COEFFICIENT DE DILATATION THERMIQUE HOMOGENE ISOTROPE
C           ( DILATA, DILATA, 0 )
            CALL REDILA( 3,NOOBSF, XD,YD,0D0,
     &                   LTDESU(LPDILA,NOOBSF),DILATA,TEMPIN )
C
            D = DILATA * (TEMPER(NUNDEL(NUELEM,L))-TEMPIN) / 6D0
            AUX(1) = AUX(1) + (ELAS(1)+ELAS(2)) * D
            AUX(2) = AUX(2) + (ELAS(2)+ELAS(3)) * D
            AUX(3) = AUX(3) + (ELAS(4)+ELAS(5)) * D
 120     CONTINUE
C
         BE(1) = BE(1) - AUX(1) * Y32 + AUX(3) * X32
         BE(2) = BE(2) + AUX(2) * X32 - AUX(3) * Y32
         BE(3) = BE(3) + AUX(1) * Y31 - AUX(3) * X31
         BE(4) = BE(4) - AUX(2) * X31 + AUX(3) * Y31
         BE(5) = BE(5) - AUX(1) * Y21 + AUX(3) * X21
         BE(6) = BE(6) + AUX(2) * X21 - AUX(3) * Y21
C
      ENDIF
C
      IF( PENALI .EQ. 0D0 ) RETURN
C
C     CONTRIBUTION DES ARETES A LA PENALISATION DES DEPLACEMENTS IMPOSES
C     ------------------------------------------------------------------
      DO 230 K = 1,3
C
C        LA LIGNE DU COTE K SUPPORTE T ELLE UNE FIXATION ?
         NL = NOOBLA(K)
         IF( NL .GT. 0 ) THEN
C
C           ADRESSE MCN DE LA FIXATION DE LA LIGNE
            MN = LTDELI(LPFIXA,NL)
            IF( MN .GT. 0 ) THEN
C
C              LE NUMERO DES 2 SOMMETS DU COTE K
               IF( K .NE. 3 ) THEN
                  KK = K+1
               ELSE
                  KK = 1
               ENDIF
C
               DO 220 M=1,2
C
C                 LE NUMERO DU SOMMET M DE L'ARETE K
                  MM = NOSOAR(M)
C
C                 CALCUL DE LA FIXATION AU SOMMET MM
                  XD = X(MM,1)
                  YD = X(MM,2)
                  CALL REFIXA( 2, NL, XD, YD, 0.D0, MN, NBCOFI, FORCE )
C                 NBCOFI LE NOMBRE DE COMPOSANTES FIXEES DU DEPLACEMENT
C
                  DO 215 J=1,NBCOFI
C                    LE NUMERO DE LA COMPOSANTE FIXEE
                     L = MCN( MN + WUCOFI - 1 + J )
C                    LA PENALISATION EST AJOUTEE POUR LE DEPLACEMENT L AU SOMMET
                     I  = MM * 2 - 2 + L
                     BE(I) = BE(I) + PENALI * FORCE(J)
 215              CONTINUE
C
 220           CONTINUE
C
            ENDIF
         ENDIF
 230  CONTINUE
C
C     CONTRIBUTION DES SOMMETS A LA PENALISATION DES DEPLACEMENTS IMPOSES
C     -------------------------------------------------------------------
      DO 250 K=1,3
C
C        NUMERO DE POINT DU SOMMET K
         NL = NOOBPS(K)
         IF( NL .GT. 0 ) THEN
C
C           LE SOMMET K EST UN POINT. EST IL SUPPORT D'UNE FIXATION?
            MN = LTDEPO(LPFIXA,NL)
            IF( MN .GT. 0 ) THEN
C
C              CALCUL DE LA FIXATION AU SOMMET K
               XD = X(K,1)
               YD = X(K,2)
               CALL REFIXA( 1, N, XD, YD, 0.D0, MN, NBCOFI, FORCE )
C              NBCOFI  LE NOMBRE DE COMPOSANTES FIXEES
C
               DO 240 J=1,NBCOFI
C
C                 LE NUMERO DE LA COMPOSANTE FIXEE
                  L = MCN( MN + WUCOFI - 1 + J )
C
C                 LA PENALISATION EST AJOUTEE
                  I  = K * 2 - 2 + L
                  BE(I) = BE(I) + PENALI * FORCE(J)
C
 240           CONTINUE
C
            ENDIF
         ENDIF
 250  CONTINUE
C
      RETURN
      END
