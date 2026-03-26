      SUBROUTINE ES3P1D( X,      PENALI,
     %                   NOOBPS, NUMIPO, NUMAPO, LTDEPO,
     %                   NOOBLA, NUMILI, NUMALI, LTDELI,
     %                   NOOBSF, NUMISU, NUMASU, LTDESU,
     %                   NOOBVC, NUMIVO, NUMAVO, LTDEVO,
     %                   ELAS,   FORCES,
     %                   NUELEM, NBELEM, NUNDEL,
     %                   MNTEMP, NTDLTE, TEMPER,
     %                   DELTA,  DP,
     %                   BE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DU SECOND MEMBRE D'UN TETRAEDRE 3P1D
C -----    INTEGRATION NUMERIQUE AUX SOMMETS DU TETRAEDRE ET DES FACES
C
C ENTREES:
C --------
C X      : COORDONNEES DES NBPOLY POINTS DE L ELEMENT FINI
C PENALI : COEFFICIENT DE PENALISATION DE LA CONDITION DE FIXATION
C          CETTE VALEUR 1/EPSILON DOIT ETRE PLUS GRANDE QUE LES
C          COEFFICIENTS DU VECTEUR ELEMENTAIRE ET GLOBAL
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
C NOOBSF : NUMERO DES SURFACES DES FACES DE L'EF
C NUMISU : NUMERO MINIMAL DES SURFACES UTILISEES
C NUMASU : NUMERO MAXIMAL DES SURFACES UTILISEES
C LTDESU : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES D'ELASTICITE
C          DES SURFACES
C
C NOOBVC : NUMERO DE L'OBJET VOLUME DE CET EF
C NUMIVO : NUMERO MINIMAL DES OBJETS VOLUMES
C NUMAVO : NUMERO MAXIMAL DES OBJETS VOLUMES
C LTDEVO : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES D'ELASTICITE
C          DES VOLUMES
C
C ELAS   : TENSEUR DE L ELASTICITE SYMETRIQUE CALCULE DANS CE SP
C FORCES : EFFORTS DE VOLUME TABLEAU REMPLI DANS CE SOUS PROGRAMME
C
C NUELEM : NUMERO DE L'ELEMENT FINI
C NBELEM : NOMBRE D'ELEMENTS FINIS
C NUNDEL : NUMERO DES NOEUDS DES ELEMENTS FINIS
C
C MNTEMP : >0 CALCUL DEMANDE DES CONTRAINTES THERMIQUES,=<0 SINON
C NTDLTE : NOMBRE DE NOEUDS THERMIQUES OU DL THERMIQUES DU MAILLAGE
C TEMPER : TEMPERATURES AUX NOEUDS DU MAILLAGE TEMPER(NTDLTE)
C
C DELTA  : JACOBIEN DE LA TRANSFORMATION EF REFERENCE -> EF
C DP     : GRADIENT (CONSTANT) DES POLYNOMES DE BASE SUR LE TETRAEDRE
C          COURANT
C
C SORTIE :
C --------
C BE     : BE(12) LE SECOND MEMBRE ELEMENTAIRE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: DENIS KYNDT WALTER PERRONNET ANALYSE NUMERIQUE UPMC  JUIN 1996
C MODIF : ALAIN PERRONNET   ANALYSE NUMERIQUE UPMC PARIS      AVRIL 1999
C23456---------------------------------------------------------------012
      include "./incl/donela.inc"
      include "./incl/a___fixation.inc"
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      INTEGER           LTDEPO(1:MXDOEL,NUMIPO:NUMAPO)
      INTEGER           LTDELI(1:MXDOEL,NUMILI:NUMALI)
      INTEGER           LTDESU(1:MXDOEL,NUMISU:NUMASU)
      INTEGER           LTDEVO(1:MXDOEL,NUMIVO:NUMAVO)
      INTEGER           NOOBPS(1:4),
     %                  NOOBLA(1:6),
     %                  NOOBSF(1:4)
      DOUBLE PRECISION  PENALI,
     %                  DELTA,
     %                  DP(3,4),
     %                  FORCES(3),
     %                  CONINI(6),
     %                  ELAS(21),
     %                  AUX(9),
     %                  TEMPER(NTDLTE),
     %                  BE(12)
C
      REAL              X(4,3)
      INTEGER           NUNDEL(NBELEM,4)
C
      INTEGER           NONOFK(3)
      DOUBLE PRECISION  DILATA, TEMPIN
      DOUBLE PRECISION  D, DGL(2,3), DGLN, DELTAK, XD, YD, ZD, VN(3)
      INTEGER           NS(2)
      EQUIVALENCE      (NS(1),N1),(NS(2),N2)
C     MISE A ZERO DE BE=LE VECTEUR ELEMENTAIRE
C     ------------------------------------------------
      CALL AZEROD( 12, BE )
C
C     CONTRIBUTION DES EFFORTS VOLUMIQUES
C     -----------------------------------
      IF( LTDEVO(LPFORC,NOOBVC) .GT. 0 ) THEN
         L = 0
         D = DELTA / 24D0
       DO 15 K=1,4
C           LA VALEUR DES EFFORTS VOLUMIQUES EN CE SOMMET K
            XD = X(K,1)
            YD = X(K,2)
            ZD = X(K,3)
C           PAS DE NORMALE EN POINT DU VOLUME
          CALL REFORC( 4, NOOBVC, 3, XD, YD, ZD, 0D0, 0D0, 0D0,
     %                   LTDEVO(LPFORC,NOOBVC), FORCES )
          DO 9 J=1,3
               L = L + 1
               BE(L) = BE(L) + D * FORCES(J)
   9        CONTINUE
   15    CONTINUE
      ENDIF
C
C     CONTRIBUTIONS DES EFFORTS SUR LES FACES
C     ---------------------------------------
      DO 50 K=1,4
         IF( NOOBSF(K) .GT. 0 ) THEN
C
C           LA FACE K EST SUR UNE SURFACE AVEC EFFORT
          IF( LTDESU(LPFORC,NOOBSF(K)) .GT. 0 ) THEN
C
C              UN TABLEAU FORCE EXISTE POUR CETTE FACE K
C              CALCUL DE LA CONTRIBUTION DE LA FACE K A BE
C              ...........................................
C              RECHERCHE DU NUMERO DES POINTS=NOEUDS DE LA FACE K
               CALL ELNOFA( 19, K, NBNOFK, NONOFK )
C              NONOFK(I) LES NUMEROS DES NOEUDS DE LA FACE K
C                        CE SONT AUSSI LES POINTS D'INTEGRATION
C
C              RECHERCHE DU JACOBIEN DE G
               N = NONOFK(1)
               DGL(1,1) = X( NONOFK(2), 1 ) - X( N, 1 )
               DGL(2,1) = X( NONOFK(3), 1 ) - X( N, 1 )
               DGL(1,2) = X( NONOFK(2), 2 ) - X( N, 2 )
               DGL(2,2) = X( NONOFK(3), 2 ) - X( N, 2 )
               DGL(1,3) = X( NONOFK(2), 3 ) - X( N, 3 )
               DGL(2,3) = X( NONOFK(3), 3 ) - X( N, 3 )
               CALL JAR2R3( DGL, DELTAK )
               D = DELTAK  / 6D0
C
C              LE VECTEUR NORMAL UNITAIRE A LA FACE
               CALL VECNOR( DGL, DGLN, VN )
C
               DO L=1,3
C
C                 CALCUL DES COORDONNEES DU SOMMET = POINT D INTEGRATION
                  II = NONOFK( L )
                  XD = X( II, 1 )
                  YD = X( II, 2 )
                  ZD = X( II, 3 )
C
C                 CALCUL DES FORCES EN CE SOMMET L DE LA FACE K
                  CALL REFORC( 3, NOOBSF(K), 3,
     %                         XD,YD,ZD, VN(1),VN(2),VN(3),
     %                         LTDESU(LPFORC,NOOBSF(K)), FORCES )
C
C                 LE NUMERO DES 3 DL DANS LA NUMEROTATION PAR SOMMETS
                  I1 = II * 3 - 2
                  I2 = I1 + 1
                  I3 = I2 + 1
C
C                 SOMMATION AVEC LE VECTEUR ELEMENTAIRE
                  BE(I1) = BE(I1) + D * FORCES(1)
                  BE(I2) = BE(I2) + D * FORCES(2)
                  BE(I3) = BE(I3) + D * FORCES(3)
C
               ENDDO
            ENDIF
         ENDIF
 50   CONTINUE
C
C     CONTRIBUTION DES CONTRAINTES INITIALES
C     --------------------------------------
      IF( LTDEVO(LPCOIN,NOOBVC) .GT. 0 ) THEN
         CALL AZEROD( 6, AUX )
         DO 70 L=1,4
C           LA VALEUR DES CONTRAINTES INITIALES EN CE SOMMET L
            XD = X(L,1)
            YD = X(L,2)
            ZD = X(L,3)
            CALL RECOIN( 4,NOOBVC,6, XD,YD,ZD,
     %                   LTDEVO(LPCOIN,NOOBVC), CONINI )
            DO 55 K=1,6
               AUX(K) = AUX(K) + CONINI(K)
 55         CONTINUE
 70      CONTINUE
C
         D = DELTA / 24D0
C
C        SOMMET 1
         BE( 1) = BE( 1)
     %   - ( DP(1,1)*AUX(1) + DP(2,1)*AUX(4) + DP(3,1)*AUX(6) )*D
         BE( 2) = BE( 2)
     %   - ( DP(1,1)*AUX(4) + DP(2,1)*AUX(2) + DP(3,1)*AUX(5) )*D
         BE( 3) = BE( 3)
     %   - ( DP(1,1)*AUX(6) + DP(2,1)*AUX(5) + DP(3,1)*AUX(3) )*D
C
C        SOMMET 2
         BE( 4) = BE( 4)
     %   - ( DP(1,2)*AUX(1) + DP(2,2)*AUX(4) + DP(3,2)*AUX(6) )*D
         BE( 5) = BE( 5)
     %   - ( DP(1,2)*AUX(4) + DP(2,2)*AUX(2) + DP(3,2)*AUX(5) )*D
         BE( 6) = BE( 6)
     %   - ( DP(1,2)*AUX(6) + DP(2,2)*AUX(5) + DP(3,2)*AUX(3) )*D
C
C        SOMMET 3
         BE( 7) = BE( 7)
     %   - ( DP(1,3)*AUX(1) + DP(2,3)*AUX(4) + DP(3,3)*AUX(6) )*D
         BE( 8) = BE( 8)
     %   - ( DP(1,3)*AUX(4) + DP(2,3)*AUX(2) + DP(3,3)*AUX(5) )*D
         BE( 9) = BE( 9)
     %   - ( DP(1,3)*AUX(6) + DP(2,3)*AUX(5) + DP(3,3)*AUX(3) )*D
C
C        SOMMET 4
         BE(10) = BE(10)
     %   - ( DP(1,4)*AUX(1) + DP(2,4)*AUX(4) + DP(3,4)*AUX(6) )*D
         BE(11) = BE(11)
     %   - ( DP(1,4)*AUX(4) + DP(2,4)*AUX(2) + DP(3,4)*AUX(5) )*D
         BE(12) = BE(12)
     %   - ( DP(1,4)*AUX(6) + DP(2,4)*AUX(5) + DP(3,4)*AUX(3) )*D
C
      ENDIF
C
C     CONTRIBUTIONS DES CONTRAINTES THERMIQUES
C     ----------------------------------------
      IF( MNTEMP .GT. 0 .AND.
     &    LTDEVO(LPDILA,NOOBVC) .GT. 0 ) THEN
C
         CALL AZEROD( 9, AUX )
C
         DO 120 L=1,4
C           LA VALEUR DES CONTRAINTES INITIALES EN CE SOMMET L
            XD = X(L,1)
            YD = X(L,2)
            ZD = X(L,3)
C
C           LE TENSEUR DE L'ELASTICITE EN CE POINT D'INTEGRATION
            CALL REELAS( 4,NOOBVC,3, XD,YD,ZD,
     &                   LTDEVO(LPYOUN,NOOBVC), ELAS )
C
C           LE COEFFICIENT DE DILATATION THERMIQUE HOMOGENE ISOTROPE
C           ( DILATA, DILATA, 0 )
            CALL REDILA( 4,NOOBVC, XD,YD,ZD,
     &                   LTDEVO(LPDILA,NOOBVC),DILATA,TEMPIN )
C
            D = DILATA*DELTA/24D0 * (TEMPER(NUNDEL(NUELEM,L))-TEMPIN)
C
            AUX(1) = AUX(1) + (ELAS( 1)+ELAS( 2)+ELAS( 4)) * D
            AUX(2) = AUX(2) + (ELAS( 7)+ELAS( 8)+ELAS( 9)) * D
            AUX(3) = AUX(3) + (ELAS(16)+ELAS(17)+ELAS(18)) * D
C
            AUX(4) = AUX(4) + (ELAS( 7)+ELAS( 8)+ELAS( 9)) * D
            AUX(5) = AUX(5) + (ELAS( 2)+ELAS( 3)+ELAS( 5)) * D
            AUX(6) = AUX(6) + (ELAS(11)+ELAS(12)+ELAS(13)) * D
C
            AUX(7) = AUX(7) + (ELAS(16)+ELAS(17)+ELAS(18)) * D
            AUX(8) = AUX(8) + (ELAS(11)+ELAS(12)+ELAS(13)) * D
            AUX(9) = AUX(9) + (ELAS( 4)+ELAS( 5)+ELAS( 6)) * D
 120     CONTINUE
C
C        SOMMET 1
         BE( 1) = BE( 1)
     %          + DP(1,1)*AUX(1) + DP(2,1)*AUX(2) + DP(3,1)*AUX(3)
         BE( 2) = BE( 2)
     %          + DP(1,1)*AUX(4) + DP(2,1)*AUX(5) + DP(3,1)*AUX(6)
         BE( 3) = BE( 3)
     %          + DP(1,1)*AUX(7) + DP(2,1)*AUX(8) + DP(3,1)*AUX(9)
C
C        SOMMET 2
         BE( 4) = BE( 4)
     %          + DP(1,2)*AUX(1) + DP(2,2)*AUX(2) + DP(3,2)*AUX(3)
         BE( 5) = BE( 5)
     %          + DP(1,2)*AUX(4) + DP(2,2)*AUX(5) + DP(3,2)*AUX(6)
         BE( 6) = BE( 6)
     %          + DP(1,2)*AUX(7) + DP(2,2)*AUX(8) + DP(3,2)*AUX(9)
C
C        SOMMET 3
         BE( 7) = BE( 7)
     %          + DP(1,3)*AUX(1) + DP(2,3)*AUX(2) + DP(3,3)*AUX(3)
         BE( 8) = BE( 8)
     %          + DP(1,3)*AUX(4) + DP(2,3)*AUX(5) + DP(3,3)*AUX(6)
         BE( 9) = BE( 9)
     %          + DP(1,3)*AUX(7) + DP(2,3)*AUX(8) + DP(3,3)*AUX(9)
C
C        SOMMET 4
         BE(10) = BE(10)
     %          + DP(1,4)*AUX(1) + DP(2,4)*AUX(2) + DP(3,4)*AUX(3)
         BE(11) = BE(11)
     %          + DP(1,4)*AUX(4) + DP(2,4)*AUX(5) + DP(3,4)*AUX(6)
         BE(12) = BE(12)
     %          + DP(1,4)*AUX(7) + DP(2,4)*AUX(8) + DP(3,4)*AUX(9)
C
      ENDIF
C
C     ============================================================
C     CONTRIBUTION DES FACES A LA CONDITION DE DIRICHLET PENALISEE
C     ============================================================
      IF( PENALI .EQ. 0D0 ) RETURN
C
      DO 140 K=1,4
C
C        LE NUMERO DE SURFACE DE LA FACE K
         NOOB = NOOBSF(K)
         IF( NOOB .GT. 0 ) THEN
C
C           LA FACE K EST SUR UNE SURFACE AVEC FIXATION
            MN = LTDESU(LPFIXA,NOOBSF(K))
          IF( MN .GT. 0 ) THEN
C
C              UN TABLEAU FORCE EXISTE POUR CETTE FACE K
C              CALCUL DE LA CONTRIBUTION DE LA FACE K A BE
C              ...........................................
C              RECHERCHE DU NUMERO DES POINTS=NOEUDS DE LA FACE K
               CALL ELNOFA( 19, K, NBNOFK, NONOFK )
C              NONOFK(I) LES NUMEROS DES NOEUDS DE LA FACE K
C
               DO 138 L=1,3
C
C                 CALCUL DES COORDONNEES DU SOMMET L
                  II = NONOFK( L )
                  XD = X( II, 1 )
                  YD = X( II, 2 )
                  ZD = X( II, 3 )
C
C                 CALCUL DES FORCES EN CE SOMMET L DE LA FACE K
                  CALL REFIXA( 3,  NOOB, XD, YD, ZD, MN,
     %                         NBCOFI, FORCES )
C                 NBCOFI LE NOMBRE DE COMPOSANTES FIXEES DU DEPLACEMENT
C
                  DO 134 J=1,NBCOFI
C
C                    LE NUMERO DE LA COMPOSANTE FIXEE
                     M = MCN( MN + WUCOFI - 1 + J )
C
C                    LE NUMERO DU DL
                     I1 = II * 3 - 3 + M
C
C                    SOMMATION AVEC LE VECTEUR ELEMENTAIRE
                     BE(I1) = BE(I1) + PENALI * FORCES(J)
 134             CONTINUE
C
 138           CONTINUE
            ENDIF
         ENDIF
 140  CONTINUE
C
C     =============================================================
C     CONTRIBUTION DES ARETES A LA CONDITION DE DIRICHLET PENALISEE
C     =============================================================
      DO 170 K=1,6
C
C           NO DE LIGNE DE L'ARETE K
            NOOB = NOOBLA(K)
            IF( NOOB .GT. 0 ) THEN
C
C              LE SOMMET K EST UN POINT. EST IL SUPPORT D'UNE FIXATION PENALISEE
               MN = LTDELI( LPFIXA, NOOB )
               IF( MN .GT. 0 ) THEN
C
C                 OUI: UN TABLEAU FIXATION EXISTE POUR CETTE LIGNE
C                 FIXATION PENALISEE = PENALI
C
C                 LE NUMERO DES 2 SOMMETS DE L'ARETE K
                  GOTO( 151, 152, 153, 154, 155, 156 ) , K
 151              N1 = 1
                  N2 = 2
                  GOTO 158
 152              N1 = 2
                  N2 = 3
                  GOTO 158
 153              N1 = 3
                  N2 = 1
                  GOTO 158
 154              N1 = 1
                  N2 = 4
                  GOTO 158
 155              N1 = 2
                  N2 = 4
                  GOTO 158
 156              N1 = 3
                  N2 = 4
C
 158              DO 168 L=1,2
C
C                    LE NUMERO DU SOMMET L DE L'ARETE K
                     N = NS(L)
C
C                    CALCUL DE LA FIXATION AU SOMMET N
                     XD = X(N,1)
                     YD = X(N,2)
                     ZD = X(N,3)
                     CALL REFIXA( 2, NOOB, XD, YD, ZD,  MN,
     %                            NBCOFI, FORCES )
C                    NBCOFI LE NOMBRE DE COMPOSANTES FIXEES DU DEPLACEMENT
C
                     DO 165 J=1,NBCOFI
C
C                       LE NUMERO DE LA COMPOSANTE FIXEE
                        M = MCN( MN + WUCOFI - 1 + J )
C
C                       LE NUMERO DU DL DU SOMMET L
                        N = N * 3 - 3 + M
C
C                       SOMMATION AVEC LE  VECTEUR ELEMENTAIRE
                        BE(N) = BE(N) + PENALI * FORCES(J)
C
 165                 CONTINUE
C
 168              CONTINUE
               ENDIF
            ENDIF
 170  CONTINUE
C
C     ==============================================================
C     CONTRIBUTION DES SOMMETS A LA CONDITION DE DIRICHLET PENALISEE
C     ==============================================================
      DO 200 K=1,4
C
C           NO DE POINT DU SOMMET K
            NOOB = NOOBPS(K)
            IF( NOOB .GT. 0 ) THEN
C
C              LE SOMMET K EST UN POINT. EST IL SUPPORT D'UNE FIXATION PENALISEE
               MN = LTDEPO(LPFIXA,NOOB)
               IF( MN .GT. 0 ) THEN
C
C                 OUI: UN TABLEAU FIXATION EXISTE POUR CE POINT
C                 FIXATION PENALISEE = PENALI
C                 CALCUL DE LA FIXATION AU SOMMET K
                  XD = X(K,1)
                  YD = X(K,2)
                  ZD = X(K,3)
                  CALL REFIXA( 1, NOOB, XD, YD, ZD, MN,
     %                         NBCOFI, FORCES )
C
                  DO 190 J=1,NBCOFI
C
C                    LE NUMERO DE LA COMPOSANTE FIXEE
                     M = MCN( MN + WUCOFI - 1 + J )
C
C                    LE NUMERO DU DL DU SOMMET K
                     N = K * 3 - 3 + M
C
C                    SOMMATION AVEC LE  VECTEUR ELEMENTAIRE
                     BE(N) = BE(N) + PENALI * FORCES(J)
C
 190              CONTINUE
C
               ENDIF
            ENDIF
 200  CONTINUE
C
      RETURN
      END
