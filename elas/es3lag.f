      SUBROUTINE ES3LAG( X,      PENALI, NDSM,
     %                   NOOBPS, NUMIPO, NUMAPO, LTDEPO,
     %                   NOOBLA, NUMILI, NUMALI, LTDELI,
     %                   NOOBSF, NUMISU, NUMASU, LTDESU,
     %                   NOOBVC, NUMIVO, NUMAVO, LTDEVO,
     %                   NBPOLY, NPI,    POLY,
     %                   NBPOTR, NPITR,  POIDTR, POLYTR, DPOLTR,
     %                   NBPOQU, NPIQU,  POIDQU, POLYQU, DPOLQU,
     %                   ELAS,   FORCES, CONINI, G1,     G2,
     %                   NUTYEL, NUELEM, NBELEM, NUNDEL,
     %                   MNTEMP, NTDLT,  TEMPER,
     %                   F,      POIDEL, DP,
     %                   BE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DES SECONDS MEMBRES DES ELEMENTS FINIS LAGRANGE 3D
C -----    DE DEGRE 1 OU 2 ISOPARAMETRIQUES SAUF TETR 3P1D
C
C ENTREES:
C --------
C X      : 3 COORDONNEES DES NBPOLY POINTS DE L ELEMENT FINI
C PENALI : COEFFICIENT DE PENALISATION DE LA CONDITION DE FIXATION
C          CETTE VALEUR 1/EPSILON DOIT ETRE PLUS GRANDE QUE LES
C          COEFFICIENTS DU VECTEUR ELEMENTAIRE ET GLOBAL
C NDSM   : NOMBRE DE SECONDS MEMBRES A CALCULER
C
C NBNSOM : NOMBRE DE SOMMETS DE L'EF
C NOOBPS : NUMERO DE POINT DES SOMMETS
C NUMIPO : NUMERO MINIMAL DES OBJETS POINTS
C NUMAPO : NUMERO MAXIMAL DES OBJETS POINTS
C LTDEPO : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES D'ELASTICITE DES POINTS
C
C NOOBLA : NUMERO DES OBJETS LIGNES DES ARETES DE L'ELEMENT
C NUMILI : NUMERO MINIMAL DES OBJETS LIGNES
C NUMALI : NUMERO MAXIMAL DES OBJETS LIGNES
C LTDELI : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES D'ELASTICITE DES LIGNES
C
C NOOBSF : NUMERO DES SURFACES DES FACES DE L'ELEMENT FINI
C NUMISU : NUMERO MINIMAL DES SURFACES DE L'OBJET
C NUMASU : NUMERO MAXIMAL DES SURFACES DE L'OBJET
C LTDESU : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES D'ELASTICITE DES SURFACES
C
C NOOBVC : NUMERO DU VOLUME DE CET ELEMENT FINI
C NUMIVO : NUMERO MINIMAL DES VOLUMES DE L'OBJET
C NUMAVO : NUMERO MAXIMAL DES VOLUMES DE L'OBJET
C LTDEVO : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES D'ELASTICITE
C          DES VOLUMES DE L'OBJET
C
C NBPOLY : NOMBRE DE POLYNOMES
C NPI    : NOMBRE DE POINTS D INTEGRATION NUMERIQUE SUR L'ELEMENT FINI
C POLY   : POLY(I,L) = PI(XL,YL)
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
C ELAS   : TENSEUR DE L ELASTICITE SYMETRIQUE CALCULE DANS CE SP
C FORCES : EFFORTS DE VOLUME TABLEAU REMPLI DANS CE SOUS PROGRAMME
C CONINI : CONTRAINTES INITIALES PAR NO DE SURFACE CONINI(NDSM,6)
C          TABLEAU REMPLI DANS CE SP
C G1,G2  : TABLEAUX AUXILIAIRES
C
C NUTYEL : NO DU TYPE DE CET ELEMENT FINI DANS LA CLASSIFICATION MEFISTO
C NUELEM : NUMERO DE L'ELEMENT FINI
C NBELEM : NOMBRE D'ELEMENTS FINIS DE CE TYPE DANS LE MAILLAGE DE L'OBJET
C NUNDEL : NUMERO DES NOEUDS DES ELEMENTS FINIS DE CE TYPE DANS L'OBJET
C
C MNTEMP : >0 CALCUL DEMANDE DES CONTRAINTES THERMIQUES, =<0 SINON
C NTDLT  : NOMBRE DE NOEUDS THERMIQUES OU DL THERMIQUES DU MAILLAGE
C TEMPER : TEMPERATURES AUX NOEUDS DU MAILLAGE TEMPER(NTDL,NDSM)
C
C F      : XYZ DES NPI POINTS D'INTEGRATION
C POIDEL : POIDS * JACOBIEN AUX NPI POINTS D'INTEGRATION DE L'EF
C DP     : DP(3,NBPOLY,NPI) GRADIENT AUX POINTS D 'INTEGRATION DES
C          FONCTIONS DE BASE ISOPARAMETRIQUES DE L'EF COURANT
C
C SORTIES:
C --------
C BE     : BE(NDSM,3*NBPOLY) LES NDSM SECONDS MEMBRES ELEMENTAIRES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE UPMC PARIS      AVRIL 1999
C-----------------------------------------------------------------------
      include"./incl/donela.inc"
      include"./incl/ponoel.inc"
      include"./incl/a___fixation.inc"
      include"./incl/pp.inc"
      COMMON           MCN(MOTMCN)
C
      DOUBLE PRECISION POIDTR(NPITR),
     %                 POLYTR(NBPOTR,NPITR),
     %                 DPOLTR(2,NBPOTR,NPITR)
      DOUBLE PRECISION POIDQU(NPIQU),
     %                 POLYQU(NBPOQU,NPIQU),
     %                 DPOLQU(2,NBPOQU,NPIQU)
      DOUBLE PRECISION PENALI
      DOUBLE PRECISION CONINI(NDSM,6),
     %                 POLY(NBPOLY,NPI),
     %                 POIDEL(NPI),
     %                 DP(3,NBPOLY,NPI),F(NPI,3),
     %                 FORCES(NDSM,3),
     %                 TEMPER(NTDLT,NDSM),
     %                 ELAS(21),
     %                 G1(NBPOLY,3),
     %                 G2(NBPOLY,3,NBPOLY),
     %                 BE(NDSM,*)
      REAL             X(NBPOLY,3)
      INTEGER          NUNDEL(1:NBELEM,1:NBPOLY)
      INTEGER          LTDEPO(1:MXDOEL,NUMIPO:NUMAPO)
      INTEGER          LTDELI(1:MXDOEL,NUMILI:NUMALI)
      INTEGER          LTDESU(1:MXDOEL,NUMISU:NUMASU)
      INTEGER          LTDEVO(1:MXDOEL,NUMIVO:NUMAVO)
      INTEGER          NOOBSF(1:*), NOOBLA(1:*), NOOBPS(1:*)
C
      DOUBLE PRECISION TEMPIN, DILATA, D, A(9), XD,YD,ZD
      INTEGER          NONOFK(8)
C
C     INITIALISATION DE BE
      CALL AZEROD( NDSM*3*NBPOLY, BE )
      NBPOL2 = 2 * NBPOLY
C
C     CONTRIBUTION DES EFFORTS VOLUMIQUES
C     -----------------------------------
      MN = LTDEVO(LPFORC,NOOBVC)
      IF( MN .GT. 0 ) THEN
C
         DO 15 L=1,NPI
C
C           LA VALEUR DES EFFORTS VOLUMIQUES EN CE POINT
C           PAS DE NORMALE EN UN POINT DU VOLUME
            CALL REFORC( 4, NOOBVC, 3,
     %                   F(L,1), F(L,2), F(L,3), 0D0, 0D0, 0D0,
     %                   MN, FORCES )
C
            DO 10 I=1,NBPOLY
C
C              CONTRIBUTION DU POIDS ET DU POLYNOME DE BASE
               D = POIDEL(L) * POLY(I,L)
C
               DO 8 K=1,3
C
C                 NUMERO DU DL
                  II = I * 3 - 3 + K
C
                  DO 5 N=1,NDSM
                     BE(N,II) = BE(N,II) + D * FORCES(N,K)
    5             CONTINUE
C
    8          CONTINUE
   10       CONTINUE
   15    CONTINUE
      ENDIF
C
C     CONTRIBUTIONS DES FORCES SUR LES FACES
C     --------------------------------------
      DO 50 K=1,NFACE
C
C        NO DE SURFACE DE LA FACE K
         NOOB = NOOBSF(K)
         IF( NOOB .GT. 0 ) THEN
C
C           LA FACE EST SUR UNE SURFACE SUPPORT DE FORCE?
            MN = LTDESU( LPFORC, NOOB )
            IF( MN .GT. 0 ) THEN
C
C              UN TABLEAU FORCE EXISTE POUR CETTE FACE K
C              CALCUL DE LA CONTRIBUTION DE LA FACE K A BE
C              ...........................................
C              RECHERCHE DU NUMERO DES POINTS=NOEUDS DE LA FACE K
               CALL ELNOFA( NUTYEL, K, NBNOFK, NONOFK )
C              NONOFK(I) LES NUMEROS DES NOEUDS DE LA FACE K
C
               IF( NBSOFA(K) .EQ. 3 ) THEN
C
C                 FACE TRIANGULAIRE
                  CALL E43LAG( X,        NDSM,   NBPOLY, NBNOFK, NONOFK,
     &                        POLYTR,    DPOLTR, NPITR,  POIDTR,
     &                        NOOBSF(K), NUMISU, NUMASU, LTDESU,
     &                        FORCES,
     &                        BE )
               ELSE
C
C                 FACE QUADRANGULAIRE
                  CALL E43LAG( X,        NDSM,   NBPOLY, NBNOFK, NONOFK,
     &                        POLYQU,    DPOLQU, NPIQU,  POIDQU,
     &                        NOOBSF(K), NUMISU, NUMASU, LTDESU,
     &                        FORCES,
     &                        BE )
               ENDIF
            ENDIF
         ENDIF
 50   CONTINUE
C
C     CONTRIBUTION DES CONTRAINTES INITIALES
C     --------------------------------------
      IF( LTDEVO(LPCOIN,NOOBVC) .GT. 0 ) THEN
C
         DO 85 L=1,NPI
C           LA VALEUR DES CONTRAINTES INITIALES EN CE POINT
            CALL RECOIN( 4, NOOBVC, 6, F(L,1),F(L,2),F(L,3),
     %                   LTDEVO(LPCOIN,NOOBVC), CONINI )
C
            DO 80 N=1,NDSM
C              CONINI =  S11 S22 S33  S12 S23 S31  DEVIENT TD * SIGMA0
C                   A =  S11 S12 S31  S12 S22 S23  S31 S23 S33
               A(1) = CONINI(N,1)
               A(2) = CONINI(N,4)
               A(3) = CONINI(N,6)
C
               A(4) = CONINI(N,4)
               A(5) = CONINI(N,2)
               A(6) = CONINI(N,5)
C
               A(7) = CONINI(N,6)
               A(8) = CONINI(N,5)
               A(9) = CONINI(N,3)
C
C              CALCUL DE G1(NBPOLY,1) = T(DP) T(A(1),A(2),A(3))
               CALL TAB0D(NBPOLY,3,1,DP(1,1,L),A(1),G1(1,1))
C
C              CALCUL DE G1(NBPOLY,2) = T(DP) T(A(4),A(5),A(6))
               CALL TAB0D(NBPOLY,3,1,DP(1,1,L),A(4),G1(1,2))
C
C              CALCUL DE G1(NBPOLY,3) = T(DP) T(A(7),A(8),A(9))
               CALL TAB0D(NBPOLY,3,1,DP(1,1,L),A(7),G1(1,3))
C
C              AJOUT AUX SECONDS MEMBRES ELEMENTAIRES
               DO 70 I=1,NBPOLY
                  I1 = I*3 - 2
                  I2 = I1 + 1
                  I3 = I2 + 1
                  BE(N,I1) = BE(N,I1) - POIDEL(L) * G1(I,1)
                  BE(N,I2) = BE(N,I2) - POIDEL(L) * G1(I,2)
                  BE(N,I3) = BE(N,I3) - POIDEL(L) * G1(I,3)
 70            CONTINUE
 80         CONTINUE
 85      CONTINUE
      ENDIF
C
C     CONTRIBUTIONS DES CONTRAINTES THERMIQUES
C     ----------------------------------------
      IF( MNTEMP .GT. 0 .AND.
     &    LTDEVO(LPDILA,NOOBVC) .GT. 0 ) THEN
C
C         INITIALISATION DE G2
          CALL AZEROD( NBPOLY*3*NBPOLY, G2 )
C
          DO 150 L=1,NPI
C
C            LE TENSEUR DE L'ELASTICITE EN CE POINT D'INTEGRATION
             CALL REELAS( 4,NOOBVC,3,F(L,1),F(L,2),F(L,3),
     &                    LTDEVO(LPYOUN,NOOBVC), ELAS )
C
C            LE COEFFICIENT DE DILATATION THERMIQUE HOMOGENE ISOTROPE
             CALL REDILA( 4,NOOBVC,F(L,1),F(L,2),F(L,3),
     &                    LTDEVO(LPDILA,NOOBVC), DILATA, TEMPIN )
C
C            CALCUL DE A=T(D)*E*DILATA*POIDS*JACOBIEN
             D    = POIDEL(L) * DILATA
             A(1) = ( ELAS( 1) + ELAS( 2) + ELAS( 4) ) * D
             A(2) = ( ELAS( 7) + ELAS( 8) + ELAS( 9) ) * D
             A(3) = ( ELAS(16) + ELAS(17) + ELAS(18) ) * D
C
             A(4) = ( ELAS( 7) + ELAS( 8) + ELAS( 9) ) * D
             A(5) = ( ELAS( 2) + ELAS( 3) + ELAS( 5) ) * D
             A(6) = ( ELAS(11) + ELAS(12) + ELAS(13) ) * D
C
             A(7) = ( ELAS(16) + ELAS(17) + ELAS(18) ) * D
             A(8) = ( ELAS(11) + ELAS(12) + ELAS(13) ) * D
             A(9) = ( ELAS( 4) + ELAS( 5) + ELAS( 6) ) * D
C
C            CALCUL DE G1(NBPOLY,1) = T(DP) T(A(1),A(2),A(3))
             CALL TAB0D(NBPOLY,3,1,DP(1,1,L),A(1),G1(1,1))
C
C            CALCUL DE G1(NBPOLY,2) = T(DP) T(A(4),A(5),A(6))
             CALL TAB0D(NBPOLY,3,1,DP(1,1,L),A(4),G1(1,2))
C
C            CALCUL DE G1(NBPOLY,3) = T(DP) T(A(7),A(8),A(9))
             CALL TAB0D(NBPOLY,3,1,DP(1,1,L),A(7),G1(1,3))
C
C            CALCUL DE G2 = (G1) (P)
             CALL AB1D(NBPOLY*3,1,NBPOLY,G1,POLY(1,L),G2)
C
 150      CONTINUE
C
C         LE SECOND MEMBRE DU AUX CONTRAINTES THERMIQUES
          DO 190 I=1,NBPOLY
             DO 180 J=1,3
                I1 = I * 3 - 3 + J
                DO 170 N=1,NDSM
                   D = 0D0
                   DO 160 L=1,NBPOLY
                      D = D + G2(I,J,L) *
     &                  ( TEMPER(NUNDEL(NUELEM,L),N) - TEMPIN )
 160               CONTINUE
                   BE(N,I1) = BE(N,I1) + D
 170            CONTINUE
 180         CONTINUE
 190      CONTINUE
C
      ENDIF
C
C     ========================================================
C     CONTRIBUTION DES FACES A LA PENALISATION DE LA FIXATION
C     ========================================================
      IF( PENALI .EQ. 0D0 ) RETURN
C
      DO 210 K=1,NFACE
C
C        NO DE SURFACE DE LA FACE K
         NOOB = NOOBSF(K)
         IF( NOOB .GT. 0 ) THEN
C
C           LA FACE EST SUR UNE SURFACE SUPPORT DE FIXATION PENALISEE?
            MN = LTDESU( LPFIXA, NOOB )
            IF( MN .LE. 0 ) GOTO 210
C
C           RECHERCHE DU NUMERO DES POINTS=NOEUDS DE LA FACE K
            CALL ELNOFA( NUTYEL, K, NBNOFK, NONOFK )
C           NONOFK(I) LES NUMEROS DES NOEUDS DE LA FACE K
C
            DO 205 I=1,NBNOFK
C
C              LE NUMERO ELEMENTAIRE DU NOEUD I DE LA FACE K
               NI = NONOFK(I)
C
C              IL EXISTE UN CONTACT SUR CETTE LIGNE
               XD = X( NI, 1 )
               YD = X( NI, 2 )
               ZD = X( NI, 3 )
               CALL REFIXA( 3, NOOB, XD, YD, ZD, MN,
     %                      NBCOFI, FORCES )
C              NBCOFI LE NOMBRE DE COMPOSANTES FIXEES DU DEPLACEMENT
C
               DO 202 J=1,NBCOFI
C
C                 LE NUMERO DE LA COMPOSANTE DE LA FIXATION
                  M = MCN( MN + WUCOFI - 1 + J )
C
C                 LE NUMERO DU DL DU SOMMET NI
                  II = NI * 3 - 3 + M
C
C                 CALCUL DE LA FIXATION PENALISEE
                  DO 200 N=1,NDSM
                     BE(N,II) = BE(N,II) + PENALI * FORCES(N,J)
 200              CONTINUE
C
 202           CONTINUE
C
 205        CONTINUE
         ENDIF
 210  CONTINUE
C
C     ========================================================
C     CONTRIBUTION DES ARETES A LA PENALISATION DE LA FIXATION
C     ========================================================
      DO 250 K=1,NARET
C
C        NO DE LIGNE DE L'ARETE K
         NOOB = NOOBLA(K)
         IF( NOOB .GT. 0 ) THEN
C
C           LE COTE K EST SUR UNE LIGNE. EST IL SUPPORT D'UN CONTACT PENALISE?
            MN = LTDELI(LPFIXA,NOOB)
            IF( MN .GT. 0 ) THEN
C
C              OUI: UN TABLEAU FIXATION PENALISE EXISTE POUR CE COTE
C              RECHERCHE DES NUMEROS LOCAUX DES NOEUDS DE L'ARETE K
               NONOFK(1) = NOSOAR(1,K)
               NONOFK(2) = NOSOAR(2,K)
               IF( NBNOAR(K) .GT. 0 ) NONOFK(3) = NONOAR(1,K)
C
               DO 240 I=1,2+NBNOAR(K)
C
C                 LE NUMERO DU I-EME NOEUD DE L'ARETE K
                  NI = NONOFK(I)
C
C                 IL EXISTE UN CONTACT SUR CETTE LIGNE
                  XD = X( NI, 1 )
                  YD = X( NI, 2 )
                  ZD = X( NI, 3 )
                  CALL REFIXA( 2, NOOB, XD, YD, ZD, MN,
     %                         NBCOFI, FORCES )
C
                  DO 230 J=1,NBCOFI
C
C                    LE NUMERO DE LA COMPOSANTE DE LA FIXATION
                     M = MCN( MN + WUCOFI - 1 + J )
C
C                    LE NUMERO DU DL DU SOMMET K
                     II = NI * 3 - 3 + M
C
C                    CALCUL DE LA FIXATION PENALISEE
                     DO 220 N=1,NDSM
                        BE(N,II) = BE(N,II) + PENALI * FORCES(N,J)
 220                 CONTINUE
C
 230              CONTINUE
C
 240           CONTINUE
C
            ENDIF
         ENDIF
 250  CONTINUE
C
C     =========================================================
C     CONTRIBUTION DES SOMMETS A LA PENALISATION DE LA FIXATION
C     =========================================================
      DO 300 K=1,NBNSOM
C
C        NO DE POINT DU SOMMET K
         NOOB = NOOBPS(K)
         IF( NOOB .GT. 0 ) THEN
C
C           LE SOMMET K EST UN POINT. EST IL SUPPORT D'UN CONTACT PENALISE?
            MN = LTDEPO( LPFIXA, NOOB )
            IF( MN .GT. 0 ) THEN
C
C              OUI: UN TABLEAU FIXATION PENALISEE EXISTE POUR CE POINT
C              LES SOMMETS SONT NUMEROTES EN PREMIER
C              PUIS, VIENNENT LES EVENTUELS MILIEUX DES COTES
               XD = X( K, 1 )
               YD = X( K, 2 )
               ZD = X( K, 3 )
               CALL REFIXA( 1, NOOB, XD, YD, ZD, MN,
     %                      NBCOFI, FORCES )
C
               DO 270 J=1,NBCOFI
C
C                 LE NUMERO DE LA COMPOSANTE DE LA FIXATION
                  M = MCN( MN + WUCOFI - 1 + J )
C
C                 LE NUMERO DU DL M DU SOMMET K
                  II = K * 3 - 3 + M
C
C                 CALCUL DE LA FIXATION PENALISEE
                  DO 260 N=1,NDSM
                     BE(N,II) = BE(N,II) + PENALI * FORCES(N,J)
 260              CONTINUE
C
 270           CONTINUE
            ENDIF
         ENDIF
 300  CONTINUE
C
      RETURN
      END
