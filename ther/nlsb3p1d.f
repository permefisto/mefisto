      SUBROUTINE NLSB3P1D( NUELEM, NONOEF, Omega,
     %                     DeltaT, X, PENALI, NBJEUX, JEU,
     %                     DELTA,  DP,
     %                     NOOBPS, NUMIPO, NUMAPO, LTDEPO,
     %                     NOOBLA, NUMILI, NUMALI, LTDELI,
     %                     NOOBSF, NUMISU, NUMASU, LTDESU,
     %                     NOOBVO, NUMIVO, NUMAVO, LTDEVO,
     %                     NBNOMA, Utn,    Utm,
     %                     BE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : NLSE: TESTNL=6 SCHEMA SEMI-IMPLICITE CALCUL DU SECOND MEMBRE BE
C ----- SUR UN EF TETRAEDRE P1 LAGRANGE 3D
C     REAL PART:
C      Fr(tn+1) +Fi(tn+1) + Rho/dt (Vn +Wn+1m-Wn)
C     +OmegaZ ( x d/dy - y d/dx )(Vn+1m-Wn+1m)
C     -N(Vn+1m**2+Wn+1m**2)(Vn+1m+Wn+1m) + N(V0**2+W0**2) Vn+1m
C     +Alfa LAPLACIEN Wn+1m
C
C     IMAG PART:
C     -Fr(tn+1) +Fi(tn+1) + Rho/dt (-Vn+1m+Vn +Wn)
C     +OmegaZ ( x d/dy - y d/dx )(Vn+1m+Wn+1m)
C     +N(Vn+1m**2+Wn+1m**2)(Vn+1m-Wn+1m) + N(V0**2+W0**2) Wn+1m
C     -Alfa LAPLACIEN Vn+1m

C ENTREES:
C --------
C NUELEM : NUMERO DE L'EF TRAITE
C NONOEF : NUMERO DES 4 NOEUDS DE L'EF NUELEM
C DELTA  : DETERMINANT DE LA MATRICE DF JACOBIENNE
C DP     : GRADIENT (CONSTANT) DES POLYNOMES DE BASE SUR LE TETRAEDRE
C          COURANT
C Omega  : VITESSE ANGULAIRE DE LA ROTATION
C DeltaT : PAS DE TEMPS DU SCHEMA IMPLICITE EN TEMPS
C X      : 3 COORDONNEES DES 4 POINTS DE L'EF
C PENALI : COEFFICIENT DE PENALISATION DE LA CONDITION DE DIRICHLET
C NBJEUX : NOMBRE DE JEUX DE DONNEES
C JEU    : NUMERO DU JEU  DE DONNEES POUR CE CALCUL DE LA MATRICE ELEMENTAIRE
C NBNOMA : NOMBRE DE NOEUDS DU MAILLAGE SUPPORT DE DL
C
C NOOBPS : NUMERO DE POINT DES SOMMETS
C NUMIPO : NUMERO MINIMAL DES OBJETS POINTS
C NUMAPO : NUMERO MAXIMAL DES OBJETS POINTS
C LTDEPO : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES THERMIQUES DES POINTS
C
C NOOBLA : NOUMERO DES OBJETS LIGNES DES ARETES DE L'ELEMENT
C NUMILI : NUMERO MINIMAL DES OBJETS LIGNES
C NUMALI : NUMERO MAXIMAL DES OBJETS LIGNES
C LTDELI : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES DE THERMIQUE
C          DES OBJETS LIGNES
C
C NOOBSF : NUMERO DE L'OBJET SURFACE DE CET ELEMENT
C NUMISU : NUMERO MINIMAL DES OBJETS SURFACES
C NUMASU : NUMERO MAXIMAL DES OBJETS SURFACES
C LTDESU : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES DE THERMIQUE
C          DES OBJETS SURFACES

C NOOBVO : NUMERO DU VOLUME DE L'EF
C NUMIVO : NUMERO MINIMAL DES OBJETS VOLUMES UTILISES
C NUMAVO : NUMERO MAXIMAL DES OBJETS VOLUMES UTILISES
C LTDEVO : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES CAPACITE
C          DES OBJETS VOLUMES

C NBNOMA : NOMBRE DE NOEUDS DU MAILLAGE
C Utn    : U(tn) (NBNOMA,2) U A L'INTANT tn
C Utm    : U(tn+1,m) (NBNOMA,2)  U A L'INTANT tn+1 iteration m

C
C SORTIE :
C --------
C BE     : BE(4,2) LE SECOND MEMBRE ELEMENTAIRE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St Pierre du Perray     Mars 2014
C23456---------------------------------------------------------------012
      include"./incl/donthe.inc"
      include"./incl/ponoel.inc"
      include"./incl/cthet.inc"
      include"./incl/ctemps.inc"
      include"./incl/cnonlin.inc"
      include"./incl/a___fixation.inc"
C
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE      (MCN(1), RMCN(1), DMCN(1))

      DOUBLE PRECISION  BE(4,2), PENALI, Omega(3)
      DOUBLE PRECISION  DP(3,4)
      REAL              X(4, 3)
      INTEGER           LTDEPO( 1:MXDOTH, 1:NBJEUX, NUMIPO:NUMAPO )
      INTEGER           LTDELI( 1:MXDOTH, 1:NBJEUX, NUMILI:NUMALI )
      INTEGER           LTDESU( 1:MXDOTH, 1:NBJEUX, NUMISU:NUMASU )
      INTEGER           LTDEVO( 1:MXDOTH, 1:NBJEUX, NUMIVO:NUMAVO )
      INTEGER           NOOBPS( 1:4 ),
     %                  NOOBLA( 1:6 ),
     %                  NOOBSF( 1:4 )
      INTEGER           NONOFK(8)
      DOUBLE PRECISION  Utn(NBNOMA,2), Utm(NBNOMA,2)
C
      DOUBLE PRECISION  S, SV, SW, D, DELTA, DeltaT, XYZ(3), FIXA(2)
      DOUBLE PRECISION  Rho, FORCE(2), DELTAK,
     %                  F1dYF2dXV, F1dYF2dXW,
     %                  VEFtm(4), WEFtm(4),
     %                  VEFtn(4), WEFtn(4),
     %                  DGL(2,3), DGLN, NLSECOEF, NLSECOEF0,
     %                  C(6), COND(6), CONDUC(10), TABAUX(12)
      INTEGER           K12(2)
      EQUIVALENCE      (K12(1),K1), (K12(2),K2)
      INTRINSIC         SQRT

C     INITIALISATION DE LA VITESSE ANGULAIRE Omega/Z SUPPOSEE CONSTANTE
C     -----------------------------------------------------------------
      IF( NUELEM .EQ. 1 ) THEN
         MN = LTDEVO(LPVIANT,JEU,NOOBVO)
         IF( MN .GT. 0 ) THEN
C           RECUPERATION DE LA VITESSE ANGULAIRE Omega SUPPOSEE CONSTANTE
C           VECTEUR(3) DE VITESSE ANGULAIRE AU 1-ER POINT D'INTEGRATION
            XYZ(1) = X(1,1)
            XYZ(2) = X(1,2)
            XYZ(3) = X(1,3)
            CALL REVIAN( 4, NOOBVO,  XYZ(1), XYZ(2),  XYZ(3),
     %                   LTDEVO(LPVIANT,JEU,NOOBVO), Omega )
          ELSE
C            PAS DE ROTATION
             Omega(1) = 0D0
             Omega(2) = 0D0
             Omega(3) = 0D0
          ENDIF
      ENDIF
C
C     RECUPERATION DE L'ONDE COMPLEXE AU TEMPS ACTUEL tn+1m et tn
C     AUX 4 DL DE L ELEMENT FINI
C     A PARTIR DU VECTEUR GLOBAL DES DL AUX NOEUDS DE L'ONDE U = V + i W
C     ------------------------------------------------------------------
      DO I=1,4
C        NO GLOBAL DU NOEUD I DE L'EF
         J = NONOEF(I)
C        AU TEMPS tn
         VEFtn(I) = Utn(J,1)
         WEFtn(I) = Utn(J,2)
C        AU TEMPS tn+1m
         VEFtm(I) = Utm(J,1)
         WEFtm(I) = Utm(J,2)
      ENDDO

C     COEFFICIENT ALFA DE - ALFA LAPLACIEN C'EST A DIRE CALCUL DU
C     TENSEUR SYMETRIQUE DE CONDUCTIVITE AU POINT D'INTEGRATION L
C     si Gross-Pitaevski: COND = 1/2 [Identite]
C     -----------------------------------------------------------
      S = DELTA / 24D0

      CALL AZEROD( 10, CONDUC )
      MN = LTDEVO(LPCOND,JEU,NOOBVO)
      IF( MN .GT. 0 ) THEN

         CALL AZEROD( 6, C )

         DO I=1,4
C           RECHERCHE DU TENSEUR SYMETRIQUE DE CONDUCTIVITE AU POINT
C           D'INTEGRATION I
            XYZ(1) = X(I,1)
            XYZ(2) = X(I,2)
            XYZ(3) = X(I,3)
            TEMPEL = VEFtm(I)
            ONDEPI = WEFtm(I)
            CALL RECOND( 4, NOOBVO, 3, XYZ, MN, COND )
            DO J=1,6
               C(J) = C(J) + COND(J)
            ENDDO
         ENDDO

C        C = CONDUCTIVITE * DELTA * POIDS
         DO J=1,6
            C(J) = C(J) * S
         ENDDO

C        CONDUC = T(DP) * C * (DP)
         CALL TABA8D( 4, 3, DP, C, CONDUC, TABAUX )

      ENDIF
C
C     ==================================
C     CONTRIBUTION DES FORCES VOLUMIQUES
C     ==================================
      DO I = 1, 4

C        COORDONNEES DU POINT D'INTEGRATION I = SOMMET I
         XYZ(1) = X(I,1)
         XYZ(2) = X(I,2)
         XYZ(3) = X(I,3)

C        DENSITE DE MASSE /DeltaTemps
C        ----------------
         MN = LTDEVO(LPMAST,JEU,NOOBVO)
         IF( MN .GT. 0 ) THEN
            TEMPEL = VEFtm(I)
            ONDEPI = WEFtm(I)
            CALL REMASS( 4, NOOBVO, 3, XYZ, MN, Rho )
         ELSE
            Rho = 1D0
         ENDIF
         Rho = Rho / DeltaT
C
C        NLSECOEF COEFFICIENT DU TERME NON LINEAIRE AU TEMPS DEMANDE
C        -----------------------------------------------------------
         MN = LTDEVO(LPCOET,JEU,NOOBVO)
         IF( MN .GT. 0 ) THEN
            TEMPEL = VEFtm(I)
            ONDEPI = WEFtm(I)
            CALL RENLSE( 4, NOOBVO, 3,XYZ, TEMPS,   TEMPEL, ONDEPI,  MN,
     %                   NLSECOEF )
            TEMPEL0 = VEFtn(I)
            ONDEPI0 = WEFtn(I)
            CALL RENLSE( 4, NOOBVO, 3,XYZ, TEMPSINI,TEMPEL0,ONDEPI0, MN,
     %                   NLSECOEF0 )
         ELSE
            NLSECOEF0 = 0D0
            NLSECOEF  = 0D0
         ENDIF
C
C        FORCE R et I AU SECOND MEMBRE
C        -----------------------------
         MN = LTDEVO(LPSOUR,JEU,NOOBVO)
         IF( MN .GT. 0 ) THEN
            TEMPEL  = VEFtm(I)
            ONDEPI  = WEFtm(I)
            CALL REFORC( 4, NOOBVO, 2, XYZ(1),XYZ(2),XYZ(3),
     %                                 TEMPEL, ONDEPI, 0D0,
     %                   MN, FORCE )
         ELSE
            FORCE(1) = 0D0
            FORCE(2) = 0D0
         ENDIF

C        CALCUL de ( x d/dy - y d/dx ) Vn+1m(bl)
C        CALCUL de ( x d/dy - y d/dx ) Wn+1m(bl)
C        ---------------------------------------
         F1dYF2dXV = 0D0
         F1dYF2dXW = 0D0
         DO J=1,4
            D = X(I,1) * DP(2,J) - X(I,2) * DP(1,J)
            F1dYF2dXV = F1dYF2dXV + D * VEFtm(J)
            F1dYF2dXW = F1dYF2dXW + D * WEFtm(J)
         ENDDO
         F1dYF2dXV = F1dYF2dXV * Omega(1)
         F1dYF2dXW = F1dYF2dXW * Omega(1)

C        LE TERME - Alfa LAPLACIEN V ou W C'EST A DIRE
C        CALCUL DU TENSEUR SYMETRIQUE DE CONDUCTIVITE AU POINT
C        D'INTEGRATION L MULTIPLIE PAR VEF ou WEF
C        -----------------------------------------------------------
         GOTO( 11, 12, 13, 14 ), I
 11      SV = CONDUC(1)*VEFtm(1) +CONDUC(2)*VEFtm(2) +CONDUC(4)*VEFtm(3)
     %       +CONDUC(7)*VEFtm(4)
         SW = CONDUC(1)*WEFtm(1) +CONDUC(2)*WEFtm(2) +CONDUC(4)*WEFtm(3)
     %       +CONDUC(7)*WEFtm(4)
         GOTO 20
 12      SV = CONDUC(2)*VEFtm(1) +CONDUC(3)*VEFtm(2) +CONDUC(5)*VEFtm(3)
     %       +CONDUC(8)*VEFtm(4)
         SW = CONDUC(2)*WEFtm(1) +CONDUC(3)*WEFtm(2) +CONDUC(5)*WEFtm(3)
     %       +CONDUC(8)*WEFtm(4)
         GOTO 20
 13      SV = CONDUC(4)*VEFtm(1) +CONDUC(5)*VEFtm(2) +CONDUC(6)*VEFtm(3)
     %       +CONDUC(9)*VEFtm(4)
         SW = CONDUC(4)*WEFtm(1) +CONDUC(5)*WEFtm(2) +CONDUC(6)*WEFtm(3)
     %       +CONDUC(9)*WEFtm(4)
         GOTO 20
 14      SV = CONDUC(7)*VEFtm(1) +CONDUC(8)*VEFtm(2) +CONDUC(9)*VEFtm(3)
     %       +CONDUC(10)*VEFtm(4)
         SW = CONDUC(7)*WEFtm(1) +CONDUC(8)*WEFtm(2) +CONDUC(9)*WEFtm(3)
     %       +CONDUC(10)*WEFtm(4)
C
C        CALCUL DU VECTEUR ELEMENTAIRE BE
C        --------------------------------
C        REAL PART:
C          +Alfa LAPLACIEN Wn+1m
C           Fr(tn+1) +Fi(tn+1) + Rho/dt (Vn +Wn+1m-Wn)
C          +OmegaZ ( x d/dy - y d/dx )(Vn+1m -Wn+1m)
C          -N(Vn+1m**2+Wn+1m**2)(Vn+1m+Wn+1m) + N(V0**2+W0**2) Vn+1m
 20      BE(I,1) = -SW + S 
     %              * ( FORCE(1) + FORCE(2)
     %                + Rho * ( VEFtn(I) + WEFtm(I) - WEFtn(I) )
     %                + F1dYF2dXV - F1dYF2dXW
     %                - NLSECOEF  * ( VEFtm(I) + WEFtm(I) )
     %                + NLSECOEF0 * VEFtm(I)  )

C        IMAG PART:
C          -Alfa LAPLACIEN Vn+1m
C          -Fr(tn+1) +Fi(tn+1) + Rho/dt (-Vn+1m+Vn +Wn)
C          +OmegaZ ( x d/dy - y d/dx )(Vn+1m +Wn+1m)
C          +N(Vn+1m**2+Wn+1m**2)(Vn+1m-Wn+1m) + N(V0**2+W0**2) Wn+1m
         BE(I,2) = SV + S
     %              * ( FORCE(2) - FORCE(1)
     %                + Rho * ( -VEFtm(I) + VEFtn(I) + WEFtn(I) )
     %                + F1dYF2dXV + F1dYF2dXW
     %                + NLSECOEF  * ( VEFtm(I) - WEFtm(I) )
     %                + NLSECOEF0 * WEFtm(I)  )

      ENDDO
C
C     ===========================================================
C     CONTRIBUTION DES FORCES OU FIXATION PENALISEE SUR LES FACES
C     ===========================================================
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
C           UN TABLEAU CONTACT EXISTE POUR CETTE FACE K
C           CALCUL DE LA CONTRIBUTION DE LA FACE K A BE
C           ..................................................
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
C              LA TEMPERATURE AU SOMMET N
               TEMPEL = VEFtm( N )
               ONDEPI = WEFtm( N )
C
C              LES 3 COORDONNEES DU SOMMET L DE LA FACE K
               XYZ(1) = X(N,1)
               XYZ(2) = X(N,2)
               XYZ(3) = X(N,3)
C
               IF( IECHAN .EQ. 1 ) THEN

C                 FORCE(2) REQUISE POUR CONDITION NEUMANN OU FOURIER
                  MN = LTDESU(LPSOUR,JEU,NOOB)
C                 ATTENTION: ICI OndeR et OndeI SONT PASSES
C                            A LA PLACE DE VN le VECTEUR NORMAL!
                  CALL REFORC( 3, NOOB, 2, XYZ(1),XYZ(2),XYZ(3),
ccc                                        VN(1), VN(2), VN(3),
     %                                     TEMPEL, ONDEPI, 0D0,
     %                         MN, FORCE )
C
               ELSE

C                 FIXATION(2) PENALISEE POUR CONDITION NEUMANN OU FOURIER
                  MN = LTDESU(LPCONT,JEU,NOOB)
                  CALL REFIXA( 3, NOOB, XYZ(1),XYZ(2),XYZ(3), MN,
     %                         NBCOFI, FIXA )
                  FORCE(1) = 0D0
                  FORCE(2) = 0D0
C                 NBCOFI LE NOMBRE DE COMPOSANTES FIXEES
                  DO I = 1, NBCOFI
C                    LE NUMERO DE LA I-EME COMPOSANTE FIXEE
                     NU = MCN( MN + WUCOFI - 1 + I )
                     FORCE(NU) = FIXA(I) * PENALI
                  ENDDO
C
               ENDIF
C
               D = DELTAK / 6D0
C              SOMMATION AVEC LE VECTEUR ELEMENTAIRE
               BE(N,1) = BE(N,1) + D * FORCE(1)
               BE(N,2) = BE(N,2) + D * FORCE(2)

 60         CONTINUE
         ENDIF
 80   CONTINUE
C
C     ======================================================
C     CONTRIBUTION DE FORCE OU FIXATION PENALISEE AUX ARETES
C     ======================================================
      IF( PENALI .NE. 0D0 ) THEN
         DO 130 K=1,6
C
C           NO DE LIGNE DE L'ARETE K
            NOOB = NOOBLA(K)
            IF( NOOB .GT. 0 ) THEN
C
C              LE SOMMET K EST UN POINT. EST IL SUPPORT
C              D'UNE SOURCE OU D'UN CONTACT PENALISE?
               IF(LTDELI(LPSOUR,JEU,NOOB).GT.0 .OR.
     %            LTDELI(LPCONT,JEU,NOOB).GT.0 .AND. PENALI.NE.0D0) THEN
C
C                 OUI: UN TABLEAU CONTACT PENALISE EXISTE POUR CE POINT
C                 CONTACT PENALISE = PENALI x TEMPERATURE
C
C                 LE NUMERO DES 2 SOMMETS DE L'ARETE K
                  GOTO( 111, 112, 113, 114, 115, 116 ) , K
 111              K1 = 1
                  K2 = 2
                  GOTO 118
 112              K1 = 2
                  K2 = 3
                  GOTO 118
 113              K1 = 3
                  K2 = 1
                  GOTO 118
 114              K1 = 1
                  K2 = 4
                  GOTO 118
 115              K1 = 2
                  K2 = 4
                  GOTO 118
 116              K1 = 3
                  K2 = 4
C
C                 LA LONGUEUR DE L'ARETE K / 2
 118              DELTAK = SQRT( (X(K2,1)-X(K1,1)) ** 2
     %                         + (X(K2,2)-X(K1,2)) ** 2 ) * 0.5D0
C
                  DO M=1,2
C
C                    NUMERO ET COORDONNEES DU SOMMET M DE L'ARETE K
                     KM = K12(M)
                     XYZ(1) = X(KM,1)
                     XYZ(2) = X(KM,2)
                     XYZ(3) = X(KM,3)
C
                     MN = LTDELI(LPSOUR,JEU,NOOB)
                     IF( MN .GT. 0 ) THEN
C
C                       FORCE(2) REQUISE POUR CONDITION NEUMANN OU FOURIER
                        MN = LTDELI(LPSOUR,JEU,NOOB)
C                       ATTENTION: ICI OndeR et OndeI SONT PASSES
C                                  A LA PLACE DE VN le VECTEUR NORMAL!
                        TEMPEL = VEFtm( KM )
                        ONDEPI = WEFtm( KM )
                        CALL REFORC( 2, NOOB, 2, XYZ(1),XYZ(2),XYZ(3),
     %                                           TEMPEL, ONDEPI, 0D0,
     %                               MN, FORCE(1) )
                     ENDIF
C
                     MN = LTDELI(LPCONT,JEU,NOOB)
                     IF( PENALI .GT. 0D0 .AND. MN .GT. 0 ) THEN
C                       FIXATION(2) PENALISEE
                        CALL REFIXA( 2, NOOB, XYZ(1),XYZ(2),XYZ(3), MN,
     %                               NBCOFI, FIXA )
C                       NBCOFI LE NOMBRE DE COMPOSANTES FIXEES
                        DO I = 1, NBCOFI
C                          LE NUMERO DE LA COMPOSANTE FIXEE
                           NU = MCN( MN + WUCOFI - 1 + I )
                           FORCE(NU) = FIXA(I) * PENALI
                        ENDDO
                     ENDIF
C
C                    LE COEFFICIENT DU SECOND MEMBRE
                     BE(KM,1) = BE(KM,1) + DELTAK * FORCE(1)
                     BE(KM,2) = BE(KM,2) + DELTAK * FORCE(2)

                  ENDDO
C
               ENDIF
            ENDIF
 130     CONTINUE
C
C        =========================================================
C        CONTRIBUTION DES FORCES OU FIXATION PENALISEE AUX SOMMETS
C        =========================================================
         DO K=1,4
C
C           NO DE POINT DU SOMMET K
            NOOB = NOOBPS(K)
            IF( NOOB .GT. 0 ) THEN
C
C              LE SOMMET K EST UN POINT. EST IL SUPPORT
C              D'UNE SOURCE OU D'UN CONTACT PENALISE?
               IF(LTDEPO(LPSOUR,JEU,NOOB).GT.0 .OR.
     %            LTDEPO(LPCONT,JEU,NOOB).GT.0 .AND. PENALI.NE.0D0) THEN
C
C                 OUI: UN TABLEAU FIXATION PENALISE EXISTE POUR CE POINT
C                 FIXATION PENALISE = PENALI x VALEUR
                  XYZ(1) = X(K,1)
                  XYZ(2) = X(K,2)
                  XYZ(3) = X(K,3)
C
                  MN = LTDEPO(LPSOUR,JEU,NOOB)
                  IF( MN .GT. 0 ) THEN
C
C                    FORCE(2) REQUISE POUR CONDITION NEUMANN OU FOURIER
C                    ATTENTION: ICI OndeR et OndeI SONT AUSSI PASSES
C                               A LA PLACE DE VN le VECTEUR NORMAL!
C                    RECUPERATION DE LA TEMPERATURE AU SOMMET K
                     TEMPEL = VEFtm( K )
                     ONDEPI = WEFtm( K )
                     CALL REFORC( 1, NOOB, 2, XYZ(1),XYZ(2),XYZ(3),
     %                                        TEMPEL, ONDEPI, 0D0,
     %                            MN, FORCE )
                  ENDIF
C
                  MN = LTDEPO(LPCONT,JEU,NOOB)
                  IF( PENALI .GT. 0D0 .AND. MN .GT. 0 ) THEN
C                    FIXATION(2) PENALISEE
                     CALL REFIXA( 1, NOOB, XYZ(1),XYZ(2),XYZ(3), MN,
     %                            NBCOFI, FIXA )
C                    NBCOFI LE NOMBRE DE COMPOSANTES FIXEES
                     DO I = 1, NBCOFI
C                       LE NUMERO DE LA COMPOSANTE FIXEE
                        NU = MCN( MN + WUCOFI - 1 + I )
                        FORCE(NU) = FIXA(I) * PENALI
                     ENDDO
                  ENDIF
C
C                 LE COEFFICIENT DU SECOND MEMBRE
                  BE(K,1) = BE(K,1) + FORCE(1)
                  BE(K,2) = BE(K,2) + FORCE(2)

               ENDIF
            ENDIF
         ENDDO
      ENDIF
C
      RETURN
      END
