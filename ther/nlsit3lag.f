      SUBROUTINE NLSIT3LAG( NUELEM, NONOEF, Omega,
     %                      DeltaT, NUTYEL, X, PENALI, NBJEUX, JEU,
     %                      NBSOMT, NOOBPS, NUMIPO, NUMAPO, LTDEPO,
     %                      NBCOTE, NOOBLA, NUMILI, NUMALI, LTDELI,
     %                      NBPOTR, NPITR,  POIDTR, POLYTR, DPOLTR,
     %                      NBPOQU, NPIQU,  POIDQU, POLYQU, DPOLQU,
     %                      NBFACE, NOOBSF, NUMISU, NUMASU, LTDESU,
     %                      NBPOLY, NPI,    POLY,
     %                      NOOBVO, NUMIVO, NUMAVO, LTDEVO,
     %                      F,      POIDEL, DP,
     %                      NBNOMA, Utn,    Utm,
     %                      BE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : NLSE i-time methode: CALCUL DU SECOND MEMBRE BG  avec TESTNL=9
C -----
C REAL PART=Rho/dt Vn + Omega ( x dWn+1m/dy - y dWn+1m/dx )
C          +( N(Vn+1m**2+Wn+1m**2)-N(V0**2-W0**2) ) Vn+1m -Fr
C IMAG PART=Rho/dt Wn - Omega ( x dVn+1m/dy - y dVn+1m/dx )
C          +( N(Vn+1m**2+Wn+1m**2)-N(V0**2-W0**2) ) Wn+1m -Fi

C ENTREES:
C --------
C NUELEM : NUMERO DE L'EF TRAITE
C NONOEF : NUMERO DES NBPOLY NOEUDS DE L'EF NUELEM
C Omega  : VITESSE ANGULAIRE DE LA ROTATION
C DeltaT : PAS DE TEMPS DU SCHEMA IMPLICITE EN TEMPS
C NUTYEL : NUMERO DU TYPE DE L'EF LAGRANGE ISOPARAMETRIQUE
C X      : 3 COORDONNEES DES NBPOLY POINTS DE L'EF
C PENALI : COEFFICIENT DE PENALISATION DE LA CONDITION DE DIRICHLET
C NBJEUX : NOMBRE DE JEUX DE DONNEES
C JEU    : NUMERO DU JEU  DE DONNEES POUR CE CALCUL DE LA MATRICE ELEMENTAIRE
C NBNOMA : NOMBRE DE NOEUDS DU MAILLAGE SUPPORT DE DL
C
C NBSOMT : NOMBRE DE SOMMETS DE L'EF
C NOOBPS : NUMERO DE POINT DES SOMMETS
C NUMIPO : NUMERO MINIMAL DES OBJETS POINTS
C NUMAPO : NUMERO MAXIMAL DES OBJETS POINTS
C LTDEPO : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES THERMIQUES DES POINTS
C
C NBCOTE : NOMBRE DES COTES DE L'EF
C NOOBLA : NOUMERO DES OBJETS LIGNES DES ARETES DE L'ELEMENT
C NUMILI : NUMERO MINIMAL DES OBJETS LIGNES
C NUMALI : NUMERO MAXIMAL DES OBJETS LIGNES
C LTDELI : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES DE THERMIQUE
C          DES OBJETS LIGNES
C
C NBPOTR : NOMBRE DE POLYNOMES DE BASE POUR UNE FACE TRIANGULAIRE DE L'EF
C NPITR  : NOMBRE DE POINTS D INTEGRATION POUR UNE FACE TRIANGULAIRE
C POIDTR : POIDS DES POINTS D INTEGRATION  POUR UNE FACE TRIANGULAIRE
C POLYTR : VALEUR DES POLYNOMES DE BASE AUX POINTS D INTEGRATION
C          D'UNE FACE TRIANGULAIRE
C DPOLTR : VALEUR DU GRADIENT DE CES MEMES POLYNOMES AUX POINTS ...
C
C NBPOQU : NOMBRE DE POLYNOMES DE BASE POUR UNE FACE QUADRANGULAIRE
C NPIQU  : NOMBRE DE POINTS D INTEGRATION POUR UNE FACE QUADRANGULAIRE
C POIDQU : POIDS DES POINTS D INTEGRATION  POUR UNE FACE QUADRANGULAIRE
C POLYQU : VALEUR DES POLYNOMES DE BASE AUX POINTS D INTEGRATION
C          D'UNE FACE QUADRANGULAIRE
C DPOLQU : VALEUR DU GRADIENT DE CES MEMES POLYNOMES AUX POINTS ...
C
C
C NBPOLY : NOMBRE DE POLYNOMES DE L'EF VOLUME
C NPI    : NOMBRE DE POINTS D INTEGRATION NUMERIQUE DE L'EF VOLUME
C POLY   : VALEUR DES POLYNOMES DE BASE AUX POINTS D'INTEGRATION
C          POLY(I,L)=P(I) (XL)
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
C F      : COORDONNEES XX YY ZZ DES NPI POINTS D INTEGRATION DE L'EF
C POIDEL : DELTA * POIDS(NPI) DES NPI POINTS D'INTEGRATION DE L'EF
C DP     : DP(3, NBPOLY, NPI) GRADIENT DES POLYNOMES DE BASE AUX POINTS
C          D'INTEGRATION SUR L'EF COURANT
C
C SORTIE :
C --------
C BE     : BE(NBPOLY,2) LE SECOND MEMBRE ELEMENTAIRE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St Pierre du Perray Novembre 2013
C23456---------------------------------------------------------------012
      PARAMETER        ( MAXPTI=27, MAXPOLY=20 )
C     MAXIMUM POUR LES TYPES D'EF 3D DES POLYNOMES,
C                                    POINTS D'INTEGRATION DU VOLUME
C                                    POINTS D'INTEGRATION D'UNE FACE
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
C
      DOUBLE PRECISION  POIDTR(NPITR),
     %                  POLYTR(NBPOTR,NPITR),
     %                  DPOLTR(2,NBPOTR,NPITR)
      DOUBLE PRECISION  POIDQU(NPIQU),
     %                  POLYQU(NBPOQU,NPIQU),
     %                  DPOLQU(2,NBPOQU,NPIQU)
      DOUBLE PRECISION  POLY(NBPOLY,NPI),
     %                  DP(3,NBPOLY,NPI),
     %                  F(NPI,3),
     %                  POIDEL(NPI)
      DOUBLE PRECISION  BE(NBPOLY,2), PENALI, Omega(3)
      REAL              X(NBPOLY, 3)
      INTEGER           LTDEPO( 1:MXDOTH, 1:NBJEUX, NUMIPO:NUMAPO )
      INTEGER           LTDELI( 1:MXDOTH, 1:NBJEUX, NUMILI:NUMALI )
      INTEGER           LTDESU( 1:MXDOTH, 1:NBJEUX, NUMISU:NUMASU )
      INTEGER           LTDEVO( 1:MXDOTH, 1:NBJEUX, NUMIVO:NUMAVO )
      INTEGER           NONOEF( 1:NBPOLY ),
     %                  NOOBPS( 1:NBSOMT ),
     %                  NOOBLA( 1:NBCOTE ),
     %                  NOOBSF( 1:NBFACE )
      INTEGER           NONOFK(8)
      DOUBLE PRECISION  Utn(NBNOMA,2), Utm(NBNOMA,2)
C
      DOUBLE PRECISION  S, SV, SW, D, DeltaT, XYZ(3), FIXA(2), PROSCD
      DOUBLE PRECISION  Rho(MAXPTI), BETA(MAXPTI), FORCE(2,MAXPTI),
     %                  Vtm(MAXPTI), Wtm(MAXPTI),
     %                  Vtn(MAXPTI), Wtn(MAXPTI),
     %                  F1dYF2dXV(MAXPTI), F1dYF2dXW(MAXPTI),
     %                  VEFtm(MAXPOLY), WEFtm(MAXPOLY),
     %                  VEFtn(MAXPOLY), WEFtn(MAXPOLY),
     %                  NLSECOEF, NLSECOEF0

C     INITIALISATION DE LA VITESSE ANGULAIRE Omega/Z SUPPOSEE CONSTANTE
C     -----------------------------------------------------------------
      IF( NUELEM .EQ. 1 ) THEN
         MN = LTDEVO(LPVIANT,JEU,NOOBVO)
         IF( MN .GT. 0 ) THEN
C            RECUPERATION DE LA VITESSE ANGULAIRE Omega SUPPOSEE CONSTANTE
C            VECTEUR(3) DE VITESSE ANGULAIRE AU 1-ER POINT D'INTEGRATION
             CALL REVIAN( 4, NOOBVO, F(1,1), F(1,2), F(1,3),
     %                    LTDEVO(LPVIANT,JEU,NOOBVO), Omega )
          ELSE
C            PAS DE ROTATION
             Omega(1) = 0D0
             Omega(2) = 0D0
             Omega(3) = 0D0
          ENDIF
      ENDIF
C
C     RECUPERATION DE L'ONDE COMPLEXE AU TEMPS ACTUEL tn+1m et tn
C     AUX NBPOLY DL DE L ELEMENT FINI
C     A PARTIR DU VECTEUR GLOBAL DES DL AUX NOEUDS DE L'ONDE U = V + i W
C     ------------------------------------------------------------------
      DO I=1,NBPOLY
C        NO GLOBAL DU NOEUD I DE L'EF
         N = NONOEF(I)
C        AU TEMPS tn
         VEFtn(I) = Utn(N,1)
         WEFtn(I) = Utn(N,2)
C        AU TEMPS tn+1m
         VEFtm(I) = Utm(N,1)
         WEFtm(I) = Utm(N,2)
      ENDDO
C
C     ==================================
C     CONTRIBUTION DES FORCES VOLUMIQUES
C     ==================================
      DO L = 1, NPI
C        COORDONNEES DU POINT D'INTEGRATION L
         XYZ(1) = F(L,1)
         XYZ(2) = F(L,2)
         XYZ(3) = F(L,3)

C        LA VALEUR INTERPOLEE DES ONDES AU POINT D'INTEGRATION bl
C        CALCUL PARTIE REELLE     a tn+1m AU POINT D'INTEGRATION bl de L'ONDE
         Vtm(L) = PROSCD( POLY(1,L), VEFtm, NBPOLY )
C        CALCUL PARTIE IMAGINAIRE a tn+1m AU POINT D'INTEGRATION bl de L'ONDE
         Wtm(L) = PROSCD( POLY(1,L), WEFtm, NBPOLY )
C        CALCUL PARTIE REELLE a tn AU POINT D'INTEGRATION bl de L'ONDE
         Vtn(L) = PROSCD( POLY(1,L), VEFtn, NBPOLY )
C        CALCUL PARTIE REELLE a tn AU POINT D'INTEGRATION bl de L'ONDE
         Wtn(L) = PROSCD( POLY(1,L), WEFtn, NBPOLY )

C        CALCUL de ( x dVn+1m/dy - y dVn+1m/dx )(bl)
C        CALCUL de ( x dWn+1m/dy - y dWn+1m/dx )(bl)
         SV = 0D0
         SW = 0D0
         IF( Omega(1) .NE. 0D0 ) THEN
            DO J=1,NBPOLY
               D  = F(L,1) * DP(2,J,L) - F(L,2) * DP(1,J,L)
               SV = SV + D * VEFtm(J)
               SW = SW + D * WEFtm(J)
            ENDDO
         ENDIF
C
C        DENSITE DE MASSE
C        ----------------
         MN = LTDEVO(LPMAST,JEU,NOOBVO)
         IF( MN .GT. 0 ) THEN
            CALL REMASS( 4, NOOBVO, 3, XYZ, MN, Rho(L) )
         ELSE
            Rho(L) = 1D0
         ENDIF
C
C        N( |U|**2 ) = r**2/2 + r**4/4 + Beta ( V**2 + W**2 )
C        TERME NON LINEAIRE AU TEMPS ACTUEL
C       -TERME NON LINEAIRE AU TEMPS INITIAL
C        ( ex: N( |U|**2 ) - N( |U0|**2 ) =
C              Beta (V**2+W**2) - Beta (V0**2+W0**2) )
C        ----------------------------------------------------
         MN = LTDEVO(LPCOET,JEU,NOOBVO)
         IF( MN .GT. 0 ) THEN
            TEMPEL  = Vtm(L)
            ONDEPI  = Wtm(L)
            CALL RENLSE( 4, NOOBVO, 3, XYZ, TEMPS, Vtm(L), Wtm(L), MN,
     %                   NLSECOEF )
            TEMPEL0 = Vtn(L)
            ONDEPI0 = Wtn(L)
            CALL RENLSE( 4, NOOBVO, 3, XYZ, TEMPSINI, Vtn(L), Wtn(L),MN,
     %                   NLSECOEF0 )
            BETA(L) = NLSECOEF - NLSECOEF0
         ELSE
            BETA(L) = 0D0
         ENDIF
C
C        FORCE R et I AU SECOND MEMBRE
C        -----------------------------
         MN = LTDEVO(LPSOUR,JEU,NOOBVO)
         IF( MN .GT. 0 ) THEN
            TEMPEL = Vtm(L)
            ONDEPI = Wtm(L)
            CALL REFORC( 4, NOOBVO, 2, XYZ(1),XYZ(2),XYZ(3),
     %                                 TEMPEL, ONDEPI, 0D0,
     %                   MN, FORCE(1,L) )
         ELSE
            FORCE(1,L) = 0D0
            FORCE(2,L) = 0D0
         ENDIF
C
C        COEF = COEF * POIDS(l) * DELTA(Bl)
C        ----------------------------------
         S          = POIDEL(L)
         Rho(L)     = Rho(L)     * S / DeltaT
         BETA(L)    = BETA(L)    * S
         FORCE(1,L) = FORCE(1,L) * S
         FORCE(2,L) = FORCE(2,L) * S

C        OmegaZ ( x dVn+1m/dy - y dVn+1m/dx )(bl)
         F1dYF2dXV(L) = Omega(1) * SV * S
C        OmegaZ ( x dWn+1m/dy - y dWn+1m/dx )(bl)
         F1dYF2dXW(L) = Omega(1) * SW * S
C
      ENDDO
C
C     CALCUL DU VECTEUR ELEMENTAIRE BE
C     --------------------------------
      DO I = 1, NBPOLY

         BE(I,1) = 0D0
         BE(I,2) = 0D0

         DO L = 1, NPI

C           Rho/dt Vn + Omega ( x dWn+1m/dy - y dWn+1m/dx )
C          +( N(Vn+1m**2+Wn+1m**2)-N(V0**2-W0**2) ) Vn+1m -Fr
            BE(I,1) = BE(I,1) + POLY(I,L)
     %              * ( Rho(L) * Vtn(L) + F1dYF2dXW(L)
     %                 +Beta(L) * VEFtm(L) - FORCE(1,L) )

C           Rho/dt Wn - Omega ( x dVn+1m/dy - y dVn+1m/dx )
C          +( N(Vn+1m**2+Wn+1m**2)-N(V0**2-W0**2) ) Wn+1m -Fi
            BE(I,2) = BE(I,2) + POLY(I,L)
     %              * ( Rho(L) * Wtn(L) - F1dYF2dXV(L)
     %                 +Beta(L) * WEFtm(L) - FORCE(2,L) )

         ENDDO

      ENDDO
C
C     =======================================
C     CONTRIBUTIONS DES SOURCES SUR LES FACES
C     =======================================
      DO 30 K=1,NBFACE
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
            IF( IECHAN .EQ. 0 ) GOTO 30
C
C           UN TABLEAU SOURCE OU CONTACT EXISTE POUR CETTE FACE K
C           CALCUL DE LA CONTRIBUTION DE LA FACE K A BE
C           .....................................................
C           RECHERCHE DU NUMERO DES POINTS=NOEUDS DE LA FACE K
            CALL ELNOFA( NUTYEL, K, NBNOFK, NONOFK )
C           NONOFK(I) LES NUMEROS DES NOEUDS DE LA FACE K

            IF( NBSOFA(K) .EQ. 3 ) THEN

C              FACE TRIANGULAIRE
               CALL N43LAG( DeltaT, NBPOLY, X,      NBNOFK, NONOFK,
     %                      POLYTR, DPOLTR, NPITR,  POIDTR,
     %                      NOOB,   NUMISU, NUMASU, NBJEUX, JEU, LTDESU,
     %                      IECHAN, PENALI,  BE )
            ELSE

C              FACE QUADRANGULAIRE
               CALL N43LAG( DeltaT, NBPOLY, X,      NBNOFK, NONOFK,
     %                      POLYQU, DPOLQU, NPIQU,  POIDQU,
     %                      NOOB,   NUMISU, NUMASU, NBJEUX, JEU, LTDESU,
     %                      IECHAN, PENALI,  BE )

            ENDIF
         ENDIF

 30   CONTINUE
C
C
C     ============================================================
C     CONTRIBUTION DES FORCES OU FIXATION PENALISEE SUR LES ARETES
C     ============================================================
      IF( PENALI .NE. 0D0 ) THEN
         DO 250 K=1,NBCOTE
C
C           NO DE LIGNE DE L'ARETE K
            NOOB = NOOBLA(K)
            IF( NOOB .GT. 0 ) THEN
C
C              LE COTE K EST SUR UNE LIGNE. EST IL SUPPORT D'UN CONTACT PENALISE
               MN = LTDELI(LPCONT,JEU,NOOB)
               IF( MN .GT. 0 ) THEN
C
C                 OUI: UN TABLEAU CONTACT PENALISE EXISTE POUR CE COTE
C                 RECHERCHE DES NUMEROS LOCAUX DES NOEUDS DE L'ARETE K
                  NONOFK(1) = NOSOAR(1,K)
                  NONOFK(2) = NOSOAR(2,K)
                  IF( NBNOAR(K) .GT. 0 ) NONOFK(3) = NONOAR(1,K)
C
                  DO 220 I=1,2+NBNOAR(K)
C
C                    LE NUMERO DU I-EME NOEUD DE L'ARETE K
                     NI = NONOFK(I)
C
C                    IL EXISTE UN CONTACT SUR CETTE LIGNE
                     XYZ(1) = X( NI, 1 )
                     XYZ(2) = X( NI, 2 )
                     XYZ(3) = X( NI, 3 )
C
C                    FORCE(2) REQUISE POUR CONDITION NEUMANN OU FOURIER
cccC                 LE VECTEUR NORMAL UNITAIRE
ccc                  VN(1) =  DGL(2) / DELTA
ccc                  VN(2) = -DGL(1) / DELTA
C
C                    ATTENTION: ICI OndeR et OndeI SONT PASSES
C                               A LA PLACE DE VN le VECTEUR NORMAL!
                     MN = LTDELI(LPSOUR,JEU,NOOB)
                     IF( MN .GT. 0 ) THEN
                        TEMPEL = VEFtm( NI )
                        ONDEPI = WEFtm( NI )
                        CALL REFORC( 2, NOOB, 2, XYZ(1),XYZ(2),XYZ(3),
     %                                           TEMPEL, ONDEPI, 0D0,
     %                               MN, FORCE(1,L) )
                     ELSE
                        FORCE(1,L) = 0D0
                        FORCE(2,L) = 0D0
                     ENDIF
C
C                    FIXATION(2) PENALISEE POUR CONDITION DE DIRICHLET
                     MN = LTDELI(LPCONT,JEU,NOOB)
                     IF( MN .GT. 0 ) THEN
                        CALL REFIXA( 2, NOOB, XYZ(1),XYZ(2),XYZ(3), MN,
     %                               NBCOFI, FIXA )
C                       NBCOFI LE NOMBRE DE COMPOSANTES FIXEES
                        DO J = 1, NBCOFI
C                          LE NUMERO DE LA COMPOSANTE FIXEE
                           NU = MCN( MN + WUCOFI - 1 + J )
                           FORCE(NU,L) = FIXA(J) * PENALI
                        ENDDO
                     ENDIF
C
C                    CALCUL DU CONTACT PENALISE = PENALI x TEMPERATURE
                     IF( TESTNL .EQ. 6 ) THEN
                        BE(NI,1) = BE(NI,1) + FORCE(1,L)
                        BE(NI,2) = BE(NI,2) + FORCE(2,L)
                     ELSE IF( TESTNL .EQ. 7 ) THEN
                        BE(NI,1) = BE(NI,1) - FORCE(1,L) * DeltaT
                        BE(NI,2) = BE(NI,2) + FORCE(2,L) * DeltaT
                     ENDIF
C
 220              CONTINUE
C
               ENDIF
            ENDIF
 250     CONTINUE
      ENDIF
C
C     ============================================
C     CONTRIBUTION DU CONTACT PENALISE AUX SOMMETS
C     ============================================
      IF( PENALI .NE. 0D0 ) THEN
         DO K = 1, NBSOMT
C
C           NO DE POINT DU SOMMET K
            NOOB = NOOBPS(K)
            IF( NOOB .GT. 0 ) THEN
C
C              LE SOMMET K EST UN POINT. EST IL SUPPORT D'UN CONTACT PENALISE?
                MN = LTDEPO(LPCONT,JEU,NOOB)
                IF( MN .GT. 0 ) THEN
C
C                 OUI: UN TABLEAU CONTACT PENALISE EXISTE POUR CE POINT
C                 CONTACT PENALISE = PENALI x ONDE
                  TEMPEL = VEFtm( K )
                  ONDEPI = WEFtm( K )
                  XYZ(1) = X(K,1)
                  XYZ(2) = X(K,2)
                  XYZ(3) = X(K,3)
C                 FIXATION(2) PENALISEE
                  CALL REFIXA( 1, NOOB, XYZ(1),XYZ(2),XYZ(3), MN,
     %                         NBCOFI, FIXA )
                  FORCE(1,L) = 0D0
                  FORCE(2,L) = 0D0
C                 NBCOFI LE NOMBRE DE COMPOSANTES FIXEES
                  DO I = 1, NBCOFI
C                    LE NUMERO DE LA COMPOSANTE FIXEE
                     NU = MCN( MN + WUCOFI - 1 + I )
                     FORCE(NU,L) = FIXA(I) * PENALI
                  ENDDO
C
C                 LE COEFFICIENT DU SECOND MEMBRE EST IMPOSE
                  IF( TESTNL .EQ. 6 ) THEN
                     BE(K,1) = BE(K,1) + FORCE(1,L)
                     BE(K,2) = BE(K,2) + FORCE(2,L)
                  ELSE IF( TESTNL .EQ. 7 ) THEN
                     BE(K,1) = BE(K,1) - FORCE(1,L) * DeltaT
                     BE(K,2) = BE(K,2) + FORCE(2,L) * DeltaT
                  ENDIF
C
               ENDIF
            ENDIF
         ENDDO
      ENDIF
C
      RETURN
      END
