      SUBROUTINE TF3P1D( X,      NDSM,  NUELEM, NBELEM, DP,
     %                   NTDL,   NUNDEL,
     %                   NOOBSF, NUMISU, NUMASU,
     %                   NOOBVO, NUMIVO, NUMAVO, LTDEVO,
     %                   TEMP,   FLUXPT,
     %                   COPNFN, FLUXNP, FLUXTO )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     CALCUL DU FLUX NORMAL DE CHALEUR AU BARYCENTRE DE CHAQUE
C -----     FACE D'UN TETRAEDRE 3P1D LAGRANGE
C
C ENTREES :
C ---------
C X      : COORDONNEES X Y Z DES 4 POINTS DE L'EF VOLUMIQUE
C NDSM   : NOMBRE DE SECONDS MEMBRES
C NOEF   : NUMERO DU TYPE DE L'EF
C NUELEM : NUMERO DE L'ELEMENT FINI PARMI CEUX DE CE TYPE(DANS NPEF"TYPE)
C NBELEM : NOMBRE D'ELEMENTS FINIS DE CE TYPE
C NTDL   : NOMBRE TOTAL DE DL DE TEMPERATURE POUR TOUT LE MAILLAGE
C NUNDEL : NUMERO DES NOEUDS DES ELEMENTS FINIS
C
C NOOBSF : NUMERO DE LA SURFACE DE CHAQUE FACE DE L'EF
C NUMISU : NUMERO MINIMAL DES OBJETS SURFACES
C NUMASU : NUMERO MAXIMAL DES OBJETS SURFACES
C
C NOOBVO : NUMERO DE L'OBJET VOLUME DE CET EF
C NUMIVO : NUMERO MINIMAL DES OBJETS VOLUMES
C NUMAVO : NUMERO MAXIMAL DES OBJETS VOLUMES
C LTDEVO : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES CONDUCTIVITE
C          DES OBJETS VOLUMES
C
C TEMP   : TEMPERATURE EN TOUS LES NOEUDS DU MAILLAGE
C FLUXPT : TABLEAU AUXILIAIRE DE 4 REELS DOUBLE PRECISION
C
C SORTIES:
C --------
C COPNFN : 2 COORDONNEES DES 4 POINTS ET DE LA NORMALE SUR LES
C          ARETES DES EF OU SONT CALCULES LES FLUX NORMAUX DE TEMPERATURE
C FLUXNP : FLUX NORMAL EN LES 4 POINTS DES FACES DES EF DE CE TYPE
C FLUXTO : FLUX NORMAL SOMME SUR CHAQUE LIGNE DE L'OBJET
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     OCTOBRE 1994
C23456---------------------------------------------------------------012
      include"./incl/donthe.inc"
      include"./incl/ponoel.inc"
      include"./incl/ctemps.inc"
      include"./incl/cthet.inc"
      include"./incl/cnonlin.inc"
C
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE      (MCN(1), RMCN(1), DMCN(1))
C
      REAL              X(4,3),
     %                  COPNFN( 3, 4, NBELEM, 2 )
      INTEGER           LTDEVO( 1:MXDOTH , NUMIVO:NUMAVO )
      INTEGER           NONOFK(3),
     %                  NOOBSF(1:4)
      DOUBLE PRECISION  FLUXPT(4),
     %                  TEMP(NTDL,NDSM),
     %                  COND(6)
      DOUBLE PRECISION  FLUXTO( NUMISU:NUMASU , 1:NDSM ),
     %                  FLUXNP( 4, NBELEM, NDSM )
      INTEGER           NUNDEL( NBELEM, 4 ), NBSM
      DOUBLE PRECISION  DP(3,4), XYZPI(3), DD, D(3), DELTAK,
     %                  DGL(2,3), C(6), DN(3)
C
C     CONTRIBUTION DES FACES
C     ======================
      DO 100 K=1,4
C
C        LE NUMERO EVENTUEL DE LA SURFACE DE CETTE FACE
         NS = NOOBSF(K)
C
C        LA FACE EST SUR UNE SURFACE UTILISATEUR
C        RECHERCHE DU NUMERO DES POINTS=NOEUDS DE LA FACE K
         CALL ELNOFA( 19, K, NBNOFK, NONOFK )
C        NONOFK(I) LES NUMEROS DES NOEUDS DE LA FACE K
C                        CE SONT AUSSI LES POINTS D'INTEGRATION
C
C        RECHERCHE DU JACOBIEN DE G
         N = NONOFK(1)
         DGL(1,1) = X( NONOFK(2), 1 ) - X( N, 1 )
         DGL(2,1) = X( NONOFK(3), 1 ) - X( N, 1 )
         DGL(1,2) = X( NONOFK(2), 2 ) - X( N, 2 )
         DGL(2,2) = X( NONOFK(3), 2 ) - X( N, 2 )
         DGL(1,3) = X( NONOFK(2), 3 ) - X( N, 3 )
         DGL(2,3) = X( NONOFK(3), 3 ) - X( N, 3 )
C
C        CALCUL DES FLUX NORMAUX AU BARYCENTRE DE LA FACE
C        LES COORDONNEES DU BARYCENTRE DE LA FACE K
         COPNFN( 1, K, NUELEM, 1 ) =
     %                (X(N,1)+X(NONOFK(2),1)+X(NONOFK(3),1)) / 3
         COPNFN( 2, K, NUELEM, 1 ) =
     %                (X(N,2)+X(NONOFK(2),2)+X(NONOFK(3),2)) / 3
         COPNFN( 3, K, NUELEM, 1 ) =
     %                (X(N,3)+X(NONOFK(2),3)+X(NONOFK(3),3)) / 3
C
C        CALCUL DE LA NORMALE: PRODUIT VECTORIEL(DG/DX1,DG/DX2)
         DN(1) = DGL(1,2) * DGL(2,3) - DGL(2,2) * DGL(1,3)
         DN(2) = DGL(1,3) * DGL(2,1) - DGL(2,3) * DGL(1,1)
         DN(3) = DGL(1,1) * DGL(2,2) - DGL(2,1) * DGL(1,2)
C
C        LES COMPOSANTES DU VECTEUR NORMAL OU EST CALCULE LE FLUX NORMAL
         DELTAK = SQRT( DN(1) ** 2 + DN(2) ** 2 + DN(3) ** 2 )
         DN(1)  = DN(1) / DELTAK
         DN(2)  = DN(2) / DELTAK
         DN(3)  = DN(3) / DELTAK
         COPNFN( 1, K, NUELEM, 2 ) = REAL( DN(1) )
         COPNFN( 2, K, NUELEM, 2 ) = REAL( DN(2) )
         COPNFN( 3, K, NUELEM, 2 ) = REAL( DN(3) )
C
C        RECHERCHE DU TENSEUR SYMETRIQUE DE CONDUCTIVITE AU BARYCENTRE DE LA FAC
         CALL AZEROD( 6, C )
         DO 50 NBSM=1,NDSM
C
            IF( MNTIME .GT. 0 ) THEN
C              LE TEMPS DU VECTEUR TEMPERATURE N EST MIS DANS LE COMMON CTEMPS
               TEMPS = RMCN(MNTIME+N)
            ENDIF
C
            DO 20 L=1,3
C
C              LE NUMERO DU SOMMET L DE LA FACE DANS LE TETRAEDRE
               N  = NONOFK( L )
               XYZPI(1) = X(N,1)
               XYZPI(2) = X(N,2)
               XYZPI(3) = X(N,3)
C
               IF( TESTNL .GE. 1 ) THEN
C                 LES COEFFICIENTS DEPENDENT DE LA TEMPERATURE
                  TEMPEL=TEMP( NUNDEL(NUELEM,N), NBSM )
               END IF
C
               CALL RECOND( 4, NOOBVO, 3, XYZPI,
     %                      LTDEVO(LPCOND,NOOBVO), COND )
C
               DO 10 I=1,6
                  C(I) = C(I) + COND(I)
 10            CONTINUE
C
 20         CONTINUE
C
C           C = CONDUCTIVITE * POIDS
            DO 30 I=1,6
               C(I) = C(I) / 3D0
 30         CONTINUE
C
C           D = TRANSPOSE(DN) * C
            CALL TABSZD( 1, 3, DN, C, D )
C
C           FLUXPT = D * DP
            CALL AB0D( 1, 3, 4, D, DP, FLUXPT )
C
C           LE FLUX NORMAL AU BARYCENTRE DE CETTE FACE K
            DD = 0
            DO 40 I=1,4
               DD = DD - FLUXPT(I) * TEMP( NUNDEL(NUELEM,I),NBSM)
 40         CONTINUE
            FLUXNP( K, NUELEM, NBSM ) = DD
C
C           LA CONTRIBUTION AUX FLUX NORMAUX DE LA SURFACE NS DE LA FACE K DE L'
            IF( NS .GT. 0 ) THEN
C              SOMMATION DES FLUX POUR LA SURFACE FRONTIERE NS
               FLUXTO( NS, NBSM ) = FLUXTO( NS, NBSM ) - DD * DELTAK
            ENDIF
 50      CONTINUE
C
 100  CONTINUE
      END
