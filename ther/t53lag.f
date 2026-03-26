      SUBROUTINE T53LAG( NUELEM, NUNDEL, X,      NBPOLY,
     &                   NS,     NUMISU, NUMASU, NBPOF,  NOPOF,
     &                   POLYF,  DPOLYF, NPIF,   POIDSF,
     &                   NOOBVO, NUMIVO, NUMAVO, LTDEVO,
     &                   COND,   DPFACE, NTDL,   TEMP,   FLUXPT, DREM,
     &                   NBPIFN, NBELEM, NDSM,   NUPTF0, COPNFN,
     &                   FLUXNP, FLUXTO  )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CONTRIBUTION D UNE FACE AU TABLEAU DU CALCUL DES FLUX NORMAUX
C -----    D UN ELEMENT TRIDIMENSIONNEL LAGRANGE DEGRE 2 ISOPARAMETRIQUE
C
C ENTREES:
C --------
C NUMELE : NUMERO DE L'EF DE FACE A TRAITER
C X      : COORDONNEES DES POINTS DE L ELEMENT FINI
C NBPOLY : NOMBRE DE POINTS DE L ELEMENT TRIDIMENSIONNEL COMPLET
C NS     : NUMERO DE SURFACE DE LA FACE
C NBPOF  : NBRE DE POINTS DEFINISSANT L APPLICATION G DE LA FACE
C NOPOF  : NO DANS L ELEMENT DES NBPOF POINTS DE LA FACE
C POLYF  : POLYF(I,J)=VALEUR DU I-EME POLYNOME AU J-EME POINT D INTEGRAT
C DPOLYF : DPOLYF(I,J,L)=VALEUR DE LA I-EME DERIVEE DU J-EME
C                      POLYNOME AU L-EME POINT D INTEGRATION DE LA FACE
C NPIF   : NOMBRE DE POINTS D INTEGRATION SUR CETTE FACE
C POIDSF : LES NPIF POIDS DE LA FORMULE D INTEGRATION DE LA FACE
C NTDL   : NOMBRE TOTAL DE DL DE TEMPERATURE POUR TOUT LE MAILLAGE
C DPFACE : DERIVEES DES POLYNOMES AUX POINTS D INTEGRATION DE LA
C          FACE K DE L EF DE REFERENCE VOLUMIQUE
C TEMP   : TEMPERATURE EN TOUS LES NOEUDS DU MAILLAGE
C FLUXPT : TABLEAU AUXILIAIRE NBPOLY * NPIF
C NBPIFN : NOMBRE TOTAL DE POINTS SUR LES FACES DE L'EF OU SONT CALCULES
C          LES FLUX NORMAUX DE CHALEUR
C( TETRA=>NPITR*NBFACE OU HEXA=>NPIQU*NBFACE OU PENTA=>NPITR*2 + NPIQU*3 )
C NBELEM : NOMBRE D'ELEMENTS FINIS DE CE TYPE
C NDSM   : NOMBRE DE SECONDS MEMBRES
C NUPTF0 : NUMERO DU DERNIER POINT AVANT LE PREMIER DE CETTE FACE
C
C SORTIES:
C --------
C COND   : TENSEUR SYMETRIQUE DE CONDUCTIVITE AU DERNIER POINT
C          D INTEGRATION DE L ELEMENT
C DREM   : LES NBPOLY VALEURS DU TABLEAU DU FLUX A TRAVERS LA FACE
C COPNFN : 2 COORDONNEES DES NBPIFN POINTS ET DE LA NORMALE SUR LES
C          ARETES DES EF OU SONT CALCULES LES FLUX NORMAUX DE TEMPERATURE
C FLUXNP : FLUX NORMAL EN LES NBPIFN POINTS DES FACES DES EF DE CE TYPE
C FLUXTO : FLUX NORMAL SOMME SUR CHAQUE SURFACE DE L'OBJET
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS     OCTOBRE 1994
C23456---------------------------------------------------------------012
      include"./incl/donthe.inc"
      include"./incl/cthet.inc"
      include"./incl/ctemps.inc"
      include"./incl/cnonlin.inc"
C
      include"./incl/pp.inc"
      COMMON            MCN(MOTMCN)
      REAL              RMCN(1)
      DOUBLE PRECISION  DMCN(1)
      EQUIVALENCE      (MCN(1), RMCN(1), DMCN(1))
C
      REAL              X(NBPOLY,3)
      INTEGER           NOPOF(NBPOF)
      DOUBLE PRECISION  POLYF(NBPOF,NPIF),
     &                  DPOLYF(2,NBPOF,NPIF),
     &                  DREM(NBPOLY),
     &                  POIDSF(NPIF),
     &                  DPFACE(3,NBPOLY,NPIF),
     %                  FLUXPT(NBPOLY,NPIF),
     %                  TEMP(NTDL,NDSM)
      DOUBLE PRECISION  FLUXTO( NUMISU:NUMASU , 1:NDSM ),
     %                  FLUXNP( NBPIFN, NBELEM, NDSM )
      DOUBLE PRECISION  DELTA,
     &                  GL(3),
     &                  DGL(2,3),
     &                  COND(6),
     &                  DF(3,3),
     &                  DF1(3,3),
     &                  DN(3),
     &                  D(3), DD(3), DGLN, PROSCD, S
      REAL              COPNFN( 3, NBPIFN, NBELEM, 2 )
      INTEGER           LTDEVO(1:MXDOTH,NUMIVO:NUMAVO)
      INTEGER           NUNDEL( NBELEM, NBPOLY )
C
      DO 90 N=1,NDSM
C
C        MISE A ZERO DU TABLEAU DREM
C        ===========================
         IF( NS .GT. 0 ) THEN
            CALL AZEROD( NBPOLY, DREM )
         ENDIF
C
         IF( MNTIME .GT. 0 ) THEN
C           LE TEMPS DU VECTEUR TEMPERATURE N EST MIS DANS LE COMMON CTEMPS
            TEMPS = RMCN(MNTIME+N)
         ENDIF
C
C        BOUCLE SUR LES POINTS D INTEGRATION DE LA FACE
C        ==============================================
C        RECUPERATION DE LA TEMPERATURE AUX POINTS D'INTEGRATION DE LA FACE
         IF( TESTNL .GT. 0 ) THEN
C           RECUPERATION DE LA TEMPERATURE AUX NOEUDS DE LA FACE
            DO 5 I=1,NBPOF
               DMCN( (MNTHDL-1)/2+I ) =
     &         TEMP( NUNDEL(NUELEM,NOPOF(I)), N )
5           CONTINUE
         ENDIF
C
         DO 20 L=1,NPIF
C
            IF( TESTNL .GT. 0 ) THEN
C              CALCUL DE LA TEMPERATURE AU POINT D'INTEGRATION L
               TEMPEL=PROSCD( POLYF(1,L), MCN(MNTHDL), NBPOF )
            END IF
C
C           CALCUL DES COORDONNEES GL DU POINT D INTEGRATION
C                  DE LA MATRICE JACOBIENNE DGL (R**2 --> R**3)
C                  DU JACOBIEN DELTA
C           ---------------------------------------------------
            CALL E23LAG( NBPOLY, NBPOF, NOPOF,
     &                  POLYF(1,L), DPOLYF(1,1,L), X,
     &                  GL, DGL, DELTA )
C           EN SORTIE GL=LES 3 COORDONNEES GL(1:3) DU POINT D'INTEGRATION L
C
C           LES COORDONNEES DU POINT OU EST CALCULE LE FLUX NORMAL
C           ------------------------------------------------------
            NUPTFN = NUPTF0 + L
            COPNFN( 1, NUPTFN, NUELEM, 1 ) = REAL( GL(1) )
            COPNFN( 2, NUPTFN, NUELEM, 1 ) = REAL( GL(2) )
            COPNFN( 3, NUPTFN, NUELEM, 1 ) = REAL( GL(3) )
C
C           CALCUL DE LA NORMALE UNITAIRE:
C           PRODUIT VECTORIEL(DG/DX1,DG/DX2) ET NORMALISATION
C           -------------------------------------------------
            CALL VECNOR( DGL, DGLN, DN )
C
C           STOCKAGE DES COMPOSANTES DE CE VECTEUR NORMAL UNITAIRE
C           ------------------------------------------------------
            COPNFN( 1, NUPTFN, NUELEM, 2 ) = REAL( DN(1) )
            COPNFN( 2, NUPTFN, NUELEM, 2 ) = REAL( DN(2) )
            COPNFN( 3, NUPTFN, NUELEM, 2 ) = REAL( DN(3) )
C
C           RECHERCHE DU TENSEUR SYMETRIQUE DE CONDUCTIVITE AU POINT GL
C           -----------------------------------------------------------
            CALL RECOND( 4, NOOBVO, 3, GL,
     %                   LTDEVO(LPCOND,NOOBVO), COND )
C
C           CALCUL DE DF(GL) AU POINT D INTEGRATION DE REFERENCE
C           ----------------------------------------------------
            CALL MATJAC( 3, NBPOLY, X, DPFACE(1,1,L), DF )
C
C           CALCUL DE DF(BL) ** -1
C           ----------------------
            CALL M33INV( DF, DELTA, DF1 )
C
C           D = TRANSPOSE(DN) * COND
C           ------------------------
            CALL TABSZD( 1, 3, DN, COND, D )
C
C           DD = D * ((DF) ** -1)
C           ---------------------
            CALL AB0D( 1, 3, 3, D, DF1, DD )
C
C           FLUXPT = DD * DP FACE AU POINT
C           ------------------------------
            CALL AB0D( 1, 3, NBPOLY, DD, DPFACE(1,1,L), FLUXPT(1,L) )
C
            IF( NS .GT. 0 ) THEN
C              CONTRIBUTION A L'INTEGRALE NUMERIQUE DU FLUX NORMAL SUR LA FACE
C              LE POIDS * NORME DU VECTEUR NORMAL A DGL
               DELTA = POIDSF(L) * DGLN
               DO 10 I=1,NBPOLY
                  DREM(I) = DREM(I) - DELTA * FLUXPT(I,L)
 10            CONTINUE
            ENDIF
 20      CONTINUE
C
C        CALCUL DU FLUX NORMAL EN LES POINTS D'INTEGRATION DE LA FACE K
C        ==============================================================
C        LE FLUX NORMAL AUX POINTS D'INTEGRATION NUMERIQUE DE CETTE FACE
         DO 40 L=1,NPIF
            S = 0
            DO 30 I=1,NBPOLY
               S = S - FLUXPT(I,L) * TEMP( NUNDEL(NUELEM,I), N )
 30         CONTINUE
            NUPTFN = NUPTF0 + L
            FLUXNP( NUPTFN, NUELEM, N ) = S
 40      CONTINUE
C
C        LA CONTRIBUTION AUX FLUX NORMAUX DE LA SURFACE NS DE CETTE FACE
         IF( NS .GT. 0 ) THEN
            S = 0
            DO 60 I=1,NBPOLY
               S = S + DREM(I) * TEMP( NUNDEL(NUELEM,I), N )
 60         CONTINUE
C           SOMMATION DES FLUX POUR LA SURFACE FRONTIERE NS
            FLUXTO( NS, N ) = FLUXTO( NS, N ) + S
         ENDIF
 90   CONTINUE
C
      RETURN
      END
