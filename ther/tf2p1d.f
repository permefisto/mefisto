      SUBROUTINE TF2P1D( X,      NDSM,   NBELEM, NUELEM, NTDL, NUNDEL,
     %                   NOOBLA, NUMILI, NUMALI,
     %                   NOOBSF, NUMISU, NUMASU, LTDESU,
     %                   TEMP,   FLUXPT,
     %                   COPNFN, FLUXNP, FLUXTO )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :     CALCUL DU FLUX NORMAL A CHAQUE ARETE DU TRIANGLE 2P1D
C -----     LAGRANGE DROIT
C           INTEGRATION NUMERIQUE AUX SOMMETS DU TRIANGLE ET DES ARETES
C
C ENTREES :
C ---------
C X      : COORDONNEES X ET Y DES 3 SOMMETS DU TRIANGLE
C NDSM   : NOMBRE DE SECONDS MEMBRES
C NBELEM : NOMBRE D'ELEMENTS FINIS DE CE TYPE
C NUELEM : NUMERO DE L'ELEMENT FINI PARMI CEUX DE CE TYPE
C NTDL   : NOMBRE TOTAL DE DL DE TEMPERATURE POUR TOUT LE MAILLAGE
C NUNDEL : NUMERO DES NOEUDS DES ELEMENTS FINIS
C
C NOOBLA : NUMERO DES LIGNES DES ARETES DE L'EF
C NUMILI : NUMERO MINIMAL DES OBJETS LIGNES
C NUMALI : NUMERO MAXIMAL DES OBJETS LIGNES
C
C NOOBSF : NUMERO DE L'OBJET SURFACE DE CET EF
C NUMISU : NUMERO MINIMAL DES OBJETS SURFACES
C NUMASU : NUMERO MAXIMAL DES OBJETS SURFACES
C LTDESU : TABLEAU DES ADRESSES DU TABLEAU DES DONNEES DE THERMIQUE
C          DES OBJETS SURFACES
C
C TEMP   : TEMPERATURE EN TOUS LES NOEUDS DU MAILLAGE
C FLUXPT : TABLEAU AUXILIAIRE DE 3 * 2 REELS DOUBLE PRECISION
C
C SORTIES:
C --------
C COPNFN : 2 COORDONNEES DES 3 POINTS ET DE LA NORMALE AU MILIEU DES
C          ARETES DES EF OU SONT CALCULES LES FLUX NORMAUX DE TEMPERATURE
C FLUXNP : FLUX NORMAL EN LES 3 POINTS DES ARETES DES EF DE CE TYPE
C FLUXTO : FLUX NORMAL SOMME SUR CHAQUE LIGNE DE L'OBJET
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS    NOVEMBRE 1994
C23456---------------------------------------------------------------012
      include"./incl/donthe.inc"
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
      DOUBLE PRECISION  FLUXPT(3),
     %                  TEMP(NTDL,NDSM),
     %                  COND(6)
      DOUBLE PRECISION  FLUXTO(NUMILI:NUMALI,1:NDSM),
     %                  FLUXNP(3,NBELEM,NDSM)
      REAL              X(3,2),
     %                  COPNFN(3,3,NBELEM,2)
      INTEGER           NOOBLA(3)
      INTEGER           LTDESU(1:MXDOTH,NUMISU:NUMASU)
      INTEGER           NUNDEL(NBELEM,3)
      DOUBLE PRECISION  DN(2), DELTAK
      DOUBLE PRECISION  X21, Y21, X31, Y31, X32, Y32, DELTA
      DOUBLE PRECISION  XYZPI(3), DD, C1, C2, C3
C
      DO 100 K=1,3
C        L'ARETE EST ELLE SUR UNE LIGNE UTILISATEUR?
         NL = NOOBLA( K )
C
C        LE NUMERO DES POINTS DU COTE DE L'ARETE K DU TRIANGLE
         IF( K .NE. 3 ) THEN
            K1 = K + 1
         ELSE
            K1 = 1
         ENDIF
C
C        LE VECTEUR NORMAL UNITAIRE A L'ARETE
         DN(1) = X(K1,2) - X(K, 2)
         DN(2) = X(K ,1) - X(K1,1)
         DELTAK  = SQRT( DN(1)**2 + DN(2)**2 )
         DN(1) = DN(1) / DELTAK
         DN(2) = DN(2) / DELTAK
C
C        LES COORDONNEES DU POINT MILIEU DE L'ARETE OU EST CALCULE LE FLUX NORMA
         COPNFN( 1, K, NUELEM, 1 ) = ( X(K,1)  + X(K1,1) ) * 0.5
         COPNFN( 2, K, NUELEM, 1 ) = ( X(K,2)  + X(K1,2) ) * 0.5
         COPNFN( 3, K, NUELEM, 1 ) = 0
C
C        LES COMPOSANTES DU VECTEUR NORMAL OU EST CALCULE LE FLUX NORMAL
         COPNFN( 1, K, NUELEM, 2 ) = REAL( DN(1) )
         COPNFN( 2, K, NUELEM, 2 ) = REAL( DN(2) )
         COPNFN( 3, K, NUELEM, 2 ) = 0
C
C        LA MATRICE JACOBIENNE DE F
         X21 = X(2,1) - X(1,1)
         X31 = X(3,1) - X(1,1)
         X32 = X(3,1) - X(2,1)
C
         Y21 = X(2,2) - X(1,2)
         Y31 = X(3,2) - X(1,2)
         Y32 = X(3,2) - X(2,2)
C
         DELTA = ( X21 * Y31 - X31 * Y21 ) * 2D0
C
C        RECHERCHE DU TENSEUR SYMETRIQUE DE CONDUCTIVITE AUX 2 SOMMETS
C        DE L'ARETE K DU TRIANGLE
         DO 90 N=1,NDSM
C
            IF( MNTIME .GT. 0 ) THEN
C              LE TEMPS DU VECTEUR TEMPERATURE N EST MIS DANS LE COMMON CTEMPS
               TEMPS = RMCN(MNTIME+N)
            ENDIF
C
            XYZPI(1) = X(K,1)
            XYZPI(2) = X(K,2)
            XYZPI(3) = 0D0
            IF( TESTNL .GE. 1 ) THEN
               TEMPEL=TEMP( NUNDEL(NUELEM,K), N )
            ENDIF
            CALL RECOND( 3, NOOBSF, 3, XYZPI,
     %                   LTDESU(LPCOND,NOOBSF), COND )
            C1 = COND(1)
            C2 = COND(2)
            C3 = COND(3)
C
            XYZPI(1) = X(K1,1)
            XYZPI(2) = X(K1,2)
            XYZPI(3) = 0D0
C
            IF( TESTNL .GE. 1 ) THEN
               TEMPEL=TEMP( NUNDEL(NUELEM,K1), N )
            ENDIF
            CALL RECOND( 3, NOOBSF, 3, XYZPI,
     %                   LTDESU(LPCOND,NOOBSF), COND )
C
C           LE TENSEUR DE CONDUCTIVITE MOYENNE SUR L'ARETE
            C1 = ( C1 + COND(1) ) / DELTA
            C2 = ( C2 + COND(2) ) / DELTA
            C3 = ( C3 + COND(3) ) / DELTA
C
            FLUXPT(1) = DN(1)*(-C1*Y32+C2*X32) + DN(2)*(-C2*Y32+C3*X32)
            FLUXPT(2) = DN(1)*( C1*Y31-C2*X31) + DN(2)*( C2*Y31-C3*X31)
            FLUXPT(3) = DN(1)*(-C1*Y21+C2*X21) + DN(2)*(-C2*Y21+C3*X21)
C
C           LE FLUX NORMAL AU POINT MILIEU DE CETTE ARETE (FLUX CONSTANT)
            DD = 0
            DO 30 I=1,3
               DD = DD - FLUXPT(I) * TEMP( NUNDEL(NUELEM,I), N )
 30         CONTINUE
            FLUXNP( K, NUELEM, N ) = DD
C
C           LA CONTRIBUTION AUX FLUX NORMAUX DE L'ARETE K DE L'EF
            IF( NL .GT. 0 ) THEN
               FLUXTO( NL, N ) = FLUXTO( NL, N ) - DD * DELTAK
            ENDIF
 90      CONTINUE
C
 100  CONTINUE
      END
