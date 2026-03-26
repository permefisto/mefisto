      SUBROUTINE NEWMBG( NTDL,   ALP,    BET,     GAM,
     &                   AMORTM, AMORTK, DT,
     &                   NORESO, LPDIAG, LPCOLO,
     &                   NCODSM, M,      NCODSK,  K,
     &                   U0,     U1,     B0,  B1, B2,   C )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DU SECOND MEMBRE DU SYSTEME A RESOUDRE PAR
C -----    LA METHODE DE NEWMARK (MATRICES M ET K CONSTANTES)
C          TRANSFERT DES VECTEURS POUR L'INSTANT SUIVANT
C
C ENTREES:
C --------
C NTDL   : DIMENSION DU SYSTEME (NOMBRE TOTAL DE DL)
C
C ALP    : PARAMETRES DE LA METHODE DE NEWMARK DEVANT [M]
C BET    : PARAMETRES DE LA METHODE DE NEWMARK DEVANT [K] et F
C GAM    : PARAMETRES DE LA METHODE DE NEWMARK DEVANT [D] AMORTISSEMENT
C
C AMORTM : COEFFICIENT DE L'AMORTISSEMENT EN MASSE
C AMORTK : COEFFICIENT DE L'AMORTISSEMENT EN RIGIDITE
C DT     : PAS DE TEMPS
C
C NORESO : 1 SI STOCKAGE PROFIL DES MATRICES
C          2 SI STOCKAGE MORSE  DES MATRICES
C LPDIAG : TABLEAU DES POSITIONS DES ELEMENTS DIAGONAUX DES MATRICES
C LPCOLO : TABLEAU DES NUMEROS DE COLONNE DES ELEMENTS DES MATRICES MORSE
C
C NCODSM : 0 SI LA MATRICE M EST DIAGONALE
C          1 SI LA MATRICE EST SYMETRIQUE
C         -1 SI LA MATRICE EST NON-SYMETRIQUE
C M      : MATRICE DE MASSE GLOBALE
C NCODSK : 0 SI LA MATRICE K EST DIAGONALE
C          1 SI LA MATRICE EST SYMETRIQUE
C         -1 SI LA MATRICE EST NON-SYMETRIQUE
C K      : MATRICE DE RIGIDITE GLOBALE
C
C U0     : LE DEPLACEMENT EN Tn-2
C U1     : LE DEPLACEMENT EN Tn-1
C B0     : SECOND MEMBRE EN Tn-2
C B1     : SECOND MEMBRE EN Tn-1
C B2     : SECOND MEMBRE EN Tn
C
C SORTIE :
C --------
C C      : SECOND  MEMBRE DU SYSTEME LINEAIRE A RESOUDRE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : M. A. FERNANDEZ VARELA ANALYSE NUMERIQUE UPMC    JANVIER 1998
C MODIFS : ALAIN PERRONNET        ANALYSE NUMERIQUE UPMC    AVRIL   1999
C-----------------------------------------------------------------------
      INTEGER           LPDIAG(0:NTDL), LPCOLO(1:*)
      DOUBLE PRECISION  AMORTM, AMORTK, DT, ALP(0:2), BET(0:2), GAM(0:2)
      DOUBLE PRECISION  M(*), K(*)
      DOUBLE PRECISION  U0(NTDL), U1(NTDL), C(NTDL)
      DOUBLE PRECISION  B0(NTDL), B1(NTDL), B2(NTDL)
      DOUBLE PRECISION  XA
C
C     CONTRIBUTION DES DEPLACEMENTS A L'INSTANT Tn
C     --------------------------------------------
C     C = XA * M * U0
      XA = -ALP(0) - DT * GAM(0) * AMORTM
      IF( NORESO .EQ. 1 ) THEN
C
C        STOCKAGE PROFIL
C        C = XA * M * U0
         CALL MAPRVE( 0,      XA,  NTDL,
     %                NCODSM, LPDIAG, M,  U0,  C )
C
         XA = - (DT**2) * BET(0) - DT * GAM(0) * AMORTK
C        C = C + XA * K * U0
         CALL MAPRVE( 1,      XA,  NTDL,
     %                NCODSK, LPDIAG, K,  U0,  C )
C
C        C = C + XA * M * U1
         XA = -ALP(1) - DT * GAM(1) * AMORTM
         CALL MAPRVE( 1,      XA,  NTDL,
     %                NCODSM, LPDIAG, M,  U1,  C )
C
C        C = C + XA * K * U1
         XA = - (DT**2) * BET(1) - DT * GAM(1) * AMORTK
         CALL MAPRVE( 1,      XA,  NTDL,
     %                NCODSK, LPDIAG, K,  U1,  C )
C
      ELSE
C
C        STOCKAGE MORSE
C        C = XA * M * U0
         CALL MAGCVE( 0,      XA,  NTDL,
     %                NCODSM, LPDIAG, LPCOLO, M,  U0,  C )
C
         XA = - (DT**2) * BET(0) - DT * GAM(0) * AMORTK
C        C = C + XA * K * U0
         CALL MAGCVE( 1,      XA,  NTDL,
     %                NCODSK, LPDIAG, LPCOLO, K,  U0,  C )
C
C        C = C + XA * M * U1
         XA = -ALP(1) - DT * GAM(1) * AMORTM
         CALL MAGCVE( 1,      XA,  NTDL,
     %                NCODSM, LPDIAG, LPCOLO, M,  U1,  C )
C
C        C = C + XA * K * U1
         XA = - (DT**2) * BET(1) - DT * GAM(1) * AMORTK
         CALL MAGCVE( 1,      XA,  NTDL,
     %                NCODSK, LPDIAG, LPCOLO, K,  U1,  C )
C
      ENDIF
C
C     CONTRIBUTION DES SECONDS MEMBRES DE L'ELASTICITE A L'INSTANT Tn
C     ---------------------------------------------------------------
      XA = DT * DT
      DO 30 I=1,NTDL
         C(I)=C(I) + (BET(0)*B0(I) + BET(1)*B1(I) + BET(2)*B2(I)) * XA
 30   CONTINUE
C
      RETURN
      END
