      SUBROUTINE BICGST( NTDL,   NTDLFX,
     S                   LPLIGN, LPCOLO, AG,  B,
     S                   LPLIGC, LPCOLC, AGC, X0,
     S                   R0, R, P, V, S, T, AUX,
     S                   X, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   BUT : RESOLUTION D'UN SYSTEME LINEAIRE SYMETRIQUE
C   ----  PAR L'ALGORITHME ''BI-CONJUGATE GRADIENT STABILISED''
C         AVEC PRECONDITIONNEMENT PAR FACTORISATION INCOMPLETE DE CROUT
C
C ENTREES :
C ---------
C NTDL   : NOMBRE D'INCONNUES ET DE LIGNES DE LA MATRICE
C NTDLFX : NTDLFX(N) =0 SI LE DL N EST LIBRE, 1 SI LE DL N EST FIXE
C          NO TEMOIN D'UN DEGRE DE LIBERTE FIXE'
C
C LPLIGN : POINTEUR SUR LE COEFFICIENT DIAGONAL DE CHAQUE LIGNE DE AG
C LPCOLO : NO COLONNE DE CHAQUE COEFFICIENT STOCKE DE AG
C AG     : MATRICE INITIALE NON FACTORISEE
C
C B      : LE TABLEAU DES NDSM SECONDS MEMBRES
C
C LPLIGC : POINTEUR SUR LE COEFFICIENT DIAGONAL DE CHAQUE LIGNE DE AGC
C LPCOLC : NO COLONNE DE CHAQUE COEFFICIENT STOCKE DE AGC
C AGC    : MATRICE FACTORISEE INCOMPLETE <=> MATRICE DE PRECONDITIONNEMENT
C
C X0     : VECTEUR INITIAL (ITERATION 0)
C R0, R, P, V, S, T, AUX: 7 VECTEURS(NTDL) AUXILIAIRES EN DOUBLE PRECISION
C
C SORTIES :
C ---------
C X      : VECTEUR SOLUTION DU SYSTEME LINEAIRE
C IERR   : 0 SI PAS D'ERREUR RENCONTREE, >0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : Alain PERRONNET LJLL UPMC & St Pierre du Perray     MARS 2009
C23456---------------------------------------------------------------012
      INTEGER           NTDLFX(NTDL)
      INTEGER           LPLIGN(1+NTDL), LPCOLO(*)
      INTEGER           LPLIGC(1+NTDL), LPCOLC(*)
      DOUBLE PRECISION  AG(*), AGC(*), X(NTDL), X0(NTDL), B(NTDL)
      DOUBLE PRECISION  R0(NTDL), R(NTDL), P(NTDL), V(NTDL), S(NTDL),
     %                  T(NTDL), AUX(NTDL)
      DOUBLE PRECISION  PROSCD
      DOUBLE PRECISION  EPS, EPSR, ALPHAK, BETAK, OMEGAK, RNORME,
     %                  TKTK, TKSK, RKR0, RK1R0
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      COMMON / EPSSSS / EPZERO, EPSXYZ
      INTRINSIC         SQRT
C
      IERR = 0
C
C     EPS SEUIL DE CONVERGENCE
      EPS = EPZERO
      WRITE (IMPRIM,10000) EPS
10000 FORMAT(/'DEBUT BICGST BI-CONJUGATE GRADIENT STABILISED avec EPS=',
     %G15.7)
C
C     INITIALISATION DE L'ALGORITHME
C     ******************************
C     (X ) = (X0)
C     (R ) = ((A)) * (X0) - (B)
C     (R0) = ((AGC))-1 * (R)
C     (P ) = (R)
C     (R0) = (R)
C     ------------------------------
      CALL MAGCVX( NTDL, NTDLFX, LPLIGN, LPCOLO, AG, X0, AUX )
      print *
      print *,'RNorme1=', PROSCD( AUX, AUX, NTDL )
C
      DO 10 I = 1,NTDL
         X(I) = X0(I)
         R(I) = B(I) - AUX(I)
 10   CONTINUE
C
ccc      print *
ccc      print *,'MATRICE AGC # L D tL debut bicgst'
ccc      DO i=1,15
ccc         print 10010,(i,LPCOLC(m),agc(m),m=LPLIGC(i)+1, LPLIGC(i+1))
ccc      enddo
ccc10010 format(5('  agc(',i3,',',i3,')=',D15.6))
C
C     R0 = AGC**-1 * R    PRECONDITIONNEMENT
      print *
      print 10011,(m,r(m),m=1,ntdl)
10011 format(5('  r(',i3,')=',D15.6))
      CALL DRCRGC( NTDL, LPLIGC, LPCOLC, AGC, R, R0 )
C
C     INITIALISATIONS DES VECTEURS P et R
      DO 20 I = 1,NTDL
         P(I) = R0(I)
         R(I) = R0(I)
 20   CONTINUE
C
C     RNORME = ( R0, R0 )
      RNORME = PROSCD( R0, R0, NTDL )
      print *
      print *,'RNORME2=',RNORME
      print 10021,(m,r0(m),m=1,ntdl)
10021 format(5('  r0(',i3,')=',D15.6))
C
C     SEUIL DE CONVERGENCE DES ITERATIONS
      EPSR = EPS * EPS * RNORME
      RKR0 = RNORME
C
C     NOMBRE MAXIMUM D'ITERATIONS
      MXITER = NTDL / 5
C
C     ITERATIONS DE L'ALGORITHME BiCG Stabilise
C     *****************************************
      K = 0
C
 30   K = K + 1
C
C     ALPHAK = ( RK, R0 ) / ( A PK, R0 )
C     ----------------------------------
      CALL MAGCVX( NTDL, NTDLFX, LPLIGN, LPCOLO, AG, P, AUX )
      CALL DRCRGC( NTDL, LPLIGC, LPCOLC, AGC, AUX, V )
      ALPHAK = RKR0 / PROSCD( V, R0, NTDL )
C
C     (SK) = (RK) - ALPHAK * VK
C     -------------------------
      DO 40 I = 1,NTDL
         S(I) = R(I) - ALPHAK * V(I)
 40   CONTINUE
C
C     (TK) = AGC**-1 AG SK
C     --------------------
      CALL MAGCVX( NTDL, NTDLFX, LPLIGN, LPCOLO, AG, S, AUX )
      CALL DRCRGC( NTDL, LPLIGC, LPCOLC, AGC, AUX, T )
C
C     OMEGAK+1 = ( TK, SK ) / ( TK, TK )
C     ----------------------------------
      TKSK = PROSCD( T, S, NTDL )
      TKTK = PROSCD( T, T, NTDL )
      OMEGAK = TKSK / TKTK
C
C     XK+1 = XK + OMEGAK SK + ALPHAK PK
C     RK+1 = SK - OMEGAK TK
C     ---------------------------------
      DO 50 I=1,NTDL
         X(I) = X(I) + OMEGAK * S(I) + ALPHAK * P(I)
         R(I) = S(I) - OMEGAK * T(I)
 50   CONTINUE
C
C     TEST DE LA CONVERGENCE
C     ----------------------
      RNORME = PROSCD( R, R, NTDL )
      IF( RNORME .GT. EPSR ) THEN
C
C        BETAK+1 = ALPHAK ( RK+1, R0 ) / ( OMEGAK ( RK, R0 ) )
C        -----------------------------------------------------
         RK1R0 = PROSCD( R, R0, NTDL )
         BETAK = ALPHAK * RK1R0 / ( OMEGAK * RKR0 )
C
C        PK+1 = RK+1 + BETAK+1 ( PK - OMEGAK VK )
C        ----------------------------------------
         DO 60 I=1,NTDL
            P(I) = R(I) + BETAK * ( P(I) - OMEGAK * V(I) )
 60      CONTINUE
C
C        TEST D'ARRET SI MAXIMUM DES ITERATIONS
C        --------------------------------------
         WRITE (IMPRIM,10060) K, RNORME, alphak, betak, omegak
10060    FORMAT(' ITERATION',I6,': RESIDU**2=',D15.6,'   alphak=',D15.6,
     %          '   betak=',D15.6, '  omegak=',D15.6 )
         IF( K .GE. MXITER ) GOTO 999
C
C        ITERATION SUIVANTE
         RKR0 = RK1R0
         GOTO 30
C
      ENDIF
C
C     CONVERGENCE CORRECTE DES ITERATIONS
      WRITE (IMPRIM,10070) K, RNORME
10070 FORMAT(' ITERATION',I6,': RESIDU**2=',D15.6,' CONVERGENCE')
C
C     CALCUL DU RESIDU SUR LE SYSTEME AG X = B
 75   CALL MAGCVX( NTDL, NTDLFX, LPLIGN, LPCOLO, AG, X, AUX )
      DO 80 I = 1,NTDL
         AUX(I) = B(I) - AUX(I)
 80   CONTINUE
      RNORME = SQRT( PROSCD( AUX, AUX, NTDL ) )
C
      WRITE (IMPRIM,10090) RNORME
10090 FORMAT(' FIN BICGST: NORME du RESIDU FINAL=',D15.6)
      RETURN
C
C     NON CONVERGENCE
C     ---------------
 999  WRITE (IMPRIM,10999) K, RNORME
10999 FORMAT(' ITERATION',I6,': NON CONVERGENCE avec RESIDU**2=',D15.6)
      IERR = 1
      GOTO 75
      END
