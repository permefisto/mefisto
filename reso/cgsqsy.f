      SUBROUTINE CGSQSY( NTDL,   NTDLFX,
     S                   LPLIGN, LPCOLO, AG,  B,
     S                   LPLIGC, LPCOLC, AGC, X0, AUX,
     S                   X,      IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   BUT : RESOLUTION D'UN SYSTEME LINEAIRE SYMETRIQUE
C   ----  PAR L'ALGORITHME ''CONJUGATE GRADIENT SQUARED''
C         AVEC PRECONDITIONNEMENT PAR FACTORISATION INCOMPLETE DE CROUT
C         MODIFICATION DE CSGPRG EN CHANGEANT GAUSS PAR CROUT
C
C   PARAMETRES D'ENTREE:
C   --------------------
C   NTDL   : NOMBRE D'INCONNUES
C   NTDLFX : NTDLFX(N) =0 SI LE DL N EST LIBRE, 1 SI LE DL N EST FIXE
C            NO TEMOIN D'UN DEGRE DE LIBERTE FIXE'
C   AG,LPLIGN,LPCOLO : MATRICE NON SYMETRIQUE DU SYSTEME LINEAIRE
C                      ET POINTEURS ASSOCIES
C   B      : SECOND MEMBRE
C   AGC,LPLIGC,LPCOLC: MATRICE DE PRECONDITIONNEMENT
C                      ET POINTEURS ASSOCIES
C   X0     : VECTEUR SOLUTION INITIALE
C   AUX    : 7 VECTEURS AUXILIAIRES
C
C   PARAMETRE DE SORTIE:
C   --------------------
C   X      : VECTEUR SOLUTION DU SYSTEME LINEAIRE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : Pascal JOLY LABORATOIRE D'ANALYSE NUMERIQUE PARIS 6  MAI 1989
C MODIFS : Alain PERRONNET LJLL UPMC & St Pierre du Perray     MARS 2009
C23456---------------------------------------------------------------012
      INTEGER           NTDLFX(NTDL)
      INTEGER           LPLIGN(1+NTDL), LPCOLO(*)
      INTEGER           LPLIGC(1+NTDL), LPCOLC(*)
      DOUBLE PRECISION  AG(*), AGC(*), X(NTDL), X0(NTDL), B(NTDL),
     %                  AUX(NTDL,7)
      DOUBLE PRECISION  DSQRT, PROSCD
      DOUBLE PRECISION  EPS, EPSR, ALPHAK, BETAK, R0R, R0AG, RNORME
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      COMMON / EPSSSS / EPZERO, EPSXYZ
C
      IERR = 0
C
C     EPS SEUIL DE CONVERGENCE
      EPS = EPZERO
      WRITE (IMPRIM,3000) EPS
 3000 FORMAT(/' DEBUT CGSQSY CONJUGATE GRADIENT SQUARED avec EPS=',
     %G15.7)
C
C     INITIALISATION DE L'ALGORITHME
C     ******************************
C     (X) = (X0)
C     (R) = ((A)) * (X0) - (B)
C     (R) = ((AGC))-1 * (R)
C     (G) = (R)
C     (U) = (R)
C     ------------------------------
      CALL MAGCVX( NTDL, NTDLFX, LPLIGN, LPCOLO, AG, X0, AUX(1,6) )
ccc      print *
ccc      print *,'RNORME1=', PROSCD( AUX(1,6), AUX(1,6), NTDL )
C
      DO 1 I = 1,NTDL
         X(I) = X0(I)
         AUX(I,6) = B(I) - AUX(I,6)
1     CONTINUE
c
ccc      print *
ccc      print *,'MATRICE AGC debut cgsqsy'
ccc      DO i=ntdl-32,ntdl
ccc         print 12345,(i,LPCOLC(m),agc(m),m=LPLIGC(i)+1, LPLIGC(i+1))
ccc      enddo
ccc12345 format(5('  agc(',i3,',',i3,')=',D15.6))
c
ccc      print *
ccc      print 99000,(m,aux(m,6),m=1,ntdl)
ccc99000 format(5('  aux(',i3,',6)=',D15.6))
      CALL DRCRGC( NTDL, LPLIGC, LPCOLC, AGC, AUX(1,6),  AUX(1,1) )
C
      DO 2 I = 1,NTDL
         AUX(I,2) = AUX(I,1)
         AUX(I,3) = AUX(I,1)
         AUX(I,4) = AUX(I,1)
2     CONTINUE
C
      RNORME = PROSCD( AUX(1,1), AUX(1,1), NTDL )
ccc      print *
ccc      print *,'RNORME2=',RNORME
ccc      print 99001,(m,aux(m,1),m=1,ntdl)
ccc99001 format(5('  aux(',i3,',1)=',D15.6))
      R0R  = RNORME
      EPSR = EPS * EPS * RNORME
C
C     ITERATIONS DE L'ALGORITHME
C     **************************
      MXITER = NTDL / 5
      K = 0
10    K = K + 1
C
C     ALPHAK = ( R0 , R ) / ( R0 , (A) * G)
C     -------------------------------------
      CALL MAGCVX( NTDL, NTDLFX, LPLIGN, LPCOLO, AG, AUX(1,2), AUX(1,7))
ccc      print 99002,(m,aux(m,7),m=1,ntdl)
ccc99002 format(5('  aux(',i3,',7)=',D15.6))
      CALL DRCRGC( NTDL, LPLIGC, LPCOLC, AGC, AUX(1,7),  AUX(1,5) )
ccc      print 99003,(m,aux(m,5),m=1,ntdl)
ccc99003 format(5('  aux(',i3,',5)=',D15.6))
      R0AG = PROSCD( AUX(1,1), AUX(1,5), NTDL )
      ALPHAK = R0R / R0AG
C
C     (H) = (U) - ALPHAK * ((A)) * (G)
C     (X) = (X) + ALPHAK * ((U) + (H))
C     --------------------------------
      DO 11 I = 1,NTDL
         AUX(I,5) = AUX(I,3) - ALPHAK * AUX(I,5)
         AUX(I,6) = AUX(I,3) + AUX(I,5)
         X(I)     = X(I)     + ALPHAK * AUX(I,6)
11    CONTINUE
C
C     (R)  = (R) - ALPHAK * ((A)) * ((U) + (H))
C     -----------------------------------------
      CALL MAGCVX( NTDL, NTDLFX, LPLIGN, LPCOLO, AG, AUX(1,6), AUX(1,7))
      CALL DRCRGC( NTDL, LPLIGC, LPCOLC, AGC, AUX(1,7),  AUX(1,6) )
      DO 12 I = 1,NTDL
         AUX(I,4) = AUX(I,4) - ALPHAK * AUX(I,6)
12    CONTINUE
C
C     RNORME = ( R , R )
C     ------------------
      RNORME = PROSCD( AUX(1,4), AUX(1,4), NTDL )
C
C     TEST DE LA CONVERGENCE
C     ----------------------
      IF( k .le. 12 .or. RNORME .GT. EPSR ) THEN
C
C        BETAK = ( R0 , RK+1 ) / ( R0 , RK )
C        -----------------------------------
         BETAK = R0R
         R0R   = PROSCD( AUX(1,1), AUX(1,4), NTDL )
         BETAK = R0R / BETAK
C
         WRITE (IMPRIM,1100) K, RNORME, alphak, betak
ccc      print 99005,(m,x(m),m=1,ntdl)
ccc99005 format(5('  x(',i3,')=',D15.6))
C
C        (U) = (R) + BETAK * (H)
C        (G) = (U) + BETAK * ( BETAK * (G) + (H))
C        ----------------------------------------
         DO 13 I = 1,NTDL
            AUX(I,3) = AUX(I,4) + BETAK * AUX(I,5)
            AUX(I,2) = AUX(I,3) + BETAK * ( AUX(I,5) + BETAK*AUX(I,2) )
 13      CONTINUE
C
C        TEST D'ARRET DES ITERATIONS
C        ---------------------------
         IF( K .GE. MXITER ) GOTO 999
         GOTO 10
C
      ELSE
C
C        CONVERGENCE
C        -----------
         WRITE (IMPRIM,1000) K, RNORME
      ENDIF
C
C     VERIFICATION SUR LE SYSTEME INITIAL
C     ***********************************
 20   CALL MAGCVX( NTDL, NTDLFX, LPLIGN, LPCOLO, AG, X, AUX(1,1) )
      DO 6 I = 1,NTDL
         AUX(I,1) = B(I) - AUX(I,1)
 6    CONTINUE
      RNORME = PROSCD( AUX(1,1), AUX(1,1), NTDL )
      RNORME = DSQRT( RNORME )
C
      WRITE (IMPRIM,2000) RNORME
      RETURN
C
C     NON CONVERGENCE
C     ---------------
 999  WRITE (IMPRIM,1500) K, RNORME
      IERR = 1
      GOTO 20
C
 1100 FORMAT(' ITERATION',I6,': RESIDU**2=',D15.6,'   alphak=',D15.6,
     %'   betak=',D15.6)
 1000 FORMAT(' ITERATION',I6,': RESIDU**2=',D15.6,' CONVERGENCE')
 1500 FORMAT(' ITERATION',I6,': NON CONVERGENCE avec RESIDU**2=',D15.6)
 2000 FORMAT(' FIN CGSQSY: NORME EUCLIDIENNE du RESIDU=',D15.6)
      END
