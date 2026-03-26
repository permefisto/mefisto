      SUBROUTINE CGSPRG( NTDL,   LPLIGN, LPCOLO, AG,
     S                   LPLIGC, LPDILU, LPCOLC, AGC,
     S                   X, X0, B, AUX )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   BUT : RESOLUTION D'UN SYSTEME LINEAIRE NON SYMETRIQUE
C   ----  PAR L'ALGORITHME ''CONJUGATE GRADIENT SQUARED''
C         AVEC PRECONDITIONNEMENT PAR FACTORISATION INCOMPLETE DE GAUSS
C
C   PARAMETRES D'ENTREE:
C   --------------------
C   NTDL   : NOMBRE D'INCONNUES
C   AG,LPLIGN,LPCOLO : MATRICE NON SYMETRIQUE DU SYSTEME LINEAIRE
C                      ET POINTEURS ASSOCIES
C   AGC,LPLIGC,LPDILU,LPCOLC: MATRICE DE PRECONDITIONNEMENT
C                             ET POINTEURS ASSOCIES
C   X0     : VECTEUR SOLUTION INITIALE
C   B      : SECOND MEMBRE
C   AUX    : VECTEUR AUXILIAIRE
C
C   PARAMETRE DE SORTIE:
C   --------------------
C   X      : VECTEUR SOLUTION DU SYSTEME LINEAIRE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : Pascal JOLY LABORATOIRE D'ANALYSE NUMERIQUE PARIS 6  MAI 1989
C23456---------------------------------------------------------------012
      DOUBLE PRECISION  AG, AGC, X, X0, B, AUX, DSQRT, PROSCD
      DOUBLE PRECISION  EPS, EPSR, ALPHAK, BETAK, R0R, R0AG, RNORME
      DIMENSION         AG(NTDL), AGC(NTDL)
      DIMENSION         X(NTDL), X0(NTDL), B(NTDL), AUX(NTDL, 7)
      DIMENSION         LPLIGN(NTDL+1), LPCOLO(NTDL)
      DIMENSION         LPLIGC(NTDL+1), LPDILU(NTDL), LPCOLC(NTDL)
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      COMMON / EPSSSS / EPZERO, EPSXYZ
C
C     EPS SEUIL DE CONVERGENCE
      EPS = EPZERO
      WRITE (IMPRIM,3000) EPS
C
C     INITIALISATION DE L'ALGORITHME
C     ******************************
C     (X) = (X0)
C     (R) = ((A)) * (X0) - (B)
C     (R) = ((AGC))-1 * (R)
C     (G) = (R)
C     (U) = (R)
C     ------------------------------
      NCODSA = 1
      CALL MAGCVE( 0, 1D0, NTDL, NCODSA, LPLIGN, LPCOLO, AG,
     S             X0 , AUX(1,6) )
C
      DO 1 I=1,NTDL
         X(I)=X0(I)
         AUX(I,6)=B(I)-AUX(I,6)
1     CONTINUE
      CALL DRREGC( NTDL, LPLIGC, LPDILU, LPCOLC,  AGC,
     S             AUX(1,6), AUX(1,1) )
      DO 2 I=1,NTDL
         AUX(I,2)=AUX(I,1)
         AUX(I,3)=AUX(I,1)
         AUX(I,4)=AUX(I,1)
2     CONTINUE
C
      RNORME = PROSCD( AUX(1,1), AUX(1,1), NTDL )
      WRITE(IMPRIM,2000) RNORME
      R0R  = RNORME
      EPSR = EPS * EPS * RNORME
C
C     ITERATIONS DE L'ALGORITHME
C     **************************
      ITEM = NTDL / 5
      K = 0
10    K = K + 1
C
C     ALPHAK = ( R0 , R ) / ( R0 , (A) * G)
C     -------------------------------------
      CALL MAGCVE( 0, 1D0, NTDL, NCODSA, LPLIGN, LPCOLO, AG,
     S             AUX(1,2) , AUX(1,7) )
C
      CALL DRREGC( NTDL, LPLIGC, LPDILU, LPCOLC, AGC,
     S             AUX(1,7), AUX(1,5) )
      R0AG=PROSCD(AUX(1,1),AUX(1,5),NTDL)
      ALPHAK=R0R/R0AG
C
C     (H) = (U) - ALPHAK * ((A)) * (G)
C     (X) = (X) + ALPHAK * ((U) + (H))
C     --------------------------------
      DO 11 I=1,NTDL
         AUX(I,5)=AUX(I,3)-ALPHAK*AUX(I,5)
         AUX(I,6)=AUX(I,3)+AUX(I,5)
         X(I)=X(I)+ALPHAK*AUX(I,6)
11    CONTINUE
C
C     (R)  = (R) - ALPHAK * ((A)) * ((U) + (H))
C     -----------------------------------------
      CALL MAGCVE( 0, 1D0, NTDL, NCODSA, LPLIGN, LPCOLO, AG,
     S             AUX(1,6) , AUX(1,7) )
C
      CALL DRREGC( NTDL, LPLIGC, LPCOLC, LPDILU, AGC,
     S             AUX(1,7), AUX(1,6) )
C
      DO 12 I=1,NTDL
         AUX(I,4)=AUX(I,4)-ALPHAK*AUX(I,6)
12    CONTINUE
C
C     RNORME = ( R , R )
C     ------------------
      RNORME=PROSCD(AUX(1,4),AUX(1,4),NTDL)
C
C     TEST DE LA CONVERGENCE
C     ----------------------
      IF( RNORME .GT. EPSR ) THEN
C
C        BETAK = ( R0 , RK+1 ) / ( R0 , RK )
C        -----------------------------------
         BETAK=1.D0/R0R
         R0R=PROSCD(AUX(1,1),AUX(1,4),NTDL)
         BETAK=BETAK*R0R
C
C        (U) = (R) + BETAK * (H)
C        (G) = (U) + BETAK * ( BETAK * (G) + (H))
C        ----------------------------------------
         DO 13 I=1,NTDL
            AUX(I,3)=AUX(I,4)+BETAK*AUX(I,5)
            AUX(I,2)=AUX(I,3)+BETAK*(BETAK*AUX(I,2)+AUX(I,5))
 13      CONTINUE
C
C        TEST D'ARRET DES ITERATIONS
C        ---------------------------
         IF(K.GE.ITEM) THEN
            GOTO 14
         ELSE
            GOTO 10
         ENDIF
C
      ELSE
C
C        CONVERGENCE
C        -----------
         WRITE (IMPRIM,1000) K
      ENDIF
C
C     VERIFICATION SUR LE SYSTEME INITIAL
C     ***********************************
20    CALL MAGCVE( 0, 1D0, NTDL, NCODSA, LPLIGN, LPCOLO, AG,
     S             X , AUX(1,1) )
      DO 6 I=1,NTDL
         AUX(I,1)=B(I)-AUX(I,1)
6     CONTINUE
      RNORME = PROSCD( AUX(1,1), AUX(1,1), NTDL )
      RNORME = DSQRT( RNORME )
C
      WRITE (IMPRIM,2000) RNORME
      RETURN
C
C     NON CONVERGENCE
C     ---------------
 14   WRITE (IMPRIM,1500) K
      GOTO 20
C
 1000 FORMAT(' ITERATION',T35,I15,' CONVERGENCE')
 1500 FORMAT(' ITERATION',T35,I15,' NON CONVERGENCE')
 2000 FORMAT(' NORME EUCLIDIENNE DU RESIDU',T35,D15.6)
 3000 FORMAT(/' ALGORITHME CONJUGATE GRADIENT SQUARED AVEC EPS='
     S        G15.7,/1X,38(1H-)//)
C
      END
