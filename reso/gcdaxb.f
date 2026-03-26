      SUBROUTINE GCDAXB( NTDL, NPDLFX, LPLIGN, LPCOLO, A, B, X0,
     %                   R, V, AV,  X, IERR )
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    METHODE DE GRADIENT CONJUGUE AVEC PRECONDITIONNEMENT PAR
C -----    L'INVERSE DE LA DIAGONALE DE LA MATRICE A POUR
C          RESOUDRE LE SYSTEME LINEAIRE SYMETRIQUE DEFINI POSITIF
C          DiagA-1 A x = DiagA-1 b
C          avec un STOCKAGE MORSE SYMETRIQUE de la MATRICE A
C
C ENTREES:
C --------
C NTDL   : NOMBRE DE LIGNES ET COLONNES DE LA MATRICE et
C          NOMBRE DE COMPOSANTES DES VECTEURS
C NPDLFX : NPDLFX(N)=0 SI LE DL N EST LIBRE
C                    1 SI LE DL N EST FIXE PAR B(N) (=X(N) FINAL)
C          NO TEMOIN D'UN DEGRE DE LIBERTE FIXE'
C          CELA REVIENT A AVOIR LA LIGNE DE A + IDENTITE POUR UN DL FIXE
C LPLIGN : POINTEURS SUR LE DERNIER COEFFICIENT (DIAGONAL) DE
C          CHAQUE LIGNE DE LA MATRICE SYMETRIQUE MORSE A
C LPCOLO : NUMEROS DES COLONNES DES COEFFICIENTS DE LA MATRICE MORSE A
C A      : COEFFICIENTS DE LA MATRICE MORSE A
C B      : VECTEUR SECOND MEMBRE DE NTDL COMPOSANTES
C X0     : VECTEUR INITIAL DE LA PREMIERE ITERATION
C
C R,V,AV : 3 VECTEURS AUXILIAIRES DE NTDL COMPOSANTES
C
C SORTIES:
C --------
C X      : VECTEUR SOLUTION DE NTDL COMPOSANTES APRES CONVERGENCE
C          ATTENTION: X DOIT ETRE DIFFERENT DE B A L'APPEL
C IERR   : 0 SI PAS D'ERREUR
C          1 SI PAS DE CONVERGENCE APRES NTDL ITERATIONS
C
C REMARQUE: IL EST POSSIBLE A L'APPEL D'UTILISER X0=X MAIS B=X INTERDIT!
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC Paris & St PIERRE du PERRAY Mai 2012
C23456---------------------------------------------------------------012
      IMPLICIT NONE
      include"./incl/langue.inc"
      COMMON / EPSSSS / EPZERO, EPSXYZ
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
      INTEGER           NTDL, NPDLFX(NTDL), LPLIGN(0:NTDL), LPCOLO(*)
      DOUBLE PRECISION  A(*), B(NTDL), X0(NTDL), X(NTDL)
      DOUBLE PRECISION  R(NTDL), V(NTDL), AV(NTDL)
      DOUBLE PRECISION  MUK, LK1, R0R0, RKRK, RK1RK1, AVKVK
      DOUBLE PRECISION  PROSCD
      INTEGER           K, I, IERR, LECTEU, IMPRIM, NUNITE
      REAL              EPZERO, EPSXYZ
C
ccc      do i=1,ntdl
ccc         if( b(i) .ne. 0d0 ) print *,'gcdaxb: b(',i,')=',b(i)
ccc      enddo
C
C     LE NOMBRE D'ITERATIONS
      K = 0

      R0R0 = 0D0
      RKRK = 1D100
      DO I=1,NTDL
         R0R0 = R0R0 + A( LPLIGN(I) ) ** 2
         RKRK = MIN( RKRK, A(LPLIGN(I)) )
      ENDDO
      print *,'gcdaxb: ||DiagA||=', SQRT(R0R0), '  MinAii =', RKRK
C
C     X = X0 SAUF DL FIXE IMPOSE PAR B(No du DL FIXE)
      DO I = 1, NTDL
         IF( NPDLFX(I) .EQ. 0 ) THEN
C           DL LIBRE
            X(I) = X0(I)
         ELSE
C           DL FIXE A B(I)
            X(I) = B(I)
         ENDIF
      ENDDO
C
C     PRODUIT A x => AV
      CALL MAGCVX( NTDL, NPDLFX, LPLIGN, LPCOLO, A, X,  AV )
C
      print *,'gcdaxb: ||B0|| =',PROSCD( B,B, NTDL ),
     %              '  ||X0|| =',PROSCD( X,X, NTDL ),
     %             '  ||AV0|| =',PROSCD( AV, AV, NTDL )
C
C     R0 = V0 = DiagA-1 (b - A x)
      DO I=1,NTDL
         R(I) = ( B(I) - AV(I) ) / A( LPLIGN(I) )
         V(I) = R(I)
      ENDDO
C
C     CARRE DU RESIDU INITIAL
      R0R0   = PROSCD( R, R, NTDL )
      print *,'gcdaxb: ||DiagA-1 (b - A x0)|| =',R0R0
      RK1RK1 = R0R0

ccc      IF( LANGAG .EQ. 0 ) THEN
ccc         WRITE(IMPRIM,*)'GCD Iteration',K,' ||DiagA-1(Ax-b)||**2=',R0R0,
ccc     %                  ' EpsZero=',EPZERO,' NTDL=',NTDL
ccc      ELSE
ccc         WRITE(IMPRIM,*)'DCG Iteration',K,' ||DiagA-1(Ax-b)||**2=',R0R0,
ccc     %                  ' EpsZero=',EPZERO,' NTDL=',NTDL
ccc      ENDIF

      IF( R0R0 .LT. 1D-10 ) GOTO 9000
C
C     LE NOMBRE D'ITERATIONS
      K = 1
      RKRK = R0R0
C
C     LES ITERATIONS DU GRADIENT CONJUGUE  ...
C     ===================================
C     PRODUIT A v => Av
 50   CALL MAGCVX( NTDL, NPDLFX, LPLIGN, LPCOLO, A, V,  AV )
C
C     DiagA-1 A v => Av
      DO I=1,NTDL
         AV(I) = AV(I) / A( LPLIGN(I) )
      ENDDO
C
      AVKVK = PROSCD( AV, V, NTDL )
C
C     RKRK = PROSCD( R, R, NTDL )
      MUK  = RKRK / AVKVK
C
C     XK+1 = XK + MUK VK
      CALL CL2VED( NTDL, 1D0, X, MUK, V, X )
C
C     RK+1 = DiagA-1 RK - MUK DiagA-1 A VK
      CALL CL2VED( NTDL, 1D0, R, -MUK, AV, R )
C
C     CARRE DU RESIDU A L'ITERATION K
      RK1RK1 = PROSCD( R, R, NTDL )

ccc      if( k .eq. 1 )  R0R0 = RK1RK1
ccc
ccc      IF( LANGAG .EQ. 0 ) THEN
ccc         WRITE(IMPRIM,*) 'GCD Iteration',K,' ||DiagA-1(Ax-b)||**2=',RK1RK1
ccc      ELSE
ccc         WRITE(IMPRIM,*) 'DCG Iteration',K,' ||DiagA-1(Ax-b)||**2=',RK1RK1
ccc      ENDIF
C
      LK1 = RK1RK1 / RKRK
C
C     VK+1 = RK+1 + LK1 VK
      CALL CL2VED( NTDL, 1D0, R, LK1, V, V )
C
C     TEST D'ARRET DES ITERATIONS
C     ---------------------------
      IF( RK1RK1 .GT. EPZERO * R0R0 ) THEN
         K = K + 1
         RKRK = RK1RK1
         IF( K .GE. NTDL ) THEN
            IERR = 1
            IF( LANGAG .EQ. 0 ) THEN
               WRITE(IMPRIM,*)'GCD Iteration',K,' ||DiagA-1(Ax-b)||**2='
     %             ,RK1RK1,' => NON CONVERGENCE avec Nb Inconnues=',NTDL
            ELSE
               WRITE(IMPRIM,*)'DCG Iteration',K,' ||DiagA-1(Ax-b)||**2='
     %          ,RK1RK1,' => NO CONVERGENCE with UNKNOWN Number=',NTDL
            ENDIF
          GOTO 9999
         ENDIF
         GOTO 50
      ENDIF
C
 9000 IF( LANGAG .EQ. 0 ) THEN
      WRITE(IMPRIM,*)'GCD Iteration',K,' de ||DiagA-1(Ax-b)||**2=',R0R0,
     %               ' a ',RK1RK1,' => CONVERGENCE avec Eps=',EPZERO
ccc         WRITE(IMPRIM,*) 'Nb Inconnues=',NTDL,' La matrice MORSE a',
ccc     %                    LPLIGN(NTDL),' coefficients'
      ELSE
         WRITE(IMPRIM,*)'DCG Iteration',K,
     %                  ' from ||DiagA-1(Ax-b)||**2=',R0R0,
     %                  ' to ',RK1RK1,' => CONVERGENCE with Eps=',EPZERO
ccc         WRITE(IMPRIM,*) 'Number of UNKNOWN=',NTDL,
ccc     %     ' The CONDENSED MATRIX has',LPLIGN(NTDL),' coefficients'
      ENDIF
      IERR = 0
C
 9999 RETURN
      END
