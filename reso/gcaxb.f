      SUBROUTINE GCAXB( NTDL, NPDLFX, LPLIGN, LPCOLO, A, B, X0,
     %                  R, V, AV,  X, IERR )
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    METHODE DE GRADIENT CONJUGUE SANS PRECONDITIONNEMENT POUR
C -----    RESOUDRE UN SYSTEME LINEAIRE SYMETRIQUE DEFINI POSITIF A x=b
C          STOCKAGE MORSE SYMETRIQUE DE LA MATRICE A
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
C REMARQUE: IL EST POSSIBLE A L'APPEL D'UTILISER X0=X MAIS PAS B=X!
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY  Decembre 2010
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      COMMON / EPSSSS / EPZERO, EPSXYZ
      COMMON / UNITES / LECTEU, IMPRIM, NUNITE(30)
      INTEGER           NPDLFX(NTDL), LPLIGN(0:NTDL), LPCOLO(*)
      DOUBLE PRECISION  A(*), B(NTDL), X0(NTDL), X(NTDL)
      DOUBLE PRECISION  R(NTDL), V(NTDL), AV(NTDL)
      DOUBLE PRECISION  MUK, LK1, R0R0, RKRK, RK1RK1, AVKVK
      DOUBLE PRECISION  PROSCD
C
C     LE NOMBRE D'ITERATIONS
      K = 0
      PRINT *
C
C     X = X0 SAUF DL FIXE IMPOSE PAR B(No du DL FIXE)
      K  = 0
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
ccc      PRINT *,'gcaxb: (B0,B0)=',PROSCD( B,B, NTDL )
ccc      PRINT *,'gcaxb: (X0,X0)=',PROSCD( X,X, NTDL )
C
C     R0 = V0 = B - A X
      CALL MAGCVX( NTDL, NPDLFX, LPLIGN, LPCOLO, A, X,  AV )
ccc      PRINT *,'gcaxb: (AV0,AV0)=',PROSCD( AV, AV, NTDL )
      DO I=1,NTDL
         R(I) = B(I) - AV(I)
         V(I) = R(I)
      ENDDO
C
C     CARRE DU RESIDU INITIAL
      R0R0   = PROSCD( R, R, NTDL )
      RK1RK1 = R0R0

ccc      PRINT *,'gcaxb: (R0,R0)=',R0R0
ccc      IF( LANGAG .EQ. 0 ) THEN
ccc         PRINT *, 'GC Iteration',K,' (Ax-b,Ax-b)=',R0R0,
ccc     %            ' EpsZero=',EPZERO,' NTDL=',NTDL
ccc      ELSE
ccc         PRINT *, 'CG Iteration',K,' (Ax-b,Ax-b)=',R0R0,
ccc     %            ' EpsZero=',EPZERO,' NTDL=',NTDL
ccc      ENDIF

      IF( R0R0 .LT. 1D-8 ) GOTO 9000
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
      AVKVK = PROSCD( AV, V, NTDL )
ccc      PRINT *,'gcaxb: (AVk,Vk)=',AVKVK
C
C     RKRK = PROSCD( R, R, NTDL )
      MUK  = RKRK / AVKVK
C
C     XK+1 = XK + MUK VK
      CALL CL2VED( NTDL, 1D0, X, MUK, V, X )
C
C     RK+1 = RK - MUK A VK
      CALL CL2VED( NTDL, 1D0, R, -MUK, AV, R )
C
C     CARRE DU RESIDU A L'ITERATION K
      RK1RK1 = PROSCD( R, R, NTDL )

CCC      IF( LANGAG .EQ. 0 ) THEN
CCC         PRINT *, 'GC Iteration',K,' (AVk,Vk)=',AVKVK,
CCC     %            ' ||Axk-b||**2=',RK1RK1
CCC      ELSE
CCC         PRINT *, 'CG Iteration',K,' (AVk,Vk)=',AVKVK,
CCC     %            ' ||Axk-b||**2=',RK1RK1
CCC      ENDIF
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
              PRINT *, 'GC Iteration',K,' de (Ax-b,Ax-b)=',R0R0,
     %                 ' a ',RK1RK1,
     %                 ' => NON CONVERGENCE avec Nb Inconnues=',NTDL
            ELSE
              PRINT *,'CG Iteration',K,' from (Ax-b,Ax-b)=',R0R0,
     %                ' to ',RK1RK1,
     %                ' => NO CONVERGENCE with UNKNOWN Nb=',NTDL
            ENDIF
            GOTO 9999
         ENDIF
         GOTO 50
      ENDIF
C
 9000 IF( LANGAG .EQ. 0 ) THEN
         PRINT *, 'GC Iteration',K,' de (Ax-b,Ax-b)=',R0R0,
     %            ' a ',RK1RK1,' => CONVERGENCE avec Eps=',EPZERO
ccc         PRINT *, 'Nb Inconnues=',NTDL,' La matrice MORSE a',
ccc     %             LPLIGN(NTDL),' coefficients'
      ELSE
         PRINT *, 'CG Iteration',K,' from (Ax-b,Ax-b)=',R0R0,
     %            ' to ',RK1RK1,' => CONVERGENCE with Eps=',EPZERO
ccc         PRINT *, 'Number of UNKNOWN=',NTDL,
ccc     %            ' The CONDENSED MATRIX has',
ccc     %            LPLIGN(NTDL),' coefficients'
      ENDIF
      IERR = 0

ccc      call affvect( 'gcaxb: Solution X=',  20,  X )
 9999 call afl1ve( 'gcaxb: Solution X=', NTDL, X )

      RETURN
      END
