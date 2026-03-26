      SUBROUTINE GCAXBK( NTDL, NPDLFX, LPLIGN, LPCOLO, A,  B,  X0,
     %                   R, V, AV,  X, KITER )
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
C                   >0 SI LE DL N EST FIXE PAR B(N) (=X(N) FINAL)
C          NO TEMOIN D'UN DEGRE DE LIBERTE FIXE'
C          CELA REVIENT A AVOIR LA LIGNE DE A + IDENTITE POUR UN DL FIXE
C LPLIGN : POINTEURS SUR LE DERNIER COEFFICIENT (DIAGONAL) DE
C          CHAQUE LIGNE DE LA MATRICE SYMETRIQUE MORSE A
C LPCOLO : NUMEROS DES COLONNES DES COEFFICIENTS DE LA MATRICE MORSE A
C A      : COEFFICIENTS DE LA MATRICE MORSE A
C B      : VECTEUR SECOND MEMBRE DE NTDL COMPOSANTES
C X0     : VECTEUR INITIAL DE LA PREMIERE ITERATION DU GC
C
C R,V,AV : 3 VECTEURS AUXILIAIRES DE NTDL COMPOSANTES DOUBLE PRECISION
C
C SORTIE :
C --------
C X      : VECTEUR SOLUTION DE NTDL COMPOSANTES APRES CONVERGENCE
C          ATTENTION: X DOIT ETRE DIFFERENT DE B A L'APPEL
C
C MODIFIE:
C --------
C KITER  : EN ENTREE: NOMBRE MAXIMAL D'ITERATIONS DE GC A FAIRE
C          EN SORTIE: NOMBRE D'ITERATIONS DE GC EFFECTUEES
C                     ET A EXPLOITER LORS DE L'APPEL SUIVANT
C
C          CETTE VALEUR EVITE LA NON CONVERGENCE DU GC LORSQUE LE
C          RESIDU INITIAL EST PETIT ET QUE LA CONVERGENCE EST DEJA ACQUISE
C
C REMARQUES:
C         . IL EST POSSIBLE A L'APPEL D'UTILISER X0=X MAIS PAS B=X!
C         . LE SEUIL DE CONVERGENCE DES ITERATIONS EST EPZERO DEFINI
C           PAR L'UTILISATEUR ET STOCKE DANS LE COMMON /EPSSSS/
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY  Decembre 2012
C MODIFS: ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY  Janvier  2013
C MODIFS: ALAIN PERRONNET             St PIERRE du PERRAY  Avril    2021
C MODIFS: ALAIN PERRONNET             St PIERRE du PERRAY  Mai      2023
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      COMMON / EPSSSS / EPZERO, EPSXYZ
      INTEGER           NPDLFX(NTDL), LPLIGN(0:NTDL), LPCOLO(*), KITER
      DOUBLE PRECISION  A(*), B(NTDL), X0(NTDL), X(NTDL)
      DOUBLE PRECISION  R(NTDL), V(NTDL), AV(NTDL)
      DOUBLE PRECISION  MUK, LK1, B0B0, R0R0, RKRK, RK1RK1, AVKVK
      DOUBLE PRECISION  PROSCD, B0R0MAX, EpsRGC

C     POUR EVITER UN ARRET IMMEDIAT ET PERMETTRE UNE AUGMENTATION
C     DU NOMBRE D'ITERATIONS POUR OBTENIR LA CONVERGENCE
      KITERM = MIN( KITER+NTDL/10, NTDL )
      KITENC = 0

C     SEUIL DE CONVERGENCE DES ITERATIONS DE GRADIENT CONJUGUE
      EpsRGC = EPZERO

C     CARRE DE LA NORME DE B INITIAL
      B0B0 = PROSCD( B, B, NTDL )
ccc      print *,'gcaxbk: Iteration 0  (b0,b0)=',B0B0

C     SI b=0 ALORS x=0 EST SOLUTION du SYSTEME LINEAIRE Ax=b
C     ======================================================
      IF( B0B0 .LE. 0D0 ) THEN
         DO I = 1, NTDL
            X(I) = 0D0
         ENDDO
         K = 0
         R0R0   = 0D0
         RK1RK1 = 0D0
         GOTO 9000
      ENDIF

C     PREPARATION DU GRADIENT CONJUGUE ITERATION=0
C     --------------------------------------------
 10   K = 0

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
ccc      print *,'gcaxbk: ( X0,X0 )=',PROSCD( X,X, NTDL )

C     R0 = V0 = B - A X
      CALL MAGCVX( NTDL, NPDLFX, LPLIGN, LPCOLO, A, X,  AV )
ccc      print *,'gcaxbk: (AV0,AV0)=',PROSCD( AV, AV, NTDL )

      DO I=1,NTDL
         R(I) = B(I) - AV(I)
         V(I) = R(I)
      ENDDO

C     CARRE DU RESIDU INITIAL (A x0 - b0, A x0 - b0)
      R0R0 = PROSCD( R, R, NTDL )
ccc      print *,'gcaxbk: ( R0,R0 )=',R0R0

      B0R0MAX = MAX( B0B0, R0R0 )
ccc      print *,'gcaxbk:  B0R0MAX =',B0R0MAX

      RKRK   = R0R0
      RK1RK1 = R0R0

      IF( R0R0 .LE. EpsRGC * B0B0 ) THEN

C        RESIDU SUFFISAMMENT PETIT. X EST SOLUTION
C        5 ITERATIONS SONT IMPOSEES AVANT ARRET
         MINITER = 5

ccc      print *,'gcaxbk: (b,b)=',B0B0,' (R0,R0)=',R0R0,
ccc     %               ' (R0,R0)/(b,b)=',R0R0/B0B0
ccc      IF( LANGAG .EQ. 0 ) THEN
ccc         PRINT *, 'gcaxbk: GC Iteration',K,' ||Ax-b||**2=',R0R0,
ccc     %   ' NTDL=',NTDL,' EpsRGC=',EpsRGC,' MaxIter=',KITERM
ccc      ELSE
ccc         PRINT *, 'gcaxbk: CG Iteration',K,' ||Ax-b||**2=',R0R0,
ccc     %   ' NTDL=',NTDL,' EpsRGC=',EpsRGC,' MaxIter=',KITERM
ccc      ENDIF
ccc      IF( R0R0 .LE. 1D-10 * B0B0 ) GOTO 9000
ccc      IF( R0R0 .LE. 1D-8  * B0B0 ) GOTO 9000

ccc      GOTO 9000

      ELSE

C        PAS DE MINIMUM DES ITERATIONS
         MINITER = NTDL

      ENDIF


C     LES ITERATIONS K DU GRADIENT CONJUGUE  ...
C     =====================================
 50   K = K + 1

C     PRODUIT A v => Av
      CALL MAGCVX( NTDL, NPDLFX, LPLIGN, LPCOLO, A, V,  AV )

      AVKVK = PROSCD( AV, V, NTDL )
ccc      print *,'gcaxbk: (AVk,Vk)=',AVKVK

C     RKRK = PROSCD( R, R, NTDL )
      MUK  = RKRK / AVKVK

C     XK+1 = XK + MUK VK
      CALL CL2VED( NTDL, 1D0, X, MUK, V, X )

C     RK+1 = RK - MUK A VK
      CALL CL2VED( NTDL, 1D0, R, -MUK, AV, R )

C     CARRE DU RESIDU A L'ITERATION K
      RK1RK1 = PROSCD( R, R, NTDL )

ccc      IF( LANGAG .EQ. 0 ) THEN
ccc         PRINT *, 'GC Iteration',K,' (AVk,Vk)=',AVKVK,
ccc     %            ' (Axk-b,Axk-b)=',RK1RK1
ccc      ELSE
ccc         PRINT *, 'CG Iteration',K,' (AVk,Vk)=',AVKVK,
ccc     %            ' (Axk-b,Axk-b)=',RK1RK1
ccc      ENDIF

      LK1 = RK1RK1 / RKRK

C     VK+1 = RK+1 + LK1 VK
      CALL CL2VED( NTDL, 1D0, R, LK1, V, V )

C     TEST D'ARRET DES ITERATIONS
C     ---------------------------
ccc      IF( K .EQ. 1 .AND. R0R0 .LE. 1D-8 * B0B0 ) GOTO 9000
ccc      IF( K .LT. 8 .OR.
ccc     %  (RK1RK1 .GT. 1D-10*B0B0 .AND. RK1RK1 .GT. EpsRGC*R0R0) )THEN

C     LE NOMBRE MINIMAL D'ITERATIONS IMPOSEES EST IL ATTEINT?
      IF( K .EQ. MINITER ) GOTO 9000

cccC     64 ITERATIONS SONT IMPOSEES AVANT LE TEST D'ARRET DES ITERATIONS
ccc      IF( K .LT. 64 .OR. RK1RK1 .GT. EpsRGC*B0R0MAX ) THEN
ccc      IF( RK1RK1 .GT. EpsRGC*B0R0MAX ) THEN

      IF( RK1RK1 .GT. EpsRGC * R0R0 ) THEN

C        A L'ITERATION K+1  PAS ENCORE DE CONVERGENCE DU RESIDU
         RKRK = RK1RK1

         IF( K .GE. KITERM  .OR.  RK1RK1 .GT. 100 * R0R0 ) THEN

C           NON CONVERGENCE ou DIVERGENCE du RESIDU
            IF( LANGAG .EQ. 0 ) THEN
               PRINT *, 'gcaxbk: Iteration',K,' de (b0,b0)=',B0B0,
     %                ' (Ax0-b0,Ax0-b0)=',R0R0,' a (Ax-b,Ax-b)=',RK1RK1,
     %                  ' => NON CONVERGENCE pour EpsRGC=',EpsRGC,
     %                  ' RK1RK1/RKRK=',RK1RK1/RKRK,' Nb DL=',NTDL
            ELSE
               PRINT *, 'gcaxbk: Iteration',K,' de (b0,b0)=',B0B0,
     %               ' (Ax0-b0,Ax0-b0)=',R0R0,' to (Ax-b,Ax-b)=',RK1RK1,
     %                 ' => NO CONVERGENCE with EpsRGC=',EpsRGC,
     %                 ' RK1RK1/RKRK=',RK1RK1/RKRK,' DoF Nb=',NTDL
            ENDIF

C           LE SEUIL D'ARRET DES ITERATIONS DU GC EST RELEVE
            KITENC = KITENC + K
            EpsRGC = EpsRGC * 5D0
            IF( EpsRGC .GT. 0.1D0 ) THEN
C              SEUIL DEVENANT TROP GRAND
               GOTO 9999
            ELSE
C              SEUIL RAISONABLE, REDEPART DU GC AVEC X0 LE VECTEUR X ACTUEL
               KITERM = NTDL
               IF( R0R0 .LT. RK1RK1 ) THEN
C                 REDEPART AVEC X0 = 0 ET EpsRGC PLUS PETIT
                  DO I = 1, NTDL
                     X0(I) = 0D0
                  ENDDO
                  GOTO 10
               ELSE
                  K = 1
                  GOTO 50
               ENDIF
            ENDIF

         ENDIF

C        UNE ITERATION K DE PLUS A FAIRE
         GOTO 50

      ENDIF

 9000 IF( LANGAG .EQ. 0 ) THEN
         PRINT *,'gcaxbk: Iteration',K,' de (b0,b0)=',B0B0,
     %           ' (Ax0-b0,Ax0-b0)=',R0R0,
     %           ' a (Ax-b,Ax-b)=',RK1RK1,
     %           ' => CONVERGENCE pour EpsRGC=',EpsRGC
ccc         PRINT *, 'Nb Inconnues=',NTDL,' La matrice MORSE a',
ccc     %             LPLIGN(NTDL),' coefficients'
      ELSE
         PRINT *,'gcaxbk: Iteration',K,' from (b0,b0)=',B0B0,
     %           ' (Ax0-b0,Ax0-b0)=',R0R0,
     %           ' to (Ax-b,Ax-b)=',RK1RK1,
     %           ' => CONVERGENCE with EpsRCG=',EpsRGC
ccc         PRINT *, 'Number of UNKNOWN=',NTDL,
ccc     %            ' The CONDENSED MATRIX has',
ccc     %            LPLIGN(NTDL),' coefficients'
      ENDIF

C     K NOMBRE DES ITERATIONS DE GC EFFECTUEES
C     MISE A JOUR DU NOMBRE D'ITERATIONS GC A FAIRE LA PROCHAINE FOIS
      KITER = MIN( K + KITENC, NTDL )

ccc      call affvect('End gcaxbk: Solution X=',   5,  X )
ccc      call afl1ve( 'End gcaxbk: Solution X=', NTDL, X )

 9999 RETURN
      END
