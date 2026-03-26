      SUBROUTINE GCPRCH( NTDL,   NDSM,   NBDIR,
     &                   LPLIGN, LPCOLO, AG, B,
     &                   LPLIGC, LPCOLC, AGC,
     &                   X0, X, R, Z, DIR, ADIR, DAD, BETA,
     &                   U, IERR )
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    METHODE DE GRADIENT CONJUGUE AVEC PRECONDITIONNEMENT PAR
C -----    FACTORISATION INCOMPLETE DE CHOLESKY (VERSION STABILISEE)

C ENTREES:
C --------
C NTDL   : NOMBRE D'INCONNUES DU SYSTEME
C NDSM   : NOMBRE DE SECONDS MEMBRES ET SOLUTIONS U
C NBDIR  : NOMBRE DE DIRECTIONS CONSERVEES POUR STABILISER

C LPLIGN : POINTEUR SUR LE COEFFICIENT DIAGONAL DE CHAQUE LIGNE DE AG
C LPCOLO : NO COLONNE DE CHAQUE COEFFICIENT STOCKE DE AG
C AG     : MATRICE INITIALE NON FACTORISEE

C B      : LE TABLEAU DES NDSM SECONDS MEMBRES

C LPLIGC : POINTEUR SUR LE COEFFICIENT DIAGONAL DE CHAQUE LIGNE DE AGC
C LPCOLC : NO COLONNE DE CHAQUE COEFFICIENT STOCKE DE AGC
C AGC    : MATRICE FACTORISEE INCOMPLETE <=> MATRICE DE PRECONDITIONNEMENT
C X0     : VECTEUR INITIAL (ITERATION 0)
C X,R,Z  : TABLEAUX AUXILIAIRES
C DIR    : LES NBDIR DIRECTIONS
C ADIR   : LES PRODUITS (A) * DIRECTIONS
C DAD    : LES PRODUITS SCALAIRES ( DIR , ADIR )
C BETA   : LES PRODUITS SCALAIRES ( R   , ADIR )

C SORTIES:
C --------
C U      : LE TABLEAU DES NDSM SOLUTIONS
C IERR   : CODE D'ERREUR, 0 SI PAS D'ERREUR, 1 SI NON CONVERGENCE DU GC
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY      ANALYSE NUMERIQUE UPMC PARIS    OCTOBRE 1989
C MODIFS : ALAIN PERRONNET  ANALYSE NUMERIQUE UPMC PARIS      AVRIL 1999
C MODIFS : ALAIN PERRONNET  LJLL UPMC & St PIERRE du PERRAY JANVIER 2009
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      DOUBLE PRECISION  AG,B,AGC,X0,X,R,Z,U,
     &                  DIR,ADIR,DAD,BETA,
     &                  PROSCD,SQRT,RNORME,RNORM,RESIDU,
     &                  ALPHA,DKADK,RKADK ,BETAK,EPS0,EPS1
      DIMENSION         LPLIGN(NTDL+1),LPCOLO(1),AG(1),B(NTDL,NDSM),
     &                  LPLIGC(NTDL+1),LPCOLC(1),AGC(1),X0(NTDL),
     &                  X(NTDL),R(NTDL),Z(NTDL),U(NTDL,NDSM),
     &                  DIR(NTDL,NBDIR),ADIR(NTDL,NBDIR),
     &                  DAD(NBDIR),BETA(NBDIR)
      COMMON / UNITES / LECTEU,IMPRIM,INTERA,NUNITE(29)
      COMMON / EPSSSS / EPZERO,EPSXYZ
      COMMON / MSIMTA / NOIMPR

10010 FORMAT(' GC: Iteration',I7,' Eps0=',1PG15.6,' Residu0=',1PG15.6,
     %' rnorme=',1PG15.6)
20010 FORMAT(' CG: Iteration',I7,' Eps0=',1PG15.6,' Residu0=',1PG15.6,
     %' rnorme=',1PG15.6)

 1000 FORMAT(' GC: Iteration',I7,' Eps1=',1PG15.6,' Residu =',1PG15.6,
     %' rnorme=',1PG15.6,'  => CONVERGENCE DU GC')
 2000 FORMAT(' CG: Iteration',I7,' Eps1=',1PG15.6,' Residu =',1PG15.6,
     %' rnorme=',1PG15.6,'  => GC CONVERGENCE')

 1100 FORMAT(/' ERREUR GCPRCH: NON CONVERGENCE apres',I10,
     &       ' ITERATIONS'/)
 2100 FORMAT(/' ERROR GCPRCH: NO CONVERGENCE after',I10,
     &       ' ITERATIONS'/)
 1300 FORMAT(' GC: NORME du RESIDU',1PG13.6)
 2300 FORMAT(' CG: NORM of RESIDUE',1PG13.6)

C     INITIALISATION
C     --------------
      ITEMAX = NTDL
      NCODSA = 1
      IERR   = 0
      EPS0   = EPZERO
      EPS1   = EPZERO

C     BOUCLE SUR LES NDSM SYSTEMES A RESOUDRE
C     ---------------------------------------
      DO 100 J=1,NDSM

C        K LE COMPTEUR D'ITERATIONS
         K = 0

C        (R) = (B) - (( A )) * (X0)
C        ---------------------------
         CALL MAGCVE( 0, 1D0, NTDL,
     %                NCODSA, LPLIGN, LPCOLO, AG,  X0,
     %                ADIR(1,1) )
         DO I=1,NTDL
            R(I) = B(I,J) - ADIR(I,1)
            X(I) = X0(I)
         ENDDO
         RESIDU = PROSCD(R,R,NTDL)

         CALL DRCHGC( NTDL, LPLIGC, LPCOLC, AGC, R, DIR(1,1) )
         RNORME = PROSCD( R, DIR(1,1), NTDL )
         IF( NOIMPR .GE. 3 ) THEN
            IF( LANGAG .EQ. 0 ) THEN
               WRITE(IMPRIM,10010) K, EPS0, RESIDU, RNORME
            ELSE
               WRITE(IMPRIM,20010) K, EPS0, RESIDU, RNORME
            ENDIF
         ENDIF
         IF( RESIDU .LE. 1D-30 ) GOTO 200

C        ITERATIONS DU GRADIENT CONJUGUE
C        -------------------------------
C        NUMDIR : LE NUMERO DE LA DIRECTION ACTUELLE
         NUMDIR = 1
C        NBDIC  : LE NOMBRE DE DIRECTIONS DE CALCUL
         NBDIC  = 1
C
 10      K = K + 1

C        (AD) = (( A )) * (D)
C        --------------------
         CALL MAGCVE( 0, 1D0, NTDL,
     %                NCODSA, LPLIGN, LPCOLO, AG,
     %                DIR(1,NUMDIR),  ADIR(1,NUMDIR) )

C        ALPHA = ( R , R ) / ( D , (A) * D )
C        -------------------------------------
         DKADK = PROSCD(DIR(1,NUMDIR),ADIR(1,NUMDIR),NTDL)
         ALPHA = RNORME / DKADK
         DAD(NUMDIR) = DKADK

C        (X) = (X) + ALPHA * (D)
C        (R) = (R) - ALPHA * ((A)) * (D)
C        ---------------------------------
         DO I=1,NTDL
            X(I) = X(I) + ALPHA *  DIR(I,NUMDIR)
            R(I) = R(I) - ALPHA * ADIR(I,NUMDIR)
         ENDDO

C        RNORME = ( R , R )
C        ------------------
         CALL DRCHGC(NTDL,LPLIGC,LPCOLC,AGC,R,Z)
         RNORM = RNORME
         RNORME = PROSCD(R,Z,NTDL)
         BETAK  = RNORME / RNORM

C        TEST D'ARRET DES ITERATIONS K AU DELA DES 4 PREMIERES
C        -----------------------------------------------------
         RESIDU = PROSCD( R, R, NTDL )

C        SEUIL POUR TEST DES ITERATIONS
         IF( K .EQ. 1 ) THEN
            IF( RESIDU .GT. 1D-30 ) THEN
               EPS1 = EPS0 * RESIDU
            ELSE
               EPS1 = EPS0
            ENDIF
         ENDIF

         IF( K .GE. 3 .AND. RESIDU .LE. EPS1 ) GOTO 200

C        BETA = ( R , ADIR ) / ( D , (A) D )
C        -----------------------------------
C        + RE-ORTHOGONALISATION EVENTUELLE
C        ---------------------------------
         DO NDI=1,NBDIC
            RKADK  = PROSCD( Z, ADIR(1,NDI), NTDL )
            BETA(NDI) = - RKADK / DAD(NDI)
         ENDDO

C        SUIVI DE LA CONVERGENCE
C        -----------------------
C        ( REDEMARRAGE )
C        ---------------
         IF( ALPHA .LT. EPS0 ) THEN
            DO NDI=1,NBDIC
               BETA(NDI) = 0.D0
            ENDDO
         END IF

C        (D) = (R) + BETAK * (D)
C        -----------------------
         DO NDI=1,NBDIC
            DO I=1,NTDL
               Z(I) = Z(I) + BETA(NDI) * DIR(I,NDI)
            ENDDO
         ENDDO
         NUMDIR = NUMDIR + 1
         NBDIC  = NBDIC  + 1
         IF (NUMDIR.GT.NBDIR) NUMDIR = 1
         IF (NBDIC .GT.NBDIR) NBDIC  = NBDIR
         DO I=1,NTDL
            DIR(I,NUMDIR) = Z(I)
         ENDDO

C-----------------------------------------------------------
C        VERIFICATION DE L'ALGORITHME
C        ----------------------------
         IF( NOIMPR .GE. 5 ) THEN
            WRITE (IMPRIM,1400) K,RESIDU,RNORME,BETAK
            WRITE (IMPRIM,1401) NUMDIR,NBDIC
            WRITE (IMPRIM,1402)
     &         (NDI,BETA(NDI),NDI=1,NBDIC-1)
            DO NDI=1,NBDIC - 1
               RKADK  = PROSCD(DIR(1,NUMDIR),ADIR(1,NDI),NTDL)
               BETA(NDI) =  RKADK
            ENDDO
            WRITE (IMPRIM,1403)
     &         (NDI,BETA(NDI),NDI=1,NBDIC-1)
         ENDIF
 1400 FORMAT(/,1X,80(1H-)/,' GC: VERIFICATION DE L''ALGORITHME',
     &       /,' ITERATION         ',I10,
     &       /,' NORME DU RESIDU   ',1PG13.6,
     &       /,' RNORME            ',1PG13.6,
     &       /,' BETA CLASSIQUE    ',1PG13.6)
 1401 FORMAT(' ORTHOGONALITE',
     &       /,' DIRECTION ACTUELLE',I10,
     &       /,' NB DE DIRECTIONS  ',I10)
 1402 FORMAT(' DIRECTION :',I10,' BETA  ',1PG13.6)
 1403 FORMAT(' DIRECTION :',I10,' (D,AD)',1PG13.6)
C-----------------------------------------------------------

C        TEST D'ARRET DU PROGRAMME
C        -------------------------
         IF ( K .GE. ITEMAX ) THEN
            IF( LANGAG .EQ. 0 ) THEN
               WRITE (IMPRIM,1100) ITEMAX
            ELSE
               WRITE (IMPRIM,2100) ITEMAX
            ENDIF
            IERR = 1
            GOTO 210
         ELSE
            GO TO 10
         END IF

C        FIN DE LA BOUCLE SUR LES ITERATIONS DU GRADIENT CONJUGUE
C        --------------------------------------------------------
 200     IF( NOIMPR .GE. 3 ) THEN
            IF( LANGAG .EQ. 0 ) THEN
               WRITE (IMPRIM,1000) K, EPS1, RESIDU, RNORME
            ELSE
               WRITE (IMPRIM,2000) K, EPS1, RESIDU, RNORME
            ENDIF
         ENDIF

C        (R) = (B) - ((A)) * (X)
C        -----------------------
 210     CALL MAGCVE( 0, 1D0, NTDL,
     %                NCODSA, LPLIGN, LPCOLO, AG,  X,  Z )
         DO I=1,NTDL
            R(I)   = B(I,J) - Z(I)
            U(I,J) = X(I)
         ENDDO

         RESIDU = PROSCD(R,R,NTDL)
         RESIDU = SQRT(RESIDU)
         IF ( NOIMPR .GE. 2 ) THEN
            IF( LANGAG .EQ. 0 ) THEN
               WRITE (IMPRIM,1300) RESIDU
            ELSE
               WRITE (IMPRIM,2300) RESIDU
            ENDIF
         ENDIF

 100  ENDDO

      RETURN
      END
