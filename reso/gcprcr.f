      SUBROUTINE GCPRCR( NTDL,   NDSM,   NBDIR, NTDLFX,
     %                   LPLIGN, LPCOLO, AG, B,
     %                   LPLIGC, LPCOLC, AGC,
     %                   X0, X, R, Z, DIR, ADIR, DAD, BETA,
     %                   U, IERR )
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    METHODE DE GRADIENT CONJUGUE AVEC PRECONDITIONNEMENT PAR
C -----    FACTORISATION INCOMPLETE DE CROUT (VERSION STABILISEE)
C
C ENTREES:
C --------
C NTDL   : NOMBRE D'INCONNUES DU SYSTEME
C NDSM   : NOMBRE DE SECONDS MEMBRES
C NBDIR  : NOMBRE DE DIRECTIONS CONSERVEES POUR STABILISER
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
C X,R,Z  : TABLEAUX AUXILIAIRES
C DIR    : LES NBDIR DIRECTIONS
C ADIR   : LES PRODUITS (A) * DIRECTIONS
C DAD    : LES PRODUITS SCALAIRES ( DIR , ADIR )
C BETA   : LES PRODUITS SCALAIRES ( R   , ADIR )
C
C SORTIES:
C --------
C U      : LE TABLEAU DES NDSM SOLUTIONS
C IERR   : CODE D'ERREUR, 0 SI PAS D'ERREUR, 1 SI NON CONVERGENCE DU GC
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PASCAL JOLY     ANALYSE NUMERIQUE UPMC PARIS     OCTOBRE 1989
C MODIFS : ALAIN PERRONNET ANALYSE NUMERIQUE UPMC PARIS       AVRIL 1999
C MODIFS : ALAIN PERRONNET Laboratoire J-L.LIONS UPMC Paris OCTOBRE 2007
C MODIFS : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY     MARS 2009
C AJOUTS : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY     MARS 2014
C23456---------------------------------------------------------------012
      include"./incl/langue.inc"
      COMMON / UNITES / LECTEU, IMPRIM, INTERA, NUNITE(29)
      COMMON / EPSSSS / EPZERO, EPSXYZ
      COMMON / MSIMTA / NOIMPR
      DOUBLE PRECISION  AG, B, AGC, X0, X, R, Z, U,
     %                  DIR, ADIR, DAD, BETA,
     %                  PROSCD, RNORME, RNORM, RESIDU, RESIDU0,
     %                  ALPHA, DKADK, RKADK , BETAK, EPS0, EPS1
      INTEGER           NTDLFX(NTDL)
      DIMENSION         LPLIGN(NTDL+1), LPCOLO(*), AG(*), B(NTDL,NDSM),
     %                  LPLIGC(NTDL+1), LPCOLC(*), AGC(*), X0(NTDL),
     %                  X(NTDL), R(NTDL), Z(NTDL), U(NTDL,NDSM),
     %                  DIR(NTDL,NBDIR), ADIR(NTDL,NBDIR),
     %                  DAD(NBDIR), BETA(NBDIR)
C
10010 FORMAT(' gcprcr: Iteration',I7,' Eps0=',G15.6,' Residu0=',G15.6,
     %' rnorme=',G15.6)
C
 1000 FORMAT(' gcprcr: Iteration',I7,' Eps1=',G15.6,' Residu =',G15.6,
     %' rnorme=',G15.6,'  => CONVERGENCE DU GC')
 2000 FORMAT(' gcprcr: Iteration',I7,' Eps1=',G15.6,' Residue=',G15.6,
     %' rnorm =',G15.6,'  => CG CONVERGENCE')
C
 1100 FORMAT(/' ERREUR GCPRCR: NON CONVERGENCE apres',I10,
     %       ' ITERATIONS'/)
 2100 FORMAT(/' ERROR GCPRCR: NO CONVERGENCE after',I10,
     %       ' ITERATIONS'/)
 1300 FORMAT(' ITERATION',I7,' DU GC  NORME du RESIDU=',G13.6)
 2300 FORMAT(' ITERATION',I7,' of CG  NORM of RESIDUE=',G13.6)
C
C     INITIALISATION
C     --------------
      ITEMAX = NTDL
      NCODSA = 1
      IERR   = 0
      EPS0   = EPZERO
      EPS1   = EPZERO
      IF( NTDL .LT. 10000 ) THEN
         KDIVGC = 20
      ELSE
         KDIVGC = 100
      ENDIF
C
C     BOUCLE SUR LES NDSM SYSTEMES A RESOUDRE
C     ---------------------------------------
      DO 100 J=1,NDSM
C
C        K COMPTEUR DES ITERATIONS DU GC
         K = 0
         RESIDU0 = 1D100
C
C        (R) = (B) - (( A )) * (X0)
C        ---------------------------
         CALL MAGCVX( NTDL, NTDLFX, LPLIGN, LPCOLO, AG,  X0,
     %                ADIR(1,1) )
C
C        R VECTEUR RESIDU INITIAL
C        ------------------------
         RESIDU = 0D0
         DO I=1,NTDL
            R(I)   = B(I,J) - ADIR(I,1)
            RESIDU = RESIDU + R(I) ** 2
            X(I)   = X0(I)
         ENDDO
         CALL DRCRGC( NTDL, LPLIGC, LPCOLC, AGC, R, DIR(1,1) )
         RNORME = PROSCD( R, DIR(1,1), NTDL )
         WRITE(IMPRIM,10010) K, EPS0, RESIDU, RNORME
         IF( RESIDU .LE. 1D-30 ) GOTO 200
C
C        ITERATIONS DU GRADIENT CONJUGUE
C        -------------------------------
C        NUMDIR: LE NUMERO DE LA DIRECTION ACTUELLE
         NUMDIR = 1
C        NBDIC : LE NOMBRE DE DIRECTIONS DE CALCUL
         NBDIC  = 1
C
 10      K = K + 1
C
C        (AD) = (( A )) * (D)
C        --------------------
         CALL MAGCVX( NTDL, NTDLFX, LPLIGN, LPCOLO, AG,
     %                DIR(1,NUMDIR),  ADIR(1,NUMDIR) )
C
C        ALPHA = ( R , R ) / ( D , (A) * D )
C        -------------------------------------
         DKADK = PROSCD( DIR(1,NUMDIR), ADIR(1,NUMDIR), NTDL )
         ALPHA = RNORME / DKADK
         DAD(NUMDIR) = DKADK
C
C        (X) = (X) + ALPHA * (D)
C        (R) = (R) - ALPHA * ((A)) * (D)
C        ---------------------------------
         DO I=1,NTDL
            X(I) = X(I) + ALPHA *  DIR(I,NUMDIR)
            R(I) = R(I) - ALPHA * ADIR(I,NUMDIR)
         ENDDO
C
C        RNORME = ( R , R )
C        ------------------
         CALL DRCRGC( NTDL, LPLIGC, LPCOLC, AGC, R, Z )
         RNORM  = RNORME
         RNORME = PROSCD( R, Z, NTDL )
         BETAK  = RNORME / RNORM
C
C        TEST D'ARRET DES ITERATIONS DU GC
C        ---------------------------------
         RESIDU0 = RESIDU
         RESIDU  = PROSCD( R, R, NTDL )

         IF( K .GT. KDIVGC .AND. RESIDU .GT. 10D0 * RESIDU0 ) THEN
C           DETECTION D'UNE DIVERGENCE DU GC
            IF( LANGAG .EQ. 0 ) THEN
               WRITE (IMPRIM,19000) K, RESIDU0, RESIDU
            ELSE
               WRITE (IMPRIM,29000) K, RESIDU0, RESIDU
            ENDIF
19000       FORMAT(' gcprcr: Iteration',I7,' Residu0=',G15.6,
     %             ' Residu=',G15.6,' => DIVERGENCE du GC'/
     %             ' UTILISER une METHODE DIRECTE pour resoudre Ax=b' )
29000       FORMAT(' gcprcr: Iteration',I7,' Residue0=',G15.6,
     %             ' Residue=',G15.6,' => DIVERGENCE of CG'/
     %             ' USE a DIRECT METHOD to solve Ax=b')
            IERR = 1
            GOTO 210
         ENDIF
C
C        SEUIL POUR TEST DES ITERATIONS
         IF( K .EQ. 1 ) THEN
            IF( RESIDU .GT. 1D-30 ) THEN
               EPS1 = EPS0 * RESIDU
            ELSE
               EPS1 = EPS0
            ENDIF
         ENDIF
C
         IF( K .GE. 3 .AND. RESIDU .LE. EPS1 ) GOTO 200
C
C        BETA = ( R , ADIR ) / ( D , (A) D )
C        -----------------------------------
C        + RE-ORTHOGONALISATION EVENTUELLE
C        ---------------------------------
         DO NDI=1,NBDIC
            RKADK  = PROSCD( Z, ADIR(1,NDI), NTDL )
            BETA(NDI) = - RKADK / DAD(NDI)
         ENDDO
C
C        SUIVI DE LA CONVERGENCE
C        -----------------------
C        ( REDEMARRAGE )
C        ---------------
         IF( ALPHA .LT. EPS0 ) THEN
            DO NDI=1,NBDIC
               BETA(NDI) = 0.D0
            ENDDO
         ENDIF
C
C        (D) = (R) + BETAK * (D)
C        -----------------------
         DO NDI=1,NBDIC
            DO I=1,NTDL
               Z(I) = Z(I) + BETA(NDI) * DIR(I,NDI)
            ENDDO
         ENDDO
C
C        AJOUT D'UNE DIRECTION POUR STABILISER
         NUMDIR = NUMDIR + 1
         NBDIC  = NBDIC  + 1
         IF ( NUMDIR .GT. NBDIR ) NUMDIR = 1
         IF ( NBDIC  .GT. NBDIR ) NBDIC  = NBDIR
         DO I=1,NTDL
            DIR(I,NUMDIR) = Z(I)
         ENDDO
C
C-----------------------------------------------------------
C        VERIFICATION DE L'ALGORITHME
C        ----------------------------
         IF( NOIMPR .GE. 5 ) THEN
            WRITE (IMPRIM,1400) K, RESIDU, RNORME, BETAK
            WRITE (IMPRIM,1401) NUMDIR, NBDIC
            WRITE (IMPRIM,1402)
     %         (NDI,BETA(NDI),NDI=1,NBDIC-1)
            DO NDI=1,NBDIC - 1
               RKADK  = PROSCD( DIR(1,NUMDIR), ADIR(1,NDI), NTDL )
               BETA(NDI) =  RKADK
            ENDDO
            WRITE (IMPRIM,1403)
     %         (NDI,BETA(NDI),NDI=1,NBDIC-1)
         ENDIF
 1400 FORMAT(/,1X,80(1H-)/,' VERIFICATION DE L''ALGORITHME',
     %       /,' ITERATION             ',I10,
     %       /,' NORME DU RESIDU (R,R) ',G13.6,
     %       /,' NORME DU RESIDU (R,Z) ',G13.6,
     %       /,' BETA CLASSIQUE        ',G13.6)
 1401 FORMAT(' ORTHOGONALITE',
     %       /,' DIRECTION ACTUELLE',I10,
     %       /,' NB DE DIRECTIONS  ',I10)
 1402 FORMAT(' DIRECTION :',I10,' BETA  ',G13.6)
 1403 FORMAT(' DIRECTION :',I10,' (D,AD)',G13.6)
C-----------------------------------------------------------
C
C        TEST D'ARRET DES ITERATIONS DU GC
C        ---------------------------------
         IF ( K .GE. ITEMAX ) THEN
C
C           MAX ATTEINT DES ITERATIONS
            IF( LANGAG .EQ. 0 ) THEN
               WRITE (IMPRIM,1100) ITEMAX
            ELSE
               WRITE (IMPRIM,2100) ITEMAX
            ENDIF
            IERR = 1
            GOTO 210
C
         ELSE
C
C           PASSAGE A L'ITERATION SUIVANTE
            GOTO 10
C
         ENDIF
C
C        FIN CONVERGEE DES ITERATIONS DU GRADIENT CONJUGUE
C        -------------------------------------------------
 200     IF( LANGAG .EQ. 0 ) THEN
            WRITE (IMPRIM,1000) K, EPS1, RESIDU, RNORME
         ELSE
            WRITE (IMPRIM,2000) K, EPS1, RESIDU, RNORME
         ENDIF
C
C        (R) = (B) - ((A)) * (X)
C        -----------------------
 210     CALL MAGCVX( NTDL, NTDLFX, LPLIGN, LPCOLO, AG,  X,  Z )
         DO I=1,NTDL
            R(I)   = B(I,J) - Z(I)
            U(I,J) = X(I)
         ENDDO
         RESIDU = PROSCD( R, R, NTDL )
ccc         RESIDU = SQRT( RESIDU )

        IF( LANGAG .EQ. 0 ) THEN
            WRITE (IMPRIM,1300) K, RESIDU
         ELSE
            WRITE (IMPRIM,2300) K, RESIDU
         ENDIF

 100  CONTINUE
C
      RETURN
      END
