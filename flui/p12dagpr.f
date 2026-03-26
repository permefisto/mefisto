      SUBROUTINE P12DAGPR( COEF,   NBSOM,  XYZSOM,
     %                     NBNOEF, NBEF,   NONOEF, NONOSO,
     %                     NBCOPG, LPDIAG, PG )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT:     CONSTRUIRE ET ASSEMBLER DANS LA MATRICE PROFIL GLOBALE LES
C ----     MATRICES ELEMENTAIRES COEF (GRAD p, GRAD q) + epsilon (p,q)
C          POUR UNE INTERPOLATION P1 EN 2D AVEC LES 3 FONCTIONS DE BASE LAMBDAi
C ENTREES:
C --------
C COEF   : COEFFICIENT MULTIPLICATEUR DE LA MATRICE
C NBSOM  : NOMBRE DE SOMMETS DU MAILLAGE
C XYZSOM : 3 COORDONNEES DES NBSOM SOMMETS DU MAILLAGE
C NBNOEF : NOMBRE DE NOEUDS D'UN EF
C         (3 POUR BREZZI-FORTIN & 6 POUR TAYLOR-HOOD)
C NBEF   : NOMBRE D'EF DU MAILLAGE
C NONOEF : NUMERO DES NBNOEF NOEUDS DES NBEF EF
C NONOSO : NONOSO( NOEUD ) = NUMERO DU SOMMET DE 1 A NBSOM
C
C NBCOPG : NOMBRE DE COEFFICIENTS DE LA MATRICE GLOBALE PG
C LPDIAG : POINTEUR SUR LES COEFFICIENTS DIAGONAUX   DE PG
C
C MODIFIES:
C ---------
C PG     : VALEUR DES COEFFICIENTS DE LA MATRICE GLOBALE ASSEMBLEE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY Decembre 2010
C23456---------------------------------------------------------------012
      IMPLICIT  NONE
      DOUBLE PRECISION EPSILON
      PARAMETER       (EPSILON=1D-7)
      include"./incl/langue.inc"
      REAL              XYZSOM(3,*)
      INTEGER           NBSOM, NBNOEF, NBEF, NONOEF(NBEF,NBNOEF),
     %                  NONOSO(1:*), LPDIAG(0:NBSOM), NBCOPG
      DOUBLE PRECISION  PG(1:NBCOPG), COEF
C
      DOUBLE PRECISION  X21, Y21, X31, Y31, X32, Y32, DELTA, AE(6),
     %                  C, EPSii, EPSij
      INTEGER           I, J, NEF, NS1, NS2, NS3, IEG, JEG, KG, KE
      INTRINSIC         ABS

C     MISE A ZERO DE LA MATRICE PROFIL GLOBALE
      DO I = 1, NBCOPG
         PG( I ) = 0D0
      ENDDO

      DO 100 NEF = 1, NBEF

C        NUMERO DES 3 SOMMETS DU TRIANGLE NEF
         NS1 = NONOEF(NEF,1)
         NS2 = NONOEF(NEF,2)
         NS3 = NONOEF(NEF,3)

C        [DFM1] [DLa]
         X21 = XYZSOM(1,NS2) - XYZSOM(1,NS1)
         X31 = XYZSOM(1,NS3) - XYZSOM(1,NS1)
         X32 = XYZSOM(1,NS3) - XYZSOM(1,NS2)

         Y21 = XYZSOM(2,NS2) - XYZSOM(2,NS1)
         Y31 = XYZSOM(2,NS3) - XYZSOM(2,NS1)
         Y32 = XYZSOM(2,NS3) - XYZSOM(2,NS2)

C        JACOBIEN de DF
         DELTA = ABS( X21 * Y31 - X31 * Y21 )

         IF( DELTA .LE. 0D0 ) THEN
            IF( LANGAG .EQ. 0 ) THEN
               PRINT*, 'p12dagpr: ATTENTION EF',NEF,
     %                 ' de SURFACE*2=',DELTA,' NON PRIS EN COMPTE'
            ELSE
               PRINT*, 'p12dagpr: ATTENTION FE',NEF,
     %                 ' of SURFACE*2=',DELTA,' is NOT COMPUTED'
            ENDIF
            GOTO 100
         ENDIF

C        COEFFICIENT MULTIPLICATEUR / ( 2 * JACOBIEN de DF )
         C = COEF / ( 2D0 * DELTA )

C        COEF EPSILON DELTA Integrale P1i P1j dX
         EPSij = COEF * DELTA * EPSILON / 24D0
         EPSii = EPSij * 2

C        AE = Integrale t([DFM1] [DLa]) ([DFM1] [DLa]) dX sur l'EF NEF
C           + EPSILON Integrale tLa La dX DELTA
         AE(1) = (  Y32 * Y32 + X32 * X32 ) * C + EPSii
         AE(2) = (- Y32 * Y31 - X32 * X31 ) * C + EPSij
         AE(3) = (  Y31 * Y31 + X31 * X31 ) * C + EPSii
         AE(4) = (  Y32 * Y21 + X32 * X21 ) * C + EPSij
         AE(5) = (- Y31 * Y21 - X31 * X21 ) * C + EPSij
         AE(6) = (  Y21 * Y21 + X21 * X21 ) * C + EPSii

C        ASSEMBLAGE DE LA MATRICE ELEMENTAIRE DANS LA MATRICE PROFIL GLOBALE PG
ccc         print *,'p12dagpr: AE=',(AE(k),k=1,6)
         KE = 0

         IF( NBNOEF .LE. 3 ) THEN
C           TRIANGLE BREZZI-FORTIN (NOEUDS=SOMMETS)
C           XYZSOM CONTIENT LES XYZ DES SOMMETS P1
C           NONOEF SONT LES NUMEROS DES SOMMETS P1
            DO I = 1, 3
C
C              LE NO GLOBAL DU I-EME DEGRE DE LIBERTE LOCAL DE L'EF NEF
               IEG = NONOEF(NEF,I)
               DO J = 1, I

C                 LE NO GLOBAL DU I-EME DEGRE DE LIBERTE LOCAL
                  JEG = NONOEF(NEF,J)

C                 ADRESSE DANS LA MATRICE GLOBALE PROFIL DE AE(I,J)
                  IF( JEG .LE. IEG ) THEN
                     KG = LPDIAG(IEG) - IEG + JEG
                  ELSE
                     KG = LPDIAG(JEG) - JEG + IEG
                  ENDIF

C                 LE NUMERO DU COEFFICIENT DE LA MATRICE ELEMENTAIRE
                  KE = KE + 1

C                 COEFFICIENT ELEMENTAIRE SOMME AU COEFFICIENT GLOBAL KG
                  PG( KG ) = PG( KG ) + AE(KE)

               ENDDO
            ENDDO

         ELSE

C           TRIANGLE TAYLOR-HOOD (NOEUDS=SOMMETS+MILIEUX ARETES)
C           XYZSOM CONTIENT LES XYZ DES NOEUDS P2
C           NONOEF SONT LES NUMEROS DES NOEUDS P2
            DO I = 1, 3

C              LE NO GLOBAL DU I-EME DEGRE DE LIBERTE LOCAL DE L'EF NEF
               IEG = NONOSO( NONOEF(NEF,I) )

               DO J = 1, I

C                 LE NO GLOBAL DU I-EME DEGRE DE LIBERTE LOCAL
                  JEG = NONOSO( NONOEF(NEF,J) )

C                 ADRESSE DANS LA MATRICE GLOBALE PROFIL DE AE(I,J)
                  IF( JEG .LE. IEG ) THEN
                     KG = LPDIAG(IEG) - IEG + JEG
                  ELSE
                     KG = LPDIAG(JEG) - JEG + IEG
                  ENDIF

C                 LE NUMERO DU COEFFICIENT DE LA MATRICE ELEMENTAIRE
                  KE = KE + 1

C                 COEFFICIENT ELEMENTAIRE SOMME AU COEFFICIENT GLOBAL KG
                  PG( KG ) = PG( KG ) + AE(KE)

               ENDDO
            ENDDO

         ENDIF

ccc         print *,'p12dagpr: PG PR=',(PG(k),k=1,20)
 100  ENDDO

      call affvect( 'p12dagpr.f: PG=', 5,      PG )
      call afl1ve(  'p12dagpr.f: PG=', NBCOPG, PG )

      RETURN
      END
