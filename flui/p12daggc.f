      SUBROUTINE P12DAGGC( COEF,   NBSOM,  XYZSOM,
     %                     NBNOEF, NBELEM, NUNOEF, NONOSO,
     %                     NBCOPG, LPLIGN, LPCOLO, PG )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT:     CONSTRUIRE ET ASSEMBLER DANS LA MATRICE MORSE GLOBALE LES
C ----     MATRICES ELEMENTAIRES COEF (GRAD p, GRAD q) + epsilon (p,q)
C          POUR UNE INTERPOLATION P1 EN 2D DE FONCTIONS DE BASE LAMBDA(1:3)
C
C ENTREES:
C --------
C COEF   : COEFFICIENT MULTIPLICATEUR DE LA MATRICE
C NBSOM  : NOMBRE DE SOMMETS DU MAILLAGE
C XYZSOM : 3 COORDONNEES DES NOEUDS DU MAILLAGE
C NBNOEF : NOMBRE DE NOEUDS D'UN EF
C         (3 POUR BREZZI-FORTIN & 6 POUR TAYLOR-HOOD)
C NBELEM : NOMBRE DE TRIANGLES DU MAILLAGE
C NUNOEF : NUMERO DES NBNOEF SOMMETS DES NBELEM EF
C NONOSO : NONOSO( NOEUD ) = NUMERO DU SOMMET DE 1 A NBSOM
C
C NBCOPG : NOMBRE DE COEFFICIENTS DE LA MATRICE GLOBALE PG
C LPLIGN : POINTEUR SUR LES COEFFICIENTS DIAGONAUX      DE PG
C LPCOLO : NUMERO DES COLONNES DES COEFFICIENTS STOCKES DE PG
C
C SORTIE :
C --------
C PG     : VALEUR DES COEFFICIENTS DE LA MATRICE GLOBALE ASSEMBLEE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY Decembre 2010
C23456---------------------------------------------------------------012
      DOUBLE PRECISION EPSILON
      PARAMETER       (EPSILON=1D-7)
      include"./incl/langue.inc"
      REAL              XYZSOM(3,*)
      INTEGER           NUNOEF(NBELEM,NBNOEF),
     %                  NONOSO(1:*),
     %                  LPLIGN(0:NBSOM),
     %                  LPCOLO(1:NBCOPG)
      DOUBLE PRECISION  COEF, PG(1:NBCOPG)
C
      DOUBLE PRECISION  X21, Y21, X31, Y31, X32, Y32, DELTA, AE(6),
     %                  C, EPSii, EPSij
      INTEGER           I, NUELEM, NOSOTR(3)
      INTRINSIC         ABS

C     MISE A ZERO DE LA MATRICE MORSE GLOBALE
      DO I = 1, NBCOPG
         PG( I ) = 0D0
      ENDDO

      DO 100 NUELEM = 1, NBELEM

C        NUMERO XYZSOM DES 3 SOMMETS DU TRIANGLE NUELEM
C        XYZSOM CONTIENT LES XYZ DES SOMMETS P1 SI BREZZI-FORTIN
C        XYZSOM CONTIENT LES XYZ DES NOEUDS  P2 SI TAYLOR-YOUNG
         DO I=1,3
            NOSOTR(I) = NUNOEF(NUELEM,I)
         ENDDO

C        [DFM1] [DLa]
         X21 = XYZSOM(1,NOSOTR(2)) - XYZSOM(1,NOSOTR(1))
         X31 = XYZSOM(1,NOSOTR(3)) - XYZSOM(1,NOSOTR(1))
         X32 = XYZSOM(1,NOSOTR(3)) - XYZSOM(1,NOSOTR(2))

         Y21 = XYZSOM(2,NOSOTR(2)) - XYZSOM(2,NOSOTR(1))
         Y31 = XYZSOM(2,NOSOTR(3)) - XYZSOM(2,NOSOTR(1))
         Y32 = XYZSOM(2,NOSOTR(3)) - XYZSOM(2,NOSOTR(2))

C        JACOBIEN DE F
         DELTA = ABS( X21 * Y31 - X31 * Y21 )

         IF( DELTA .LE. 0D0 ) THEN
            IF( LANGAG .EQ. 0 ) THEN
               PRINT*, 'p12daggc: ATTENTION EF',NUELEM,
     %                 ' de SURFACE*2=',DELTA,' NON PRIS EN COMPTE'
            ELSE
               PRINT*, 'p12daggc: ATTENTION FE',NUELEM,
     %                 ' of SURFACE*2=',DELTA,' is NOT COMPUTED'
            ENDIF
            GOTO 100
         ENDIF

C        COEFFICIENT MULTIPLICATEUR / ( 2 * JACOBIEN de DF )
         C = COEF / ( 2D0 * DELTA )

C        COEF EPSILON DELTA Integrale P1i P1j dX
         EPSij = COEF * DELTA * EPSILON / 24D0
         EPSii = EPSij * 2

C        AE = Integrale t([DFM1] [DLa]) ([DFM1] [DLa]) dX sur l'EF NUELEM
C           + EPSILON Integrale tLa La dX DELTA
         AE(1) = (  Y32 * Y32 + X32 * X32 ) * C + EPSii
         AE(2) = (- Y32 * Y31 - X32 * X31 ) * C + EPSij
         AE(3) = (  Y31 * Y31 + X31 * X31 ) * C + EPSii
         AE(4) = (  Y32 * Y21 + X32 * X21 ) * C + EPSij
         AE(5) = (- Y31 * Y21 - X31 * X21 ) * C + EPSij
         AE(6) = (  Y21 * Y21 + X21 * X21 ) * C + EPSii

ccc         print *,'p12daggc: NOSOTR=',NOSOTR
ccc         print *,'p12daggc: AE=',(AE(k),k=1,6)

         IF( NBNOEF .GT. 3 ) THEN
C           EF TAYLOR-YOUNG
C           XYZSOM CONTIENT LES XYZ DES NOEUDS P2
C           NUNOEF SONT LES NUMEROS DES NOEUDS P2
C           NONOSO SONT LES NUMEROS DE SOMMET DES NOEUDS P2
C           MAIS LPLIGN, LPCOLO SONT DEFINIS SUR UN MAILLAGE P1
            DO I=1,3
               NOSOTR(I) = NONOSO( NUNOEF(NUELEM,I) )
            ENDDO
         ENDIF

C        ASSEMBLAGE DE LA MATRICE ELEMENTAIRE DANS LA MATRICE MORSE PG
         CALL ASMEGC( 3, NOSOTR, 1, AE, AE, 1, LPLIGN, LPCOLO, PG )

ccc         print *,'p12daggc: PG GC=',(PG(k),k=1,20)
 100  ENDDO

ccc      call affvect( 'p12daggc.f: PG=', 20,     PG )
ccc      call afl1ve(  'p12daggc.f: PG=', NBCOPG, PG )

      RETURN
      END
