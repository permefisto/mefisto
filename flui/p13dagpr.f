      SUBROUTINE P13DAGPR( COEF,   NBSOM,  XYZSOM,
     %                     NBNOEF, NBEF,   NONOEF, NONOSO,
     %                     NBCOPG, LPDIAG, PG )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT: CONSTRUIRE ET ASSEMBLER DANS LA MATRICE PROFIL GLOBALE LES
C ---- MATRICES ELEMENTAIRES COEF (GRAD p, GRAD q) + epsilon (p,q)
C      POUR UNE INTERPOLATION P1 EN 3D DE FONCTIONS DE BASE COORDONNEES
C      BARYCENTRIQUES LAMBDA(1:4)
C
C ENTREES:
C --------
C COEF   : COEFFICIENT MULTIPLICATEUR DE LA MATRICE
C NBSOM  : NOMBRE DE SOMMETS DU MAILLAGE
C XYZSOM : 3 COORDONNEES DES NBSOM SOMMETS DU MAILLAGE
C NBNOEF : NOMBRE DE SOMMETS (4) D'UN EF
C NBEF   : NOMBRE D'EF DU MAILLAGE
C NONOEF : NUMERO DES NBNOEF SOMMETS DES NBEF EF
C NONOSO : NONOSO( NOEUD ) = NUMERO DU SOMMET DE 1 A NBSOM
C
C NBCOPG : NOMBRE DE COEFFICIENTS DE LA MATRICE GLOBALE PG
C LPDIAG : POINTEUR SUR LES COEFFICIENTS DIAGONAUX   DE PG
C
C MODIFIES:
C ---------
C PG     : VALEUR DES COEFFICIENTS DE LA MATRICE GLOBALE ASSEMBLEE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY Novembre 2010
C23456---------------------------------------------------------------012
      IMPLICIT  NONE
      DOUBLE PRECISION EPSILON
      PARAMETER       (EPSILON=1D-7)
      include"./incl/langue.inc"
      REAL              XYZSOM(3,NBSOM)
      INTEGER           NBSOM, NBNOEF, NBEF, NONOEF(NBEF,NBNOEF),
     %                  NONOSO(1:*), LPDIAG(0:NBSOM), NBCOPG
      DOUBLE PRECISION  PG(1:NBCOPG), COEF
C
      DOUBLE PRECISION  DF(3,3), DFM1(3,3), DELTA, DFM1DLa(3,4), S, D,
     %                  X, Y, Z, EPSii, EPSij
      INTEGER           I, J, K, NEF, NS1, NS2, NS3, NS4, IEG, JEG, KG
      INTRINSIC         ABS
C
C     MISE A ZERO DE LA MATRICE PROFIL GLOBALE
      DO I = 1, NBCOPG
         PG( I ) = 0D0
      ENDDO
C
      DO 100 NEF = 1, NBEF
C
C        NUMERO DES 4 SOMMETS DU TETRAEDRE NEF
         NS1 = NONOEF(NEF,1)
         NS2 = NONOEF(NEF,2)
         NS3 = NONOEF(NEF,3)
         NS4 = NONOEF(NEF,4)
C
C        CONSTRUCTION DE LA MATRICE DF
         X = XYZSOM(1,NS1)
         DF(1,1) = XYZSOM(1,NS2) - X
         DF(2,1) = XYZSOM(1,NS3) - X
         DF(3,1) = XYZSOM(1,NS4) - X
C
         Y = XYZSOM(2,NS1)
         DF(1,2) = XYZSOM(2,NS2) - Y
         DF(2,2) = XYZSOM(2,NS3) - Y
         DF(3,2) = XYZSOM(2,NS4) - Y
C
         Z = XYZSOM(3,NS1)
         DF(1,3) = XYZSOM(3,NS2) - Z
         DF(2,3) = XYZSOM(3,NS3) - Z
         DF(3,3) = XYZSOM(3,NS4) - Z
C
C        VALEUR ABSOLUE DU DETERMINANT DE DF
         DELTA = ABS(DF(1,1) * (DF(2,2) * DF(3,3) - DF(3,2) * DF(2,3))
     %              +DF(2,1) * (DF(3,2) * DF(1,3) - DF(1,2) * DF(3,3))
     %              +DF(3,1) * (DF(1,2) * DF(2,3) - DF(2,2) * DF(1,3)) )
C
         IF( DELTA .LE. 0D0 ) THEN
            IF( LANGAG .EQ. 0 ) THEN
               PRINT*,'p13dagpr: ATTENTION EF',NEF,
     %                ' de VOLUME*6=',DELTA,' NON PRIS EN COMPTE'
            ELSE
               PRINT*,'p13dagpr: ATTENTION FE',NEF,
     %                ' of VOLUME*6=',DELTA,' is NOT COMPUTED'
            ENDIF
            GOTO 100
         ENDIF
C
C        LE TETRAEDRE EST SUPPOSE DE VOLUME NON NUL
C        LES 9 COEFFICIENTS DE LA MATRICE INVERSE DFM1 SANS / DELTA
         DFM1(1,1) = DF(2,2) * DF(3,3) - DF(3,2) * DF(2,3)
         DFM1(2,1) = DF(2,3) * DF(3,1) - DF(3,3) * DF(2,1)
         DFM1(3,1) = DF(2,1) * DF(3,2) - DF(3,1) * DF(2,2)
C
         DFM1(1,2) = DF(1,3) * DF(3,2) - DF(1,2) * DF(3,3)
         DFM1(2,2) = DF(1,1) * DF(3,3) - DF(1,3) * DF(3,1)
         DFM1(3,2) = DF(1,2) * DF(3,1) - DF(1,1) * DF(3,2)
C
         DFM1(1,3) = DF(1,2) * DF(2,3) - DF(2,2) * DF(1,3)
         DFM1(2,3) = DF(1,3) * DF(2,1) - DF(2,3) * DF(1,1)
         DFM1(3,3) = DF(1,1) * DF(2,2) - DF(2,1) * DF(1,2)
C
C        [DFM1] [DLa]
         DFM1DLa(1,1) = -DFM1(1,1) - DFM1(1,2) - DFM1(1,3)
         DFM1DLa(1,2) =  DFM1(1,1)
         DFM1DLa(1,3) =  DFM1(1,2)
         DFM1DLa(1,4) =  DFM1(1,3)
C
         DFM1DLa(2,1) = -DFM1(2,1) - DFM1(2,2) - DFM1(2,3)
         DFM1DLa(2,2) =  DFM1(2,1)
         DFM1DLa(2,3) =  DFM1(2,2)
         DFM1DLa(2,4) =  DFM1(2,3)
C
         DFM1DLa(3,1) = -DFM1(3,1) - DFM1(3,2) - DFM1(3,3)
         DFM1DLa(3,2) =  DFM1(3,1)
         DFM1DLa(3,3) =  DFM1(3,2)
         DFM1DLa(3,4) =  DFM1(3,3)
C
C        COEFFICIENT MULTIPLICATEUR
         D = COEF / ( 6D0 * DELTA )

C        COEF EPSILON DELTA Integrale P1i P1j dX
         EPSij = COEF * DELTA * EPSILON / 120D0
         EPSii = EPSij * 2
C
         IF( NBNOEF .LE. 4 ) THEN
C
C           TETRAEDRE BREZZI-FORTIN (NOEUDS=SOMMETS)
            DO I = 1, 4
C
C              ASSEMBLAGE DE LA MATRICE ELEMENTAIRE DANS LA MATRICE PROFIL GLOBA
C              LE NO GLOBAL DU I-EME DEGRE DE LIBERTE LOCAL DE L'EF NEF
               IEG = NONOEF(NEF,I)
C
C              Integrale t([DFM1] [DLa]) ([DFM1] [DLa]) sur l'EF NEF
               DO J = 1, I
C
                  S = 0D0
                  DO K = 1, 3
                     S = S + DFM1DLa(K,I) * DFM1DLa(K,J)
                  ENDDO
C
C                 LE NO GLOBAL DU I-EME DEGRE DE LIBERTE LOCAL
                  JEG = NONOEF(NEF,J)
C
C                 ADRESSE DANS LA MATRICE GLOBALE PROFIL DE AE(I,J)
                  IF( JEG .LE. IEG ) THEN
                     KG = LPDIAG(IEG) - IEG + JEG
                  ELSE
                     KG = LPDIAG(JEG) - JEG + IEG
                  ENDIF
C
C                 COEFFICIENT ELEMENTAIRE SOMME AU COEFFICIENT GLOBAL KG
                  IF( I .EQ. J ) THEN
                     S = S * D + EPSii
                  ELSE
                     S = S * D + EPSij
                  ENDIF
                  PG( KG ) = PG( KG ) + S
C
               ENDDO
            ENDDO
C
         ELSE
C
C           TETRAEDRE TAYLOR-HOOD (NOEUDS=SOMMETS+MILIEUX des ARETES)
            DO I = 1, 4
C
C              ASSEMBLAGE DE LA MATRICE ELEMENTAIRE DANS LA MATRICE PROFIL GLOBA
C              LE NO GLOBAL DU I-EME DEGRE DE LIBERTE LOCAL DE L'EF NEF
               IEG = NONOSO( NONOEF(NEF,I) )
C
C              Integrale t([DFM1] [DLa]) ([DFM1] [DLa]) sur l'EF NEF
               DO J = 1, I
C
                  S = 0D0
                  DO K = 1, 3
                     S = S + DFM1DLa(K,I) * DFM1DLa(K,J)
                  ENDDO
C
C                 LE NO GLOBAL DU I-EME DEGRE DE LIBERTE LOCAL
                  JEG = NONOSO( NONOEF(NEF,J) )
C
C                 ADRESSE DANS LA MATRICE GLOBALE PROFIL DE AE(I,J)
                  IF( JEG .LE. IEG ) THEN
                     KG = LPDIAG(IEG) - IEG + JEG
                  ELSE
                     KG = LPDIAG(JEG) - JEG + IEG
                  ENDIF
C
C                 COEFFICIENT ELEMENTAIRE SOMME AU COEFFICIENT GLOBAL KG
                  IF( I .EQ. J ) THEN
                     S = S * D + EPSii
                  ELSE
                     S = S * D + EPSij
                  ENDIF
                  PG( KG ) = PG( KG ) + S
C
               ENDDO
            ENDDO
C
         ENDIF
C
 100  CONTINUE

c     call affvect( 'p13dagpr: PG=', 20, PG )
c     call afl1ve(  'p13dagpr: PG=', NBCOPG, PG )

      RETURN
      END
