      SUBROUTINE INTSOLFA( NOSURF, NBNOFA, NONOFA, NBCOOR, XYZNOE,
     %                     NCAS,   solution,
     %                     NUMISU, NUMASU, INTSOL )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER L'INTEGRALE DE LA SOLUTION SUR UNE FACE D'UN EF 3D
C -----    APPARTENANT A LA SURFACE NOSURF DE R3

C ENTREES:
C --------
C NOSURF : NUMERO DE LA SURFACE (de NUMISU a NUMASU)
C NBNOFA : NOMBRE DE  NOEUDS DE LA FACE
C NONOFA : NUMERO DES NOEUDS DE LA FACE
C NBCOOR : NOMBRE DE COORDONNEES DES NOEUDS (3 ou 6)
C XYZNOE : COORDONNEES DES NOEUDS DU MAILLAGE DE L'OBJET
C NCAS   : NUMERO DU VECTEUR SOLUTION d'INTEGRALE A CALCULER
C solution: TABLEAU de pointeurs sur le tableau solution
C NUMISU : NUMERO MINIMAL DES NUMEROS DE SURFACE DE L'OBJET
C NUMASU : NUMERO MAXIMAL DES NUMEROS DE SURFACE DE L'OBJET

C MODIFIE:
C --------
C INTSOL : INTEGRALE DU CAS NCAS DE LA SOLUTION SUR LES FACES DES SURFACES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY SEPTEMBRE 2010
C2345X7..............................................................012
      DOUBLE PRECISION  INTSOL(NUMISU:NUMASU)
      INTEGER           NONOFA(NBNOFA)
      REAL              XYZNOE(NBCOOR,*)
      DOUBLE PRECISION  DELTAG, DGL(2,3)

      type typ_dptab
         DOUBLE PRECISION, dimension(:), pointer :: dptab
      end type typ_dptab
      type( typ_dptab ), dimension(:), allocatable :: solution

      IF( NBNOFA .NE. 3 .AND. NBNOFA .NE. 6 ) THEN
         print *,'INTSOLFA: NBNOFA=',NBNOFA,
     %   ' NON 3 (Interpolation P1) et NON 6 (P2) => NON PROGRAMME'
         RETURN
      ENDIF

C     TRIANGLE P1 DE NUMEROS DE SOMMETS
      N1 = NONOFA(1)
      N2 = NONOFA(2)
      N3 = NONOFA(3)

C     RECHERCHE DU JACOBIEN DE G: TRIANGLE REFERENCE --> FACE
      DGL(1,1) = XYZNOE( 1, N2 ) - XYZNOE( 1, N1 )
      DGL(1,2) = XYZNOE( 2, N2 ) - XYZNOE( 2, N1 )
      DGL(1,3) = XYZNOE( 3, N2 ) - XYZNOE( 3, N1 )

      DGL(2,1) = XYZNOE( 1, N3 ) - XYZNOE( 1, N1 )
      DGL(2,2) = XYZNOE( 2, N3 ) - XYZNOE( 2, N1 )
      DGL(2,3) = XYZNOE( 3, N3 ) - XYZNOE( 3, N1 )
      CALL JAR2R3( DGL, DELTAG )

C     L'INTEGRALE DE LA SOLUTION SUR LE TRIANGLE
      IF( NBNOFA .EQ. 6 ) THEN

C        TRIANGLE P2 DE NUMEROS DE MILIEUX DES 3 COTES
         N1 = NONOFA(1)
         N2 = NONOFA(2)
         N3 = NONOFA(3)

      ENDIF

C     INTEGRATION NUMERIQUE EXACTE P1 et P2
      DELTAG = DELTAG / 6D0 * ( solution( NCAS )%dptab( N1 )
     %                        + solution( NCAS )%dptab( N2 )
     %                        + solution( NCAS )%dptab( N3 ) )

C     L'INTEGRALE DE LA SOLUTION SUR LA SURFACE NOSURF
      INTSOL( NOSURF ) = INTSOL( NOSURF ) + DELTAG

      RETURN
      END
