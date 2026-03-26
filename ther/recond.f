      SUBROUTINE RECOND( NYOBJT, NUOBJT, NBCOOR, XYZPI, MNCOND,
     &                   COND  )
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++012
C BUT :    REMPLIR LE TABLEAU COND COEFFICIENTS DU TENSEUR
C -----    DE CONDUCTIVITE SOUS FORME D'UN VECTEUR
C          EN UN POINT D'INTEGRATION NUMERIQUE A L'INSTANT TEMPS
C
C          FILL the VECTOR COND with the SYMMETRIC CONDUCTIVITY
C          TENSOR COEFFICIENTS COMPUTED AT A QUADRATURE FORMULA POINT
C          AT a TIME

C ENTREES:
C --------
C NYOBJT : NO DU TYPE D'OBJET (1:POINT, 2:LIGNE, ... )
C NUOBJT : NO DE L'OBJET DANS SON LEXIQUE
C NBCOOR : NOMBRE (3 ou 6) DE COORDONNEES DES POINTS D'INTEGRATION
C          EN 2D FOURNIR NBCOOR=3 ET LA 3-EME VALEUR=0D0
C XYZPI  : LES 3 ou 6 COORDONNEES DU POINT D'INTEGRATION
C          EN 2D FOURNIR 3 COORDONNEES AVEC LA 3-EME VALEUR=0D0
C MNCOND : ADRESSE MCN DU TABLEAU 'CONDUCTIVITE'

C SORTIES:
C --------
C COND : LE TENSEUR SYMETRIQUE DE CONDUCTIVITE RANGE SOUS LA FORME
C        1: C11  2: C21  3:C22  4:C31, 5:C32, 6:C33
C        ET SI NBCOOR=6 ALORS 7:C41, ... 21:C66
C        REMARQUE: SI PB 2D, C11 a C33 SONT FOURNIS BIEN QUE C31 C32 C33
C                  NE SONT PAS UTILISES ENSUITE!...
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET Texas A & M University           Juillet 2005
C23456---------------------------------------------------------------012
      include"./incl/donthe.inc"
      include"./incl/a___conductivite.inc"
      include"./incl/ctemps.inc"
      include"./incl/cthet.inc"
      include"./incl/cnonlin.inc"

      DOUBLE PRECISION XYZPI(NBCOOR), COND(1:21), PXYZ(10)
      DOUBLE PRECISION CONDUC
      include"./incl/pp.inc"
      COMMON           MCN(MOTMCN)
      REAL             RMCN(1)
      EQUIVALENCE     (MCN(1),RMCN(1))

C     LE TYPE DES DONNEES DE LA CONDUCTIVITE
      LTCOND = MCN( MNCOND + WTCOND )

C     REMPLISSAGE SELON LE TYPE
      IF( LTCOND .EQ. 10 ) THEN

C        MATERIAU ANISOTROPE CONSTANT
C        ============================
         DO I=1,MCN(MNCOND+WVCOND)
            COND(I) = RMCN( MNCOND + WECOND + I - 1 )
         ENDDO
         RETURN

      ELSE IF( LTCOND .GE. 1 .AND. LTCOND .LE. 3 ) THEN

C        MATERIAU HOMOGENE ISOTROPE
C        ==========================
         CONDUC = RMCN( MNCOND + WONDUC )

      ELSE IF( LTCOND .EQ. -1 ) THEN

C        FONCTION UTILISATEUR POUR CONDUCTIVITE ISOTROPE
C        F(temps, x,y,z, no_type_object, materiau, temperature)
C        ======================================================
         PXYZ(1) = TEMPS
         N = 1
         DO I=1,NBCOOR
            N = N + 1
            PXYZ( N ) = XYZPI(I)
         ENDDO
         PXYZ(N+1) = NYOBJT
         PXYZ(N+2) = NUOBJT
         N = N + 3
         PXYZ(N) = TEMPEL
         CALL FONVAL( MCN(MNCOND+WFCOND), N, PXYZ,
     %                NCODEV, CONDUC )
C        NCODEV: 0 CONDUC N'EST PAS INITIALISEE EN SORTIE
C                1 CONDUC   EST     INITIALISEE EN SORTIE
C
      ELSE IF( LTCOND .EQ. -2 ) THEN

C        FONCTION UTILISATEUR POUR CONDUCTIVITE ANISOTROPE
C        =================================================
C        NOMBRE DE COMPOSANTES DU TENSOR
         NVCOAN = MCN(MNCOND+WVCOAN)
         PXYZ(1) = TEMPS
         N = 1
         DO I=1,NBCOOR
            N = N + 1
            PXYZ( N ) = XYZPI(I)
         ENDDO
         PXYZ(N+1) = NYOBJT
         PXYZ(N+2) = NUOBJT
         N = N + 3
         PXYZ(N) = TEMPEL
         N = N + 1
         DO K=1,NVCOAN
            COND(K) = 0D0
C           NO DU COEFFICIENT
C           1:C11, 2:C12=C21, 3:C22, 4:C13=C31, 5:C23=C32, 6:C33
            PXYZ(N) = K
            CALL FONVAL( MCN(MNCOND+WFCOAN), N, PXYZ,
     %                   NCODEV, COND(K) )
            IF( NCODEV .EQ. 0 ) THEN
C              NCODEV: 0 COND(K) N'EST PAS INITIALISEE EN SORTIE
C                      1 COND(K)   EST     INITIALISEE EN SORTIE
C              ERREUR DE DONNEES
               CONDUC = 0D0
               GOTO 50
            ENDIF
         ENDDO
         RETURN

      ELSE

C        ERREUR DE DONNEES
         CONDUC = 0.D0

      ENDIF

C     CREATION DU TENSEUR DE CONDUCTIVITE ISOTROPE HOMOGENE
C     =====================================================
C     MATERIAU 3D ou 2D HOMOGENE ISOTROPE      => C11 a C33
C     MATERIAU AXISYMETRIQUE HOMOGENE ISOTROPE => C11 a C33
C     MATERIAU 6D HOMOGENE ISOTROPE            => C11 a C66

C     LES 3 PREMIERS COEFFICIENTS DIAGONAUX
 50   COND( 1 ) = CONDUC
      COND( 3 ) = CONDUC
      COND( 6 ) = CONDUC
C     LES 3 COEFFICIENTS NON DIAGONAUX
      COND( 2 ) = 0.D0
      COND( 4 ) = 0.D0
      COND( 5 ) = 0.D0
C
      IF( NBCOOR .EQ. 6 ) THEN
         N = 6
         DO I=4,6
            DO J=1,I-1
               N = N + 1
               COND( N ) = 0D0
            ENDDO
            N = N+ 1
            COND( N ) = CONDUC
         ENDDO
      ENDIF

      RETURN
      END
