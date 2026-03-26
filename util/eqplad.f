      SUBROUTINE EQPLAD( COORPO , A , IERR )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT : GENERER LES COEFFICIENTS A DE L EQUATION DU PLAN PASSANT PAR
C ----- 3 POINTS DE COORDONNEES COORPO
C
C PARAMETRE D ENTREE :
C --------------------
C COORPO : COORPO(I,J)=I-EME COORDONNEE DU J-EME POINT
C
C PARAMETRES RESULTATS :
C ----------------------
C A      : COEFFICIENTS DE L EQUATION DU PLAN
C          A(1) * X + A(2) * Y + A(3) * Z + A(4) = 0
C IERR   : 0 PAS D ERREUR
C          1 LES 3 POINTS SONT ALIGNES OU CONFONDUS
C         -1 LES 3 POINTS SONT PROCHES DE L'ALIGNEMENT
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C PROGRAMMEUR : ALAIN PERRONNET ANALYSE NUMERIQUE PARIS       AVRIL 1982
C ......................................................................
      DOUBLE PRECISION  COSMAX
      PARAMETER       ( COSMAX=0.9962D0 )
C     COSMAX : SEUIL DU COSINUS DE L'ANGLE FORME PAR 3 POINTS
C              ( 0.99     => 8.11 DEGRES )
C              ( 0.9962   => 5    DEGRES )
C              ( 0.99756  => 4    DEGRES )
C              ( 0.99863  => 3    DEGRES )
C              ( 0.999    => 2.56 DEGRES )
C              ( 0.9999   => 0.8  DEGRES )
      DOUBLE PRECISION  COORPO(3,3), A(4), COS3PD, D, D1, DMAX, DMAX1
C
C     LES 3 POINTS SONT ILS PROCHES DE L'ALIGNEMENT ?
C     ===============================================
      IERR = 0
      D    = COS3PD( COORPO(1,1), COORPO(1,2), COORPO(1,3) )
      IF( ABS( D ) .GT. 2D0 ) THEN
C        P1 CONFONDU AVEC P2 OU P3
         IERR = 1
         RETURN
      ELSE IF( ABS( D ) .GT. COSMAX ) THEN
         IERR = -1
      ENDIF
C
C     DETERMINATION DU VECTEUR A ORTHOGONAL AU PLAN FORME PAR
C     LES 2 DROITES P1P2  P1P3
C     =======================================================
      IMAX = 0
      DMAX1= -1
      DMAX = -1
      DO 10 I=1,3
         IF( I .NE. 3 ) THEN
            I1 = I + 1
         ELSE
            I1 = 1
         ENDIF
C
C        LE DETERMINANT DES INCONNUES A(I) A(I1)
         D1 = (COORPO(I ,2) - COORPO(I ,1)) *
     &        (COORPO(I1,3) - COORPO(I1,1)) -
     &        (COORPO(I1,2) - COORPO(I1,1)) *
     &        (COORPO(I ,3) - COORPO(I ,1))
C
         D  = ABS( D1 )
         IF( D .GT. DMAX ) THEN
            DMAX  = D
            DMAX1 = D1
            IMAX  = I
         ENDIF
   10 CONTINUE
C
      IF( DMAX .LE. 0D0 ) THEN
C
C        LES 3 POINTS SONT ALIGNES ERREUR
C        ================================
         IERR = 1
C
      ELSE
C
C        LES COEFFICIENTS DU VECTEUR ORTHOGONAL (INCONNUES IMAX,IMAX+1)
C        ==============================================================
         I1 = IMAX + 1
         IF( I1 .GT. 3 ) I1 = 1
         I2 = IMAX + 2
         IF( I2 .GT. 3 ) I2 = I2 - 3
         D = (COORPO(I1,2) - COORPO(I1,1)) *
     &       (COORPO(I2,3) - COORPO(I2,1)) -
     &       (COORPO(I2,2) - COORPO(I2,1)) *
     &       (COORPO(I1,3) - COORPO(I1,1))
         A( IMAX ) = D / DMAX1
         D = (COORPO(I2  ,2) - COORPO(I2  ,1)) *
     &       (COORPO(IMAX,3) - COORPO(IMAX,1)) -
     &       (COORPO(IMAX,2) - COORPO(IMAX,1)) *
     &       (COORPO(I2  ,3) - COORPO(I2  ,1))
         A( I1   ) = D / DMAX1
         A( I2   ) = 1.
C
C        RECHERCHE PARMI TOUS LES PLANS ORTHOGONAUX AU VECTEUR
C        A(1),A(2),A(3) DE CELUI PASSANT PAR LE POINT 1
C        =====================================================
         A(4) = -( A(1) * COORPO(1,1) + A(2) * COORPO(2,1) +
     &             A(3) * COORPO(3,1) )
      ENDIF
      END
