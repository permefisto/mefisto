      SUBROUTINE PLDIM2( ST, NBVECT, VECT, PLAN, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LES 4 COEFFICIENTS DU PLAN AX+BY+CZ+D
C -----    PASSANT PAR LE SOMMET ST ET QUI MINIMISE LA SOMME
C          DES CARRES DES DISTANCES AXi+BYi+D - Zi  D'UN NUAGE
C          DE VECTEURS DE R3
C
C ENTREES:
C --------
C ST     : LE SOMMET IMPOSE DU PLAN
C NBVECT : LE NOMBRE DE POINTS DU NUAGE
C VECT   : LES 3 COORDONNEES DES NBVECT POINTS
C
C SORTIES:
C --------
C PLAN   : LES 4 COEFFICIENTS A,B,C,D DU PLAN AX + BY + CZ + D = 0
C IERR   : 0 SI PAS D'ERREUR, 1 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS        JUIN 1998
C2345X7..............................................................012
      REAL              ST(3), VECT(3,NBVECT)
      DOUBLE PRECISION  PLAN(4)
      DOUBLE PRECISION  SYST3(3,3),S11,S21,S22,B1,B2,D,DNORME
      EQUIVALENCE      (S11,SYST3(1,1)),(S21,SYST3(2,1)),
     %                 (S22,SYST3(2,2)),
     %                 (B1,SYST3(1,3)), (B2,SYST3(2,3))
C
      IF( NBVECT .LE. 1 ) THEN
C        PAS ASSEZ DE POINTS POUR DEFINIR UN PLAN
         IERR = 1
         RETURN
      ENDIF
C
      IF( NBVECT .EQ. 2 ) THEN
C        LE PLAN EST FORME DES 3 POINTS ( PASSAGE EN DOUBLE PRECISION)
         DO 4 J=1,2
            DO 2 I=1,3
               SYST3(I,J) = ST(I) + VECT(I,J)
 2          CONTINUE
 4       CONTINUE
         SYST3(1,3) = ST(1)
         SYST3(2,3) = ST(2)
         SYST3(3,3) = ST(3)
         CALL EQPLAD( SYST3, PLAN, IERR )
         IF( IERR .NE. 0 ) RETURN
         GOTO 9000
      ENDIF
C
C     AU MOINS 3 POINTS DE R**3 + LE SOMMET ST
C     CALCUL DES COEFFICIENTS DU SYSTEME LINEAIRE POUR Z = AX + BY + D
      S11 = 0D0
      S21 = 0D0
      S22 = 0D0
      B1  = 0D0
      B2  = 0D0
C
      DO 30 I=1,NBVECT
         S22 = S22 + VECT(1,I) ** 2
         S21 = S21 - VECT(1,I) * VECT(2,I)
         S11 = S11 + VECT(2,I) ** 2
         B1  = B1  + VECT(1,I) * VECT(3,I)
         B2  = B2  + VECT(2,I) * VECT(3,I)
 30   CONTINUE
C
C     LE DETERMINANT
      D      = S11 * S22 - S21 ** 2
      DNORME = S11**2 + S22**2 + S21**2
      IF( ABS(D) .GT. 1D-5 * DNORME ) THEN
C
C        PLAN  Z = A X + B Y + D  =>  PLAN A X + B Y - 1 Z + D = 0
         PLAN(1) = (S11 * B1 + S21 * B2) / D
         PLAN(2) = (S21 * B1 + S22 * B2) / D
         PLAN(3) = -1D0
         PLAN(4) = ST(3) - PLAN(1) * ST(1) - PLAN(2) * ST(2)
         GOTO 9000
C
      ENDIF
C
C     LE NUAGE DE POINTS EST DANS UN PLAN ORTHOGONAL AU PLAN Z=Cte
C     CALCUL DES COEFFICIENTS DU SYSTEME LINEAIRE POUR X = AY + BZ + C
      S11 = 0D0
      S21 = 0D0
      S22 = 0D0
      B1  = 0D0
      B2  = 0D0
C
      DO 60 I=1,NBVECT
         S22 = S22 + VECT(2,I) ** 2
         S21 = S21 - VECT(2,I) * VECT(3,I)
         S11 = S11 + VECT(3,I) ** 2
         B1  = B1  + VECT(2,I) * VECT(1,I)
         B2  = B2  + VECT(3,I) * VECT(1,I)
 60   CONTINUE
C
C     LE DETERMINANT
      D      = S11 * S22 - S21 ** 2
      DNORME = S11**2 + S22**2 + S21**2
      IF( ABS(D) .GT. 1D-5 * DNORME ) THEN
C
C        PLAN  X = A Y + B Z + D  =>  PLAN -1 X + A Y + B Z + D = 0
         PLAN(1) = -1D0
         PLAN(2) = (S11 * B1 + S21 * B2) / D
         PLAN(3) = (S21 * B1 + S22 * B2) / D
         PLAN(4) = ST(1) - PLAN(2) * ST(2) - PLAN(3) * ST(3)
         GOTO 9000
C
      ENDIF
C
C     LE NUAGE DE POINTS EST DANS UN PLAN ORTHOGONAL AU PLAN Z=Cte et X=Cte
C     => DANS UN PLAN Y=Cte  =>  PLAN  0 X - 1 Y + 0 Z + D = 0
      PLAN(1) =  0D0
      PLAN(2) = -1D0
      PLAN(3) =  0D0
      PLAN(4) =  ST(2)
C
 9000 IERR = 0
      RETURN
      END
