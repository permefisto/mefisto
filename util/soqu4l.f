      SUBROUTINE SOQU4L( LONGCOTE, X3, Y3, X4, Y4, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :  CALCUL DES 2 COORDONNEES DES 2 SOMMETS SUPERIEURS DU QUADRANGLE
C -----  DEFINI PAR LA LONGUEUR DE SES 4 COTES
C
C ENTREES:
C --------
C LONGCOTE: LONGUEUR DES 3 ARETES DES 4 COTES DU QUADRANGLE
C           LONGCOTE(1) >= LONGCOTE(2),LONGCOTE(3),LONGCOTE(4)
C           LONGCOTE(1) <  LONGCOTE(2)+LONGCOTE(3)+LONGCOTE(4)
C
C SORTIES:
C --------
C X3,Y3, X4,Y4 : COORDONNEES DU SOMMET SUPERIEUR DU QUADRANGLE
C                DE SOMMET 1 L'ORIGINE, DE SOMMET 2 (LONGCOTE(1),0)
C IERR   : 0 SI X3,Y3, X4,Y4 SONT CALCULES NORMALEMENT; >0 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET  ANALYSE NUMERIQUE PARIS         OCTOBRE 1993
C MODIFS : ALAIN PERRONNET LABORATOIRE J.-L. LIONS UPMC PARIS  MARS 2006
C234567..............................................................012
      REAL  LONGCOTE( 4 )
C
C     VERIFICATION DE L'EXISTENCE D'UNE SOLUTION
      A = 1E-3 * LONGCOTE(1)
      IF( LONGCOTE(1) .GE. LONGCOTE(2)+LONGCOTE(3)+LONGCOTE(4) ) THEN
C        COTE 1 TROP LONG POUR LES 3 AUTRES COTES
C        LE QUADRANGLE EST LE RECTANGLE DE COTES 1 3 LONGCOTE(1)
C        ET DE COTES 2 4 DE MAX(LONGCOTE(2), LONGCOTE(4))
         IERR = 1
         X3 = LONGCOTE(1)
         X4 = 0D0
         D = MAX( LONGCOTE(2), LONGCOTE(4) )
         Y3 = D
         Y4 = D
         RETURN
      ELSE IF( LONGCOTE(2) .LE. A ) THEN
         IERR = 2
         X3 = LONGCOTE(1)
         Y3 = A
         X4 = 0D0
         Y4 = LONGCOTE(4)
         RETURN
      ELSE IF( LONGCOTE(3) .LE. A ) THEN
         C3 = LONGCOTE(3)
         LONGCOTE(3) = LONGCOTE(4)
         CALL SOTR3L( LONGCOTE, X3, Y3, IERR )
         LONGCOTE(3) = C3
         IF( IERR .EQ. 0 ) THEN
C           TRIANGLE S1 S2 S3
            IERR = 3
            X4 = X3 - A
            Y4 = Y3 - A
         ELSE
C           LE QUADRANGLE EST LE RECTANGLE DE COTES 1 3 LONGCOTE(1)
C           ET DE COTES 2 4 DE MAX(LONGCOTE(2), LONGCOTE(4))
            IERR = 5
            X3 = LONGCOTE(1)
            X4 = 0D0
            D = MAX( LONGCOTE(2), LONGCOTE(4) )
            Y3 = D
            Y4 = D
         ENDIF
         RETURN
      ELSE IF( LONGCOTE(4) .LE. A ) THEN
         IERR = 4
         X3 = LONGCOTE(1)
         Y3 = LONGCOTE(2)
         X4 = 0D0
         Y4 = A
         RETURN
      ENDIF
C
C     RECHERCHE DE 2 POINTS DES CERCLES CENTRES
C     EN S2 DE RAYON L2, EN S1 DE RAYON L4
      C3     = LONGCOTE(3) ** 2
      NBP    = 10
C
C     DEPART AVEC PI/2
 5    A      = ATAN(1.) * 2.
      ANGINC = A / NBP
C
C     ITERATIONS POUR AFFINER L'ANGLE AUX SOMMETS 1 ET 2
      DO 20 ITER = 1,5
C        BOUCLE SUR LES ANGLES
 10      A = A - ANGINC
         IF( A .LE. 0 ) THEN
C           PAS DE SOLUTION AVEC CET INCREMENT D'ANGLE
            NBP = NBP * 2
            IF( NBP .GT. 1024 ) THEN
C              PAS DE CONVERGENCE
               GOTO 80
            ELSE
               GOTO 5
            ENDIF
         ENDIF
C
         COSA = COS( A )
         SINA = SIN( A )
         X4 = LONGCOTE(4) * COSA
         Y4 = LONGCOTE(4) * SINA
         X3 = LONGCOTE(1) - LONGCOTE(2) * COSA
         Y3 =               LONGCOTE(2) * SINA
         D34 = (X4-X3)**2 + (Y4-Y3)**2
         IF( D34 .GT. C3 ) GOTO 10
C
C        PASSAGE AU DELA. RETOUR A L'ANGLE PRECEDENT
         A = A + ANGINC
C        DIMINUTION DE L'INCREMENT DE L'ANGLE
         ANGINC = ANGINC / NBP
 20   CONTINUE
C
      IERR = 0
      RETURN
C
C     PAS DE CONVERGENCE
 80   IERR = 4
      X3 = LONGCOTE(1)
      Y3 = LONGCOTE(2)
      X4 = 0D0
      Y4 = LONGCOTE(4)
      RETURN
      END
