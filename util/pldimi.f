      SUBROUTINE PLDIMI( NBPT, PT, PLAN, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LES 4 COEFFICIENTS DU PLAN AX+BY+CZ+D
C -----    QUI MINIMISE LA SOMME DES CARRES DES DISTANCES
C          AXi+BYi+D - Zi  D'UN NUAGE DE POINTS DE R3
C
C ENTREES:
C --------
C NBPT   : LE NOMBRE DE POINTS DU NUAGE
C PT     : LES 3 COORDONNEES DES NBPT POINTS
C
C SORTIES:
C --------
C PLAN   : LES 4 COEFFICIENTS A,B,C,D DU PLAN AX + BY + CZ + D = 0
C IERR   : 0 SI PAS D'ERREUR, 1 SINON
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS        JUIN 1998
C2345X7..............................................................012
      REAL              PT(3,NBPT)
      DOUBLE PRECISION  PLAN(4)
      DOUBLE PRECISION  SYST(3,4)
      INTEGER           NBPT,I,J,IERR
C
      IERR = 0
      IF( NBPT .LE. 2 ) THEN
C        PAS ASSEZ DE POINTS POUR DEFINIR UN PLAN
         IERR = 1
         RETURN
      ENDIF
C
      IF( NBPT .EQ. 3 ) THEN
C        LE PLAN FORME DES 3 POINTS ( PASSAGE EN DOUBLE PRECISION)
         DO 4 J=1,3
            DO 2 I=1,3
               SYST(I,J) = PT(I,J)
 2          CONTINUE
 4       CONTINUE
         CALL EQPLAD( SYST, PLAN, IERR )
         RETURN
      ENDIF
C
C     AU MOINS 4 POINTS DE R**3
C     CALCUL DES COEFFICIENTS DU SYSTEME LINEAIRE POUR Z = AX + BY + D
      DO 20 J=1,4
         DO 10 I=1,3
            SYST(I,J) = 0D0
 10      CONTINUE
 20   CONTINUE
C
      DO 30 I=1,NBPT
         SYST(1,1) = SYST(1,1) + PT(1,I) ** 2
         SYST(2,1) = SYST(2,1) + PT(1,I) * PT(2,I)
         SYST(2,2) = SYST(2,2) + PT(2,I) ** 2
         SYST(3,1) = SYST(3,1) + PT(1,I)
         SYST(3,2) = SYST(3,2) + PT(2,I)
         SYST(1,4) = SYST(1,4) + PT(1,I) * PT(3,I)
         SYST(2,4) = SYST(2,4) + PT(2,I) * PT(3,I)
         SYST(3,4) = SYST(3,4) + PT(3,I)
 30   CONTINUE
      SYST(3,3) = NBPT
C     SYMETRISATION
      SYST(1,2) = SYST(2,1)
      SYST(1,3) = SYST(3,1)
      SYST(2,3) = SYST(3,2)
C
C     RESOLUTION DU SYSTEME LINEAIRE
      CALL GAUSPT( 3, 1 , SYST , IERR )
C
      IF( IERR .EQ. 0 ) THEN
C
C        PLAN Z = A X + B Y + D  =>  PLAN A X + B Y - 1 Z + D = 0
         PLAN(1) = SYST(1,4)
         PLAN(2) = SYST(2,4)
         PLAN(3) = -1D0
         PLAN(4) = SYST(3,4)
         RETURN
      ENDIF
C
C     LE NUAGE DE POINTS EST DANS UN PLAN ORTHOGONAL AU PLAN Z=Cte
C     CALCUL DES COEFFICIENTS DU SYSTEME LINEAIRE POUR X = AY + BZ + C
      DO 50 J=1,4
         DO 40 I=1,3
            SYST(I,J) = 0D0
 40      CONTINUE
 50   CONTINUE
C
      DO 60 I=1,NBPT
         SYST(1,1) = SYST(1,1) + PT(2,I) ** 2
         SYST(2,1) = SYST(2,1) + PT(2,I) * PT(3,I)
         SYST(2,2) = SYST(2,2) + PT(3,I) ** 2
         SYST(3,1) = SYST(3,1) + PT(2,I)
         SYST(3,2) = SYST(3,2) + PT(3,I)
         SYST(1,4) = SYST(1,4) + PT(2,I) * PT(1,I)
         SYST(2,4) = SYST(2,4) + PT(3,I) * PT(1,I)
         SYST(3,4) = SYST(3,4) + PT(1,I)
 60   CONTINUE
      SYST(3,3) = NBPT
C     SYMETRISATION
      SYST(1,2) = SYST(2,1)
      SYST(1,3) = SYST(3,1)
      SYST(2,3) = SYST(3,2)
C
C     RESOLUTION DU SYSTEME LINEAIRE
      CALL GAUSPT( 3, 1 , SYST , IERR )
C
      IF( IERR .EQ. 0 ) THEN
C
C        PLAN X = A Y + B Z + D  =>  PLAN -1 X + A Y + B Z + D = 0
         PLAN(1) = -1D0
         PLAN(2) = SYST(1,4)
         PLAN(3) = SYST(2,4)
         PLAN(4) = SYST(3,4)
         RETURN
      ENDIF
C
C     LE NUAGE DE POINTS EST DANS UN PLAN ORTHOGONAL AU PLAN Z=Cte et X=Cte
C     => DANS UN PLAN Y=Cte  =>  PLAN 0 X - 1 Y + 0 Z + D = 0
      PLAN(1) =  0D0
      PLAN(2) = -1D0
      PLAN(3) =  0D0
      PLAN(4) =  0D0
      DO 70 I=1,NBPT
         PLAN(4) = PLAN(4) + PT(2,I)
 70   CONTINUE
      PLAN(4) = PLAN(4) / NBPT
      IERR = 0
      RETURN
      END
