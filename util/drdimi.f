      SUBROUTINE DRDIMI( NBPT, PT, DROITE, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCULER LES 3 COEFFICIENTS DE LA DROITE AX+BY+C
C -----    QUI MINIMISE LA SOMME DES CARRES DES DISTANCES
C          AXi+BYi+D - Zi  D'UN NUAGE DE POINTS DE R2
C
C ENTREES:
C --------
C NBPT   : LE NOMBRE DE POINTS DU NUAGE
C PT     : LES 2 COORDONNEES DES NBPT POINTS
C
C SORTIES:
C --------
C DROITE : LES 3 COEFFICIENTS A,B,C DE LA DROITE AX + BY + C = 0
C IERR   : 0 SI PAS D'ERREUR, 1 SI MOINS DE 2 POINTS DANS LE NUAGE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN LJLL UPMC & VEULETTES SUR MER    JUILLET 2008
C2345X7..............................................................012
      REAL              PT(2,NBPT), X21
      DOUBLE PRECISION  DROITE(3), SYST(2,3)
      INTEGER           NBPT, I, IERR
C
      IF( NBPT .LT. 2 ) THEN
C        PAS ASSEZ DE POINTS POUR DEFINIR UNE DROITE
         IERR = 1
         RETURN
      ENDIF
C
      IF( NBPT .EQ. 2 ) THEN
C
C        LE DROITE JOIGNANT LES 2 POINTS ( PASSAGE EN DOUBLE PRECISION)
         X21 = PT(1,2) - PT(1,1)
         IF( X21 .NE. 0.0 ) THEN
C
C           DROITE AX + Y + C = 0
            DROITE(1) = ( PT(2,1) - PT(2,2) ) / X21
            DROITE(2) = 1D0
            DROITE(3) = ( PT(2,2) * PT(1,1) - PT(2,1) * PT(1,2) ) / X21
            GOTO 9900
C
         ELSE
C
C           DROITE X = Cte = X1  =>  DROITE 1 X + 0 Y - X1 = 0
            DROITE(1) =  1D0
            DROITE(2) =  0D0
            DROITE(3) = -PT(1,1)
            GOTO 9900
C
         ENDIF
      ENDIF
C
C     NUAGE AVEC AU MOINS 3 POINTS DE R**2
C     CALCUL DES COEFFICIENTS DU SYSTEME LINEAIRE POUR Y = AX + B
      SYST(1,1) = 0D0
      SYST(2,1) = 0D0
      SYST(1,3) = 0D0
      SYST(2,3) = 0D0
      DO 30 I=1,NBPT
         SYST(1,1) = SYST(1,1) + PT(1,I) ** 2
         SYST(2,1) = SYST(2,1) + PT(1,I)
         SYST(1,3) = SYST(1,3) + PT(1,I) * PT(2,I)
         SYST(2,3) = SYST(2,3) + PT(2,I)
 30   CONTINUE
      SYST(2,2) = NBPT
C     SYMETRISATION
      SYST(1,2) = SYST(2,1)
C
C     RESOLUTION DU SYSTEME LINEAIRE
      CALL GAUSPT( 2, 1, SYST, IERR )
C
      IF( IERR .NE. 0 ) THEN
C
C        DROITE X = Cte = X1  =>  DROITE 1 X + 0 Y - X1 = 0
         DROITE(1) =  1D0
         DROITE(2) =  0D0
         DROITE(3) = -PT(1,1)
C
      ELSE
C
C        DROITE Y=AX+B => DROITE A X - Y + B = 0
         DROITE(1) = SYST(1,3)
         DROITE(2) =-1D0
         DROITE(3) = SYST(2,3)
C
      ENDIF
C
 9900 IERR = 0
      RETURN
      END
