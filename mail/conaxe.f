      SUBROUTINE CONAXE( PAX1, PAX2, RPCONE, HAUTEU )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DES 3 VECTEURS UNITAIRES D'UN TRONC DE CONE
C -----    A SECTIONS ORTHOGONALES AUX POINTS PAX1 ET PAX2 DE L'AXE
C          LES COMPOSANTES DES 3 VECTEURS SONT EXPRIMEES DANS LE
C          REPERE GLOBAL ANCIEN
C
C ENTREES:
C --------
C PAX1   : POINT DE L'AXE DU CONE
C PAX2   : POINT DE L'AXE DU CONE (DIFFERENT DE AX1)
C
C SORTIES:
C --------
C RPCONE : LES 3 VECTEURS UNITAIRES DU REPERE LIE AU CONE
C          OU MATRICE DE PASSAGE DES COORDONNEES D'UN POINT
C          XX, YY, ZZ DANS LE REPERE DU CONE
C          X , Y , Z  DANS LE REPERE GLOBAL
C
C  (X       )    (        ) ( XX )    ( XX ) ( T        )(X       )
C  (Y - PAX1) =  ( RPCONE ) ( YY ) ;  ( YY )=(   RPCONE )(Y - PAX1)
C  (Z       )    (        ) ( ZZ )    ( ZZ ) (          )(Z       )
C
C HAUTEU : DISTANCE ENTRE LES POINTS AX1 ET AX2 (0 => ERREUR)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : PERRONNET ALAIN UPMC ANALYSE NUMERIQUE PARIS    DECEMBRE 1997
C2345X7..............................................................012
      DOUBLE PRECISION  PAX1(3), PAX2(3)
      DOUBLE PRECISION  RPCONE(3,3), HAUTEU
C
C     L'AXE 3 "VERTICAL" EST PORTE PAR LA DIRECTION AX1 VERS AX2
      HAUTEU = 0D0
      DO 10 I=1,3
         RPCONE(I,3) = PAX2(I) - PAX1(I)
         HAUTEU = HAUTEU + RPCONE(I,3) ** 2
 10   CONTINUE
      IF( HAUTEU .LE. 0D0 ) THEN
         RETURN
      ENDIF
      HAUTEU = SQRT( HAUTEU )
      DO 20 I=1,3
         RPCONE(I,3) = RPCONE(I,3) / HAUTEU
 20   CONTINUE
C
C     LE PREMIER AXE EST X SAUF SI LE VECTEUR VERTICAL EST DEJA L'ANCIEN AXE X
      IF( RPCONE(2,3) .NE. 0D0 .OR. RPCONE(3,3) .NE. 0D0 ) THEN
C
C        L'AXE 1 EST L'ANCIEN AXE X
         RPCONE(1,1) = 1D0
         RPCONE(2,1) = 0D0
         RPCONE(3,1) = 0D0
C
C        L'AXE 2 EST AXE3 VECTORIEL AXE1
         CALL PROVEC( RPCONE(1,3), RPCONE(1,1), RPCONE(1,2) )
         CALL NORME1( 3, RPCONE(1,2), I )
C
C        L'AXE 1 EST ORTHONORMALISE
         CALL PROVEC( RPCONE(1,2), RPCONE(1,3), RPCONE(1,1) )
         CALL NORME1( 3, RPCONE(1,1), I )
C
      ELSE
C
C        L'AXE 1 EST Y , L'AXE 2 EST Z , L'AXE 3 EST X
         RPCONE(1,1) = 0D0
         RPCONE(2,1) = 1D0
         RPCONE(3,1) = 0D0
C
         RPCONE(1,2) = 0D0
         RPCONE(2,2) = 0D0
         RPCONE(3,2) = 1D0
C
         RPCONE(1,3) = 1D0
         RPCONE(2,3) = 0D0
         RPCONE(3,3) = 0D0
C
      ENDIF
      RETURN
      END
