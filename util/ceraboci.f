      SUBROUTINE CERABOCI( P1, P2, P3, P4,  CENTRE, VOLUTE, IERR )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :   CALCUL DES 3 COORDONNEES DU CENTRE DE LA BOULE CIRCONSCRITE
C -----   AU TETRAEDRE P1 P2 P3 P4 ET DU CARRE DE SON RAYON
C         CALCUL DU VOLUME ORIENTE DU TETRAEDRE P1P2P3P4
C
C ENTREES:
C --------
C P1 P2 P3 P4 : LES 3 COORDONNEES DES 4 SOMMETS
C
C SORTIE :
C --------
C CENTRE : 3 COORDONNEES DU CENTRE ET CARRE DU RAYON DE LA BOULE CIRCONSCRITE
C          0,0,0,-1 SI TETRAEDRE INCORRECT
C VOLUTE : LE VOLUME DU TETRAEDRE (>0 SI ORIENTE COMME UN REPERE
C                                  =0 SI TETRAEDRE DEGENERE
C                                  <0 SI TETRAEDRE MAL ORIENTE)
C IERR   : =0 SI CALCUL CORRECT SANS PROBLEME RENCONTRE
C          =1 SI TETRAEDRE PLAT (4 SOMMETS DANS UN PLAN)
C          =2 SI TETRAEDRE MAL ORIENTE DE VOLUME STRICTEMENT NEGATIF
C          =3 SI RAYON**2 DE LA BOULE CIRCONSCRITE PAR LA FORMULE est <0
C          =4 SI TROP GRANDE DIFFERENCE SUR LE CALCUL DU RAYON AUX 4 SOMMETS
C                SOUVENT DU A UN VOLUME TROP FAIBLE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC Paris & St PIERRE du PERRAY Mai 2017
C2345X7..............................................................012
      DOUBLE PRECISION  P1(3), P2(3), P3(3), P4(3), CENTRE(4), VOLUTE,
     %                  D, VOLTET, VOLU6T, DP1, DP2, DP3, DP4, GAMA
ccc      DOUBLE PRECISION  C4

      IERR = 0

cccC     VOLU6T EST LE DETERMINANT DE LA MATRICE [P2-P1, P3-P1, P4-P1]
cccC     VOLU6T EST AUSSI -6 FOIS LE VOLUME DU TETRAEDRE
ccc      CALL DETMAT44( P1(1), P1(2), P1(3), 1D0,
ccc     %               P2(1), P2(2), P2(3), 1D0,
ccc     %               P3(1), P3(2), P3(3), 1D0,
ccc     %               P4(1), P4(2), P4(3), 1D0,  VOLU6T )
cccC     LE VOLUME DU TETRAEDRE
ccc      VOLUTE = -VOLU6T / 6D0
ccc      PRINT*,'ceraboci 1: VOLU6T=',VOLU6T,' VOLUTE=',VOLUTE


C     LE VOLUME DU TETRAEDRE
      VOLUTE = VOLTET( P1, P2, P3, P4 )

C     VOLU6T EST LE DETERMINANT DE LA MATRICE [P2-P1, P3-P1, P4-P1]
C     VOLU6T EST AUSSI -6 FOIS LE VOLUME DU TETRAEDRE
      VOLU6T = - VOLUTE * 6D0

ccc      PRINT*,'ceraboci 2: VOLU6T=',VOLU6T,' VOLUTE=',VOLUTE

      IF( VOLUTE .LT. 0.D0 ) THEN

ccc       PRINT 10010, VOLUTE, P1, P2, P3, P4
ccc10010  FORMAT('ceraboci: VOLUME du TETRAEDRE =',G25.16,'<0 de SOMMETS'/
ccc     %         ('  X=',G25.16,'   Y=',G25.16,'   Z=',G25.16))

C        VALEURS IMPOSEES POUR UN TETRAEDRE de VOLUME<0
C        CENTRE et CARRE DU RAYON SONT NULS
         IERR = 2
         GOTO 9000

ccc      ELSE IF( VOLU6T .LE. SQRT(A*B*C)*1D-3 ) THEN   19/8/2014

      ELSE IF( VOLUTE .EQ. 0D0 ) THEN

ccc10020  FORMAT('ceraboci: VOLUME=0 du TETRAEDRE PLAT de SOMMETS'/
ccc     %         ('  X=',G25.16,'   Y=',G25.16,'   Z=',G25.16))
ccc         PRINT 10020, VOLUTE, P1, P2, P3, P4

         IERR = 1
         GOTO 9000

      ENDIF

C     VOLUME POSITF
C     || P1 || ** 2
      DP1 = P1(1) ** 2 + P1(2) ** 2 + P1(3) ** 2

C     || P2 || ** 2
      DP2 = P2(1) ** 2 + P2(2) ** 2 + P2(3) ** 2

C     || P3 || ** 2
      DP3 = P3(1) ** 2 + P3(2) ** 2 + P3(3) ** 2

C     || P4 || ** 2
      DP4 = P4(1) ** 2 + P4(2) ** 2 + P4(3) ** 2

C     COEFFICIENT GAMA
      CALL DETMAT44( DP1, P1(1), P1(2), P1(3),
     %               DP2, P2(1), P2(2), P2(3),
     %               DP3, P3(1), P3(2), P3(3),
     %               DP4, P4(1), P4(2), P4(3),  GAMA )

C     ABSCISSE X DU CENTRE DE LA BOULE CIRCONSCRITE
      CALL DETMAT44( DP1, P1(2), P1(3), 1D0,
     %               DP2, P2(2), P2(3), 1D0,
     %               DP3, P3(2), P3(3), 1D0,
     %               DP4, P4(2), P4(3), 1D0,  CENTRE(1) )

C     -ORDONNEE Y DU CENTRE DE LA BOULE CIRCONSCRITE
      CALL DETMAT44( DP1, P1(1), P1(3), 1D0,
     %               DP2, P2(1), P2(3), 1D0,
     %               DP3, P3(1), P3(3), 1D0,
     %               DP4, P4(1), P4(3), 1D0,  CENTRE(2) )
      CENTRE(2) = -CENTRE(2)

C     COTE Z DU CENTRE DE LA BOULE CIRCONSCRITE
      CALL DETMAT44( DP1, P1(1), P1(2), 1D0,
     %               DP2, P2(1), P2(2), 1D0,
     %               DP3, P3(1), P3(2), 1D0,
     %               DP4, P4(1), P4(2), 1D0,  CENTRE(3) )

cccC     RAYON**2 DE LA BOULE CIRCONSCRITE PAR LA FORMULE
ccc      CENTRE(4) = ( CENTRE(1)**2 + CENTRE(2)**2 + CENTRE(3)**2
ccc     %            - 4D0 * VOLU6T * GAMA ) / ( 4D0 * VOLU6T * VOLU6T )

ccc      IF( CENTRE(4) .LE. 0D0 ) THEN
cccccc         PRINT *,'ceraboci: ERREUR d''ARRONDI? RAYON**2=',CENTRE(4),
cccccc     %           ' EST <0 . VOLUME du tetraedre=',VOLUTE
ccc         IERR = 3
ccc      ENDIF

C     XYZ DU CENTRE DE LA BOULE CIRCONSCRITE
      D = 1D0 / ( 2D0 * VOLU6T )
      CENTRE( 1 ) = CENTRE( 1 ) * D
      CENTRE( 2 ) = CENTRE( 2 ) * D
      CENTRE( 3 ) = CENTRE( 3 ) * D

C     VERIFICATION DE L'EGALITE DES 5 RAYONS**2
      DP1 = ( CENTRE(1) - P1(1) ) ** 2
     %    + ( CENTRE(2) - P1(2) ) ** 2
     %    + ( CENTRE(3) - P1(3) ) ** 2

      DP2 = ( CENTRE(1) - P2(1) ) ** 2
     %    + ( CENTRE(2) - P2(2) ) ** 2
     %    + ( CENTRE(3) - P2(3) ) ** 2

      DP3 = ( CENTRE(1) - P3(1) ) ** 2
     %    + ( CENTRE(2) - P3(2) ) ** 2
     %    + ( CENTRE(3) - P3(3) ) ** 2

      DP4 = ( CENTRE(1) - P4(1) ) ** 2
     %    + ( CENTRE(2) - P4(2) ) ** 2
     %    + ( CENTRE(3) - P4(3) ) ** 2

ccc      IF( IERR .EQ. 0 ) THEN

cccC        CARRE DU RAYON = MOYENNE DES 5 DIFFERENTES VALEURS
ccc         CENTRE(4) = ( CENTRE(4) + DP1 + DP2 + DP3 + DP4 ) * 0.2D0

ccc      ELSE

ccc         C4 = CENTRE(4)

C     FINALEMENT: RAYON**2 FINAL=MOYENNE DES 4 DISTANCES AU CARRE
      CENTRE(4) = ( DP1 + DP2 + DP3 + DP4 ) * 0.25D0

C     VERIFICATION et IMPRESSION
      D = MAX( CENTRE(4), DP1, DP2, DP3, DP4 )
     %  - MIN( CENTRE(4), DP1, DP2, DP3, DP4 )

      IF( D .GT. CENTRE(4) * 5D-3 ) THEN

         IERR = 4
C        PROVOQUE PAR UN VOLUME TROP PETIT DU TETRAEDRE

C        AFFICHAGE POUR VOIR LES DIFFERENCES
         PRINT *
         PRINT *,'ceraboci: Volume du TETRAEDRE =',VOLUTE,
     %           ' XYZ du CENTRE de la BOULE CIRCONSCRITE=',
     %           (CENTRE(I),I=1,3)
         PRINT *,'RAYON**2 FINAL    =',CENTRE(4)
         PRINT *,'||CENTRE - P1||**2=',DP1
         PRINT *,'||CENTRE - P2||**2=',DP2
         PRINT *,'||CENTRE - P3||**2=',DP3
         PRINT *,'||CENTRE - P4||**2=',DP4
         PRINT *,'ECART MAX-MIN RELATIF des RAYONS**2=',D/CENTRE(4),
     %           ' IERR=',IERR
ccc            PRINT *,'RAYON**2  FORMULE =',C4
         PRINT *
         GOTO 9000

      ELSE
         IERR = 0
      ENDIF

ccc      ENDIF

      GOTO 9999


C     VALEURS IMPOSEES POUR UN TETRAEDRE de VOLUME<=0 ou PIRE
C     CENTRE XYZ NULS et CARRE DU RAYON =-1 POUR EVITER DIVISION/R=0D0
 9000 CENTRE(1) = 0D0
      CENTRE(2) = 0D0
      CENTRE(3) = 0D0
      CENTRE(4) =-1D0


 9999 RETURN
      END
