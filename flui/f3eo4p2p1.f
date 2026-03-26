      SUBROUTINE F3EO4P2P1( Omega,  XYZEF, NONOTE, TP2P2,
     %                      NBNOVI, WM,    WTN,    VE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DU SECOND MEMBRE DU TETRAEDRE TAYLOR-HOOD DU A LA ROTATION
C -----    P2 CONTINU POUR LES 3 COMPOSANTES DE LA VITESSE
C          Integrale tP2 [ - 2 Omega(tn+1) x (Wn+1,m - W(tn)) ] dX
C          FORCE DE CORIOLIS
C
C ENTREES:
C --------
C Omega  : VITESSE ANGULAIRE SUPPOSEE CONSTANTE DE LA ROTATION DU FLUIDE
C XYZEF  : 3 COORDONNEES DES 10 NOEUDS DU TETRAEDRE
C NONOTE : NONOTE(I) NO GLOBAL DU I-EME NOEUD DU TETRAEDRE I=1,...,10
C TP2P2  : Integrale P2i P2j dX SUR LE TETRAEDRE UNITE
C NBNOVI : NOMBRE DE NOEUDS VITESSE SOMMETS et BARYCENTRES DES TETRAEDRES
C WM     : 3 COMPOSANTES DE LA VITESSE WM    EN LES NBNOVI NOEUDS VITESSE
C WTN    : 3 COMPOSANTES DE LA VITESSE W(tn) EN LES NBNOVI NOEUDS VITESSE
C
C MODIFIE:
C --------
C VE     : SECOND MEMBRE ELEMENTAIRE DU A LA ROTATION
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC & St Pierre du Perray      Juin 2011
C23456---------------------------------------------------------------012
      IMPLICIT   NONE
      DOUBLE PRECISION  Omega(3), WM(NBNOVI,3), WTN(NBNOVI,3),
     %                  TP2P2(10,10)
      REAL              XYZEF(10,3)
      INTEGER           NBNOVI
      INTEGER           NONOTE(10), I, J, K, K1, K2, NOEUDJ
      DOUBLE PRECISION  VE(10,3)
      DOUBLE PRECISION  DELTAe, DF(3,3)
      DOUBLE PRECISION  X1, Y1, Z1
      INTRINSIC         ABS
C
C     CALCUL DE LA MATRICE JACOBIENNE ET DE SON DETERMINANT
C     -----------------------------------------------------
      X1 = XYZEF(1,1)
      Y1 = XYZEF(1,2)
      Z1 = XYZEF(1,3)
      DF(1,1) = XYZEF(2,1) - X1
      DF(1,2) = XYZEF(2,2) - Y1
      DF(1,3) = XYZEF(2,3) - Z1
C
      DF(2,1) = XYZEF(3,1) - X1
      DF(2,2) = XYZEF(3,2) - Y1
      DF(2,3) = XYZEF(3,3) - Z1
C
      DF(3,1) = XYZEF(4,1) - X1
      DF(3,2) = XYZEF(4,2) - Y1
      DF(3,3) = XYZEF(4,3) - Z1
C
C     LE DETERMINANT DE DF
      DELTAe= ABS( DF(1,1) * ( DF(2,2) * DF(3,3) - DF(3,2) * DF(2,3) )
     %           + DF(2,1) * ( DF(3,2) * DF(1,3) - DF(1,2) * DF(3,3) )
     %           + DF(3,1) * ( DF(1,2) * DF(2,3) - DF(2,2) * DF(1,3) ) )
C
C     COEFFICIENTS DU SECOND MEMBRE LIE AUX EFFORTS DE ROTATION
C     ---------------------------------------------------------
      DO K=1,3
C
C        COMPOSANTE K DE LA VITESSE
         IF( K .NE. 3 ) THEN
            K1 = K + 1
         ELSE
            K1 = 1
         ENDIF
C
         IF( K1 .NE. 3 ) THEN
            K2 = K1 + 1
         ELSE
            K2 = 1
         ENDIF
C
         DO I=1,10
C
C           POLYNOME P2 TERME I
C
            DO J=1,10
C              POLYNOME P2 TERME J
C              NUMERO DU NOEUD J DU TETRAEDRE
               NOEUDJ = NONOTE(J)
C
C              TP2P2(i,j) = integrale P2i P2j dX
               VE(I,K) = VE(I,K) - 2D0 * TP2P2(I,J) * DELTAe *
     %                 ( Omega(K1) * (WM(NOEUDJ,K2)-WTN(NOEUDJ,K2))
     %                 - Omega(K2) * (WM(NOEUDJ,K1)-WTN(NOEUDJ,K1)) )
C
            ENDDO
C
         ENDDO
C
      ENDDO
C
      RETURN
      END
