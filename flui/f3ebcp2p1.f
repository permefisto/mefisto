      SUBROUTINE F3EBCP2P1( RhoDt2, Omega, NBNOVI, XYZNOE, NONOTE,
     %                      Wm, Wtn,   BE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DU SECOND MEMBRE DU TETRAEDRE TAYLOR-HOOD DU A LA ROTATION
C -----    P2 CONTINU POUR LES 3 COMPOSANTES DE LA VITESSE
C          Integrale tP2 RhoDt2 [ Omega(tn+1) x (Wn+1,m - W(tn)) ] dX
C          FORCE DE CORIOLIS

C ENTREES:
C --------
C RhoDt2 : -2 * DENSITE VOLUMIQUE DE MASSE DU FLUIDE
C             * PAS DE TEMPS DU SCHEMA EN TEMPS
C Omega  : VITESSE ANGULAIRE SUPPOSEE CONSTANTE DE LA ROTATION DU FLUIDE
C XYZNOE : 3 COORDONNEES DES NOEUDS VITESSE DES TETRAEDRES DU MAILLAGE
C NONOTE : NONOTE(I) NO GLOBAL DU I-EME NOEUD DU TETRAEDRE I=1,...,10
C NBNOVI : NOMBRE DE NOEUDS VITESSE SOMMETS et BARYCENTRES DES TETRAEDRES
C Wm     : 3 COMPOSANTES DE LA VITESSE Wn+1,m EN LES NBNOVI NOEUDS VITESSE
C Wtn    : 3 COMPOSANTES DE LA VITESSE W(tn)  EN LES NBNOVI NOEUDS VITESSE

C MODIFIE:
C --------
C BE     : SECOND MEMBRE ELEMENTAIRE DU A LA ROTATION (FORCE CORIOLIS)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR : ALAIN PERRONNET LJLL UPMC & St PIERRE du PERRAY  Fevrier 2013
C23456---------------------------------------------------------------012
      IMPLICIT   NONE
      include"./incl/p2p23d.inc"
C     DOUBLE PRECISION  P2P23D(10,10) P2P23D(i,j) = integrale P2i P2j dX
      DOUBLE PRECISION  RhoDt2, Omega(3),
     %                  Wm(NBNOVI,3), Wtn(NBNOVI,3)
      REAL              XYZNOE(3,NBNOVI)
      INTEGER           NBNOVI
      INTEGER           NONOTE(10), I, J, K, K1, K2, NOEUDJ
      DOUBLE PRECISION  BE(10,3)
      DOUBLE PRECISION  DELTAe, DF(3,3), X1, Y1, Z1
      INTRINSIC         ABS

C     CALCUL DE LA MATRICE JACOBIENNE ET DE SON DETERMINANT
C     -----------------------------------------------------
      I = NONOTE(1)
      X1 = XYZNOE(1,I)
      Y1 = XYZNOE(2,I)
      Z1 = XYZNOE(3,I)
C
      NOEUDJ = NONOTE(2)
      DF(1,1) = XYZNOE(1,NOEUDJ) - X1
      DF(1,2) = XYZNOE(2,NOEUDJ) - Y1
      DF(1,3) = XYZNOE(3,NOEUDJ) - Z1
C
      NOEUDJ = NONOTE(3)
      DF(2,1) = XYZNOE(1,NOEUDJ) - X1
      DF(2,2) = XYZNOE(2,NOEUDJ) - Y1
      DF(2,3) = XYZNOE(3,NOEUDJ) - Z1
C
      NOEUDJ = NONOTE(4)
      DF(3,1) = XYZNOE(1,NOEUDJ) - X1
      DF(3,2) = XYZNOE(2,NOEUDJ) - Y1
      DF(3,3) = XYZNOE(3,NOEUDJ) - Z1

C     LE DETERMINANT DE DF
      DELTAe= ABS( DF(1,1) * ( DF(2,2) * DF(3,3) - DF(3,2) * DF(2,3) )
     %           + DF(2,1) * ( DF(3,2) * DF(1,3) - DF(1,2) * DF(3,3) )
     %           + DF(3,1) * ( DF(1,2) * DF(2,3) - DF(2,2) * DF(1,3) ) )

C     COEFFICIENTS DU SECOND MEMBRE LIE AUX EFFORTS DE ROTATION
C     ---------------------------------------------------------
      Z1 = RhoDt2 * DELTAe
      DO K=1,3

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

C           POLYNOME P2 TERME I

            DO J=1,10

C              POLYNOME P2 TERME J
C              NUMERO DU NOEUD J DU TETRAEDRE
               NOEUDJ = NONOTE(J)

C              P2P23D(i,j) = integrale P2i P2j dX
               BE(I,K) = BE(I,K) + P2P23D(I,J) * Z1 *
     %                 ( Omega(K1) * (Wm(NOEUDJ,K2)-Wtn(NOEUDJ,K2))
     %                 - Omega(K2) * (Wm(NOEUDJ,K1)-Wtn(NOEUDJ,K1)) )

            ENDDO

         ENDDO

      ENDDO

      RETURN
      END
