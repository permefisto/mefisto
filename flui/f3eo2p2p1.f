      SUBROUTINE F3EO2P2P1( DELTAT, Omega0, Omega,  XYZEF,
     %                      NONOTE, NBNOVI, VITES0, TP2P2,
     %                      VE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DU SECOND MEMBRE DU TETRAEDRE TAYLOR-HOOD DU A LA ROTATION
C -----    P2 CONTINU POUR LES 3 COMPOSANTES DE LA VITESSE
C          Integrale tP2 [ - 2 Omega(tn+1) x W(tn)
C                          - ( Omega(tn+1)-Omega(tn) )/dt  x r
C                          -   Omega(tn+1) x ( Omega(tn+1) x r ) ] dX
C          CES FORCES SONT INTERPOLEES P2
C ENTREES:
C --------
C Omega  : VITESSE ANGULAIRE SUPPOSEE CONSTANTE DE LA ROTATION DU FLUIDE
C XYZEF  : 3 COORDONNEES DES 10 NOEUDS DU TETRAEDRE
C NONOTE : NONOTE(I) NO GLOBAL DU I-EME NOEUD DU TETRAEDRE I=1,...,10
C NBNOVI : NOMBRE DE NOEUDS VITESSE SOMMETS et BARYCENTRES DES TETRAEDRES
C VITES0 : DL DES 3 COMPOSANTES DE LA VITESSE A L'INSTANT temps
C TP2P2  : Integrale P2i P2j dX SUR LE TETRAEDRE UNITE
C
C MODIFIE:
C --------
C VE     : SECOND MEMBRE ELEMENTAIRE DU A LA ROTATION
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC & St Pierre du Perray      Juin 2011
C AJOUTg: ALAIN PERRONNET TEXAS A & M University at QATAR      Mars 2012
C23456---------------------------------------------------------------012
      IMPLICIT   NONE
      DOUBLE PRECISION  DELTAT, Omega0(3), Omega(3), VITES0(NBNOVI,3),
     %                  TP2P2(10,10)
      REAL              XYZEF(10,3)
      INTEGER           NBNOVI
      INTEGER           NONOTE(10), I, J, NOEUDJ
      DOUBLE PRECISION  VE(10,3), FORCE(3,10)
      DOUBLE PRECISION  DELTAe, DF(3,3), X1, Y1, Z1
      INTRINSIC         ABS
C
C     CALCUL DE LA MATRICE JACOBIENNE ET DE SON DETERMINANT
C     -----------------------------------------------------
      X1 = XYZEF(1,1)
      Y1 = XYZEF(1,2)
      Z1 = XYZEF(1,3)
C
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
C     LE TETRAEDRE EST SUPPOSE DE VOLUME NON NUL
C
C     CONTRIBUTION VITESSE ANGULAIRE PAR INTERPOLATION DES FORCES
C     DE ROTATION  AUX 10 NOEUDS DU TETRAEDRE
C     -----------------------------------------------------------
      DO J=1,10
C
C        NUMERO DU NOEUD J DU TETRAEDRE
         NOEUDJ = NONOTE(J)
C
C        COORDONNEES DU NOEUD J  et r = i 0 + j y + k z
ccc         X1 = XYZEF(J,1)
         X1 = 0D0
         Y1 = XYZEF(J,2)
         Z1 = XYZEF(J,3)
C
C        CALCUL DES FORCES AU NOEUD J DUES A LA ROTATION
C        Integrale tP2 [ -   Omega(tn+1) x ( Omega(tn+1) x r )
C                        - ( Omega(tn+1)-Omega(tn) ) /dt x r
C                        - 2 Omega(tn+1) x W(tn)  ] dX
cccc
ccccC                    + Grad( -g z k )        ] dX
C
         FORCE(1,J) =
     %             Omega(3) * ( Omega(3) * X1 - Omega(1) * Z1 )
     %           - Omega(2) * ( Omega(1) * Y1 - Omega(2) * X1 )
     %           +(( Omega(3)-Omega0(3) )* Y1
     %           - ( Omega(2)-Omega0(2) )* Z1 ) / DELTAT
     %           + 2D0 * ( Omega(3) * VITES0(NOEUDJ,2)
     %                   - Omega(2) * VITES0(NOEUDJ,3) )
C
         FORCE(2,J) =
     %             Omega(1) * ( Omega(1) * Y1 - Omega(2) * X1 )
     %           - Omega(3) * ( Omega(2) * Z1 - Omega(3) * Y1 )
     %           +(( Omega(1)-Omega0(1) )* Z1
     %           - ( Omega(3)-Omega0(3) )* X1 ) / DELTAT
     %           + 2D0 * ( Omega(1) * VITES0(NOEUDJ,3)
     %                   - Omega(3) * VITES0(NOEUDJ,1) )
C
         FORCE(3,J) =
     %             Omega(2) * ( Omega(2) * Z1 - Omega(3) * Y1 )
     %           - Omega(1) * ( Omega(3) * X1 - Omega(1) * Z1 )
     %           +(( Omega(2)-Omega0(2) )* X1
     %           - ( Omega(1)-Omega0(1) )* Y1 ) / DELTAT
     %           + 2D0 * ( Omega(2) * VITES0(NOEUDJ,1)
     %                   - Omega(1) * VITES0(NOEUDJ,2) )
C
ccc
ccc    Gravite a prendre en compte comme une Force=-9.8 * Rho * Z-Direction
ccc     %           - 9.8D0
C
      ENDDO
C
C     COEFFICIENTS DU SECOND MEMBRE LIE AUX EFFORTS DE ROTATION
C     ---------------------------------------------------------
      DO I=1,10
         DO J=1,10
C           TP2P2(i,j) = integrale P2i P2j dX
            X1 = TP2P2(I,J) * DELTAe
            VE(I,1) = VE(I,1) + X1 * FORCE(1,J)
            VE(I,2) = VE(I,2) + X1 * FORCE(2,J)
            VE(I,3) = VE(I,3) + X1 * FORCE(3,J)
         ENDDO
      ENDDO
C
      RETURN
      END
