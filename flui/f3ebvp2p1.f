      SUBROUTINE F3EBVP2P1( Rho,    DELTAT, Omega0, Omega,
     %                      NBNOVI, XYZNOE, NONOTE, NUDDL, Wtn,
     %                      VE )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DU SECOND MEMBRE DU TETRAEDRE TAYLOR-HOOD DU A LA ROTATION
C -----    P2 CONTINU POUR LES 3 COMPOSANTES DE LA VITESSE
C          Integrale tP2 Rho dt [ - 2 Omega(tn+1) x W(tn)
C                                 - ( Omega(tn+1)-Omega(tn) )/dt  x r
C                                 -   Omega(tn+1) x ( Omega(tn+1) x r ) ] dX
C          CES FORCES SONT INTERPOLEES P2
C ENTREES:
C --------
C Rho    : DENSITE VOLUMIQUE DE MASSE DU FLUIDE
C DELTAT : PAS DE TEMPS DU SCHEMA D'INTEGRATION EN TEMPS
C Omega0 : VITESSE ANGULAIRE SUPPOSEE CONSTANTE DE LA ROTATION DU FLUIDE a tn
C Omega  : VITESSE ANGULAIRE SUPPOSEE CONSTANTE DE LA ROTATION DU FLUIDE a tn+1
C NBNOVI : NOMBRE DE NOEUDS VITESSE DES TETRAEDRES DU MAILLAGE
C XYZNOE : 3 COORDONNEES DES NBNOVI NOEUDS VITESSE DU MAILLAGE
C NONOTE : NONOTE(I) NO GLOBAL DU I-EME NOEUD DU TETRAEDRE I=1,...,10
C NUDDL  : NUMERO DU DERNIER DL DE CHAQUE NOEUD VITESSE (0:NBNOVI)
C          A UTILISER POUR RETROUVER LES DL PRESSION DU VECTEUR Wtn
C Wtn    : DL DES 3 COMPOSANTES DE LA VITESSE W(tn) VXVYVZPR(Noeud), ...
C
C MODIFIE:
C --------
C VE     : SECOND MEMBRE ELEMENTAIRE DU A LA ROTATION DU DOMAINE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET LJLL UPMC & St Pierre du Perray   Fevrier 2013
C23456---------------------------------------------------------------012
      IMPLICIT   NONE
      include "./incl/p2p23d.inc"
C     INTEGRALES P2 P2 SUR LE TETRAEDRE P2 de REFERENCE
C     DOUBLE PRECISION  P2P23D(10,10)
C                       P2P23D(i,j) = integrale P2i P2j dX

      DOUBLE PRECISION  Rho, DELTAT, Omega0(3), Omega(3), Wtn(*)
      REAL              XYZNOE(3,NBNOVI)
      INTEGER           NBNOVI, NUDDL(0:NBNOVI)
      INTEGER           NONOTE(10), I, J, NOEUDJ, NDL
      DOUBLE PRECISION  VE(10,3), FORCE(3,10), RhoDt
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
C        COORDONNEES DU NOEUD J telles que r = i X + j Y + k Z
         X1 = XYZNOE(1,NOEUDJ)
         Y1 = XYZNOE(2,NOEUDJ)
         Z1 = XYZNOE(3,NOEUDJ)
C
C        CALCUL DES FORCES AU NOEUD J DUES A LA ROTATION
C        Integrale tP2 [- 2 Omega(tn+1) x W(tn) 
C                       -   Omega(tn+1) x ( Omega(tn+1) x r )
C                       - ( Omega(tn+1)-Omega(tn) ) /dt x r   ] dX
ccccC                   + dt Grad( -Rho g z k )    ] dX

C        NUMERO DU DL PRECEDANT LES DL DU NOEUDJ
         NDL = NUDDL( NOEUDJ - 1 )

         FORCE(1,J) = 2D0 * ( Omega(3) * Wtn(NDL+2)
     %                      - Omega(2) * Wtn(NDL+3) )
     %              + Omega(3) * ( Omega(3) * X1 - Omega(1) * Z1 )
     %              - Omega(2) * ( Omega(1) * Y1 - Omega(2) * X1 )
     %              +(( Omega(3)-Omega0(3) )* Y1
     %              - ( Omega(2)-Omega0(2) )* Z1 ) / DELTAT

         FORCE(2,J) = 2D0 * ( Omega(1) * Wtn(NDL+3)
     %                      - Omega(3) * Wtn(NDL+1) )
     %              + Omega(1) * ( Omega(1) * Y1 - Omega(2) * X1 )
     %              - Omega(3) * ( Omega(2) * Z1 - Omega(3) * Y1 )
     %              +(( Omega(1)-Omega0(1) )* Z1
     %              - ( Omega(3)-Omega0(3) )* X1 ) / DELTAT

         FORCE(3,J) = 2D0 * ( Omega(2) * Wtn(NDL+1)
     %                      - Omega(1) * Wtn(NDL+2) )
     %              + Omega(2) * ( Omega(2) * Z1 - Omega(3) * Y1 )
     %              - Omega(1) * ( Omega(3) * X1 - Omega(1) * Z1 )
     %              +(( Omega(2)-Omega0(2) )* X1
     %              - ( Omega(1)-Omega0(1) )* Y1 ) / DELTAT

ccc
ccc    Gravite a prendre en compte comme une Force=-9.8 * Rho * Z-Direction
ccc     %           - 9.8D0 * Rho * DELTAT
C
      ENDDO
C
C     COEFFICIENTS DU SECOND MEMBRE LIE AUX EFFORTS DE ROTATION
C     ---------------------------------------------------------
      RhoDt = Rho * DELTAT * DELTAe
      DO I=1,10
         DO J=1,10
C           P2P23D(i,j) = integrale P2i P2j dX
            X1 = P2P23D(I,J) * RhoDt
            VE(I,1) = VE(I,1) + X1 * FORCE(1,J)
            VE(I,2) = VE(I,2) + X1 * FORCE(2,J)
            VE(I,3) = VE(I,3) + X1 * FORCE(3,J)
         ENDDO
      ENDDO
C
      RETURN
      END
