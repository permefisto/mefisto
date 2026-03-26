      SUBROUTINE FL2VITBF3( NBPOE,  X,      VOLTET, NONOEF,
     %                      NTDLVP, NBVECT, VXYZPR, NDDLNV,
     %                      INTVIT, NOFOVI, NBNOVI, VITEXA,  INTVER )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DE L'INTEGRALE L2 DES NBVECT VITESSE
C -----    SUR UN TETRAEDRE DE BREZZI-FORTIN DE VITESSE AVEC 3 COMPOSANTES
C          POLYNOME DE LAGRANGE DE DEGRE 1 + BULLE AU BARYCENTRE
C
C ENTREES:
C --------
C NBPOE  : NOMBRE DE NOEUDS=POINTS DE L'EF
C X      : 3 COORDONNEES DES NBPOE NOEUDS DU TETRAEDRE
C VOLTET : VOLUME DU TETRAEDRE
C NONOEF : NUMEROS DES 4? ou 5? NOEUDS DU TETRAEDRE COURANT
C NTDLVP : NOMBRE TOTAL DE DEGRES DE LIBERTE DU MAILLAGE DE L'OBJET
C NBVECT : NOMBRE DE VECTEURS SOLUTION VITESSEPRESSION
C VXYZPR  : NBVECT VECTEURS VITESSEPRESSION SUR LE MAILLAGE
C NDDLNV : TABLEAU DES POINTEURS SUR LE DERNIER DL POUR CHAQUE NOEUD VITESSE
C          NDDLNV(I)= NO DERNIER DL DU NOEUD (SOMMET ou MILIEU ou BARYCENTRE)
C          CE TABLEAU EST DIMENSIONNE A 0:NBNOVI
C NOFOVI : NUMERO DE LA FONCTION UTILISATEUR VITESSE_EXACTE(t,x,y,z,nocomp)
C          0 SI ELLE N'EXISTE PAS => INTVER N'EST PAS CALCULE
C NBNOVI : NOMBRE DE NOEUDS SUPPORT DE LA VITESSE SUR LE MAILLAGE
C VITEXA : VITEXA(NBNOVI,NBVECT,NDIM) VITESSE EXACTE EN CHAQUE NOEUD TEMPS
C
C SORTIE :
C --------                                                      (Vke1)
C INTVIT : Som  Jacobien Som (Vke1,...,Vke5) [int P1Bi P1Bj DX] (  . )
C          e de E        k=1 ... 5                              (Vke5)
C                                                                     (Vk1E-Vk1C
C INTVER : Som  Som Jacobien (Vk1E-Vk1C,...,Vk5E-Vk5C) [int Pi Pj dX] (   ...
C       e de E  k=1,...,2                                             (Vk5E-Vk5C
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET TEXAS A & M UNIVERSITY at QATAR   Janvier 2012
C23456---------------------------------------------------------------012
C     LES 4 COEFFICIENTS GENERIQUES DE LA MATRICE Integrale tP1B P1B dX
      DOUBLE PRECISION  A, B, C, D
      PARAMETER        (A=29836D0/2494800D0, B=9046D0/2494800D0,
     %                  C=  956D0/ 155925D0, D=4096D0/ 155925D0 )
C
      INTEGER           NONOEF(NBPOE), NDDLNV(0:*)
      REAL              X(NBPOE,3)
C           NBPOE=4? et pour X?
      DOUBLE PRECISION  VXYZPR(NTDLVP,NBVECT),
     %                  VITEXA(NBNOVI,NBVECT,3),
     %                  INTVIT(NBVECT), INTVER(NBVECT),
     %                  VOLTET, INTVEF, INTVERR,
     %                  VitExkI, VitExkJ(3), XBAR, YBAR, ZBAR, DELTA
C
C     INTEGRALES P1B P1B SUR LE TETRAEDRE P1B de REFERENCE
      DOUBLE PRECISION  P1BP1B3D(5,5)
C                       P1BP1B3D(i,j) = integrale P1Bi P1Bj dX
      DATA   P1BP1B3D / A, B, B, B, C,
     %                  B, A, B, B, C,
     %                  B, B, A, B, C,
     %                  B, B, B, A, C,
     %                  C, C, C, C, D /
C
C     JACOBIEN = 6 x (VOLUME DU TETRAEDRE)
      DELTA = 6D0 * VOLTET
C
C     COORDONNEES DU BARYCENTRE
      IF( NOFOVI .GT. 0 ) THEN
         XBAR = ( X(1,1) + X(2,1) + X(3,1) + X(4,1) ) * 0.25D0
         YBAR = ( X(1,2) + X(2,2) + X(3,2) + X(4,2) ) * 0.25D0
         ZBAR = ( X(1,3) + X(2,3) + X(3,3) + X(4,3) ) * 0.25D0
      ENDIF
C
      DO NV=1, NBVECT
C
         INTVEF  = 0D0
         INTVERR = 0D0
C
         DO J=1,5
C
C           NUMERO DU NOEUD J DU TETRAEDRE
            NSJ = NONOEF( J )
C
C           NUMERO DU DL DE VITESSE AU NOEUD AVANT NSJ
            NDLJ = NDDLNV( NSJ-1 )
C
            IF( NOFOVI .GT. 0 ) THEN
C              LA VITESSE_EXACTE(temps,x,y,z,k) AU NOEUD J
               DO K=1,3
                  VitExkJ(K) = VITEXA(NSJ,NV,K)
               ENDDO
            ENDIF
C
            DO I=1,5
C
C              NUMERO DU NOEUD I DU TETRAEDRE
               NSI = NONOEF( I )
C
C              NUMERO DU DL DE VITESSE AU NOEUD AVANT NSI
               NDLI = NDDLNV( NSI-1 )
C
C              NORME L1 DE VITESSE**2 INTEGREE EXACTEMENT
               DO K=1,3
C                 COMPOSANTE K DE LA VITESSE
                  INTVEF = INTVEF + P1BP1B3D(I,J) * VXYZPR(NDLI+K,NV)
     %                                            * VXYZPR(NDLJ+K,NV)
C
                  IF( NOFOVI .GT. 0 ) THEN
C
C                    LA VITESSE_EXACTE(temps,x,y,z,k) AU NOEUD I
                     VitExkI = VITEXA(NSI,NV,K)
C
C                    COMPOSANTE K DE LA VITESSE
                     INTVERR = INTVERR+ P1BP1B3D(I,J)
     %                                * ( VitExkI   -VXYZPR(NDLI+K,NV) )
     %                                * ( VitExkJ(K)-VXYZPR(NDLJ+K,NV) )
                  ENDIF
C
               ENDDO
C
            ENDDO
C
         ENDDO
C
C        CARRE DE LA NORME L2 DE LA VITESSE INTEGREE EXACTEMENT
         INTVIT(NV) = INTVIT(NV) + INTVEF  * DELTA
         INTVER(NV) = INTVER(NV) + INTVERR * DELTA
C
      ENDDO
C
      RETURN
      END
