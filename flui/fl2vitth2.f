      SUBROUTINE FL2VITTH2( NBPOE,  SURFTR, NONOEF,
     %                      NTDLVP, NBVECT, VXYPR,  NDDLNV,
     %                      INTVIT, NOFOVI, NBNOVI, VITEXA,  INTVER )
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C BUT :    CALCUL DE L'INTEGRALE DES NBVECT  VITESSE**2
C -----    SUR UN TRIANGLE DE TAYLOR-HOOD DE VITESSE AVEC 2 COMPOSANTES
C          POLYNOME DE LAGRANGE DE DEGRE 2
C
C ENTREES:
C --------
C NBPOE  : NOMBRE DE NOEUDS=POINTS DE L'EF
C SURFTR : SURFACE DU TRIANGLE
C NONOEF : NUMEROS DES NBPOE NOEUDS DU TRIANGLE COURANT
C NTDLVP : NOMBRE TOTAL DE DEGRES DE LIBERTE DU MAILLAGE DE L'OBJET
C NBVECT : NOMBRE DE VECTEURS SOLUTION VITESSEPRESSION
C VXYPR  : NBVECT VECTEURS VITESSEPRESSION SUR LE MAILLAGE
C NDDLNV : TABLEAU DES POINTEURS SUR LE DERNIER DL POUR CHAQUE NOEUD VITESSE
C          NDDLNV(I)= NO DERNIER DL DU NOEUD (SOMMET ou MILIEU ou BARYCENTRE)
C          CE TABLEAU EST DIMENSIONNE A 0:NBNOVI
C NOFOVI : NUMERO DE LA FONCTION UTILISATEUR VITESSE_EXACTE(t,x,y,z,nocomp)
C          0 SI ELLE N'EXISTE PAS => INTVER N'EST PAS CALCULE
C NBNOVI : NOMBRE DE NOEUDS SUPPORT DE LA VITESSE SUR LE MAILLAGE
C VITEXA : VITEXA(NBNOVI,NBVECT,NDIM) VITESSE EXACTE EN CHAQUE NOEUD TEMPS
C
C SORTIE :
C --------                                                    (Vke1)
C INTVIT : Som  Jacobien Som (Vke1,...,Vke6) [int P2i P2j DX] (  . )
C          e de E        k=1 ... 6                            (Vke6)
C                                                                     (Vk1E-Vk1C
C INTVER : Som  Som Jacobien (Vk1E-Vk1C,...,Vk6E-Vk6C) [int Pi Pj dX] (   ...
C       e de E  k=1,...,2                                             (Vk6E-Vk6C
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C AUTEUR: ALAIN PERRONNET TEXAS A & M UNIVERSITY at QATAR   Janvier 2012
C23456---------------------------------------------------------------012
C     INTEGRALES P2 P2 SUR LE TRIANGLE P2 de REFERENCE
      include"./incl/p2p22d.inc"
CCC   DOUBLE PRECISION  P2P22D(6,6)
CCCC                    P2P22D(i,j) = integrale P2i P2j dX
C
      INTEGER           NONOEF(NBPOE), NDDLNV(0:NBNOVI)
      DOUBLE PRECISION  VXYPR(NTDLVP,NBVECT),
     %                  VITEXA(NBNOVI,NBVECT,2),
     %                  INTVIT(NBVECT), INTVER(NBVECT),
     %                  SURFTR, INTVEF, INTVERR, DELTA,
     %                  VitExkI, VitExkJ(2)
C
C     JACOBIEN = 2 x (SURFACE DU TRIANGLE)
      DELTA = 2D0 * SURFTR
C
      DO NV=1, NBVECT
C
         INTVEF  = 0D0
         INTVERR = 0D0
C
         DO J=1,6
C
C           NUMERO DU NOEUD J DU TRIANGLE
            NSJ = NONOEF( J )
C
C           NUMERO DU DL DE VITESSE AU NOEUD AVANT NSJ
            NDLJ = NDDLNV( NSJ-1 )
C
            IF( NOFOVI .GT. 0 ) THEN
C              LA VITESSE_EXACTE(temps,x,y,z,k) AU NOEUD J
               DO K=1,2
                  VitExkJ(K) = VITEXA(NSJ,NV,K)
               ENDDO
            ENDIF
C
            DO I=1,6
C
C              NUMERO DU NOEUD I DU TRIANGLE
               NSI = NONOEF( I )
C
C              NUMERO DU DL DE VITESSE AU NOEUD AVANT NSI
               NDLI = NDDLNV( NSI-1 )
C
C              CARRE DE LA NORME L2 DE |VITESSE| INTEGREE EXACTEMENT
               DO K=1,2
C
C                 COMPOSANTE K DE LA VITESSE
                  INTVEF = INTVEF + P2P22D(I,J) * VXYPR(NDLI+K,NV)
     %                                          * VXYPR(NDLJ+K,NV)
                  IF( NOFOVI .GT. 0 ) THEN
C
C                    LA VITESSE_EXACTE(temps,x,y,z,k) AU NOEUD I
                     VitExkI = VITEXA(NSI,NV,K)
C
C                    COMPOSANTE K DE LA VITESSE
                     INTVERR = INTVERR + P2P22D(I,J)
     %                                 * ( VitExkI   -VXYPR(NDLI+K,NV) )
     %                                 * ( VitExkJ(K)-VXYPR(NDLJ+K,NV) )
                  ENDIF
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
